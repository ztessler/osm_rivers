import csv
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import palettable
import geopandas
import rasterio
import rasterio.features as rfeatures
import cartopy.crs as ccrs
import networkx as nx
import skimage.morphology as morph
from scipy.ndimage.filters import generic_filter
from netCDF4 import Dataset
import pyproj
import itertools


def project_and_clip_osm_rivers(source, target, env):
    rivers_ll = geopandas.read_file(str(source[0]))
    delta_ll = geopandas.read_file(str(source[1]))

    deltahull_ll = geopandas.GeoDataFrame(delta_ll.dissolve(by='Delta').convex_hull, columns=['geometry'])
    deltahull_ll.crs = delta_ll.crs
    lon0, lat0 = np.array(deltahull_ll.centroid.squeeze())

    laea = ccrs.LambertAzimuthalEqualArea(central_longitude = lon0,
                                          central_latitude = lat0)

    rivers = rivers_ll.to_crs(laea.proj4_params)
    delta = delta_ll.to_crs(laea.proj4_params)
    deltahull = deltahull_ll.to_crs(laea.proj4_params)

    rivers_clip = geopandas.overlay(rivers, deltahull, how='intersection') #slow

    rivers_clip.to_file(str(target[0]))
    with open(str(target[1]), 'w') as fout:
        fout.write(laea.proj4_init + '\n')
    return 0


def clip_osm_rivers(source, target, env):
    rivers_ll = geopandas.read_file(str(source[0]))
    delta_ll = geopandas.read_file(str(source[1]))

    deltahull_ll = geopandas.GeoDataFrame(delta_ll.dissolve(by='Delta').convex_hull, columns=['geometry'])
    deltahull_ll.crs = delta_ll.crs

    rivers_clip = geopandas.overlay(rivers_ll, deltahull_ll, how='intersection') #slow

    rivers_clip.to_file(str(target[0]))
    return 0


def thin_vec(source, target, env):
    rivers = geopandas.read_file(str(source[0]))
    thresh = env.get('thresh', 100)

    rivers_eroded = rivers.buffer(-thresh)
    rivers_clean = rivers_eroded[rivers_eroded.area>0].buffer(thresh)
    rivers_merge = geopandas.overlay(geopandas.GeoDataFrame(rivers_clean, columns=['geometry']), geopandas.GeoDataFrame(rivers_clean, columns=['geometry']), how='union')

    rivers_merge.to_file(str(target[0]))
    return 0


def plot_vec_rivs(source, target, env):
    rivers = geopandas.read_file(str(source[0]))

    mpl.style.use('ggplot')

    fig, ax = plt.subplots(1, 1, figsize=(8,6))
    rivers.plot(ax=ax)
    ax.set_aspect('equal')
    fig.savefig(str(target[0]))

    return 0


def rasterize_riv(source, target, env):
    rivers = geopandas.read_file(str(source[0]))
    imsize = env.get('imsize', 1000)

    bounds = rivers.bounds
    width = bounds['maxx'].max() - bounds['minx'].min()
    height = bounds['maxy'].max() - bounds['miny'].min()
    dx = (width + height) / 2. / imsize
    #dx = width / imshape[0]
    #dy = height / imshape[1]
    imshape = (int(np.ceil(height / dx)), int(np.ceil(width / dx)))
    buffx = np.ceil(width / dx) - (width / dx)
    buffy = np.ceil(height / dx) - (height / dx)
    x0 = bounds['minx'].min() - (buffx * dx)/2
    x1 = bounds['maxx'].max() + (buffx * dx)/2
    y0 = bounds['miny'].min() - (buffy * dx)/2
    y1 = bounds['maxy'].max() + (buffy * dx)/2


    a=dx; b=0; c=x0; d=0; e=-dx; f=y1
    affine = rasterio.Affine(a,b,c,d,e,f)

    burned = rfeatures.rasterize(
            ((f,1) for f in rivers['geometry']),
            out_shape=imshape,
            transform=affine,
            all_touched=True,
            dtype=np.uint8)

    with rasterio.open(str(target[0]), 'w',
            driver='GTiff',
            width=imshape[1], height=imshape[0],
            count=1, dtype=str(burned.dtype),
            crs=rivers.crs,
            transform=affine,
            nodata=0) as rast:
        rast.write(burned, 1)
    return 0


def georef_nc(env, target, source):
    nc = Dataset(str(source[0]))
    var = nc.variables[nc.subject]
    data = var[:].squeeze().data.astype(np.int32)
    nodata = var.missing_value
    lat_bnds = nc.variables['latitude_bnds'][:]
    lon_bnds = nc.variables['longitude_bnds'][:]

    yoff, xoff = data.shape
    sx = np.diff(lon_bnds).mean()
    sy = np.diff(lat_bnds).mean()

    affine = rasterio.Affine.translation(lon_bnds.min(), lat_bnds.max()) * rasterio.Affine.scale(sx, -sy)
    with rasterio.open(str(target[0]), 'w',
            driver='GTiff', width=xoff, height=yoff,
            crs={'init':'epsg:4326'}, transform=affine,
            count=1, nodata=nodata, dtype=str(data.dtype)) as dst:
        dst.write(np.flipud(data), 1)

    nc.close()
    return 0


def skeleton_riv(source, target, env):
    with rasterio.open(str(source[0])) as rast:
        rivers = rast.read(1)
        meta = rast.meta.copy()

    closing = env.get('closing', 5)
    rivers = morph.binary_closing(rivers, morph.square(closing))

    holethresh = env.get('holethresh', 1000)
    rivers = morph.remove_small_holes(rivers, min_size=holethresh, connectivity=2)

    rivers, distance = morph.medial_axis(rivers, return_distance=True)
    rivers = rivers.astype(np.uint8)

    # another closing to fix weird checkerboards and stuff
    rivers = morph.binary_closing(rivers, morph.square(3))
    skeleton, distance = morph.medial_axis(rivers, return_distance=True)
    skeleton = skeleton.astype(np.uint8)
    skeleton[skeleton>0] = 1

    with rasterio.open(str(target[0]), 'w', **meta) as out:
        out.write(skeleton, 1)
    return 0


def keep_n_rivers(source, target, env):
    with rasterio.open(str(source[0])) as rast:
        rivers = rast.read(1)
        meta = rast.meta.copy()
    n = env.get('n', 1)

    labels = morph.label(rivers)
    rivsizes = [(i, (labels==i).sum()) for i in range(1, labels.max()+1)]
    maxlen = 0
    longest = 0
    for rivsize in rivsizes:
        if rivsize[1] > maxlen:
            longest = rivsize[0]
            maxlen = rivsize[1]
    rivers[labels != longest] = 0

    with rasterio.open(str(target[0]), 'w', **meta) as out:
        out.write(rivers, 1)
    return 0


def _count_and_trim_segment(j, i, skip, n, minlen, rivers, rivval):
    if n > minlen:
        return [(j,i)]
    downstream = []
    for dj in [-1, 0, 1]:
        for di in [-1, 0, 1]:
            if dj == di == 0:
                continue
            j2 = j + dj
            i2 = i + di
            if (j2<rivers.shape[0]) and (i2<rivers.shape[1]) and ((j2, i2) not in skip) and (rivers[j2, i2]==rivval):
                downstream.extend(_count_and_trim_segment(j2, i2, skip+[(j,i)], n+1, minlen, rivers, rivval))
    return (downstream + [(j,i)])

def _notnextto(a,b):
    if a[0] == b[0]:
        return abs(a[1]-b[1]) > 1
    if a[1] == b[1]:
        return abs(a[0]-b[0]) > 1
    return True

def trim_short_rivs(source, target, env):
    with rasterio.open(str(source[0])) as rast:
        rivers = rast.read(1)
        meta = rast.meta.copy()
    minlen = env.get('minlen', 10)

    rivval = rivers.max()
    wet = np.where(rivers==rivval)
    todelete = []
    for j,i in zip(*wet):
            branches = []
            for dj in [-1, 0, 1]:
                for di in [-1, 0, 1]:
                    if dj == di == 0:
                        continue
                    j2 = j + dj
                    i2 = i + di
                    if (j2<rivers.shape[0]) and (i2<rivers.shape[1]) and (rivers[j2, i2]==rivval):
                        branches.append((j2,i2))
            if len(branches) == 3 and np.all([_notnextto(a,b) for (a,b) in itertools.combinations(branches, 2)]): # self and 3 neighbors, not right next to each other. if they share faces then other configs that aren't splits can still have 3 neighbors
                for (j2, i2) in branches:
                    segment = _count_and_trim_segment(j2, i2, [(j,i)]+branches, 1, minlen, rivers, rivval)
                    if len(segment) < minlen:
                        todelete.extend(segment)
    for j,i in todelete:
        rivers[j,i] = 0

    with rasterio.open(str(target[0]), 'w', **meta) as out:
        out.write(rivers, 1)
    return 0


def import_rgis_network(source, target, env):
    with rasterio.open(str(source[0]), 'r') as rast:
        cellids = rast.read(1)
        nodata = rast.nodata
        affine = rast.affine
        assert rast.crs['init'] == 'epsg:4326'
    with rasterio.open(str(source[1]), 'r') as rast:
        basins = rast.read(1)
    with rasterio.open(str(source[2]), 'r') as rast:
        flowdir = rast.read(1)
    with rasterio.open(str(source[3]), 'r') as rast:
        rivers = rast.read(1)
        affine_riv = rast.affine
    with open(str(source[4]), 'r') as fin:
        proj4str = fin.read().strip()
    proj = pyproj.Proj(proj4str)

    neighbors = {
            1: (1, 0),
	    2: (1, 1),
	    4: (0, 1),
	    8: (-1, 1),
	    16: (-1, 0),
	    32: (-1, -1),
	    64: (0, -1),
	    128: (1, -1),
	    }

    minx, maxy = affine_riv * (0,0)
    maxx, miny = affine_riv * (rivers.shape[1], rivers.shape[0])

        #if not ((minx < px < maxx) and (miny < py < maxy)):

    G = nx.DiGraph()
    Gclip = nx.DiGraph()
    for j in range(cellids.shape[0]):
        for i in range(cellids.shape[1]):
            if cellids[j,i] != nodata:
                ll = affine * (i, j)
                xy = proj(*ll)
                G.add_node((i,j), **dict(ll=ll, xy=xy, basin=basins[j,i], cellid=cellids[j,i]))
                if (minx < xy[0] < maxx) and (miny < xy[1] < maxy):
                    Gclip.add_node((i,j), **dict(ll=ll, xy=xy, basin=basins[j,i], cellid=cellids[j,i]))
                tocell = flowdir[j,i]
                if tocell != 0:
                    di, dj = neighbors[tocell]
                    i2 = i + di
                    j2 = j + dj
                    if cellids[j2, i2] != nodata:
                        ll2 = affine * (i2, j2)
                        xy2 = proj(*ll2)
                        G.add_node((i2,j2), **dict(ll=ll2, xy=xy2, basin=basins[j2,i2], cellid=cellids[j2,i2]))
                        G.add_edge((i, j), (i2, j2))
                        if (minx < xy2[0] < maxx) and (miny < xy2[1] < maxy):
                            Gclip.add_node((i2,j2), **dict(ll=ll2, xy=xy2, basin=basins[j2,i2], cellid=cellids[j2,i2]))
                            if (i,j) in Gclip:
                                Gclip.add_edge((i, j), (i2, j2))
    for node in G.nodes():
        G.node[node]['upstream'] = len(nx.ancestors(G, node))
        G.node[node]['downstream'] = len(nx.descendants(G, node))
    nx.write_yaml(G, str(target[0]))

    for node in Gclip.nodes():
        Gclip.node[node]['upstream'] = len(nx.ancestors(G, node)) # use numbers from full network
        Gclip.node[node]['downstream'] = len(nx.descendants(G, node)) # use numbers from full network
    nx.write_yaml(Gclip, str(target[1]))
    return 0


def find_bifurs(source, target, env):
    with rasterio.open(str(source[0]), 'r') as rast:
        rivers = rast.read(1)
        meta = rast.meta.copy()

    # count all neighbor river cells for each river cell
    neighbor_elem = np.array([[1,1,1],
                              [1,0,1],
                              [1,1,1]])
    neighbors = generic_filter(rivers, np.sum, footprint=neighbor_elem)
    neighbors[rivers==0] = 0

    # bifurcations have 3 river neighbors. some curves can also have 3, but in these cases two of
    # the neighbors are adjacent to each other. so find those adjecent cells and remove them
    pair_elems = []
    pair_elems.append(np.array([[1,1,0],[0,0,0],[0,0,0]]))
    pair_elems.append(np.array([[0,1,1],[0,0,0],[0,0,0]]))
    pair_elems.append(np.array([[0,0,1],[0,0,1],[0,0,0]]))
    pair_elems.append(np.array([[0,0,0],[0,0,1],[0,0,1]]))
    pair_elems.append(np.array([[0,0,0],[0,0,0],[0,1,1]]))
    pair_elems.append(np.array([[0,0,0],[0,0,0],[1,1,0]]))
    pair_elems.append(np.array([[0,0,0],[1,0,0],[1,0,0]]))
    pair_elems.append(np.array([[1,0,0],[1,0,0],[0,0,0]]))
    pair_arrs = []
    for pair_elem in pair_elems:
        p = (generic_filter(rivers, np.sum, footprint=pair_elem)==2).astype(np.int)
        p[rivers==0] = 0
        pair_arrs.append(p)

    pairs = np.array(pair_arrs).sum(axis=0)

    bifurs = (neighbors - pairs).astype(np.uint8)
    # bifurs values:
    # 0: dry land
    # 1: river head or mouth
    # 2: river
    # 3: river bifurcation

    with rasterio.open(str(target[0]), 'w', **meta) as out:
        out.write(bifurs, 1)
    return 0


def plot_network_map(source, target, env):
    G = nx.read_yaml(str(source[0]))
    with rasterio.open(str(source[1]), 'r') as rast:
        bifurs = rast.read(1)
        affine = rast.affine

    mpl.style.use('ggplot')
    fig, ax = plt.subplots(1,1, figsize=(8, 12), dpi=300)

    pos = {node: G.node[node]['xy'] for node in G.node}
    basin = np.array([G.node[node]['basin'] for node in G.node])
    upstream = np.array([np.log(G.node[node]['upstream']+1) for node in G.node])

    labels = {node: G.node[node]['cellid'] for node in G.node}
    #labels = {node: node for node in G.node}
    nx.draw_networkx(G, pos, node_size=(upstream*10), node_color=basin,
            with_labels=True, labels=labels, font_size=6,
            arrowsize=8, edge_color='.5',
            cmap=palettable.cartocolors.qualitative.Bold_10.mpl_colormap, ax=ax)
    for t in ax.texts:
        t.set_clip_on(False)

    I, J = np.meshgrid(np.arange(bifurs.shape[1]), np.arange(bifurs.shape[0]))
    xs, ys = affine * (I.flatten(), J.flatten())
    X = xs.reshape(I.shape)
    Y = ys.reshape(J.shape)

    bifurs_mask = np.ma.masked_equal(bifurs, 0)
    ax.pcolormesh(X, Y, bifurs_mask, cmap=mpl.cm.Reds)
    ax.axis([X.min(), X.max(), Y.min(), Y.max()])
    ax.set_aspect('equal')

    fig.savefig(str(target[0]))
    return 0


def simple_bifurcations_OLD(source, target, env):
    G = nx.read_yaml(str(source[0]))
    with rasterio.open(str(source[1]), 'r') as rast:
        rivers = rast.read(1)
        affine = rast.affine
    with rasterio.open(str(source[2]), 'r') as rast:
        basins = rast.read(1)

    nodes = [node for node in G.nodes()]
    positions = [G.node[node]['xy'] for node in nodes]

    minx, maxy = affine * (0,0)
    maxx, miny = affine * (rivers.shape[1], rivers.shape[0])

    # remove network nodes outside delta domain
    #for n, (px, py) in list(zip(nodes, positions)):
        #if not ((minx < px < maxx) and (miny < py < maxy)):
            #nodes.remove(n)
            #positions.remove((px,py))

    # loop over river pixels, find nearest node
    wet = np.where(rivers>0)
    nearestnode = {}
    for j,i in zip(*wet):
        x, y = affine * (i,j)
        mindist = np.inf
        for idx, (px, py) in enumerate(positions):
            dist = np.sqrt((px-x)**2 + (py-y)**2)
            if dist < mindist:
                mindist = dist
                minidx = idx
        nearestnode[(j,i)] = nodes[minidx] #(px,py)

    nearestriv = {}
    for idx, (px, py) in enumerate(positions):
        mindist = np.inf
        for j,i in zip(*wet):
            x, y = affine * (i,j)
            dist = np.sqrt((px-x)**2 + (py-y)**2)
            if dist < mindist:
                mindist = dist
                minji = (j,i)
        nearestriv[nodes[idx]] = minji

    nodesonriv = set(nearestnodes.values())

# TODO - walk river from top to bottom, find where river meets new basin. bifurcate there!
# also TODO - check in WBM/MFlib for cycles? pretty sure a cycle would mess the sorting up since the recursive travel calc would never complete. maybe check if 100% of the flow is diverted, and if so delete the original downlink? Just avoid cycles here for now.

    # find node on river with largest "upstream"
    maxnupstream = -np.inf
    #for node in nodes:
    for node in nearestnode.values():
        #if node in nearestriv:
        nupcells = G.node[node]['upstream']
        if nupcells > maxnupstream:
            maxnupstream = nupcells
            maxnupnode = node
    curnode = maxnupnode
    # curnode should now be downstream-most node on OSM river

    #import ipdb;ipdb.set_trace() # TODO code here can diverge from OSM river
    # TODO: also an issue where nearest node isn't actually on the mainstem, and a close but larger node should be considered as "nearest". need some scoring function to asses how far away is still ok.
    while True:
        maxnupstream = -np.inf
        for upcell in G.predecessors(curnode):
            upcellupstream = G.nodes[upcell]['upstream']
            if (ipcellupstream in upcellupstream) > maxnupstream:
                maxnupstream = upcellupstream
                nextnode = upcell
        if maxnupstream == -np.inf:
            break
        curnode = nextnode
    # curnode should now be upstream-most node on OSM river.
    topnode = curnode

    (j,i) = nearestriv[topnode] # ASSUME THIS IS ON MAIN BASIN. other way?
    mainbasin = basins[j,i] # ASSUMING THIS IS MAIN BASIN

    offbasin = []
    for (j,i) in zip(*wet):
        if basins[j,i] != mainbasin:
            offbasin.append((j,i))

    # nodes for those river pixels not on main basin
    offnodes = [(nearestnode[(j,i)], basins[j,i]) for (j,i) in offbasin]

    otherbasins = {n[1] for n in offnodes}
    for otherbasin in otherbasins:
        othernodes = [offnode[0] for offnode in offbasin if offnode[1] == otherbasin]
        maxndown = -np.inf
        for othernode in othernodes:
            ndown = G.nodes[othernode]['downstream']
            if ndown > maxndown:
                maxndown = ndown
                topbasinnode = othernode
        # find nearest source node. start at top
        curnode = topnode
        mindist = np.inf
        while True:
            curnode = G.descendants(curnode)[0]
            topbasinxy = G.nodes[topbasinnode]['xy']
            curnodexy = G.nodes[curnode]['xy']
            dist = np.sqrt((topbasinxy[0] - curnodexy[0])**2 + (topbasinxy[1] - curnodexy[1])**2)
            if dist < mindist:
                mindist = dist
                nearestmainstemnode = curnode
        print("Fill basin", otherbasin, ": divert from", nearestmainstemnode, "to", topbasinnode)
        with open(str(target[0]), 'a') as fout:
            fout.write('{0}, {1}, 0.5\n'.format(G.nodes[nearestmainstemnode]['cellid']+1, G.nodes[topbasinnode]['cellid']+1))

    return 0


def _find_largest_node_i(rivxy, positions, nupstream, nodemask, resolution, searchradius):
    rivx, rivy = rivxy
    dists = [np.sqrt((rivx - nodex)**2 + (rivy - nodey)**2) for (nodex, nodey) in positions]
    score = []
    for d, n, oknode in zip(dists, nupstream, nodemask):
        if oknode and d<(resolution * searchradius):
            score.append(n)
        else:
            score.append(-1)
    if np.max(score) == -1:
        return None
    return np.argmax(score)

def _find_nearest_node_i(rivxy, positions, nodemask):
    rivx, rivy = rivxy
    dists = [np.sqrt((rivx - nodex)**2 + (rivy - nodey)**2) for (nodex, nodey) in positions]
    score = []
    for i, (d, oknode) in enumerate(zip(dists, nodemask)):
        if oknode:
            score.append(d)
        else:
            score.append(np.inf)
    return np.argmin(score)

def simple_bifurcations(source, target, env):
    G = nx.read_yaml(str(source[0]))
    with rasterio.open(str(source[1]), 'r') as rast:
        rivers = rast.read(1)
        affine = rast.affine
    with rasterio.open(str(source[2]), 'r') as rast:
        basins = rast.read(1)

    nodes = [node for node in G.nodes()]
    positions = [G.node[node]['xy'] for node in nodes]
    nupstream = [G.node[node]['upstream'] for node in nodes]
    ndownstream = [G.node[node]['downstream'] for node in nodes]
    nodebasins = [G.node[node]['basin'] for node in G.node]
    cellid = [G.node[node]['cellid'] for node in G.node]

    p0 = positions[0]
    resolution = np.min([np.sqrt((p0[0]-p[0])**2 + (p0[1]-p[1])**2) for p in positions[1:]])

    minx, maxy = affine * (0,0)
    maxx, miny = affine * (rivers.shape[1], rivers.shape[0])

    mouthnodei = np.argmax(nupstream)
    mouthnode = nodes[mouthnodei]
    mainbasin = basins[mouthnode[1], mouthnode[0]]
    mainbasinmask = [b == mainbasin for b in nodebasins]
    allbasinmask = [1 for b in nodebasins]

    #nearestmainnode = {}
    #wet = np.where(rivers>0)
    #for j, i in zip(*wet):
        #xy = affine * (i,j)
        #nearest_node_i = _find_nearest_node_i(xy, positions, nupstream, mainbasinmask, resolution, 2)
        #nearestmainnode[(j,i)] = nodes[nearest_node_i]

    # find upstream most river point
    endpoints = np.where(rivers==1)
    maxdownstream = -np.inf
    headnode_i = None
    head_endpoint_i = None
    for ii, (j, i) in enumerate(zip(*endpoints)):
        xy = affine * (i,j)
        largest_near_node_i = _find_largest_node_i(xy, positions, nupstream, mainbasinmask, resolution, 2)
        if largest_near_node_i is not None:
            if ndownstream[largest_near_node_i] > maxdownstream:
                maxdownstream = ndownstream[largest_near_node_i]
                headnode_i = largest_near_node_i
                head_endpoint_i = ii
    headnode = nodes[headnode_i]

    # walk down nodes to identify mainstem nodes
    mainstem_nodes = [headnode]
    to_visit = list(G.successors(headnode))
    while to_visit:
        node = to_visit.pop()
        mainstem_nodes.append(node)
        downstream = list(G.successors(node))
        to_visit.extend(downstream)
    mainstem_mask = [n in mainstem_nodes for n in nodes]

    # walk down river
    bifurcations = {}
    visited = set()
    to_visit = [((endpoints[0][head_endpoint_i], endpoints[1][head_endpoint_i]), headnode_i)] # riv point, last node index
    while to_visit:
        (rivj, rivi), last_node_i = to_visit.pop(0)
        visited.add((rivj, rivi))

        xy = affine * (rivi,rivj)
        this_node_i = _find_largest_node_i(xy, positions, nupstream, allbasinmask, resolution, 2)
        this_node_basin = nodebasins[this_node_i]
        if (this_node_basin != nodebasins[last_node_i]) and (nodebasins[last_node_i] == mainbasin): # only allows bifurs from mainstem
            mainstem_source_i = _find_nearest_node_i(xy, positions, mainstem_mask)
            from_cell = cellid[mainstem_source_i]
            to_cell = cellid[this_node_i]
            if (from_cell, to_cell) not in bifurcations:
                bifurcations[(from_cell, to_cell)] = .5 #assumes from_cell is currently at 1. need to track?
                print('Bifurcation: cellid {0} to {1}'.format(from_cell, to_cell))

        for dj in [-1,0,1]:
            for di in [-1,0,1]:
                rivj2 = rivj + dj
                rivi2 = rivi + di
                if rivers[rivj2,rivi2]>0 and (rivj2, rivi2) not in visited:
                    to_visit.append(((rivj2, rivi2), this_node_i))

    with open(str(target[0]), 'w', newline='') as fout:
        csvwriter = csv.writer(fout)
        for (from_cell, to_cell), frac in bifurcations.items():
            csvwriter.writerow([from_cell, to_cell, frac])
    return 0

