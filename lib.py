import csv
import pickle
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
from collections import defaultdict

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
    deltahull = geopandas.GeoDataFrame(deltahull.buffer(25000), columns=['geometry'])

    rivers_clip = geopandas.overlay(rivers, deltahull, how='intersection') #slow

    rivers_clip.to_file(str(target[0]))
    with open(str(target[1]), 'w') as fout:
        fout.write(laea.proj4_init + '\n')
    return 0


def project_and_clip_osm_waterways(source, target, env):
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
    deltahull = geopandas.GeoDataFrame(deltahull.buffer(25000), columns=['geometry'])

    rivgeom = rivers['geometry']
    inside = rivgeom.intersects(deltahull['geometry'].item())
    rivers_clip = rivers[inside]

    rivers_clip.to_file(str(target[0]))
    with open(str(target[1]), 'w') as fout:
        fout.write(laea.proj4_init + '\n')
    return 0


def clip_osm_rivers(source, target, env):
    rivers_ll = geopandas.read_file(str(source[0]))
    delta_ll = geopandas.read_file(str(source[1]))

    deltahull_ll = geopandas.GeoDataFrame(delta_ll.dissolve(by='Delta').convex_hull, columns=['geometry'])
    deltahull_ll.crs = delta_ll.crs

    rivers_clip = geopandas.overlay(rivers_ll, deltahull_ll, how='intersection')

    rivers_clip.to_file(str(target[0]))
    return 0


def filter_waterway_types(source, target, env):
    rivers = geopandas.read_file(str(source[0]))
    rivers = rivers[rivers['fclass'] == 'river']

    rivers.to_file(str(target[0]))
    return 0


def thin_vec(source, target, env):
    rivers = geopandas.read_file(str(source[0]))

    if not env.get('wetlands', True):
        # keep unspecified, reservoirs, and rivers and river classes. lakes are often unspec.
        rivers = rivers[(rivers['code']==8200) | (rivers['code']==8201) | (rivers['code']==8202)]
    # drop small area stuff
    minarea = env.get('minarea', 0)
    rivers = rivers[rivers.area > minarea]

    # drop thin stuff
    thinning = env.get('thinning', 100)
    rivers_eroded = rivers.buffer(-thinning)
    rivers_clean = rivers_eroded[rivers_eroded.area>0].buffer(thinning)

    rivers_merge = geopandas.overlay(geopandas.GeoDataFrame(rivers_clean, columns=['geometry']), geopandas.GeoDataFrame(rivers_clean, columns=['geometry']), how='union')
    rivers_merge.crs = rivers.crs

    rivers_merge.to_file(str(target[0]))
    return 0


def plot_vec_rivs(source, target, env):
    mpl.style.use('ggplot')

    rivers = geopandas.read_file(str(source[0]))

    fig, ax = plt.subplots(1, 1, figsize=(10,10))
    ax.set_facecolor('w')
    rivers.plot(color='k', ax=ax)
    ax.set_aspect('equal')
    ax.xaxis.set_ticks([])
    ax.yaxis.set_ticks([])

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

    rivers = morph.skeletonize(rivers)
    rivers = rivers.astype(np.uint8)

    # another closing to fix weird checkerboards and stuff
    rivers = morph.binary_closing(rivers, morph.square(3))
    skeleton = morph.skeletonize(rivers)
    skeleton = skeleton.astype(np.uint8)
    skeleton[skeleton>0] = 1


    # at this point can still find non-single-pixel rivers. If multiple rivers intersect
    # can have "clumps" that aren't identified as bifurs
    # look for 2x2 river clumps, expand to 3x3 region, remove pixels in clump one at a time
    # and find a removal that doesn't cut off a branch
    clump = np.ones((2,2))
    p = generic_filter(skeleton, np.sum, footprint=clump)
    clumps = np.where(p == 4)
    fixed = []
    for j,i in zip(*clumps):
        fixed.append(None)
        region4x4 = skeleton[j-2:j+2, i-2:i+2] # 4x4 around 2x2 clump
        for dj, di in [(1,1), (1,2), (2,1), (2,2)]:
            region = region4x4.copy()
            region[dj,di] = 0
            labels = morph.label(region, neighbors=8)
            if labels.max() == 1:
                fixed[-1] = (dj, di)
                break
        if fixed[-1] is None:
            raise NotImplementedError('changing single pixel is not enough to break 2x2 clump')
        else:
            dj, di = fixed[-1]
            skeleton[j + (dj - 2), i + (di - 2)] = 0

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


def _walk_to_branch_end(j, i, skip, n, minlen, rivers):
    if n > minlen:
        return [(j,i)]
    if rivers[j,i] >= 3:
        # found another bifurcation. return a long list of fake points so this segment
        #isn't deleted. last point can't have rivers == 1, so use this one, which is == 3
        return [(j,i)] * (minlen+1)
    downstream = []
    for (dj, di) in [(-1,0),(0,1),(1,0),(0,-1),(-1,-1),(1,-1),(1,1),(-1,1)]: # list with horizontal and vert before diag, in case of ambig route: example: bifur to right and a branch from that to down. down could be jumped to diagonally, so check hor/vert first and break if bifur is found
        if dj == di == 0:
            continue
        j2 = j + dj
        i2 = i + di
        if (j2<rivers.shape[0]) and (i2<rivers.shape[1]) and ((j2, i2) not in skip) and (rivers[j2, i2]>0):
            downstream.extend(_walk_to_branch_end(j2, i2, skip+[(j,i)], n+1, minlen, rivers))
            if len(downstream) > minlen: # need this if bifur point was found, dont want to keep looping around looking for more branches
                break
    return [(j,i)] + downstream

def trim_short_rivs(source, target, env):
    with rasterio.open(str(source[0])) as rast:
        rivers = rast.read(1)
        meta = rast.meta.copy()
    minlen = env.get('minlen', 10)

    # need to loop since after a pass we could still be left with short branches.
    # if we have a short segment between two bifur points, with two short branches on one end,
    # cleaning one short branch can still leave the bifur->bifur segment + one branch
    # too short. trim again and loop as long as rivers keeps changing
    while True:
        prevrivers = rivers.copy()

        bifurs = np.where(rivers>=3)
        for j,i in zip(*bifurs):
            branches = []
            todelete = []
            for dj in [-1, 0, 1]:
                for di in [-1, 0, 1]:
                    if dj == di == 0:
                        continue
                    j2 = j + dj
                    i2 = i + di
                    if (j2<rivers.shape[0]) and (i2<rivers.shape[1]) and (rivers[j2, i2]>0):
                        branches.append((j2,i2))
            for (j2, i2) in branches:
                segment = _walk_to_branch_end(j2, i2, [(j,i)]+branches, 1, minlen, rivers)
                if rivers[segment[-1]] == 1: # found terminating segment shorter than minlen
                    todelete.append(segment)
            #if len(segments) == 1, just a stub, gets deleted (other branches are longer, or non-terminating)
            if todelete:
                if len(todelete) == 2:
                    # two terminating branches
                    # keep longer segment, will get appended to last segment
                    seglens = [len(segment) for segment in todelete]
                    longest_i = np.argmax(seglens)
                    todelete.pop(longest_i)
                rivers[j,i] -= 1 # old bifur point becomes normal river, or a four-way becomes three-way
                for segment in todelete:
                    for rivpt in segment: # bifur point not included on segment
                        rivers[rivpt] = 0

        # run skeleton to clean corners where short segments were clipped off
        # dont want to overwrite bifur info though
        justrivers = (rivers>0).astype(np.uint8)
        justrivers = morph.skeletonize(justrivers)
        rivers[justrivers==0] = 0

        if np.all(rivers == prevrivers):
            break

    with rasterio.open(str(target[0]), 'w', **meta) as out:
        out.write(justrivers.astype(np.uint8), 1)
    return 0


def import_rgis_network(source, target, env):
    with rasterio.open(str(source[0]), 'r') as rast:
        cellids = rast.read(1)
        nodata = rast.nodata
        affine = rast.transform
        assert rast.crs['init'] == 'epsg:4326'
    with rasterio.open(str(source[1]), 'r') as rast:
        basins = rast.read(1)
    with rasterio.open(str(source[2]), 'r') as rast:
        flowdir = rast.read(1)
    with rasterio.open(str(source[3]), 'r') as rast:
        rivers = rast.read(1)
        affine_riv = rast.transform
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

    G = nx.DiGraph()
    Gclip = nx.DiGraph()
    nj, ni = cellids.shape
    for j in range(nj):
        for i in range(ni):
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
                    if (0<=i2<ni and 0<=j2<nj) and cellids[j2, i2] != nodata:
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
        affine = rast.transform
    with rasterio.open(str(source[2]), 'r') as rast:
        extended_bifurs = rast.read(1)
    labeltype = env['labels']

    # reset upstream counts
    for node in G.nodes():
        G.node[node]['upstream'] = len(nx.ancestors(G, node))
        G.node[node]['downstream'] = len(nx.descendants(G, node))

    mpl.style.use('ggplot')
    fig, ax = plt.subplots(1,1, figsize=(8, 12))#, dpi=300)

    pos = {node: G.node[node]['xy'] for node in G.node}
    basin = np.array([G.node[node]['basin'] for node in G.node])
    upstream = np.array([np.log(G.node[node]['upstream']+1) for node in G.node])

    if labeltype == 'cells':
        with_labels = True
        labels = {node: G.node[node]['cellid'] for node in G.node}
    elif labeltype == 'nodes':
        with_labels = True
        labels = {node: node for node in G.node}
    elif labeltype == 'none':
        with_labels = False
        labels = {}
    nx.draw_networkx_nodes(G, pos, node_size=(upstream*20), node_color=basin,
            with_labels=with_labels, labels=labels, font_size=6,
            alpha=.5, cmap=palettable.cartocolors.qualitative.Bold_10.mpl_colormap, ax=ax)
    nx.draw_networkx_edges(G, pos, arrowstyle='->', arrowsize=10, node_size=20, edge_color='k', ax=ax)
    if with_labels:
        nx.draw_networkx_labels(G, pos, labels=labels, font_size=6, ax=ax)
        for t in ax.texts:
            t.set_clip_on(False)
            t.set_rotation(30)

    I, J = np.meshgrid(np.arange(bifurs.shape[1]), np.arange(bifurs.shape[0]))
    xs, ys = affine * (I.flatten(), J.flatten())
    X = xs.reshape(I.shape)
    Y = ys.reshape(J.shape)

    ext_bifurs_mask = np.ma.masked_equal(extended_bifurs, 0)
    ax.pcolormesh(X, Y, ext_bifurs_mask, cmap=mpl.cm.Blues)

    bifurs_mask = np.ma.masked_equal(bifurs, 0)
    ax.pcolormesh(X, Y, bifurs_mask, cmap=mpl.cm.Reds)
    ax.axis([X.min(), X.max(), Y.min(), Y.max()])
    ax.set_aspect('equal')


    ax.xaxis.set_ticks([])
    ax.yaxis.set_ticks([])

    fig.savefig(str(target[0]))
    if env['inspect'] is not None:
        plt.show()
    return 0


def plot_flowdirs_map(source, target, env):
    with rasterio.open(str(source[0]), 'r') as rast:
        bifurs = rast.read(1)
        affine = rast.transform

    with rasterio.open(str(source[1]), 'r') as rast:
        extended_bifurs = rast.read(1)

    with open(str(source[2]), 'rb') as fin:
        segments = pickle.load(fin)

    with open(str(source[3]), 'rb') as fin:
        next_rivpts = pickle.load(fin)

    mpl.style.use('ggplot')
    fig, ax = plt.subplots(1,1, figsize=(8, 12))#, dpi=300)

    for segi, segment in segments.items():
        x1, y1 = affine * segment[0][::-1]
        x2, y2 = affine * segment[-1][::-1]
        if segment[1] in next_rivpts[segment[0]]:
            # flows from 0 to end
            ax.annotate("", xy=(x2,y2), xytext=(x1,y1),
                    arrowprops=dict(facecolor='k', edgecolor='k', arrowstyle='-|>'))
            ax.text((x2+x1)/2, (y2+y1)/2, segi)
        elif segment[0] in next_rivpts[segment[1]]:
            # flows from end to 0
            ax.annotate("", xy=(x1,y1), xytext=(x2,y2),
                    arrowprops=dict(facecolor='k', edgecolor='k', arrowstyle='-|>'))
            ax.text((x2+x1)/2, (y2+y1)/2, segi)

    I, J = np.meshgrid(np.arange(bifurs.shape[1]), np.arange(bifurs.shape[0]))
    xs, ys = affine * (I.flatten(), J.flatten())
    X = xs.reshape(I.shape)
    Y = ys.reshape(J.shape)

    ext_bifurs_mask = np.ma.masked_equal(extended_bifurs, 0)
    ax.pcolormesh(X, Y, ext_bifurs_mask, cmap=mpl.cm.Blues)

    bifurs_mask = np.ma.masked_equal(bifurs, 0)
    ax.pcolormesh(X, Y, bifurs_mask, cmap=mpl.cm.Reds)
    ax.axis([X.min(), X.max(), Y.min(), Y.max()])
    ax.set_aspect('equal')

    ax.xaxis.set_ticks([])
    ax.yaxis.set_ticks([])

    fig.savefig(str(target[0]))
    if env['inspect'] is not None:
        plt.show()
    return 0


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


def remap_riv_network(source, target, env):
    import skimage.draw

    def _merge_riv_path_to_mainstem(rivers, head_rivpt, initial_riv_pts, nodes, nupstream, ndownstream, positions, nearestnode_to_riv, next_rivpt, prev_rivpt, affine):

        def _trim_river(rivers, rivpt, initial_riv_pts, next_rivpt, prev_rivpt, nearestnode_to_riv):
            # erase river above rivpt, and adjust next_rivpt and prev_rivpt dicts
            thispt = rivpt
            rivers[thispt] = 1
            initial_riv_pts.append(thispt)
            to_visit = [(thispt, uppt) for uppt in prev_rivpt[thispt]]
            while to_visit:
            #while thispt in prev_rivpt:
                downpt, thispt = to_visit.pop()
                #downpt = thispt
                #thispt = prev_rivpt[downpt]
                rivers[thispt] = 0
                try:
                    del prev_rivpt[downpt]
                except KeyError: # could have already been deleted if on branch upstream of a convergence
                    pass
                del next_rivpt[thispt]
                del nearestnode_to_riv[thispt]
                if thispt in initial_riv_pts:
                    initial_riv_pts.remove(thispt)
                if thispt in prev_rivpt:
                    to_visit.extend([(thispt, uppt) for uppt in prev_rivpt[thispt]])
            return rivers, initial_riv_pts, next_rivpt, prev_rivpt, nearestnode_to_riv

        def _extend_river(rivers, mindist_rivpt, head_node_ij, initial_riv_pts, next_rivpt, prev_rivpt, nearestnode_to_riv, positions):
            # draw new river in straight line from closest river approach to head_node
            head_node_ji = tuple(head_node_ij[::-1])
            rowidx, colidx = skimage.draw.line(*mindist_rivpt, *head_node_ji)
            if (rowidx[-1] == mindist_rivpt[0]) and (colidx[-1] == mindist_rivpt[0]):
                # reverse so we walk through from rivpt toward node (upstream)
                rowidx = rowidx[::-1]
                colidx = colidx[::-1]
            rivers[rowidx, colidx] = 2
            rivers[head_node_ji] = 1
            downj = None; downi = None
            for j, i in zip(rowidx, colidx):
                # walk new line, setup prev, next, and nearest dicts
                if (downj is not None) and (downi is not None):
                    if (j,i) not in prev_rivpt[downj,downi]:
                        prev_rivpt[downj,downi].append((j, i))
                    if (downj,downi) not in next_rivpt[j,i]:
                        next_rivpt[j,i].append((downj, downi))
                    xy = affine * (i, j)
                    allbasinmask = [1 for pos in positions]
                    nearest_node_i = _find_nearest_node_i(xy, positions, allbasinmask)
                    nearestnode_to_riv[j, i] = nodes[nearest_node_i]
                downj = j
                downi = i
            initial_riv_pts.remove(mindist_rivpt)
            initial_riv_pts.append(head_node_ji)
            return rivers, initial_riv_pts, next_rivpt, prev_rivpt, nearestnode_to_riv

        score = np.array(nupstream) * np.array(ndownstream) # maximized at "center" of mainstem, which we assume is upstream of the delta boundary
        head_node_i = np.argmax(score)
        head_node = nodes[head_node_i]
        nearest = nearestnode_to_riv[head_rivpt]

        if head_node == nearest:
            return rivers, head_rivpt, initial_riv_pts, head_node, next_rivpt, prev_rivpt, nearestnode_to_riv

        found_head_node_on_riv = False
        for rivpt, node in nearestnode_to_riv.items():
            if node == head_node:
                found_head_node_on_riv = True
                break
        if found_head_node_on_riv:
            if rivpt not in initial_riv_pts:
                rivers, initial_riv_pts, next_rivpt, prev_rivpt, nearestnode_to_riv = _trim_river(rivers, rivpt, initial_riv_pts, next_rivpt, prev_rivpt, nearestnode_to_riv)
            return rivers, rivpt, initial_riv_pts, node, next_rivpt, prev_rivpt, nearestnode_to_riv
        else:
            to_visit = initial_riv_pts.copy()
            head_node_xy = positions[head_node_i]
            head_node_x, head_node_y = positions[head_node_i]
            head_node_ij = [int(val) for val in ~affine * (head_node_x, head_node_y)]
            mindist = np.inf
            while to_visit:
                rivpt = to_visit.pop()
                rivj, rivi = rivpt
                rivx, rivy = affine * (rivi,rivj)
                dist = np.sqrt((head_node_x - rivx)**2 + (head_node_y - rivy)**2)
                if dist < mindist:
                    mindist = dist
                    mindist_rivpt = rivpt
                to_visit.extend(next_rivpt[rivpt])
            if mindist_rivpt not in initial_riv_pts:
                rivers, initial_riv_pts, next_rivpt, prev_rivpt, nearestnode_to_riv = _trim_river(rivers, mindist_rivpt, initial_riv_pts, next_rivpt, prev_rivpt, nearestnode_to_riv)
            rivers, initial_riv_pts, next_rivpt, prev_rivpt, nearestnode_to_riv = _extend_river(rivers, mindist_rivpt, head_node_ij, initial_riv_pts, next_rivpt, prev_rivpt, nearestnode_to_riv, positions)
            return rivers, head_rivpt, initial_riv_pts, head_node, next_rivpt, prev_rivpt, nearestnode_to_riv


    G = nx.read_yaml(str(source[0]))
    Gorig = G.copy()
    with rasterio.open(str(source[1]), 'r') as rast:
        rivers = rast.read(1)
        affine = rast.transform
        meta = rast.meta.copy()
    with rasterio.open(str(source[2]), 'r') as rast:
        basins = rast.read(1)

    flowdir_weights = env.get('flowdir_weights', (1,2))

    nodes = [node for node in G.nodes()]
    positions = [G.node[node]['xy'] for node in nodes]
    nupstream = [G.node[node]['upstream'] for node in nodes]
    ndownstream = [G.node[node]['downstream'] for node in nodes]
    nodebasins = [G.node[node]['basin'] for node in nodes]
    cellid = [G.node[node]['cellid'] for node in nodes]

    p0 = positions[0]
    resolution = np.min([np.sqrt((p0[0]-p[0])**2 + (p0[1]-p[1])**2) for p in positions[1:]])

    minx, maxy = affine * (0,0)
    maxx, miny = affine * (rivers.shape[1], rivers.shape[0])

    mouthnodei = np.argmax(nupstream)
    mouthnode = nodes[mouthnodei]
    mainbasin = basins[mouthnode[1], mouthnode[0]]
    mainbasinmask = [b == mainbasin for b in nodebasins]
    allbasinmask = [1 for b in nodebasins]

    nearestnode_to_riv = {}
    nearestnode_ndownstream = {} # dict on riv points, lower ndownstream means closer to coast. use to assess flow dir on river
    wet = np.where(rivers>0)
    for j, i in zip(*wet):
        xy = affine * (i,j)
        nearest_node_i = _find_nearest_node_i(xy, positions, allbasinmask)
        nearestnode_to_riv[(j,i)] = nodes[nearest_node_i]
        nearestnode_ndownstream[(j,i)] = ndownstream[nearest_node_i]

    endpoints = np.where(rivers==1)
    endpoints_ndown = [nearestnode_ndownstream[(j,i)] for (j,i) in zip(*endpoints)]
    ndown_mean = np.mean(endpoints_ndown)
    ndown_std = np.std(endpoints_ndown)
    ndown_z = (np.array(endpoints_ndown) - ndown_mean) / ndown_std
    upstream_endpoints_ind = np.where(ndown_z > 0)
    upstream_endpoints = [list(zip(*endpoints))[i] for i in upstream_endpoints_ind[0]]
    head_endpoint_i = np.argmax(ndown_z[upstream_endpoints_ind])
    head_rivpt = (endpoints[0][head_endpoint_i], endpoints[1][head_endpoint_i])

    # Find "coastal" nodes by seeing where a node has a missing neighbor
    # (ignore border nodes since those probably on clipped delta boundary - ok since dist will find a neighbor)
    coastal = [] # true or false
    # Use that to calc dist-to-coast for each node
    dist_to_coast = []
    # That determines flowdir

    # find coastal nodes
    minj = min([node[0] for node in nodes])
    maxj = max([node[0] for node in nodes])
    mini = min([node[1] for node in nodes])
    maxi = max([node[1] for node in nodes])
    for node in nodes:
        if (not (minj < node[0] < maxj) or not (mini < node[1] < maxi)):
            coastal.append(False)
            continue
        foundedge = False
        for dj in [-1,0,1]:
            for di in [-1,0,1]:
                if (dj == di == 0):
                    continue
                othernode = (node[0] + dj, node[1] + di)
                if othernode not in nodes:
                    foundedge = True
                    break
            if foundedge:
                break
        coastal.append(foundedge)

    # calc dist-to-coast
    for node in nodes:
        mindist = np.inf
        for i, iscoastal in enumerate(coastal):
            if not iscoastal:
                continue
            coastalnode = nodes[i]
            dist = np.sqrt((node[0]-coastalnode[0])**2 + (node[1]-coastalnode[1])**2)
            if dist < mindist:
                mindist = dist
        dist_to_coast.append(mindist)


    # Determine flow direction of each river segment. Given branching and convergence, some segs might be ambiguoous
    # starting from highest endpoint, walk downstream on riv following branches
    # find segments as all riv points between branches
    # calc mean(diff([nearestnode_ndownstream(rivpt) for each riv pt in segment in order from bifur point])
    # if value is positive, it means the branch is mostly moving upstream. set all riv points from bifur point to end to zero (remove river)
    (j,i) = head_rivpt
    cursegi = 0
    tovisit = [(j,i,cursegi)]
    visited = set()
    segments = defaultdict(list)
    next_rivpt = defaultdict(list) # dictionary keyed on rivpt (j,i) indicating next pt flow-dir-wise (will be multiple pts at bifur
    prev_rivpt = defaultdict(list)
    maxsegi = cursegi
    while tovisit:
        j, i, cursegi = tovisit.pop(0)
        # check if we're starting along an already traveled branch (could happend with bifur -> convergence)
        if ((rivers[j,i] == 2) and # on regular river point
                (len(segments[cursegi])==1) and # first point is bifur, so check if this is second
                np.any([(j,i) in seg for seg in segments.values()])): # if point is on another segment
            # ignore this branch, already traveled
            earlier = np.where([(j,i) in seg for seg in segments.values()])[0].squeeze()
            print('Deleting segment {0}, already seen in segment(s) {1}'.format(cursegi, earlier))
            del segments[cursegi]
            continue
        segments[cursegi].append((j,i))
        if rivers[j,i] in [1,3]: # on start/endpoint or bifur point
            segpts = segments[cursegi]
            if len(segpts) > 1: # at end of segment
                # use both ndownstream metric (which uses initial non-bifur network), and dist-to-coast by summimg them. double weight on dist_to_coast. both individually have problems
                nodescores = [flowdir_weights[0]*nearestnode_ndownstream[pt] + flowdir_weights[1]*dist_to_coast[nodes.index(nearestnode_to_riv[pt])] for pt in segpts]
                n_nodes_on_seg = len({nearestnode_to_riv[pt] for pt in segpts})

                diffscores = np.diff(nodescores).tolist()
                while 0 in diffscores:
                    diffscores.remove(0)
                diffscores = np.array(diffscores)
                diffscores[diffscores<0] = -1
                diffscores[diffscores>0] = 1
                if len(diffscores) == 0:
                    diffscores = [0]
                dir_metric = np.mean(diffscores) # collapse to 0,1 to remove influence of a few outlier nodes
                # mark next riv pts
                if ((dir_metric <= 0) or (n_nodes_on_seg <= 1) or ((dir_metric <= .03) and (segpts[0] in upstream_endpoints))) and not ((dir_metric >= -.03) and segpts[-1] in upstream_endpoints): # downstream, and dont set very short segs to upstream, and downstream if segment has upstream endpoint (but if dir_metric is very positive, then ignore upstream_endpoint
                    print(cursegi, segpts[0], segpts[-1], dir_metric, 'downstream')
                    for thisji, nextji in zip(segpts[:-1], segpts[1:]):
                        if nextji not in next_rivpt[thisji]:
                            next_rivpt[thisji].append(nextji)
                        if thisji not in prev_rivpt[nextji]:
                            prev_rivpt[nextji].append(thisji)
                else: # upstream
                    print(cursegi, segpts[0], segpts[-1], dir_metric, 'upstream')
                    for thisji, nextji in zip(segpts[1:], segpts[:-1]):
                        if nextji not in next_rivpt[thisji]:
                            next_rivpt[thisji].append(nextji)
                        if thisji not in prev_rivpt[nextji]:
                            prev_rivpt[nextji].append(thisji)

        visited.add((j,i))
        for (dj, di) in [(-1,0),(0,1),(1,0),(0,-1),(-1,-1),(1,-1),(1,1),(-1,1)]: # list with horizontal and vert before diag, in case of ambig route: example: bifur to right and a branch from that to down. down could be jumped to diagonally, so check hor/vert first and break if bifur is found
                if dj == di == 0:
                    continue
                j2 = j + dj
                i2 = i + di
                if (rivers[j2,i2]>0) and ((j2, i2) not in segments[cursegi]) and (((j2,i2) not in visited) or rivers[j2,i2]==3): # dont go to a point already on segment, or one already visited on another segment (but you can repeat if its a bifur point - convergences can be repeated)
                    if rivers[j,i] == 3: # bifur point
                        maxsegi += 1 # each branch gets new cursegi
                        nextsegi = maxsegi
                        segments[nextsegi].append((j,i)) # add current bifur point as start of new segment
                        tovisit.append((j2, i2, nextsegi)) # add new segments to end of queue
                    else:
                        tovisit.insert(0, (j2, i2, cursegi)) # on current segment, add to front of queue to stay on this branch
                    if rivers[j,i] == 2:
                        # regular path, not bifur, only one neighbor. break so as not to find incorrect diag branch
                        break
# by here, next_rivpt gives us flowdirection

    ## Go back over segments: for short segments, use direction of longest neighbor. effectively attaches it
    # do all first, then another pass to check/change small ones. look at segments that share same endpoints
    # upstream endpoints marked good
    #import ipdb;ipdb.set_trace()
    seglens = {}
    for segi, segment in segments.items():
        seglens[segi] = len({nearestnode_to_riv[pt] for pt in segment})
    goodsegs = {segi: ((seglens[segi] > 2) or (seg[0] in upstream_endpoints) or (seg[-1] in upstream_endpoints)) for segi, seg in segments.items()}
    # set direction of small segments to match longest neighboring segment that flows INTO segment
    # track "good" (long or fixed) segments, iterate until all good
    prevsegs = None
    while not np.all(list(goodsegs.values())):
        if goodsegs == prevsegs:
            with open(str(target[3]), 'wb') as fout:
                pickle.dump(segments, fout)
            with open(str(target[4]), 'wb') as fout:
                pickle.dump(next_rivpt, fout)
            raise NotImplementedError('Flow direction can not be determined for some segments. Goodsegs: {}'.format(str(goodsegs)))
        prevsegs = goodsegs

        for segi, segment in segments.items():
            if not goodsegs[segi]:
                neighbor_segments = []
                for segj, otherseg in segments.items():
                    if otherseg == segment:
                        continue
                    if len(otherseg) < 3:
                        continue
                    if (((otherseg[0] == segment[0]) and (otherseg[1] in next_rivpt[otherseg[2]])) or
                        ((otherseg[0] == segment[-1]) and (otherseg[1] in next_rivpt[otherseg[2]])) or
                        ((otherseg[-1] == segment[0]) and (otherseg[-2] in next_rivpt[otherseg[-3]])) or
                        ((otherseg[-1] == segment[-1]) and (otherseg[-2] in next_rivpt[otherseg[-3]]))):
                        # otherseg flows INTO segment
                        if goodsegs[segj]:
                            neighbor_segments.append((segj, otherseg))
                if not neighbor_segments:
                    continue
                longest_i = np.argmax([len(seg) for segj, seg in neighbor_segments])
                segj, otherseg = neighbor_segments[longest_i]
                # check how segments are aligned, and set segment based on flow dir in otherseg
                if otherseg[0] == segment[0]:
                    if otherseg[2] in next_rivpt[otherseg[1]]: # go inside otherseg a bit, endpoint may have been overwritten by this short segment in code above
                        pairs = zip(segment[1:], segment[:-1])
                    else:
                        pairs = zip(segment[:-1], segment[1:])
                if otherseg[0] == segment[-1]:
                    if otherseg[2] in next_rivpt[otherseg[1]]:
                        pairs = zip(segment[:-1], segment[1:])
                    else:
                        pairs = zip(segment[1:], segment[:-1])
                if otherseg[-1] == segment[0]:
                    if otherseg[-2] in next_rivpt[otherseg[-3]]:
                        pairs = zip(segment[:-1], segment[1:])
                    else:
                        pairs = zip(segment[1:], segment[:-1])
                if otherseg[-1] == segment[-1]:
                    if otherseg[-2] in next_rivpt[otherseg[-3]]:
                        pairs = zip(segment[1:], segment[:-1])
                    else:
                        pairs = zip(segment[:-1], segment[1:])
                for i, (thisji, nextji) in enumerate(pairs):
                    # dont just overwrite. bifur points need multiple next_rivpts, dont want to lose other branch
                    if i==0:
                        if nextji not in next_rivpt[thisji]:
                            next_rivpt[thisji].append(nextji)
                    else:
                        next_rivpt[thisji] = [nextji]
                    if i==(len(list(pairs)) - 1):
                        if thisji not in prev_rivpt[nextji]:
                            prev_rivpt[nextji].append(thisji)
                    else:
                        prev_rivpt[nextji] = [thisji]
                goodsegs[segi] = True
                print('Fixed direction on segment {0}, matched to segment {1}'.format(segi, segj))



    initial_riv_pts = [(endpoints[0][endpoint_i], endpoints[1][endpoint_i]) for endpoint_i in upstream_endpoints_ind[0]]
    # Force river path to start at largest upstream mainstem node.
    # Trim river so it starts at mainstem, and if necessary extend river to mainstem node
    # and adjust network data to reflect changes
    rivers, head_rivpt, initial_riv_pts, head_node, next_rivpt, prev_rivpt, nearestnode_to_riv = _merge_riv_path_to_mainstem(rivers, head_rivpt, initial_riv_pts, nodes, nupstream, ndownstream, positions, nearestnode_to_riv, next_rivpt, prev_rivpt, affine)

    branch = 1
    next_branch = 2
    to_visit = []
    for riv_pt in initial_riv_pts:
        node = nearestnode_to_riv[riv_pt]
        node_i = nodes.index(node)
        G.nodes[node]['branches'] = {branch}
        to_visit.append((riv_pt, node_i, branch))
        branch = next_branch
        next_branch += 1

    # walk down river
    edits = defaultdict(lambda: defaultdict(dict))
    visited = set()
    outlets = set()
    while to_visit:
        (rivj, rivi), last_node_i, branch = to_visit.pop(0)
        xy = affine * (rivi,rivj)
        last_node = nodes[last_node_i]
        last_cell = cellid[last_node_i]
        visited.add(((rivj, rivi), last_node_i)) # include last_node_i so that different branches coming from different nodes can re-visit a node. but branches that are on the same node and same rivpts are dropped, since they will trace the same route

        next_node = nearestnode_to_riv[(rivj,rivi)]
        next_node_i = nodes.index(next_node)
        next_cell = cellid[next_node_i]

        xriv, yriv = xy
        xnode, ynode = positions[next_node_i]
        xp = xriv - xnode
        yp = yriv - ynode
        dist = np.sqrt(xp**2 + yp**2)
        halfres = .5 * resolution
        widehalfres = .6 * resolution
        # nearness determined by star-shape region around node. half resolution dist ensures all
        # gridlines get captured, but second terms chops out region in middle. makes diagonal
        # connections easier
        #riv_near_node = ((dist <= halfres) and
                         #(np.abs(xp/resolution * yp/resolution) <= 0.02))
        # or use a diamond shape, slighlty expand halfres to widen target at corners (but clip to halfres dist)
        #riv_near_node = ((dist <= halfres) and
                         #(yp + xp < widehalfres) and (yp - xp > -widehalfres) and
                         #(yp + xp > -widehalfres) and (yp - xp < widehalfres))
        # or a star with straight edges, gathers a bit more of the center than first star method
        riv_near_node = ((dist <= halfres) and
                         ((yp + 2*xp <  widehalfres) or (yp + xp/2 <  widehalfres/2)) and
                         ((yp - 2*xp > -widehalfres) or (yp - xp/2 > -widehalfres/2)) and
                         ((yp + 2*xp > -widehalfres) or (yp + xp/2 > -widehalfres/2)) and
                         ((yp - 2*xp <  widehalfres) or (yp - xp/2 <  widehalfres/2)))

        if ((next_node != last_node) and # moving to new node
                (riv_near_node)):        # helps reduce zig-zag
            if ((next_cell in edits) and
                   (last_cell in edits[next_cell]) and
                   (branch in edits[next_cell][last_cell]) and
                   (edits[next_cell][last_cell][branch] is not None)):
                # if water already rerouted from next_cell to last_cell, undo that and dont add new edit
                # fixes issue where meandering river goes back and forth between two nodes and second
                # rerouting overwrites the first. remove both to leave water where it was.
                print('Undo:', branch, ': cellid {0} to {1}'.format(next_cell, last_cell))
                G.remove_edge(next_node, last_node)
                G.nodes[last_node]['branches'].discard(branch)
                del edits[next_cell][last_cell][branch]
            else:
                # regular re-routing
                # but dont make link if another branch already goes from next_cell TO last_cell
                # this means(or, can happend when?) a new branch backtracks to previous cell
                # dont want to delete link since that branch should stay. just dont make new link
                # and branch will continue on elsewhere
                if not ((next_cell in edits) and
                        (last_cell in edits[next_cell]) and
                        (len(edits[next_cell][last_cell].values()) > 0) and
                        (branch not in edits[next_cell][last_cell]) and
                        (None not in edits[next_cell][last_cell].values())):

                    # create a new link
                    if 'branches' not in G.nodes[next_node]:
                        G.nodes[next_node]['branches'] = {branch}
                    else:
                        G.nodes[next_node]['branches'].add(branch)
                    print('New:', branch, ': cellid {0} to {1}'.format(last_cell, next_cell))
                    edits[last_cell][next_cell][branch] = True
                    for downstream_node in list(G.successors(last_node)):
                        # remove existing downstream links
                        downstream_cell = cellid[nodes.index(downstream_node)]
                        if (downstream_cell not in edits[last_cell]):
                            print('Remove(a):', branch, ': cellid {0} to {1}'.format(last_cell, downstream_cell))
                            G.remove_edge(last_node, downstream_node)
                            edits[last_cell][downstream_cell][branch] = None
                    if last_node in list(G.successors(next_node)):
                        # if changing direction, remove previous wrong-direction link
                        print('Remove(b):', branch, ': cellid {0} to {1}'.format(next_cell, last_cell))
                        G.remove_edge(next_node, last_node)
                        edits[next_cell][last_cell][branch] = None
                    G.add_edge(last_node, next_node)

                    # check to see if we just crossed a connection. if so, move existing crossed connection to next_cell as well. only an issue with diagonal fluxes
                    if ((abs(last_node[0]-next_node[0]) == 1) and
                        (abs(last_node[1]-next_node[1]) == 1)):
                        corner1_node = (next_node[0], last_node[1])
                        corner2_node = (last_node[0], next_node[1])
                        corner1_cell = cellid[nodes.index(corner1_node)]
                        corner2_cell = cellid[nodes.index(corner2_node)]
                        if (corner1_node in G.succ[corner2_node]):
                            G.remove_edge(corner2_node, corner1_node)
                            edits[corner2_cell][corner1_cell][branch] = None
                            G.add_edge(corner2_node, next_node)
                            edits[corner2_cell][next_cell][branch] = True # anything other than None
                        if (corner2_node in G.succ[corner1_node]):
                            G.remove_edge(corner1_node, corner2_node)
                            edits[corner1_cell][corner2_cell][branch] = None
                            G.add_edge(corner1_node, next_node)
                            edits[corner1_cell][next_cell][branch] = True # anything other than None

        elif (not riv_near_node): #or
              #(next_node in nx.ancestors(G, last_node))):
            # not stepping to next node, skip it by waiting until a different node is closest to river
            # next visited rivpt will still be attached to last_node
            next_node = last_node
            next_node_i = last_node_i
            next_cell = last_cell

        next_rivpts = next_rivpt[rivj,rivi]
        for rivj2,rivi2 in next_rivpts:
            if ((rivj2, rivi2), last_node_i) not in visited:
                # when rivers converg, keep both branches, could loose a convergence here if last node on one branch is too far away to jump before the branch is dropped
                if len(next_rivpts) > 1: # new branches
                    newbranch = next_branch
                    next_branch += 1
                else:
                    newbranch = branch # really same branch
                to_visit.append(((rivj2, rivi2), next_node_i, newbranch)) # first branch stays the same
        if ((len(next_rivpt[rivj,rivi]) == 0) and
            (dist_to_coast[nodes.index(nearestnode_to_riv[rivj,rivi])] <= 3)): # node units
            # no downstream points, AND CLOSE TO COAST, remove downstream flow from node
            # dont do this for likely upstream points, just leave existing connections
            # helps with errors in osm_river, dont want to strand water upstream
            # next_node is next if we just moved to new one, or last_node if we didn't
            outlets.add(next_cell)
            for node2 in list(G.successors(next_node)):
                G.remove_edge(next_node, node2)
                node2_cell = cellid[nodes.index(node2)]
                edits[next_cell][node2_cell][branch] = None

    with open(str(target[0]), 'w', newline='') as fout:
        csvwriter = csv.writer(fout)
        for (from_cell, to_cells) in edits.items():
            downstream = 0
            for to_cell, branches in to_cells.items():
                branch_found = False
                for branch, val in branches.items():
                    if val: # either None or True
                        branch_found = True
                if branch_found:
                    downstream += 1
            if downstream:
                frac = 1/downstream
            for to_cell, branches in to_cells.items():
                from_node = nodes[cellid.index(from_cell)]
                to_node = nodes[cellid.index(to_cell)]
                successors = list(Gorig.successors(from_node))
                if list(branches.values()) == [None]: # zero out for removed links
                    csvwriter.writerow([from_cell, to_cell, 0])
                else:
                    if ((to_node in successors) and (frac != 1)): # zero out if single downlink changing to multiple
                        csvwriter.writerow([from_cell, to_cell, 0])
                    if ((to_node not in successors) or # record new link
                        (frac != 1)): # record if fractional flux
                        csvwriter.writerow([from_cell, to_cell, frac])

    for node in G.nodes():
        G.node[node]['upstream'] = len(nx.ancestors(G, node))
        G.node[node]['downstream'] = len(nx.descendants(G, node))
    nx.write_yaml(G, str(target[1]))

    with rasterio.open(str(target[2]), 'w', **meta) as rast:
        rast.write(rivers, 1)

    with open(str(target[3]), 'w') as fout:
        for outlet in sorted(outlets):
            fout.write(str(outlet)+'\n')

    with open(str(target[4]), 'wb') as fout:
        pickle.dump(segments, fout)
    with open(str(target[5]), 'wb') as fout:
        pickle.dump(next_rivpt, fout)

    return 0
