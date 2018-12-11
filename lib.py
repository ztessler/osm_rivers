import csv
import pickle
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import palettable
import geopandas
import shapely.geometry as sgeom
import rasterio
import rasterio.features as rfeatures
import cartopy.crs as ccrs
import networkx as nx
import skimage.morphology as morph
import skimage.draw
from scipy.ndimage.filters import generic_filter
from netCDF4 import Dataset
import pyproj
import itertools
from collections import defaultdict, Counter


def set_projection(source, target, env):
    delta_ll = geopandas.read_file(str(source[0]))
    lon0, lat0 = np.array(delta_ll.dissolve(by='Delta').convex_hull.centroid.squeeze())
    laea = ccrs.LambertAzimuthalEqualArea(central_longitude = lon0,
                                          central_latitude = lat0)
    with open(str(target[0]), 'w') as fout:
        fout.write(laea.proj4_init)
    return 0


def project_and_clip_osm_rivers(source, target, env):
    rivers_ll = geopandas.read_file(str(source[0]))
    delta_ll = geopandas.read_file(str(source[1]))
    with open(str(source[2]), 'r') as fin:
        proj4_init = fin.read()

    deltahull_ll = geopandas.GeoDataFrame(delta_ll.dissolve(by='Delta').convex_hull, columns=['geometry'])
    deltahull_ll.crs = delta_ll.crs

    rivers = rivers_ll.to_crs(proj4_init)
    delta = delta_ll.to_crs(proj4_init)
    deltahull = deltahull_ll.to_crs(proj4_init)
    deltahull = geopandas.GeoDataFrame(deltahull.buffer(25000), columns=['geometry'], crs=deltahull.crs)

    rivers_clip = geopandas.overlay(rivers, deltahull, how='intersection') #slow
    rivers_clip.to_file(str(target[0]))
    return 0


def project_and_clip_osm_waterways(source, target, env):
    rivers_ll = geopandas.read_file(str(source[0]))
    delta_ll = geopandas.read_file(str(source[1]))
    with open(str(source[2]), 'r') as fin:
        proj4_init = fin.read()

    deltahull_ll = geopandas.GeoDataFrame(delta_ll.dissolve(by='Delta').convex_hull, columns=['geometry'])
    deltahull_ll.crs = delta_ll.crs

    rivers = rivers_ll.to_crs(proj4_init)
    delta = delta_ll.to_crs(proj4_init)
    deltahull = deltahull_ll.to_crs(proj4_init)
    deltahull = geopandas.GeoDataFrame(deltahull.buffer(25000), columns=['geometry'], crs=deltahull.crs)

    rivgeom = rivers['geometry']
    rivers_clip = rivgeom.intersection(deltahull['geometry'].item())
    rivers['geometry'] = rivers_clip
    rivers = rivers[rivers.length > 0]

    rivers.to_file(str(target[0]), encoding='utf-8')
    return 0


def project_and_clip_coastline(source, target, env):
    coast_ll = geopandas.read_file(str(source[0]))
    delta_ll = geopandas.read_file(str(source[1]))
    with open(str(source[2]), 'r') as fin:
        proj4_init = fin.read()

    deltahull_ll = geopandas.GeoDataFrame(delta_ll.dissolve(by='Delta').convex_hull, columns=['geometry'], crs=delta_ll.crs)

    deltahull3_ll = geopandas.GeoDataFrame(deltahull_ll.buffer(3), columns=['geometry'], crs=deltahull_ll.crs)
    coast_ll_clip = geopandas.overlay(coast_ll, deltahull3_ll, how='intersection')
    coast = coast_ll_clip.to_crs(proj4_init)
    deltahull = deltahull_ll.to_crs(proj4_init)

    coastbuff_poly = coast.buffer(2000).unary_union
    #if isinstance(coastbuff_poly, sgeom.MultiPolygon):
        #i = np.argmax([p.area for p in coastbuff_poly])
        #coastbuff_poly = coastbuff_poly[i]
    coastline = coastbuff_poly.buffer(-2000).boundary

    deltahull_buff = deltahull.buffer(50000).unary_union
    coastline_clip = coastline.intersection(deltahull_buff)

    geopandas.GeoSeries(coastline_clip).to_file(str(target[0]))
    return 0


def clip_osm_rivers(source, target, env):
    rivers_ll = geopandas.read_file(str(source[0]))
    delta_ll = geopandas.read_file(str(source[1]))

    deltahull_ll = geopandas.GeoDataFrame(delta_ll.dissolve(by='Delta').convex_hull, columns=['geometry'])
    deltahull_ll.crs = delta_ll.crs

    rivers_clip = geopandas.overlay(rivers_ll, deltahull_ll, how='intersection')

    rivers_clip.to_file(str(target[0]))
    return 0


def select_waterway_rivers(source, target, env):
    rivers = geopandas.read_file(str(source[0]))
    rivers = rivers[rivers['fclass'] == 'river']

    rivers.to_file(str(target[0]), encoding='utf-8')
    return 0


def filter_waterway_rivers(source, target, env):
    rivers = geopandas.read_file(str(source[0]))

    minwaterway = env.get('minwaterway', 0)
    rivers = rivers[~(rivers['name'].isnull()) & (rivers.length > minwaterway)]

    rivers.to_file(str(target[0]), encoding='utf-8')
    return 0


def get_river_widths(source, target, env):
    rivers = geopandas.read_file(str(source[0]))
    rivers = rivers[rivers['fclass'] == 'river']

    # estimate width of each segment using area=w*l, perimeter=2w+2l, solve with quadratic formula
    simplified = rivers['geometry'].simplify(0.001)
    a = simplified.area
    p = simplified.boundary.length
    #x1 = (a + np.sqrt(p**2 - 16*a)) / 4 # larger, will be length
    width = (a - np.sqrt(p**2 - 16*a)) / 4 # smaller, will be width
    #= rivers[rivers['width'] > 0]
    rivers['width_est'] = width
    #rivers = rivers[width > 100]

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
    rivers['geometry'] = rivers.buffer(-thinning).buffer(thinning)
    rivers = rivers[rivers.area>0]


    # fill holes, braided rivers
    #minhole = env.get('minhole', 0)
    # just fill all holes. might need to leave large ones...
    newgeoms = []
    for i, poly in rivers['geometry'].items():
        if isinstance(poly, sgeom.MultiPolygon):
            poly2 = sgeom.MultiPolygon([sgeom.Polygon(p.exterior) for p in poly])
        else:
            poly2 = sgeom.Polygon(poly.exterior)
        newgeoms.append(poly2)
    #rivers.loc[i, 'geometry'] = poly2
    rivers['geometry'] = geopandas.GeoSeries(newgeoms, index=rivers.index)

    #rivers_merge = geopandas.overlay(geopandas.GeoDataFrame(rivers, columns=['geometry']), geopandas.GeoDataFrame(rivers, columns=['geometry']), how='union')
    #rivers_merge.crs = rivers.crs

    rivers.to_file(str(target[0]))
    return 0


def merge_water_waterway_vecs(source, target, env):
    water = geopandas.read_file(str(source[0]))
    waterways = geopandas.read_file(str(source[1]))

    buff = env.get('buff', 100)
    waterway_polys = geopandas.GeoDataFrame(waterways.buffer(buff), columns=['geometry'], crs=waterways.crs)

    merged = geopandas.overlay(water, waterway_polys, how='union')
    merged.to_file(str(target[0]))
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
                newj = None
                if ((len(todelete) == 2 and rivers[j,i] == 3) or
                        (len(todelete) == 3 and rivers[j,i] == 4)):
                    #two (or three) terminating branches
                    #replace with single line to average location
                    ## keep longer segment, will get appended to last segment
                    #seglens = [len(segment) for segment in todelete]
                    #longest_i = np.argmax(seglens)
                    #todelete.pop(longest_i)
                    newj = int(round(sum([segtodelete[-1][0] for segtodelete in todelete]) / len(todelete)))
                    newi = int(round(sum([segtodelete[-1][1] for segtodelete in todelete]) / len(todelete)))
                rivers[j,i] -= len(todelete) # old bifur point becomes normal river, or a four-way becomes three-way
                for segment in todelete:
                    for rivpt in segment: # bifur point not included on segment
                        rivers[rivpt] = 0
                if newj is not None:
                    rowidx, colidx = skimage.draw.line(j, i, newj, newi)
                    rivers[rowidx, colidx] = 2
                    rivers[newj, newi] = 1

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
                ll = affine * (i+.5, j+.5)
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
                        ll2 = affine * (i2+.5, j2+.5)
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

    upstream = {node: Gclip.node[node]['upstream'] for node in Gclip.nodes()}
    downstream = {node: Gclip.node[node]['downstream'] for node in Gclip.nodes()}
    positions = {node: Gclip.node[node]['xy'] for node in Gclip.nodes()}
    with open(str(target[2]), 'wb') as fout:
        pickle.dump(upstream, fout)
    with open(str(target[3]), 'wb') as fout:
        pickle.dump(downstream, fout)
    with open(str(target[4]), 'wb') as fout:
        pickle.dump(positions, fout)
    return 0


def find_bifurs(source, target, env):
    with rasterio.open(str(source[0]), 'r') as rast:
        rivers = rast.read(1)
        meta = rast.meta.copy()

    # just in case this is run on a bifur_grid
    rivers[rivers != 0] = 1

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

    with open(str(source[1]), 'rb') as fin:
        segments = pickle.load(fin)

    mpl.style.use('ggplot')
    fig, ax = plt.subplots(1,1, figsize=(8, 12))#, dpi=300)

    for segi, segment in segments.items():
        x1, y1 = affine * segment[0][::-1]
        x2, y2 = affine * segment[-1][::-1]
        # flows from segment[0] to segment[-1]
        ax.annotate("", xy=(x2,y2), xytext=(x1,y1),
                arrowprops=dict(facecolor='k', edgecolor='k', arrowstyle='-|>'))
        ax.text((x2+x1)/2, (y2+y1)/2, segi)

    I, J = np.meshgrid(np.arange(bifurs.shape[1]), np.arange(bifurs.shape[0]))
    xs, ys = affine * (I.flatten(), J.flatten())
    X = xs.reshape(I.shape)
    Y = ys.reshape(J.shape)

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


def find_nearest_nodes_to_riv(source, target, env):
    G = nx.read_yaml(str(source[0]))
    with rasterio.open(str(source[1]), 'r') as rast:
        rivers = rast.read(1)
        affine = rast.transform

    nodes = [node for node in G.nodes()]
    positions = [G.node[node]['xy'] for node in nodes]
    #nupstream = [G.node[node]['upstream'] for node in nodes]
    #ndownstream = [G.node[node]['downstream'] for node in nodes]

    allbasinmask = [1 for n in nodes]

    nearestnode = {}
    #nearestnode_ndownstream = {} # dict on riv points, lower ndownstream means closer to coast. use to assess flow dir on river
    #nearestnode_nupstream = {} # dict on riv points, lower ndownstream means closer to coast. use to assess flow dir on river
    wet = np.where(rivers>0)
    for j, i in zip(*wet):
        xy = affine * (i,j)
        nearest_node_i = _find_nearest_node_i(xy, positions, allbasinmask)
        nearestnode[(j,i)] = nodes[nearest_node_i]
        #nearestnode_ndownstream[(j,i)] = ndownstream[nearest_node_i]
        #nearestnode_nupstream[(j,i)] = nupstream[nearest_node_i]
        # if these are needed, just get with ndownstream[nodes.index(node)]

    with open(str(target[0]), 'wb') as fout:
        pickle.dump(nearestnode, fout)
    return 0


def calc_dist_to_coast(source, target, env):
    G = nx.read_yaml(str(source[0]))
    with rasterio.open(str(source[1]), 'r') as rast:
        rivers = rast.read(1)
        affine = rast.transform
    coastline = geopandas.read_file(str(source[2])).loc[0,'geometry']

    # calc node_dist_to_coast
    nodes = list(G.nodes())
    positions = geopandas.GeoSeries([sgeom.Point(G.node[node]['xy']) for node in nodes], index=nodes)
    node_dist_to_coast = positions.distance(coastline)

    # calc riv_dist_to_coast
    wet = np.where(rivers > 0)
    positions = geopandas.GeoSeries([sgeom.Point(affine * (i,j)) for j,i in zip(*wet)], index=zip(*wet))
    riv_dist_to_coast = positions.distance(coastline)

    node_dist_to_coast.to_pickle(str(target[0]))
    riv_dist_to_coast.to_pickle(str(target[1]))
    return 0


def find_river_segments(source, target, env):
    with rasterio.open(str(source[0]), 'r') as rast:
        rivers = rast.read(1)
        affine = rast.transform

    endpoints = np.where(rivers == 1)
    # just start at the first one. all rivers are connected, so will get everywhere

    branchpoints = set()
    for j,i in zip(*np.where(rivers==3)):
        for dj in [-1,0,1]:
            for di in [-1,0,1]:
                if dj == di == 0:
                    continue
                if rivers[j+dj, i+di] > 0:
                    branchpoints.add((j+dj,i+di))

    (j,i) = endpoints[0][0], endpoints[1][0]
    cursegi = 0
    tovisit = [(j,i,cursegi)]
    onfullsegment = set()
    segments = defaultdict(list)
    fullsegments = {}
    maxsegi = cursegi
    while tovisit:
        j, i, cursegi = tovisit.pop(0)
        segments[cursegi].append((j,i))
        if (rivers[j,i] == 3) or ((rivers[j,i] == 1) and (len(segments[cursegi])>1)):
            # found segment end
            # mark points on this segment as complete
            for (_j, _i) in segments[cursegi]:
                onfullsegment.add((_j,_i))
            fullsegments[cursegi] = segments[cursegi]

        for (dj, di) in [(-1,0),(0,1),(1,0),(0,-1),(-1,-1),(1,-1),(1,1),(-1,1)]: # list with horizontal and vert before diag, in case of ambig route: example: bifur to right and a branch from that to down. down could be jumped to diagonally, so check hor/vert first and break if bifur is found
                if dj == di == 0:
                    continue
                j2 = j + dj
                i2 = i + di

                if rivers[j2,i2] > 0:
                    if (rivers[j,i]<3) and ((j,i) in branchpoints) and ((j2,i2) in branchpoints):
                        # don't leap to another branch by accident
                        continue
                    if (j2, i2) in segments[cursegi]:
                        # dont go backwards up river
                        continue
                    if (rivers[j2,i2]<3) and ((j2, i2) in onfullsegment):
                        # point has already been assigned to a full segment (bifurs can be on multiple)
                        continue

                    # found next river point
                    if rivers[j,i] == 3: # currently on bifur point
                        maxsegi += 1 # each branch gets new cursegi
                        nextsegi = maxsegi
                        segments[nextsegi].append((j,i)) # add current bifur point as start of new segment
                        tovisit.append((j2, i2, nextsegi)) # add new segments to end of queue
                    else:
                        tovisit.insert(0, (j2, i2, cursegi)) # on current segment, add to front of queue to stay on this branch
                    if rivers[j,i] == 2:
                        # regular path, not bifur, only one neighbor. break so as not to find incorrect diag branch
                        break

    with open(str(target[0]), 'wb') as fout:
        pickle.dump(fullsegments, fout)
    return 0


def set_segment_flowdir(source, target, env):
    with open(str(source[0]), 'rb') as fin:
        segments = pickle.load(fin)
    with rasterio.open(str(source[1])) as rast:
        rivers = rast.read(1)
        affine = rast.transform
    waterways = geopandas.read_file(str(source[2]))
    with open(str(source[3]), 'rb') as fin:
        riv_dist_to_coast = pickle.load(fin)
    lines = waterways['geometry']

    #for each segment
    counts = Counter()
    scores = Counter()
    directed_segments = {}
    x1, y1 = affine * (rivers.shape[1]//2, rivers.shape[0]//2)
    x2, y2 = affine * (rivers.shape[1]//2 + 1, rivers.shape[0]//2 + 1)
    dx = abs(x2-x1)
    dy = abs(y2-y1)
    pixelsize = np.sqrt(dx**2 + dy**2)
    for segi in sorted(segments):
        segment = segments[segi]
        counts.clear()
        scores.clear()
        # for each point on segment
        for (j,i) in segment:
            x, y = affine * (i+.5, j+.5)
            # find nearest waterway
            pt = sgeom.Point(x, y)
            dists = [pt.distance(line) for line in lines]
            ind = np.argmin(dists)
            scores[ind] += 1/max(dists[ind], pixelsize) # inverse distance weight so closest points count most, but dont go closer than nominal resolution since one very very close line could blow up comparison
            counts[ind] += (dists[ind] < (2*pixelsize))
            #TODO clip distances so that past a certain dist (pixelsize*3??) just add zero. require total score to be 3 of these, if not then use dist_to_coast metric
        # get direction of most common waterway
        ind, score = scores.most_common(1)[0]
        count = counts[ind]
        use_dist_to_coast = False
        if count >= min(7, len(segments)):
            # have enough overlap between osm river waterway and segment, use it's direction
            line = lines.iloc[ind]
            # check if segment direction is correct or revered
            # for each verticies from start of line to end
            j_seg_start, i_seg_start = segment[0]
            j_seg_end, i_seg_end = segment[-1]
            x_seg_start, y_seg_start = affine * (i_seg_start, j_seg_start)
            x_seg_end, y_seg_end = affine * (i_seg_end, j_seg_end)
            mindist_to_start = np.inf
            mindist_to_end = np.inf
            nearest_verti_to_start = None
            nearest_verti_to_end = None
            for verti, (x, y) in enumerate(line.coords):
                # calc dists, track mins, to start and end of segment
                start_dist = np.sqrt((x - x_seg_start)**2 + (y - y_seg_start)**2)
                end_dist = np.sqrt((x - x_seg_end)**2 + (y - y_seg_end)**2)
                if start_dist < mindist_to_start:
                    mindist_to_start = start_dist
                    nearest_verti_to_start = verti
                if end_dist < mindist_to_end:
                    mindist_to_end = end_dist
                    nearest_verti_to_end = verti
            if nearest_verti_to_start < nearest_verti_to_end:
                # mindist to start comes earlier, then directions match, do nothing
                directed_segments[segi] = segment
                print('Segment {0}: {1} to {2} (no change)'.format(segi, segment[0], segment[-1]))
            elif nearest_verti_to_end < nearest_verti_to_start:
                # mindist to end comes earlier, then swap segment direction
                directed_segments[segi] = segment[::-1]
                print('Segment {0}: {1} to {2} (reversed)'.format(segi, segment[-1], segment[0]))
            else:
                print('  Segment {0}, not enough line/segment overlap, use dist to coast'.format(segi))
                use_dist_to_coast = True
        if (count < min(7, len(segments))) or use_dist_to_coast:
            # use dist_to_coast since not enough segment pixels align with an osm waterway
            if ((riv_dist_to_coast[segment[0]] < riv_dist_to_coast[segment[-1]])): # and
                    #(riv_dist_to_coast[segment[0]] < 10000) and
                    #(riv_dist_to_coast[segment[-1]] > 10000)): # replace with node resolution?
                directed_segments[segi] = segment[::-1]
                print('Segment {0}: {1} to {2} (reversed (b))'.format(segi, segment[-1], segment[0]))
            else:
                directed_segments[segi] = segment
                print('Segment {0}: {1} to {2} (no change (b))'.format(segi, segment[0], segment[-1]))

    with open(str(target[0]), 'wb') as fout:
        pickle.dump(directed_segments, fout)
    return 0


def remove_small_loops(source, target, env):
    with rasterio.open(str(source[0])) as rast:
        rivers = rast.read(1)
        meta = rast.meta
    with open(str(source[1]), 'rb') as fin:
        segments = pickle.load(fin)
    minlen = env.get('minlen', 10)

    for segi, segj in itertools.combinations(segments.keys(), 2):
        seg1 = segments[segi]
        seg2 = segments[segj]
        if (((seg1[0] == seg2[0]) and (seg1[-1] == seg2[-1]) and (seg1 != seg2)) or
                ((seg1[0] == seg2[-1]) and (seg1[-1] == seg2[0]) and (seg1 != seg2[::-1]))):
            # parallel paths
            if (seg1[0] == seg2[0]):
                samedir = True
            else:
                samedir = False
            if (len(seg1) < minlen) and (len(seg2) < minlen):
                # both segs short, delete both and replace with straight line
                for seg in [seg1, seg2]:
                    for j,i in seg[1:-1]:
                        rivers[j,i] = 0
                rowidx, colidx = skimage.draw.line(*seg1[0], *seg1[-1])
                rivers[rowidx, colidx] = 1
            elif len(seg1) < minlen:
                for j,i in seg1[1:-1]:
                    rivers[j,i] = 0
            elif len(seg2) < minlen:
                for j,i in seg2[1:-1]:
                    rivers[j,i] = 0

    with rasterio.open(str(target[0]), 'w', **meta) as out:
        out.write(rivers, 1)
    return 0


def next_prev_pts(source, target, env):
    with open(str(source[0]), 'rb') as fin:
        segments = pickle.load(fin)
    next_rivpts = defaultdict(list)
    prev_rivpts = defaultdict(list)
    for segment in segments.values():
        prevj = None
        previ = None
        for thisj, thisi in segment:
            if prevj is not None:
                next_rivpts[prevj, previ].append((thisj, thisi))
                prev_rivpts[thisj, thisi].append((prevj, previ))
            prevj = thisj
            previ = thisi

    with open(str(target[0]), 'wb') as fout:
        pickle.dump(next_rivpts, fout)
    with open(str(target[1]), 'wb') as fout:
        pickle.dump(prev_rivpts, fout)
    return 0


def find_head_rivpt(source, target, env):
    with rasterio.open(str(source[0])) as rast:
        rivers = rast.read(1)
        affine = rast.transform
    with open(str(source[1]), 'rb') as fin:
        nearestnode = pickle.load(fin)
    with open(str(source[2]), 'rb') as fin:
        ndownstream = pickle.load(fin)
    with open(str(source[3]), 'rb') as fin:
        positions = pickle.load(fin)

    maxdown = 0
    for node, ndown in ndownstream.items():
        if ndown > maxdown:
            maxdown = ndown
            maxdownnode = node
    nodex, nodey = positions[maxdownnode]

    endpoints = np.where(rivers==1)
    mindist = np.inf
    for (j,i) in zip(*endpoints):
        x, y = affine * (i, j)
        dist = np.sqrt((nodex-x)**2 + (nodey-y)**2)
        if dist < mindist:
            mindist = dist
            head_rivpt = (j,i)

    with open(str(target[0]), 'wb') as fout:
        pickle.dump(head_rivpt, fout)
    return 0


def merge_riv_path_to_mainstem(source, target, env):
    with rasterio.open(str(source[0])) as rast:
        rivers = rast.read(1)
        affine = rast.transform
        meta = rast.meta
    with open(str(source[1]), 'rb') as fin:
        head_rivpt = pickle.load(fin)
    with open(str(source[2]), 'rb') as fin:
        next_rivpt = pickle.load(fin)
    with open(str(source[3]), 'rb') as fin:
        prev_rivpt = pickle.load(fin)
    with open(str(source[4]), 'rb') as fin:
        nearestnode = pickle.load(fin)
    with open(str(source[5]), 'rb') as fin:
        nupstream = pickle.load(fin)
    with open(str(source[6]), 'rb') as fin:
        ndownstream = pickle.load(fin)
    with open(str(source[7]), 'rb') as fin:
        positions = pickle.load(fin)

    nodes = sorted(nupstream) # sorting keys of dict
    nupstream = [nupstream[node] for node in nodes]
    ndownstream = [ndownstream[node] for node in nodes]
    positions = [positions[node] for node in nodes]

    def _trim_river(rivers, rivpt, next_rivpt, prev_rivpt):
        # erase river above rivpt, and adjust next_rivpt and prev_rivpt dicts
        thispt = rivpt
        rivers[thispt] = 1
        to_visit = [(thispt, uppt) for uppt in prev_rivpt[thispt]]
        while to_visit:
            downpt, thispt = to_visit.pop()
            rivers[thispt] = 0
            if thispt in prev_rivpt:
                to_visit.extend([(thispt, uppt) for uppt in prev_rivpt[thispt]])
        return rivers

    def _extend_river(rivers, mindist_rivpt, head_node_ij):
        # draw new river in straight line from closest river approach to head_node
        head_node_ji = tuple(head_node_ij[::-1])
        rowidx, colidx = skimage.draw.line(*mindist_rivpt, *head_node_ji)
        rivers[rowidx, colidx] = 2
        rivers[head_node_ji] = 1
        return rivers

    score = np.array(nupstream) * np.array(ndownstream) # maximized at "center" of mainstem, which we assume is upstream of the delta boundary
    head_node_i = np.argmax(score)
    head_node = nodes[head_node_i]
    nearest = nearestnode[head_rivpt]

    if head_node == nearest:
        with rasterio.open(str(target[0]), 'w', **meta) as rast:
            rast.write(rivers, 1)
        return 0

    found_head_node_on_riv = False
    for rivpt, node in nearestnode.items():
        if node == head_node:
            found_head_node_on_riv = True
            break
    if found_head_node_on_riv:
        rivers = _trim_river(rivers, rivpt, next_rivpt, prev_rivpt)
    else:
        to_visit = [head_rivpt]
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
        if mindist_rivpt != head_rivpt:
            rivers = _trim_river(rivers, mindist_rivpt, next_rivpt, prev_rivpt)
        rivers = _extend_river(rivers, mindist_rivpt, head_node_ij)

    with rasterio.open(str(target[0]), 'w', **meta) as rast:
        rast.write(rivers, 1)
    return 0


def remap_riv_network(source, target, env):
    G = nx.read_yaml(str(source[0]))
    Gorig = G.copy()
    with rasterio.open(str(source[1]), 'r') as rast:
        rivers = rast.read(1)
        affine = rast.transform
        meta = rast.meta.copy()
    with open(str(source[2]), 'rb') as fin:
        head_rivpt = pickle.load(fin)
    with open(str(source[3]), 'rb') as fin:
        next_rivpt = pickle.load(fin)
    with open(str(source[4]), 'rb') as fin:
        prev_rivpt = pickle.load(fin)
    with open(str(source[5]), 'rb') as fin:
        nearestnode = pickle.load(fin)
    with open(str(source[6]), 'rb') as fin:
        dist_to_coast = pickle.load(fin)

    nodes = [node for node in G.nodes()]
    positions = [G.node[node]['xy'] for node in nodes]
    nupstream = [G.node[node]['upstream'] for node in nodes]
    ndownstream = [G.node[node]['downstream'] for node in nodes]
    cellid = [G.node[node]['cellid'] for node in nodes]

    p0 = positions[0]
    resolution = np.min([np.sqrt((p0[0]-p[0])**2 + (p0[1]-p[1])**2) for p in positions[1:]])

    minx, maxy = affine * (0,0)
    maxx, miny = affine * (rivers.shape[1], rivers.shape[0])

    endpoints = np.where(rivers == 1)
    upstream_endpoints = []
    for j,i in zip(*endpoints):
        if not prev_rivpt[j,i]:
            upstream_endpoints.append((j,i))

    branch = 1
    next_branch = 2
    to_visit = []
    for riv_pt in upstream_endpoints:
        # dont start on a node, wait to snap to one
        to_visit.append((riv_pt, None, branch))
        branch = next_branch
        next_branch += 1

    # walk down river
    edits = defaultdict(lambda: defaultdict(dict))
    visited = set()
    outlets = set()
    while to_visit:
        (rivj, rivi), last_node_i, branch = to_visit.pop(0)
        xy = affine * (rivi,rivj)
        if last_node_i is not None:
            last_node = nodes[last_node_i]
            last_cell = cellid[last_node_i]
            visited.add(((rivj, rivi), last_node_i)) # include last_node_i so that different branches coming from different nodes can re-visit a node. but branches that are on the same node and same rivpts are dropped, since they will trace the same route
        else:
            last_node = None
            last_cell = None

        next_node = nearestnode[(rivj,rivi)]
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
                if not ((last_node is None) or
                        ((next_cell in edits) and
                         (last_cell in edits[next_cell]) and
                         (len(edits[next_cell][last_cell].values()) > 0) and
                         (branch not in edits[next_cell][last_cell]) and
                         (None not in edits[next_cell][last_cell].values()))):

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
                        (abs(last_node[1]-next_node[1]) == 1) and
                        ((next_node[0], last_node[1]) in nodes) and
                        ((last_node[0], next_node[1]) in nodes)):
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
            (dist_to_coast[nearestnode[rivj,rivi]] <= 1) and # node units
            (next_node is not None)):
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

    with open(str(target[2]), 'w') as fout:
        for outlet in sorted(outlets):
            fout.write(str(outlet)+'\n')

    return 0
