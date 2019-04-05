import csv
import pickle
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import palettable
import pandas
import geopandas
import shapely.geometry as sgeom
import shapely.ops as sops
import rasterio
import rasterio.features as rfeatures
import fiona
import cartopy.crs as ccrs
import networkx as nx
import skimage.morphology as morph
import skimage.draw
from scipy.ndimage.filters import generic_filter, gaussian_filter
from scipy.interpolate import griddata
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


def find_grwl_list(source, target, env):
    delta = geopandas.read_file(str(source[0]))
    boundary = delta.buffer(1).unary_union.boundary
    inbounds = []
    for grwlfile in [str(s) for s in source[1:env['nsources']]]:
        grwl = fiona.open(grwlfile)
        gminx, gminy, gmaxx, gmaxy = grwl.bounds
        footprint = sgeom.Polygon([(gminx, gminy), (gminx, gmaxy), (gmaxx, gmaxy), (gmaxx, gminy), (gminx, gminy)])
        if footprint.intersects(boundary):
            inbounds.append(grwlfile)
    with open(str(target[0]), 'w') as fout:
        for fname in inbounds:
            fout.write(fname + '\n')
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
    coastpoly = coastbuff_poly.buffer(-2000)
    coastline = coastpoly.boundary

    deltahull_buff = deltahull.buffer(50000).unary_union
    coast_clip = coastpoly.intersection(deltahull_buff)
    coastline_clip = coastline.intersection(deltahull_buff)

    geopandas.GeoSeries(coastline_clip, crs=proj4_init).to_file(str(target[0]))
    geopandas.GeoSeries(coast_clip, crs=proj4_init).to_file(str(target[1]))
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


def filter_river_types(source, target, env):
    rivers = geopandas.read_file(str(source[0]))
    if not env.get('wetlands', True):
        # keep unspecified, reservoirs, and rivers and river classes. lakes are often unspec.
        rivers = rivers[(rivers['code']==8200) | (rivers['code']==8201) | (rivers['code']==8202)]
    rivers.to_file(str(target[0]))
    return 0


def get_river_widths(source, target, env):
    rivers = geopandas.read_file(str(source[0]))

    def est_width(s):
        shps = s['geometry']
        if not isinstance(shps, sgeom.MultiPolygon):
            shps = [shps]
        widths = []
        for shp in shps:
            a = shp.area
            p = shp.length
            for interior in shp.interiors:
                hole = sgeom.Polygon(interior)
                p -= hole.length
            widths.append((p - np.sqrt(p**2 - 16*a)) / 4) # - term is smaller, is width
        return np.mean(widths)

    # estimate width of each segment using area=w*l, perimeter=2w+2l, solve with quadratic formula
    #simplified = rivers['geometry'].simplify(0.001)
    #a = simplified.area
    #p = simplified.boundary.length
    #x1 = (p + np.sqrt(p**2 - 16*a)) / 4 # larger, will be length
    #width = (p - np.sqrt(p**2 - 16*a)) / 4 # smaller, will be width
    #= rivers[rivers['width'] > 0]
    #rivers['width_est'] = width
    #rivers = rivers[width > 100]

    rivers['width_est'] = rivers.apply(est_width, axis=1)

    rivers.to_file(str(target[0]))
    return 0


def thin_vec(source, target, env):
    rivers = geopandas.read_file(str(source[0]))


    # drop thin stuff
    thinning = env.get('thinning', 100) # removes thin portions of geoms
    minwidth = env.get('minwidth', 0)   # removes whole rivers geoms
    #rivers = rivers[rivers.width_est > minwidth]
    rivers['geometry'] = rivers.buffer(-thinning/2).buffer(thinning/2)
    rivers = rivers[(rivers.area>0) & (rivers.width_est > minwidth)]

    # drop remaining small area stuff
    minarea = env.get('minarea', 0)
    rivers = rivers[rivers.area > minarea]

    # fill holes, braided rivers
    #minhole = env.get('minhole', 0)
    # just fill all holes. might need to leave large ones...
    #newgeoms = []
    #for i, poly in rivers['geometry'].items():
        #if isinstance(poly, sgeom.MultiPolygon):
            #poly2 = sgeom.MultiPolygon([sgeom.Polygon(p.exterior) for p in poly])
        #else:
            #poly2 = sgeom.Polygon(poly.exterior)
        #newgeoms.append(poly2)
    ##rivers.loc[i, 'geometry'] = poly2
    #rivers['geometry'] = geopandas.GeoSeries(newgeoms, index=rivers.index)

    #rivers_merge = geopandas.overlay(geopandas.GeoDataFrame(rivers, columns=['geometry']), geopandas.GeoDataFrame(rivers, columns=['geometry']), how='union')
    #rivers_merge.crs = rivers.crs

    rivers.to_file(str(target[0]))
    return 0


def merge_water_waterway_vecs(source, target, env):
    water = geopandas.read_file(str(source[0]))
    unthinned_water = geopandas.read_file(str(source[1]))
    waterways = geopandas.read_file(str(source[2]))

    # remove any waterways that mostly overlap a polygon that is in unthinned_water but not water
    # these waterways have already been removed in the polygon layer (due to width), so don't add them back
    #import ipdb;ipdb.set_trace()
    # probably very slow. think about how to speed up
    #import time
    #t = time.time()
    #print('Running...')
    #removed_polys = geopandas.overlay(unthinned_water, water, how='difference').unary_union()
    #waterways = waterways.difference(removed_polys)
    #print('Took {} seconds'.format(t-time.time()))
    unthinned_index = unthinned_water.sindex
    thinned_index = water.sindex
    goodlines = []
    filtered_waterways = waterways.copy()
    for i, line in waterways.iterrows():
        unthinned_maybe_ind = list(unthinned_index.intersection(line['geometry'].bounds))
        if not unthinned_maybe_ind:
            # line is not in original water bodies layer, so keep
            continue
        unthinned_maybe = unthinned_water.iloc[unthinned_maybe_ind]
        unthinned_matches = unthinned_maybe[unthinned_maybe['geometry'].intersects(line['geometry'])]

        thinned_maybe_ind = list(thinned_index.intersection(line['geometry'].bounds))
        if not thinned_maybe_ind:
            # line was in original, but none left in thinned. remove from waterways too
            filtered_waterways.drop(i, inplace=True)
            continue
        thinned_maybe = water.iloc[thinned_maybe_ind]
        thinned_matches = thinned_maybe[thinned_maybe['geometry'].intersects(line['geometry'])]
        if thinned_matches.empty:
            # line was in original, but none left in thinned. remove from waterways too
            filtered_waterways.drop(i, inplace=True)
            continue

        # remove segments that are inside original water polygons. either they were removed during
        # thinning, in which we remove them here too, or they were kept, in which case we already
        # either way, we dont need the overlap
        # but we DO need to keep segments that are outside of unthinned water bodies, since those
        # are new river segs
        segments_to_remove = []
        for j, unthinned_match in unthinned_matches.iterrows():
            overlap = line['geometry'].intersection(unthinned_match['geometry'])
            if isinstance(overlap, sgeom.Point):
                overlap = sgeom.LineString([overlap, overlap])
            if isinstance(overlap, sgeom.MultiLineString):
                for linestring in overlap:
                    segments_to_remove.append(linestring)
            else:
                segments_to_remove.append(overlap)

            #keep = True
            #for k, thinned_match in thinned_matches.iterrows():
                #if not overlap.intersection(thinned_match['geometry']).almost_equals(overlap):
                    #keep = False
            #if keep:
                #newlines.append(overlap)
        if segments_to_remove:
            if len(segments_to_remove)>1:
                toremove = sgeom.MultiLineString(segments_to_remove)
                # clean it up, segments might be out of order. line could double back.
                toremove = sops.linemerge(toremove.intersection(toremove))
            else:
                toremove = segments_to_remove[0]
            newlines = line['geometry'].difference(toremove)
            filtered_waterways.drop(i, inplace=True)
            if newlines.length > 0:
                if not isinstance(newlines, sgeom.MultiLineString):
                    newlines = [newlines]
                for newline in newlines:
                    line['geometry'] = sgeom.LineString(newline)
                    filtered_waterways = filtered_waterways.append(line)
    waterways = filtered_waterways[filtered_waterways.length > 0].reset_index()

    buff = env.get('buff', 100)*2
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

    with open(str(target[1]), 'wb') as fout:
        pickle.dump(affine, fout)
    with open(str(target[2]), 'wb') as fout:
        pickle.dump(imshape, fout)
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

    #holethresh = env.get('holethresh', 1000)
    #rivers = morph.remove_small_holes(rivers, min_size=holethresh, connectivity=2)

    rivers = morph.skeletonize(rivers)
    rivers = rivers.astype(np.uint8)

    # another closing to fix weird checkerboards and stuff
    rivers = morph.binary_closing(rivers, morph.square(3))
    #skeleton = morph.skeletonize(rivers)
    skeleton = morph.thin(rivers)
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
    #if n > minlen:
        #return [(j,i)]
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
            tokeep = []
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
                if (len(segment)<minlen) and (rivers[segment[-1]] == 1): # found terminating segment shorter than minlen
                    todelete.append(segment)
                else:
                    tokeep.append(segment)
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
                # dont delete segment if it will leave a circular loop
                # that is, if one branch will be deleted, but the other two branches directly connect
                if (len(todelete)==1 and (len(tokeep)==2) and (tokeep[0][-1]==tokeep[1][-1])):
                    tokeep.append(todelete.pop())
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
    with open(str(source[3]), 'rb') as fin:
        affine_riv = pickle.load(fin)
    with open(str(source[4]), 'rb') as fin:
        riv_shape = pickle.load(fin)
    with open(str(source[5]), 'r') as fin:
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
    maxx, miny = affine_riv * (riv_shape[1], riv_shape[0])

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
    nx.write_gpickle(G, str(target[0]))

    for node in Gclip.nodes():
        Gclip.node[node]['upstream'] = len(nx.ancestors(G, node)) # use numbers from full network
        Gclip.node[node]['downstream'] = len(nx.descendants(G, node)) # use numbers from full network
    nx.write_gpickle(Gclip, str(target[1]))

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


def convert_network_to_graphml(source, target, env):
    G = nx.read_gpickle(str(source[0]))
    nodes = list(G.nodes)
    for node in nodes:
        G.nodes[node]['basin'] = int(G.nodes[node]['basin'])
        G.nodes[node]['cellid'] = int(G.nodes[node]['cellid'])
        G.nodes[node]['downstream'] = int(G.nodes[node]['downstream'])
        G.nodes[node]['upstream'] = int(G.nodes[node]['upstream'])
        G.nodes[node]['lon'] = float(G.nodes[node]['ll'][0])
        G.nodes[node]['lat'] = float(G.nodes[node]['ll'][1])
        del G.nodes[node]['ll']
        G.nodes[node]['x'] = float(G.nodes[node]['xy'][0])
        G.nodes[node]['y'] = float(G.nodes[node]['xy'][1])
        del G.nodes[node]['xy']
        if 'branches' in G.nodes[node]:
            G.nodes[node]['branches'] = ' '.join([str(b) for b in G.nodes[node]['branches']])

    edges = list(G.edges)
    for edge in edges:
        if G.edges[edge]:
            branches = []
            weights = []
            for branch, weight in G.edges[edge]['branches'].items():
                branches.append(branch)
                weights.append(weight)
            branches = ' '.join([str(b) for b in branches])
            weights = ' '.join([str(w) for w in weights])
            G.edges[edge]['branches'] = branches
            G.edges[edge]['weights'] = weights

    nx.write_graphml(G, str(target[0]))
    return 0



def find_bifurs(source, target, env):
    with rasterio.open(str(source[0]), 'r') as rast:
        rivers = rast.read(1)
        meta = rast.meta.copy()

    # just in case this is run on a bifur_grid
    rivers[rivers != 0] = 1

    # add one pixel border to make sure any endpoints on edge are picked up
    bordered = np.zeros((rivers.shape[0]+2, rivers.shape[1]+2))
    bordered[1:-1,1:-1] = rivers
    rivers = bordered

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

    # remove one-pixel border
    bifurs = bifurs[1:-1, 1:-1]

    with rasterio.open(str(target[0]), 'w', **meta) as out:
        out.write(bifurs, 1)
    return 0


def plot_network_map(source, target, env):
    G = nx.read_gpickle(str(source[0]))
    with rasterio.open(str(source[1]), 'r') as rast:
        bifurs_pre = (rast.read(1) > 0).astype(np.int)
        affine = rast.transform
    with rasterio.open(str(source[2]), 'r') as rast:
        bifurs = (rast.read(1) > 0).astype(np.int)
    with open(str(source[3]), 'r') as fin:
        proj4_init = fin.read()
    labeltype = env['labels']

    proj4_params = {}
    for kv in proj4_init.split():
        if '=' in kv:
            k,v = kv.split('=')
            k = k.replace('+','')
            try:
                proj4_params[k] = float(v)
            except ValueError:
                proj4_params[k] = v

    assert proj4_params['proj'] == 'laea'
    laea = ccrs.LambertAzimuthalEqualArea(central_longitude = proj4_params['lon_0'],
                                          central_latitude = proj4_params['lat_0'])

    # reset upstream counts
    for node in G.nodes():
        G.node[node]['upstream'] = len(nx.ancestors(G, node))
        G.node[node]['downstream'] = len(nx.descendants(G, node))

    mpl.style.use('ggplot')
    fig, ax = plt.subplots(1,1, figsize=(10, 15), subplot_kw={'projection':laea})#, dpi=300)
    ax.coastlines('10m')

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
            with_labels=with_labels, labels=labels, font_size=8,
            alpha=.5, cmap=palettable.cartocolors.qualitative.Bold_10.mpl_colormap, ax=ax)
    nx.draw_networkx_edges(G, pos, arrowstyle='->', arrowsize=10, node_size=20, edge_color='k', ax=ax)
    if with_labels:
        nx.draw_networkx_labels(G, pos, labels=labels, font_size=8, ax=ax)
        for t in ax.texts:
            t.set_clip_on(False)
            t.set_rotation(30)

    I, J = np.meshgrid(np.arange(bifurs.shape[1]), np.arange(bifurs.shape[0]))
    xs, ys = affine * (I.flatten(), J.flatten())
    X = xs.reshape(I.shape)
    Y = ys.reshape(J.shape)

    norm = mpl.colors.Normalize(vmin=0, vmax=2)
    mpl.cm.Reds.set_bad(alpha=0)
    bifurs_mask = np.ma.masked_equal(bifurs, 0)
    ax.pcolormesh(X, Y, bifurs_mask, cmap=mpl.cm.Reds, norm=norm)

    mpl.cm.Greens.set_bad(alpha=0)
    old = bifurs_pre - bifurs
    old_mask = np.ma.masked_less_equal(old, 0)
    ax.pcolormesh(X, Y, old_mask, cmap=mpl.cm.Greens, norm=norm)

    mpl.cm.Blues.set_bad(alpha=0)
    new = bifurs - bifurs_pre
    new_mask = np.ma.masked_less_equal(new, 0)
    ax.pcolormesh(X, Y, new_mask, cmap=mpl.cm.Blues, norm=norm)

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
    G = nx.read_gpickle(str(source[0]))
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
    G = nx.read_gpickle(str(source[0]))
    with rasterio.open(str(source[1]), 'r') as rast:
        rivers = rast.read(1)
        affine = rast.transform
    coastline = geopandas.read_file(str(source[2])).loc[0,'geometry']
    coast = geopandas.read_file(str(source[3])).loc[0,'geometry']

    # calc node_dist_to_coast
    nodes = list(G.nodes())
    positions = geopandas.GeoSeries([sgeom.Point(G.node[node]['xy']) for node in nodes], index=nodes)
    node_dist_to_coast = positions.distance(coastline)

    # calc riv_dist_to_coast
    wet = np.where(rivers > 0)
    positions = geopandas.GeoSeries([sgeom.Point(affine * (i,j)) for j,i in zip(*wet)], index=zip(*wet))
    riv_dist_to_coast = positions.distance(coastline)
    riv_dist_to_coast[~positions.intersects(coast)] *= -1 # if river goes past coastline into ocean, make those distances negative so direction-finding still works

    node_dist_to_coast.to_pickle(str(target[0]))
    riv_dist_to_coast.to_pickle(str(target[1]))
    return 0

def calc_riv_flowdist_to_coast(source, target, env):
    with rasterio.open(str(source[0]), 'r') as rast:
        rivers = rast.read(1)
        affine = rast.transform
        meta = rast.meta.copy()
    with open(str(source[1]), 'rb') as fin:
        riv_dist_to_coast = pickle.load(fin)

    endpoints = np.where(rivers==1)
    outlets = [endpt for endpt in zip(*endpoints) if riv_dist_to_coast[endpt]<10000]

    mindist = np.zeros_like(rivers) * np.nan
    for outlet in outlets:
        tovisit = [(outlet, 0)]
        visited = set(outlet)
        while tovisit:
            rivpt, dist = tovisit.pop()
            mindist[rivpt] = np.nanmin([mindist[rivpt], dist])
            for dj in [-1,0,1]:
                for di in [-1,0,1]:
                    if dj==di==0:
                        continue
                    nextpt = (rivpt[0]+dj, rivpt[1]+di)
                    if (nextpt[0]<0 or nextpt[0]==rivers.shape[0] or
                            nextpt[1]<0 or nextpt[1]==rivers.shape[1]):
                        continue
                    if rivers[nextpt]>0 and nextpt not in visited:
                        tovisit.append((nextpt, dist+1))
            visited.add(rivpt)

    meta['dtype'] = mindist.dtype
    with rasterio.open(str(target[0]), 'w', **meta) as rast:
        rast.write(mindist, 1)
    return 0


def find_river_segments(source, target, env):
    with rasterio.open(str(source[0]), 'r') as rast:
        rivers = rast.read(1)
        affine = rast.transform

    endpoints = np.where(rivers == 1)
    # just start at the first one. all rivers are connected, so will get everywhere

    branchpoints = set()
    for j,i in zip(*np.where(rivers>=3)):
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
    nj, ni = rivers.shape
    while tovisit:
        j, i, cursegi = tovisit.pop(0)
        segments[cursegi].append((j,i))
        if (rivers[j,i] >= 3) or ((rivers[j,i] == 1) and (len(segments[cursegi])>1)):
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
                if j2 < 0 or j2 >= nj or i2 < 0 or i2 >= ni:
                    continue

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
                    if rivers[j,i] >= 3: # currently on bifur point
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
    #with open(str(source[3]), 'rb') as fin:
        #riv_dist_to_coast = pickle.load(fin)
    with rasterio.open(str(source[3])) as rast:
        flowdist = rast.read(1)

    lines = waterways['geometry']
    sindex = lines.sindex

    #for each segment
    counts = Counter()
    scores = Counter()
    directed_segments = {}
    x1, y1 = affine * (rivers.shape[1]//2, rivers.shape[0]//2)
    x2, y2 = affine * (rivers.shape[1]//2 + 1, rivers.shape[0]//2 + 1)
    dx = abs(x2-x1)
    dy = abs(y2-y1)
    pixelsize = np.sqrt(dx**2 + dy**2)
    pixel10 = pixelsize * 10
    used_dist_to_coast = []
    for segi in sorted(segments):
        segment = segments[segi]
        counts.clear()
        scores.clear()
        # for each point on segment
        for (j,i) in segment:
            x, y = affine * (i+.5, j+.5)
            # find nearest waterway
            pt = sgeom.Point(x, y)
            nearbylines_idx = list(sindex.intersection((x-pixel10, y-pixel10, x+pixel10, y+pixel10)))
            if nearbylines_idx:
                dists = [pt.distance(line) for line in lines.iloc[nearbylines_idx]]
                mindist = np.min(dists)
                ind = nearbylines_idx[np.argmin(dists)]
                scores[ind] += 1/max(mindist, pixelsize) # inverse distance weight so closest points count most, but dont go closer than nominal resolution since one very very close line could blow up comparison
                counts[ind] += (mindist < (2*pixelsize))
        # get direction of most common waterway
        ind, score = scores.most_common(1)[0]
        count = counts[ind]
        use_dist_to_coast = False
        if count >= min(7, len(segment)):
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
        if (count < min(7, len(segment))) or use_dist_to_coast:
            # use dist_to_coast since not enough segment pixels align with an osm waterway
            #if ((riv_dist_to_coast[segment[0]] < riv_dist_to_coast[segment[-1]])): # and
            if ((flowdist[segment[0]] < flowdist[segment[-1]])): # and
                    #(riv_dist_to_coast[segment[0]] < 10000) and
                    #(riv_dist_to_coast[segment[-1]] > 10000)): # replace with node resolution?
                directed_segments[segi] = segment[::-1]
                print('Segment {0}: {1} to {2} (reversed (b))'.format(segi, segment[-1], segment[0]))
            else:
                directed_segments[segi] = segment
                print('Segment {0}: {1} to {2} (no change (b))'.format(segi, segment[0], segment[-1]))
            used_dist_to_coast.append(segi)

    # check for bifurs which are only sinks or sources. suggests an incorrect direction. change one of the dist_to_coast method segments if there is one. otherwise change shortest segment. repeat until no changes
    bifurs = np.where(rivers >= 3)
    while True:
        reversed_segment = False
        for j,i in zip(*bifurs):
            sinks = []
            sources = []
            for segi in list(directed_segments.keys()):
                segment = directed_segments[segi]
                if (j,i) == segment[0]:
                    sources.append(segi)
                if (j,i) == segment[-1]:
                    sinks.append(segi)
            if (not sources) or (not sinks):
                segis = sources + sinks
                # water can't get to or from this bifur
                candidates = [segi for segi in segis if segi in used_dist_to_coast]
                if len(candidates) == 1:
                    segi = candidates[0]
                    print('Found source-only or sink-only bifur - reversing segment {}'.format(segi))
                    directed_segments[segi] = directed_segments[segi][::-1]
                else:
                    lengths = [len(directed_segments[segi]) for segi in segis]
                    segi = segis[np.argmin(lengths)]
                    print('Found source-only or sink-only bifur - reversing shortest segment {}'.format(candidates))
                    directed_segments[segi] = directed_segments[segi][::-1]
                reversed_segment = True
        if not reversed_segment:
            break

    # check for pairs of segments with same start/end points, but flowing in different
    # directions. should be in same direction. change to ensure no source-only or sink-only
    # endpoints
    while True:
        reversed_segment = False
        for segi, segj in itertools.combinations(directed_segments.keys(), 2):
            seg1 = directed_segments[segi]
            seg2 = directed_segments[segj]
            if (seg1[0] == seg2[-1] and seg1[-1] == seg2[0]):
                for segk, segl in itertools.combinations(directed_segments.keys(), 2):
                    if (segk in [segi, segj]) or (segl in [segi, segj]):
                        continue
                    seg3 = directed_segments[segk]
                    seg4 = directed_segments[segl]
                    if seg3[-1] == seg1[0] and seg4[0] == seg1[-1]: #seg3 upstream, seg4 down
                        # change seg2
                        directed_segments[segj] = directed_segments[segj][::-1]
                        print('Mismatched loop, reversing segment {}'.format(segi))
                        reversed_segment = True
                    elif seg3[-1] == seg2[0] and seg4[0] == seg2[-1]: #seg3 upstream, seg4 down
                        # change seg1
                        directed_segments[segi] = directed_segments[segi][::-1]
                        print('Mismatched loop, reversing segment {}'.format(segi))
                        reversed_segment = True
        if not reversed_segment:
            break

    with open(str(target[0]), 'wb') as fout:
        pickle.dump(directed_segments, fout)
    return 0


def set_segment_widths(source, target, env):
    with open(str(source[0]), 'rb') as fin:
        segments = pickle.load(fin)
    with rasterio.open(str(source[1])) as rast:
        rivers = rast.read(1)
        affine = rast.transform
        meta = rast.meta.copy()
    widths = geopandas.read_file(str(source[2]))
    widths_rast = np.zeros_like(rivers) * np.nan
    lines = widths['geometry']
    sindex = lines.sindex

    segment_widths = {}
    scores = Counter()
    x1, y1 = affine * (rivers.shape[1]//2, rivers.shape[0]//2)
    x2, y2 = affine * (rivers.shape[1]//2 + 1, rivers.shape[0]//2 + 1)
    dx = abs(x2-x1)
    dy = abs(y2-y1)
    pixelsize = np.sqrt(dx**2 + dy**2)
    pixel10 = pixelsize * 10
    for segi in sorted(segments):
        segment = segments[segi]
        # find nearest GRWL river. each segment gets a starting width and ending width as average along some distance of nearest GRWL river segments
        scores.clear()
        nearby_widths = defaultdict(list) # keep flowdir-ordered list of widths for each line that is close enough to river
        # for each point on segment
        # segment is directed, so widths are added to list in order
        for (j,i) in segment:
            x, y = affine * (i+.5, j+.5)
            # find nearest waterway
            pt = sgeom.Point(x, y)
            nearbylines_idx = list(sindex.intersection((x-pixel10, y-pixel10, x+pixel10, y+pixel10)))
            #dists = [pt.distance(line) for line in lines]
            if nearbylines_idx:
                dists = [pt.distance(line) for line in lines.iloc[nearbylines_idx]]
                mindist = np.min(dists)
                #ind = np.argmin(dists)
                ind = nearbylines_idx[np.argmin(dists)]
                GRWL_segmentID = widths.iloc[ind]['segmentID']
                scores[GRWL_segmentID] += 1/max(mindist, pixelsize) # inverse distance weight so closest points count most, but dont go closer than nominal resolution since one very very close line could blow up comparison
                if (mindist < 2*pixelsize):
                    nearby_widths[GRWL_segmentID].append(widths.iloc[ind]['width_m'])

        # get widths of most common waterway
        ws = []
        if nearby_widths:
            GRWL_segmentID, score = scores.most_common(1)[0]
            ws = [w for w in nearby_widths[GRWL_segmentID] if w>1] # GRWL seems to use width==1 as missing data placeholder
        if len(ws) <= 7:
            n = len(ws)
        elif len(ws) <= 14:
            n = 7
        else:
            n = len(ws)//2
        width_start = np.mean(ws[:n]) if n>0 else None
        width_end = np.mean(ws[-n:]) if n>0 else None
        segment_widths[segi] = (width_start, width_end)
        print('Segment {0}: {1}, {2}'.format(segi, GRWL_segmentID, segment_widths[segi]))

    # give each point on each segment an interpolated width. will make branching easier
    river_widths = defaultdict(list)
    for segi, segment in segments.items():
        n = len(segment)
        widths = segment_widths[segi]
        fracs = np.linspace(1, 0, n)
        for i, (rivj,rivi) in enumerate(segment):
            if widths == (None, None):
                river_widths[rivj,rivi].append(None)
            else:
                frac = fracs[i]
                river_widths[rivj,rivi].append(frac*widths[0] + (1-frac)*widths[1])
    # bifur points will have multiple widths.
    for (rivj,rivi), widths in list(river_widths.items()):
        widths = [w for w in widths if w is not None]
        if widths:
            river_widths[rivj,rivi] = np.mean(widths)
            widths_rast[rivj,rivi] = np.mean(widths)
        else:
            river_widths[rivj,rivi] = None

    with open(str(target[0]), 'wb') as fout:
        pickle.dump(river_widths, fout)
    meta['dtype'] = widths_rast.dtype
    with rasterio.open(str(target[1]), 'w', **meta) as rast:
        rast.write(widths_rast, 1)
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
        nupstream = pickle.load(fin)
    with open(str(source[2]), 'rb') as fin:
        ndownstream = pickle.load(fin)
    with open(str(source[3]), 'rb') as fin:
        positions = pickle.load(fin)

    nodes = sorted(nupstream) # sorting keys of dict
    nupstream = [nupstream[node] for node in nodes]
    ndownstream = [ndownstream[node] for node in nodes]
    positions = [positions[node] for node in nodes]

    def _extend_river(rivers, mindist_rivpt, head_node_ij):
        # draw new river in straight line from closest river approach to head_node
        head_node_ji = tuple(head_node_ij[::-1])
        rowidx, colidx = skimage.draw.line(*mindist_rivpt, *head_node_ji)
        rivers[rowidx, colidx] = 1
        return rivers

    score = np.array(nupstream) * np.array(ndownstream) # maximized at "center" of mainstem, which we assume is upstream of the delta boundary
    head_node_i = np.argmax(score)
    head_node = nodes[head_node_i]

    wet = np.where(rivers)
    head_node_xy = positions[head_node_i]
    head_node_x, head_node_y = positions[head_node_i]
    head_node_ij = [int(val) for val in ~affine * (head_node_x, head_node_y)]
    mindist = np.inf
    for rivpt in zip(*wet):
        rivj, rivi = rivpt
        rivx, rivy = affine * (rivi,rivj)
        dist = np.sqrt((head_node_x - rivx)**2 + (head_node_y - rivy)**2)
        if dist < mindist:
            mindist = dist
            mindist_rivpt = rivpt
    rivers = _extend_river(rivers, mindist_rivpt, head_node_ij)

    with rasterio.open(str(target[0]), 'w', **meta) as rast:
        rast.write(rivers, 1)
    return 0


def extend_rivers_to_coast(source, target, env):
    with rasterio.open(str(source[0]), 'r') as rast:
        bifurs = rast.read(1)
        affine = rast.transform
        meta = rast.meta.copy()
    with open(str(source[1]), 'rb') as fin:
        next_rivpt = pickle.load(fin)
    with open(str(source[2]), 'rb') as fin:
        prev_rivpt = pickle.load(fin)
    with open(str(source[3]), 'rb') as fin:
        rivdisttocoast = pickle.load(fin)
    coastline = geopandas.read_file(str(source[4]))['geometry']
    coast = geopandas.read_file(str(source[5]))

    polys = coast.explode()
    ind = polys.area.values.argmax()
    coast = polys.iloc[ind:ind+1]

    lines = coastline.explode()
    ind = lines.length.values.argmax()
    coastline = lines.iloc[ind:ind+1]

    coastbuff = coastline.buffer(20000)
    coastbuff_boundary = coastbuff.boundary
    inside = coastbuff_boundary.intersection(coast)
    outside = coastbuff_boundary.difference(inside)

    validzone = np.zeros_like(bifurs) * np.nan
    rfeatures.rasterize(((geom, 1) for geom in coastbuff), out_shape=validzone.shape,
            out=validzone, transform=affine, all_touched=True)

    maxdepth = 20
    depths = np.zeros_like(bifurs) * np.nan
    rfeatures.rasterize(((geom, -maxdepth) for geom in inside), out_shape=depths.shape,
            out=depths, transform=affine, all_touched=True)
    rfeatures.rasterize(((geom, 0) for geom in coastline), out_shape=depths.shape,
            out=depths, transform=affine, all_touched=True)
    rfeatures.rasterize(((geom, maxdepth) for geom in outside), out_shape=depths.shape,
            out=depths, transform=affine, all_touched=True)

    xs = np.arange(depths.shape[1])
    ys = np.arange(depths.shape[0])
    X, Y = np.meshgrid(xs, ys)
    valid = np.isfinite(depths)
    depthsgood = depths[valid]
    Xgood = X[valid]
    Ygood = Y[valid]
    depths = griddata(np.c_[Xgood, Ygood], depthsgood, (X, Y), method='linear')
    depths[np.isnan(depths)] = maxdepth
    depths = gaussian_filter(depths, 5, mode='nearest')
    depths[np.isnan(validzone)] = np.nan

    grad_y, grad_x = np.gradient(depths)
    grad_y[np.isclose(grad_y, 0)] = 0
    grad_x[np.isclose(grad_x, 0)] = 0
    angle = np.arctan2(grad_y, grad_x)
    angle[np.isnan(depths)] = np.nan

    endpoints = np.where(bifurs==1)
    outlets = [p for p in zip(*endpoints) if ((not next_rivpt[p]) and (rivdisttocoast[p] < 20000))]
    rivers = (bifurs > 0).astype(bifurs.dtype)

    for rivj, rivi in outlets:
        # in case no gradient at start, start with initial dj,di leading away from prev
        prevjs = [prev_rivpt[rivj, rivi][0][0]]
        previs = [prev_rivpt[rivj, rivi][0][1]]
        dj = rivj - prevjs[-1]
        di = rivi - previs[-1]
        rivj_i = rivj
        rivi_i = rivi
        rivj, rivi = rivj+.5, rivi+.5
        nj, ni = angle.shape
        while (0<=rivj_i<nj and
               0<=rivi_i<ni and
               np.isfinite(angle[rivj_i,rivi_i])):
            # check all neighbors of potential new point
            # if any (other than this river's prevs) is already a river, break
            # dont want to connect to any other branchs
            collision = False
            for dj in [-1,0,1]:
                for di in [-1,0,1]:
                    if dj == di == 0:
                        continue
                    rivj2 = rivj_i+dj
                    rivi2 = rivi_i+di
                    if rivj2 in prevjs and rivi2 in previs:
                        continue
                    if rivj2 < 0 or rivj2 >= nj or rivi2 < 0 or rivi2 >= ni:
                        continue
                    if rivers[rivj2, rivi2]:
                        collision = True
            if collision:
                break
            rivers[rivj_i,rivi_i] = 1
            rivxy = affine * (rivi, rivj)
            a = angle[rivj_i,rivi_i]
            if (grad_y[rivj_i, rivi_i]!= 0) or (grad_x[rivj_i, rivi_i]!= 0):
                # only update if gradient is not zero
                # otherwise just continue in last direction
                dj = np.sin(a)
                di = np.cos(a)
            prevjs.append(rivj_i)
            previs.append(rivi_i)
            rivj += dj
            rivi += di
            rivj_i = int(round(rivj))
            rivi_i = int(round(rivi))

    # clean up any multipixel corners, turns
    rivers = (rivers>0)
    rivers = morph.skeletonize(rivers).astype(np.uint8)

    with rasterio.open(str(target[0]), 'w', **meta) as out:
        out.write(rivers, 1)
    return 0


def remap_riv_network(source, target, env):
    G = nx.read_gpickle(str(source[0]))
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
        riv_dist_to_coast = pickle.load(fin)
    with open(str(source[7]), 'rb') as fin:
        river_widths = pickle.load(fin)

    nodes = [node for node in G.nodes()]
    positions = [G.node[node]['xy'] for node in nodes]
    nupstream = [G.node[node]['upstream'] for node in nodes]
    ndownstream = [G.node[node]['downstream'] for node in nodes]
    cellid = [G.node[node]['cellid'] for node in nodes]

    maxresolution = -np.inf
    for node1 in nodes:
        x1, y1 = positions[nodes.index(node1)]
        # check both diagonals
        node2 = (node1[0]+1, node1[1]+1)
        if node2 in nodes:
            x2, y2 = positions[nodes.index(node2)]
            dist = np.sqrt((x2-x1)**2 + (y2-y1)**2)
            if dist > maxresolution:
                maxresolution = dist
        node3 = (node1[0]+1, node1[1]-1)
        if node3 in nodes:
            x3, y3 = positions[nodes.index(node3)]
            dist = np.sqrt((x3-x1)**2 + (y3-y1)**2)
            if dist > maxresolution:
                maxresolution = dist

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
        # dont start on a node, wait to snap to one. should be first step since river was extended to hit the head node
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
        xy2 = affine * (rivi+1,rivj)
        xy3 = affine * (rivi,rivj+1)
        pixelres = max(abs(xy2[0]-xy[0]), abs(xy3[1]-xy[1]))
        if last_node_i is not None:
            last_node = nodes[last_node_i]
            last_cell = cellid[last_node_i]
            visited.add(((rivj, rivi), last_node_i)) # include last_node_i so that different branches coming from different nodes can re-visit a node. but branches that are on the same node and same rivpts are dropped, since they will trace the same route
            if 'branches' not in G.nodes[last_node]:
                # first time the river snaps to a node, that node has not yet been added to a branch. do that
                G.nodes[last_node]['branches'] = {branch}
        else:
            last_node = None
            last_cell = None

        next_node = nearestnode[(rivj,rivi)]
        next_node_i = nodes.index(next_node)
        next_cell = cellid[next_node_i]
        if last_node is not None:
            isdiag = int(abs(next_node[0]-last_node[0])==1 and abs(next_node[1]-last_node[1])==1)
        else:
            isdiag = 0

        xriv, yriv = xy
        xnode, ynode = positions[next_node_i]
        xp = xriv - xnode
        yp = yriv - ynode
        dist = np.sqrt(xp**2 + yp**2)
        #halfres = .5 * maxresolution
        widehalfres = .55 * maxresolution
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
        #riv_near_node = ((dist <= halfres) and
                         #((yp + 2*xp <  widehalfres) or (yp + xp/2 <  widehalfres/2)) and
        # squished star
        #riv_near_node = (((yp + 2*xp <  widehalfres) or (yp + xp/2 <  widehalfres/2)) and
                         #((yp - 2*xp > -widehalfres) or (yp - xp/2 > -widehalfres/2)) and
                         #((yp + 2*xp > -widehalfres) or (yp + xp/2 > -widehalfres/2)) and
                         #((yp - 2*xp <  widehalfres) or (yp - xp/2 <  widehalfres/2)))
        # small circles with legs along boundary, with expanded diagonal circles
        #riv_near_node = ((dist < ((.3 * maxresolution) * (1 + isdiag * (np.sqrt(2) - 1)))) or # within small circle around next_node (larger for diagonal nextnode)
                         #((abs(yp) < (2 * pixelres)) and (abs(xp) < (.5 * maxresolution))) or # within 2 river pixels of nodeline
                         #((abs(xp) < (2 * pixelres)) and (abs(yp) < (.5 * maxresolution))))   # within 2 river pixels of nodeline
        # small circles with legs along node-cell boundary
        riv_near_node = ((dist < (.2 * maxresolution)) or # within small circle around next_node (larger for diagonal nextnode)
                         ((abs(yp) < (2 * pixelres)) and (abs(xp) < (.5 * maxresolution))) or # within 2 river pixels of nodeline
                         ((abs(xp) < (2 * pixelres)) and (abs(yp) < (.5 * maxresolution))))   # within 2 river pixels of nodeline
        # just thin legs along boundary
        #riv_near_node = (((abs(yp) < (2 * pixelres)) and (abs(xp) < (.5 * maxresolution))) or # within 2 river pixels of nodeline
                         #((abs(xp) < (2 * pixelres)) and (abs(yp) < (.5 * maxresolution))))    # within 2 river pixels of nodeline
        valid_next_node = (riv_near_node and # helps reduce zig-zag
                           ((last_node is None) or
                               abs(next_node[0]-last_node[0]) <= 1) and # dont jump past neighbor cell
                           ((last_node is None) or
                               abs(next_node[1]-last_node[1]) <= 1))   #  (possible along coastline)


        if ((next_node != last_node) and # moving to new node
                valid_next_node):

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
                    print('New:', branch, ': cellid {0} to {1}'.format(last_cell, next_cell))
                    G.add_edge(last_node, next_node)
                    edits[last_cell][next_cell][branch] = True
                    if 'branches' not in G.nodes[next_node]:
                        G.nodes[next_node]['branches'] = {branch}
                    else:
                        G.nodes[next_node]['branches'].add(branch)
                    if 'branches' not in G.edges[last_node, next_node]: # can we ever even have multiple???
                        G.edges[last_node, next_node]['branches'] = {branch: river_widths[rivj,rivi]}
                    else:
                        G.edges[last_node, next_node]['branches'][branch] = river_widths[rivj,rivi]
                    for downstream_node in list(G.successors(last_node)):
                        # remove existing downstream links (other than this or other new ones)
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
                            print('Swap overlap: from {} to {}'.format(corner2_cell, corner1_cell))
                            G.remove_edge(corner2_node, corner1_node)
                            edits[corner2_cell][corner1_cell][branch] = None
                            G.add_edge(corner2_node, next_node)
                            if 'branches' not in G.edges[corner2_node, next_node]: # can we ever even have multiple???
                                G.edges[corner2_node, next_node]['branches'] = {branch: river_widths[rivj,rivi]}
                            else:
                                G.edges[corner2_node, next_node]['branches'][branch] = river_widths[rivj,rivi]
                            edits[corner2_cell][next_cell][branch] = True # anything other than None
                        if (corner2_node in G.succ[corner1_node]):
                            print('Swap overlap: from {} to {}'.format(corner1_cell, corner2_cell))
                            G.remove_edge(corner1_node, corner2_node)
                            edits[corner1_cell][corner2_cell][branch] = None
                            G.add_edge(corner1_node, next_node)
                            if 'branches' not in G.edges[corner1_node, next_node]: # can we ever even have multiple???
                                G.edges[corner1_node, next_node]['branches'] = {branch: river_widths[rivj,rivi]}
                            else:
                                G.edges[corner1_node, next_node]['branches'][branch] = river_widths[rivj,rivi]
                            edits[corner1_cell][next_cell][branch] = True # anything other than None

        elif (not valid_next_node): #or
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
            (riv_dist_to_coast[rivj,rivi] <= 10000) and # meters
            (next_node is not None)):
            # no downstream points, AND CLOSE TO COAST, remove downstream flow from node
            # dont do this for likely upstream points, just leave existing connections
            # helps with errors in osm_river, dont want to strand water upstream
            # next_node is next if we just moved to new one, or last_node if we didn't
            outlets.add(next_cell)
            for node2 in list(G.successors(next_node)):
                node2_cell = cellid[nodes.index(node2)]
                print('Removing original outflows from outlet cell {} to {}'.format(next_cell, node2_cell))
                G.remove_edge(next_node, node2)
                edits[next_cell][node2_cell][branch] = None

    with open(str(target[0]), 'w', newline='') as fout:
        csvwriter = csv.writer(fout)
        wrote = []
        fracs = {}
        for oldedge in Gorig.edges:
            # zero out removed links
            if oldedge not in G.edges:
                from_cell, to_cell = [cellid[nodes.index(n)] for n in oldedge]
                csvwriter.writerow([from_cell, to_cell, 0])
                wrote.append((from_cell, to_cell, 0))
        for newedge in G.edges:
            if newedge not in Gorig.edges:
                # make new link
                from_node, to_node = newedge
                from_cell, to_cell = [cellid[nodes.index(n)] for n in newedge]
                discharges = {}
                total_discharge = 0
                for neighbor_node in G.succ[from_node]:
                    discharge = 0
                    for branch, width in G.edges[(from_node, neighbor_node)]['branches'].items():
                        # power law relationship between river width and discharge
                        # kellerhals and church 1989, bray 1973,1975,
                        # also http://publications.gc.ca/collections/collection_2007/dfo-mpo/Fs97-6-2637E.pdf page 10
                        # W ~ a * Q**.5
                        # a ~ 4-4.5 for bankfull flow, ~9-10 for mean flow
                        # exponent ~.4 for very small rivers, ~.55 for very large, ~.5 for medium
                        # constant a cancels out in these ratios. use exponent=.5 for all (Q ~ w**2)
                        if width is not None:
                            discharge += width**2
                        else:
                            discharge += 0
                    total_discharge += discharge
                    if discharge == 0:
                        discharge = None
                    discharges[neighbor_node] = discharge
                if total_discharge == 0:
                    # all widths are none. just set each value to a constant, will end up with balanced branches
                    for neighbor_node in G.succ[from_node]:
                        discharges[neighbor_node] = 1
                else:
                    # check for any missing widths (discharges) set each missing to 10% of total we do have
                    for neighbor_node in G.succ[from_node]:
                        # adjust neighbor links to account for fractional flow
                        if discharges[neighbor_node] is None:
                            discharges[neighbor_node] = .1 * total_discharge
                total_discharge = sum(list(discharges.values()))
                dis_fracs = [np.round(discharges[dnode]/total_discharge * 1000)/1000 for dnode in G.succ[from_node]]
                if sum(dis_fracs) != 1.0:
                    print('Discharge fractions from cell {0}: {1} (sum: {2})'.format(from_cell, dis_fracs, sum(dis_fracs)))
                    ind = np.argmin(dis_fracs)
                    smallest = dis_fracs.pop(ind)
                    replacement = np.round((1 - sum(dis_fracs)) * 1000) / 1000
                    dis_fracs.insert(ind, replacement)
                    print(' Adjusted frac {0} to {1}'.format(smallest, replacement))
                for neighbor_node, dis_frac in zip(G.succ[from_node], dis_fracs):
                    #dis_frac = discharges[neighbor_node] / total_discharge
                    from_cell = cellid[nodes.index(from_node)]
                    neighbor_cell = cellid[nodes.index(neighbor_node)]
                    if (from_node, neighbor_node) in Gorig.edges:
                        # first zero out a link before adjusting for fractional flow
                        if (from_cell, neighbor_cell, 0) not in wrote:
                            csvwriter.writerow([from_cell, neighbor_cell, 0])
                            wrote.append((from_cell, neighbor_cell, 0))
                    if (from_cell, neighbor_cell, 1) not in wrote: # use 1 as sentinel to avoid round-off errors. will revisit some edges if both bifur branches are new. neighbor_cell will be computed twice, but dont want to write twice
                        csvwriter.writerow([from_cell, neighbor_cell, dis_frac])
                        wrote.append((from_cell, neighbor_cell, 1))
                        fracs[(from_cell, neighbor_cell)] = dis_frac

    # starting at head_rivpt (which cooresponds to upstream discharge source), walk down and
    # track fraction of original flow at each node. for estimating outlet fraction
    tovisit = [(nearestnode[head_rivpt], 1)]
    while tovisit:
        node, flow = tovisit.pop(0)
        if 'flow' in G.node[node]:
            G.node[node]['flow'] += flow
        else:
            G.node[node]['flow'] = flow
        for to_node in G.succ[node]:
            from_cell = cellid[nodes.index(node)]
            to_cell = cellid[nodes.index(to_node)]
            if (from_cell, to_cell) in fracs:
                flow_frac = flow * fracs[(from_cell, to_cell)]
            else:
                flow_frac = flow
            tovisit.append((to_node, flow_frac))

    for node in G.nodes():
        G.node[node]['upstream'] = len(nx.ancestors(G, node))
        G.node[node]['downstream'] = len(nx.descendants(G, node))
    nx.write_gpickle(G, str(target[1]))

    with open(str(target[2]), 'w') as fout:
        for outlet in sorted(outlets):
            node = nodes[cellid.index(outlet)]
            fout.write('{0},{1:0.3f}\n'.format(outlet, G.node[node]['flow']))

    return 0


def convert_bifur_cellids_to_SSEA(source, target, env):
    stn = pandas.read_csv(str(source[0]), sep='\t', header=0)[['CellXCoord', 'CellYCoord', 'CellID']]
    ssea = pandas.read_csv(str(source[1]), sep='\t', header=0)[['CellXCoord', 'CellYCoord', 'CellID']]
    stnbifurs = pandas.read_csv(str(source[2]), header=None, names=['FromCell_stn', 'ToCell_stn', 'frac'])

    cells = stn.set_index(['CellXCoord', 'CellYCoord']).join(ssea.set_index(['CellXCoord', 'CellYCoord']), how='left', lsuffix='_stn', rsuffix='_ssea')
    stnbifurs = stnbifurs.join(cells.set_index('CellID_stn'), on='FromCell_stn', how='left')
    stnbifurs = stnbifurs.join(cells.set_index('CellID_stn'), on='ToCell_stn', how='left', lsuffix='_from', rsuffix='_to')

    sseabifurs = stnbifurs[['CellID_ssea_from', 'CellID_ssea_to', 'frac']]
    sseabifurs.to_csv(str(target[0]),index=False,header=False)
    return 0
