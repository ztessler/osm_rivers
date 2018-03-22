import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import geopandas
import rasterio
import rasterio.features as rfeatures
import cartopy.crs as ccrs
import networkx as nx
import skimage.morphology as morph
from netCDF4 import Dataset


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


def thin_vec(source, target, env):
    rivers = geopandas.read_file(str(source[0]))
    thresh = env.get('thresh', 100)

    rivers_eroded = rivers.buffer(-thresh)
    rivers_clean = rivers_eroded[rivers_eroded.area>0].buffer(thresh)
    rivers_merge = geopandas.overlay(geopandas.GeoDataFrame(rivers_clean, columns=['geometry']), geopandas.GeoDataFrame(rivers_clean, columns=['geometry']), how='union')

    rivers_merge.to_file(str(target[0]))
    return 0


def rasterize_riv(source, target, env):
    rivers = geopandas.read_file(str(source[0]))
    imshape = env.get('imshape', (1000, 1000))

    bounds = rivers.bounds
    width = bounds['maxx'].max() - bounds['minx'].min()
    height = bounds['maxy'].max() - bounds['miny'].min()
    dx = width / imshape[0]
    dy = height / imshape[1]

    a=dx; b=0; c=bounds['minx'].min(); d=0; e=-dy; f=bounds['maxy'].max()
    affine = rasterio.Affine(a,b,c,d,e,f)

    burned = rfeatures.rasterize(
            ((f,255) for f in rivers['geometry']),
            out_shape=imshape,
            transform=affine,
            all_touched=True,
            dtype=np.uint8)

    with rasterio.open(str(target[0]), 'w',
            driver='GTiff',
            width=imshape[0], height=imshape[1],
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
    rivers = morph.closing(rivers, morph.square(closing))

    holethresh = env.get('holethresh', 1000)
    rivers = morph.remove_small_holes(rivers, min_size=holethresh, connectivity=2)

    skeleton, distance = morph.medial_axis(rivers, return_distance=True)
    skeleton = skeleton.astype(np.uint8)
    skeleton[skeleton>0] = 255

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


def import_rgis_network(source, target, env):
    with rasterio.open(str(source[0]), 'r') as rast:
        basins = rast.read(1)
        nodata = rast.nodata
    with rasterio.open(str(source[1]), 'r') as rast:
        flowdir = rast.read(1)

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
    G = nx.DiGraph()
    for y in range(basins.shape[0]):
        for x in range(basins.shape[1]):
            if basins[y,x] != nodata:
                G.add_node((x,y))
                tocell = flowdir[y,x]
                if tocell != nodata:
                    dx, dy = neighbors[tocell]
                    if basins[y+dy, x+dx] != nodata:
                        G.add_edge((x, y), (x+dx, y+dy))
    nx.write_yaml(G, str(target[0]))
    return 0


def plot_vec_rivs(source, target, env):
    rivers = geopandas.read_file(str(source[0]))

    mpl.style.use('ggplot')

    fig, ax = plt.subplots(1, 1, figsize=(8,6))
    rivers.plot(ax=ax)
    ax.set_aspect('equal')
    fig.savefig(str(target[0]))

    return 0



