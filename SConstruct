# vim: fileencoding=UTF-8
# vim: filetype=python

### TO RUN:

import os
import sys
import lib

SetOption('max_drift', 1)

env = Environment(ENV = {'PATH' : os.environ['PATH'],
                         'GDAL_DATA': os.environ['GDAL_DATA'],
                         'GHAASBIN': os.environ.get('GHAASDIR', '/Users/ecr/ztessler/opt/ghaas_bifur/ghaas/bin/'),
                         })
env.Decider('MD5-timestamp')

GHAASBIN = env['ENV']['GHAASBIN']
work = 'work'
output = 'output'
figures = 'figures'

delta = os.environ.get('DELTA', 'Mekong')
OSMriver = os.environ.get('OSMriver', 'vietnam')
STNres = os.environ.get('STNres', '06min')

STNnetwork = '/Users/ecr/ztessler/projects/CHART/WBM/tools/buildNetwork/output/{delta}_Network_{res}.gdbn'.format(delta=delta, res=STNres)
OSMrivers = '/Users/ecr/ztessler/projects/CHART/WBM/tools/osm_rivers/osm_data/{0}/gis.osm_water_a_free_1.shp'.format(OSMriver)
deltashp = '/Users/ecr/ztessler/data/deltas_LCLUC/maps/{0}_shp/{0}.shp'.format(delta)

work = 'work'
output = 'output'

# project and clip river vectors to delta
clipped_vec = os.path.join(work, '{0}_riv_clipped/{0}_riv_clipped.shp'.format(delta))
proj4str = os.path.join(work, '{}_proj4.txt'.format(delta))
env.Command(
        source=[OSMrivers, deltashp],
        target=[clipped_vec, proj4str],
        action=lib.project_and_clip_osm_rivers)

thinned_vec = os.path.join(work, '{0}_riv_thinned/{0}_riv_thinned.shp'.format(delta))
env.Command(
        source=clipped_vec,
        target=thinned_vec,
        action=lib.thin_vec,
        thresh=100)

# plot
env.Command(
        source=thinned_vec,
        target=os.path.join(figures, '{}_thinned_rivs.png'.format(delta)),
        action=lib.plot_vec_rivs,
        imshape=(1000, 1000))


# rasterize
riv_rast = os.path.join(work, '{0}_riv_rast.tif'.format(delta))
env.Command(
        source=thinned_vec,
        target=riv_rast,
        action=lib.rasterize_riv)

# skeletonize raster
riv_skel = os.path.join(work, '{0}_riv_skeleton.tif'.format(delta))
env.Command(
        source=riv_rast,
        target=riv_skel,
        action=lib.skeleton_riv,
        closing=5,
        holethresh=1000)

# drop small rivers
riv_clean = os.path.join(work, '{0}_riv_cleaned.tif'.format(delta))
env.Command(
        source=riv_skel,
        target=riv_clean,
        action=lib.keep_n_rivers,
        n=1)

# add record id column to network cell table
network = os.path.join(work, '{0}_{1}_network.gdbn'.format(delta, STNres))
env.Command(
        source=STNnetwork,
        target=network,
        action=GHAASBIN+'tblAddIdXY $SOURCE $TARGET')

# import RGIS network
cellid = os.path.join(work, '{0}_{1}_cellid.{{ext}}'.format(delta, STNres))
basins = os.path.join(work, '{0}_{1}_basins.{{ext}}'.format(delta, STNres))
flowdir = os.path.join(work, '{0}_{1}_flowdir.{{ext}}'.format(delta, STNres))
for (path, varname) in [(cellid, 'CellID'),
                        (basins, 'BasinID'),
                        (flowdir, 'ToCell')]:
    env.Command(
            source=network,
            target=path.format(ext='nc'),
            action=[
                'netCells2Grid -f {0} -t {0} -u {0} -d {1} $SOURCE ${{TARGET}}.1'.format(varname, delta),
                'grdRenameLayers -r 1 XXXX ${TARGET}.1 ${TARGET}.2',
                'grdDateLayers -y 1 -e day ${TARGET}.2 ${TARGET}.3',
                'rgis2netcdf ${TARGET}.3 $TARGET'])
    env.Command(
            source=path.format(ext='nc'),
            target=path.format(ext='tif'),
            action=lib.georef_nc)

network = os.path.join(work, '{0}_{1}_network.nx.yaml'.format(delta, STNres))
env.Command(
        source=[cellid.format(ext='tif'),
                basins.format(ext='tif'),
                flowdir.format(ext='tif')],
        target=network,
        action=lib.import_rgis_network)

