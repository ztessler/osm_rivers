# vim: fileencoding=UTF-8
# vim: filetype=python

### TO RUN:

import os
import sys
import hashlib

import lib

SetOption('max_drift', 1)

env = Environment(ENV = {'PATH' : os.environ['PATH'],
                         'GDAL_DATA': os.environ['GDAL_DATA'],
                         'GHAASBIN': os.environ.get('GHAASDIR', '/Users/ecr/ztessler/opt/ghaas_bifur/ghaas/bin/'),
                         })
env.Decider('MD5-timestamp')

def myCommand(target, source, action, **kwargs):
    '''
    env.Command wrapper that forces env override arguments to be sconsign
    signature database. Wraps all extra kwargs in env.Value nodes and adds
    them to the source list, after the existing sources. Changing the extra
    arguments will cause the target to be rebuilt, as long as the data's string
    representation changes.
    '''
    def hash(v):
        # if this is changed then all targets with env overrides will be rebuilt
        return hashlib.md5(repr(v).encode('utf-8')).hexdigest()
    if not isinstance(source, list):
        source = [source]
    if None in source:
        source.remove(None)
    kwargs['nsources'] = len(source)
    source.extend([env.Value('{}={}'.format(k,hash(v))) for k,v in kwargs.items()])
    return env.Command(target=target, source=source, action=action, **kwargs)

GHAASBIN = env['ENV']['GHAASBIN']
work = 'work'
output = 'output'
figures = 'figures'

delta = os.environ.get('DELTA', 'Mekong')
OSMriver = os.environ.get('OSMriver', 'vietnam')
STNres = os.environ.get('STNres', '06min')

STNnetwork = '/Users/ecr/ztessler/projects/CHART/WBM/tools/buildNetwork/output/{delta}/{delta}_Network_{res}.gdbn'.format(delta=delta, res=STNres)
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
#env.Command(
        #source=[OSMrivers, deltashp],
        #target=clipped_vec,
        #action=lib.clip_osm_rivers)

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
        action=lib.plot_vec_rivs)


# rasterize
riv_rast = os.path.join(work, '{0}_riv_rast.tif'.format(delta))
env.Command(
        source=thinned_vec,
        target=riv_rast,
        action=lib.rasterize_riv,
        imsize=1000)

# skeletonize raster
riv_skel = os.path.join(work, '{0}_riv_skeleton.tif'.format(delta))
env.Command(
        source=riv_rast,
        target=riv_skel,
        action=lib.skeleton_riv,
        closing=5,
        holethresh=1000)

# drop small rivers
riv_dropped_small = os.path.join(work, '{0}_riv_dropped_pieces.tif'.format(delta))
env.Command(
        source=riv_skel,
        target=riv_dropped_small,
        action=lib.keep_n_rivers,
        n=1)

riv_clean = os.path.join(work, '{0}_riv_cleaned.tif'.format(delta))
myCommand(
        source=riv_dropped_small,
        target=riv_clean,
        action=lib.trim_short_rivs,
        minlen=40)

bifur_grid = os.path.join(work,'{0}_bifurs.tif'.format(delta))
env.Command(
        source=riv_clean,
        target=bifur_grid,
        action=lib.find_bifurs)


# add record id column to network cell table
network = os.path.join(work, '{0}_{1}_network.gdbn'.format(delta, STNres))
env.Command(
        source=STNnetwork,
        target=network,
        action=os.path.join(GHAASBIN, 'tblAddIdXY') + ' $SOURCE $TARGET')

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
                os.path.join(GHAASBIN,'netCells2Grid') + ' -f {0} -t {0} -u {0} -d {1} $SOURCE ${{TARGET}}.1'.format(varname, delta),
                os.path.join(GHAASBIN, 'grdRenameLayers') + ' -r 1 XXXX ${TARGET}.1 ${TARGET}.2',
                os.path.join(GHAASBIN, 'grdDateLayers') + ' -y 1 -e day ${TARGET}.2 ${TARGET}.3',
                os.path.join(GHAASBIN, 'rgis2netcdf') + ' ${TARGET}.3 $TARGET'])
    env.Command(
            source=path.format(ext='nc'),
            target=path.format(ext='tif'),
            action=lib.georef_nc)

network = os.path.join(work, '{0}_{1}_network.nx.yaml'.format(delta, STNres))
networkdelta = os.path.join(work, '{0}_{1}_network_delta.nx.yaml'.format(delta, STNres))
env.Command(
        source=[cellid.format(ext='tif'),
                basins.format(ext='tif'),
                flowdir.format(ext='tif'),
                bifur_grid,
                proj4str],
        target=[network, networkdelta],
        action=lib.import_rgis_network)

env.Command(
        source=[networkdelta, bifur_grid],
        target='output/{0}_{1}_map.png'.format(delta, STNres),
        action=[lib.plot_network_map,
                'convert -trim $TARGET $TARGET'],
            )


bifurs = os.path.join(output, 'bifurcations.csv')
env.Command(
        source=[networkdelta, bifur_grid, basins.format(ext='tif')],
        target=bifurs,
        action=lib.simple_bifurcations) # finds where osm river hits a subbasin that isn't on the main basin, and creates a bifurcation there

bifurs = os.path.join(output, 'bifurcations2.csv')
bifurnetwork = os.path.join(work, '{0}_{1}_network_delta_bifur.nx.yaml'.format(delta, STNres))
env.Command(
        source=[networkdelta, bifur_grid, basins.format(ext='tif')],
        target=[bifurs, bifurnetwork],
        action=lib.remap_riv_network) # more complete remapping of network to match osm rivers

env.Command(
        source=[bifurnetwork, bifur_grid],
        target='output/{0}_{1}_bifur_map.png'.format(delta, STNres),
        action=[lib.plot_network_map,
                'convert -trim $TARGET $TARGET'],
            )
