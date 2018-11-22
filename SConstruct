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

domain = os.environ.get('DOMAIN', 'Mekong')
delta = os.environ.get('DELTA', 'Mekong')
#OSMrivers = os.environ.get('OSMriver', 'vietnam:cambodia').split(':')
OSMrivers = os.environ.get('OSMriver', 'vietnam').split(':')
STNres = os.environ.get('STNres', '06min')

deltawork = os.path.join('work', delta) # SSEA domain can reuse some single-domain files
domainwork = os.path.join('work', domain, delta, STNres)
domain_nores_work = os.path.join('work', domain, delta)
output = os.path.join('output', domain, STNres)
deltafigures = os.path.join('figures', delta)
domainfigures = os.path.join('figures', domain, delta, STNres)

STNnetwork = '/Users/ecr/ztessler/projects/CHART/WBM/tools/buildNetwork/output/{domain}/{res}/{domain}_Network_{res}.gdbn'.format(domain=domain, res=STNres)
#STNnetwork = '/Users/ecr/ztessler/projects/CHART/WBM/tools/buildNetwork/output/SSEA/{res}/SSEA_Network_{res}.gdbn'.format(res=STNres)
OSMshps = ['/Users/ecr/ztessler/projects/CHART/WBM/tools/osm_rivers/osm_data/{0}/gis.osm_water_a_free_1.shp'.format(OSMriver) for OSMriver in OSMrivers]
deltashp = '/Users/ecr/ztessler/data/deltas_LCLUC/maps/{0}_shp/{0}.shp'.format(delta)

thumbnail_size = 300

# merge multiple country-level data if necessary
if len(OSMshps) > 1:
    merged_shps = os.path.join(domain_nores_work, 'merged_rivs.shp')
    env.Command(
            source=OSMshps,
            target=merged_shps,
            action='ogrmerge.py -single -o $TARGET $SOURCES')
else:
    merged_shps = OSMshps[0]

# project and clip river vectors to delta
clipped_vec = os.path.join(deltawork, '{0}_riv_clipped/{0}_riv_clipped.shp'.format(delta))
proj4str = os.path.join(deltawork, '{}_proj4.txt'.format(delta))
myCommand(
        source=[merged_shps, deltashp],
        target=[clipped_vec, proj4str],
        action=lib.project_and_clip_osm_rivers)
p = myCommand(
        source=clipped_vec,
        target=os.path.join(deltafigures, '{}_vec_rivs.png'.format(delta)),
        action=[lib.plot_vec_rivs,
                'convert -fuzz 40% -trim -trim -resize {0} $TARGET $TARGET'.format(thumbnail_size)])
env.Default(p)
p = myCommand(
        source=clipped_vec,
        target=os.path.join(deltafigures, '{}_vec_rivs_full.png'.format(delta)),
        action=lib.plot_vec_rivs)
env.Default(p)

thinned_vec = os.path.join(deltawork, '{0}_riv_thinned/{0}_riv_thinned.shp'.format(delta))
myCommand(
        source=clipped_vec,
        target=thinned_vec,
        action=lib.thin_vec,
        thresh=100)
p = myCommand(
        source=thinned_vec,
        target=os.path.join(deltafigures, '{}_vec_thinned_rivs.png'.format(delta)),
        action=[lib.plot_vec_rivs,
                'convert -fuzz 40% -trim -trim -resize {0} $TARGET $TARGET'.format(thumbnail_size)])
env.Default(p)


# rasterize
riv_rast = os.path.join(deltawork, '{0}_riv_rast.tif'.format(delta))
myCommand(
        source=thinned_vec,
        target=riv_rast,
        action=lib.rasterize_riv,
        imsize=1000)
p = env.Command(
        source=riv_rast,
        target=os.path.join(deltafigures, '{}_riv_rast.0.png').format(delta),
        action='convert -resize {0} -negate -normalize $SOURCE $TARGET'.format(thumbnail_size))
env.Default(p)

# skeletonize raster
riv_skel = os.path.join(deltawork, '{0}_riv_skeleton.tif'.format(delta))
myCommand(
        source=riv_rast,
        target=riv_skel,
        action=lib.skeleton_riv,
        closing=5,
        holethresh=1000)
p = env.Command(
        source=riv_skel,
        target=os.path.join(deltafigures, '{}_riv_skel.1.png').format(delta),
        action='convert -resize {0} -negate -normalize $SOURCE $TARGET'.format(thumbnail_size))
env.Default(p)

# drop small rivers
riv_dropped_small = os.path.join(deltawork, '{0}_riv_dropped_pieces.tif'.format(delta))
myCommand(
        source=riv_skel,
        target=riv_dropped_small,
        action=lib.keep_n_rivers,
        n=1)
p = env.Command(
        source=riv_dropped_small,
        target=os.path.join(deltafigures, '{}_riv_dropped_small.2.png').format(delta),
        action='convert -resize {0} -negate -normalize $SOURCE $TARGET'.format(thumbnail_size))
env.Default(p)

bifur_precleaned = os.path.join(deltawork,'{0}_bifurs_precleaned.tif'.format(delta))
env.Command(
        source=riv_dropped_small,
        target=bifur_precleaned,
        action=lib.find_bifurs)

riv_clean = os.path.join(deltawork, '{0}_riv_cleaned.tif'.format(delta))
myCommand(
        source=bifur_precleaned,
        target=riv_clean,
        action=lib.trim_short_rivs,
        minlen=40)
p = env.Command(
        source=riv_clean,
        target=os.path.join(deltafigures, '{}_riv_clean.3.png').format(delta),
        action='convert -resize {0} -negate -normalize $SOURCE $TARGET'.format(thumbnail_size))
env.Default(p)

bifur_grid = os.path.join(deltawork,'{0}_bifurs.tif'.format(delta))
env.Command(
        source=riv_clean,
        target=bifur_grid,
        action=lib.find_bifurs)


# add record id column to network cell table
network = os.path.join(domainwork, '{0}_{1}_network.gdbn'.format(domain, STNres))
env.Command(
        source=STNnetwork,
        target=network,
        action=os.path.join(GHAASBIN, 'tblAddIdXY') + ' $SOURCE $TARGET')

# import RGIS network
cellid = os.path.join(domainwork, '{0}_{1}_cellid.{{ext}}'.format(domain, STNres))
basins = os.path.join(domainwork, '{0}_{1}_basins.{{ext}}'.format(domain, STNres))
flowdir = os.path.join(domainwork, '{0}_{1}_flowdir.{{ext}}'.format(domain, STNres))
for (path, varname) in [(cellid, 'CellID'),
                        (basins, 'BasinID'),
                        (flowdir, 'ToCell')]:
    env.Command(
            source=network,
            target=path.format(ext='nc'),
            action=[
                os.path.join(GHAASBIN,'netCells2Grid') + ' -f {0} -t {0} -u {0} -d {1} $SOURCE ${{TARGET}}.1'.format(varname, domain),
                os.path.join(GHAASBIN, 'grdRenameLayers') + ' -r 1 XXXX ${TARGET}.1 ${TARGET}.2',
                os.path.join(GHAASBIN, 'grdDateLayers') + ' -y 1 -e day ${TARGET}.2 ${TARGET}.3',
                os.path.join(GHAASBIN, 'rgis2netcdf') + ' ${TARGET}.3 $TARGET'])
    env.Command(
            source=path.format(ext='nc'),
            target=path.format(ext='tif'),
            action=lib.georef_nc)

network = os.path.join(domainwork, '{0}_{1}_network.nx.yaml'.format(domain, STNres))
networkdelta = os.path.join(domainwork, '{0}_{1}_{2}_network_delta.nx.yaml'.format(domain, delta, STNres))
env.Command(
        source=[cellid.format(ext='tif'),
                basins.format(ext='tif'),
                flowdir.format(ext='tif'),
                bifur_grid,
                proj4str],
        target=[network, networkdelta],
        action=lib.import_rgis_network)

#bifurs = os.path.join(output, '{0}_{1}_simple_bifurcations.csv'.format(delta, STNres))
#b = env.Command(
        #source=[networkdelta, bifur_grid, basins.format(ext='tif')],
        #target=bifurs,
        #action=lib.simple_bifurcations) # finds where osm river hits a subbasin that isn't on the main basin, and creates a bifurcation there
#env.Default(b)

bifurs = os.path.join(output, '{0}_{1}_{2}_bifurcations.csv'.format(domain, delta, STNres))
bifurnetwork = os.path.join(domainwork, '{0}_{1}_{2}_network_delta_bifur.nx.yaml'.format(domain, delta, STNres))
bifuroutlets = os.path.join(output, '{0}_{1}_{2}_bifur_outlet_cellids.csv'.format(domain, delta, STNres))
riversegments = os.path.join(output, '{0}_{1}_{2}_river_segments.pkl'.format(domain, delta, STNres))
flowdir = os.path.join(output, '{0}_{1}_{2}_river_flowdirs.pkl'.format(domain, delta, STNres))
b = env.Command(
        source=[networkdelta, bifur_grid, basins.format(ext='tif')],
        target=[bifurs, bifurnetwork, bifuroutlets, riversegments, flowdir],
        action=lib.remap_riv_network) # more complete remapping of network to match osm rivers
env.Default(b)

for networkversion, network_name in [(bifurnetwork, 'bifur_'), (networkdelta, '')]:
    for labels, label_name in [('none', ''), ('nodes', '_nodes'), ('cells', '_cells')]:
        p = myCommand(
                source=[networkversion, bifur_grid],
                target=os.path.join(domainfigures, '{0}_{1}_{2}_{3}map{4}.png'.format(domain, delta, STNres, network_name, label_name)),
                action=[lib.plot_network_map,
                        'convert -trim $TARGET $TARGET'],
                labels=labels)
        env.Default(p)
p = myCommand(
        source=[bifur_grid, riversegments, flowdir],
        target=os.path.join(domainfigures, '{0}_{1}_{2}_river_flowdirs.png'.format(domain, delta, STNres)),
        action=[lib.plot_flowdirs_map,
                'convert -trim $TARGET $TARGET'])
env.Default(p)
