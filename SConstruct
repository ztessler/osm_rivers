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
OSMrivers = os.environ.get('OSMriver', 'vietnam:cambodia').split(':')
STNres = os.environ.get('STNres', '06min')
INSPECTfig = os.environ.get('INSPECTfig', None)

deltawork = os.path.join('work', delta) # SSEA domain can reuse some single-domain files
domainwork = os.path.join('work', domain, delta, STNres)
domain_nores_work = os.path.join('work', domain, delta)
output = os.path.join('output', domain, delta, STNres)
deltafigures = os.path.join('figures', delta)
domainfigures = os.path.join('figures', domain, delta, STNres)

STNnetwork = '/Users/ecr/ztessler/projects/CHART/WBM/tools/buildNetwork/output/{domain}/{res}/{domain}_Network_{res}.gdbn'.format(domain=domain, res=STNres)
OSMshps = ['/Users/ecr/ztessler/data/OpenStreetMaps/Geofabrik/{0}/gis_osm_water_a_free_1.shp'.format(OSMriver) for OSMriver in OSMrivers]
OSMwaterways = ['/Users/ecr/ztessler/data/OpenStreetMaps/Geofabrik/{0}/gis_osm_waterways_free_1.shp'.format(OSMriver) for OSMriver in OSMrivers]
deltashp = '/Users/ecr/ztessler/data/deltas_LCLUC/maps/{0}_shp/{0}.shp'.format(delta)

thumbnail_size = 300


params = {
        'Mekong': {'wetlands': False, 'thinning': 100, 'minarea': 0, 'minlen': 40, 'flowdir_weights': (1,2)},
        'Chao_Phraya': {'wetlands': False, 'thinning': 20, 'minarea': 10000, 'minlen': 80, 'flowdir_weights': (5,1)},
        }

# merge multiple country-level data if necessary
if len(OSMshps) > 1:
    merged_shps = os.path.join(domain_nores_work, 'merged_rivs.shp')
    env.Command(
            source=OSMshps,
            target=merged_shps,
            action='ogrmerge.py -single -lco ENCODING=UTF-8 -o $TARGET $SOURCES')
    merged_waterways = os.path.join(domain_nores_work, 'merged_waterways.shp')
    env.Command(
            source=OSMwaterways,
            target=merged_waterways,
            action='ogrmerge.py -single -lco ENCODING=UTF-8 -o $TARGET $SOURCES')
else:
    merged_shps = OSMshps[0]
    merged_waterways = OSMwaterways[0]

# project and clip river vectors to delta
clipped_vec = os.path.join(deltawork, '{0}_riv_clipped/{0}_riv_clipped.shp'.format(delta))
proj4str = os.path.join(deltawork, '{}_proj4.txt'.format(delta))
myCommand(
        source=[merged_shps, deltashp],
        target=[clipped_vec, proj4str],
        action=lib.project_and_clip_osm_rivers)
clipped_ww_vec = os.path.join(deltawork, '{0}_ww_clipped/{0}_ww_clipped.shp'.format(delta))
ww_proj4str = os.path.join(deltawork, '{}_ww_proj4.txt'.format(delta))
myCommand(
        source=[merged_waterways, deltashp],
        target=[clipped_ww_vec, ww_proj4str],
        action=lib.project_and_clip_osm_waterways)

p = myCommand(
        source=clipped_vec,
        target=os.path.join(deltafigures, '{}_vec_rivs_full.png'.format(delta)),
        action=[lib.plot_vec_rivs,
                'convert -trim $TARGET $TARGET'])
env.Default(p)
p = myCommand(
        source=os.path.join(deltafigures, '{}_vec_rivs_full.png'.format(delta)),
        target=os.path.join(deltafigures, '{}_vec_rivs.png'.format(delta)),
        action='convert -fuzz 40% -trim -trim -resize {0} $SOURCE $TARGET'.format(thumbnail_size))
env.Default(p)

p = myCommand(
        source=clipped_ww_vec,
        target=os.path.join(deltafigures, '{}_ww_vec_rivs_full.png'.format(delta)),
        action=[lib.plot_vec_rivs,
                'convert -trim $TARGET $TARGET'])
env.Default(p)
p = myCommand(
        source=os.path.join(deltafigures, '{}_ww_vec_rivs_full.png'.format(delta)),
        target=os.path.join(deltafigures, '{}_ww_vec_rivs.png'.format(delta)),
        action='convert -fuzz 40% -trim -trim -resize {0} $SOURCE $TARGET'.format(thumbnail_size))
env.Default(p)

filtered_ww_vec = os.path.join(deltawork, '{0}_ww_filtered/{0}_ww_filtered.shp'.format(delta))
myCommand(
        source=clipped_ww_vec,
        target=filtered_ww_vec,
        action=lib.filter_waterway_types)
p = myCommand(
        source=filtered_ww_vec,
        target=os.path.join(deltafigures, '{}_ww_filtered_vec_rivs_full.png'.format(delta)),
        action=[lib.plot_vec_rivs,
                'convert -trim $TARGET $TARGET'])
env.Default(p)
p = myCommand(
        source=os.path.join(deltafigures, '{}_ww_filtered_vec_rivs_full.png'.format(delta)),
        target=os.path.join(deltafigures, '{}_ww_filtered_vec_rivs.png'.format(delta)),
        action='convert -fuzz 40% -trim -trim -resize {0} $SOURCE $TARGET'.format(thumbnail_size))
env.Default(p)


thinned_vec = os.path.join(deltawork, '{0}_riv_thinned/{0}_riv_thinned.shp'.format(delta))
myCommand(
        source=clipped_vec,
        target=thinned_vec,
        action=lib.thin_vec,
        wetlands=params[delta]['wetlands'],
        thinning=params[delta]['thinning'],
        minarea=params[delta]['minarea'])
        #minhole=params[delta]['minhole'])
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

bifur_grid1 = os.path.join(deltawork,'{0}_bifurs.1.tif'.format(delta))
env.Command(
        source=riv_dropped_small,
        target=bifur_grid1,
        action=lib.find_bifurs)

riv_clean1 = os.path.join(deltawork, '{0}_riv_cleaned.1.tif'.format(delta))
myCommand(
        source=bifur_grid1,
        target=riv_clean1,
        action=lib.trim_short_rivs,
        minlen=params[delta]['minlen'])

bifur_grid = os.path.join(deltawork,'{0}_bifurs.2.tif'.format(delta))
env.Command(
        source=riv_clean1,
        target=bifur_grid,
        action=lib.find_bifurs)

segments_1 = os.path.join(deltawork, 'river_segments.1.pkl')
myCommand(
        source=bifur_grid,
        target=segments_1,
        action=lib.find_river_segments)

riv_clean = os.path.join(deltawork, '{0}_riv_cleaned.tif'.format(delta))
myCommand(
        source=[riv_clean1, segments_1],
        target=riv_clean,
        action=lib.remove_small_loops,
        minlen=params[delta]['minlen'])
#bifur_grid = os.path.join(deltawork, '{0}_bifurs.tif'

segments_2 = os.path.join(deltawork, 'river_segments.2.pkl')
myCommand(
        source=[segments_1, bifur_grid, filtered_ww_vec],
        target=segments_2,
        action=lib.set_segment_flowdir)

p = env.Command(
        source=riv_clean,
        target=os.path.join(deltafigures, '{}_riv_clean.3.png').format(delta),
        action='convert -resize {0} -negate -normalize $SOURCE $TARGET'.format(thumbnail_size))
env.Default(p)

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


nearestnodes_to_riv_1 = os.path.join(domainwork, 'nearestnode_to_riv.1.pkl')
myCommand(
        source=[networkdelta, bifur_grid],
        target=nearestnodes_to_riv_1,
        action=lib.find_nearest_nodes_to_riv)

node_dist_to_coast = os.path.join(domainwork, 'node_dist_to_coast.pkl')
myCommand(
        source=networkdelta,
        target=node_dist_to_coast,
        action=lib.calc_dist_to_coast)

bifurs = os.path.join(output, '{0}_{1}_{2}_bifurcations.csv'.format(domain, delta, STNres))
bifurnetwork = os.path.join(domainwork, '{0}_{1}_{2}_network_delta_bifur.nx.yaml'.format(domain, delta, STNres))
extended_bifurgrid = os.path.join(deltawork, '{0}_bifurs_extended.tif'.format(delta))
bifuroutlets = os.path.join(output, '{0}_{1}_{2}_bifur_outlet_cellids.csv'.format(domain, delta, STNres))
riversegments = os.path.join(output, '{0}_{1}_{2}_river_segments.pkl'.format(domain, delta, STNres))
flowdir = os.path.join(output, '{0}_{1}_{2}_river_flowdirs.pkl'.format(domain, delta, STNres))
b = myCommand(
        source=[networkdelta, bifur_grid, basins.format(ext='tif')],
        #source=[networkdelta, bifur_grid, nearestnodes_to_riv_1, node_dist_to_coast],
        target=[bifurs, bifurnetwork, extended_bifurgrid, bifuroutlets, riversegments, flowdir],
        action=lib.remap_riv_network,
        flowdir_weights=params[delta]['flowdir_weights']) # more complete remapping of network to match osm rivers
env.Default(b)

for networkversion, network_name in [(bifurnetwork, 'bifur_'), (networkdelta, '')]:
    for labels, label_name in [('none', ''), ('nodes', '_nodes'), ('cells', '_cells')]:
        p = myCommand(
                source=[networkversion, bifur_grid, extended_bifurgrid],
                target=os.path.join(domainfigures, '{0}_{1}_{2}_{3}map{4}.png'.format(domain, delta, STNres, network_name, label_name)),
                action=[lib.plot_network_map,
                        'convert -trim $TARGET $TARGET'],
                labels=labels,
                inspect=INSPECTfig)
        env.Default(p)
p = myCommand(
        source=[bifur_grid, extended_bifurgrid, riversegments, flowdir],
        target=os.path.join(domainfigures, '{0}_{1}_{2}_river_flowdirs.png'.format(domain, delta, STNres)),
        action=[lib.plot_flowdirs_map,
                'convert -trim $TARGET $TARGET'],
        inspect=INSPECTfig)
env.Default(p)
