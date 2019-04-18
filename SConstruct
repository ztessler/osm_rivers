# vim: fileencoding=UTF-8
# vim: filetype=python

### TO RUN:

import os
import sys
import hashlib
import glob
from collections import defaultdict

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


alldeltas = ['Mekong', 'Chao_Phraya', 'Godavari']
allOSMrivers = {'Mekong': 'vietnam:cambodia',
                'Chao_Phraya': 'thailand',
                'Godavari': 'india',
                }
params = defaultdict(lambda:{'wetlands': False, 'thinning': 40, 'minwaterway_len': 0, 'minwidth':
    40, 'minarea': 50000}) #, 'minlen': 40})
params['Mekong'] = {'wetlands': False, 'thinning': 300, 'minwaterway_len': 30000, 'minwidth': 300,
        'minarea': 50000} #, 'minlen': 40}
        #'Chao_Phraya': {'wetlands': False, 'thinning': 20, 'minwaterway_len': 0, 'minarea': 0},# 'minlen': 40},
        #'Godavari': {'wetlands': False, 'thinning': 20, 'minwaterway_len': 0, 'minarea': 0},# 'minlen': 40},
        #}

GHAASBIN = env['ENV']['GHAASBIN']
STNres = os.environ.get('STNres', '06min')
INSPECTfig = os.environ.get('INSPECTfig', None)

deltas = os.environ.get('DELTA', alldeltas)
if not isinstance(deltas, list):
    deltas = [deltas]

for delta in deltas:
    work = os.path.join('work', delta)
    reswork = os.path.join('work', delta, STNres)
    output = os.path.join('output', delta, STNres)
    figures = os.path.join('figures', delta)
    resfigures = os.path.join('figures', delta, STNres)

    OSMrivers = allOSMrivers[delta].split(':')
    STNnetwork = '/Users/ecr/ztessler/projects/CHART/WBM/tools/buildNetwork/output/{delta}/{res}/{delta}_Network_{res}.gdbn'.format(delta=delta, res=STNres)
    SSEAnetwork = '/Users/ecr/ztessler/projects/CHART/WBM/tools/buildNetwork/output/SSEA/{res}/SSEA_Network_{res}.gdbn'.format(res=STNres)
    OSMshps = ['/Users/ecr/ztessler/data/OpenStreetMaps/Geofabrik/{0}/gis_osm_water_a_free_1.shp'.format(OSMriver) for OSMriver in OSMrivers]
    OSMwaterways = ['/Users/ecr/ztessler/data/OpenStreetMaps/Geofabrik/{0}/gis_osm_waterways_free_1.shp'.format(OSMriver) for OSMriver in OSMrivers]
    deltashp = '/Users/ecr/ztessler/data/deltas_LCLUC/maps/{0}_shp/{0}.shp'.format(delta)
    GSHHSshp = '/Users/ecr/ztessler/data/Coastline/GSHHG/gshhg-shp-2.3.6/GSHHS_shp/f/GSHHS_f_L1.shp'
    GRWLshps = '/Users/ecr/ztessler/data/Rivers/GRWL/GRWL_vector_V01.01/*.shp'

    thumbnail_size = 300

    # merge multiple country-level data if necessary
    if len(OSMshps) > 1:
        merged_shps = os.path.join(work, 'merged_rivs', 'merged_rivs.shp')
        env.Command(
                source=OSMshps,
                target=merged_shps,
                action='ogrmerge.py -single -lco ENCODING=UTF-8 -o $TARGET $SOURCES')
        merged_waterways = os.path.join(work, 'merged_rivs', 'merged_waterways.shp')
        env.Command(
                source=OSMwaterways,
                target=merged_waterways,
                action='ogrmerge.py -single -lco ENCODING=UTF-8 -o $TARGET $SOURCES')
    else:
        merged_shps = OSMshps[0]
        merged_waterways = OSMwaterways[0]

    proj4str = os.path.join(work, '{}_proj4.txt'.format(delta))
    myCommand(
            source=deltashp,
            target=proj4str,
            action=lib.set_projection)

    grwl_shp_list = os.path.join(work, 'grwl_shp_list.txt')
    grwl_shps = glob.glob(GRWLshps)
    myCommand(
            source=[deltashp]+grwl_shps,
            target=grwl_shp_list,
            action=lib.find_grwl_list)
    merged_GRWL_ll = os.path.join(work, '{0}_grwl_ll'.format(delta), '{0}_grwl_ll.shp'.format(delta))
    myCommand(
            source=[grwl_shp_list]+grwl_shps,
            target=merged_GRWL_ll,
            action='cat ${SOURCES[0]} | xargs ogrmerge.py -single -lco ENCODING=UTF-8 -o $TARGET')
    merged_GRWL = os.path.join(work, '{0}_grwl'.format(delta), '{0}_grwl.shp'.format(delta))
    myCommand(
            source=[proj4str, merged_GRWL_ll],
            target=merged_GRWL,
            action='cat ${SOURCES[0]} | xargs -I {} ogr2ogr -t_srs {} $TARGET ${SOURCES[1]}')

    # project and clip river vectors to delta
    clipped_vec = os.path.join(work, '{0}_riv_clipped/{0}_riv_clipped.shp'.format(delta))
    myCommand(
            source=[merged_shps, deltashp, proj4str],
            target=clipped_vec,
            action=lib.project_and_clip_osm_rivers)
    clipped_ww_vec = os.path.join(work, '{0}_ww_clipped/{0}_ww_clipped.shp'.format(delta))
    ww_proj4str = os.path.join(work, '{}_ww_proj4.txt'.format(delta))
    myCommand(
            source=[merged_waterways, deltashp, proj4str],
            target=clipped_ww_vec,
            action=lib.project_and_clip_osm_waterways)

    cleaned_ww_vec0 = os.path.join(work, '{0}_cleaned_ww_vec0/{0}_cleaned_ww_vec0.shp'.format(delta))
    env.Command(
            source=clipped_ww_vec,
            target=cleaned_ww_vec0,
            action='ogr2ogr -nlt LINESTRING -explodecollections -lco ENCODING=UTF-8 $TARGET $SOURCE')

    cleaned_ww_vec = os.path.join(work, '{0}_cleaned_ww_vec/{0}_cleaned_ww_vec.shp'.format(delta))
    env.Command(
            source=cleaned_ww_vec0,
            target=cleaned_ww_vec,
            action=lib.remove_ww_duplicates)

    clipped_coastline = os.path.join(work,
            '{0}_coastline_clipped/{0}_coastline_clipped.shp'.format(delta))
    clipped_coast = os.path.join(work,
            '{0}_coast_clipped/{0}_coast_clipped.shp'.format(delta))
    myCommand(
            source=[GSHHSshp, deltashp, proj4str],
            target=[clipped_coastline, clipped_coast],
            action=lib.project_and_clip_coastline)

    p = myCommand(
            source=clipped_vec,
            target=os.path.join(figures, '{}_vec_rivs_full.png'.format(delta)),
            action=[lib.plot_vec_rivs,
                    'convert -trim $TARGET $TARGET'])
    env.Default(p)
    p = myCommand(
            source=os.path.join(figures, '{}_vec_rivs_full.png'.format(delta)),
            target=os.path.join(figures, '{}_vec_rivs.png'.format(delta)),
            action='convert -fuzz 40% -trim -trim -resize {0} $SOURCE $TARGET'.format(thumbnail_size))
    env.Default(p)

    p = myCommand(
            source=cleaned_ww_vec,
            target=os.path.join(figures, '{}_ww_vec_rivs_full.png'.format(delta)),
            action=[lib.plot_vec_rivs,
                    'convert -trim $TARGET $TARGET'])
    env.Default(p)
    p = myCommand(
            source=os.path.join(figures, '{}_ww_vec_rivs_full.png'.format(delta)),
            target=os.path.join(figures, '{}_ww_vec_rivs.png'.format(delta)),
            action='convert -fuzz 40% -trim -trim -resize {0} $SOURCE $TARGET'.format(thumbnail_size))
    env.Default(p)

    rivers_ww_vec = os.path.join(work, '{0}_ww_rivers/{0}_ww_rivers.shp'.format(delta))
    myCommand(
            source=cleaned_ww_vec,
            target=rivers_ww_vec,
            action=lib.select_waterway_rivers)
    p = myCommand(
            source=rivers_ww_vec,
            target=os.path.join(figures, '{}_ww_rivers_vec_rivs_full.png'.format(delta)),
            action=[lib.plot_vec_rivs,
                    'convert -trim $TARGET $TARGET'])
    env.Default(p)
    p = myCommand(
            source=os.path.join(figures, '{}_ww_rivers_vec_rivs_full.png'.format(delta)),
            target=os.path.join(figures, '{}_ww_rivers_vec_rivs.png'.format(delta)),
            action='convert -fuzz 40% -trim -trim -resize {0} $SOURCE $TARGET'.format(thumbnail_size))
    env.Default(p)

    filtered_rivers_ww_vec = os.path.join(work,'{0}_filtered_ww_rivers/{0}_filtered_ww_rivers.shp'.format(delta))
    p = myCommand(
            source=rivers_ww_vec,
            target=filtered_rivers_ww_vec,
            action=lib.filter_waterway_rivers,
            minwaterway_len=params[delta]['minwaterway_len'])
    p = myCommand(
            source=filtered_rivers_ww_vec,
            target=os.path.join(figures, '{}_ww_filtered_rivers_full.png'.format(delta)),
            action=[lib.plot_vec_rivs,
                    'convert -trim $TARGET $TARGET'])
    env.Default(p)
    p = myCommand(
            source=os.path.join(figures, '{}_ww_filtered_rivers_full.png'.format(delta)),
            target=os.path.join(figures, '{}_ww_filtered_rivers.png'.format(delta)),
            action='convert -fuzz 40% -trim -trim -resize {0} $SOURCE $TARGET'.format(thumbnail_size))
    env.Default(p)

    # use the filtered waterway lines to correct the direction
    # apply corrections to full unfiltered set also so the dense version can be used in direction-setting
    corrected_ww_vec = os.path.join(work,'{0}_corrected_ww_vec/{0}_corrected_ww_vec.shp'.format(delta))
    filt_corrected_ww_vec = os.path.join(work,'{0}_filt_corrected_ww_vec/{0}_filt_corrected_ww_vec.shp'.format(delta))
    p = myCommand(
            source=[filtered_rivers_ww_vec, cleaned_ww_vec],
            target=[filt_corrected_ww_vec, corrected_ww_vec],
            action=lib.correct_waterway_flowdir)


    filtered_vec = os.path.join(work, '{0}_filtered_vec/{0}_filtered_vec.shp'.format(delta))
    p = myCommand(
            source=clipped_vec,
            target=filtered_vec,
            action=lib.filter_river_types,
            wetlands=params[delta]['wetlands'])
    p = myCommand(
            source=filtered_vec,
            target=os.path.join(figures, '{}_filtered_vec_rivs_full.png'.format(delta)),
            action=[lib.plot_vec_rivs,
                    'convert -trim $TARGET $TARGET'])
    env.Default(p)
    p = myCommand(
            source=os.path.join(figures, '{}_filtered_vec_rivs_full.png'.format(delta)),
            target=os.path.join(figures, '{}_filtered_vec_rivs.png'.format(delta)),
            action='convert -fuzz 40% -trim -trim -resize {0} $SOURCE $TARGET'.format(thumbnail_size))
    env.Default(p)

    width_est_vec = os.path.join(work, '{0}_riv_width_est/{0}_riv_width_est.shp'.format(delta))
    myCommand(
            source=filtered_vec,
            target=width_est_vec,
            action=lib.get_river_widths)

    thinned_vec = os.path.join(work, '{0}_riv_thinned/{0}_riv_thinned.shp'.format(delta))
    myCommand(
            #source=filtered_vec,
            source=width_est_vec,
            target=thinned_vec,
            action=lib.thin_vec,
            thinning=params[delta]['thinning'],
            minwidth=params[delta]['minwidth'],
            minarea=params[delta]['minarea'])
            #minhole=params[delta]['minhole'])
    p = myCommand(
            source=thinned_vec,
            target=os.path.join(figures, '{}_vec_thinned_rivs_full.png'.format(delta)),
            action=[lib.plot_vec_rivs,
                    'convert -trim $TARGET $TARGET'])
    env.Default(p)
    p = myCommand(
            source=os.path.join(figures, '{}_vec_thinned_rivs_full.png'.format(delta)),
            target=os.path.join(figures, '{}_vec_thinned_rivs.png'.format(delta)),
            action='convert -fuzz 40% -trim -trim -resize {0} $SOURCE $TARGET'.format(thumbnail_size))
    env.Default(p)

    merged_vec = os.path.join(work, '{0}_riv_merged/{0}_riv_merged.shp'.format(delta))
    myCommand(
            source=[thinned_vec, filtered_vec, filt_corrected_ww_vec],
            target=merged_vec,
            action=lib.merge_water_waterway_vecs,
            buff=params[delta]['thinning'])
    p = myCommand(
            source=merged_vec,
            target=os.path.join(figures, '{}_vec_merged_rivs_full.png'.format(delta)),
            action=[lib.plot_vec_rivs,
                    'convert -trim $TARGET $TARGET'])
    env.Default(p)
    p = myCommand(
            source=os.path.join(figures, '{}_vec_merged_rivs_full.png'.format(delta)),
            target=os.path.join(figures, '{}_vec_merged_rivs.png'.format(delta)),
            action='convert -fuzz 40% -trim -trim -resize {0} $SOURCE $TARGET'.format(thumbnail_size))
    env.Default(p)



    # rasterize
    riv_rast = os.path.join(work, '{0}_riv_rast.tif'.format(delta))
    riv_affine = os.path.join(work, '{0}_riv_affine.pkl'.format(delta))
    riv_shape = os.path.join(work, '{0}_riv_shape.pkl'.format(delta))
    myCommand(
            source=merged_vec,
            target=[riv_rast, riv_affine, riv_shape],
            action=lib.rasterize_riv,
            imsize=1000)
    p = env.Command(
            source=riv_rast,
            target=os.path.join(figures, '{}_riv_rast_full.0.png').format(delta),
            action='convert -negate -normalize $SOURCE $TARGET')
    env.Default(p)
    p = env.Command(
            source=riv_rast,
            target=os.path.join(figures, '{}_riv_rast.0.png').format(delta),
            action='convert -resize {0} -negate -normalize $SOURCE $TARGET'.format(thumbnail_size))
    env.Default(p)

    # skeletonize raster
    riv_skel = os.path.join(work, '{0}_riv_skeleton.tif'.format(delta))
    myCommand(
            source=riv_rast,
            target=riv_skel,
            action=lib.skeleton_riv,
            closing=5,
            holethresh=1000)
    p = env.Command(
            source=riv_skel,
            target=os.path.join(figures, '{}_riv_skel.1.png').format(delta),
            action='convert -resize {0} -negate -normalize $SOURCE $TARGET'.format(thumbnail_size))
    env.Default(p)

    # drop small rivers
    riv_clean = os.path.join(work, '{0}_riv_cleaned.tif'.format(delta))
    myCommand(
            source=riv_skel,
            target=riv_clean,
            action=lib.keep_n_rivers,
            n=1)
    p = env.Command(
            source=riv_clean,
            target=os.path.join(figures, '{}_riv_clean.png').format(delta),
            action='convert -resize {0} -negate -normalize $SOURCE $TARGET'.format(thumbnail_size))
    env.Default(p)

    bifur_grid = os.path.join(work,'{0}_bifurs.tif'.format(delta))
    env.Command(
            source=riv_clean,
            target=bifur_grid,
            action=lib.find_bifurs)

    segments1 = os.path.join(work, 'river_segments.1.pkl')
    myCommand(
            source=bifur_grid,
            target=segments1,
            action=lib.find_river_segments)
    myCommand(
            source=[bifur_grid, segments1],
            target=os.path.join(resfigures, '{0}_{1}_river_flowdirs.0.png'.format(delta, STNres)),
            action=[lib.plot_flowdirs_map,
                    'convert -trim $TARGET $TARGET'],
            inspect=INSPECTfig)

    # add record id column to network cell table
    network = os.path.join(reswork, '{0}_network.gdbn'.format(STNres))
    env.Command(
            source=STNnetwork,
            target=network,
            action=os.path.join(GHAASBIN, 'tblAddIdXY') + ' $SOURCE $TARGET')

    # import RGIS network
    cellid = os.path.join(reswork, '{0}_cellid.{{ext}}'.format(STNres))
    basins = os.path.join(reswork, '{0}_basins.{{ext}}'.format(STNres))
    flowdir = os.path.join(reswork, '{0}_flowdir.{{ext}}'.format(STNres))
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

    network = os.path.join(reswork, '{0}_{1}_network_fullbasins.nx.pkl'.format(delta, STNres))
    networkdelta = os.path.join(reswork, '{0}_{1}_network_delta.nx.pkl'.format(delta, STNres))
    nupstream = os.path.join(reswork, 'nupstream.pkl')
    ndownstream = os.path.join(reswork, 'ndownstream.pkl')
    nodepositions = os.path.join(reswork, 'nodepositions.pkl')
    env.Command(
            source=[cellid.format(ext='tif'),
                    basins.format(ext='tif'),
                    flowdir.format(ext='tif'),
                    riv_affine, riv_shape,
                    proj4str],
            target=[network, networkdelta, nupstream, ndownstream, nodepositions],
            action=lib.import_rgis_network)

    nearestnodes1 = os.path.join(reswork, 'nearestnodes.1.pkl')
    myCommand(
            source=[networkdelta, bifur_grid],
            target=nearestnodes1,
            action=lib.find_nearest_nodes_to_riv)

    node_dist_to_coast = os.path.join(reswork, 'node_dist_to_coast.pkl')
    riv_dist_to_coast = os.path.join(reswork, 'riv_dist_to_coast.pkl')
    myCommand(
            source=[networkdelta, bifur_grid, clipped_coastline, clipped_coast],
            target=[node_dist_to_coast, riv_dist_to_coast],
            action=lib.calc_dist_to_coast)

    riv_flowdist_to_coast = os.path.join(reswork, 'riv_flowdist_to_coast.tif')
    myCommand(
            source=[bifur_grid, riv_dist_to_coast],
            target=riv_flowdist_to_coast,
            action=lib.calc_riv_flowdist_to_coast)

    head_rivpt = os.path.join(reswork, 'head_rivpt.1.pkl')
    myCommand(
            source=[bifur_grid, nearestnodes1, ndownstream, nodepositions],
            target=head_rivpt,
            action=lib.find_head_rivpt)

    segments2 = os.path.join(reswork, 'river_segments.2.pkl')
    myCommand(
            source=[segments1, bifur_grid, corrected_ww_vec, riv_flowdist_to_coast],
            target=segments2,
            action=lib.set_segment_flowdir)

    next_rivpts = os.path.join(reswork, 'next_rivpts.pkl')
    prev_rivpts = os.path.join(reswork, 'prev_rivpts.pkl')
    myCommand(
            source=segments2,
            target=[next_rivpts, prev_rivpts],
            action=lib.next_prev_pts)

    # extend rivers (at coastal end) out toward/past coastline, some initial rgis networks extend
    # farther out to direct estuary flow.
    river_adj1 = os.path.join(reswork, '{0}_river_adj_to_network.1.tif'.format(delta))
    myCommand(
            source=[bifur_grid, next_rivpts, prev_rivpts, riv_dist_to_coast, clipped_coastline, clipped_coast],
            target=river_adj1,
            action=lib.extend_rivers_to_coast)

    # final river version, cleaned and merged network. put in reswork dir since depends on rgis res
    river_adj = os.path.join(reswork, '{0}_river_adj_to_network.tif'.format(delta))
    myCommand(
            source=[river_adj1, nupstream, ndownstream, nodepositions],
            target=river_adj,
            action=lib.merge_riv_path_to_mainstem)

    bifur_adj = os.path.join(reswork,'{0}_adj_bifurs.tif'.format(delta))
    env.Command(
            source=river_adj,
            target=bifur_adj,
            action=lib.find_bifurs)

    node_dist_to_coast = os.path.join(reswork, 'node_dist_to_coast.1.pkl') # same as other
    riv_dist_to_coast = os.path.join(reswork, 'riv_adj_dist_to_coast.pkl')
    myCommand(
            source=[networkdelta, bifur_adj, clipped_coastline, clipped_coast],
            target=[node_dist_to_coast, riv_dist_to_coast],
            action=lib.calc_dist_to_coast)

    segments3 = os.path.join(reswork, '{0}_river_segments.3.pkl'.format(delta))
    myCommand(
            source=bifur_adj,
            target=segments3,
            action=lib.find_river_segments)

    segments = os.path.join(reswork, '{0}_river_segments.pkl'.format(delta))
    myCommand(
            source=[segments3, bifur_adj, corrected_ww_vec, riv_flowdist_to_coast],
            target=segments,
            action=lib.set_segment_flowdir)

    river_widths = os.path.join(reswork, '{0}_river_widths.pkl'.format(delta))
    river_widths_rast = os.path.join(reswork, '{0}_river_widths.tif'.format(delta))
    myCommand(
            source=[segments, bifur_adj, merged_GRWL],
            target=[river_widths, river_widths_rast],
            action=lib.set_segment_widths)

    next_rivpts = os.path.join(reswork, '{0}_next_rivpts.pkl'.format(delta))
    prev_rivpts = os.path.join(reswork, '{0}_prev_rivpts.pkl'.format(delta))
    myCommand(
            source=segments,
            target=[next_rivpts, prev_rivpts],
            action=lib.next_prev_pts)

    nearestnodes = os.path.join(reswork, '{0}_nearestnodes.pkl'.format(delta))
    myCommand(
            source=[networkdelta, bifur_adj],
            target=nearestnodes,
            action=lib.find_nearest_nodes_to_riv)


    head_rivpt = os.path.join(reswork, '{0}_head_rivpt.pkl'.format(delta))
    myCommand(
            source=[bifur_adj, nearestnodes, ndownstream, nodepositions],
            target=head_rivpt,
            action=lib.find_head_rivpt)

    bifurs = os.path.join(output, '{0}_{1}_bifurcations.csv'.format(delta, STNres))
    bifurnetwork = os.path.join(reswork, '{0}_{1}_network_delta_bifur.nx.pkl'.format(delta, STNres))
    bifuroutlets = os.path.join(output, '{0}_{1}_bifur_outlet_cellids.csv'.format(delta, STNres))
    b = myCommand(
            source=[networkdelta, network, bifur_adj, head_rivpt, next_rivpts, prev_rivpts, nearestnodes,
                    riv_dist_to_coast, river_widths],
            target=[bifurs, bifurnetwork, bifuroutlets],
            action=lib.remap_riv_network)
    env.Default(b)

    for networkversion, network_name in [(bifurnetwork, 'bifur_'), (networkdelta, '')]:
        for labels, label_name in [('none', ''), ('nodes', '_nodes'), ('cells', '_cells'), ('outlets', '_outlets')]:
            p = myCommand(
                    source=[networkversion, bifur_grid, bifur_adj, bifuroutlets, proj4str],
                    target=os.path.join(resfigures, '{0}_{1}_{2}map{3}.png'.format(delta, STNres, network_name, label_name)),
                    action=[lib.plot_network_map,
                            'convert -trim $TARGET $TARGET'],
                    labels=labels,
                    inspect=INSPECTfig)
            env.Default(p)
    for labels, label_name in [('none', ''), ('nodes', '_nodes'), ('cells', '_cells'), ('outlets', '_outlets')]:
        p = myCommand(
                source=[networkdelta, bifurnetwork, bifur_grid, bifur_adj, bifuroutlets, proj4str],
                target=os.path.join(resfigures, '{0}_{1}_bifur_map{2}_diff.png'.format(delta, STNres, label_name)),
                action=[lib.plot_network_diff_map,
                        'convert -trim $TARGET $TARGET'],
                labels=labels,
                inspect=INSPECTfig)
        env.Default(p)
    p = myCommand(
            source=[bifur_adj, segments],
            target=os.path.join(resfigures, '{0}_{1}_river_flowdirs.png'.format(delta, STNres)),
            action=[lib.plot_flowdirs_map,
                    'convert -trim $TARGET $TARGET'],
            inspect=INSPECTfig)
    env.Default(p)

    bifurnetworkgraphml = os.path.join(reswork, '{0}_{1}_network_delta_bifur.nx.graphml'.format(delta, STNres))
    t = myCommand(
            source=bifurnetwork,
            target=bifurnetworkgraphml,
            action=lib.convert_network_to_graphml)
    env.Default(t)

    STNcells = os.path.join(reswork, 'STNnetwork_dbcells.txt')
    env.Command(
            source=STNnetwork,
            target=STNcells,
            action='rgis2table -a DBCells $SOURCE > $TARGET')
    SSEAcells = os.path.join(reswork, 'SSEAnetwork_dbcells.txt')
    env.Command(
            source=SSEAnetwork,
            target=SSEAcells,
            action='rgis2table -a DBCells $SOURCE > $TARGET')
    SSEA_bifurs = os.path.join(output, '{0}_{1}_SSEA_bifurcations.csv'.format(delta, STNres))
    b = env.Command(
            source=[STNcells, SSEAcells, bifurs],
            target=SSEA_bifurs,
            action=lib.convert_bifur_cellids_to_SSEA)
    env.Default(b)
