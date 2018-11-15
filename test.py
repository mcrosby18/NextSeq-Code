#!/usr/bin/env python

import argparse
import csv
import fnmatch
import glob
import json
import numpy
import os
import pandas
import re
import sys

from bitstring import BitString

#from boto.s3.connection    import S3Connection
#from boto.s3.key           import Key
#from boto.sqs.connection   import SQSConnection
#from boto.sqs.message      import Message
from os                    import environ
from xml.etree.ElementTree import ElementTree

#from ggplot import *


def valid_runfolder(path):

    bn = os.path.basename(path.rstrip(os.sep))

    if (re.search('^[0-9]{6}_[A-Za-z0-9\-]+_[0-9]{4}_[A-Za-z0-9\-]+(_[0-9]{4})?$', bn)):

        if (os.path.exists(path)):
            return True

    return False

def lookup_runinfo(path):

    dict = { 'folder_location' : path }

    runparameters_xml      = os.path.join(path, 'runParameters.xml')
    runparameters_xml_rta3 = os.path.join(path, 'RunParameters.xml')

    tree = ElementTree()

    if (os.path.exists(runparameters_xml_rta3)):
        tree.parse(runparameters_xml_rta3)
    else:
        tree.parse(runparameters_xml)

    exp_elem = tree.find("Setup/ExperimentName")

    if (exp_elem is None):
        exp_elem = tree.find("ExperimentName")

    dict['run_number'] = exp_elem.text.replace(" ", "_")

    control_lane_elem = tree.find("Setup/ControlLane")

    if (control_lane_elem is not None):
        dict['control_lane'] = int(control_lane_elem.text)

    rta_version_elem = tree.find("Setup/RTAVersion")

    if (rta_version_elem is None):
        rta_version_elem = tree.find("RTAVersion")

    if (rta_version_elem is None):
        rta_version_elem = tree.find("RtaVersion")

    dict['rta_version'] = rta_version_elem.text

    runinfo_xml = os.path.join(path, 'RunInfo.xml')

    tree = ElementTree()
    tree.parse(runinfo_xml)
    run_elem = tree.find("Run")
    dict['run_id'] = run_elem.attrib['Id']
    flowcell_elem = tree.find("Run/Flowcell")
    dict['flowcell_id']  = re.sub('0{9}-', '', flowcell_elem.text) 
    read_elem = tree.find("Run/Reads")
    reads = read_elem.findall("Read")

    inst_elem = tree.find("Run/Instrument")
    dict['instrument'] = inst_elem.text

    flowcell_layout_elem = tree.find("Run/FlowcellLayout")
    num_lanes = int(flowcell_layout_elem.attrib['LaneCount'])
   
    dict['lanes'] = { }
    for i in (range(1,num_lanes+1)):
        dict['lanes'][str(i)] = { }
        

    for lane in dict['lanes']:
        dict['lanes'][lane]['reads'] = [ ]        

    dict['total_cycles'] = 0
    dict['read_cycles']  = 0
    dict['index_cycles'] = 0

    dict['index_reads'] = set()

    for read in reads:

        read_num        = int(read.attrib['Number'])
        num_cycles      = int(read.attrib['NumCycles'])
        is_indexed_read = read.attrib['IsIndexedRead']
   
        dict['total_cycles'] += num_cycles

        for lane in dict['lanes']:
            dict['lanes'][lane]['reads'].append({'read_number':read_num})
 
        if (is_indexed_read == 'N'):
            dict['read_cycles'] += num_cycles 
        elif(is_indexed_read == 'Y'):
            dict['index_cycles'] += num_cycles
            dict['index_reads'].add(read_num)
        else:
            sys.exit('unknown IsIndexedRead value %s' % (is_indexed_read)) 

    read_num_lookup = [ ]
    last_read_num   = None

    for read in reads:
        for i in range(int(read.attrib['NumCycles'])):
            read_num_lookup.append(int(read.attrib['Number']))
        if (last_read_num is not None):
            if (last_read_num >= int(read.attrib['Number'])):
                sys.exit('reads out of order') 
        last_read_num = int(read.attrib['Number']) 

    dict['read_num_lookup'] = read_num_lookup

    return dict

#def transmit_stats(runinfo):
#
#    json_stats = json.dumps(runinfo)
#
#    conn = SQSConnection(aws_access_key, aws_secret_key)
#
#    q = conn.get_queue('lims_publish_stats')
#
#    m = Message()
#       
#    m.set_body(json_stats)
#
#    status = q.write(m)
#  
#    if (not status):
#        sys.exit('SQS failed to accept message')    

#def transmit_plots(folder_name, png_files):
#
#    conn = S3Connection(aws_access_key, aws_secret_key)
#
#    for b in conn.get_all_buckets():
#
#        if (b.name == 'gtac_lims_ivc_plots'):
#
#            for png_file in png_files:
#
#                png_file_basename = os.path.basename(png_file)
#
#                k = Key(b)
#                k.key = folder_name + '_' + png_file_basename 
#                k.set_contents_from_filename(png_file)

def parse_interop(path, read_cycle_lookup):

    interop_path = os.path.join(path, 'InterOp')
   
    tm = parse_tilemetrics(interop_path)
    qm = parse_qualitymetrics(interop_path, read_cycle_lookup)
    em = parse_errormetrics(interop_path, read_cycle_lookup)
    cm = parse_correctedintensitymetrics(interop_path)
     
    return dict(tm.items() + qm.items() + em.items() + cm.items())

def parse_tilemetrics(path):

    fn = 'TileMetricsOut.bin'
    fn = os.path.join(path, fn)

    dict = { }

    a = BitString(bytes=open(fn, 'rb').read())

    file_version       = a.read('uintle:8')
    record_byte_length = a.read('uintle:8')
    
    if not (2 <= file_version <= 3):
        sys.exit('tile metrics file version %d not supported' % (file_version)) 

    header_size = 16
    tile_area   = 0

    if (file_version == 3):
       tile_area = a.read('floatle:32')
       header_size = 48

    cluster_density = { } 
    percent_aligned = { }
    num_cluster     = { }
    num_cluster_pf  = { }

    # header_size bits of header, 8 bits to a byte 
    for i in range(0,((a.len - header_size) / (record_byte_length * 8))):
 
        lane_number  = a.read('uintle:16')

        tile_number = 0

        if (file_version <= 2):
            tile_number = a.read('uintle:16')
        else:
            tile_number = a.read('uintle:32')

        if (file_version <= 2):
        
            metric_code  = a.read('uintle:16')
            metric_value = a.read('floatle:32')

            # only grab desired metrics
            if (metric_code == 100):
                if (not cluster_density.has_key(lane_number)):
                    cluster_density[lane_number] = [ ]
                cluster_density[lane_number].append(metric_value)
            elif (metric_code == 102):
                if (not num_cluster.has_key(lane_number)):
                    num_cluster[lane_number] = 0
                num_cluster[lane_number]+= metric_value
            elif (metric_code == 103):
                if (not num_cluster_pf.has_key(lane_number)):
                    num_cluster_pf[lane_number] = 0
                num_cluster_pf[lane_number] += metric_value
            elif (metric_code == 300):
                if (not percent_aligned.has_key(lane_number)):
                    percent_aligned[lane_number] = { } 
                if (not percent_aligned[lane_number].has_key(1)):
                    percent_aligned[lane_number][1] = [ ]
                percent_aligned[lane_number][1].append(metric_value)
            elif (metric_code == 301):
                if (not percent_aligned.has_key(lane_number)):
                    percent_aligned[lane_number] = { } 
                if (not percent_aligned[lane_number].has_key(2)):
                    percent_aligned[lane_number][2] = [ ]
                percent_aligned[lane_number][2].append(metric_value)
            elif (metric_code == 302):
                if (not percent_aligned.has_key(lane_number)):
                    percent_aligned[lane_number] = { } 
                if (not percent_aligned[lane_number].has_key(3)):
                    percent_aligned[lane_number][3] = [ ]
                percent_aligned[lane_number][3].append(metric_value)
        else:

            metric_code  = a.read('bytes:1')

            if (metric_code == 't'):
                
                cluster_count    = a.read('floatle:32')            
                pf_cluster_count = a.read('floatle:32')
                
                calc_cluster_density = cluster_count / tile_area

                if (not cluster_density.has_key(lane_number)):
                    cluster_density[lane_number] = [ ]

                cluster_density[lane_number].append(calc_cluster_density)

                if (not num_cluster.has_key(lane_number)):
                    num_cluster[lane_number] = 0
                num_cluster[lane_number]+= cluster_count 

                if (not num_cluster_pf.has_key(lane_number)):
                    num_cluster_pf[lane_number] = 0
                num_cluster_pf[lane_number] += pf_cluster_count 

            elif (metric_code == 'r'):

                read_number = a.read('uintle:32') 
                pct_aligned = a.read('floatle:32')
                
                if (not percent_aligned.has_key(lane_number)):
                    percent_aligned[lane_number] = { }
                if (not percent_aligned[lane_number].has_key(read_number)):
                    percent_aligned[lane_number][read_number] = [ ]
                percent_aligned[lane_number][read_number].append(pct_aligned)

            elif (metric_code == '\00'):
                empty = a.read('bytes:8')
 
            else:
                sys.exit('illegal metric_code %s encountered' % (metric_code))
                
 
    dict['cluster_density']        = { }
    dict['cluster_density_stddev'] = { }

    for lane in cluster_density:

        narray = numpy.array(cluster_density[lane])

        dict['cluster_density'][lane]        = round(numpy.mean(narray) / 1000, 2)
        dict['cluster_density_stddev'][lane] = round(numpy.std(narray) / 1000, 2)

    dict['num_cluster']            = { }
    dict['num_cluster_pf']         = { }
    dict['pct_cluster_pf']         = { }

    for lane in num_cluster:

        dict['num_cluster'][lane]    = num_cluster[lane]
        dict['num_cluster_pf'][lane] = num_cluster_pf[lane]
        dict['pct_cluster_pf'][lane] = round((float(num_cluster_pf[lane]) / float(num_cluster[lane])) * 100, 2)
   
    dict['percent_aligned'] = { }

    for lane in percent_aligned:
  
        for read in percent_aligned[lane]:

            narray = numpy.array(percent_aligned[lane][read])
            
            if (not dict['percent_aligned'].has_key(lane)):
                dict['percent_aligned'][lane] = { }
            if (not dict['percent_aligned'][lane].has_key(read)):
                dict['percent_aligned'][lane][read] = { }
            dict['percent_aligned'][lane][read] = round(numpy.mean(narray), 2)
    
    return dict

def parse_qualitymetrics(path, read_cycle_lookup):

    fn = 'QMetricsOut.bin'
    fn = os.path.join(path, fn)

    dict = { }

    a = BitString(bytes=open(fn, 'rb').read())

    file_version       = a.read('uintle:8')
    record_byte_length = a.read('uintle:8')

    header_len = 16

    if not (4 <= file_version <= 7):
        sys.exit('Q metrics file version %d not supported' % (file_version))

    qscore_bin_flag = 0

    if ( 5 <= file_version):
        qscore_bin_flag = a.read('uintle:8')
        header_len      += 8
   
    num_qscore_bins  = None
    remapped_qscores = None 

    if (qscore_bin_flag == 1):
        num_qscore_bins   = a.read('uintle:8')
        remapped_qscores  = [0] * num_qscore_bins
        # lower boundary of quality score bins
        for i in range(num_qscore_bins):
            a.read('uintle:8')
        # upper boundary of quality score bins
        for i in range(num_qscore_bins):
            a.read('uintle:8')
        # remapped scores of quality score bins
        for i in range(num_qscore_bins):
            remapped_qscores[i] = a.read('uintle:8')

    total_clusters = { }
    q30_clusters   = { }
    percent_q30    = { }
   
    for i in range(0,((a.len - header_len) / (record_byte_length * 8))):

        lane_number  = a.read('uintle:16')

        tile_number = 0
 
        if (file_version <= 6):
            tile_number = a.read('uintle:16')
        else:
            tile_number = a.read('uintle:32')

        cycle_number = a.read('uintle:16')
        
        read_number  = read_cycle_lookup[cycle_number - 1] 

        num_buckets = 50
    
        if (qscore_bin_flag):
            num_buckets = num_qscore_bins
        
        for i in range(num_buckets):

           count = a.read('uintle:32')

           if (not total_clusters.has_key(lane_number)):
               total_clusters[lane_number] = { }
           
           if (not total_clusters[lane_number].has_key(read_number)):
               total_clusters[lane_number][read_number] = 0

           total_clusters[lane_number][read_number] += count
           
           bucket_qscore = None 

           if (qscore_bin_flag == 0):
               bucket_qscore = i + 1 
           else:
               bucket_qscore = remapped_qscores[i]

           if (bucket_qscore >= 30):

               if (not q30_clusters.has_key(lane_number)):
                   q30_clusters[lane_number] = { }
 
               if (not q30_clusters[lane_number].has_key(read_number)):
                   q30_clusters[lane_number][read_number] = 0

               q30_clusters[lane_number][read_number] += count   

    dict['percent_q30'] = { }

    for lane in (total_clusters):

        dict['percent_q30'][lane] = { }

        for read in (total_clusters[lane]):
            percent_q30 = None
            if (q30_clusters[lane][read] == 0):
                percent_q30 = 0.00
            else:
                percent_q30 = (q30_clusters[lane][read] / float(total_clusters[lane][read])) * 100
            dict['percent_q30'][lane][read] = round(percent_q30, 2)
    
    return dict

def parse_errormetrics(path, read_cycle_lookup):

    fn = 'ErrorMetricsOut.bin'
    fn = os.path.join(path, fn)

    dict = { }

    if (not os.path.exists(fn)):
        return dict

    a = BitString(bytes=open(fn, 'rb').read())

    file_version       = a.read('uintle:8')
    record_byte_length = a.read('uintle:8')

    if not (3 <= file_version <= 4):
        sys.exit('error metrics file version %d not supported' % (file_version)) 

    error_rates = { }

    # 16 bits of header, 8 bits to a byte 
    for i in range(0,((a.len - 16) / (record_byte_length * 8))):
 
        lane_number       = a.read('uintle:16')

        tile_number = 0
 
        if (file_version < 4):
            tile_number       = a.read('uintle:16')
        else:
            tile_number       = a.read('uintle:32')

        cycle_number      = a.read('uintle:16')
        error_rate        = a.read('floatle:32')

        if (file_version < 4):

            num_perfect_reads = a.read('uintle:32') 
            reads_one_error   = a.read('uintle:32') 
            reads_two_error   = a.read('uintle:32') 
            reads_three_error = a.read('uintle:32') 
            reads_four_error  = a.read('uintle:32') 
    
        read_number = read_cycle_lookup[cycle_number - 1]
        
        if (not error_rates.has_key(lane_number)):
            error_rates[lane_number] = { } 
        
        if (not error_rates[lane_number].has_key(read_number)):
            error_rates[lane_number][read_number] = [ ]

        error_rates[lane_number][read_number].append(error_rate)            

    dict['error_rate']        = { }
    dict['error_rate_stddev'] = { }

    for lane in error_rates:
       
        dict['error_rate'][lane]        = { }
        dict['error_rate_stddev'][lane] = { }
 
        for read in error_rates[lane]:

            narray = numpy.array(error_rates[lane][read])
            
            dict['error_rate'][lane][read]        = round(numpy.mean(narray), 2)
            dict['error_rate_stddev'][lane][read] = round(numpy.std(narray), 2) 
   
    return dict

def parse_correctedintensitymetrics(path):

    fn = 'CorrectedIntMetricsOut.bin'
    fn = os.path.join(path, fn)

    dict       = { }
    num_called = { }

    a = BitString(bytes=open(fn, 'rb').read())

    file_version       = a.read('uintle:8')
    record_byte_length = a.read('uintle:8')

    if not (2 <= file_version <= 4):
        sys.exit('extraction metrics file version %d not supported' % (file_version))

    # 16 bits of header, 8 bits to a byte 
    for i in range(0,((a.len - 16) / (record_byte_length * 8))):

        tile_number = 0

        lane_number  = a.read('uintle:16')
        
        if (file_version <= 3):
            tile_number  = a.read('uintle:16')
        else:
            tile_number  = a.read('uintle:32')

        cycle_number = a.read('uintle:16') 
  
        if (file_version == 2):

            avg_int             = a.read('uintle:16')
            avg_corrected_int_A = a.read('uintle:16')
            avg_corrected_int_C = a.read('uintle:16')
            avg_corrected_int_G = a.read('uintle:16')
            avg_corrected_int_T = a.read('uintle:16')

        if (file_version <= 3):

            avg_corrected_int_called_A = a.read('uintle:16')
            avg_corrected_int_called_C = a.read('uintle:16')
            avg_corrected_int_called_G = a.read('uintle:16')
            avg_corrected_int_called_T = a.read('uintle:16')

        num_N = 0
        num_A = 0
        num_C = 0
        num_G = 0
        num_T = 0
        
        if (file_version == 2):

            num_N = a.read('floatle:32')
            num_A = a.read('floatle:32')
            num_C = a.read('floatle:32')
            num_G = a.read('floatle:32')
            num_T = a.read('floatle:32')
    
        if (file_version >= 3):

            num_N = float(a.read('uintle:32'))
            num_A = float(a.read('uintle:32'))
            num_C = float(a.read('uintle:32'))
            num_G = float(a.read('uintle:32'))
            num_T = float(a.read('uintle:32'))

        if (file_version == 2):

            snr = a.read('floatle:32')

        if (not num_called.has_key(lane_number)):
            num_called[lane_number] = { }

        if (not num_called[lane_number].has_key(cycle_number)):
            num_called[lane_number][cycle_number] = { }
            num_called[lane_number][cycle_number]['N'] = 0 
            num_called[lane_number][cycle_number]['A'] = 0 
            num_called[lane_number][cycle_number]['C'] = 0  
            num_called[lane_number][cycle_number]['G'] = 0  
            num_called[lane_number][cycle_number]['T'] = 0 
     
        num_called[lane_number][cycle_number]['N'] += num_N
        num_called[lane_number][cycle_number]['A'] += num_A
        num_called[lane_number][cycle_number]['C'] += num_C
        num_called[lane_number][cycle_number]['G'] += num_G
        num_called[lane_number][cycle_number]['T'] += num_T

    pct_called = { } 
     
    for l in num_called:
 
        if (not pct_called.has_key(l)):
            pct_called[l] = { }

        for c in num_called[l]:

            if (not pct_called[l].has_key(c)):
                pct_called[l][c] = {}

            total_calls = num_called[l][c]['N'] + num_called[l][c]['A'] + num_called[l][c]['C'] + num_called[l][c]['G'] + num_called[l][c]['T']
           
            try: 
                pct_called[l][c]['A'] = (num_called[l][c]['A'] / total_calls) * 100
                pct_called[l][c]['C'] = (num_called[l][c]['C'] / total_calls) * 100
                pct_called[l][c]['G'] = (num_called[l][c]['G'] / total_calls) * 100
                pct_called[l][c]['T'] = (num_called[l][c]['T'] / total_calls) * 100
            except ZeroDivisionError:
                pct_called[l][c]['A'] = 0.0
                pct_called[l][c]['C'] = 0.0
                pct_called[l][c]['G'] = 0.0
                pct_called[l][c]['T'] = 0.0

    dict['percent_called'] = pct_called
            
    return dict

#def ivc_plot(percent_called, output_fn):
    
#    da = { 'cycle': [ ], 'percent': [ ] }
#    dc = { 'cycle': [ ], 'percent': [ ] }
#    dg = { 'cycle': [ ], 'percent': [ ] }
#    dt = { 'cycle': [ ], 'percent': [ ] }

#    for cycle in percent_called:

#        da['cycle'].append(cycle)
#        dc['cycle'].append(cycle)
#        dg['cycle'].append(cycle)
#        dt['cycle'].append(cycle)
#        da['percent'].append(percent_called[cycle]['A'])
#        dc['percent'].append(percent_called[cycle]['C'])
#        dg['percent'].append(percent_called[cycle]['G'])
#        dt['percent'].append(percent_called[cycle]['T'])

#    dfa = pandas.DataFrame(data=da)
#    dfc = pandas.DataFrame(data=dc)
#    dfg = pandas.DataFrame(data=dg)
#    dft = pandas.DataFrame(data=dt)
    
#    plt = ggplot(dfa, aes(x='cycle', y='percent')) + geom_line(color='green') + geom_line(data=dfc, color='red') +  geom_line(data=dfg, color='blue') +  geom_line(data=dft, color='purple')

#    ggsave(plot=plt, filename=output_fn)


#aws_access_key = environ['AWS_ACCESS_KEY'] 
#aws_secret_key = environ['AWS_SECRET_KEY'] 

parser = argparse.ArgumentParser(description='Illumina LIMS Run Publisher')

parser.add_argument(
                    '--run-folder',
                    '-run-folder',
                    '-r',
                    required=True,
                    help='full path to run folder',
                   )

args = parser.parse_args()

run_folder_full_path = args.run_folder

if (not valid_runfolder(run_folder_full_path)):
    sys.exit('invalid run folder %s' % (run_folder_full_path))

os.umask(0007)

runinfo    = lookup_runinfo(run_folder_full_path)
run_number = runinfo['run_number']
instrument = runinfo['instrument']

basecalls_full_path = os.path.join(run_folder_full_path, 'Data', 'Intensities', 'BaseCalls')

interop = parse_interop(run_folder_full_path, runinfo['read_num_lookup'])

plot_files = [ ]

for ln in runinfo['lanes']:

    lane_number = ln

    #plot_fn = os.path.join(run_folder_full_path, 's_%d_percent_all.png' % (int(lane_number)))

    #ivc_plot(interop['percent_called'][int(lane_number)], plot_fn)

    #plot_files.append(plot_fn)

    runinfo['lanes'][lane_number]['cluster_density']        = interop['cluster_density'][int(lane_number)]
    runinfo['lanes'][lane_number]['cluster_density_stddev'] = interop['cluster_density_stddev'][int(lane_number)]

    num_reads = interop['num_cluster_pf'][int(lane_number)]
    pct_pf    = interop['pct_cluster_pf'][int(lane_number)]

    for rd in runinfo['lanes'][lane_number]['reads']:

        read_number = rd['read_number']

        rd['num_cycles'] = runinfo['read_num_lookup'].count(read_number)

        rd['num_reads'] = num_reads 
        rd['percent_pass_filter'] = pct_pf
        rd['percent_q30'] = interop['percent_q30'][int(lane_number)][read_number]

        if (read_number in runinfo['index_reads']):
            rd['error_rate']        = 0
            rd['error_rate_stddev'] = 0
            rd['percent_aligned']   = 0
        else:
            if (interop.has_key('error_rate')):
                if (interop['error_rate'].has_key(int(lane_number))):
                    if (interop['error_rate'][int(lane_number)].has_key(read_number)): 
                        rd['error_rate'] = interop['error_rate'][int(lane_number)][read_number]
                    else:
                        rd['error_rate'] = 0
                else:
                    rd['error_rate'] = 0
            else:
                rd['error_rate'] = 0
 
            if (interop.has_key('error_rate_stddev')):
                if (interop['error_rate_stddev'].has_key((int(lane_number)))):
                    if (interop['error_rate_stddev'][int(lane_number)].has_key(read_number)):
                        rd['error_rate_stddev'] = interop['error_rate_stddev'][int(lane_number)][read_number]
                    else:
                        rd['error_rate_stddev'] = 0    
                else:
                    rd['error_rate_stddev'] = 0
            else:
                rd['error_rate_stddev'] = 0

            if (interop['percent_aligned'].has_key((int(lane_number)))):
                if (interop['percent_aligned'][int(lane_number)].has_key(read_number)):
                    rd['percent_aligned'] = interop['percent_aligned'][int(lane_number)][read_number]
                else:
                    rd['percent_aligned'] = 0
            else:
                rd['percent_aligned'] = 0       

del runinfo['index_reads'] 

print runinfo

#transmit_stats(runinfo)
#transmit_plots(runinfo['run_id'], plot_files)
