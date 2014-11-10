'''
Created on Feb 13, 2014
@author yniknafs
'''


import os
import sys
import argparse
import logging
import numpy as np
import subprocess
import shutil
import math

header = [
             'ss_name',
             'ss_cat',
             'ss_type',
             'frac_cutoff',
             'es_sign',
             'fdr_cutoff',
             'metric', 
             'threshold',
             'category',
             'count'
            ]

def chunks(l, n):
    """ Yield successive n-sized chunks from l.
    """
    n = int(n)
    for i in xrange(0, len(l), n):
        yield l[i:i+n]

    
def main():
    # parse command line
    parser = argparse.ArgumentParser()
#     parser.add_argument("metadata_json")
    parser.add_argument("sets_file")
#     parser.add_argument("root_dir")
    parser.add_argument("-p", dest = 'num_proc')
    args = parser.parse_args()
    
    num_sets = len(open(args.sets_file).readlines())
    num_chunks = math.ceil(num_sets/float(args.num_proc))
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    if not os.path.isdir('STATS_TMPS'):
        os.mkdir('STATS_TMPS')
    prefix = 'STATS_TMPS'    
    sets_list = [x.split() for x in  open(args.sets_file).readlines()]
    i = 0
    chunk_file_list = []
    chunk_list =  list(chunks(sets_list, num_chunks))
    for chunk in chunk_list:
        filename = 'chunk' + str(i)
        i+=1
        file_path = os.path.join(prefix, filename)
        chunk_file_list.append(file_path)
        with open(file_path, 'w') as f: 
            for set in chunk:
                print >>f, set[0]
    i = 0
    out_file_list = []
    procs = []
    for fchunk in chunk_file_list: 
        out_file = 'out' + str(i)
        i+=1
        out_file_path = os.path.join(prefix, out_file)
        out_file_list.append(out_file_path)
        pargs = [
                 'python',
                 '/mctp/users/yniknafs/scripts/workspace/local_scripts/compendia/ssea_set_stats/set_stats_for_parallel.py',
                 '/mctp/projects/mitranscriptome/ssea/metadata_tucp.json',
                 fchunk,
                 '/mctp/projects/mitranscriptome/ssea/results',
                 out_file_path
                 ]
        
        p = subprocess.Popen(pargs)
        procs.append(p)
    
    for proc in procs: 
        proc.wait()
    
    with open('set_stats_all.txt', 'w') as outfile:
        outfile.write('\t'.join(header)+'\n')
        for fname in out_file_list:
            with open(fname) as infile:
                outfile.write(infile.read())
    
    shutil.rmtree(prefix)
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())
