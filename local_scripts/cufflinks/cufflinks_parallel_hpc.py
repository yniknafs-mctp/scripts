    
'''
Created on May 10, 2014
@author yniknafs
'''


import os
import sys
import argparse
import logging
import subprocess
import glob
import time

COMPUTE_GENE_EXPR_SCRIPT = '/mctp/users/yniknafs/scripts/workspace/assemblyline/assemblyline/pipeline/compute_gene_expression.py'

def worker(args):
    (sample,tot, i) = args 
    logging.info('Running cufflinks for %s; worker %d/%d' % (os.path.basename(sample), i, tot))
    sh_args = [
               'sh',
               os.path.join(sample, 'cufflinks.sh')
              ]
    with open(os.path.join(sample, 'log.txt'), 'w') as f:
        subprocess.call(sh_args, stderr=f)
    
def sleeper(max_procs):
    job_count = int(subprocess.check_output('qstat -u yniknafs| grep " R " | wc -l', shell=True))
    while job_count >= max_procs: 
        time.sleep(30)
        job_count = int(subprocess.check_output('qstat -u yniknafs| grep " R " | wc -l', shell=True)) 
        
    

def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("config_xml")
    parser.add_argument("library_table")
    parser.add_argument("gtf_file")
    parser.add_argument("out_dir")
    parser.add_argument("-p", dest = 'proc',
                    default = 4,
                    help = 'number of processors to use')
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
     
    args.proc = int(args.proc)
    #run the script to make the PBS scripts
    script_gen_args = [
                       'python',
                       COMPUTE_GENE_EXPR_SCRIPT,
                       '-o',
                       args.out_dir, 
                       args.config_xml,
                       args.library_table,
                       args.gtf_file,
                       '--mode',
                       'cufflinks' 
                       ]
     
    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)
     
    logging.debug('Generating PBS scripts...')
    subprocess.call(script_gen_args)
     
    sample_dirs = glob.glob(os.path.join(args.out_dir, '*'))
     
    #kick off jobs on cluster run cufflinks
    logging.debug('Running %d jobs on cluster using %d cores' % (len(sample_dirs), args.proc))
    tasks = []
    i=0
    for sample in sample_dirs:
        i+=1
        script = os.path.join(sample, 'cufflinks.sh')
        sleeper(args.proc)
        logging.debug('Running job %d/%d: %s' % (i, len(sample_dirs), os.path.basename(sample)))
        subprocess.call('qsub ' + script, shell=True)
        
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())    

