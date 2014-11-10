'''
Created on Aug 29, 2014
@author yniknafs
'''


import os
import sys
import argparse
import logging
import subprocess
import glob
import time
import re

GSEA_CORR_PATH = '/mctp/users/yniknafs/scripts/workspace/local_scripts/gsea/gsea_corr.py'
SLEEP_TIME = 30
PROC_PERCENT = 5
USER = 'yniknafs'

TISSUES = ['breast','cervical','colorectal',
           'gbm','head_neck','liver',
           'lgg','luad','lusc','pancreatic',
           'prostate','kich','kirc',
           'kirp','stomach','uterine',
           'thyroid', 'melanoma', 'ovarian']
TISSUES = set(TISSUES)

def hpc_load():
    used_nodes = 0
    #check hpc
    HOST="mctp-hpc-login3"
    COMMAND='qstat -u %s| grep " R " | wc -l' % USER
    ssh = subprocess.Popen(["ssh", "%s" % HOST, COMMAND],
                           shell=False,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE)
    qstat = int(ssh.stdout.read().strip())
    used_nodes += qstat
    COMMAND='qstat -u %s| grep " Q " | wc -l' % USER
    ssh = subprocess.Popen(["ssh", "%s" % HOST, COMMAND],
                           shell=False,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE)
    qstat = int(ssh.stdout.read().strip())
    used_nodes += qstat
    
#     for line in qstat: 
#         fields = re.split('\s+', line)
#         if fields[7] == 'R':
#             used_nodes+=int(fields[4])
    return used_nodes


def pb8_load():
    #check pb8
    used_nodes =0
    HOST="pathbio-8"
    COMMAND="ps -u %s|grep python" % USER
    ssh = subprocess.Popen(["ssh", "%s" % HOST, COMMAND],
                           shell=False,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE)
    ps = ssh.stdout.readlines()
    used_nodes +=len(ps)
#     COMMAND="ps -u %s|grep java" % USER
#     ssh = subprocess.Popen(["ssh", "%s" % HOST, COMMAND],
#                            shell=False,
#                            stdout=subprocess.PIPE,
#                            stderr=subprocess.PIPE)
#     ps = ssh.stdout.readlines()
#     used_nodes += len(ps)
    return used_nodes


def pb9_load():
    #check pb9
    used_nodes =0
    HOST="pathbio-9"
    COMMAND="ps -u %s|grep python" % USER
    ssh = subprocess.Popen(["ssh", "%s" % HOST, COMMAND],
                           shell=False,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE)
    ps = ssh.stdout.readlines()
    used_nodes +=len(ps)
#     COMMAND="ps -u %s|grep java" % USER
#     ssh = subprocess.Popen(["ssh", "%s" % HOST, COMMAND],
#                            shell=False,
#                            stdout=subprocess.PIPE,
#                            stderr=subprocess.PIPE)
#     ps = ssh.stdout.readlines()
#     used_nodes += len(ps)
    return used_nodes



def sleeper(hpc_procs, pb8_procs, pb9_procs):
    hpc_load_test = hpc_load()
    pb8_load_test = pb8_load()
    pb9_load_test = pb9_load()
    while hpc_load()>= hpc_procs and pb8_load()>=pb8_procs and pb9_load()>=pb9_procs: 
#         logging.debug('sleeping...')
#         hpc_load_test = hpc_load()
#         pb8_load_test = pb8_load()
#         pb9_load_test = pb9_load()
#         logging.debug(hpc_load_test)
#         logging.debug(pb8_load_test)
#         logging.debug(pb9_load_test)
        time.sleep(15)
    free_agents = []
    
    if hpc_load() < hpc_procs: 
        free_agents.append('hpc')
    if pb8_load() < pb8_procs: 
        free_agents.append('pb8')
    if pb9_load() < pb9_procs: 
        free_agents.append('pb9')
    return free_agents

def get_pbs_header(job_name,
                   node_processors=1,
                   node_memory=4096,
                   pbs_script_lines=None, 
                   working_dir=None, 
                   deps=None,
                   stdout_filename=None,
                   stderr_filename=None,
                   modules=None):
    '''
    job_name: string name of job
    node_processors: number of cores available per node
    node_memory: amount of memory per node (MB)
    pbs_script_lines: list of PBS directives to be added to the script
    working_dir: the "working directory" of the job (allows scripts to access files using relative pathnames)
    deps: 'None' if no dependencies, or a python list of job ids
    stdout_filename: string filename for storing stdout
    stderr_filename: string filename for storing stderr    
    '''    
    if pbs_script_lines is None:
        pbs_script_lines = []
    if isinstance(deps, basestring):
        deps = [deps]
    # add PBS parameters
    lines = ["#PBS -N %s" % job_name]
    lines.extend(pbs_script_lines)
    lines.append('#PBS -q batch')
    lines.append('#PBS -m abe')
#     lines.append('#PBS -M yniknafs@med.umich.edu')
    lines.append('#PBS -V')
    lines.append('#PBS -l nodes=1:ppn=1,mem=4000mb,walltime=240:00:00')
    
    
    if stdout_filename is None:
        stdout_filename = "/dev/null"
    lines.append("#PBS -o %s" % (stdout_filename))
    if stderr_filename is None:
        stderr_filename = "/dev/null"        
    lines.append("#PBS -e %s" % (stderr_filename))
    if deps is not None:
        lines.append("#PBS -W depend=afterok:%s" % (":".join([d for d in deps])))    
    if working_dir is not None: 
        lines.append("cd %s" % (working_dir))
    lines.append('source /mctp/wkgrps/bioinfo/sw/rhel6/init.sh')
    lines.append('module purge')
    lines.append('module add epd')
    if modules is not None: 
        for module in modules: 
            lines.append('module add %s' % module)
    return lines


def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("associations")
    parser.add_argument("trans_meta")
    parser.add_argument("lib_meta")
    parser.add_argument("matrix")
    parser.add_argument("--perms", dest = "perms")
    parser.add_argument("--gmt", dest = "gmt")
    parser.add_argument("-o", dest = "out", default = None)
    parser.add_argument("--expr_mat_dir", dest = "expr_mat_dir", default = None)
    parser.add_argument("--hpc", dest = 'hpc',
                    default = 0,
                    help = 'number of HPC processors to use')
    parser.add_argument("--pb8", dest = 'pb8',
                    default = 0,
                    help = 'number of pb8 processors to use')
    parser.add_argument("--pb9", dest = 'pb9',
                    default = 0,
                    help = 'number of pb9 processors to use')
    
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    args.hpc = int(args.hpc)
    args.pb8 = int(args.pb8)
    args.pb9 = int(args.pb9)
    args.gmt = os.path.abspath(args.gmt)
    args.matrix = os.path.abspath(args.matrix)
    args.trans_meta = os.path.abspath(args.trans_meta)
    args.lib_meta = os.path.abspath(args.lib_meta)
    args.out = os.path.abspath(args.out)
    logging.debug(args.out)
    if not os.path.exists(args.out):
        os.mkdir(args.out)
    
    def hpc_job(tid, assoc, i, out_dir):
        #kick off job on hpc
        HOST="mctp-hpc-login3"
        gsea_args = [
                'python',
                GSEA_CORR_PATH,
                tid,
                args.trans_meta,
                args.lib_meta,
                assoc,
                args.matrix,
                '--gmt',
                args.gmt,
                '--perms',
                args.perms,
                '--job_id',
                i,
                '-o',
                out_dir,
                '--expr_mat_dir',
                args.expr_mat_dir
                ]
        
        script_file = os.path.join(out_dir, 'job_%s_%s_%s.sh' % (i, tid, assoc))
        out_file = os.path.join(out_dir, 'job_%s_%s_%s.out' % (i, tid, assoc))
        err_file = os.path.join(out_dir, 'job_%s_%s_%s.err' % (i, tid, assoc))
        shell_commands = get_pbs_header(job_name='%s_%s' % (tid, assoc),
                                        stdout_filename=os.path.abspath(out_file),
                                        stderr_filename=os.path.abspath(err_file))
        gsea_command = ' '.join(map(str,gsea_args))
        shell_commands.append(gsea_command)
        
        f = open(script_file, "w")
        for command in shell_commands:
            print >>f, command
        f.close()
        COMMAND = 'qsub %s' % os.path.abspath(script_file)
#         logging.debug(COMMAND)
        p = subprocess.Popen(["ssh", "%s" % HOST, COMMAND],
                               shell=False,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
#         logging.debug(p.stdout.read())
#         logging.debug(p.stderr.read())
    def pb8_job(tid, assoc, i, out_dir):
        #kick off job on hpc
        HOST="pathbio-8"
        
        gsea_args = [
                'python',
                GSEA_CORR_PATH,
                tid,
                args.trans_meta,
                args.lib_meta,
                assoc,
                args.matrix,
                '--gmt',
                args.gmt,
                '--perms',
                args.perms,
                '--job_id',
                i,
                '-o',
                out_dir,
                '--expr_mat_dir',
                args.expr_mat_dir
                ]
        
        COMMAND= ' '.join(map(str,gsea_args))
        subprocess.Popen(["ssh", "%s" % HOST, COMMAND],
                               shell=False,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
        
    def pb9_job(tid, assoc, i, out_dir):
        #kick off job on hpc
        HOST="pathbio-9"
        gsea_args = [
                'python',
                GSEA_CORR_PATH,
                tid,
                args.trans_meta,
                args.lib_meta,
                assoc,
                args.matrix,
                '--gmt',
                args.gmt,
                '--perms',
                args.perms,
                '--job_id',
                i,
                '-o',
                out_dir,
                '--expr_mat_dir',
                args.expr_mat_dir
                ]
        
        COMMAND= ' '.join(map(str,gsea_args))
        subprocess.Popen(["ssh", "%s" % HOST, COMMAND],
                               shell=False,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    
    
    tot = 0
    fileh = open(args.associations)
    for line in fileh: 
        line = line.strip().split('\t')
        tid = line[0]
        associations = line[1].split(',')
        for assoc in associations: 
            if assoc not in TISSUES:
                continue
            tot+=1
    logging.debug('Running %d jobs total...' % tot)
    logging.debug('Using %d cores on HPC' % args.hpc)
    logging.debug('Using %d cores on PB8' % args.pb8)
    logging.debug('Using %d cores on PB9' % args.pb9)
    i=0
    
    #read associations to get list of jobs to do
    fileh = open(args.associations)
    for line in fileh: 
        line = line.strip().split('\t')
        tid = line[0]
        associations = line[1].split(',')
        for assoc in associations: 
            if assoc not in TISSUES:
#                 logging.error("Tissue: %s does not have enough samples for correlation analysis" % assoc)
                continue
            i+=1
            JOB_RUNNING = os.path.join(args.out, '%s_%s.running' % (tid, assoc))
            JOB_DONE = os.path.join(args.out, '%s_%s.done' % (tid, assoc))
            if os.path.exists(JOB_RUNNING):
                logging.info("SKIPPING %s because it is running..." % os.path.basename(JOB_RUNNING).split('.')[0])
                continue
            elif os.path.exists(JOB_DONE):
                logging.info("SKIPPING %s because it is done..." % os.path.basename(JOB_DONE).split('.')[0])
                continue
                
            free_agents = sleeper(args.hpc, args.pb8, args.pb9)
            if 'hpc' in free_agents:
                logging.debug('Running job %d/%d on HPC: %s' % (i, tot, tid))
                hpc_job(tid, assoc, i, args.out)
            elif 'pb8' in free_agents:
                logging.debug('Running job %d/%d on PB8: %s' % (i, tot, tid))
                pb8_job(tid, assoc, i, args.out)
            elif 'pb9' in free_agents:
                logging.debug('Running job %d/%d on PB9: %s' % (i, tot, tid))
                pb9_job(tid, assoc, i, args.out)
             
            with open(JOB_RUNNING, 'w'):
                pass
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())

