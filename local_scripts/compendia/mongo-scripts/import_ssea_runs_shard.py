'''
Created on Dec 17, 2013

@author: yniknafs

When given a directory containing ssea runs, will import all 
ssea runs within the directory

'''

import sys
import argparse
import os
import logging
import glob
import subprocess
import ssea
import shutil 

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("runs_dir")
    parser.add_argument("-m", dest = "matrix_dir",
                        default = '/mctp/projects/ssea/isoform_count_matrix_v7/')
    parser.add_argument("-n", "--name", dest = 'name',
                        default = 'compendia',
                        help = 'name for ssea run (will be name of database)')
    parser.add_argument("--host", dest = 'host',
                        default = '172.20.54.118:27017',
                        help = 'name of mongodb server to connect to')
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    _ssea_path = ssea.__path__[0]
    _merge_path = os.path.join(_ssea_path, 'utils/mongo_import_ssea_v3.py')
    
    dirs = glob.glob(args.runs_dir + '/*')
    for dir in dirs:
        logging.debug("Checking dir: \'%s\'" % dir)
        tsv = os.path.join(dir, 'sample_sets.tsv')
        if not os.path.exists(tsv):
            print tsv
            logging.debug('Directory: \'%s\'  is not an SSEA run directory' % args.runs_dir)
            continue
        imported = os.path.join(dir, 'imported')
#         if not os.path.exists(imported):
#             os.mkdir(imported)
        ls = glob.glob(dir + '/*')
        for item in ls: 
            done = os.path.join(item, 'job.done')
            importing = os.path.join(item, 'job.importing')
            imported = os.path.join(item, 'job.imported')
            if item != tsv and item != imported and os.path.exists(done) and (not os.path.exists(importing)):
                with open(importing, "w"): pass
                import_args = ['python', _merge_path, item, args.matrix_dir, '-n', args.name, '--host', args.host]
                p1 = subprocess.call(import_args)
#                 shutil.move(item, imported)
                os.rename(importing, imported)
    
    


if __name__ == '__main__':
    sys.exit(main())