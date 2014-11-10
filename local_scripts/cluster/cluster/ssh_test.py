'''
Created on Aug 21, 2014
@author yniknafs
'''


import os
import sys
import argparse
import logging
import paramiko
import subprocess




    
def main():
    # parse command line
    parser = argparse.ArgumentParser()
    
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    
#     logging.info('Testing paramiko')
#     ssh = paramiko.SSHClient()
#     ssh.connect('pathbio-8', username='yniknafs', password='Rd^5CmaF5')
#     ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
#     stin,stout,sterr = ssh.exec_command('echo hi')
    

    
    HOST="pathbio-8"
    # Ports are handled in ~/.ssh/config since we use OpenSSH
    COMMAND="sleep 15"
    
    ssh = subprocess.Popen(["ssh", "%s" % HOST, COMMAND],
                           shell=False,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE)
    result = ssh.stdout.readlines()
    if result == []:
        error = ssh.stderr.readlines()
        print >>sys.stderr, "ERROR: %s" % error
    else:
        print result
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())
