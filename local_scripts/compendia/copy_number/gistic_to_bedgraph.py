'''
Created on Apr 7, 2014
@author yniknafs
'''


import os
import sys
import argparse
import logging
import numpy as np

def track_info(name,desc,cnv):
    if cnv=='Amp': 
        color='255,0,0'
    if cnv=='Del':
        color='0,0,255'
    info = [
              'track',
              'type=bedGraph',
              'name='+name,
              'description='+desc,
              'graphType=bar',
              'visibility=full',
              'color='+color,            
           ]
    return ' '.join(info)

def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("gistic")
    parser.add_argument("out_prefix")
    
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    gistic_fh = open(args.gistic)
    gistic_head = gistic_fh.next().strip().split('\t')
    amp_lines = []
    del_lines = []
    for line in gistic_fh: 
                line = line.strip().split('\t')
                type = line[gistic_head.index('Type')]
                chr = line[gistic_head.index('Chromosome')].strip()
                chr = 'chr'+chr
                start = line[gistic_head.index('Start')]
                end = line[gistic_head.index('End')]
                score = line[gistic_head.index('G-score')]
                lineo = [chr, start, end, score]
                if type == "Amp":
                    amp_lines.append(lineo)
                if type == 'Del':
                    del_lines.append(lineo)   
    
    new_amps = []
    print track_info('pancan_cnv_amp', 'pancan_cnv_amp', 'Amp')
    for x in xrange(len(amp_lines)): 
        chr, start, end, score = amp_lines[x]
        new_amps.append(amp_lines[x])
        if x == len(amp_lines)-1: 
            continue
        chrx, startx, endx, scorex = amp_lines[x+1]
        if chr != chrx: 
            continue
        if end == startx:
            continue
        scoren = np.mean([float(score), float(scorex)])
        startn = end
        endn = startx
        linen = [chr, startn, endn, scoren]
        new_amps.append(linen)
    
    for line in new_amps:
        print '\t'.join(map(str,line))
    
    return 0    
        
    
    
    with open(args.out_prefix+'.amp.bedGraph','w') as famp:
        print >>famp, track_info('pancan_cnv_amp', 'pancan_cnv_amp', 'Amp')
        with open(args.out_prefix+'.del.bedGraph','w') as fdel:
            print >>fdel, track_info('pancan_cnv_del', 'pancan_cnv_del', 'Del')
            for line in gistic_fh: 
                line = line.strip().split('\t')
                type = line[gistic_head.index('Type')]
                chr = line[gistic_head.index('Chromosome')].strip()
                chr = 'chr'+chr
                start = line[gistic_head.index('Start')]
                end = line[gistic_head.index('End')]
                score = line[gistic_head.index('G-score')]
                lineo = [chr, start, end, score]
                if type == "Amp":
                    print >>famp, '\t'.join(lineo)
                if type == 'Del':
                    print >>fdel, '\t'.join(lineo)    
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())
