import argparse
import logging
import os
import sys
import csv
import time
import collections
import bisect
import operator
import subprocess
import glob


def main():

    
    logging.basicConfig(level=logging.DEBUG, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    parser = argparse.ArgumentParser()
    parser.add_argument("lesion_file")
    parser.add_argument("bed_file")
    parser.add_argument("output_dir")
    args = parser.parse_args()
    
    a = os.path.abspath(args.lesion_file)
    b = os.path.dirname(a)
    c = os.path.basename(b)
    
    
    def makebeds():
        logging.info("Creating bed files for intersection")
        with open(args.lesion_file, 'r') as f:
            labels = f.readline().replace(" ", "").strip().split('\t')
    
            label_dict = collections.defaultdict(lambda: [])
            for x in xrange(len(labels)): 
                label_dict[labels[x]]=x 
            
            narrow_peaks = []
            wide_peaks = []
            
            reader=csv.reader(f,delimiter='\t')
            for labels in reader: 
                labels = map(lambda s:s.replace(" ",""), labels)
                if labels[0].endswith("CNvalues") == False:
                    np_loc = labels[label_dict['PeakLimits']].replace('(',':').replace(')', '').replace('-',':').split(':')
                    wp_loc = labels[label_dict['WidePeakLimits']].replace('(',':').replace(')', '').replace('-',':').split(':')
                    narrow_peaks.append((np_loc[0], np_loc[1], np_loc[2], labels[label_dict['UniqueName']],labels[label_dict['qvalues']]))
                    wide_peaks.append((wp_loc[0], wp_loc[1], wp_loc[2], labels[label_dict['UniqueName']], labels[label_dict['qvalues']]))
            outfile_n = open('narrow.tmp.bed', 'w')
            outfile_w = open('wide.tmp.bed', 'w')
            for x in narrow_peaks: 
                outfile_n.write('\t'.join(map(str, x)))
                outfile_n.write('\n')
            for x in wide_peaks: 
                outfile_w.write('\t'.join(map(str, x)))
                outfile_w.write('\n')
            outfile_n.close()
            outfile_w.close() 

    peak_file_narrow = open(args.output_dir + "/"+ c + "_narrowpeaks_by_transcript.txt", 'w')
    peak_file_wide = open(args.output_dir + "/"+ c + "_widepeaks_by_transcript.txt", 'w')
    
    def bed_intersect():
        logging.info("Creating intersection")
        p1 = subprocess.check_output(["bedtools", 'intersect', '-a', 'narrow.tmp.bed', '-b', args.bed_file, '-loj'])
        p2 = subprocess.check_output(["bedtools", 'intersect', '-a', 'wide.tmp.bed', '-b', args.bed_file, '-loj'])
        
        for line in p1.splitlines():
            items = line.split('\t')
            if items[8] != '.':
                a, b = (items[8].split('|')[1], items[8].split('|')[0])
                outlist = [a, b, items[3], items[0], items[1], items[2]]
                peak_file_narrow.write('\t'.join(map(str, outlist)))
                peak_file_narrow.write('\n')
        peak_file_narrow.close()
        
        for line in p2.splitlines():
            items = line.split('\t')
            if items[8] != '.':
                a, b = (items[8].split('|')[1], items[8].split('|')[0])
                outlist = [a, b, items[3], items[0], items[1], items[2]]
                peak_file_wide.write('\t'.join(map(str, outlist)))
                peak_file_wide.write('\n')
        peak_file_wide.close()
    
    
    
    
    makebeds()
    bed_intersect()
    #p2 = subprocess.Popen(['rm','*.tmp.bed'])            
    #print os.path.abspath(os.path.curdir)
    #print os.listdir(".")
    os.remove("narrow.tmp.bed")
    os.remove("wide.tmp.bed")
    
    logging.info("Finished")   
    


    return 0



if __name__ == '__main__':

    sys.exit(main())