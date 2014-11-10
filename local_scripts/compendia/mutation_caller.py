'''
Created on Dec 4, 2013

@author: yniknafs

testestest
'''
import os 
import sys
import argparse
import logging
import pymongo


PEAK_TYPE_US_BLOCK = 1

def db_connect(name, host):
    logging.info('connecting to %s database on mongo server: %s' % (name, host))
    client = pymongo.MongoClient(host)
    db = client[name]
    row_metadata = db['metadata']
    col_metadata = db['samples']
    sample_sets = db['sample_sets']
    config = db['config']
    reports = db['reports']
    merged = db['merged']
    colls = {'row_meta':row_metadata, 'col_meta':col_metadata, 'ss':sample_sets, 'config':config, 'reports':reports, 'merged': merged}
    return colls

def tip(db, host, out, tip_file):
    colls = db_connect(db, host)
    sample_sets = colls['ss']
    reports = colls['reports']
    row_meta = colls['row_meta']
    
    
    file_out = open(out, 'w')
    #print header for the output file 
    headers = ['t_name',
               'category',
               'peak',
               'peak_type',
               'nes',
               'fdr']
    print >> file_out, '\t'.join(map(str, headers)) 
    
    
    tip_file = open(tip_file,'r')
    tip_list = tip_file.readlines()
    tip_len = len(tip_list)
    i = 0
    for line in tip_list: 
        line = line.strip()
        t_name, peak_id = line.split('\t')
        peak_type = peak_id.split('_')[PEAK_TYPE_US_BLOCK]
        #query the row_meta collection to obtain transcript category
        meta = row_meta.find_one({'name':t_name}, ['category', '_id'])
        if meta == None:
            logging.debug('Skipping transcript %s, not included in SSEA run' % t_name)
            continue
        category = meta['category']
        t_id = meta['_id']
        
        #will need the ss_id for getting the NES and FDR from the reports collection
        ss_id_cur = sample_sets.find_one({'name':peak_id}, {'_id':1})
        if ss_id_cur == None:
            logging.debug('Skipping sample set %s, not included in SSEA run' % peak_id)
            continue
        ss_id = ss_id_cur['_id']
        
        report = reports.find_one({'t_id': t_id, 'ss_id':ss_id}, ['nes', 'ss_fdr_q_value'])
        nes = report['nes']
        fdr = report['ss_fdr_q_value']
        
        list_out = [t_name,
                    category, 
                    peak_id, 
                    peak_type,
                    nes,
                    fdr]
        
        print >> file_out, '\t'.join(map(str,list_out))
        
    file_out.close()

def peak(db, host, out, peak, tip_file):
    colls = db_connect(db, host)
    sample_sets = colls['ss']
    merged = colls['merged']
    row_meta = colls['row_meta']
    
    
    file_out = open(out, 'w')
    #print header for the output file 
    headers = ['t_id',
               'category',
               'peak',
               'chr',
               'peak_status',
               'nes',
               'fdr']
    print >> file_out, '\t'.join(map(str, headers)) 
    
    tip_set = set([line.strip().split('\t')[0] for line in open(tip_file)])
    peak_id = sample_sets.find_one({'name':peak})['_id']
    loc = peak.split('_')[2]
    loc_rep = loc.replace('q','p')
    chr = loc_rep.split('p')[0]
    
    tot = merged.find({'ss_id': int(peak_id)}).count()
    i = 0
    for line in merged.find({'ss_id': int(peak_id)}):
        if (i % 1000) == 0:
            logging.debug('Finished %d/%d' % (i, tot))
            
        i +=1
        t_id = line['t_id']
        
        #query the row_meta collection to obtain transcript category
        meta = row_meta.find_one({'_id':int(t_id)}, ['category', 'locus', 'name'])
        t_name = meta['name']
        locus = meta['locus']
        category = meta['category']
        if t_name in tip_set:
            peak_status = 'in_peak'
        else: 
            chr_str = locus.split(':')[0]
            chr_int = chr_str.split('r')[1]
            if chr_int == chr:
                peak_status = 'same_chr'
            else: 
                peak_status = 'diff_chr'
                
        
        
        
        #will need the ss_id for getting the NES and FDR from the reports collection
        report = merged.find_one({'t_id': t_id, 'ss_id':peak_id}, ['nes', 'ss_fdr_q_value'])
        nes = report['nes']
        fdr = report['ss_fdr_q_value']
        
        list_out = [t_name,
                    category,
                    peak,
                    chr,
                    peak_status,
                    nes,
                    fdr]
        
        print >> file_out, '\t'.join(map(str,list_out))
        
    file_out.close()
    
    
    
def main():
    '''Command line options.'''    
    parser = argparse.ArgumentParser()
    # Add command line parameters
    parser.add_argument("tip_file", 
                        help="CNV transcripts in peaks file")
    parser.add_argument('db',
                        help = 'name of mongodb database to connect to')
    parser.add_argument('out',
                        help = 'output file')
    parser.add_argument("--host", dest = 'host',
                        default = 'localhost:27017',
                        help = 'name of mongodb server to connect to')
    parser.add_argument("--tip", dest="tip", 
                        action="store_true", default=False, 
                        help="process for all transcripts in peaks")
    parser.add_argument("--peak", dest="peak", 
                        default=False, 
                        help="process for one peak")
    args = parser.parse_args()
    level = logging.DEBUG
    logging.basicConfig(level=level,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")


        
    if args.tip != False:
        tip(args.db, args.host, args.out, args.tip_file)
    
    if args.peak != False: 
        peak(args.db, args.host, args.out, args.peak, args.tip_file)
        
        
        
        
        
        
        
        
        


if __name__ == '__main__':
    sys.exit(main())