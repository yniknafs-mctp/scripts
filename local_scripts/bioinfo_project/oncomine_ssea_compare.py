'''
Created on Jan 30, 2014

@author: mkiyer
'''
import logging
import argparse
import os
import sys
import operator
import re
import random
import collections
import xlrd

# local imports
from ssea.lib.config import Config 
from ssea.lib.countdata import BigCountMatrix
from ssea.lib.base import Result, JobStatus, computerize_name
import ssea.lib.cfisher as fisher

def parse_results(filename):
    with open(filename, 'r') as fp:
        for line in fp:
            result = Result.from_json(line.strip())
            yield result

def get_direction(name):
    if name.find('Over-expressed') != -1:
        return 'up'
    elif name.find('Under-expressed') != -1:
        return 'dn'
    
def get_study(name):
    p = re.compile(r'.+expressed\s\((.+)\)')
    m = p.match(name)
    return m.groups()[0]

class SSEAGeneSet(object):
    def __init__(self):
        pass
    @staticmethod
    def parse(input_path, fracs, bm, meta):
        results_file = os.path.join(input_path, Config.RESULTS_JSON_FILE)
        ss_compname = os.path.basename(input_path)
        ss_up = {}
        ss_dn = {}
        i = 0
        # collapse to genes by taking best performing transcript
        for res in parse_results(results_file):
            # logging
            i += 1
            if (i % 10000) == 0:
                logging.debug('Parsed %d results' % (i))
            # skip if zero
            if res.ss_frac == 0:
                continue       
            # skip if transcript id not found
            transcript_id = bm.rownames[res.t_id]
            if transcript_id not in meta:
                continue
            # convert to ref gene symbols
            symbols = meta[transcript_id]
            ss = ss_up if res.ss_frac > 0 else ss_dn
            for s in symbols:
                ss[s] = (transcript_id, res.nes, res.ss_frac)
        logging.debug("Found %d up and %d dn symbols for set '%s'" % (len(ss_up), len(ss_dn), ss_compname))
        # sort by ss frac
        ss_up = sorted([(k,) + v for k,v in ss_up.iteritems()], key=operator.itemgetter(3), reverse=True)
        ss_dn = sorted([(k,) + v for k,v in ss_dn.iteritems()], key=operator.itemgetter(3), reverse=False)
        # null list is all gene symbols
        null_list = set([x[0] for x in ss_up])
        null_list.update(set([x[0] for x in ss_dn]))
        for frac in fracs:
            for direction in ('up', 'dn'):
                ss = ss_up if direction == 'up' else ss_dn
                n = int(len(ss) * frac)
                symbs = [x[0] for x in ss[:n]]
                self = SSEAGeneSet()
                self.name = '%s_%s_%s' % (ss_compname, str(frac), direction)
                self.direction = direction
                self.study = 'SSEA'
                self.null_name = 'mitranscriptome'
                self.null_list = null_list
                self.gene_list = set(symbs)
                self.missing = set()
                yield self

class DESeqGeneSet(object):
    def __init__(self):
        pass
    @staticmethod
    def parse(input_path, fracs, meta):
        results_file = open(input_path)
        ss_compname = 'deseq_'+os.path.basename(input_path)
        ss_up = collections.defaultdict(lambda: (1, 1, 1))
        ss_dn = collections.defaultdict(lambda: (1, 1, 1))
        i = 0
        # collapse to genes by taking best performing transcript
        deseq_header = results_file.next().strip().split('\t')
        for line in results_file:
            # logging
            i += 1
            if (i % 10000) == 0:
                logging.debug('Parsed %d results' % (i))
            # skip if zero
            line = line.strip().split('\t')
            fold_change = line[deseq_header.index('foldChange')]
            if fold_change=='0' or fold_change=="NA":
                continue
            log_fold = float(line[deseq_header.index('log2FoldChange')])
            pval = float(line[deseq_header.index('pval')])
            transcript_id = line[deseq_header.index('id')].replace('\"','')
            if fold_change == 0:
                continue       
            # skip if transcript id not found
            if transcript_id not in meta:
                continue
            # convert to ref gene symbols
            symbols = meta[transcript_id]
            ss = ss_up if log_fold < 0 else ss_dn
            for s in symbols:
                if pval < ss[s][1]:
                    ss[s] = (transcript_id, pval, log_fold)
        logging.debug("Found %d up and %d dn symbols for set '%s'" % (len(ss_up), len(ss_dn), ss_compname))
        # sort by ss frac
        ss_up = sorted([(k,) + v for k,v in ss_up.iteritems()], key=operator.itemgetter(2), reverse=False)
        ss_dn = sorted([(k,) + v for k,v in ss_dn.iteritems()], key=operator.itemgetter(2), reverse=False)
        # null list is all gene symbols
        null_list = set([x[0] for x in ss_up])
        null_list.update(set([x[0] for x in ss_dn]))
        for frac in fracs:
            for direction in ('up', 'dn'):
                ss = ss_up if direction == 'up' else ss_dn
                n = int(len(ss) * frac)
                symbs = [x[0] for x in ss[:n]]
                self = DESeqGeneSet()
                self.name = '%s_%s_%s' % (ss_compname, str(frac), direction)
                self.direction = direction
                self.study = 'DESeq'
                self.null_name = 'mitranscriptome'
                self.null_list = null_list
                self.gene_list = set(symbs)
                self.missing = set()
                yield self

class OncomineGeneSet(object):
    def __init__(self):
        pass
    @staticmethod
    def parse(filename, platforms):
        wkbook = xlrd.open_workbook(filename)
        for i in xrange(wkbook.nsheets):
            wksheet = wkbook.sheet_by_index(i)
            # count sets
            nsets = 0
            fields = wksheet.row_values(1)
            for j in xrange(0, len(fields), 2):
                if fields[j] == 'Concept Information':
                    nsets += 1
            # names
            studies = []
            names = []
            directions = []
            fields = wksheet.row_values(3)
            for j in xrange(0, len(fields), 2):
                assert fields[j] == 'Name:'
                name = fields[j+1]
                study = get_study(name)
                direction = get_direction(name)
                names.append(name)
                studies.append(study)
                directions.append(direction)
            # null list
            null_names = []
            null_lists = []
            fields = wksheet.row_values(6)
            for j in xrange(0, len(fields), 2):
                assert fields[j] == 'Null List:'
                null_name = fields[j+1].strip()
                gene_list = platforms[null_name]
                null_names.append(null_name)
                null_lists.append(gene_list)
            # gene lists
            for j in xrange(0, wksheet.ncols, 2):
                n = (j >> 1)
                null_list = null_lists[n]
                null_name = null_names[n]
                gene_list = set()
                missing = set()
                values = wksheet.col_values(j, 9, wksheet.nrows)
                for v in values:
                    if not v:
                        continue
                    if v in null_list:
                        gene_list.add(v)
                    else:
                        missing.add(v)
                self = OncomineGeneSet()
                self.name = names[n]
                self.direction = directions[n]
                self.study = studies[n]
                self.null_name = null_name
                self.null_list = set(null_list)
                self.gene_list = set(gene_list)
                self.missing = missing
                yield self

def ssea_vs_oncomine(aset, gset, anull, gnull):
    # subset both sets
    nullset = anull.intersection(gnull)
    aset = aset.intersection(nullset)
    gset = gset.intersection(nullset)
    # find overlap
    tp = len(aset.intersection(gset))
    fp = len(aset.difference(gset))
    fn = len(gset.difference(aset))
    tn = len(nullset.difference(aset).difference(gset))
    # fisher exact test (one-sided hypothesis that LE is enricheD)
    fisher_p_value = fisher.pvalue(tp, fp, fn, tn).right_tail
    if (fp == 0) or (fn == 0):
        odds_ratio = float('inf')
    else:
        odds_ratio = (tp * tn) / float(fp * fn)
    return tp, fp, fn, tn, odds_ratio, fisher_p_value

def oncomine_vs_oncomine(a, b, anull, bnull):
    # subset both sets
    null = anull.intersection(bnull)
    a = a.intersection(null)
    b = b.intersection(null)
    # find overlap
    tp = len(a.intersection(b))
    fp = len(a.difference(b))
    fn = len(b.difference(a))
    tn = len(null.difference(a).difference(b))
    # fisher exact test (one-sided hypothesis that LE is enricheD)
    fisher_p_value = fisher.pvalue(tp, fp, fn, tn).right_tail
    if (fp == 0) or (fn == 0):
        odds_ratio = float('inf')
    else:
        odds_ratio = (tp * tn) / float(fp * fn)
    return tp, fp, fn, tn, odds_ratio, fisher_p_value

def compare_two_sets(a, b, anull, bnull):
    # subset both sets
    null = anull.intersection(bnull)
    a = a.intersection(null)
    b = b.intersection(null)
    # find overlap
    tp = len(a.intersection(b))
    fp = len(a.difference(b))
    fn = len(b.difference(a))
    tn = len(null.difference(a).difference(b))
    # fisher exact test (one-sided hypothesis that LE is enricheD)
    #fisher_p_value = fisher.pvalue(tp, fp, fn, tn).right_tail
    # fisher exact test (two-tailed hypothesis)
    fisher_p_value = fisher.pvalue(tp, fp, fn, tn).two_tail
    if (fp == 0) or (fn == 0):
        odds_ratio = float('inf')
    else:
        odds_ratio = (tp * tn) / float(fp * fn)
    return tp, fp, fn, tn, odds_ratio, fisher_p_value


def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("--platforms", dest="platforms_dir", default=None)
    parser.add_argument("--oncomine", dest="oncomine_files", action='append')
    parser.add_argument('--matrix', dest='matrix_dir', default=None)
    parser.add_argument('--meta', dest='metadata_file', default=None)
    parser.add_argument("--ssea", dest="ssea_paths", action='append')
    parser.add_argument("--deseq", dest="deseq_paths", action='append')
    parser.add_argument('-f', '--frac', dest='fracs', type=float, action='append')
    parser.add_argument('--output-sets', dest='output_sets_prefix', default=None)
    args = parser.parse_args()
    # get args
    matrix_dir = args.matrix_dir
    ssea_paths = args.ssea_paths
    deseq_paths = args.deseq_paths
    metadata_file = args.metadata_file
    oncomine_files = args.oncomine_files
    platforms_dir = args.platforms_dir
    output_sets_prefix = args.output_sets_prefix
    fracs = args.fracs
    if fracs is None:
        fracs = [0.01, 0.05, 0.1]
    else:
        fracs = sorted(set(fracs))
    # check args
    if oncomine_files is None:
        parser.error("Specify oncomine gene set file(s) using --oncomine")
    else:
        for filename in oncomine_files:
            if not os.path.exists(filename):
                parser.error("Oncomine file '%s' not found" % (filename))
    if platforms_dir is None:
        parser.error("Specify platforms dir using --platforms")
    elif not os.path.exists(platforms_dir):
        parser.error("Platforms directory '%s' not found" % (platforms_dir))
    platforms_file = os.path.join(platforms_dir, 'platforms.txt')
    if not os.path.exists(platforms_file):
        parser.error('Platforms file "%s" not found' % (platforms_file))
    if ssea_paths is None:
        parser.error('Specify SSEA result paths(s) using --ssea')
    else:
        for path in ssea_paths:
            if not os.path.exists(path):
                parser.error('SSEA path "%s" not found' % (path))
    if metadata_file is None:
        parser.error("Specify transcript metadata using --meta")
    elif not os.path.exists(metadata_file):
        parser.error("Transcript metadata file '%s' not found" % (metadata_file))
    if matrix_dir is None:
        parser.error('Specify SSEA matrix file using --matrix')
    elif not os.path.exists(matrix_dir):
        parser.error("Matrix directory '%s' not found" % (matrix_dir))
    # parse platforms
    logging.info('Parsing platforms')
    platforms = {}
    with open(platforms_file) as f:
        for line in f:
            name, filename = line.strip().split('\t')
            path = os.path.join(platforms_dir, os.path.basename(filename))
            genelist = set([line.strip() for line in open(path)])
            platforms[name] = genelist
            logging.debug('Platform name="%s" has %d genes' % (name, len(genelist)))
    # setup oncomine
    gsets = []
    logging.info('Parsing Oncomine gene sets')
    for filename in oncomine_files:    
        for gset in OncomineGeneSet.parse(filename, platforms):
            logging.debug('Name=%s Size=%d Null Size=%d Missing=%d' % (gset.name, len(gset.gene_list), len(gset.null_list), len(gset.missing)))
            gsets.append(gset)
    # setup metadata
    assembly_gene_set = set()
    logging.debug('Parsing transcript metadata')
    meta = {}
    with open(metadata_file) as f:
        header_fields = f.next().strip().split()
        ref_gene_name_col = header_fields.index('ref_gene_name')
        category_col = header_fields.index('category')
        for line in f:
            fields = line.strip().split('\t')
            t_id = fields[0]
            category = fields[category_col]
            if ((category == 'same_strand') or
                (category == 'read_through')):
                ref_gene_names = fields[ref_gene_name_col].split(',')
                meta[t_id] = ref_gene_names
                assembly_gene_set.update(ref_gene_names)
    logging.debug('Found metadata for %d transcripts' % (len(meta)))
    logging.debug('Assembly has %d gene mappings' % (len(assembly_gene_set)))
    # parse ssea sets
    bm = BigCountMatrix.open(matrix_dir)
    for ssea_path in ssea_paths:
        logging.debug('Parsing SSEA path %s' % (ssea_path))
        for gset in SSEAGeneSet.parse(ssea_path, fracs, bm, meta):
            logging.debug('Name=%s Size=%d Null Size=%d Missing=%d' % (gset.name, len(gset.gene_list), len(gset.null_list), len(gset.missing)))
            gsets.append(gset)
    # parse deseq sets
    for deseq_path in deseq_paths:
        logging.debug('Parsing DESeq path %s' % (deseq_path))
        for gset in DESeqGeneSet.parse(deseq_path, fracs, meta):
            logging.debug('Name=%s Size=%d Null Size=%d Missing=%d' % (gset.name, len(gset.gene_list), len(gset.null_list), len(gset.missing)))
            gsets.append(gset)
    # output gene sets    
    if output_sets_prefix is not None:
        logging.info('Writing gene sets to files using prefix %s' % (output_sets_prefix))
        for gset in gsets:
            name = computerize_name(gset.name)
            filename = os.path.join(output_sets_prefix + '.' + name)
            with open(filename, 'w') as f:
                for g in gset.gene_list:
                    print >>f, g
    # compare gene sets
    header_fields = ['study1', 'dir1', 'name1', 'study2', 'dir2', 'name2', 
                     'tp', 'fp', 'fn', 'tn', 'oddsratio',' pvalue']
    print '\t'.join(header_fields)
    for i in xrange(len(gsets)):
        for j in xrange(len(gsets)):
            a = gsets[i]
            b = gsets[j]
            fields = [a.study, a.direction, a.name, b.study, b.direction, b.name]
            res = compare_two_sets(a.gene_list, b.gene_list, a.null_list, b.null_list)
            fields.extend(map(str, res))
            print '\t'.join(fields)
    bm.close()
    return 0



    # if no assembly path specified compare oncomine to itself
    if (input_path is None) or (metadata_file is None):
        logging.debug('Comparing Oncomine gene sets to each other')
        header_fields = ['study1', 'dir1', 'name1', 'study2', 'dir2', 'name2', 
                         'tp', 'fp', 'fn', 'tn', 'oddsratio',' pvalue']
        print '\t'.join(header_fields)
        for i in xrange(len(gsets)-1):
            for j in xrange(i + 1, len(gsets)):
                a = gsets[i]
                b = gsets[j]
                if a.study == b.study:
                    continue
                if a.direction != b.direction:
                    continue
                fields = [a.study, a.direction, a.name, b.study, b.direction, b.name]
                fields.extend(map(str, oncomine_vs_oncomine(a.gene_list, b.gene_list, a.null_list, b.null_list)))
                print '\t'.join(fields)        
        return 1

    # parse ssea results
    bm = BigCountMatrix.open(matrix_dir)
    logging.debug('Parsing path %s' % (input_path))
    results_file = os.path.join(input_path, Config.RESULTS_JSON_FILE)
    ss_compname = os.path.basename(input_path)
    ss_up = {}
    ss_dn = {}
    i = 0    
    # collapse to genes by taking best performing transcript
    for res in parse_results(results_file):
        # logging
        i += 1
        if (i % 10000) == 0:
            logging.debug('Parsed %d results' % (i))
        # skip if zero
        if res.ss_frac == 0:
            continue       
        # skip if transcript id not found
        transcript_id = bm.rownames[res.t_id]
        if transcript_id not in meta:
            continue
        # convert to ref gene symbols
        symbols = meta[transcript_id]
        ss = ss_up if res.ss_frac > 0 else ss_dn
        for s in symbols:
            ss[s] = (transcript_id, res.nes, res.ss_frac)
    logging.debug("Found %d up and %d dn symbols for set '%s'" % (len(ss_up), len(ss_dn), ss_compname))
    # sort by ss frac
    ss_up = sorted([(k,) + v for k,v in ss_up.iteritems()], key=operator.itemgetter(3), reverse=True)
    ss_dn = sorted([(k,) + v for k,v in ss_dn.iteritems()], key=operator.itemgetter(3), reverse=False)
    # now do comparison with gene set
    if output_sets_file is not None:
        output_sets_fileh = open(output_sets_file, 'w')
    header_fields = ['study1', 'dir1', 'frac', 'type', 'name1', 'study2', 'dir2', 
                     'name2', 'tp', 'fp', 'fn', 'tn', 'oddsratio',' pvalue']
    print '\t'.join(header_fields)
    for frac in fracs:
        # output gene sets to separate file for later use        
        if output_sets_file is not None:
            for direction in ('up', 'dn'):
                ss = ss_up if direction == 'up' else ss_dn
                n = int(len(ss) * frac)
                null = [x[0] for x in ss]
                symbs = null[:n]
                fields = [ss_compname, str(frac), direction]
                fields.extend(symbs)
                print >>output_sets_fileh, '\t'.join(symbs)
        # compare to oncomine sets
        for gset in gsets:
            ss = ss_up if gset.direction == 'up' else ss_dn
            n = int(len(ss) * frac)
            null = [x[0] for x in ss]
            symbs = null[:n]
            fields = ['ssea_' + ss_compname, gset.direction, str(frac), 'observed', 
                      'SSEA %s top %f %s' % (ss_compname, frac, gset.direction),
                      gset.study, gset.direction, gset.name]
            fields.extend(map(str, ssea_vs_oncomine(set(symbs), gset.gene_list, set(null), gset.null_list)))
            print '\t'.join(fields)
            # negative control is shuffled gene ranks
            for i in xrange(rand_iters):
                random.shuffle(null)
                symbs = null[:n]
                fields = ['ssea_' + ss_compname, gset.direction, str(frac), 'random',
                          'Random %s top %f %s' % (ss_compname, frac, gset.direction),
                          gset.study, gset.direction, gset.name]
                fields.extend(map(str, ssea_vs_oncomine(set(symbs), gset.gene_list, set(null), gset.null_list)))
                print '\t'.join(fields)            
    if output_sets_file is not None:
        output_sets_fileh.close()
    bm.close()

if __name__ == '__main__':
    sys.exit(main())