'''
Created on Jul 17, 2014
@author yniknafs
'''


import os
import sys
import argparse
import logging
import collections


def add_meta(args):
    new_meta = args.new_meta
    if new_meta is None:
        logging.error('Please specify new meta file using "-n" or "--new_meta"')
        return 1
    fields = args.fields
    if fields is None:
        logging.error('Please specify fields using "-f" or "--fields"')
        return 1
    lib_field = args.lib_field
    if lib_field is None:
        logging.error('Please specify library field using "-lf" or "--lib_field"')
        return 1
    
    #make dict for new fields 
    fileh = open(new_meta)
    new_header = fileh.next().strip().split('\t')
    
    #make sure provided fields are in file
    field_break = 0
    for field in fields: 
        if field not in new_header:
            field_break+=1
            logging.error("Field: %s not in provided library metadata file" % field)
    if field_break >0:
        return 1
    
    #make sure library field is in provided file
    if lib_field not in new_header:
        logging.error("Provided library field (%s) is not in meta file (the new one)" % lib_field)
        return 1
    fields_dict = collections.defaultdict(lambda: collections.defaultdict(lambda: ''))
    for line in fileh:
        line = line.strip().split('\t')
        lib = line[new_header.index(lib_field)]
        for field in fields:
            field_val = line[new_header.index(field)]
            fields_dict[field][lib] = field_val
    
    old_meta_fh = open(args.library_meta_file)
    old_meta_header = old_meta_fh.next().strip().split('\t')
    for field in fields: 
        old_meta_header.append(field)
    print '\t'.join(old_meta_header)
    for line in old_meta_fh:
        line = line.strip().split('\t')
        lib = line[old_meta_header.index('library_id')]
        for field in fields:
            field_val = fields_dict[field][lib]
            line.append(field_val)
        print '\t'.join(map(str, line))
               
    
    
    
    
def main():
    # parse command line
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()
    # add to meta
    subparser = subparsers.add_parser('add')
    subparser.add_argument("library_meta_file")
    subparser.add_argument("-n", "--new_meta", dest="new_meta", default=None)
    subparser.add_argument("-lf", "--lib_field", dest="lib_field", default=None)
    subparser.add_argument("-f", "--fields", dest="fields", action = 'append', default=None)
    subparser.set_defaults(func=add_meta)
    
    args = parser.parse_args()
    args.func(args)
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())
