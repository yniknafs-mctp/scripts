'''
Created on Jan 21, 2014
@author yniknafs
'''


import os
import sys
import argparse
import logging
import numpy as np
import subprocess




    
def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("pathbio")
    
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")

    assembly = ([
                'python', '/mctp/users/yniknafs/scripts/workspace/local_scripts/GWAS/GWAS_locus_shuffle_OR_parallel.py',
                '--gtf', '/mctp/projects/mitranscriptome/annotation/assembly.gtf',
                '--shuffs', '100',
                '-p', '60'
                ], 'gwas_assembly.txt')
    
    assembly_exon = ([
                'python', '/mctp/users/yniknafs/scripts/workspace/local_scripts/GWAS/GWAS_locus_shuffle_OR_parallel.py',
                '--gtf', '/mctp/projects/mitranscriptome/annotation/assembly.gtf',
                '--shuffs', '100',
                '-p', '60',
                '--exon',
                ], 'gwas_assembly_exon.txt')
    
    lncrna = ([
               'python', '/mctp/users/yniknafs/scripts/workspace/local_scripts/GWAS/GWAS_locus_shuffle_OR_parallel.py',
                '--gtf', '/mctp/projects/mitranscriptome/gwas/gtfs/tcat/lncrna.gtf',
                '--shuffs', '100',
                '-p', '60'
               ], 'gwas_lncrna.txt')
    
    lncrna_exon = ([
               'python', '/mctp/users/yniknafs/scripts/workspace/local_scripts/GWAS/GWAS_locus_shuffle_OR_parallel.py',
                '--gtf', '/mctp/projects/mitranscriptome/gwas/gtfs/tcat/lncrna.gtf',
                '--shuffs', '100',
                '--exon',
                '-p', '60'
               ], 'gwas_lncrna_exon.txt')
    
    tucp = ([
               'python', '/mctp/users/yniknafs/scripts/workspace/local_scripts/GWAS/GWAS_locus_shuffle_OR_parallel.py',
                '--gtf', '/mctp/projects/mitranscriptome/gwas/gtfs/tcat/tucp.gtf',
                '--shuffs', '100',
                '-p', '60'
               ], 'gwas_tucp.txt')
    
    tucp_exon = ([
               'python', '/mctp/users/yniknafs/scripts/workspace/local_scripts/GWAS/GWAS_locus_shuffle_OR_parallel.py',
                '--gtf', '/mctp/projects/mitranscriptome/gwas/gtfs/tcat/tucp.gtf',
                '--shuffs', '100',
                '--exon',
                '-p', '60'
               ], 'gwas_tucp_exon.txt')
    protein_coding = ([
               'python', '/mctp/users/yniknafs/scripts/workspace/local_scripts/GWAS/GWAS_locus_shuffle_OR_parallel.py',
                '--gtf', '/mctp/projects/mitranscriptome/gwas/gtfs/tcat/protein_coding.gtf',
                '--shuffs', '100',
                '-p', '60'
               ], 'gwas_protein_coding.txt')
    
    protein_coding_exon = ([
               'python', '/mctp/users/yniknafs/scripts/workspace/local_scripts/GWAS/GWAS_locus_shuffle_OR_parallel.py',
                '--gtf', '/mctp/projects/mitranscriptome/gwas/gtfs/tcat/protein_coding.gtf',
                '--shuffs', '100',
                '--exon',
                '-p', '60'
               ], 'gwas_protein_coding_exon.txt')
    pseudogene = ([
               'python', '/mctp/users/yniknafs/scripts/workspace/local_scripts/GWAS/GWAS_locus_shuffle_OR_parallel.py',
                '--gtf', '/mctp/projects/mitranscriptome/gwas/gtfs/tcat/pseudogene.gtf',
                '--shuffs', '100',
                '-p', '60'
               ], 'gwas_pseudogene.txt')
    
    pseudogene_exon = ([
               'python', '/mctp/users/yniknafs/scripts/workspace/local_scripts/GWAS/GWAS_locus_shuffle_OR_parallel.py',
                '--gtf', '/mctp/projects/mitranscriptome/gwas/gtfs/tcat/pseudogene.gtf',
                '--shuffs', '100',
                '--exon',
                '-p', '60'
               ], 'gwas_pseudogene_exon.txt')
    mixed_read_through = ([
               'python', '/mctp/users/yniknafs/scripts/workspace/local_scripts/GWAS/GWAS_locus_shuffle_OR_parallel.py',
                '--gtf', '/mctp/projects/mitranscriptome/gwas/gtfs/tcat/mixed_read_through.gtf',
                '--shuffs', '100',
                '-p', '60'
               ], 'gwas_mixed_read_through.txt')
    
    mixed_read_through_exon = ([
               'python', '/mctp/users/yniknafs/scripts/workspace/local_scripts/GWAS/GWAS_locus_shuffle_OR_parallel.py',
                '--gtf', '/mctp/projects/mitranscriptome/gwas/gtfs/tcat/mixed_read_through.gtf',
                '--shuffs', '100',
                '--exon',
                '-p', '60'
               ], 'gwas_mixed_read_through_exon.txt')
    lncrna_intergenic = ([
               'python', '/mctp/users/yniknafs/scripts/workspace/local_scripts/GWAS/GWAS_locus_shuffle_OR_parallel.py',
                '--gtf', '/mctp/projects/mitranscriptome/gwas/gtfs/tcat/tgenic/lncrna_intergenic.gtf',
                '--shuffs', '100',
                '-p', '60'
               ], 'gwas_lncrna_intergenic.txt')
    
    lncrna_intergenic_exon = ([
               'python', '/mctp/users/yniknafs/scripts/workspace/local_scripts/GWAS/GWAS_locus_shuffle_OR_parallel.py',
                '--gtf', '/mctp/projects/mitranscriptome/gwas/gtfs/tcat/tgenic/lncrna_intergenic.gtf',
                '--shuffs', '100',
                '--exon',
                '-p', '60'
               ], 'gwas_lncrna_exon_intergenic.txt')
    lncrna_intragenic = ([
               'python', '/mctp/users/yniknafs/scripts/workspace/local_scripts/GWAS/GWAS_locus_shuffle_OR_parallel.py',
                '--gtf', '/mctp/projects/mitranscriptome/gwas/gtfs/tcat/tgenic/lncrna_intragenic.gtf',
                '--shuffs', '100',
                '-p', '60'
               ], 'gwas_lncrna_intragenic.txt')
    
    lncrna_intragenic_exon = ([
               'python', '/mctp/users/yniknafs/scripts/workspace/local_scripts/GWAS/GWAS_locus_shuffle_OR_parallel.py',
                '--gtf', '/mctp/projects/mitranscriptome/gwas/gtfs/tcat/tgenic/lncrna_intragenic.gtf',
                '--shuffs', '100',
                '--exon',
                '-p', '60'
               ], 'gwas_lncrna_exon_intragenic.txt')
    tucp_intergenic = ([
               'python', '/mctp/users/yniknafs/scripts/workspace/local_scripts/GWAS/GWAS_locus_shuffle_OR_parallel.py',
                '--gtf', '/mctp/projects/mitranscriptome/gwas/gtfs/tcat/tgenic/tucp_intergenic.gtf',
                '--shuffs', '100',
                '-p', '60'
               ], 'gwas_tucp_intergenic.txt')
    
    tucp_intergenic_exon = ([
               'python', '/mctp/users/yniknafs/scripts/workspace/local_scripts/GWAS/GWAS_locus_shuffle_OR_parallel.py',
                '--gtf', '/mctp/projects/mitranscriptome/gwas/gtfs/tcat/tgenic/tucp_intergenic.gtf',
                '--shuffs', '100',
                '--exon',
                '-p', '60'
               ], 'gwas_tucp_exon_intergenic.txt')
    tucp_intragenic = ([
               'python', '/mctp/users/yniknafs/scripts/workspace/local_scripts/GWAS/GWAS_locus_shuffle_OR_parallel.py',
                '--gtf', '/mctp/projects/mitranscriptome/gwas/gtfs/tcat/tgenic/tucp_intragenic.gtf',
                '--shuffs', '100',
                '-p', '60'
               ], 'gwas_tucp_intragenic.txt')
    
    tucp_intragenic_exon = ([
               'python', '/mctp/users/yniknafs/scripts/workspace/local_scripts/GWAS/GWAS_locus_shuffle_OR_parallel.py',
                '--gtf', '/mctp/projects/mitranscriptome/gwas/gtfs/tcat/tgenic/tucp_intragenic.gtf',
                '--shuffs', '100',
                '--exon',
                '-p', '60'
               ], 'gwas_tucp_exon_intragenic.txt')
    
    
    
#     lncrna_exon = ([
#                'python', '/mctp/users/yniknafs/scripts/workspace/local_scripts/GWAS/GWAS_locus_shuffle_OR_parallel.py',
#                 '--gwas', '/mctp/projects/mitranscriptome/gwas/bed_files/gwas_xconv.bed',
#                 '--excl', '/mctp/projects/mitranscriptome/gwas/bed_files/gap.bed',
#                 '--gtf', '/mctp/projects/mitranscriptome/annotation/transcript_categories/assembly.annotation.lincRNA.gtf',
#                 '--snps', '/mctp/projects/mitranscriptome/gwas/bed_files/snp_arrays.sorted.bed',
#                 '--shuffs', '100',
#                 '-p', '60',
#                 '--exon'
#                ], 'gwas_lincRNA_exon.txt')
# 
#     lncRNA_genic = ([
#                'python', '/mctp/users/yniknafs/scripts/workspace/local_scripts/GWAS/GWAS_locus_shuffle_OR_parallel.py',
#                 '--gwas', '/mctp/projects/mitranscriptome/gwas/bed_files/gwas_xconv.bed',
#                 '--excl', '/mctp/projects/mitranscriptome/gwas/bed_files/gap.bed',
#                 '--gtf', '/mctp/projects/mitranscriptome/annotation/transcript_categories/assembly.annotation.lncRNA_genic.gtf',
#                 '--snps', '/mctp/projects/mitranscriptome/gwas/bed_files/snp_arrays.sorted.bed',
#                 '--shuffs', '100',
#                 '-p', '50'
#                ], 'gwas_lncRNA_genic.txt')
#     
#     lncRNA_genic_exon = ([
#                'python', '/mctp/users/yniknafs/scripts/workspace/local_scripts/GWAS/GWAS_locus_shuffle_OR_parallel.py',
#                 '--gwas', '/mctp/projects/mitranscriptome/gwas/bed_files/gwas_xconv.bed',
#                 '--excl', '/mctp/projects/mitranscriptome/gwas/bed_files/gap.bed',
#                 '--gtf', '/mctp/projects/mitranscriptome/annotation/transcript_categories/assembly.annotation.lncRNA_genic.gtf',
#                 '--snps', '/mctp/projects/mitranscriptome/gwas/bed_files/snp_arrays.sorted.bed',
#                 '--shuffs', '100',
#                 '-p', '50',
#                 '--exon'
#                ], 'gwas_lncRNA_genic_exon.txt')
#     
#     protein = ([
#                'python', '/mctp/users/yniknafs/scripts/workspace/local_scripts/GWAS/GWAS_locus_shuffle_OR_parallel.py',
#                 '--gwas', '/mctp/projects/mitranscriptome/gwas/bed_files/gwas_xconv.bed',
#                 '--excl', '/mctp/projects/mitranscriptome/gwas/bed_files/gap.bed',
#                 '--gtf', '/mctp/projects/mitranscriptome/annotation/transcript_categories/assembly.annotation.protein_coding.gtf',
#                 '--snps', '/mctp/projects/mitranscriptome/gwas/bed_files/snp_arrays.sorted.bed',
#                 '--shuffs', '100',
#                 '-p', '50'
#                ], 'gwas_protein.txt')    
#     
#     protein_exon = ([
#                'python', '/mctp/users/yniknafs/scripts/workspace/local_scripts/GWAS/GWAS_locus_shuffle_OR_parallel.py',
#                 '--gwas', '/mctp/projects/mitranscriptome/gwas/bed_files/gwas_xconv.bed',
#                 '--excl', '/mctp/projects/mitranscriptome/gwas/bed_files/gap.bed',
#                 '--gtf', '/mctp/projects/mitranscriptome/annotation/transcript_categories/assembly.annotation.protein_coding.gtf',
#                 '--snps', '/mctp/projects/mitranscriptome/gwas/bed_files/snp_arrays.sorted.bed',
#                 '--shuffs', '100',
#                 '-p', '50',
#                 '--exon'
#                ], 'gwas_protein_exon.txt')
#     
#     intergenic = ([
#                'python', '/mctp/users/yniknafs/scripts/workspace/local_scripts/GWAS/GWAS_locus_shuffle_OR_parallel.py',
#                 '--gwas', '/mctp/projects/mitranscriptome/gwas/bed_files/gwas_xconv.bed',
#                 '--excl', '/mctp/projects/mitranscriptome/gwas/bed_files/ref.merged.gaps.sorted.bed',
#                 '--gtf', '/mctp/projects/mitranscriptome/gwas/gtfs/intergenic.gtf',
#                 '--snps', '/mctp/projects/mitranscriptome/gwas/bed_files/snp_arrays.sorted.bed',
#                 '--shuffs', '1000',
#                 '-p', '80'
#                ], 'gwas_intergenic.txt')
#     
#     intergenic_exon = ([
#                'python', '/mctp/users/yniknafs/scripts/workspace/local_scripts/GWAS/GWAS_locus_shuffle_OR_parallel.py',
#                 '--gwas', '/mctp/projects/mitranscriptome/gwas/bed_files/gwas_xconv.bed',
#                 '--excl', '/mctp/projects/mitranscriptome/gwas/bed_files/ref.merged.gaps.sorted.bed',
#                 '--gtf', '/mctp/projects/mitranscriptome/gwas/gtfs/intergenic.gtf',
#                 '--snps', '/mctp/projects/mitranscriptome/gwas/bed_files/snp_arrays.sorted.bed',
#                 '--shuffs', '1000',
#                 '-p', '80',
#                 '--exon'
#                ], 'gwas_intergenic_exon.txt')
#     ref = ([
#                 'python', '/mctp/users/yniknafs/scripts/workspace/local_scripts/GWAS/GWAS_locus_shuffle_OR_parallel.py',
#                 '--gwas', '/mctp/projects/mitranscriptome/gwas/bed_files/gwas_xconv.bed',
#                 '--excl', '/mctp/projects/mitranscriptome/gwas/bed_files/gap.bed',
#                 '--gtf', '/mctp/projects/mitranscriptome/annotation/ref.merged.filterchrom.gtf',
#                 '--snps', '/mctp/projects/mitranscriptome/gwas/bed_files/snp_arrays.sorted.bed',
#                 '--shuffs', '100',
#                 '-p', '50'
#                 ], 'gwas_ref.txt')
#     
#     ref_exon = ([
#                 'python', '/mctp/users/yniknafs/scripts/workspace/local_scripts/GWAS/GWAS_locus_shuffle_OR_parallel.py',
#                 '--gwas', '/mctp/projects/mitranscriptome/gwas/bed_files/gwas_xconv.bed',
#                 '--excl', '/mctp/projects/mitranscriptome/gwas/bed_files/gap.bed',
#                 '--gtf', '/mctp/projects/mitranscriptome/annotation/ref.merged.filterchrom.gtf',
#                 '--snps', '/mctp/projects/mitranscriptome/gwas/bed_files/snp_arrays.sorted.bed',
#                 '--shuffs', '100',
#                 '-p', '50',
#                 '--exon'
#                 ], 'gwas_ref_exon.txt')
#     
#     
#     
#     assembly_nonprot = ([
#                 'python', '/mctp/users/yniknafs/scripts/workspace/local_scripts/GWAS/GWAS_locus_shuffle_OR_parallel.py',
#                 '--gwas', '/mctp/projects/mitranscriptome/gwas/bed_files/gwas_xconv.bed',
#                 '--excl', '/mctp/projects/mitranscriptome/gwas/bed_files/assembly.annotation.protein_coding.gap.bed',
#                 '--gtf', '/mctp/projects/mitranscriptome/annotation/assembly.gtf',
#                 '--snps', '/mctp/projects/mitranscriptome/gwas/bed_files/snp_arrays.sorted.bed',
#                 '--shuffs', '100',
#                 '-p', '50'
#                 ], 'gwas_assembly_nonprot.txt')
#     
#     assembly_exon_nonprot = ([
#                 'python', '/mctp/users/yniknafs/scripts/workspace/local_scripts/GWAS/GWAS_locus_shuffle_OR_parallel.py',
#                 '--gwas', '/mctp/projects/mitranscriptome/gwas/bed_files/gwas_xconv.bed',
#                 '--excl', '/mctp/projects/mitranscriptome/gwas/bed_files/assembly.annotation.protein_coding.gap.bed',
#                 '--gtf', '/mctp/projects/mitranscriptome/annotation/assembly.gtf',
#                 '--snps', '/mctp/projects/mitranscriptome/gwas/bed_files/snp_arrays.sorted.bed',
#                 '--shuffs', '100',
#                 '-p', '50',
#                 '--exon'
#                 ], 'gwas_assembly_exon_nonprot.txt')
#     
#     lincRNA_nonprot = ([
#                'python', '/mctp/users/yniknafs/scripts/workspace/local_scripts/GWAS/GWAS_locus_shuffle_OR_parallel.py',
#                 '--gwas', '/mctp/projects/mitranscriptome/gwas/bed_files/gwas_xconv.bed',
#                 '--excl', '/mctp/projects/mitranscriptome/gwas/bed_files/assembly.annotation.protein_coding.gap.bed',
#                 '--gtf', '/mctp/projects/mitranscriptome/annotation/transcript_categories/assembly.annotation.lincRNA.gtf',
#                 '--snps', '/mctp/projects/mitranscriptome/gwas/bed_files/snp_arrays.sorted.bed',
#                 '--shuffs', '100',
#                 '-p', '50'
#                ], 'gwas_lincRNA_nonprot.txt')
#     
#     lincRNA_exon_nonprot = ([
#                'python', '/mctp/users/yniknafs/scripts/workspace/local_scripts/GWAS/GWAS_locus_shuffle_OR_parallel.py',
#                 '--gwas', '/mctp/projects/mitranscriptome/gwas/bed_files/gwas_xconv.bed',
#                 '--excl', '/mctp/projects/mitranscriptome/gwas/bed_files/assembly.annotation.protein_coding.gap.bed',
#                 '--gtf', '/mctp/projects/mitranscriptome/annotation/transcript_categories/assembly.annotation.lincRNA.gtf',
#                 '--snps', '/mctp/projects/mitranscriptome/gwas/bed_files/snp_arrays.sorted.bed',
#                 '--shuffs', '100',
#                 '-p', '50',
#                 '--exon'
#                ], 'gwas_lincRNA_exon_nonprot.txt')
# 
#     lncRNA_genic_nonprot = ([
#                'python', '/mctp/users/yniknafs/scripts/workspace/local_scripts/GWAS/GWAS_locus_shuffle_OR_parallel.py',
#                 '--gwas', '/mctp/projects/mitranscriptome/gwas/bed_files/gwas_xconv.bed',
#                 '--excl', '/mctp/projects/mitranscriptome/gwas/bed_files/assembly.annotation.protein_coding.gap.bed',
#                 '--gtf', '/mctp/projects/mitranscriptome/annotation/transcript_categories/assembly.annotation.lncRNA_genic.gtf',
#                 '--snps', '/mctp/projects/mitranscriptome/gwas/bed_files/snp_arrays.sorted.bed',
#                 '--shuffs', '100',
#                 '-p', '50'
#                ], 'gwas_lncRNA_genic_nonprot.txt')
#     
#     lncRNA_genic_exon_nonprot = ([
#                'python', '/mctp/users/yniknafs/scripts/workspace/local_scripts/GWAS/GWAS_locus_shuffle_OR_parallel.py',
#                 '--gwas', '/mctp/projects/mitranscriptome/gwas/bed_files/gwas_xconv.bed',
#                 '--excl', '/mctp/projects/mitranscriptome/gwas/bed_files/assembly.annotation.protein_coding.gap.bed',
#                 '--gtf', '/mctp/projects/mitranscriptome/annotation/transcript_categories/assembly.annotation.lncRNA_genic.gtf',
#                 '--snps', '/mctp/projects/mitranscriptome/gwas/bed_files/snp_arrays.sorted.bed',
#                 '--shuffs', '100',
#                 '-p', '50',
#                 '--exon'
#                ], 'gwas_lncRNA_genic_exon_nonprot.txt')
#     
#     protein_nonprot = ([
#                'python', '/mctp/users/yniknafs/scripts/workspace/local_scripts/GWAS/GWAS_locus_shuffle_OR_parallel.py',
#                 '--gwas', '/mctp/projects/mitranscriptome/gwas/bed_files/gwas_xconv.bed',
#                 '--excl', '/mctp/projects/mitranscriptome/gwas/bed_files/assembly.annotation.protein_coding.gap.bed',
#                 '--gtf', '/mctp/projects/mitranscriptome/annotation/transcript_categories/assembly.annotation.protein_coding.gtf',
#                 '--snps', '/mctp/projects/mitranscriptome/gwas/bed_files/snp_arrays.sorted.bed',
#                 '--shuffs', '100',
#                 '-p', '50'
#                ], 'gwas_protein_nonprot.txt')    
#     
#     protein_exon_nonprot = ([
#                'python', '/mctp/users/yniknafs/scripts/workspace/local_scripts/GWAS/GWAS_locus_shuffle_OR_parallel.py',
#                 '--gwas', '/mctp/projects/mitranscriptome/gwas/bed_files/gwas_xconv.bed',
#                 '--excl', '/mctp/projects/mitranscriptome/gwas/bed_files/assembly.annotation.protein_coding.gap.bed',
#                 '--gtf', '/mctp/projects/mitranscriptome/annotation/transcript_categories/assembly.annotation.protein_coding.gtf',
#                 '--snps', '/mctp/projects/mitranscriptome/gwas/bed_files/snp_arrays.sorted.bed',
#                 '--shuffs', '100',
#                 '-p', '50',
#                 '--exon'
#                ], 'gwas_protein_exon_nonprot.txt')
    
    
    
    pathbio8 = [assembly, 
                lncrna, 
                lncrna_intergenic, 
                lncrna_intragenic,
                tucp,
                tucp_intergenic,
                tucp_intragenic,
                protein_coding,
                mixed_read_through,
                pseudogene
                ]
    
    pathbio9 = [assembly_exon, 
                lncrna_exon, 
                lncrna_intergenic_exon, 
                lncrna_intragenic_exon,
                tucp_exon,
                tucp_intergenic_exon,
                tucp_intragenic_exon,
                protein_coding_exon,
                mixed_read_through_exon,
                pseudogene_exon
                ]
    
    
    
    if args.pathbio == '8':
        for run_args, file_name in pathbio8: 
            logging.info("PATHBIO8 run")
            logging.info('Running %s' % file_name)
            with open(file_name, 'w') as fileh:
                    subprocess.call(run_args, stdout=fileh)
    if args.pathbio == '9':
        for run_args, file_name in pathbio9: 
            logging.info("PATHBIO9 run")
            logging.info('Running %s' % file_name)
            with open(file_name, 'w') as fileh:
                    subprocess.call(run_args, stdout=fileh)
    
    
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())
