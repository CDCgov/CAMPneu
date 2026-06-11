#!/usr/bin/env python3

import argparse
import re
import sys

### Read input TSV into dictionary ###
def parse_tsv(file_path,key):
    '''
    input:
        file_path - path to input TSV file (str)
        key - column to use as keys (str)
    output:
        contents - dictionary with specified key as keys and other columns as values (dict)
    '''
    contents = {}
    with open(file_path) as file:
        first = True
        for line in file:
            if first:
                header = line.strip().split('\t')
                first = False
                try:
                    key_index = header.index(key)
                except ValueError:
                    print(f"{key} column not found in {file_path}")
                    sys.exit()
            else:
                fields = line.strip().split('\t')
                contents[fields[key_index]] = {}
                for col in range(len(header)):
                    if col != key_index:
                        contents[fields[key_index]][header[col]] = fields[col]
    return contents

### Parse info column in annotated VCF ###
def parse_info(info_string, ann_info_header):
    '''
    input:
    output:
    '''
    relevant_info = {"AO":"ALT_COUNT","RO":"REF_COUNT","TYPE":"TYPE"}
    info_list = info_string.split(";")
    ann_info = info_list[-1].split("|")
    info = {}
    for col in info_list[0:-1]:
        for key in relevant_info:
            if col.startswith(key):
                info[relevant_info[key]] = col.split("=")[1]
    for i in range(len(ann_info_header)):
        info[ann_info_header[i]] = ann_info[i]
    return info

### Read input annotated VCF into dictionary of variants ###
def parse_vcf(input_path):
    '''
    input:
    output:
    '''
    variants = {}
    prev = ""
    ann_info_header = []
    count = 0
    with open(input_path) as file:
        for line in file:
            if not line.startswith("#"):
                if count == 0:
                    cols = prev.strip().split('\t')
                fields = line.strip().split('\t')
                variants[count] = {}
                for i in range(len(cols)):
                    if cols[i] == "INFO":
                        info = parse_info(fields[i],ann_info_header)
                        variants[count][cols[i]] = info
                    else:
                        variants[count][cols[i]] = fields[i]
                count += 1
            else:
                if "Functional annotations" in line: #extract the header for the functional annotations
                    half = line.strip().split(":")
                    ann_info_string = re.sub(r"[\s,',>]", "", half[1])
                    ann_info_header = ann_info_string.split("|")
                prev = line
    return variants

### format variant for addition to dict containing variants to report ###
def format_var(pos,var,type):
    '''
    '''
    key = pos['MPN']+"-"+pos[type] #same key format as initial amr_pos dict
    value = {}
    value['GENE'] = pos['Gene']
    value['POS'] = var['POS']
    value['REF'] = var['REF']
    value['ALT'] = var['ALT']
    value['ALT_COUNT'] = var['INFO']['ALT_COUNT']
    value['QUAL'] = var['QUAL']
    if type == 'Nucleotide Position':
        if pos['MPN'] == 'MPN_RS00530': #23S
            value['EFFECT'] = value['REF'] + str(int(var['POS'])-120056) + value['ALT']
    else:
        value['EFFECT'] = var['INFO']['HGVS.p']
    value['DEPTH'] = int(var['INFO']['REF_COUNT'])+int(value['ALT_COUNT'])
    alt_frac = int(value['ALT_COUNT'])/value['DEPTH']
    if value['DEPTH'] < 10 or alt_frac < 0.9:
        value['PASS/FAIL'] = 'FAIL'
    else:
        value['PASS/FAIL'] = 'PASS'
    
    return key,value

### check annotated VCF for known AMR variants ###
def find_amr(amr_positions,variants,nucl_pos):
    '''
    '''
    amr_vars = {}
    for gene in amr_positions:
        amr_vars[list(amr_positions[gene].values())[0]['Drug_Class']] = {} #keep track of all of our drug classes
    for var in variants:
        var_pos = variants[var]['POS']
        if variants[var]['INFO']['Gene_ID'] in amr_positions or var_pos in nucl_pos: #for rRNA the MPN nums won't match
            if var_pos in nucl_pos:
                id = nucl_pos[var_pos]
            else:
                id = variants[var]['INFO']['Gene_ID']
            snp = (variants[var]['INFO']['TYPE'] == 'snp')
            if var_pos in amr_positions[id]: #this should only happen for nucl variants but we'll double check
                if amr_positions[id][var_pos]['Nucleotide Position'] != "NA":
                    if snp and variants[var]['ALT'] in amr_positions[id][var_pos]['Alt']:
                        key,value = format_var(amr_positions[id][var_pos],variants[var],'Nucleotide Position')
                        amr_vars[amr_positions[id][var_pos]['Drug_Class']][key] = value #group by drug class
    print(amr_vars)

                    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
                        prog='python3 amr_variant_parsing.py',
                        description='This program performs the variant analysis for Mp relative to M129 using a SnpEff annotated VCF')
    parser.add_argument('-i','--input',type=str,required=True,help='Full path to annotated VCF file from SnpEff (.ann.vcf)')
    parser.add_argument('-g','--genes',type=str,required=True,help='Full path to TSV file containing AMR gene information')
    parser.add_argument('-s','--snps',type=str,required=True,help='Full path to TSV file containing AMR SNP information')
    parser.add_argument('-o','--out_path', type=str, required=True, help='path of output directory')
    
        
    args = parser.parse_args()

    input_path = args.input
    genes_path = args.genes
    snps_path = args.snps
    out = args.out_path

    amr_genes = parse_tsv(genes_path,"MPN")
    amr_pos = parse_tsv(snps_path,"Name")
    
    # reformat dictionary to make it easier to search
    amr_positions = {}
    nucl_pos = {} #nucleotide positions for rRNA variants
    for pos in amr_pos.keys():
        keys = pos.split('-')
        if keys[0] not in amr_positions:
            amr_positions[keys[0]] = { keys[1]: amr_pos[pos] }
        else:
            amr_positions[keys[0]][keys[1]] = amr_pos[pos]
        if amr_pos[pos]['Nucleotide Position'] != "NA":
            nucl_pos[keys[1]] = keys[0]

    #print(amr_positions)
    variants = parse_vcf(input_path)
    #l = list(variants.keys())
    #print(variants[l[0]]['INFO'])
    find_amr(amr_positions,variants,nucl_pos)