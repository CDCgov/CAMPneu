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
    relevant_info = {"AO":"ALT_COUNT","RO":"REF_COUNT"}
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

# ### Read AMR_defaults.tsv into a dictionary ###
# def get_AMR_defaults(defaults):
#     '''
#     input:
#         defaults - path of file with AMR gene position default AA/Nucs 
#     output:
#         hgt_defaults - locuses of AMR genes that are (mostly) NOT present in FA19 / horizontally transferred genes with positions of interest as keys (nested dictionary)
#         WG_Defaults - contains defaults for positions in AMR genes that are present in FA19 with positions of interest as keys (nested dictionary)
#         mtrR_promoter - contains nucleotide positions as keys and default nucleotides as values for the mtrR promoter (dictionary)
#     '''
#     hgt_defaults = {}
#     WG_defaults = {}
#     mtrR_promoter = {}
#     with open(defaults) as file: 
#         count = 0
#         for line in file:
#             if count > 0:
#                 fields = line.strip().split('\t')
#                 if fields[2] == "CP012026": #Locus for FA19
#                     if fields[0] == "mtrR promoter": #this promoter contains 5 positions, won't work well in WG dictionary
#                         positions=fields[3].split(',')
#                         for i in range(len(positions)):
#                             mtrR_promoter[positions[i]] = fields[5][i]
#                     else:
#                         WG_defaults[fields[0]] = {"Gene":fields[1], "Locus":fields[2], "Nucleotide Position":fields[3], "AA Position":fields[4], "Default":fields[5]}
#                 else:
#                     hgt_defaults[fields[0]] = {"Gene":fields[1], "Locus":fields[2], "Nucleotide Position":fields[3], "AA Position":fields[4], "Default":fields[5]}
#             count += 1
#     return hgt_defaults, WG_defaults, mtrR_promoter

# ### Check all whole genome positions with default values to see if mutations appear at those positions ###
# def get_FA19_calls(WG_defaults,file,AA):
#     '''
#     input:
#         WG_defaults - contains defaults for positions in AMR genes that are present in FA19 with positions of interest as keys (nested dictionary)
#         file - name of tab separated vcf file with just relevant AMR genes (string)
#         AA - amino acid multi letter code as keys and single letter code as values (dictionary)
#     output:
#         results - stores variants for sample with positions as keys and found variants or defaults as values (dictionary)
#     '''
#     results = {}    
#     for field in WG_defaults.keys():
#         if field != 'penA D345ins':
#             if WG_defaults[field]["Nucleotide Position"] != 'NA':
#                 pattern = '.*[^\S]+' + WG_defaults[field]["Nucleotide Position"] + '[^\S]+.*'
#                 found,variant = run_grep(pattern,file)
#                 if found and len(variant[4]) == 1:
#                     results[field] = variant[4]
#                 else: 
#                     results[field] = WG_defaults[field]["Default"]
#             else:
#                 if WG_defaults[field]["Gene"] == 'ponA':
#                     pattern = '.*p(.)[a-zA-Z]+' + WG_defaults[field]["AA Position"] + '[a-zA-Z]+.*' + 'mrcA' + '.*' #ponA is annotated as mrcA in FA19
#                 elif WG_defaults[field]["Gene"] == 'mtrD':
#                     pattern = '.*p(.)[a-zA-Z]+' + WG_defaults[field]["AA Position"] + '[a-zA-Z]+.*' + 'mexB' + '.*' #mtrD is annotated as mexB in FA19
#                 else: pattern = '.*p(.)[a-zA-Z]+' + WG_defaults[field]["AA Position"] + '[a-zA-Z]+.*' + WG_defaults[field]["Gene"] + '.*' 
#                 found,variant = run_grep(pattern,file)
#                 if found:
#                     result = variant[10].split(' ')[-1].split(WG_defaults[field]["AA Position"]) #parse resulting AA from "EFFECT" column
#                     if len(result) == 2 and len(result[1]) == 3:
#                         aa = result[1]
#                         results[field] = AA[aa]
#                     elif len(result) == 2 and len(result[1]) > 3: #deal with complex mutations that START at the AA position of interest (takes the first AA change in the complex mutation)
#                         aa = ''.join(list(result[1])[0:3])
#                         results[field] = AA[aa]
#                     else:
#                         results[field] = WG_defaults[field]["Default"] 
#                 else: 
#                     results[field] = WG_defaults[field]["Default"]
#         else: results[field] = penA_D345ins(WG_defaults,WG_defaults[field]["Gene"],field,file)
#     return results

# ### Check vcf as far as 5 AA positions or 15 nucleotide positions back for complex mutations ###
# def check_vcf(poi,file,lookback,AA):
#     '''
#     input: 
#         poi - position of interest dictionary with "Gene", "Locus", "Nucleotide Position", "AA Position", and "Default" as keys with the corresponding information as values (dictionary)
#         file - name of tab separated vcf file with just relevant AMR genes (string)
#         lookback - number of positions to look back for considering complex mutations (int)
#         AA - amino acid multi letter code as keys and single letter code as values (dictionary)
#     output:
#         not_found - True if no complex variants impacting the point of interest are found (bool) 
#         variant - amino acid or nucleotide found at position of interest within complex mutation
#     '''
#     not_found = True 
#     variant = ''
#     for i in range(1,lookback+1):

#         if poi['AA Position'] != 'NA':
#             position = str(int(poi['AA Position']) - i)
#             if poi['Gene'] == 'ponA':
#                 pattern = '.*complex.*.*p(.)[a-zA-Z]+' + position + '[a-zA-Z]+.*' + 'mrcA' + '.*' #ponA is annotated as mrcA in FA19
#             elif poi['Gene'] == 'mtrD':
#                  pattern = '.*complex.*.*p(.)[a-zA-Z]+' + position + '[a-zA-Z]+.*' + 'mexB' + '.*' #mtrD is annotated as mexB in FA19
#             else: pattern = '.*complex.*.*p(.)[a-zA-Z]+' + position + '[a-zA-Z]+.*' + poi['Gene'] + '.*' 
#         else: 
#             position = str(int(poi['Nucleotide Position']) - i)
#             pattern = '^' + poi["Locus"] + '[^\S]+' + position + '[^\S]+complex.*' 

#         found,variant = run_grep(pattern,file)
#         if found:
#             if poi['AA Position'] != 'NA':
#                 AA_list = variant[10].split(' ')[-1].split(position) #parse resulting AA from "EFFECT" column 
#                 num_aa = i*3
#                 if len(AA_list[1]) > num_aa: #check if complex variant reaches position of interest
#                     aa = ''.join(list(AA_list[1])[num_aa:(num_aa+3)])
#                     variant = AA[aa]
#                     not_found = False 
#             else:
#                 nuc_list = list(variant[4])
#                 if len(nuc_list) > i: #check if complex variant reaches position of interest
#                     variant = nuc_list[i]
#                     not_found = False
#             break #if the variant does not reach the point of interest no further away variants will and the point of interest will stay set to wild type default, so we exit either way

#     return not_found,variant

# ### Check if snippy called complex mutations impacting positions of interest ###
# def check_complex(results,file,hgt_defaults,WG_defaults,AA):
#     '''
#     input:
#         results - stores variants for sample with positions as keys and found variants or defaults as values (dictionary)
#         file - name of tab separated vcf file with just relevant AMR genes (string)
#         hgt_defaults - locuses of AMR genes that are (mostly) NOT present in FA19 / horizontally transferred genes with positions of interest as keys (nested dictionary)
#         WG_defaults - contains defaults for positions in AMR genes that are present in FA19 with positions of interest as keys (nested dictionary)
#         AA - amino acid multi letter code as keys and single letter code as values (dictionary)
#     output:
#         results - results dictionary updated to include mutations at positions of interest due to complex mutations
#     '''
#     not_found = True

#     for result in results:

#         if result in hgt_defaults and re.search('freq',result) == None:
#             if results[result] == hgt_defaults[result]['Default']:
#                 if hgt_defaults[result]['Nucleotide Position'] != 'NA':
#                     not_found,variant = check_vcf(hgt_defaults[result],file,30,AA)
#                 else: not_found,variant = check_vcf(hgt_defaults[result],file,10,AA)

#         elif result in WG_defaults and result != 'penA D345ins':
#             if results[result] == WG_defaults[result]['Default']:
#                 if WG_defaults[result]['Nucleotide Position'] != 'NA':
#                     not_found,variant = check_vcf(WG_defaults[result],file,30,AA)
#                 else: not_found,variant = check_vcf(WG_defaults[result],file,10,AA)

#         if not not_found: #sorry for the double neg
#             results[result] = variant
#             not_found = True #reset this for the next position of interest

#     return results

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
    amr_positions = parse_tsv(snps_path,"Name")

    variants = parse_vcf(input_path)