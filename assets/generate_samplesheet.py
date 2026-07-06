#!/usr/bin/env python

import os
import sys
import glob
import argparse
import re

# adapted from viralrecon

def csv_file(value):
    if not value.lower().endswith(".csv"):
        raise argparse.ArgumentTypeError("SAMPLESHEET_FILE must be a .csv file")
    return value

def parse_args(args=None):
    Description = (
        "Generate samplesheet from a directory of FastQ/FastA files and/or a list of SRA accessions."
    )
    Epilog = "Example usage: python3 generate_samplesheet.py -i <SAMPLE_DIR> -s <SRA_ACCESSIONS_FILE> <SAMPLESHEET_FILE>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("SAMPLESHEET_FILE", type=csv_file, help="Output samplesheet file (.csv only).")
    parser.add_argument(
        "-i",
        "--input_dir",
        type=str,
        dest="SAMPLE_DIR",
        default="./",
        help="Folder containing raw FastQ and/or FastA files. (default: './')",
    )
    parser.add_argument(
        "-s",
        "--sra_accessions",
        type=str,
        dest="SRA_ACCESSIONS_FILE",
        default="sra_accessions.txt",
        help="Text file containing SRA accessions, 1 accession per line (default: sra_accessions.txt)",
    )
    parser.add_argument(
        "-sn",
        "--sanitise_name",
        dest="SANITISE_NAME",
        action="store_true",
        help="Whether to further sanitise FastQ file name to get sample id. Used in conjunction with --sanitise_name_delimiter and --sanitise_name_index.",
    )
    parser.add_argument(
        "-sd",
        "--sanitise_name_delimiter",
        type=str,
        dest="SANITISE_NAME_DELIMITER",
        default="_",
        help="Delimiter to use to sanitise sample name. (default: '_')",
    )
    parser.add_argument(
        "-si",
        "--sanitise_name_index",
        type=int,
        dest="SANITISE_NAME_INDEX",
        default=1,
        help="After splitting FastQ file name by --sanitise_name_delimiter all elements before this index (1-based) will be joined to create final sample name. (default: 1)",
    )
    parser.add_argument(
        "-sna",
        "--sanitise_name_fa",
        dest="SANITISE_NAME_FA",
        action="store_true",
        help="Whether to further sanitise FastA file name to get sample id. Used in conjunction with --sanitise_name_delimiter_fa and --sanitise_name_index_fa.",
    )
    parser.add_argument(
        "-sda",
        "--sanitise_name_delimiter_fa",
        type=str,
        dest="SANITISE_NAME_DELIMITER_FA",
        default="_",
        help="Delimiter to use to sanitise sample name. (default: '_')",
    )
    parser.add_argument(
        "-sia",
        "--sanitise_name_index_fa",
        type=int,
        dest="SANITISE_NAME_INDEX_FA",
        default=1,
        help="After splitting FastA file name by --sanitise_name_delimiter_fa all elements before this index (1-based) will be joined to create final sample name. (default: 1)",
    )
    return parser.parse_args(args)

def make_samplesheet(samplesheet_file):
    out_dir = os.path.dirname(samplesheet_file)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir)

    with open(samplesheet_file, "w") as fout:
        header = ["sample", "fastq_1", "fastq_2","fasta","sra_accession"]
        fout.write(",".join(header) + "\n")

def get_ext(dir):
    read1_regex = r'(_R1_001\.fastq\.gz$)|(_1\.fastq\.gz$)|(_R1\.fastq\.gz$)|(_R1_001\.fq\.gz$)|(_1\.fq\.gz$)|(_R1\.fq\.gz$)'
    read1_ext = { "_R1_001.fastq.gz": 0, "_1.fastq.gz": 0, "_R1.fastq.gz": 0, "_R1_001.fq.gz": 0, "_1.fq.gz": 0, "_R1.fq.gz": 0}

    fasta_regex = r'(\.fasta$)|(\.fa$)|(\.fas$)|(\.fna$)'
    fasta_ext = {".fasta": 0, ".fa": 0, ".fas": 0, ".fna": 0 }

    for filename in os.listdir(dir):
        read_match = re.search(read1_regex,filename)
        fasta_match = re.search(fasta_regex,filename)
        if read_match:
            if read_match.group(1):
                read1_ext[read_match.group(1)] += 1
            elif read_match.group(2):
                read1_ext[read_match.group(2)] += 1
            elif read_match.group(3):
                read1_ext[read_match.group(3)] += 1
            elif read_match.group(4):
                read1_ext[read_match.group(4)] += 1
            elif read_match.group(5):
               read1_ext[read_match.group(5)] += 1
            elif read_match.group(6):
                read1_ext[read_match.group(6)] += 1
        elif fasta_match:
            if fasta_match.group(1):
                fasta_ext[fasta_match.group(1)] += 1
            elif fasta_match.group(2):
                fasta_ext[fasta_match.group(2)] += 1
            elif fasta_match.group(3):
                fasta_ext[fasta_match.group(3)] += 1
            elif fasta_match.group(4):
                fasta_ext[fasta_match.group(4)] += 1
    return read1_ext,fasta_ext

def fastq_dir_to_samplesheet(
    fastq_dir,
    samplesheet_file,
    read1_extension="_R1_001.fastq.gz",
    read2_extension="_R2_001.fastq.gz",
    sanitise_name=False,
    sanitise_name_delimiter="_",
    sanitise_name_index=1,
):
    def sanitize_sample(path, extension):
        """Retrieve sample id from filename"""
        sample = os.path.basename(path).replace(extension, "")
        if sanitise_name:
            sample = sanitise_name_delimiter.join(
                os.path.basename(path).split(sanitise_name_delimiter)[
                    :sanitise_name_index
                ]
            )
        return sample

    def get_fastqs(extension):
        """
        Needs to be sorted to ensure R1 and R2 are in the same order
        when merging technical replicates. Glob is not guaranteed to produce
        sorted results.
        See also https://stackoverflow.com/questions/6773584/how-is-pythons-glob-glob-ordered
        """
        return sorted(
            glob.glob(os.path.join(fastq_dir, f"*{extension}"), recursive=False)
        )

    if (os.path.isdir(fastq_dir)):

        read_dict = {}

        ## Get read 1 files
        for read1_file in get_fastqs(read1_extension):
            sample = sanitize_sample(read1_file, read1_extension)
            if sample not in read_dict:
                read_dict[sample] = {"R1": [], "R2": []}
            read_dict[sample]["R1"].append(read1_file)

        ## Get read 2 files
        for read2_file in get_fastqs(read2_extension):
            sample = sanitize_sample(read2_file, read2_extension)
            read_dict[sample]["R2"].append(read2_file)

        ## Write to file
        if len(read_dict) > 0:
            with open(samplesheet_file, "a") as fout:
                for sample, reads in sorted(read_dict.items()):
                    for idx, read_1 in enumerate(reads["R1"]):
                        read_2 = ""
                        if idx < len(reads["R2"]):
                            read_2 = reads["R2"][idx]
                        sample_info = ",".join([sample, read_1, read_2,"",""])
                        fout.write(f"{sample_info}\n")
        else:
            error_str = (
                "\nERROR No paired FastQ files found with " + read1_extension + " and " + read2_extension + "samplesheet may be empty!\n\n"
            )
            error_str += "Check the contents of the provided directory for unpaired files:\n"
            print(error_str)
            exit(1)
    else:
        error_str = (
                "\nERROR input directory" + fastq_dir + " does not exist!\n\n"
            )
        error_str += "Check the values provided for the:\n"
        error_str += "  - '--input_dir' Path to the directory containing the FastQ files\n"
        print(error_str)
        exit(1)

def fasta_dir_to_samplesheet(
    fasta_dir,
    samplesheet_file,
    fasta_extension=".fasta",
    sanitise_name_fa=False,
    sanitise_name_delimiter_fa="_",
    sanitise_name_index_fa=1,
):
    def sanitize_sample(path, extension):
        """Retrieve sample id from filename"""
        sample = os.path.basename(path).replace(extension, "")
        if sanitise_name_fa:
            sample = sanitise_name_delimiter_fa.join(
                os.path.basename(path).split(sanitise_name_delimiter_fa)[
                    :sanitise_name_index_fa
                ]
            )
        return sample

    def get_fastas(extension):
        return sorted(
            glob.glob(os.path.join(fasta_dir, f"*{extension}"), recursive=False)
        )

    if (os.path.isdir(fasta_dir)):
        assembly_dict = {}

        ## Get assembly files
        for file in get_fastas(fasta_extension):
            sample = sanitize_sample(file, fasta_extension)
            assembly_dict[sample] = file

        ## Write to file
        if len(assembly_dict) > 0:
            out_dir = os.path.dirname(samplesheet_file)
            if out_dir and not os.path.exists(out_dir):
                os.makedirs(out_dir)

            with open(samplesheet_file, "a") as fout:
                for sample in sorted(assembly_dict.keys()):
                    sample_info = ",".join([sample,"","",assembly_dict[sample],""])
                    fout.write(f"{sample_info}\n")
    else:
        error_str = (
                "\nERROR input directory" + fasta_dir + " does not exist!\n\n"
            )
        error_str += "Check the values provided for the:\n"
        error_str += "  - '--input_dir' Path to the directory containing the FastQ files\n"
        print(error_str)
        exit(1)

def sra_list_to_samplesheet(
        sra_list="sra_accessions.txt",
        samplesheet_file=""
        
):
    if (os.path.isfile(sra_list)):
        sra_acc = []
        with open(sra_list, "r") as list:
            for acc in list:
                sra_acc.append(acc.strip())

        with open(samplesheet_file, "a") as fout:
                for sample in sra_acc:
                    sample_info = ",".join([sample,"","","",sample])
                    fout.write(f"{sample_info}\n")

    else:
        error_str = (
            "\nNo SRA accession list file found so none have been added to samplesheet!\n\n"
        )
        error_str += "If you intended to add SRA accessions please check the values provided for the:\n"
        error_str += "  - '--sra_accessions' Path to the file containing the list of SRA accessions\n"
        print(error_str)
    

def main(args=None):
    args = parse_args(args)

    make_samplesheet(args.SAMPLESHEET_FILE)

    read1_ext, fasta_ext = get_ext(args.SAMPLE_DIR)

    for ext1 in read1_ext:
        if read1_ext[ext1] > 0:
            ext2 = ext1.replace("1","2",1)
            fastq_dir_to_samplesheet(
                fastq_dir=args.SAMPLE_DIR,
                samplesheet_file=args.SAMPLESHEET_FILE,
                read1_extension=ext1,
                read2_extension=ext2,
                sanitise_name=args.SANITISE_NAME,
                sanitise_name_delimiter=args.SANITISE_NAME_DELIMITER,
                sanitise_name_index=args.SANITISE_NAME_INDEX,
            )
            log_str = (
                "\nSUCCESS! FastQ files found with " + ext1 + " extension have been added to samplesheet!\n"
            )
            print(log_str)
        else:
            error_str = (
                "No FastQ files found with " + ext1 + " extension so none have been added to samplesheet!"
            )
            print(error_str)

    for ext in fasta_ext:
        if fasta_ext[ext] > 0:
            fasta_dir_to_samplesheet(
                fasta_dir=args.SAMPLE_DIR,
                fasta_extension=ext,
                samplesheet_file=args.SAMPLESHEET_FILE,
                sanitise_name_fa=args.SANITISE_NAME_FA,
                sanitise_name_delimiter_fa=args.SANITISE_NAME_DELIMITER_FA,
                sanitise_name_index_fa=args.SANITISE_NAME_INDEX_FA,
            )
            log_str = (
                "\nSUCCESS! FastA files found with " + ext + " extension have been added to samplesheet!\n"
            )
            print(log_str)
        else:
            error_str = (
                "No FastA files found with " + ext + " extension so none have been added to samplesheet!"
            )
            print(error_str)

    sra_list_to_samplesheet(
        sra_list=args.SRA_ACCESSIONS_FILE,
        samplesheet_file=args.SAMPLESHEET_FILE
    )


if __name__ == "__main__":
    sys.exit(main())
