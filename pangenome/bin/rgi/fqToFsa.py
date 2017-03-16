import sys
import os
import filepaths
from Bio import SeqIO

working_directory = os.getcwd()
script_path = filepaths.determine_path()

path = script_path

def convertFqToFsa(fq_path,fa_path):
	SeqIO.convert(fq_path, "fastq", fa_path, "fasta")

def main(argvfile):
	file_name = os.path.basename(argvfile)
	convertFqToFsa(argvfile,working_directory+'/'+file_name+'.read.fsa')


if __name__ == "__main__":
    main(sys.argv[1])
