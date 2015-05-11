from subprocess import Popen, PIPE
from pickle import dumps
from os import path

def findbps(reads, output, bowtie_options, motif, length, threshold, strand):
    """
    Input:
    reads: str of name of file where single-end, stranded 
           RNA-seq reads in fastq format are located
    output:str of desired basename of output files
    bowtie_options: str of bowtie options you wish to 
           be used for alignment of reads after splitting.
           See the bowtie manual. 
           Recommend "-y -p 2 -v 0 -X 5000 -m 1 <index>"
    motif: list of dictionaries representing 5'ss motif
           position weight matrix. Each dictionary has a
           key for each nucleotide, with a float of the 
           probability as keys. 
    length:int of the lowest acceptable number of bases
           used to align a fragment of a read.
    threshold: float of the lowest acceptable probability 
           that a sequence would be sampled from the 
           given martrix in order to attempt mapping.
           Recommend 0.0 unless many false positives
    strand:str either 'first' if reads are first-stranded
           or 'second' if reads are second-stranded
    Output:
    output + '.bed': 
           A file in paired-end bed format with
           information about the reads with a valid
           alignment. 
    output + '_no_alignment.fastq':
           Reads with no valid alignment in the 
           paired-end tab-delimited format 
           described in the bowtie manual split 
           as they were attempted to be aligned.
    """
           
    #gets the name of the directory of this file
    directory =  path.dirname(path.realpath(__file__))

    #make these arguments into strings so they can be passed to fp_checker.py
    motif = '"' + dumps(motif) + '"'
    length = str(length)
    threshold = str(threshold)

    #this process splits each read at the most likely 5'SS based on the 
    # given weight matrix and sends them to bowtie to be mapped
    # see fp_checker.py for further details
    fp_checker = Popen('python ' + directory + '/fp_checker.py ' + 
                       motif +' '+ length +' '+ threshold +' '+ strand,
                       stdin = open(reads,'r'), stdout = PIPE, shell = True)

    #this process maps each split read to the given genome 
    bowtie = Popen('bowtie --ff ' + bowtie_options + ' --12 - --un ' +
                   output+'_no_alignment.fastq',
                   stdin = fp_checker.stdout, stdout = PIPE, shell = True)

    fp_checker.stdout.close()

    #this process converts the bowtie output into a bed file
    # see make_bed.py for further details
    make_bed = Popen('python ' + directory + '/make_bed.py',
                     stdin = bowtie.stdout,
                     stdout = open(output + ".bed",'w'), shell = True)

    bowtie.stdout.close()
    make_bed.wait()
    return 0

if __name__ == '__main__':
    from sys import argv
    reads = argv[1]
    output = argv[2]
    bowtie_options = argv[3]
    motif = eval(argv[4])
    length = int(argv[5])
    threshold = float(argv[6])
    strand = argv[7]
    findbps(reads, output, bowtie_options, motif, length, threshold, strand)
