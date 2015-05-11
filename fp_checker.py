# This module processes reads given in fastq format
# Finds the subsequence that most closely matches the given 5'ss
# returns paired end reads broken just upstream of 5'ss
# returned in the tab delimited format as described in the bowtie manual

#call on command line with 4 arguments:
#motif: the string of a pickled dictionary representing weight matrix
#length: the shortest acceptable read length after splitting
#threshold: the minimum acceptable probability that the best motif would
#           be sampled from the given weight matrix
#strand: 'first' if first-stranded 'second' if second-stranded
#also need to set stdin to fastq file and stdout to bowtie or a file

def reverse_seq(seq):
    out = ''
    for letter in seq[::-1]:
        if letter == 'A':
            out += 'T'
        elif letter == 'T':
            out += 'A'
        elif letter == 'G':
            out += 'C'
        elif letter == 'C':
            out += 'G'
        else:
            out += letter
    return out

def interpret_fastq(fastq):
    """
    Given:
    fastq: an opened fastq file for which the next
    line is the beginning of an entry
    Return:
    ID: str of read sequence ID
    seq: str of read sequence
    quality: str of quality scores
    Note: importantly fastq is left
    such that the next line is the 
    beginning of the next entry
    Returns all empty strings at EOF
    """
    try:
        ID=fastq.readline()[1:-1].split(' ')[0]
    except IndexError:
        return '','',''
    seq=fastq.readline()[:-1]
    fastq.readline()
    quality=fastq.readline()[:-1]
    return ID, seq, quality

def score_seq(motif, seq):
    """
    Given:
    motif: list of dictionaries of nucleotide --> probability
    seq: str of part of read sequence must be same length as motif
    Return:
    score: float of probability that it would be sampled from motif
    """
    score= 1
    for x in range(len(motif)):
        #no N's allowed in motif
        try:
            #bayes rule transition
            score = score*motif[x][seq[x]]
        except KeyError:
            return 0
    return score

def find_break(seq, length, best_score, motif):
    """
    Search for most likely location of motif in read
    Given:
    seq: str of read sequence
    motif: list of dictionaries of nucleotide --> probability
    length: int minimum length of each fragment
    best_score: float threshold probability
    Return: int of most likely if exceeds threshold
    otherwise return int 0
    """
    break_point = 0
    for x in range(length, len(seq) - length + 1):
        score = score_seq(motif,seq[x:x+len(motif)])
        if score > best_score:
            best_score = score
            break_point = x
    return break_point

def fp_checker(reads, motif, length, threshold, output, first_stranded):
    """
    This is the only function you should call from outside
    this module.
    It runs the whole process and prints tab delimited reads
    to standard out.
    """
    ID, seq, quality = interpret_fastq(reads)
    while ID:
        if first_stranded:
            seq = reverse_seq(seq)
            quality = quality[::-1]
        break_point = find_break(seq, length, threshold,  motif)
        if break_point:
            #add bp nucleotide and quality score to begining of ID
            #this lets us maintain information without aligning it
            x = break_point
            output.write(seq[x-1] + quality[x-1] + ID + '\t' +
                         seq[x:] + '\t' + quality[x:] + '\t' + 
                         seq[:x-1] + '\t' + quality[:x-1] + '\n')
        ID, seq, quality = interpret_fastq(reads)
    output.close()
    return 0

##############################################################

from sys import stdin, stdout, argv
from pickle import loads

motif = loads(argv[1])
length = int(argv[2])
threshold = float(argv[3])
first_stranded = (argv[4] == 'first')

fp_checker(stdin, motif, length, threshold, stdout, first_stranded)
