#this module takes standard bowtie paired end output and makes
#    it into a paired end bed file
#this code is really quite boring and tedious, and shouldnt need to be changed
#output is as follows: chromosome, start_upstream, end_upstream, 
#     chromosome, start_downstream, end_downstream, read ID
#     quality scores, strand, strand, upstream seq, downstream seq

# sequences are given for the + strand regardless of strand of alignment
# the branchpoint nucleotide is re-added in its respective position


def get_first_data(first):
    bp_nucleotide = first[0]
    bp_quality = first[1]
    first = first[2:]
    (ID, strand, chromosome, first_position, first_seq, first_quality,
     ceiling, first_mismatches) = first.split('\t')
    return (bp_nucleotide, bp_quality, ID, strand, chromosome,
            int(first_position), first_seq, first_quality, first_mismatches)

def get_second_data(second):
    (ID, strand, chromosome, second_position,second_seq, second_quality,
     ceiling, second_mismatches) = second.split('\t')
    return int(second_position), second_seq, second_quality, second_mismatches

def complement(n):
    if n == 'A':
        out = 'T'
    elif n == 'T':
        out = 'A'
    elif n == 'G':
        out = 'C'
    elif n == 'C':
        out = 'G'
    else:
        out = n 
    return out

def read_pair(inp):
#reads two lines and breaks up the first one to get information as described
    first=inp.readline()
    second=inp.readline()
    if first == second:
        return '','','','','','','','','',''
    (bp_nucleotide, bp_quality, ID, strand, chromosome, first_position, first_seq,
     first_quality, first_mismatches) = get_first_data(first)
    second_position, second_seq, second_quality, second_mismatches = get_second_data(second)

#when strand is plus first corresponds to the 5'SS part
    if strand=='+':
        fp_start = first_position
        fp_end = first_position + len(first_quality)
        fp_seq = first_seq
        
        bp_start = second_position
        bp_end = second_position + len(second_quality)
        bp_seq = second_seq + bp_nucleotide
        
        quality = second_quality + bp_quality + first_quality

#when strand is minus second corresponds to the 5'SS part
    else:
        fp_start = second_position - 1
        fp_end = second_position + len(second_quality) - 1
        fp_seq = second_seq
        
        bp_start = first_position - 1
        bp_end = first_position + len(first_quality)
        bp_seq = complement(bp_nucleotide) + first_seq 

        quality = first_quality + bp_quality + second_quality
    
    return (ID,strand,chromosome,quality,
            fp_start, fp_end, fp_seq,
            bp_start, bp_end, bp_seq)

def make_bed(reads, out):
    (ID,strand,chromosome,quality,
     fp_start, fp_end, fp_seq,
     bp_start, bp_end, bp_seq) = read_pair(reads)
    total = 0
    while ID:
        total += 1
        fp_start = str(fp_start) 
        fp_end = str(fp_end)
        bp_start = str(bp_start) 
        bp_end = str(bp_end)
        if strand == '+':
            start = chromosome + '\t' + fp_start + '\t' + fp_end + '\t'
            stop = chromosome + '\t' + bp_start + '\t' + bp_end + '\t'
            seq = fp_seq + '\t' + bp_seq + '\n' 
        else:
            start = chromosome + '\t' + bp_start + '\t' + bp_end + '\t'
            stop = chromosome + '\t' + fp_start + '\t' + fp_end + '\t'
            seq = bp_seq + '\t' + fp_seq + '\n' 
        line = start+stop+ID+'\t'+quality+'\t'+strand+'\t'+strand+'\t'+seq 
        out.write(line)
        (ID,strand,chromosome,quality,
         fp_start, fp_end, fp_seq,
         bp_start, bp_end, bp_seq) = read_pair(reads)
    return total

from sys import stdin, stdout

make_bed(stdin, stdout)
