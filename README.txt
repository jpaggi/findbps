The MIT License (MIT)

Copyright (c) 20015 Joseph Paggi

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

#########################################################

This module can be invoked using either the command line or by importing this 
module directly.

#########################################################

To use via command line, use the following command:

python <reads> <output> <bowtie_options> <motif> <length> <threshold> <strand>

Where arguments are defined as follows:
reads:  path to file where single-end, stranded
        RNA-seq reads in fastq format are located
output: desired basename of output files
bowtie_options: bowtie options you wish to 
     	be used for alignment of reads after splitting.
        See the bowtie manual for details. 
        Recommend "-y -p 2 -v 0 -X 5000 -m 1 <index>"
	MUST BE PUT IN QUOTES
motif: 	list of dictionaries representing 5'ss motif
        position weight matrix. Each dictionary has a
        key for each nucleotide, with a float of the 
        probability as keys. Write in PYTHON format
	and PUT IN QUOTES. 
length: lowest acceptable number of bases
        used to align a fragment of a read.
threshold: float of the lowest acceptable probability 
        that a sequence would be sampled from the 
        given martrix in order to attempt mapping.
        Recommend 0.0 unless many false positives
strand: str either 'first' if reads are first-stranded
       	or 'second' if reads are second-stranded

Will right reads with valid alignment to <output>.bed.
Information recorded in paired-end bed format.

Reads with no valid alignment are written to
<output>_no_alignment.fastq in the paired-end tab-
deliminated format described in the bowtie manual.
Reads are split as they were attempted to be aligned.

Below is an example command:

python findbps/findbps.py /path/to/reads /path/to/output_files "-y -p 2 -v 0 -X 5000 -m 1 index/path" '[{"A":0,"T":0,"C":0,"G":1},{"A":0,"T":.95,"C":.05,"G":0},{"A":1,"T":0,"C":0,"G":0},{"A":.02,"T":.8,"C":.18,"G":0},{"A":0,"T":0,"C":0,"G":1}]' 10 0 second

#######################################################

In order to use this module directly as a python library,
you must add this directory to your python path.

This can easily be done temporarily by using the command:
sys.path.append(<path to this directory>).

You can then import the function findbps() using:
from findbps import findbps

You can then use the findbps function directly using:
findbps(reads, output, bowtie_options, motif, length, 
	threshold, strand)

Arguments are the same as described above, except that 
arguments should be entered as thier respective types
instead of as strings. 

Below is an example:

from sys import path
path.append('~/findbps')
from findbps import findbps
#specify reads you want to process
reads = '/path/to/reads'

#specify basename of output file here
output = 'path/to/output_files'

#specify bowtie options
bowtie_options = "-y -p 2 -v 0 -X 5000 -m 1 index/path"

#specify 5'SS weight matrix
motif = [{"A":0,"T":0,"C":0,"G":1},
         {"A":0,"T":.95,"C":.05,"G":0},
         {"A":1,"T":0,"C":0,"G":0},
         {"A":.02,"T":.8,"C":.18,"G":0},
         {"A":0,"T":0,"C":0,"G":1},
         {"A":.05,"T":.9,"C":.05,"G":0}]

#specify shortest acceptable mapping length
length = 10

#specify lowest probability that a reads best
# 5'SS was sampled from the given matrix
threshold = 0

strand = 'second'

findbps(reads, output, bowtie_options, motif, length, threshold, strand)
