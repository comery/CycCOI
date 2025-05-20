# Cyclone-DNAbarcode-COI


### DESCRIPTION
A set of tootkit for dealing with COI amplicons using Cyclone sequencing platform

### INSTALLATION
- Clone from github
```bash
$ git clone https://github.com/comery/CycCOI.git
```
### Requirements 
#### (1) python modules

```
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
from collections import defaultdict
import argparse
import logging
import sys
import os
from multiprocessing import Pool
import multiprocessing
from icecream import ic
```



### DATA requirements:

#### (1) CycloneSEQ reads
#### (2) primers list
primer.list
```text
for	TAAACTTCTGGATGTCCAAAAAATCA
rev	TTTCAACAAATCATAAAGATATTGG
```



#### (3) index(barcodes for identifying samples in which plate) list

plate.index.tsv

```
P1	TCGGTCTTAGACG
P2	TGTGAAGTTGCCA
P3	AGATTCTACACAA
P4	ATGCGATTAATTG
P5	GGCTGTTACAACA
```


#### (4) index(barcodes for identifying samples in a plate) list

cell.index.tsv

```text
001     AAAGC 
002     AACAG 
003     AACCT 
004     AACTC 
005     AAGCA  
```


### Get start

1. QC, filtering sequencing by length, gc and generate report figures

```bash
python3 bin/CycFqFilter.py -q 7 -l 700 -L 770 -g 0.2 -G 0.6 -o test.clean test.fastq.gz
```

2. assign sequencing by plate index and well index

```bash
$ python3 ../bin/pcr_demultiplex.py -p primer.txt --plate-index plate.index.tsv --well-index cell.index.tsv -f test.fa  -o out  --primer_max_mismatch 3 --primer_max_indel 3 --index_max_mismatch 2 --index_max_indel 1
```





