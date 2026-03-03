# Deduper
This script (along with a sorted sam and a text file containing the known UMIs) remove all PCR duplicates, while retaining only a single copy of each read, from a sorted SAM file of uniquely mapped reads. Python 3.12 compatible code.
    
    `samtools sort` needs to be use outside of the python script
   


## Usage
UMI text file should have all known UMIS with one UMI per line 

sam files need to be sorted (chromosome tracking is reset when chromosome change is detected)

Files should not be gzipped

Script accounts of soft clipping and only single-end reads


Example of how to run script:

python3 deduper.py -f input_file.sam -o output_file.sam -u umi_file.txt

    - -f, --file: designates absolute file path to sorted sam file
    
    - -o, --outfile: designates absolute file path to deduplicated sam file
    
    - -u, --umi: designates file containing the list of UMIs



Example Input of a sorted SAM:
```
@HD     VN:1.4
@SQ     SN:1    LN:195154279
NS500451:154:HWKTMBGXX:4:23611:3821:5431:ATCGAACC       0       19      5848248 255     1S71M   *       0       0       ACCTGCCTTAAA
NS500451:154:HWKTMBGXX:4:23611:3821:5431:ATCGAACC       0       19      5848248 255     1S71M   *       0       0       ACCTGCCTTAAA
NS500451:154:HWKTMBGXX:4:23611:17006:5432:GAGAAGTC      0       17      56414979        255     1S71M   *       0       0       TCCCCGCTCTTCC
```

Output:
```
@HD     VN:1.4
@SQ     SN:1    LN:195154279
NS500451:154:HWKTMBGXX:4:23611:3821:5431:ATCGAACC       0       19      5848248 255     1S71M   *       0       0       ACCTGCCTTAAA
NS500451:154:HWKTMBGXX:4:23611:17006:5432:GAGAAGTC      0       17      56414979        255     1S71M   *       0       0       TCCCCGCTCTTCC
```




- Account for: 
    - all possible CIGAR strings (including adjusting for soft clipping, etc.)
    - Single-end reads
    - Known UMIs
- Considerations:
    - Millions of reads – avoid loading everything into memory!
    - Be sure to utilize functions appropriately
    - Appropriately comment code and include doc strings
- **CHALLENGE**: In a **separate branch**, implement options for
    - Single-end vs paired-end
    - Known UMIs vs randomers
    - Error correction of known UMIs
    - Choice of duplicate written to file
    
