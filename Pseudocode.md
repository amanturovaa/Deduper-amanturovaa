{\rtf1\ansi\ansicpg1252\cocoartf2822
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fmodern\fcharset0 Courier;}
{\colortbl;\red255\green255\blue255;\red0\green0\blue0;}
{\*\expandedcolortbl;;\cssrgb\c0\c0\c0;}
\margl1440\margr1440\vieww17160\viewh13200\viewkind0
\deftab720
\pard\pardeftab720\partightenfactor0

\f0\fs26 \cf0 \expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 Process sorted SAM file, which contains uniquely mapped single end reads, and remove PCR duplicated\
PCR duplicated contain the same RNAME, POS, FLAG, and UMI\
Keep the first unique read encountered and remove any subsequent reads.\
\
\pard\pardeftab720\partightenfactor0
\cf0 input\
    sam file \
    file with valid UMIs (STL96)\
output\
    deduplicated sam\
\
\
setup\
  load STL96.txt to set umi_set\
  open input sam for reading\
  open output sam for writing\
  write all header lines, starting with @, from in to out\
  make empty set seen_keys\
\
for each line in\
  parse fields: qname (col1), flag (col2), rname (col3), pos (col4), cigar (col6)\
umi\
  umi = extract_umi(qname) and split on :\
\
  if umi is not in umi_set:\
    continue and discard and move to next read\
\
Find strand information\
  strand = extract_strand(flag)  - if flag & 0x10 else +\
\
Find true position to account for any possible soft clipping\
  true_pos = calculate_pos(pos, cigar, strand); \outl0\strokewidth0  calculate_pos() function to \outl0\strokewidth0 \strokec2 extract pos\
    if strand +, subtract leading S length from pos\
    if strand -, true_pos = pos + sum of M,D,N,=,X + trailing S length \
\
Extract chromosome information from rname\
  chrom = rname\
\
\
Keep or discard duplicates\
  check if umi is valid\
	if umi is already in list of know barcode, move on to following steps\
	if umi is not in known barcodes, discard\
\pard\pardeftab720\partightenfactor0
\cf0 \outl0\strokewidth0   make key = (umi, true_pos, strand, chrom)\outl0\strokewidth0 \strokec2 \
\pard\pardeftab720\partightenfactor0
\cf0   populate with keys \
	if key is not in the set \
		add key to set\
		add to output file\
	if key already in set\
		discard read and break to next line\
\
\
close in and out files\
\pard\pardeftab720\partightenfactor0
\cf0 \
}