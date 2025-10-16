Process sorted SAM file, which contains uniquely mapped single end reads, and remove PCR duplicated
PCR duplicated contain the same RNAME, POS, FLAG, and UMI
Keep the first unique read encountered and remove any subsequent reads.

input
    sam file 
    file with valid UMIs (STL96)
output
    deduplicated sam


setup
  load STL96.txt to set umi_set
  open input sam for reading
  open output sam for writing
  write all header lines, starting with @, from in to out
  make empty set seen_keys

for each line in
  parse fields: qname (col1), flag (col2), rname (col3), pos (col4), cigar (col6)
umi
  umi = extract_umi(qname) and split on :

  if umi is not in umi_set:
    continue and discard and move to next read

Find strand information
  strand = extract_strand(flag)  - if flag & 0x10 else +

Find true position to account for any possible soft clipping
  true_pos = calculate_pos(pos, cigar, strand);  calculate_pos() function to extract pos
    if strand +, subtract leading S length from pos
    if strand -, true_pos = pos + sum of M,D,N,=,X + trailing S length 

Extract chromosome information from rname
  chrom = rname


Keep or discard duplicates
  check if umi is valid
	if umi is already in list of know barcode, move on to following steps
	if umi is not in known barcodes, discard
  make key = (umi, true_pos, strand, chrom)
  populate with keys 
	if key is not in the set 
		add key to set
		add to output file
	if key already in set
		discard read and break to next line


close in and out files

