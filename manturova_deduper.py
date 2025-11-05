#!/usr/bin/env python
import argparse 
import re 

parser = argparse.ArgumentParser(description="SAM file with uniquely mapped reads, remove all PCR duplicates (keep only a single copy of each read). Needs to be sorted SAM file.") 
parser.add_argument("-f", "--file", help="file path to sorted sam file", required=True)
parser.add_argument("-o", "--outfile", help="file path to output sam file", required=True)
parser.add_argument("-u", "--umi", help="file containing the list of UMIs",required=True)
args = parser.parse_args()


def five_prime_pos_funct(lmp:int, cigar:str, strand:str) -> int:
    '''Given the leftmost +1 starting position, cigar, and strand (+/-), will return the 5' starting position. Accounts for soft clipping.'''
    #check if cigar is valid. "*"" indicates unmapped read. lmp (left most position) of 0 is also unmapped, since reads are 1-based.
    if cigar == '*' or lmp == 0: 
        return 0

    if strand == "+": # account for soft clipping at start of + strand
       length = re.match(r"^([0-9]+)S", cigar)
       if length:
           soft_clip_length = int(length.group(1))
           return lmp - soft_clip_length
       else:
           return lmp

    elif strand == "-":
    # need to calculate the end position of reverse strand read
        cigar_parts = re.findall(r'(\d+)([MDNS])', cigar)
        total_length = 0

        for length_str, cigar_code in cigar_parts: 
            if cigar_code in ('M', 'D', 'N', 'S'):  # consumes reference
                total_length += int(length_str) 

        return lmp + total_length - 1  #subtract 1 to get 5' position since lmp is inclusive

    return 0  #return 0 if strand is not recognized

if __name__ == "__main__":


# Load UMIs into set
    umi_set = set()
    with open(args.umi, "r") as f:
        for line in f:
            umi = line.strip()
            if umi:
                umi_set.add(umi)

# Initialize variables
    seen_keys = set()
    chromosome_counts = {}
    current_chromosome = None
    reverse_strand_flag = 16 # SAM flag for reverse strand
    read_count = 0
    header_count = 0
    unique_count = 0
    wrong_umi_count = 0
    duplicate_count = 0

# Process the SAM file
    with open(args.file, "r") as fh, open(args.outfile, "w") as out: 
        for line in fh:
            line_stripped = line.strip()


            if line_stripped.startswith("@"):
                out.write(line) 
                header_count += 1
                continue


            fields = line_stripped.split("\t")
            if len(fields) < 11: # skip any incomplete/broken lines
                continue
            
            read_count += 1

            flag = int(fields[1])
            rname = fields[2]
            pos = int(fields[3])
            cigar = fields[5] 
            qname = fields[0]

            umi = qname.split(":")[-1]

            # skip reads with wrong UMIs
            if umi not in umi_set: 
                wrong_umi_count += 1
                continue

            # reset seen keys when chromosome changes
            if rname != current_chromosome:
                current_chromosome = rname
                seen_keys.clear() 

            # skip unmapped reads
            if pos == 0 or cigar == '*': 
                continue
                
            strand = "-" if (flag & reverse_strand_flag) else "+" #determine strand, if true then strand is "-", otherwise "+"

            five_prime_pos = five_prime_pos_funct(pos, cigar, strand) #get 5' position

            dedup_key = (five_prime_pos, strand, umi)  #create deduplication key

            if dedup_key in seen_keys: 
                duplicate_count += 1
                continue
            else:
                out.write(line) 
                seen_keys.add(dedup_key)
                unique_count += 1

                if rname in chromosome_counts:
                    chromosome_counts[rname] += 1
                else:
                    chromosome_counts[rname] = 1

    # Print out stats
    print(f"Total reads {read_count}")
    print(f"Total headers {header_count}")
    print(f"Duplicate reads removed {duplicate_count}")
    print(f"Unique reads {unique_count}")
    print(f"Wrong UMI reads {wrong_umi_count}")

    for chrom in sorted(chromosome_counts):
        print(f"{chrom}\t{chromosome_counts[chrom]}")