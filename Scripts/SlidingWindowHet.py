# Script to count number of called genotypes and number of heterozygotes per sample in 
# sliding windows.
# Input file is a single- or multi-sample VCF file that has been filtered (passing sites 
# have "PASS" in the FILTER column) and compressed with gzip/bgzip.
#
# Usage: 
# python ./SlidingWindowHet.py [vcf] [window size] [step size] [chromosome/scaffold name] [output]
#
# Windows will be non-overlapping if step size == window size.
#
# Example: 
# python ./SlidingWindowHet.py input.vcf.gz 100000 10000 chr01

import sys
import pysam
import os
import gzip


# Open input file and make sure the VCF file is indexed (if not, create index)
filename = sys.argv[1]
VCF = gzip.open(filename, 'rt')

if not os.path.exists("%s.tbi" % filename):
    pysam.tabix_index(filename, preset="vcf")
parsevcf = pysam.Tabixfile(filename)


# Set variables
window_size = int(sys.argv[2])
step_size = int(sys.argv[3])
chrom = sys.argv[4]
output = sys.argv[5]


# Generate a dictionary with chromosomes and chromosome lengths
# File "chrom_lengths.txt" is a two-column tab-delimited list of chromosomes and 
# chromosome lengths in the reference genome
cc=open("chrom_lengths.txt", 'r')
chrom_size={line.strip().split('\t')[0]:line.strip().split('\t')[1] for line in cc}
cc.close()


# Get list of samples from VCF file header
samples=[]
for line in VCF:
    if line.startswith('##'):
        pass
    else:
        for i in line.split()[9:]: samples.append(i)
        break


# Get start and end positions of chromosome
#for line in VCF:
 #   if line[0] != '#':
  #      start_pos = int(line.strip().split()[1])
   #     end_pos = int(chrom_size[chrom])
    #    break
start_pos = 1
end_pos = int(chrom_size[chrom])

# Create output file
output = open(output + '_' + chrom + '_het_%swin_%sstep.txt' % (window_size, step_size), 'w')
output.write('chrom\twindow_start\twindow_end\tsites_total\tcalls_%s\thets_%s\n' % ('\tcalls_'.join(samples), '\thets_'.join(samples)) )


# Fetch a region, ignore sites that fail filters, tally genotype calls and heterozygotes        
def snp_cal(chrom,window_start,window_end):
    #print("%s:%s" % (chrom,window_start))
    rows = tuple(parsevcf.fetch(region="%s:%s-%s" % (chrom, window_start, window_end), parser=pysam.asTuple()))    
    sites_total=0
    calls=[0]*len(samples)
    hets=[0]*len(samples)
    hetRate=[0]*len(samples)
    for line in rows:
        if line[6]!="PASS": continue
        sites_total+=1
        for i in range(0,len(samples)):
            if line[i+9][:1]=='.': continue
            calls[i]+=1
            GT=line[i+9].split(':')[0]
            if '/' in GT: sp='/'
            if '|' in GT: sp='|'
            if GT.split(sp)[0]!=GT.split(sp)[1]: hets[i]+=1
    output.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (chrom,window_start,window_end,sites_total,'\t'.join(map(str,calls)),'\t'.join(map(str,hets))) )


# Initialize window start and end coordinates
window_start = start_pos
window_end = start_pos+window_size-1


# Calculate stats for window, update window start and end positions, 
# repeat to end of chromosome
while window_end <= end_pos:    
    if window_end < end_pos:
        snp_cal(chrom,window_start,window_end)
        window_start = window_start + step_size
        window_end = window_start + window_size - 1
    else:
        snp_cal(chrom,window_start,window_end)
        break    
else:
    window_end = window_start + window_size - 1
    snp_cal(chrom,window_start,window_end)


# Close files and exit
VCF.close()
output.close()

exit()



