#Estimate genetic diversity using ANGSD
~/angsd/angsd -bam pop_bam.list \
-ref SMQ.Chrplus.fasta \
-out pop -doSaf 1 -gl 1 -minInd `wc -l pop_bam.list | awk '{print $1}'` -minMapQ 30 -minQ 20 

~/angsd/misc/realSFS pop.saf.idx -maxIter 100 > pop.sfs 

~/angsd/misc/realSFS saf2theta pop.saf.idx -outname pop -sfs pop.sfs

~/angsd/misc/thetaStat do_stat pop.thetas.idx -win 20000 -step 20000 -outnames pop.window.theta

#Estimate genome-wide heterozygosity using script from Kyriazis et al. (2023), https://github.com/ckyriazis/moose_WGS_project
python SlidingWindowHet.py Tree_sparrow.filter2.pop.vcf.gz 20000 20000 chr_name pop

#Estimate fst using vcftools
vcftools --gzvcf Tree_sparrow.filter2.vcf.gz --weir-fst-pop pop1.list --weir-fst-pop pop2.list --fst-window-size 20000 --fst-window-step 20000 --out pop1_2

#Estimate runs of homozygosity (ROH)
plink --vcf Tree_sparrow.filter2.vcf.gz --double-id --allow-extra-chr -chr-set 41 --make-bed --out Tree_sparrow.ROH

plink --bfile Tree_sparrow.ROH --homozyg --homozyg-window-snp 50 --homozyg-density 50 --homozyg-gap 200 \
 --homozyg-window-het 2 --homozyg-window-missing 5 --homozyg-kb 350 --allow-extra-chr --chr-set 41 --out Tree_sparrow.ROH.350kb

plink --bfile Tree_sparrow.ROH --homozyg --homozyg-window-snp 50 --homozyg-density 50 --homozyg-gap 200 \
 --homozyg-window-het 2 --homozyg-window-missing 5 --homozyg-kb 1000 --allow-extra-chr --chr-set 41 --out Tree_sparrow.ROH.1mb

#Estimate genetic load
java -jar ~/snpEff/snpEff.jar -v \
-lof \
-c ~/snpEff/snpEff.config Tree_sparrow_newGenome \
Tree_sparrow.filter2.ancestry.vcf.gz > Tree_sparrow.anc.ann.vcf

cat Tree_sparrow.anc.ann.vcf | java -jar ~/snpEff/SnpSift.jar filter "(exists LOF)" > Tree_sparrow.anc.ann_lof.vcf 

java -jar ~/snpEff/SnpSift.jar filter "ANN[*].EFFECT = 'missense_variant'"  Tree_sparrow.anc.ann.vcf > Tree_sparrow.anc.ann_missense.vcf 

java -jar ~/snpEff/SnpSift.jar filter "ANN[*].EFFECT = 'intergenic_region'"  Tree_sparrow.anc.ann.vcf.gz > Tree_sparrow.anc.ann_intergenic.vcf 

#Infer balance selection using BetaScan2
#For more details, please refer to https://github.com/ksiewert/BetaScan
~/glactools/glactools vcfm2acf --onlyGT --fai SMQ.Chrplus.fasta.fai Tree_sparrow.filter2.pop.vcf.gz > acf/pop.acf.gz 

~/glactools/glactools acf2betascan --fold acf/pop.acf.gz | gzip > beta/pop.beta.txt.gz

python BetaScan.py -i beta/pop.beta.txt.gz -fold -w 20000 -m 0.05 -o output/pop.betascores.txt

