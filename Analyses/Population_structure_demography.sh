#Eastimate linkage disequilibrium (LD) using PopLDdecay
PopLDdecay \
-i Tree_sparrow.filter2.vcf.gz \
-s pop.list -o LD.stat.gz

#Remove SNPs at high LD
plink --vcf Tree_sparrow.filter2.vcf.gz \
--double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--indep-pairwise 20kb 1 0.1 --out Tree_sparrow.filter_LD

sed -i 's/:/\t/g' Tree_sparrow.filter_LD.prune.in

vcftools --gzvcf Tree_sparrow.filter2.vcf.gz \
--positions Tree_sparrow.filter_LD.prune.in --recode --out Tree_sparrow.filter_LD

#Principal component analysis (PCA)
plink --vcf Tree_sparrow.filter_LD.vcf \
--double-id --recode --out Tree_sparrow.filter_LD --allow-extra-chr -chr-set 41

plink --file Tree_sparrow.filter_LD --make-bed --allow-extra-chr -chr-set 41 --out Tree_sparrow.filter_LD

plink --allow-extra-chr -chr-set 41 -bfile Tree_sparrow.filter_LD --pca --out Tree_sparrow.filter_LD.pca

#Population genetic structure using R package LEA
#First, convert vcf to genotype matrix using vcftools
vcftools --vcf Tree_sparrow.filter_LD.vcf --012 --out Tree_sparrow.filter_LD
library(LEA)
library(data.table)
Tree_sparrow.filter_LD.geno <- fread("Tree_sparrow.filter_LD.012")
write.lfmm(Tree_sparrow.filter_LD.geno[,2:ncol(Tree_sparrow.filter_LD.geno)], "Tree_sparrow.filter_LD.lfmm")#the first column is ID
project = NULL
project = snmf("Tree_sparrow.filter_LD.lfmm",K = 2:4,entropy = TRUE,repetitions = 5,iterations = 200,project = "new")

#Neighbor-joining (NJ) tree using PHYLIPNEW
VCF2Dis -InPut c.vcf -OutPut Tree_sparrow.filter_LD.mat
~/PHYLIPNEW-3.69.650/src/fneighbor -datafile Tree_sparrow.filter_LD.mat \
-outfile Tree_sparrow.filter_LD.mat.txt \
-outgrno 1 \
-matrixtype s -treetype n \
-outtreefile Tree_sparrow.filter_LD.mat.nj.tre

#Infer source population for introduced population using f3-statistic
#First, convert vcf to ADMIXTOOLS input format using convert.sh written by Joana Meier
bash convert.sh Tree_sparrow.filter_LD.vcf --renameScaff
#Then, running f3-statistic using in R
library(admixtools)
prefix = 'Tree_sparrow.filter_LD'
my_f2_diR = 'Source_population'
extract_f2(prefix, my_f2_dir,overwrite = T)
f2_blocks = f2_from_precomp(my_f2_dir,afprod = TRUE)
outgroup = c("Outgroup")
introduced_pop = c("introduced_pop")#target pop, e.g., Australia or USA pop in my study
test_group=c("XX1","XX2")#list all populations you want to test
f3_run <- qp3pop(f2_blocks, outgroup, introduced_pop, test_group)#f3

#Infer recombination rate using ReLERNN
~/ReLERNN/ReLERNN_SIMULATE \
-v Tree_sparrow.filter2.vcf \
--phased \
-g SMQ.Chrplus.fasta.fai \
-d ./output \
-l 1 
~/ReLERNN/ReLERNN_TRAIN \
-d ./output 
~/ReLERNN/ReLERNN_PREDICT \
-v  Tree_sparrow.filter2.vcf \
-d ./output \
--phased

#Infer recent population history using GONE
#For more details, e.g., setting of INPUT_PARAMETERS_FILE, please refer to https://github.com/esrud/GONE
bash script_GONE.sh Tree_sparrow.gone.input