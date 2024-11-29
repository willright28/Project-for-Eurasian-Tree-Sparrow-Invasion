#We used weighted-Z analysis (WZA) to identify climate-adaptive genetic loci
#First the WZA need a rank transform summary statistics for individual SNPs as input
#We used RDA to get the statistics, script can be found in https://github.com/willright28/Tibet-mammals-and-birds or refer to https://popgen.nescent.org/2018-03-27_RDA_GEA.html
#After this, we can compute the window-based statistics by running commond below, for more details, please refer to https://github.com/TBooker/WZA

python general_WZA_script.py --correlations pop.input.csv \
--summary_stat q.values \
--window class \
--output pop_WZA_output.csv \
--sep "\t" \
--min_snps 5 \
--MAF ALT_FREQS
