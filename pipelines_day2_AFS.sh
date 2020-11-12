# note: --reservation argument is now genomics_day2 (it will be day3 tomorrow!)
# first, we need to identify sites to work with and estimate SFS (site frequency spectrum) for each population.

# Running ANGSD on all bams with filters that do not disturb AFS, to find sites to work on.
FILTERS="-uniqueOnly 1 -skipTriallelic 1 -minMapQ 30 -minQ 30 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -maxHetFreq 0.5 -minInd 38 "
TODO="-doMajorMinor 1 -doMaf 1 -dosnpstat 1 -doPost 2 -doGeno 11"
echo "angsd -b bams.nc -GL 1 -P 12 $FILTERS $TODO -out oksites">ac
ls5_launcher_creator.py -j ac -n ac -t 0:30:00 -e matz@utexas.edu -w 1 -a mega2014 -q normal
cat ac.slurm | perl -pe 's/module/#SBATCH --reservation=genomics_day2\nmodule/' > ac.R.slurm
sbatch ac.R.slurm

# saving and indexing sites that pass the filters (excluding mitochondrion)
zcat oksites.mafs.gz | cut -f 1,2 | tail -n +2 | grep "chr" > goodsites
angsd sites index goodsites

# saving bam lists corresponding to each population
grep O bams.nc >o
grep K bams.nc >k

# generating per-population SAF (site allele frequency likelihoods)
export GENOME_REF=/work/01211/cmonstr/pipelines/amilV2_chroms.fasta
TODO="-doSaf 1 -doMajorMinor 1 -doMaf 1 -doPost 1 -anc $GENOME_REF -ref $GENOME_REF"
echo "angsd -sites goodsites -b o -GL 1 -P 4 $TODO -out o
angsd -sites goodsites -b k -GL 1 -P 4 $TODO -out k">sfsa
ls5_launcher_creator.py -j sfsa -n sfsa -t 0:30:00 -e matz@utexas.edu -w 2 -a mega2014
cat sfsa.slurm | perl -pe 's/module/#SBATCH --reservation=genomics_day2\nmodule/' > sfsa.R.slurm
sbatch sfsa.R.slurm

# ----- generating bootstrapped SFS data (resampling chromosomes)

# first generate 100 series of 5 bootstraps for each population:
export GENOME_REF=/work/01211/cmonstr/pipelines/amilV2_chroms.fasta
>b100
for B in `seq 1 100`; do
echo "sleep $B && realSFS o.saf.idx k.saf.idx -ref $GENOME_REF -anc $GENOME_REF -bootstrap 5 -P 1 -resample_chr 1 >ok_$B">>b100;
done
ls5_launcher_creator.py -j b100 -n b100 -t 1:30:00 -e matz@utexas.edu -w 24 -a mega2014
cat b100.slurm | perl -pe 's/module/#SBATCH --reservation=genomics_day2\nmodule/' > b100.R.slurm
sbatch b100.R.slurm

# then, do "bagging" (averaging of 5 bootstrap replicates within each of the 100 series)
SFSIZE="73 53" # 2N+1 for each population.
for B in `seq 1 100`; do
echo $SFSIZE >ok_${B}.sfs;
cat ok_${B} | awk '{for (i=1;i<=NF;i++){a[i]+=$i;}} END {for (i=1;i<=NF;i++){printf "%.3f", a[i]/NR; printf "\t"};printf "\n"}' >> ok_${B}.sfs;
done

# ------ get pre-computed bootstrapped sfs (only if you are stuck with the above):
 
cp /work/01211/cmonstr/pipelines/sfs/* .

# ----- plotting 2dAFS (try different projections)

# full AFS
2dAFS.py ok_1.sfs O K 72 52
# to 0.9
2dAFS.py ok_1.sfs O K 66 48
# to 0.8
2dAFS.py ok_1.sfs O K 58 42

# transfer 2dAFS_ok_1.sfs_O_K* files to laptop for viewing


# ----- searching for best-fitting model

# (see README at https://github.com/z0on/AFS-analysis-with-moments for explanations)

# cloning Misha's Moments model collection and accessory scripts
cd
rm -rf AFS-analysis-with-moments/
git clone https://github.com/z0on/AFS-analysis-with-moments.git
cp ~/AFS-analysis-with-moments/multimodel_inference/py2/* ~/bin/
cd -

# writing a huge list of model runs
module load Rstats/3.5.1
Rscript ~/AFS-analysis-with-moments/modSel_write.R contrast=ok nreps=3 args="o k 58 42 0.02 0.005"
# NOTE: do not run this! (it would take too much computer time for all of us). Use pre-computed result ok.modsel

# summarizing results, writing commands to bootstrap the winning model 
module load Rstats/3.5.1
Rscript ~/AFS-analysis-with-moments/modSel_summary.R modselResult=ok.modsel args="o k 58 42 0.02 0.005"

# ----- bootstrapping the winning model

# The previous script, modSel_summary.R, produced the file ok.winboots.runs containing a list of commands to run - same model on 'nboots' bootstrap replicates. We need to launch these commands all in parallel. 
# (how do we do it? write it yourself!)

ls5_launcher_creator.py -j ok.winboots.runs -n ok.winboots.runs -t 1:00:00 -e matz@utexas.edu -w 24 -a mega2014 -q normal
cat ok.winboots.runs.slurm | perl -pe 's/module/#SBATCH --reservation=genomics_day2\nmodule/' > ok.winboots.runs.R.slurm
sbatch ok.winboots.runs.slurm


# the result will be collected in ok.winboots. 
# we must concatenate all these results - please email your ok.winboots to Misha, matz@utexas.edu, with subject "ok.winboots". Wait for Misha to email you back the concatenated file.
# To summarize and plot results, do
module load Rstats/3.5.1
Rscript ~/AFS-analysis-with-moments/bestBoot_summary.R bootRes=ok.winboots


#--------------- GADMA

#installing dadi
cd
git clone https://bitbucket.org/gutenkunstlab/dadi.git
cd dadi
python setup.py install --user

# installing Pillow
python -m pip install Pillow

# installing GADMA
cd
rm -rf GADMA
git clone https://github.com/ctlab/GADMA.git
cd GADMA
python setup.py install --user

# writing GADMA parameters file

mv ok_1.sfs_20_20.projected ok_1_20_20.projected.sfs

echo "# param_file
Input file : ok_1_20_20.projected.sfs
Output directory : gadma_ok
Population labels : o , k
Initial structure : 1,1" >gadma_params

echo "gadma -p gadma_params">gad
#echo "gadma --resume gadma_12dp">gad
ls5_launcher_creator.py -j gad -n gad -t 14:00:00 -e matz@utexas.edu -w 1 -a tagmap -q normal
sbatch gad.slurm


