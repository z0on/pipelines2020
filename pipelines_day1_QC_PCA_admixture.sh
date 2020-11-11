Acropora millepora: Orpheus-Keppel comparison

BAMS=/work/01211/cmonstr/pipelines/*.bam
GENOME_REF=/work/01211/cmonstr/pipelines/amilV2_chroms.fasta

SAFs, ibsMat (just in case, dont download them yet):
/work/01211/cmonstr/pipelines/angsdResults/*

BEFORE STARTING, replace, in this whole file:
	- matz@utexas.edu by your actual email;
	- yourusername with your TACC user name.

# ============== installations ================================

# ------- ANGSD: 

# install xz first from https://tukaani.org/xz/

cd
wget https://tukaani.org/xz/xz-5.2.4.tar.gz --no-check-certificate
tar vxf xz-5.2.4.tar.gz 
cd xz-5.2.4/
./configure --prefix=$HOME/xz-5.2.4/
make
make install

# edit .bashrc:
cd
nano .bashrc
   export LD_LIBRARY_PATH=$HOME/xz-5.2.4/lib:$LD_LIBRARY_PATH
   export LIBRARY_PATH=$HOME/xz-5.2.4/lib:$LIBRARY_PATH
   export C_INCLUDE_PATH=$HOME/xz-5.2.4/include:$C_INCLUDE_PATH
logout
# re-login

# now, install htslib:
cd
git clone https://github.com/samtools/htslib.git
cd htslib
make CFLAGS=" -g -Wall -O2 -D_GNU_SOURCE -I$HOME/xz-5.2.4/include"

# install ANGSD
cd
git clone https://github.com/ANGSD/angsd.git 
cd angsd
make HTSSRC=../htslib

# now adding ANGSD to $PATH
cd
nano .bashrc
# section 2:
   export PATH=$HOME/angsd:$PATH
   export PATH=$HOME/angsd/misc:$PATH
# save (Ctl-O, Ctl-X)

# ---- PCAngsd

cd
module load python2
git clone https://github.com/Rosemeis/pcangsd.git
cd pcangsd/
python setup.py build_ext --inplace


#-------  stairwayPlot (skip this for now)

# project page: https://sites.google.com/site/jpopgen/stairway-plot
cdw
# get version from June 2016 (v2beta2)
wget https://www.dropbox.com/s/toxnlvk8rhe1p5h/stairway_plot_v2beta2.zip
unzip stairway_plot_v2beta2.zip
mv stairway_plot_v2beta2 stairway_plot_v2beta

# ------- Moments: 
cd
git clone https://bitbucket.org/simongravel/moments.git 
cd moments
module load python2
python setup.py build_ext --inplace

# add this to .bashrc, section 2:
  export PYTHONPATH=$PYTHONPATH:$HOME/moments
# re-login

# to see if it worked:
module load python2
python
import moments
# if you get an error message something is wrong, if you just see >>> it is all fine
quit()

# ------ installing 2bRAD scripts in $HOME/bin 

cd
mkdir bin 
cd ~/bin 
# cloning github repositories
git clone https://github.com/z0on/2bRAD_denovo.git
# move scripts to ~/bin from sub-directories
mv 2bRAD_denovo/* . 
# remove now-empty directories
rm -rf 2bRAD_denovo 


#===================== A  N  G  S  D =====================

# "FUZZY genotyping" with ANGSD - without calling actual genotypes but working with genotype likelihoods at each SNP. Optimal for low-coverage data (<10x).
# if your coverage is >10x, go to GATK section below

# install ANGSD (see InstallingSoftware.txt file... this can be a challenge, so let me know when/if you get stuck)

#----------- assessing base qualities and coverage depth

# This is our data:
BAMS=/work/01211/cmonstr/pipelines/*.bam

# switching to scratch, creating working directory
cds
mkdir pipes
cd pipes

# creating a list of bam files to work on
ls $BAMS > bams

# entering interactive session, giving all node's memory to one process:
idev -tpn 1 -N 1

FILTERS="-uniqueOnly 1 -minMapQ 20 -maxDepth 10000"
# if only looking at high-confidence  SNPs
# FILTERS="-uniqueOnly 1 -minMapQ 20 -maxDepth 10000 -snp_pval 1e-5"

# T O   D O : 
 TODO="-doQsDist 1 -doDepth 1 -doCounts 1 -dumpCounts 2"
# if only looking at high-confidence SNPs:
# TODO="-doQsDist 1 -doDepth 1 -doCounts 1 -dumpCounts 2 -doMajorMinor 1 -doMaf 1"

# in the following line, -r argument is ~1 Mb (no need to do this for whole genome)
# (look up lengths of your contigs in the header of *.sam files if you need)
angsd -b bams -r chr10:1-1000000 -GL 1 $FILTERS $TODO -P 12 -out dd

# summarizing results (using modified script by Matteo Fumagalli)
Rscript ~/bin/plotQC.R prefix=dd

# scp dd.pdf to laptop to see distribution of base quality scores and fraction of sites in each sample depending on coverage threshold

# note the new file bams.qc: it lists only the bam files that are 3 SDs less than mean quality (quality = fraction of sites with >5x coverage, written in file quality.txt)

#--------------- population structure (based on common polymorphisms, allele freq >0.05)

# --- step1: identifying weird and clonal samples:

FILTERS="-uniqueOnly 1 -minMapQ 20 -minQ 20 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -minInd 34 -snp_pval 1e-5 -minMaf 0.05"
TODO="-doMajorMinor 1 -skipTriallelic 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doPost 1 -doGlf 2"

# Starting angsd with -P the number of parallel processesors. 
echo "angsd -b bams.qc -GL 1 $FILTERS $TODO -P 12 -out OKall" >aa
ls5_launcher_creator.py -j aa -n aa -t 0:30:00 -a mega2014 -e matz@utexas.edu -w 1
cat aa.slurm | perl -pe 's/module/#SBATCH --reservation=genomics_day1\nmodule/' > aa.R.slurm
sbatch aa.R.slurm

#if you are stuck, get pre-made OKall.ibsMat here: /work/01211/cmonstr/pipelines/angsdResults/

# analyze OKall.ibsMat using OKall_ibs.R to identify clones to remove
# (also need to scp file bams.qc to laptop to do this)

# --- step 2: RERUNNING angsd after cleaning up bams list

# scp bams.nc (without clones and weirdos) from laptop to here, or get it from /work/01211/cmonstr/pipelines/angsdResults/

FILTERS="-uniqueOnly 1 -skipTriallelic 1 -minMapQ 20 -minQ 20 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -minInd 52 -snp_pval 1e-5 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doPost 1 -doGlf 2"
echo "angsd -b bams.nc -GL 1 $FILTERS $TODO -P 12 -out OK" >ab
ls5_launcher_creator.py -j ab -n ab -t 0:30:00 -a mega2014 -e matz@utexas.edu -w 1
cat ab.slurm | perl -pe 's/module/#SBATCH --reservation=genomics_day1\nmodule/' > ab.R.slurm
sbatch ab.R.slurm

#if you are stuck, get pre-made OK.ibsMat and OK.beagle.gz here: /work/01211/cmonstr/pipelines/angsdResults/

# ------- PCAngsd: estimating admixture, kinship, and SNP covariance

module load python2
python ~/pcangsd/pcangsd.py -beagle OK.beagle.gz -admix -o pcangsd -inbreed 2 -kinship -selection -threads 12

# making a table of bams : population correspondence
cat bams.nc | perl -pe 's/(.+)([OK])(\d+)(.+)/$2$3\t$2/' > inds2pops 

# transfer inds2pops, pcangsd* and OK.ibsMat files to laptop, proceed with admixturePlotting_pcangsd.R and OK_ibs.R








