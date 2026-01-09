# ZPeaks_U

The program ZPeaks_U calls peaks from a genomic experiment, comparing experiment and control over optimized windows and taking into account the biases associated with nucleotide analogs such as BrdU or EdU

Author: Ugo Bastolla <ubastolla@cbm.csic.es>

Centro de Biologia Molecular Severo Ochoa (CSIC-UAM), Madrid Spain

Associated paper: Correcting the biases of nascent strand detection with thymidine analogues sheds light on the replication time of the Arabidopsis thaliana genome and its coordination with gene function. by Ugo Bastolla and Crisanto Gutierrez (submitted)

## OVERVIEW:

ZPeaks_U smooths the experiment and the control over windows whose size is optimized in such a way to maximize the discrimination power of the method, and it calculates their ratio, which constitutes the ZPeaks score.

ZPeaks_U is particularly designed for detecting genome replication origins (ORIs) with Short Nascent Strand (SNS) labelled with Thymidine analogoues such as BrdU or EdU, since it addresses two biases that are frequently overlooked with this method.

Firstly, T-analogues are incorporated only where a Thymidine is present, which implies that the method is biased towards T-rich regions. Zpeaks_U corrects this bias upon demand by weighting the control with the local T content. Secondly, we expect that the DNA polymerase has preferential affinity for the natural nucleotide than for the analogue, which implies that the rate of T-analogue incorporation increases with time as the replication process consumes the natural nucleotide. To take into account this factor, Zpeaks_U considers local maxima of the ZPeaks score and compares their value with the local average of the score. In this way, it classifies the candidate ORIs into five classes that we correlate with the replication time.

## Installation (linux):
====================

Download ZPeaks.zip to your computer and execute the following commands:

```sh
>mkdir <dir-name> (create a directory where ZPeaks_U will be stored)
>mv ZPeaks_U.zip <dir-name>
>cd <dir-name>
>unzip ZPeaks_U.zip
>make 
>cp ZPeaks <your-path-directory>
```

## Configuration file:

You have to build a configuration file structured as the sample file "Input_ZPeaks_example.in" provided in the package.

### Input:
The configuration file must specify the files that contain the experimental data in WIG format (EXPER record in configuration file), the data used for control (CONTROL record in configuration file; if the control is absent the program will use averaged experimental reads as a control), the genome of the organism (GENOME record in configuration file) if the user requires to weight the control with its A+T content (option -at_norm), an optional BED format file containing the peaks that you want to compare with reference peaks used for comparison (PREDICTION record in configuration file) and profiles (PROF) of genomic or epigenomic data measured across the chromosomes in GR or WIG formats (PROF parameter in configuration file).
All types of data can be contained in a single file or in a list of files, one for each chromosome. Both individual and groups of files must be terminated with the record END at the beginning of a line.

GR format:      https://stackoverflow.com/questions/28880086/what-is-gr-file-format

BED format:     https://genome.ucsc.edu/FAQ/FAQformat.html#format1

WIG format:     https://genome.ucsc.edu/FAQ/FAQformat.html#format6

### Example of input

```sh
EXPER:
### BrdU:
DIR=/home/ubastolla/RESEARCH/COLLABORATIONS/ARABIDOPSIS/ZPEAKS_U/INPUT/
BrdU_exp_3mismatch_uniq_nosibs_chr1.gr
BrdU_exp_3mismatch_uniq_nosibs_chr2.gr
BrdU_exp_3mismatch_uniq_nosibs_chr3.gr
BrdU_exp_3mismatch_uniq_nosibs_chr4.gr
BrdU_exp_3mismatch_uniq_nosibs_chr5.gr
END
CONTROL:
DIR=/home/ubastolla/RESEARCH/COLLABORATIONS/ARABIDOPSIS/ZPEAKS_U/INPUT/
# Not normalized:
BrdU_contr_3mismatch_uniq_nosibs_chr1.gr
BrdU_contr_3mismatch_uniq_nosibs_chr2.gr
BrdU_contr_3mismatch_uniq_nosibs_chr3.gr
BrdU_contr_3mismatch_uniq_nosibs_chr4.gr
BrdU_contr_3mismatch_uniq_nosibs_chr5.gr
# AT normalized:
#BrdU_contr_ATnorm_3mismatch_uniq_nosibs_chr1.gr
#BrdU_contr_ATnorm_3mismatch_uniq_nosibs_chr2.gr
#BrdU_contr_ATnorm_3mismatch_uniq_nosibs_chr3.gr
#BrdU_contr_ATnorm_3mismatch_uniq_nosibs_chr4.gr
#BrdU_contr_ATnorm_3mismatch_uniq_nosibs_chr5.gr
END
GENOME:
DIR=/home/ubastolla/RESEARCH/COLLABORATIONS/ARABIDOPSIS/ZPEAKS_U/INPUT/
#DIR=/ngs/evo/Replication_Oris/TAIR9_TRANSLATED/FASTA_SEQUENCES/
TAIR10_chr1.fas
TAIR10_chr2.fas
TAIR10_chr3.fas
TAIR10_chr4.fas
TAIR10_chr5.fas
END
NORM_AT=1 ! Normalize control by T count
PREDICTION:
DIR=/home/ubastolla/RESEARCH/COLLABORATIONS/ARABIDOPSIS/ZPEAKS_U/INPUT/
Zpeaks_4d_10d_J200_S200.bed
#Zpeaks_BrdU_Norm_AT_T2_W1500_J0.bed
END
PROF CDC6
DIR=/home/ubastolla/RESEARCH/COLLABORATIONS/ARABIDOPSIS/ZPEAKS_U/INPUT/
CDC6_34_vs_mpTiling_COL_123_original_CEL_Unique_PM_hmm_chr1.gr
CDC6_34_vs_mpTiling_COL_123_original_CEL_Unique_PM_hmm_chr2.gr
CDC6_34_vs_mpTiling_COL_123_original_CEL_Unique_PM_hmm_chr3.gr
CDC6_34_vs_mpTiling_COL_123_original_CEL_Unique_PM_hmm_chr4.gr
CDC6_34_vs_mpTiling_COL_123_original_CEL_Unique_PM_hmm_chr5.gr
END
```

## Detailed description of the algorithm

The program Zpeaks_U calls putative ORIs by comparing the number of reads in each genomic bin for an experiment of short nascent strand (SNS) sequencing and its control, which are input as wig files. As its predecessor program Zpeaks (Sequeira-Mendes et al. 2019), Zpeaks_U detects ORIs as local maxima of the difference between the normalized experiment and control and it smooths this difference across neighboring windows according to parameters that it optimizes internally. Zpeaks_U adds new features that are specific to the case of thymidine analogues:

(I) Under demand (option -at_norm from command line or NORM_AT=1 in configuration file), ZPeaks_U weights the control of each bin multiplying it times the number of Thymine nucleotides that it contains, in order to take into account that sequences with more Ts are more likely to be selected in the experiment.

(II) ZPeaks_U considers the possibility that the thymidine analogue, which we indicate with TA, and the natural thymine nucleotide may be incoporated by the DNA polymerases with different rates. In particular, we expect that the incorporation rate is larger for the natural T nucleotide and that the rate of incorporation of the TA increase as the DNA synthesis proceeds, since the pool of T decreases and the relative weight of TA/T increases. The program compares the Zpeaks score at a local maximum with the Zpeaks score of a region of range $R$ outside the maximum, excluding regions that we attribute to peaks, as we explain in detail below.

### Flow of the algorithm

ZPeaks_U works with the following steps, which can be modified by the user:

1) Under demand (option -at_norm from command line or NORM_AT=1 in configuration file), ZPeaks_U weights the control such that the control score of each bin is proportional to the number of Thymine nucleotides that it contains, in order to take into account that sequences with more Ts are more likely to be selected in the experiment. This is done by reading the genome sequence, determining the maximum number of Thymidine bases of the leading and lagging strand of each bin, and multiplying the number of reads in the control times this number.

2) ZPeaks_U normalizes experiment and control separately, either chromosome by chromosome (SEPARATE_NORM=1) or for the whole genome, so that their mean value across each chromosome equals 1, and it computes their difference $y(i)=(\mathrm{experiment}(i)-\mathrm{control}(i)$ at each bin $i$.

3) ZPeaks_U smooths the difference signal $y(i)$ over windows of $2l+1$ genomic bins, with weight 1 for the central bin, reduced by a factor $\exp(-\mathrm{B/DAMP})$ at each neighboring bin, until a threshold weight EPS is reached. Here $B$ is the bin size and the smoothing window size is calculated as WIN=EPS*DAMP/B. The parameters DAMP and EPS are internally optimized as explained below. The smoothed score is transformed into a Z-score, using the total standard deviation excluding outliers that deviate more than OUTL total standard deviations from the mean. If OUTL<0, outliers are not removed. The smoothing is performed either over all the genome or chromosome by chromosome (SEP_SMOOTH=1). The Z score of the smoothed difference score is called the Zpeaks score and indicated below as $\tilde{y}(i)$ .

4) Subsequently, ZPeaks_U identifies the local maxima $k$ of the Zpeaks score such that $\tilde{y}(i)$ monotonically decreases for $i>k$ and for $i<k$ for at least two bins (i.e., it has a bell shape over at least five bins). The size at which the decrease is monotonic in both directions defines the width of the peak, $W$.

5) Then ZPeaks_U computes the local mean and standard deviation of the Zpeaks score across a genomic region of range $R$ around the peak $k$, excluding bins with scores $\tilde{y}(i)$ similar to those within the candidate peak, since these bins can correspond to neighboring ORIs, i.e. it discards bins whose value is larger than the threshold $t \tilde{y}_\mathrm{max}(k) +(1-t)\tilde{y}_\mathrm{min}(k)$, where $t\in [0,1]$ is a parameter that combines the minimum and maximum value of $\tilde{y}(i)$ at the candidate peak.
To avoid being too permissive, the standard deviation $\sigma_k$ associated to the peak $k$ is not smaller than the standard deviation across the whole genome.

6) Finally, ZPeaks_U computes the local Z score around the peak $k$ using the local mean and standard deviation of the Zpeaks score and it calls the peak if the local Z score is larger than a threshold.
For improving the detection of peaks in difficult cases, the threshold $T(W)$ decreases with the peak width $W$ from the value $T_\mathrm{max}$ to $T_\mathrm{min}$, i.e. the program is more tolerant for wider peaks: $T(W)=\mathrm{max}(T_\mathrm{max}-(W-2)T_\mathrm{step},T_\mathrm{min})$
We include into the candidate ORI all the bins that surround the peak k that fulfil  $Z_\mathrm{local}(i)=(\tilde{y}(i)-<\tilde{y}>(k,R,t))/\sigma_k > T(W)$.

7) The sum of the local Z score over all bins in the peak constitutes the score of the peak Score$(k)$.

### Parameter determination

The program Zpeaks_U determines internally the smoothing parameters $DAMP$ and $EPS$ by maximizing the sum of the peak scores described at point 7) through iterative quadratic fits. For all cases that we considered these computations converge rapidly and robustly.

Another crucial parameter is the range $R$ across which the program computes the local mean and standard deviation of the score. The program Zpeaks_U determines $R$ by maximizing the mean score of the bins that belong to peaks multiplied times the square root of the number of called peaks $N_p$, i.e.

Discriminative_score=$\sqrt{N_p}\sum_k \mathrm{Score}(k)}/(N_p <m_p>$

where $<m_p>$ is the mean number of bins in a peak. We call this quantity discriminative score since it quantifies the discriminative power for distinguishing peaks from non-peaks, whose score is zero. Since the standard deviation of the Z score is one and the called peaks are independent, the discriminative score estimates the standard error of the mean difference of the scores between peaks and non-peaks.

The program also allows, upon request, to determine $R$ by maximizing the total score or maximizing the number of peaks, but the previous criterion is preferrable because it trades-off large mean score and large number of peaks.

We use as default the thresholds T_max=2.0, T_min=1.5, T_step=0.5. These parameters can be modified by the users in the configuration file. We also use $t=0.15$ for computing the local mean and standard deviation of the Zpeaks score around the candidate peaks, omitting values comparable with the cadidate peak.

### BrdU (EdU) classes

Finally, the program discretizes the mean local score $<\tilde{y}>(k, R, t)$ in five discrete values and associates each peak to one of the five classes. We expect that these classes represent different replication times, in the order -1,0,1,2,-2 from earlier to later replication time.

### FORMAT of the configuration file:

#### Input and output
NAME Zpeaks_U_BrdU_Norm_AT  ! Name of output file
EXPER:
DIR=<path/to/exper/files>
<exper_file_1> (either one file or one file for each chromosome in WIG format, MANDATORY)
<exper_file_2> (OPTIONAL)
... 
<exper_file_n> (OPTIONAL)
END
CONTROL:
DIR=<path/to/control/files>
<control_file> (either one file or one file for each chromosome, OPTIONAL)
END
DIR=<path/to/genome>
<genome_file_1> (OPTIONAL)
... 
<genome_file_n> (OPTIONAL)
END
PREDICTION:
DIR=<path/to/reference/peaks>
<reference/peaks> (just one file, bed format, OPTIONAL)
END
PROF1 <profile_name_1> (chromatin profile such as histones, histone modification...) 
DIR=<path/to/prof_1>
<prof_1_file_1> (either single file or one for each chromosome, OPTIONAL)
... 
<prof_1_file_n> (OPTIONAL)
END
PROF2 <profile_name_2>
DIR=<path/to/prof_2>
<prof_2_file_1> (either single file or one for each chromosome, OPTIONAL)
... 
<prof_2_file_n> (OPTIONAL)
END
....
#### Parameters
PRINT_SCORE=1    ! Print Peakscore as wig file
PRINT_CLASS=1    ! Print bed file for each class?
#### Parameters
NORM_AT=1        ! Weight the control with AT content
WP_MIN=2 	 ! Min.tested value of w2 def: 4
THR=2.0          ! Max. threshold on score
THR_STEP=0.0	! Decrease of threshold per unit of WP
THR_MIN=1.5	! Minimum value of Thr for WP>=WP_MIN
SCORE_RANGE=s  ! Score to optimize range: S (score) N (n.peaks) D (S/sqrt(N)) s (S/N)
RANGE_MIN=20	 ! Minimum range for local average of score
RANGE_MAX=100          ! Range for local average of score
N_RANGE=20       ! Trials for optimizing RANGE
Y_MAX=0.0 	! Omit bin if y>y_min+Y_MAX*(y_max-y_min)
LMIN_e=3         ! Min length for excluding putative peak
LMIN_p=3         ! Min length for accepting putative peak
SCORE=1          ! Optimize window to maximize number of selected bins (SCORE=0)
SEPARATE_NORM=1  ! Normalize experiment and control separately for each chr
SEP_SMOOTH=1     ! Smooth (exper-control) separately for each chr?
OUTL=3           ! Normalize smoothed score excluding outliers such that
                 ! (score-ave)/sd>OUTL (if OUTL<=0, no outlier is removed)
LOCAL=1		! Compute standard deviation locally
EXP_MIN= 1.5	 ! Min. reads per bp at peaks def: 2

## Execution

```sh
>ZPeaks_U -file <config_file>
Options:
          -at_norm  ! Normalize the control weighting it with AT content 
```

## Output:


The Zscore output comprises several files:

<name>_score.wig                            Zscore in wig format
<name>_ALL_<Parameters>.bed     Called peaks in bed format, where <Parameters>=WIN<optimized_window_size>_T<THR>
Properties_<name>.dat                       Table in text format with the input genomic and epigenomic properties of every called peak, if these properties are provided

Optional output: 

If a reference set of peaks is input as PREDICTION, the program outputs the list of peaks that overlap with the reference (<name>_OLD_<Parameters>.bed),
that do not overlap with it (<name>_NEW_<Parameters>.bed) and the reference peaks that are not found in the current experiment (<name>_notfound_<Parameters>.bed)

