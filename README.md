# RCADE
Recognition Code-Assisted Discovery of regulatory Elements (http://rcade.ccbr.utoronto.ca/)

Motif discovery from ChIP-seq data is often limited by presence of non-targeted transcription factor motifs, as well as similarity of peak sequences due to common ancestry rather than common binding factors. The latter aspect particularly affects a large number of proteins from the Cys2His2 zinc finger (C2H2-ZF) class of transcription factors, as their binding sites are often dominated by endogenous repeat elements (EREs) that have highly similar sequences. To overcome these limits, RCADE combines predictions from a DNA recognition code of C2H2-ZFs with ChIP-seq data to identify models that represent the genuine DNA binding prefer-ences of C2H2-ZF proteins.


## Requirements
- Unix-compatible OS
- R version 3.0.1 or later (http://www.r-project.org/)
- R “randomForest” library (http://cran.r-project.org/web/packages/randomForest/index.html)
- GNU-compatible MAKE software (https://gcc.gnu.org/)
- MEME Suite (http://meme.nbcr.net/meme/downloads.html)
 
## Installation
- Step 1. To install the program, extract the package, and run the "make" command.
- Step 2. Change the value of line 7 of the “RCOpt.sh” script to where the executable MEME files are located on your computer.
 
To test the pipeline, execute this command:

> bash RCOpt.sh MyTestJob examples/CTCF/CTCF.fasta examples/CTCF/GSM1407629.500bp.fasta

This should create a “./out/MyTestJob” folder, with the RCADE output files described below.
 
## Usage
Use the RCOpt.sh script to run RCADE on your dataset:

> bash RCOpt.sh jobName fastaC2H2 fastaChIP

The argument _jobName_ is a unique identifier for your job. The output files of RCADE will be placed in "./out/_jobName_". These files will include:

- **results.ps**: A postscript file that visualizes a summary of the optimization results. RCADE identifies several motifs from the ChIP-seq data, which are sorted in this file based on their AUROC values for distinguishing ChIP-seq peaks from dinucleotide-shuffled sequences. For each motif, the corresponding zinc fingers are shown on the top (for example, CTCF:3-7 means that zinc fingers 3-7 of the CTCF protein are used for predicting the initial seed motif that is then optimized). The seed motif that is directly predicted from protein sequence is then shown, followed by the motif that is optimized based on ChIP-seq data. The AUROC value for each motif, the associated p-value, as well as the Pearson similarity of the seed and optimized motifs are also shown.
- **results.opt.ps**: Same as the above output, except that it only includes the top-scoring optimized motif.
- **results.opt.PFM.txt**: A text file containing the PFM of the top-scoring optimized motif, in a format similar to what is used in the CisBP database (http://cisbp.ccbr.utoronto.ca/).
- **results.opt.PFM.meme.txt**: A text file containing the PFM of the top-scoring optimized motif, in a format suitable for the MEME suite (http://meme.nbcr.net/meme/).
- **results.PFM.txt**: A text file containing all seed motifs and their optimized versions (the optimized motif names end with the phrase “opt”). The motifs are in CisBP format.
- **results.report.txt**: A report table, summarizing the optimization results for the motifs.
- **log.info.txt**: A short summary of warning/error/info messages.
