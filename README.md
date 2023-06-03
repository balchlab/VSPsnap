# VSPsnap  [![DOI](https://zenodo.org/badge/472567724.svg)](https://zenodo.org/badge/latestdoi/472567724)
VSPsnap is a collection of R and Python code for Gaussian Process regression in a kriging-like setting (i.e. two features (X,Y) and a target (Z)). Includes functions for:
- GP hyperparameter search and tuning;
- Optimization;
- GP map and tabular output and plot (GP map) generation;
- Parallelization to multiple compute nodes in a HPC setting.  

Primitive kriging functions are provided by the _gstat_ R library.

VSPsnap was used to build the early warning system for SARS-CoV2 variants of concern described in “Understanding the Host-Pathogen Evolutionary Balance through Gaussian Process Modelling of SARS-CoV-2” (Loguercio et al., 2023).

### Installation

Minimal requirements:
- R (3.5 and above)
- *(Optional, for co-occurrence)*: Python3 (3.7 and above)
- *(Optional, for HPC runs)*: SLURM scheduler (20 and above)

1. Clone the repository; unzip it;
2. Generate / update the input data (see `Data/IR_FR_data_update.docx`);
3. Run VSPsnap in batch mode (see below). 

### Quick Start

Generate a GP regression map (annotated with a single VOC (Omicron), for a specific date):

*Script to use:* `VSPsnap_one/custom/VSPsnap_IR_FR_log_bestModel_one_Omicron.R`
```R
R CMD BATCH --no-restore --no-save VSPsnap_IR_FR_log_bestModel_one_Omicron.R
```

### Batch Mode

The single *VSPsnap* function (as in the **Quick Start**) processes a single date and a specific set of variant labels at the time (i.e. Omicron only; Delta only, all VOCs etc.). To generate maps for multiple dates and/or multiple sets of variant labels, a convenient way to proceed is generate batches of VSPsnap scripts, each looping through a small interval of dates. An example for generating a batch of VSPsnap scripts is given below, where a file of date intervals (with 3 intervals and 32 days in total) is used to generate three batch VSPsnap scripts, each processing ~10 days). 
1. Create a directory for the batch scripts to be generated:
``` bash
mkdir update_01_16_22
```
2. Start R and generate the batch VSPsnap scripts:
``` R
library("readtext")
setwd("update_01_16_22/")
batches<-read.table("VSPsnap/data/batchesCMT_update_01_16_22.tsv",header=F,quote="",sep=" ")
# read in template VSPsnap function to be used
tmp<-readtext("../VSPsnap/VSPsnap_IR_FR_log_bestModel_one_update_01_16_22_template.R")

# split the text in two chunks of 5K characters. This workaround is used b/c sprintf (used below for substitution) has a max character limit of 5K. 
tmp0<-substr(tmp$text,1,5000)
tmp1<-substr(tmp$text,5001,nchar(tmp$text))

for(i in 1:nrow(batches)){

  out0<-sprintf(tmp0,batches[i,1],batches[i,2])
  out<-paste0(out0,tmp1)

  cat(out,file=paste0("VSPsnap_ir_fr_log_one_011622_",i,".R"))
}
```

3. Run the batch scripts.

If a cluster of compute nodes is available, a convenient way to run the batch VSPsnap scripts thus generated is through a simple job spawner. The line below assembles the shell command to run in order to launch the parallel jobs (assuming SLURM as scheduler - job parameters can be inspected/modified through the `run_Rscript_slurm.sh` script.
```R
spawner_cmd_log_one<-paste(c("./foo_spawner_slurm.sh","run_Rscript_slurm.sh",dir(pattern = "VSPsnap_ir_fr_log_one_011622_[^(_)]")),sep="",collapse=" ")
```

```bash
./foo_spawner_slurm.sh run_Rscript_slurm.sh VSPsnap_ir_fr_log_one_011622_1.R VSPsnap_ir_fr_log_one_011622_2.R VSPsnap_ir_fr_log_one_011622_3.R
```

If the compute nodes are not available, the individual batch scripts can be simply run via `R CMD BATCH`, e.g.:
```R
R CMD BATCH --no-restore --no-save VSPsnap_ir_fr_log_one_011622_1.R
```
And so on.

### Variant Co-occurrence
Variant co-occurrence is computed on each gff3 file available in GISAID (each corresponding to a viral isolate) - each time a pair of variants is observed to occur together in a viral isolate, the corresponding CC is incremented by 1. The CC computation is iterated across all days of the Covid pandemic (720+ days for this project, considering the time interval March'20 - March'22).  
CC computation is implemented in Python and, similarly to VSPsnap, processes single days at the time to be easily batched/parallelized. Cumulative CC is then obtained by summing up CC matrices across the time interval.  

Example (compute CC matrix for a single day):
```python
python coOccurrenceMiner_singleDay2.py <day, in mm/dd/yy format)>
```


