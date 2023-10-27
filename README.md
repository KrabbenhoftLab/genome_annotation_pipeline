# Krabbenhoft Lab genome annotation pipeline using BRAKER and GeMoMa
## WORK IN PROGRESS

Authors: Dan MacGuigan*, Nate Backenstose, Christopher Osborne

[*dmacguig@buffalo.edu](dmacguig@buffalo.edu)

## Annotation pipeline flowchart
![flowchart](/annotation_flowchart.png)

## Dependencies

> [!NOTE] 
> This pipeline has been built work specifically for the Krabbenhoft Lab's servers and the University at Buffalo HPC cluster. 
> We are in the process of revising this pipeline to work on any Linux system. 
> We plan to distribute this pipeline with a Docker image containing all dependencies in the future.
> Please stay tuned for updates.

- [RepeatModeler](https://www.repeatmasker.org/RepeatModeler/)
- [RepeatMasker](https://www.repeatmasker.org/)
- [NCBI BLAST 2.4.0](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.4.0/)(for compatibility with ProtExcluder)
- [ProtExcluder](https://github.com/NBISweden/ProtExcluder)
- perl
- [bbmap](https://sourceforge.net/projects/bbmap/)
- [bedtools](https://bedtools.readthedocs.io/en/latest/)
- [samtools](https://www.htslib.org/)
- [HISAT2](https://github.com/DaehwanKimLab/hisat2)
- [NCBI Datasets command line tools](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/)
- [BRAKER3](https://github.com/Gaius-Augustus/BRAKER)
- [AGAT](https://agat.readthedocs.io/en/latest/agat_how_does_it_work.html)
- [GeMoMa](https://www.jstacs.de/index.php/GeMoMa)
- [GenomeTools](http://genometools.org/)
- [EVidenceModeler](https://github.com/EVidenceModeler/EVidenceModeler/wiki)
- [R](https://cran.r-project.org/)
	- ggplot2 package
	- cowplot package

## Usage

First create a new directory that will contain all of your annotation data and results. 
Then clone this repository to that directory. 
You will need to make `genome-annotation` executable. 
To do this, run `chmod 755 genome-annotation` in your annotation directory.
Then run `./genome-annotation -h` to make sure it worked.

Before running the pipeline, be sure to set all of the variables in the `config.txt` file.

To see help options, run `./genome-annotation -h`.

To perform a step of the pipeline, run `./genome-annotation -s 1 -c config.txt`. 
Pipeline steps should be performed sequentially, except for steps 5 and 6, which can run simultaneously.
