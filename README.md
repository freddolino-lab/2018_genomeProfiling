# High-resolution mapping of transcriptional propensity in the *Escherichia coli* chromosome

## Description  

This project of high-resolution mapping of transcriptional propensity in *E. coli* genome is to profile the effect of position on chromosome on gene expression, using a multiplex strategy to randomly integrate standardized barcoded reporters in a *E. coli* population, combined with sequencing to pair transcriptional propensity, defined as amount of RNA produced per unit DNA, with positions of integration sites on genome.  

This code repository is a collection of source codes used in data analysis process, especially for high-throughput sequencing data. For biological interference of the project results, please see publication.  

## Installation

This projected presented as a repository is a collection of source codes, which means there is not a standardized installation process. Each source code file severs a certain analysis goal, can be run independently, and is documented in README files in the same directory. Certain dependencies are needed, e.g. GNU bash, python, Biopython, etc. See documentation for each source file for dependency details.  

## Usage  

The project, as a collection of source codes, can be divided into three parts, each representing a stage of data analysis and is collected in one of the sub-directories of the repository.  

For each part:
 + `raw_transcriptional_propensity_footprinting`: for high-throughput data analysis for mapping transcriptional propensity to genomic positions, supplemented data re-analysis, and statistics exploration on gene or annotation level;  
 + `spline_normalization_module_source_code`: for spline normalization, performed on raw RNA/DNA ratios as needed;  
 + `windower_snapshot_genome_profiling`: for windowing average of transcriptional propensity across genome.  

Detail usage, examples, input, and output documentations can be found in READMEs of each sub-directory, and the comments or docstrings of source code files.  

## Credits  

This project is credited to: Scott A. Scholz, Rucheng Diao, Michael B. Wolfe, Elayne M. Fivenson, Xiaoxia Nina Lin, and Peter L. Freddolino. See publication for details.  

The authors of the codes (in order of the listing of sub-directories):  
Rucheng Diao: diaorch@umich.edu;   
Peter L. Freddolino: petefred@umich.edu;  
Michael B. Wolfe: mbwolfe@umich.edu.  

Special thanks to Shweta Ramdas, for providing [the tool to calculate auto-correlation](https://github.com/shwetaramdas/autocorrelation).  
## Contacts

For questions regarding this repository, please contact Rucheng Diao (diaorch@umich.edu). For questions regarding the publication, please see the author listed ans contact the corresponding author(s).  

## License  

For license information, please see `raw_transcriptional_propensity_footprinting/LICENSE.md`, and in codes under `windower_snapshot_genome_profiling`, for licensing from different authors.
