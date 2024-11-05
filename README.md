## Human VCP mutant microglia display immune and lysosomal phenotypes independently of GPNMB
Benjamin E. Clarke1,2†*, Oliver J. Ziff1,2,3†*, Giulia Tyzack1,2, Marija Petrić Howe1,2, Yiran Wang1,2, Pierre Klein1,2, Claudia A. Smith4, Cameron A. Hall4, Adel Helmy4, Michael Howell2, Gavin Kelly2, Rickie Patani1,2,3*
2024
1 Department of Neuromuscular Diseases, Queen Square Institute of Neurology, University College London, London WC1N 3BG, UK.
2 The Francis Crick Institute, 1 Midland Road, London NW1 1AT, UK.
3 National Hospital for Neurology and Neurosurgery, University College London NHS Foundation Trust, London WC1N 3BG, UK.  
4 Division of Neurosurgery and Wolfson Brain Imaging Centre, Department of Clinical Neurosciences, University of Cambridge, Cambridge, UK.
†Co-first authors. *Corresponding authors: Benjamin Clarke ben.clarke@crick.ac.uk, Oliver Ziff oliver.ziff@crick.ac.uk and Rickie Patani: rickie.patani@ucl.ac.uk 

This repository contains an R Markdown file that includes the scripts necessary to generate each figure in the manuscript. The Rmd file is organised into sections corresponding to the figures in the manuscript. Please refer to the instructions provided within the Rmd file to reproduce the figures. We have included a separate script called "microglia_deseq_objects.R" that generates the DESeq2 and DEP objects required for the analysis. Additionally, the script entitled "R_workspace.R" loads all the necessary packages and functions required for the figures. Please run these scripts before executing the R Markdown file to ensure all dependencies are properly set up.

#### [Manuscript link](https://molecularneurodegeneration.biomedcentral.com/)
DOI : 10.1186/s13024-024-00773-1
MOND-D-24-00068R2

## Data availability

RNA sequencing data have been deposited in the GEO under accession number [GSE***](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE***). 

Mass Spectrometry data have been deposited in the ProteomeXchange Consortium via the PRIDE partner repository under identifier PXD***.

Post-mortem RNA sequencing data is accessible at [GSE137810](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137810) and https://collaborators.nygenome.org/. 

Raw fastq files were processed with [nf-core/rnaseq](https://nf-co.re/rnaseq) v3.9 utilising alignment with STAR and read quantification with salmon. 

Differential transcript expression was performed using DESeq2. Differential protein expression was performed with MaxQuant followed by DEP, which is a wrapper around limma.
