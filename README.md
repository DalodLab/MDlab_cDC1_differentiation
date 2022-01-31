# Harnessing single cell RNA sequencing to identify dendritic cell types, characterize their biological states and infer their activation trajectory

Article Information:

Authors: Ammar Sabir Cheema 1, Kaibo Duan 2, Marc Dalod 1, Thien-Phong Vu Manh 1

1 Aix-Marseille Univ, Centre National de la Recherche Scientifique (CNRS), Institut National de la Santé et de la Recherche Médicale (INSERM), Centre d'Immunologie de Marseille-Luminy (CIML), Marseille, France.

2 Singapore Immunology Network (SIgN), A*STAR, 8A Biomedical Grove, Singapore 138648, Singapore

**Summary:** 

Dendritic cells (DCs) orchestrate innate and adaptive immunity, by translating the sensing of distinct danger signals into the induction of different effector lymphocyte responses, to induce different defense mechanisms suited to face distinct types of threats. Hence, DCs are very plastic, which results from two key characteristics. First, DCs encompass distinct cell types specialized in different functions. Second, each DC type can undergo different activation states, fine-tuning its functions depending on its tissue microenvironment and the pathophysiological context, by adapting the output signals it delivers to the input signals it receives. Hence, to better understand DC biology and harness it in the clinic, we must determine which combinations of DC types and activation states mediate which functions, and how.
To decipher the nature, functions and regulation of DC types and their physiological activation states, one of the methods that can be harnessed most successfully is ex vivo single cell RNA sequencing (scRNAseq). However, for new users of this approach, determining which analytics strategy and computational tools to choose can be quite challenging, considering the rapid evolution and broad burgeoning of the field. In addition, awareness must be raised on the need for specific, robust and tractable strategies to annotate cells for cell type identity and activation states. It is also important to emphasize the necessity of examining whether similar cell activation trajectories are inferred by using different, complementary methods. In this chapter, we take these issues into account for providing a pipeline for scRNAseq analysis and illustrating it with a tutorial reanalyzing a public dataset of mononuclear phagocytes isolated from the lungs of naïve or tumor-bearing mice. We describe this pipeline step-by-step, including data quality controls, dimensionality reduction, cell clustering, cell cluster annotation, inference of the cell activation trajectories and investigation of the underpinning molecular regulation. It is accompanied with a more complete tutorial on Github. We anticipate that this method will be helpful for both wet lab and bioinformatics researchers interested in harnessing scRNAseq data for deciphering the biology of DCs or other cell types, and that it will contribute to establishing high standards in the field.


# Goal of the Github:

This github project contains the instructions and material to reproduce the analysis reported in the book chapter. Source code is available in the github repository. Required data and builded Docker images are available for download from [zenodo](https://zenodo.org/ "Google's Homepage"). Instructions to reproduce the analysis are provided below. 

To reproduce the analysis, you have to first, prepare the environments (see "Prepare the Environments" section below), then execute the analysis step by step (see "Run the analysis" section below).

# Datasets used in the analysis

In this analysis 2 pre-processed raw counts of scRNAseq data (one file for naïve and one file for tumor-bearing lungs) were used which can be downloaded from 

https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3832nnn/GSM3832735/suppl/GSM3832735_wt_naive_gex.csv.gz
https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3832nnn/GSM3832737/suppl/GSM3832737_wt_tumor_gex.csv.gz

and the 2 metadata files containing the Antibody-Derived Tags (ADT) information, coming from to the original publication [6].

https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3832nnn/GSM3832736/suppl/GSM3832736_wt_naive_adt.csv.gz
https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3832nnn/GSM3832738/suppl/GSM3832738_wt_tumor_adt.csv.gz


# Prepare the Environments

1. Download “Immgen phase 1” files  that are both necessary for the generation of the signature files used for the Connectivity Map (CMAP) analysis [12], one expression file containing the normalized gene expression data (.gct), one class file providing the cell type identity for each sample/microarray (.cls).

2. Download two signature files used as inputs for the CMAP analysis, (one for positive and one for negative signatures), in case one wants to skip the signature generation step.

3. Download the two .R scripts necessary to run single cell CMAP analyses, in order to assess the enrichment of transcriptomic signatures on single cells, for cell type identification.
    
4. Although not mandatory, we provide a Docker image in order to simplify the reproducibility of our analyses. 
