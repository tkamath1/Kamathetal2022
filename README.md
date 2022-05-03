# Kamathetal2022
## Code used to generate analyses and reproduce main figures for Kamath et al., 2022 Nature Neuro

### Single-cell genomic profiling of human dopamine neurons identifies a population that selectively degenerates in Parkinson’s disease (_Nature Neuroscience_, 2022)
#### Tushar V. Kamath<sup>1,2</sup>, Abdulraouf Abdulraouf, SJ Burris, Jonah Langlieb<sup>1</sup>, Vahid Gazestani<sup>1</sup>, Naeem Nadaf, Karol Balderrama, Charles Vanderburg, Evan Z. Macosko

<sup>1 -  Performed analysis</sup>

<sup>2 - Analysis and project lead (contact: tkamath@broadinstitute.org)</sup> 

## Code

Analyses/[^1] - Grouped by primary analyses performed for this paper

Auxiliary_scripts/ - A couple of extra scripts for helper functions

Main_Figs/ - code to generate main figure panels

ExtendedData_Figs/ - code to generate extended data/supp figure panels

## Data
All raw, unaligned sequence-level data from human samples are available via dbGAP (requires formal request for access to ensure federal privacy concerns compliance): phs002879.v1.p1[^2]



All digital gene expression matrices generated from this study have been made publicly available via GEO at: [GSE178265](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178265)

In addition, all processed gene expression matrices, UMAP coordinates, sample metadata, and cell-type annotations have been made publicly available via Broad's Single Cell Portal:

SCP link 1 (single-nuc RNA-seq data): https://singlecell.broadinstitute.org/single_cell/study/SCP1768/

SCP link 2 (Slide-seq data): https://singlecell.broadinstitute.org/single_cell/study/SCP1769/

## FAQs  [^3]
Q: Where is the metadata?

A: The easiest place to access the metadata associated with each postmortem human donor and animal is in the metadata .tsv file hosted on Broad SCP. GEO also has metadata associated with each sample in the individual sample landing pages (all listed within the main GEO link, but you will need to navigate/scrape across each individual sample page to gather them together).

Q: Why are there two macaque datasets?

A: We generated Nurr-enriched DA neuron profiles from macaque (M. mulatta) for our preprint. Since our preprint release, we received higher-quality freshly perfused and well-preserved macaque (M. fascicularis) midbrain specimens. These specimens were used to generate snRNA-seq profiles of DA neurons in the final version of the paper. The M. fascicularis snRNA-seq data was used for species integration, defining subtypes in the macaque, and Slide-seq in the final version of the paper.

## Acknowledgements
This work was supported by the following funding sources: National Institutes of Health F30AG069446-01 (TK), National Institutes of Health DP2AG058488 (E.Z.M.), National Institutes of Health U01MH124602 (E.Z.M.), Chan Zuckerberg Initiative 2017-175259 (E.Z.M.).

In addition, we would like to thank Djordje Gveric (Multiple Sclerosis and Parkinson's Tissue Bank), Randall Woltjer (Oregon Health and Sciences University), and Rashed Nagra (Human Brain and Spinal Fluid Resource Center) for their contributions of postmortem tissue. We thank the NIH NeuroBioBank for facilitating the acquisition of postmortem brain tissue samples. We also thank Nicole Shultz and David Fitzpatrick (Max Planck Florida Institute for Neuroscience) for their donation of tree shrew brains. 

[^1]: There are many intermediate files that are used within these scripts. I've made an attempt at linking how to create each one in the scripts but please contact if you are having a hard time figuring out where each intermediate was generated from.
[^2]: NB: will update this with the actual link when we have the final landing page as there is an extended process we are actively working on to ensure full compliance with federal laws associated with sequence-level human data deposition
[^3]: This section will be periodically updated based on the correspondences received.

