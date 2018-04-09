# Abstract

The JCVI Pan-Genome Pipeline is a collection of programs designed to allow users to run PanOCT along with numerous other tools that support and extend the capabilities of PanOCT.  PanOCT, or Pan-genome Ortholog Clustering Tool, was created as a tool for pan-genomic analysis of closely related prokaryotic species or strains.  The JCVI Pan-Genome Pipeline consists of a series of command-line utilities, invoked by a single wrapper, that: prepare input for analysis, invoke third-party analysis tools such as NCBI Blast+, run PanOCT, annotate features of the resulting pan-genome, and generate figures and tables to visualize the results of the analysis.  The pipeline is capable of running in a hierarchical mode, lowering the RAM used by the pipeline, by analyzing subdivisions of the input before running a final PanOCT run that incorporates the results in those subdivisions.

# Publications
  1. PanOCT: automated clustering of orthologs using conserved gene neighborhood for pan-genomic analysis of bacterial strains and closely related species.
   Fouts DE, Brinkac L, Beck E, Inman J, Sutton G.
   http://www.ncbi.nlm.nih.gov/pubmed/22904089
