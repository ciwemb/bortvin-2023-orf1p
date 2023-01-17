### This repository contains the pipelines used for sequencing analyses published in De Luca C, Gupta A, and Bortvin A. (2023) Ribonucleoprotein condensation driven by retrotransposon LINE-1 sustains RNA integrity and translation in mouse spermatocytes. bioRxiv https://doi.org/10.1101/2023.01.09.523313

![Figure 1](https://github.com/ciwemb/bortvin-2023-orf1p/blob/main/banner.png)
Colocalization of L1 RNA (green) and ORF1p (magenta) in the cytoplasm of *Mael*-mutant mouse spermatocytes.
---
### Characterization of the RNA present in L1 ORF1p macromolecular complexes isolated from Maelstrom-mutant mouse testes

piRNA-deficient Maelstrom (*Mael*) null mice are characterized by a strong upregulation of retrotransposon LINE-1 (L1) in meiotic spermatocytes. This defect turns out in the accumulation of L1 RNA and ORF1p in their cytoplasm and the formation of prominent ribonucleoprotein aggregates. 
We used 3-months-old *Mael-/-* male mice to characterize the RNA present in those ORF1p aggregates. To favor the isolation of complexed versus free ORF1p protein, we first fractionated *Mael-/-* testis extracts (that we refer to as TOTAL) by sucrose gradient ultracentrifugation, in the presence of EDTA. We then pooled the sucrose fractions where ORF1p macromolecular complexes sediment (fractions 5-8) and used this pool as the INPUT for an anti-ORF1p co-immunoprecipitation (IP) followed by RNA-seq.

The folder "/RNA-seq" contains the pipeline used for the analysis of the generated RNA-seq data. Raw data have been deposited in NCBI’s Gene Expression Omnibus and are accessible through GEO Series accession number GSE222416 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE222416).

Code contributors: Anuj Gupta, Chiara De Luca

---
### Effect of the upregulation of retrotransposon LINE-1 on translation in the mouse testis

We used *Mael* null mice to investigate the effect of LINE-1 retrotransposon upregulation on translation efficiency (TE) in the germline. We collected *Mael+/-* and *Mael-/-* testes at postnatal day 16 (P16), when LINE-1 upregulation in *Mael-/-* spermatocytes is manifest, and performed ribosome profiling. RNA-seq analysis of the same samples was performed in parallel to evaluate the contribution of transcriptional changes to variations in footprint levels.

The folder "/ribosome_profiling" contains the pipeline used for the ribosome profiling data analysis. Raw data have been deposited in NCBI’s Gene Expression Omnibus and are accessible through GEO Series accession number GSE222415 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE222415).

Code contributor: Chiara De Luca
