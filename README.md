RNA silencing safeguards plant fertility during viral infection and
decreases Turnip rosette virus vertical transmission
================
Aimer Gutiérrez-Díaz, Sanjana Holla, Inês Moura and Anders Hafrén

# Comparative viral transcriptomic analysis

RNA-seq expression data for uninfected and infected Arabidopsis thaliana
were obtained from NCBI Bioprojects: TuMV PRJNA788379
\[[1](#ref-gyula2022ecotype)\], TuYV and CaMV PRJEB49403
\[[2](#ref-chesnais2022comparative)\], TyMV PRJNA1103879
\[[3](#ref-clavel2024metabolic)\], TCV PRJNA336058
\[[4](#ref-wu2016analyses)\], CMV PRJNA1124548
\[[5](#ref-liu2025mutually)\] and ArLV1 PRJNA863409
\[[6](#ref-jiang2024deciphering)\]. Read processing, alignment, and
gene-level quantification was addressed by mapping with HISAT2 v2.2.1
(Kim et al., 2019) to TAIR10 reference genome and quantifying with
featureCounts (Liao et al., 2014). Differential expression was computed
separately within each using DESeq2 DEGs were calculated using Deseq2
(Love et al., 2014).

# Viral RdRps structural prediction and phylogenetics

Protein sequences corresponding to CMV 2a, TRV 134K, TyMV 206K, TRoV
P2ab, ALV1 P1, TuMV NIb, and PLrV, TuYV, RYMV, AhPV1, TMV, YoMV, TCV,
P1AMV RdRps were curated prior to structure prediction. For viruses
where the replication protein is polyprotein-derived (e.g., TRoV P2ab),
sequences were processed to extract the annotated mature peptide
corresponding to the RdRP-containing product, final sequences are
available in this repository. Each curated protein was then folded using
AlphaFold2 (Jumper et al., 2021), while CaMV P5 PBD was the only RdRp
experimentally elucidated (PDB: 8R0S) (Prabaharan et al., 2024). To
build a structure-based phylogeny, an initial structural reconstruction
was performed using the predicted RdRP models, and a non-LTR
retrotransposon reverse transcriptase structure (PDB: 8GH6) was included
as an outgroup, in a similar way to Wolf et al. (2018). Finally, a
consensus topology was obtained by performing an agreement analysis
between trees generated from DALI (Holm, 2022) and Foldtree (Moi et al.,
2025) outputs, retaining and scoring clades supported by both
approaches.

# References

<div id="refs" class="references csl-bib-body">

<div id="ref-gyula2022ecotype" class="csl-entry">

1\. Gyula P, Tóth T, Gorcsa T, Nyikó T, Sós-Hegedűs A, Szittya G.
Ecotype-specific blockage of tasiARF production by two different RNA
viruses in arabidopsis. Plos one. 2022;17:e0275588.

</div>

<div id="ref-chesnais2022comparative" class="csl-entry">

2\. Chesnais Q, Golyaev V, Velt A, Rustenholz C, Brault V, Pooggin MM,
et al. Comparative plant transcriptome profiling of arabidopsis thaliana
col-0 and camelina sativa var. Celine infested with myzus persicae
aphids acquiring circulative and noncirculative viruses reveals
virus-and plant-specific alterations relevant to aphid feeding behavior
and transmission. Microbiology Spectrum. 2022;10:e00136–22.

</div>

<div id="ref-clavel2024metabolic" class="csl-entry">

3\. Clavel M, Bianchi A, Kobylinska R, Groh R, Ma J, Papareddy RK, et
al. Metabolic enzymes moonlight as selective autophagy receptors to
protect plants against viral-induced cellular damage. bioRxiv.
2024;2024–05.

</div>

<div id="ref-wu2016analyses" class="csl-entry">

4\. Wu C, Li X, Guo S, Wong S-M. Analyses of RNA-seq and sRNA-seq data
reveal a complex network of anti-viral defense in TCV-infected
arabidopsis thaliana. Scientific reports. 2016;6:36007.

</div>

<div id="ref-liu2025mutually" class="csl-entry">

5\. Liu J-H, Lin Y, Li Y-X, Lang Z, Zhang Z, Duan C-G. A mutually
antagonistic mechanism mediated by RNA m6A modification in plant-virus
interactions. Nature Communications. 2025;16:10378.

</div>

<div id="ref-jiang2024deciphering" class="csl-entry">

6\. Jiang Z, Verhoeven A, Li Y, Geertsma R, Sasidharan R, Zanten M van.
Deciphering acclimation to sublethal combined and sequential abiotic
stresses in arabidopsis thaliana. Plant Physiology. 2024;kiae581.

</div>

</div>
