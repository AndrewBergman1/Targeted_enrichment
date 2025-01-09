**Author: Andrew Bergman**
**Course: BINP50 30 CTS**

**Description** 
This thesis covers two master's degrees, one in molecular genetics and biotechnology and another in bioinformatics, both from Lunds University, Sweden. The project has been conducted at the Swedish Defence Research Institue (FOI) in Umeå, supervised by Andreas Sjödin. The aim of this project is to generate a TWIST panel for targeted enrichment, then  assess the efficacy of said TWIST panel.

**Generating in-solution hybridisation baits**: Generating TWIST baits. The entire workflow to generate TWIST baits is found in the Snakefile at ~/bait_gen/Snakefile. In short, the process is:

	a) Generate probes using Syotti (bait-gen tool) 
	b) Filter probes on thermodynamic bases (GC%, Gibbs free energy, Entropy)
	c) Filter probes on alignment to unwanted sequences (mitochondrial, chloroplast, human)

**Performing in-solution hybdridisation**: Targeted enrichment takes place in the lab and follows the protocol available in the appendix, refered to as: TWIST PROTOCOL - TARGET ENRICHMENT, authored by Olivia Wesula Luande, Aron  X and Y.

**Bioinformatic processing**: Bioinformatic assessment of targeted enrichment. The code for assembly, classification and alignments are be fond in the Snakefile at ~/analysis/*.snakefile. 

	a) Assembly of reads into contigs using Trinity.
	b) Classification of contigs using Kraken2.
	c) Alignment of reads to NT using BLAST. Code for plots are found in ~/analysis/*.R

 Plotting code is found in ~/analysis/*.R
