
The script is used to produce the needed DNA fragments for merging into the full-length gene with PCR. All of the resultant DNA fragments have similar melting temperatures in the overlap regions, which is important for the successful assembly of target gene via PCR.

Here we use the example file ‘beta_original.fasta’, which is the original gene sequence before sequence optimization, to show the running process of our script. We show the steps in the macOS. The steps in the Windows are largely the same, except that the Terminal prompt is different.

(1) Download the test gene sequence file 'beta_original.fasta' from the S2 directory of this project.
(2) Open the 'beta_original.fasta' with TextEdit to get the gene sequence. Perform gene sequence optimization with the Codon Optimization Tool (https://sg.idtdna.com/pages/tools/codon-optimization-tool) from Integrated DNA Technologies, Inc. (IA, USA), and obtain the optimized gene sequence ‘beta_optimized.fasta’. This step can not only improve gene expression but also simplify the DNA structure to help PCR synthesis of the target gene. There are other commercial companies also provide sequence optimization services, for example, the GenSmart Codon Optimization (https://www.genscript.com/gensmart-free-gene-codon-optimization.html) of GenScript (NJ, USA).
(3) Open Terminal and make a new directory (for example, 'gene_fragments') on the Desktop and copy ‘gene_synthesis.py’ (in directory S1), and ‘beta_optimized.fasta’ (in directory S3) to the directory. The copy can be done with mouse drag and drop, or with terminal command cp (in Windows, the command is copy). Go the directory 'gene_fragments' in the terminal.
```zsh
% cd
% cd Desktop
% mkdir gene_fragments
% cd gene_fragments
```
(4) Run the script with the following command.
```zsh
python3.10 gene_synthesis.py beta_optimized.fasta
```
(5) A figure showing the melting temperature of each DNA fragment will show up, and the DNA fragments are stored in a folder in the directory ‘gene_fragments’.
