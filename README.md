# mpox_surveillance

1) Create consensus sequences using the Viralrecon pipeline. See the "Mpox" directory for the parameters and output from viralrecon.
2) Use the [Nextstrain Mpox pipeline](https://github.com/nextstrain/mpox) for downloading reference sequences and create phylogenetic trees. I cloned the mpox repo separately and run the analysis from there. 

Followed the "phylogenetic" pipeline in the mpox repo:

- Downloaded the references sequences from NCBI and metadata file using these links: https://data.nextstrain.org/files/workflows/mpox/sequences.fasta.xz and https://data.nextstrain.org/files/workflows/mpox/metadata.tsv.gz
- Created a similar metadata file for our own sequences. Use the consensus sequences from the viralrecon pipeline.
- I also ran the norwegian sequences throught nextclade (reference clade IIb) and downloaded the results. Add some of this info to the metadata file.
- Entered submission date as the same as collection date. Otherwise got error in "recency" step...
- Added all norwegian strains to "include.txt".
- Used the default config.yaml file for the hmpxv1 build. 