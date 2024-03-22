# mpox_surveillance

This is a repo for the first test of sequencing and analysing mpox samples at NIPH. This is still work under development.  

Steps taken:
1) Created consensus sequences using the Viralrecon pipeline. Used the following parameters:
```
   {
    "input": "samplesheet.csv",
    "platform": "illumina",
    "protocol": "metagenomic",
    "outdir": "viralrecon",
    "genome": "NC_063383.1",
    "fasta": "https://github.com/nf-core/test-datasets/raw/viralrecon/genome/NC_063383.1/GCF_014621545.1_ASM1462154v1_genomic.220824.fna.gz",
    "gff": "https://github.com/nf-core/test-datasets/raw/viralrecon/genome/NC_063383.1/GCF_014621545.1_ASM1462154v1_genomic.220824.gff.gz",
    "nextclade_dataset": "https://github.com/nf-core/test-datasets/raw/viralrecon/genome/NC_063383.1/nextclade_hMPXV_NC_063383.1_2022-08-19T12_00_00Z.tar.gz",
    "nextclade_dataset_name": "hMPXV",
    "nextclade_dataset_reference": "NC_063383.1",
    "nextclade_dataset_tag": "2022-08-19T12:00:00Z",
    "skip_pangolin": true,
    "skip_asciigenome": true,
    "skip_variants_quast": true,
    "skip_variants_long_table": true,
    "kraken2_variants_host_filter": true,
    "variant_caller": "bcftools",
    "skip_ivar_trim": true,
    "skip_assembly": true
}
```

2) Used the [Nextstrain Mpox pipeline](https://github.com/nextstrain/mpox) for downloading reference sequences and create a phylogenetic trees. The analyses were run from within a clone of the repo and inside the `phylogenetic` directory using these steps:
   
- Downloaded the references sequences from NCBI and metadata file using these links: https://data.nextstrain.org/files/workflows/mpox/sequences.fasta.xz and https://data.nextstrain.org/files/workflows/mpox/metadata.tsv.gz
- Created a metadata file for our own sequences with the same structure as the downloaded data. See the R script `create_metadata.R`.   
- I also ran the norwegian sequences throught nextclade (reference clade IIb) and downloaded the results as a tsv file. Some of the nextclade results were added to the metadata (see the R script).  
- I entered submission date as the same as collection date. Otherwise got error in a "recency" step... But the sequences are not submitted anywhere.
- Added all norwegian strains to the file `include.txt`.  
- Used the default `config.yaml` file for the hmpxv1 build. This should probably be tweaked in the subsampling.
