# isoforms_project
Looking for isoforms in single cell RNAseq data (smartSeq 3 data)

# repo structure
- **processOnDardel** folder contains the Snakefile with analyses pertaining to processing the raw fastq files per plate. The reads undergo special adapter trimming, UMI extraction, mapping, cell-barcode correction, de-duplication, de-multiplexing and counting to produce an exon by cell count matrix.

- **Dania** folder contains (test) and initial scripts which were written to process the data, and which were then polished and re-written in the Snakemake workflow. This folder can be ignored.

- **cellbarcodes.py** is the script which does cell barcode correction and (optional) de-multiplexing.

