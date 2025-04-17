#!/bin/bash
#
#' @author Kalaivani Raju, Minkyung Sung, Julian Gough
#
# Download and prepare gene annotation data, if it does not exist already
# Checks for existence of gene annotation files for hg19, hg38, mm9 and mm10
# If they do not exist, downloads and prepares gene annotation files to [result_dir]/ext_annotation
#
# Args:
#       results_dir (str) is the results directory to store Mogrify outputs and intermediate result files
#                       e.g., ../../results
#                                   /projects/ci/C001_mogrify/results
# Returns:
#   None
# Notes:
#   This script is stored so we know exactly what version of gene annotation was used.
#   We have had problems in the past with changing gene annotations with different versions.
#   This script is primarily for the reproducibility.

results_dir='data/'
######################### Human hg19 ######################################################
if [ ! -d "${results_dir}" ]
then
  mkdir ${results_dir}
fi
file="${results_dir}/GRCh37_geneannot.txt"
if [ ! -f "$file" ]
then
  wget -O $file 'http://grch37.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query>
  <Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
    <Dataset name = "hsapiens_gene_ensembl" interface = "default" >
      <Attribute name = "ensembl_gene_id"	/>
      <Attribute name = "ensembl_transcript_id"	/>
      <Attribute name = "description"	/>
      <Attribute name = "external_gene_name"	/>
      <Attribute name = "external_gene_source"	/>
      <Attribute name = "gene_biotype"	/>
      <Attribute name = "transcript_biotype"	/>
      <Attribute name = "entrezgene_id"	/>
      <Attribute name = "hgnc_symbol"	/>
      <Attribute name = "hgnc_id"	/>
      <Attribute name = "ensembl_gene_id_version"	/>
      <Attribute name = "ensembl_transcript_id_version"	/>
      <Attribute name = "ensembl_peptide_id"	/>
    </Dataset>
  </Query>'
  sed -i -e '1iensembl_gene_id	ensembl_transcript_id	description	gene_name	external_gene_source	gene_biotype	transcript_biotype	entrezgene_id	hgnc_symbol	hgnc_id	ensembl_gene_id	version	ensembl_transcript_id_version	ensembl_peptide_id\' $file
fi

######################### Human hg38 ######################################################
file="${results_dir}/GRCh38_geneannot.txt"
if [ ! -f "$file" ]
then
  wget -O $file 'http://Sep2019.archive.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query>
  <Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
    <Dataset name = "hsapiens_gene_ensembl" interface = "default" >
      <Attribute name = "ensembl_gene_id"	/>
      <Attribute name = "ensembl_transcript_id"	/>
      <Attribute name = "description"	/>
      <Attribute name = "external_gene_name"	/>
      <Attribute name = "external_gene_source"	/>
      <Attribute name = "gene_biotype"	/>
      <Attribute name = "transcript_biotype"	/>
      <Attribute name = "entrezgene_id"	/>
      <Attribute name = "hgnc_symbol"	/>
      <Attribute name = "hgnc_id"	/>
      <Attribute name = "ensembl_gene_id_version"	/>
      <Attribute name = "ensembl_transcript_id_version"	/>
      <Attribute name = "ensembl_peptide_id"	/>
    </Dataset>
  </Query>'
  sed -i -e '1iensembl_gene_id	ensembl_transcript_id	description	gene_name	external_gene_source	gene_biotype	transcript_biotype	entrezgene_id	hgnc_symbol	hgnc_id	ensembl_gene_id_version	ensembl_transcript_id_version	ensembl_peptide_id\' $file
fi

############################### Mouse mm10 ####################################################
# Sep.2019.Ensembl98
file="${results_dir}/GRCm38_geneannot.txt"
if [ ! -f "$file" ]
then
  wget -O $file 'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?>
  <!DOCTYPE Query>
  <Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >

        <Dataset name = "mmusculus_gene_ensembl" interface = "default" >
                <Attribute name = "ensembl_gene_id"	/>
                <Attribute name = "ensembl_transcript_id"	/>
                <Attribute name = "description"	/>
                <Attribute name = "external_gene_name"	/>
                <Attribute name = "gene_biotype"	/>
                <Attribute name = "transcript_biotype"	/>
                <Attribute name = "entrezgene_trans_name"	/>
                <Attribute name = "mgi_symbol"	/>
                <Attribute name = "mgi_id"	/>
                <Attribute name = "ensembl_peptide_id"	/>
        </Dataset>
  </Query>'
  sed -i -e '1iensembl_gene_id	ensembl_transcript_id	description	gene_name	gene_biotype	transcript_biotype	entrezgene_trans_name	mgi_symbol	mgi_id	ensembl_peptide_id\' $file
fi
