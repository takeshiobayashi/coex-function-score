# About

This script caculates a partial AUC (false positive rate from 0% to 1%) for a givin coexpression matrix and pathway data. This score is used as a function score in [ATTED-II](https://atted.jp).

# Usage

./function_score.pl -d Ath-r.v21-01.G18957-S14741.combat_pca_subagging.ls.d -k dummy_pathway_data.txt -g dummy_paralog_data.txt

## Options
 -d: Directory of coexpression data provided in the bulk download page in ATTED-II.
     Example of data source
     https://zenodo.org/record/4961962/files/Ath-r.v21-01.G18957-S14741.combat_pca_subagging.ls.d.zip

 -k: Annotation file (tab-separated: Pathway ID, Entrez Gene IDs)
     Example of data source (KEGG FTP license is required):
     ftp://ftp.bioinformatics.jp/kegg/genes/organisms/ath/ath_link.tar.gz

 -g: (option) Paralog genes (tab-separated: Paralog ID, Entrez Gene IDs)
     Example of data source (KEGG FTP license is required):
     ftp://ftp.kegg.net/kegg/genes/links/genes_ncbi-proteinid.list.gz
     ftp://ftp.kegg.net/kegg/genes/links/genes_ko.list.gz

-M: use a smaller-is-better coexpression index. (Default: a larger-is-better index)

(c) 2022 Takeshi Obayashi