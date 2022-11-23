# GeneCoeus
A Plotly dashboard to facilitate visual exploration of gene neighborhood data. 

Genecoeus's primary use case is for the analysis of the genomic contexts of antimicrobial resistance (AMR) genes in bacterial genomes. It allows for easy comparison of gene order alongside neighborhood clustering results showcasing similarities and differences between a gene's neighborhood representations across genomes. 

For generating clustering outputs, see Gene-Order-Workflow for a ready-to-use Nextflow pipeline that extracts AMR gene neighborhoods from assembly and annotation files and performs unsupervised learning to assign neighborhood clusters. 

### Acknowledgements 
This project uses code developed by Cameron Gilchrist (gamcil) for clustermap.js, which is used and modified here under the MIT license.

Citation: 
```
clinker & clustermap.js: Automatic generation of gene cluster comparison figures.
Gilchrist, C.L.M., Chooi, Y.-H., 2020.
Bioinformatics. doi: https://doi.org/10.1093/bioinformatics/btab007
```

