![Coeus Logo](coeus/assets/coeus.png)

A Plotly dashboard to facilitate the visual exploration of gene neighborhood data. 

Coeus's primary use case is for the analysis of the genomic contexts of antimicrobial resistance (AMR) genes in bacterial genomes. It allows for easy visual comparison of gene order alongside gene neighborhood similarity-based clustering results showcasing differences between a gene's neighborhood across genomes. 

Coeus requires output files generated from [Gene-Order-Workflow](https://github.com/JTL-lab/Gene-Order-Workflow). Gene-Order-Workflow is a Nextflow based workflow to extract AMR gene neighborhoods from assembly and annotation files, perform unsupervised learning to assign neighborhood clusters, and generate ML classifier predictions for candidate divergent AMR genes (WIP). 

### Installation 
Coeus is set up with Poetry for dependencies management. 

To install Poetry: 
```
python3 -m pip install pipx
python3 -m pipx ensurepath
pipx install poetry
```

To then install all project dependencies, within the Coeus repository run: 
```
poetry install
```

If you're interested in viewing details on the package dependencies, within the Coeus repository you can run: 
```
poetry show --help 
```

### Using the dashboard
1. Ensure the `JSON` directory containing all gene order visualizations generated from `Gene-Order-Workflow` is present in `assets/clustermap/JSON` (sample data is included as a placeholder). 

2. If you do not intend to generate clustering figures dynamically using the available hyperparameter controls in the dashboard and want to use precomputed figures instead (generally recommended for large datasets), ensure the `clustering` directory containing all clustering ouputs generated from `Gene-Order-Workflow` is present in `assets/clustering` (sample data is included as a placeholder).

3.Once you've installed the required dependencies, you can use the dashboard on your local machine by running the following within the src directory (/Coeus/coeus): 
```
python coeus.py
```

### Acknowledgements 
To render gene order visualizations within the dashboard, this project uses code developed by Cameron Gilchrist (gamcil) for the D3 chart clustermap.js, which is used and modified here under the MIT license.

Citation: 
```
clinker & clustermap.js: Automatic generation of gene cluster comparison figures.
Gilchrist, C.L.M., Chooi, Y.-H., 2020.
Bioinformatics. doi: https://doi.org/10.1093/bioinformatics/btab007
```

