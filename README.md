## Data Processing for CVD Associated Genes Analysis

This project utilizes Python for preprocessing a large dataset involving genes and their associations with cardiovascular diseases (CVDs). Relations to CVDs and genes were acquired from the Ensembl database. Only markings from the Human Phenotype Ontology (HPO) project were used.

1. **Gene Presence Across Diseases:**
   - A data table is created to record the presence of each gene across different CVDs, based on raw data files from Ensembl.
   - If a gene is present for a particular disease, it is assigned a score of 1; otherwise, it gets a score of 0.
   - The total number of diseases associated with each gene is calculated.

2. **Gene Grouping:**
   - Each gene is classified into a gene group using the classification provided by the HUGO Gene Nomenclature Committee (HGNC).
   - A new data table is created for gene groups, where each gene's presence in diseases is scored and summed up.


 ### Requirements
- Python 3.9, [link](https://www.python.org/downloads/)
- Pandas

### Installation
Run the following commands in your terminal:

```bash
# Use git to clone or download the repository as a zip file
git clone https://github.com/LauriSarap/CVDAssociatedGenesAnalysis.git
cd CVDAssociatedGenesAnalysis
pip install pandas requests
python main.py
```

