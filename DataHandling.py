import os
import pandas as pd


class DataHandling:
    def __init__(self, data_dir='data/', data_analysis_dir='data_analysis/'):
        self.data_dir = data_dir
        self.data_analysis_dir = data_analysis_dir

    def read_data(self, disease_name: str):
        data_file = os.path.join(self.data_dir, f'{disease_name}.csv')
        df = pd.read_csv(data_file)
        return df

    def handle_missing_values(self, df: pd.DataFrame):
        # Handle empty 'Reported gene(s)
        df['Reported gene(s)'].replace(['null', ''], pd.NA, inplace=True)
        df['Reported gene(s)'] = df.apply(
            lambda row: row['Name(s)'] if pd.isnull(row['Reported gene(s)']) and row['Type'] == 'Gene' else (
                'Other' if pd.isnull(row['Reported gene(s)']) and row['Type'] != 'Gene' else row['Reported gene(s)']),
            axis=1
        )
        return df

    def split_gene_fields(self, df: pd.DataFrame):
        # Splitting gene fields by various newline characters and explode into separate rows
        df['Reported gene(s)'] = df['Reported gene(s)'].astype(str).str.split(r'[\r\n]+')
        df = df.explode('Reported gene(s)')
        df['Reported gene(s)'] = df['Reported gene(s)'].str.strip()
        df = df[df['Reported gene(s)'].str.lower().replace('null', '').str.strip().astype(bool)]
        return df

    def calculate_presence(self, df: pd.DataFrame, disease_name: str):
        gene_freq_df = df['Reported gene(s)'].value_counts().reset_index()
        gene_freq_df.columns = ['Gene', 'Frequency']

        df['p_desc'] = df['p_desc'].astype(str).str.strip()
        phenotype_freq_df = df['p_desc'].value_counts().reset_index()
        phenotype_freq_df.columns = ['Phenotype', 'Frequency']

        if not os.path.exists(self.data_analysis_dir):
            os.makedirs(self.data_analysis_dir)

        if not os.path.exists(os.path.join(self.data_analysis_dir, disease_name)):
            os.makedirs(os.path.join(self.data_analysis_dir, disease_name))

        gene_freq_file = os.path.join(self.data_analysis_dir, f'{disease_name}/{disease_name}_gene_frequency.csv')
        phenotype_freq_file = os.path.join(self.data_analysis_dir, f'{disease_name}/{disease_name}_phenotype_frequency.csv')

        gene_freq_df.to_csv(gene_freq_file, index=False)
        phenotype_freq_df.to_csv(phenotype_freq_file, index=False)

        return gene_freq_file, phenotype_freq_file

