import requests
import time
import os
import pandas as pd


class GeneGrouping:
    def __init__(self, groups_file_path):
        self.group_cache = {}
        self.load_grouped_genes(groups_file_path)

    def load_grouped_genes(self, groups_file_path):
        if os.path.exists(groups_file_path):
            df = pd.read_csv(groups_file_path)
            for index, row in df.iterrows():
                self.group_cache[row['Gene']] = row['Group']

    def get_gene_group(self, gene):
        if gene in self.group_cache:
            print(f"Gene {gene} is in cache")
            return self.group_cache[gene]

        def fetch_gene_data(url):
            time.sleep(0.11) # To avoid overloading the API (max 10 requests per second
            response = requests.get(url, headers={"Accept": "application/json"})
            if response.status_code == 200:
                return response.json()
            else:
                print(f"Error {response.status_code} with gene {gene}")
                return None

        def get_largest_group(gene_data):
            gene_groups = gene_data.get('gene_group', [])
            gene_group_ids = gene_data.get('gene_group_id', [])
            largest_group, max_members = None, 0

            for group_id in gene_group_ids:
                group_data = fetch_gene_data(f"https://rest.genenames.org/search/gene_group_id/{group_id}")
                if group_data and group_data['response']['numFound'] > max_members:
                    max_members = group_data['response']['numFound']
                    largest_group = gene_groups[gene_group_ids.index(group_id)]

            self.group_cache[gene] = largest_group
            return largest_group

        gene_data = fetch_gene_data(f"https://rest.genenames.org/fetch/symbol/{gene}")
        if not gene_data or not gene_data['response']['docs']:
            # If the initial gene symbol doesn't return data, check previous symbols
            prev_symbol_data = fetch_gene_data(f"https://rest.genenames.org/search/prev_symbol/{gene}")
            if prev_symbol_data and prev_symbol_data['response']['docs']:
                new_gene = prev_symbol_data['response']['docs'][0]['symbol']
                print(f"Gene {gene} is a previous symbol for {new_gene}")
                gene_data = fetch_gene_data(f"https://rest.genenames.org/fetch/symbol/{new_gene}")
            else:
                return "Unknown"

        if gene_data and gene_data['response']['docs'][0].get('gene_group', []):
            return get_largest_group(gene_data['response']['docs'][0])
        else:
            print(f"Gene {gene} has no group")
            return "No group"
