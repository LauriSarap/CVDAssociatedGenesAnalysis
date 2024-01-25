import os
import pandas as pd
import time
from DataHandling import DataHandling
from GeneGrouping import GeneGrouping

# Settings
DO_GENE_PRESENCE_ACROSS_DISEASES = False
DO_GENE_VARIANTS_ACROSS_DISEASES = False
DO_GENE_GROUPING = False
DO_GENE_GROUP_PRESENCE = True

dh = DataHandling('data/', 'data_analysis/')

if DO_GENE_PRESENCE_ACROSS_DISEASES:

    gene_presence = pd.DataFrame()
    for filename in os.listdir(dh.data_dir):
        if filename.endswith(".csv"):
            disease_name = filename[:-4]

            # Read and process data
            df = dh.read_data(disease_name)
            df = dh.handle_missing_values(df)
            df = dh.split_gene_fields(df)
            gene_freq_file, phenotype_freq_file = dh.calculate_presence(df, disease_name)

            gene_freq_df = pd.read_csv(gene_freq_file)

            # Add a column for the disease presence
            gene_freq_df[disease_name] = 1

            # Initialize gene_presence DataFrame if it's the first file
            if gene_presence.empty:
                gene_presence = gene_freq_df[['Gene', disease_name]].copy()
            else:
                # Merge the new data with the existing gene_presence DataFrame
                gene_presence = pd.merge(gene_presence, gene_freq_df[['Gene', disease_name]],
                                         on='Gene', how='outer')

            print(f'Processed {disease_name}.')


    ### GENE PRESENCE TABLE ###
    # Replace NaN with 0 and add a column to count the number of diseases per gene
    gene_presence.fillna(0, inplace=True)
    gene_presence['Number of Diseases'] = gene_presence.drop('Gene', axis=1).sum(axis=1)

    # Sort the columns and rows
    cols = ['Gene', 'Number of Diseases'] + [col for col in gene_presence.columns if col not in ['Gene', 'Number of Diseases']]
    gene_presence = gene_presence[cols]
    gene_presence.sort_values(by='Number of Diseases', ascending=False, inplace=True)

    # Save the gene presence table
    gene_presence_file = os.path.join(dh.data_analysis_dir, 'gene_presence_across_diseases.csv')
    gene_presence.to_csv(gene_presence_file, index=False)

    print(f"Gene presence table created.")

if DO_GENE_VARIANTS_ACROSS_DISEASES:

    gene_variants = pd.DataFrame()
    for filename in os.listdir(dh.data_dir):
        if filename.endswith(".csv"):
            disease_name = filename[:-4]

            # Read and process data
            df = dh.read_data(disease_name)
            df = dh.handle_missing_values(df)
            df = dh.split_gene_fields(df)
            gene_variant_count = df['Reported gene(s)'].value_counts().reset_index()
            gene_variant_count.columns = ['Gene', disease_name]

            # Initialize gene_variants DataFrame if it's the first file
            if gene_variants.empty:
                gene_variants = gene_variant_count.copy()
            else:
                # Merge the new data with the existing gene_variants DataFrame
                gene_variants = pd.merge(gene_variants, gene_variant_count,
                                         on='Gene', how='outer')

            print(f'Processed {disease_name}.')

    gene_variants.fillna(0, inplace=True)
    gene_variants['Total Variants Across Diseases'] = gene_variants.drop('Gene', axis=1).sum(axis=1)
    cols = ['Gene', 'Total Variants'] + [col for col in gene_variants.columns if col not in ['Gene', 'Total Variants']]
    gene_variants = gene_variants[cols]
    gene_variants.sort_values(by='Total Variants', ascending=False, inplace=True)

    # Save the gene variants table
    gene_variants_file = os.path.join(dh.data_analysis_dir, 'gene_variants_across_diseases.csv')
    gene_variants.to_csv(gene_variants_file, index=False)

if DO_GENE_GROUPING:

    # Read the gene presence table
    gene_presence_file = os.path.join(dh.data_analysis_dir, 'gene_presence_across_diseases.csv')
    df_genes = pd.read_csv(gene_presence_file)

    # Cache the grouped genes
    gene_groups_table = os.path.join(dh.data_analysis_dir, 'gene_groups.csv')

    gp = GeneGrouping(groups_file_path=gene_groups_table)

    if os.path.exists(gene_groups_table):
        gene_groups = pd.read_csv(gene_groups_table)
    else:
        gene_groups = pd.DataFrame(columns=['Gene', 'Group'])

    gene_groups_presence = []

    if len(gene_groups_table) == 940:
        print("Grouping genes...")
        processed_genes = set(gene_groups['Gene'])
        processed_count = 0
        start_time = time.time()
        for gene in df_genes['Gene']:
            if gene != 'Other':
                group = gp.get_gene_group(gene)
            else:
                group = 'Other'
            gene_groups_presence.append(group)
            processed_count += 1
            gene_groups = gene_groups.append({'Gene': gene, 'Group': group}, ignore_index=True)
            gene_groups.to_csv(gene_groups_table, index=False)
            processed_genes.add(gene)

            if processed_count % 10 == 0:
                print(f'Processed {processed_count} genes in {time.time() - start_time} seconds.')
                start_time = time.time()
                print(f'Remaining: {len(df_genes) - processed_count} genes.')
    else:
        print("All genes are already grouped!")

if DO_GENE_GROUP_PRESENCE:

    gene_presence_file = os.path.join(dh.data_analysis_dir, 'gene_presence_across_diseases.csv')
    gene_groups_file = os.path.join(dh.data_analysis_dir, 'gene_groups.csv')
    df_gene_presence = pd.read_csv(gene_presence_file)
    df_gene_groups = pd.read_csv(gene_groups_file)

    gene_to_group = dict(zip(df_gene_groups['Gene'], df_gene_groups['Group']))

    group_presence_data = {}

    for index, row in df_gene_presence.iterrows():
        gene = row['Gene']
        group = gene_to_group.get(gene, 'No group')
        if group not in group_presence_data:
            group_presence_data[group] = row[2:]  # Initialize with the first occurrence
        else:
            group_presence_data[group] += row[2:]  # Accumulate the counts

    # Convert the accumulated data into a DataFrame
    df_group_presence = pd.DataFrame.from_dict(group_presence_data, orient='index')
    df_group_presence.reset_index(inplace=True)
    df_group_presence.columns = ['Gene Group'] + list(df_gene_presence.columns[2:])

    # Total columns
    df_group_presence['Number of Diseases'] = (df_group_presence.iloc[:, 1:] > 0).sum(axis=1)
    df_group_presence['Total Presence Strength'] = df_group_presence.iloc[:, 1:-1].sum(axis=1)
    cols = ['Gene Group', 'Number of Diseases', 'Total Presence Strength'] + list(df_group_presence.columns[1:-2])
    df_group_presence = df_group_presence[cols]

    # Save the new DataFrame
    output_file = 'gene_group_presence_across_diseases.csv'
    df_group_presence.to_csv(os.path.join(dh.data_analysis_dir, output_file), index=False)
