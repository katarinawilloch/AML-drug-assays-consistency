import pandas as pd 
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from hungarian_algorithm import algorithm


fimm_oslo_clinical_info = pd.read_csv('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/FIMM_oslo_clinical_info.csv')

fimm_oslo_clinical_info.loc[fimm_oslo_clinical_info['lab'] == 'Enserink','X'] = 'Patient ' + fimm_oslo_clinical_info['X'].astype(str)
fimm_oslo_clinical_info.loc[fimm_oslo_clinical_info['FAB_class'].isna(),'FAB_class'] = 'None'
fimm_oslo_clinical_info.loc[fimm_oslo_clinical_info['mutations'].isna(), 'mutations'] = 'None'
fimm_oslo_clinical_info.loc[fimm_oslo_clinical_info['karyotype_class'].isna(), 'karyotype_class'] = 'None'

fimm_oslo_clinical_info['mutations'] = fimm_oslo_clinical_info['mutations'].str.split(',')
fimm_oslo_clinical_info['FAB_class'] = fimm_oslo_clinical_info['FAB_class'].str.split(',')

#fimm_clinical_info1 = fimm_clinical_info.explode('mutations')
#fimm_clinical_info['m_value'] = 1

fimm_oslo_clinical_info.to_csv('/Users/katarinawilloch/Desktop/UiO/Project 1/Code/test.csv')

category_1 = fimm_oslo_clinical_info[fimm_oslo_clinical_info['lab'] == 'FIMM']
category_2 = fimm_oslo_clinical_info[fimm_oslo_clinical_info['lab'] == 'Enserink']
print(fimm_oslo_clinical_info.columns)

# Create an empty list to collect results
results_list = []

for i, row_1 in category_1.iterrows():
    genes_1 = set(row_1['mutations'])
    fab_class_1 = set(row_1['FAB_class'])
    for j, row_2 in category_2.iterrows():
        genes_2 = set(row_2['mutations'])
        fab_class_2 = set(row_2['FAB_class'])

        common_values = genes_1.intersection(genes_2)
        common_fab = fab_class_1.intersection(fab_class_2)
        
        common_values = common_values | common_fab
        if common_values: 
            print(common_values)
            num_common_genes = len(list(common_values))
            print(len(list(common_values)))
            start = row_1['X']
            end = row_2['X']

            results_list.append({
                'start': start, 
                'end': end, 
                'weight': num_common_genes,
                'common_values':common_values
            })

# Create a DataFrame from the results list
result_df = pd.DataFrame(results_list)
result_df.to_csv('/Users/katarinawilloch/Desktop/UiO/Project 1/Code/bipartite_df.csv')
print(result_df)

testing = {}
for _, row in result_df.iterrows():
    if row['start'] not in testing:
        testing[row['start']] = {}
    testing[row['start']][row['end']] = row['weight']
print(testing)

G = {
	'Ann': {'RB': 3, 'CAM': 2, 'GK': 1},
	'Ben': {'LW': 3, 'S': 2, 'CM': 1},
	'Cal': {'CAM': 3, 'RW': 2, 'SWP': 1},
	'Dan': {'S': 3, 'LW': 2, 'GK': 1},
	'Ela': {'GK': 3, 'LW': 2, 'F': 1},
	'Fae': {'CM': 3, 'GK': 2, 'CAM': 1},
	'Gio': {'GK': 3, 'CM': 2, 'S': 1},
	'Hol': {'CAM': 3, 'F': 2, 'SWP': 1},
	'Ian': {'S': 3, 'RW': 2, 'RB': 1},
	'Jon': {'F': 3, 'LW': 2, 'CB': 1},
	'Kay': {'GK': 3, 'RW': 2, 'LW': 1, 'LB': 0}
}
print(G)


# Create a bipartite graph
B = nx.Graph()

# Add nodes with the bipartite attribute
category_1_nodes = result_df['start'].unique()
category_2_nodes = result_df['end'].unique()

B.add_nodes_from(category_1_nodes, bipartite=0)  # Set category 1 nodes
B.add_nodes_from(category_2_nodes, bipartite=1)  # Set category 2 nodes

# Add edges with weights
for _, row in result_df.iterrows():
    B.add_edge(u_of_edge = row['start'],  v_of_edge = row['end'], weight=row['weight'])

# Draw the bipartite graph
# Define the position of nodes based on bipartite sets
pos = nx.bipartite_layout(B, nodes=category_1_nodes)

# Increase figure size
plt.figure(figsize=(100, 100))

# Draw nodes with increased size
nx.draw_networkx_nodes(B, pos, nodelist=category_1_nodes, node_color='lightblue', node_size=800, label='FIMM')
nx.draw_networkx_nodes(B, pos, nodelist=category_2_nodes, node_color='lightgreen', node_size=800, label='Enserink')

# Draw edges with increased width
edges = B.edges(data=True)
weights = [data['weight'] for _, _, data in edges]
nx.draw_networkx_edges(B, pos, edgelist=edges, width=[weight * 2 for weight in weights], alpha=0.6)

# Draw labels with increased font size
nx.draw_networkx_labels(B, pos, font_size=40, font_family='sans-serif')

# Add legend and title with larger font size
plt.legend(["FIMM", "Enserink"], fontsize=50)
plt.title("Bipartite Graph", fontsize=50)

# Display the plot
plt.show()



# Check if each node in category 1 has edges to all nodes in category 2
missing_edges = []

for start in category_1_nodes:
    start_edges = set(B.neighbors(start))
    required_edges = set(category_2_nodes)
    
    if not required_edges.issubset(start_edges):
        missing_for_start = required_edges - start_edges
        missing_edges.append((start, list(missing_for_start)))

# Print out missing edges
for start, missing in missing_edges:
    print(f"Patient {start} is missing edges to: {', '.join(missing)}")

# If no missing edges were found
if not missing_edges:
    print("Each start node has edges to every end node.")

print(nx.get_edge_attributes(B,'weight'))

# Compute the maximum bipartite matching
matching = nx.algorithms.bipartite.matching.maximum_matching(B, top_nodes=category_1_nodes)
#matching = algorithm.find_matching(testing, matching_type = 'max', return_type = 'list')

# Separate matches into two dictionaries for clear display
matches_category_1 = {k: v for k, v in matching.items() if k in category_1_nodes}
matches_category_2 = {v: k for k, v in matching.items() if v in category_2_nodes}

# Print the results
print("Matches:")
for start in matches_category_1:
    end = matches_category_1[start]
    print(f"{start} -> {end}")

results_list = []
common_id_counter = 1
for start, end in matches_category_1.items():
    results_list.append({
        'start': start,
        'end': end,
        'common_id': f"ID_{common_id_counter}"
    })
    common_id_counter += 1  # Increment the counter for each new ID

result_df = pd.DataFrame(results_list)
print(result_df)
result_df.to_csv('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/matching/FIMM_Oslo_matching.csv')

# Print the number of matches
print(f"\nTotal number of matches: {len(matches_category_1)}")
print(matching)