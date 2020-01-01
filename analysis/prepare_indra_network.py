import os
import pathlib
import pickle

import networkx
import numpy as np
import pandas as pd
import synapseclient

cwd = pathlib.Path(os.getcwd()) / "clustering" / "indra"
pathlib.Path.mkdir(cwd)

syn = synapseclient.Synapse()
syn.login()

syn_indra_graph = syn.get("syn21479217")
with open(syn_indra_graph.path, "rb") as f:
    indra_graph = pickle.load(f)

syn_function_clustering = syn.get("syn21478380")
function_clustering = pd.read_csv(syn_function_clustering.path)


# Modify edge attributes so that statements are removed
for edge_name in indra_graph.edges:
    edge = indra_graph.edges[edge_name]
    statement = edge["statements"][0]
    edge["stmt_type"] = statement["stmt_type"]
    edge["evidence_count"] = statement["evidence_count"]
    edge["curated"] = statement["curated"]
    del edge["statements"]
    for k, v in edge.items():
        if type(v) == np.float128:
            edge[k] = float(v)

# Add info about clusters to nodes
for row in function_clustering.itertuples():
    try:
        node = indra_graph.nodes[row.gene_name]
    except KeyError:
        print(row.gene_name, "not found")
        continue
    node["class_combined"] = row.class_combined
    node["class"] = row._3
    node["direction"] = row.direction
    for k, v in node.items():
        if type(v) == np.float128:
            node[k] = float(v)

# Remove node metadata dict
del indra_graph.graph["node_by_ns_id"]

with (cwd / "indra_graph_cleaned.pickle").open("wb") as f:
    pickle.dump(indra_graph, f)

networkx.write_graphml(indra_graph, cwd / "indra_graph_cleaned.graphml")


syn_indra = syn.get("syn21435525")

syn.store(
    synapseclient.File(str(cwd / "indra_graph_cleaned.graphml"), parent=syn_indra),
    activityName="Preprocess python INDRA network pickle file for R",
    used=["syn21479217", "syn21478380"],
    executed="https://github.com/clemenshug/erk_senescence/blob/master/analysis/prepare_indra_network.py"
)
