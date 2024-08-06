import os
import matplotlib.pyplot as plt
import networkx as nx

def build_tree_graph(root_dir):
    G = nx.DiGraph()
    for root, dirs, files in os.walk(root_dir):
        for d in dirs:
            parent = os.path.relpath(root, root_dir)
            child = os.path.join(parent, d)
            G.add_edge(parent, child)
        for f in files:
            parent = os.path.relpath(root, root_dir)
            child = os.path.join(parent, f)
            G.add_edge(parent, child)
    return G

def plot_tree_graph(G, output_file):
    pos = nx.multipartite_layout(G, subset_key=lambda n: n.count(os.sep))
    plt.figure(figsize=(20, 10))
    nx.draw(G, pos, with_labels=True, node_size=50, font_size=8, arrows=False)
    plt.savefig(output_file, format='PNG')
    plt.show()

# Automatically determine the path to the current script's directory
current_script_dir = os.path.dirname(os.path.abspath(__file__))

# Find the project root by looking for the root marker directory
def find_project_root(starting_dir, root_marker='FFSP_code_revised_2024'):
    current_dir = os.path.abspath(starting_dir)
    while True:
        if root_marker in os.listdir(current_dir):
            return os.path.join(current_dir, root_marker)
        parent_dir = os.path.abspath(os.path.join(current_dir, os.pardir))
        if parent_dir == current_dir:
            raise FileNotFoundError(f"Project root directory '{root_marker}' not found.")
        current_dir = parent_dir

# Determine the project root directory
try:
    project_root = find_project_root(current_script_dir)
    print(f"Project root found: {project_root}")
except FileNotFoundError as e:
    print(e)
    project_root = None

if project_root:
    # Build and plot the tree graph
    tree_graph = build_tree_graph(project_root)
    output_file = os.path.join(current_script_dir, 'project_structure.png')
    plot_tree_graph(tree_graph, output_file)
    print(f"Tree plot saved to {output_file}")
else:
    print("Project root not found. Tree plot not generated.")
