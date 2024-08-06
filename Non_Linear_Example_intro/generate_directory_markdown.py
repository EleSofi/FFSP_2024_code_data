import os

def generate_directory_structure(root_dir, output_file):
    with open(output_file, 'w') as f:
        f.write("# Project Structure\n\n")
        for root, dirs, files in os.walk(root_dir):
            level = root.replace(root_dir, '').count(os.sep)
            indent = ' ' * 4 * level
            f.write(f"{indent}- {os.path.basename(root)}/\n")
            sub_indent = ' ' * 4 * (level + 1)
            for file in files:
                f.write(f"{sub_indent}- {file}\n")

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
    output_file = os.path.join(current_script_dir, 'project_structure.md')
    generate_directory_structure(project_root, output_file)
    print(f"Project structure saved to {output_file}")
else:
    print("Project root not found. Project structure not generated.")
