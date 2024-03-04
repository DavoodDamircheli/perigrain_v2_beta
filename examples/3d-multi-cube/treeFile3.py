import ast

def list_classes_and_methods(tree):
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef):
            if hasattr(node, 'parent') and isinstance(node.parent, ast.Module):
                print(f'Top-level function: {node.name}')

        if isinstance(node, ast.ClassDef):
            print(f'Class: {node.name}')
            for sub_node in ast.walk(node):
                if isinstance(sub_node, ast.FunctionDef):
                    print(f'    Method: {sub_node.name}')

def add_parent(tree):
    for node in ast.walk(tree):
        for child in ast.iter_child_nodes(node):
            child.parent = node

if __name__ == "__main__":
    filename = 'setup.py'  # Replace with the actual file name
    with open(filename, 'r') as f:
        tree = ast.parse(f.read())

    add_parent(tree)  # This will add the 'parent' attribute
    list_classes_and_methods(tree)

