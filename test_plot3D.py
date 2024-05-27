import pytest
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from io import StringIO
import sys

def plot_molecule_3D(smiles):
    """
    Plot a molecule in 3D with different colors for different types of atoms and bonds between atoms.
    Args:
    - smiles (str): SMILES representation of the molecule.
    """
    # Generate RDKit molecule from SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print("Error: Unable to generate molecule from SMILES.")
        return
    
    # Add explicit hydrogens
    mol = Chem.AddHs(mol)

    # Generate 3D coordinates
    AllChem.EmbedMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol)

    # Extract atom coordinates and symbols
    coords = mol.GetConformer().GetPositions()
    symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]

    # Create 3D plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Define colors for different atom types
    atom_colors = {'H': 'gray', 'C': 'black', 'O': 'red', 'N': 'blue', 'F': 'green', 'Br': 'purple', 'Cl': 'yellow'}

    # Plot atoms with different colors
    for atom_idx, symbol in enumerate(symbols):
        x, y, z = coords[atom_idx]
        color = atom_colors.get(symbol, 'brown')  # Default to brown for unknown atoms
        ax.scatter(x, y, z, c=color, label=symbol, s=100)

    # Helper function to plot bonds
    def plot_bond(ax, start_pos, end_pos, offset=np.array([0, 0, 0])):
        ax.plot([start_pos[0] + offset[0], end_pos[0] + offset[0]], 
                [start_pos[1] + offset[1], end_pos[1] + offset[1]], 
                [start_pos[2] + offset[2], end_pos[2] + offset[2]], c='gray')

    # Plot bonds between atoms with different styles for single and double bonds
    for bond in mol.GetBonds():
        start_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        start_pos = coords[start_idx]
        end_pos = coords[end_idx]
        bond_type = bond.GetBondType()
        
        if bond_type == Chem.rdchem.BondType.SINGLE:
            plot_bond(ax, start_pos, end_pos)
        elif bond_type == Chem.rdchem.BondType.DOUBLE:
            # Calculate an offset direction perpendicular to the bond
            bond_vector = end_pos - start_pos
            normal_vector = np.cross(bond_vector, np.array([1, 0, 0]))
            if np.linalg.norm(normal_vector) < 1e-3:
                normal_vector = np.cross(bond_vector, np.array([0, 1, 0]))
            normal_vector /= np.linalg.norm(normal_vector)  # Normalize the vector
            offset = normal_vector * 0.125
            plot_bond(ax, start_pos, end_pos, offset)
            plot_bond(ax, start_pos, end_pos, -offset)

    # Set labels
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('Molecule in 3D')

    # Show legend
    handles = [plt.Line2D([0], [0], marker='o', color='w', markersize=10, markerfacecolor=color, label=symbol)
               for symbol, color in atom_colors.items()]
    ax.legend(handles=handles, title='Atom types', loc='best')

    # Show plot
    plt.show()


valid_smiles = "CCO"  
invalid_smiles = "invalid_smiles_string"

def test_valid_smiles():
    mol = Chem.MolFromSmiles(valid_smiles)
    assert mol is not None, "Failed to generate molecule from valid SMILES."

def test_add_hydrogens():
    mol = Chem.MolFromSmiles(valid_smiles)
    mol = Chem.AddHs(mol)
    num_hydrogens = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'H')
    assert num_hydrogens > 0, "No hydrogens were added to the molecule."

def test_generate_3D_coordinates():
    mol = Chem.MolFromSmiles(valid_smiles)
    mol = Chem.AddHs(mol)
    result = AllChem.EmbedMolecule(mol)
    assert result == 0, "3D coordinate generation failed."
    coords = mol.GetConformer().GetPositions()
    assert coords is not None and len(coords) > 0, "No 3D coordinates were generated."

def test_plot_molecule_3D_valid_smiles():
    try:
        plot_molecule_3D(valid_smiles)
    except Exception as e:
        pytest.fail(f"plot_molecule_3D raised an exception for valid SMILES: {e}")

def test_plot_molecule_3D_invalid_smiles():
    captured_output = StringIO()
    sys.stdout = captured_output
    plot_molecule_3D(invalid_smiles)
    sys.stdout = sys.__stdout__
    assert "Error: Unable to generate molecule from SMILES." in captured_output.getvalue()


plt.show = lambda: None