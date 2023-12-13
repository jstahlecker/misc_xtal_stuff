import argparse


def fetch_and_init(PDB):
    """
    Fetches the PDB structure and initializes the PyMOL visualization.
    
    Parameters:
    - PDB (str): The PDB code of the structure to fetch and visualize.
    """
    cmd.fetch(PDB)
    cmd.hide("everything")
    cmd.color("white")

def find_interactions_and_color(PDB, cutoff):
    """
    Finds interactions between chains and colors them based on a cutoff distance.
    
    Parameters:
    - PDB (str): The PDB code of the structure.
    - cutoff (float): The distance cutoff for interactions.
    """
    contact_distance = cutoff  
    chains = cmd.get_chains(PDB)

    # Create a set to store unique chain pairs
    unique_pairs = set()

    # Iterate over all possible pairs of chains
    for chain1 in chains:
        for chain2 in chains:
            if chain1 != chain2:
                # Sort chains for consistent unique pairs
                sorted_pair = tuple(sorted([chain1, chain2]))
                unique_pairs.add(sorted_pair)

    # Calculate and color interactions for unique pairs
    for chain1, chain2 in unique_pairs:
        # Select residues from each chain
        selection1 = f"{PDB} and chain {chain1}"
        selection2 = f"{PDB} and chain {chain2}"

        # Find residues within the updated distance threshold and color them red
        interaction_name = f"contacts_{chain1}_{chain2}"
        cmd.select(
            interaction_name,
            f"({selection1} within {contact_distance} of {selection2}) or ({selection2} within {contact_distance} of {selection1})"
        )
        cmd.color("red", f"byres {interaction_name}")


def symmetry_mates(PDB, cutoff):
    """
    Calculates symmetry contacts and colors them.
    
    Parameters:
    - PDB (str): The PDB code of the structure.
    - cutoff (float): The distance cutoff for symmetry contacts.
    """

    contact_distance = cutoff
    cmd.symexp("sym", f"{PDB}", selection=f"{PDB}", cutoff=contact_distance)

    cmd.select("sym_contacts", f"{PDB} within {contact_distance} of sym*")
    cmd.color("blue", "sym_contacts")

    cmd.hide("everything", "sym*")
    cmd.show("surface", f"{PDB}")
    cmd.deselect()




def main(pdb_code, cutoff):
    """
    Main function to fetch the structure, find interactions, and calculate symmetry contacts.
    
    Parameters:
    - pdb_code (str): The PDB code of the structure.
    - cutoff (float): The distance cutoff for interactions and symmetry contacts.
    """
    fetch_and_init(pdb_code)
    find_interactions_and_color(pdb_code, cutoff)
    symmetry_mates(pdb_code, cutoff)




if __name__ == "pymol":
    parser = argparse.ArgumentParser(description="Find and color interactions with symmetry contacts.")

    parser.add_argument("--pdb", type=str, required=True, help="PDB code for the structure")
    parser.add_argument("--cutoff", type=float, default=4.0, help="Cutoff distance for interactions and symmetry contacts")   

    args = parser.parse_args()

    main(args.pdb, args.cutoff)