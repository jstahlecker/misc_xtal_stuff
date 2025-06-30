import json
import requests


def get_list_of_pdbs(uniprot, result_type="polymer_instance"):
    """Get list of PDBs based on uniprot id"""
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    query = {
            "query": {
                "type": "group",
                "logical_operator": "and",
                "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                    "operator": "exact_match",
                    "value": uniprot,
                    "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession"
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                    "operator": "exact_match",
                    "value": "UniProt",
                    "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_name"
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                    "operator": "exact_match",
                    "value": "X-RAY DIFFRACTION",
                    "attribute": "exptl.method"
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                    "operator": "less_or_equal",
                    "value": 3.0, #2.5
                    "attribute": "rcsb_entry_info.resolution_combined"
                    }
                }
                ]
            },
            "request_options": {
                "return_all_hits": True,
            },
            "return_type": result_type
            }   
    # Convert query to json object
    try:
        r = requests.post(url, json=query)
        results = r.json()
    except json.decoder.JSONDecodeError:
        return []
    # Get the list of PDBs
    pdb_list = []
    for result in results["result_set"]:
        pdb_list.append(result["identifier"])
    # Convert pdb_list (pdb_id.chain_id) to dictionary where key is the pdb_id and values are a list of chain_ids
    try:
        pdb_dict = {}
        for pdb in pdb_list:
            pdb_id = pdb.split(".")[0]
            chain_id = pdb.split(".")[1]
            if pdb_id in pdb_dict:
                pdb_dict[pdb_id].append(chain_id)
            else:
                pdb_dict[pdb_id] = [chain_id]

        return pdb_dict
    except IndexError:
        return pdb_list



def fetch_and_align_pdbs_create_pocket(list_of_pdbs, uniprot, reference_structure, reference_chain, rmsd_cutoff=3):
    pdb_obj = []
    for pdb in list_of_pdbs:
        print("fetching pdbs", pdb)
        cmd.fetch(pdb)
        cmd.split_chains(pdb)
        cmd.delete(pdb)

      
    align_errors = []
    obj_list = cmd.get_object_list()

    if reference_structure is None:
        reference_structure = obj_list[0]

    else:
        try:
            cmd.fetch(f"{reference_structure}_{reference_chain}")
        except Exception:
            pass
        
        reference_structure = f"{reference_structure}_{reference_chain}"
    

    for pdb in obj_list:
        if pdb == reference_structure:
            continue
        print("Aligning", pdb, reference_structure)
        try:
            align_score = cmd.align(pdb + " and backbone", reference_structure + " and backbone")[0]

            if align_score > rmsd_cutoff:
                print(f"Alignment score {align_score} > {rmsd_cutoff}, skipping")
                align_errors.append(pdb)

        except Exception:
            print("Error aligning", pdb, reference_structure)
            align_errors.append(pdb)


    cmd.center(reference_structure)
    
    print("The following pdbs could not be aligned, check chain names and align manually")
    print(align_errors)
    print("Please Keep in mind that the PDBs were split into chains, so contacts are lost or badly aligned structures could be different proteins!")

    cmd.center(obj_list[0])

    # Save pse with Error PDBs

    error_selection = " or ".join(align_errors)
    cmd.save(f"{uniprot}_error.pse", error_selection)

    # Delete error selection
    cmd.delete(error_selection)
    
    cmd.save(f"{uniprot}.pse")





def main(uniprot, result_type, reference_structure, reference_chain):
    pdbs = list(get_list_of_pdbs(uniprot, result_type))

    fetch_and_align_pdbs_create_pocket(pdbs, uniprot, reference_structure, reference_chain, rmsd_cutoff=3)


INPUT_UNIPROT = "XXXXXXX" # Replace with your UniProt ID
REFERENCE_STRUCTURE = None # If you want to use a specific reference structure, provide its PDB ID here
REFERENCE_CHAIN = None # If you want to use a specific reference chain, provide its chain ID here


main(INPUT_UNIPROT, result_type="polymer_instance", reference_structure=REFERENCE_STRUCTURE, reference_chain=REFERENCE_CHAIN)