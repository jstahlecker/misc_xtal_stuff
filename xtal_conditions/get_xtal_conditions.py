#!/usr/bin/env python
import urllib
import json
import pandas as pd
import argparse
import requests
#from get_fasta import get_fasta

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

def split_fasta_string(fasta_string):
    """Split a FASTA string into a list of FASTA strings."""
    fasta_list = fasta_string.split(">")[1].split("\n")[:-1][1]    

    return fasta_list


def get_fasta(pdb_id):
    """Get the FASTA sequence for a PDB ID."""
    url = f"https://www.rcsb.org/fasta/entry/{pdb_id}/display"
    with urllib.request.urlopen(url) as response:
        fastas = response.read().decode("utf-8")
    
    fasta_list = split_fasta_string(fastas)

    return fasta_list

def make_xtal_csv(pdb_list, filename):
    """Make a csv file with the crystallographic conditions for each PDB in the list"""
    RCSB_BASE = "https://data.rcsb.org/rest/v1/core/entry/"

    error_list = []
    xtal_data = []
    for pdb_id in pdb_list:
        print(pdb_id)
        page = urllib.request.urlopen(RCSB_BASE + pdb_id)

        data = json.load(page)
        try:
            try:
                # Xtal details
                xtal_details = data["exptl_crystal_grow"][0]["pdbx_details"]
            except:
                xtal_details = "-"
            try:
                xtal_temp = data["exptl_crystal_grow"][0]["temp"]
            except:
                xtal_temp = "-"
            try:
                xtal_method = data["exptl_crystal_grow"][0]["method"]
            except:
                xtal_method = "-"
            # Author Information
            author_list = []
            for author in data["audit_author"]:
                author_list.append(author["name"])

            # Citation
            try:
                citation = data["citation"][0]["pdbx_database_id_doi"]
            except:
                try:
                    citation = data["citation"][0]["journal_abbrev"]
                except:
                    citation = "-"
            try:
                #resolution = data["pdbx_vrpt_summary"]["pdbresolution"]
                resolution = data["rcsb_entry_info"]["diffrn_resolution_high"]["value"]
            except:
                resolution = "-"
            
            # Cell data
            try:
                alpha = data["cell"]["angle_alpha"]
            except:
                alpha = "-"

            try:
                beta = data["cell"]["angle_beta"]
            except:
                beta = "-"

            try:
                gamma = data["cell"]["angle_gamma"]
            except:
                gamma = "-"

            try:
                a = data["cell"]["length_a"]
            except:
                a = "-"

            try:
                b = data["cell"]["length_b"]
            except:
                b = "-"

            try:
                c = data["cell"]["length_c"]
            except:
                c = "-"

            try:
                symmetry = data["symmetry"]["space_group_name_hm"]
            except:
                symmetry = "-"

            try:
                sg = data["symmetry"]["int_tables_number"]
            except:
                sg = "-"
            
            try:
                fasta = get_fasta(pdb_id)
            except:
                fasta = "-"

        except Exception:
            error_list.append(pdb_id)

        xtal_data.append([pdb_id, resolution, symmetry, sg, ",".join([str(x) for x in [alpha, beta, gamma]]), ",".join([str(x) for x in [a,b,c]]), xtal_details,
                          xtal_temp, xtal_method, citation, ", ".join([author for author in author_list]), fasta])

    print("Errors in following IDs", error_list)

    df = pd.DataFrame(xtal_data, columns=[
                      "PDB_ID", "RESOLUTION", "SYMMETRY", "SG", "ANGLE", "LENGTH", "XTAL_DETAILS", "XTAL_TEMP", "XTAL_METHOD", "CITATION", "AUTHOR_LIST", "FASTA"])

    df.to_csv(filename, index=False)


def main(uniprot_id, filename):

    pdb_list = get_list_of_pdbs(uniprot_id, result_type="entry")
    make_xtal_csv(pdb_list, filename)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("-uniprot", help="UNIPROT ID",
                        type=str, action="store")
    parser.add_argument("-filename", help="filename",
                        type=str, action="store")

    args = parser.parse_args()
    main(args.uniprot, args.filename)
