#!/usr/bin/env python
import json
import pandas as pd
import argparse
import requests
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import partial

# Constants
RCSB_BASE_ENTRY = "https://data.rcsb.org/rest/v1/core/entry/"
RCSB_BASE_POLYMER = "https://data.rcsb.org/rest/v1/core/polymer_entity"
FASTA_URL_TEMPLATE = "https://www.rcsb.org/fasta/entry/{pdb_id}/display"
SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"

# Initialize a single session for all requests
session = requests.Session()

def get_list_of_pdbs(uniprot, result_type="entry"):
    """Get list of PDBs based on UniProt ID."""
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
                        "value": 3.0,  # 2.5
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

    try:
        response = session.post(SEARCH_URL, json=query)
        response.raise_for_status()
        results = response.json()
    except (requests.RequestException, json.JSONDecodeError) as e:
        print(f"Error fetching PDB list: {e}")
        return {}

    pdb_list = [result["identifier"] for result in results.get("result_set", [])]

    pdb_dict = {}
    for pdb in pdb_list:
        try:
            pdb_id, chain_id = pdb.split(".")
            pdb_dict.setdefault(pdb_id, []).append(chain_id)
        except ValueError:
            pdb_dict.setdefault(pdb, [])

    return pdb_dict

def get_fasta(pdb_id):
    """Get the FASTA sequence for a PDB ID."""
    url = FASTA_URL_TEMPLATE.format(pdb_id=pdb_id)
    try:
        response = session.get(url)
        response.raise_for_status()
        fastas = response.text
        fasta_list = [line for line in fastas.split('\n') if line and not line.startswith('>')]
        return ''.join(fasta_list)
    except requests.RequestException:
        return "-"

def get_entity_names(pdb_id):
    """Gets the names of the entities for a PDB ID."""
    url = f"{RCSB_BASE_ENTRY}{pdb_id}"
    try:
        response = session.get(url)
        response.raise_for_status()
        data = response.json()
        return data["rcsb_entry_container_identifiers"]["polymer_entity_ids"]
    except (requests.RequestException, KeyError):
        return []

def get_expression_system(pdb_id):
    """Gets the expression system for a PDB ID."""
    entity_names = get_entity_names(pdb_id)
    if not entity_names:
        return "-"

    expression_systems = set()
    for entity in entity_names:
        url = f"{RCSB_BASE_POLYMER}/{pdb_id}/{entity}"
        try:
            response = session.get(url)
            response.raise_for_status()
            data = response.json()
            expression = data.get("rcsb_entity_host_organism", [{}])[0].get("ncbi_scientific_name", "-")
            expression_systems.add(expression)
        except (requests.RequestException, IndexError, KeyError):
            expression_systems.add("-")
    return ", ".join(expression_systems) if expression_systems else "-"

def make_xtal_csv(pdb_dict, filename, max_workers=10):
    """Make a CSV file with the crystallographic conditions for each PDB in the list."""
    RCSB_BASE_ENTRY = "https://data.rcsb.org/rest/v1/core/entry/"
    xtal_data = []
    error_list = []

    def process_pdb(pdb_id):
        try:
            response = session.get(f"{RCSB_BASE_ENTRY}{pdb_id}")
            response.raise_for_status()
            data = response.json()

            # Xtal details
            xtal_details = data.get("exptl_crystal_grow", [{}])[0].get("pdbx_details", "-")
            xtal_temp = data.get("exptl_crystal_grow", [{}])[0].get("temp", "-")
            xtal_method = data.get("exptl_crystal_grow", [{}])[0].get("method", "-")

            # Author Information
            author_list = [author.get("name", "") for author in data.get("audit_author", [])]

            # Citation
            citation = data.get("citation", [{}])[0].get("pdbx_database_id_doi") or \
                       data.get("citation", [{}])[0].get("journal_abbrev", "-")

            # Resolution
            resolution = data.get("rcsb_entry_info", {}).get("diffrn_resolution_high", {}).get("value", "-")

            # Cell data
            cell = data.get("cell", {})
            alpha = cell.get("angle_alpha", "-")
            beta = cell.get("angle_beta", "-")
            gamma = cell.get("angle_gamma", "-")
            a = cell.get("length_a", "-")
            b = cell.get("length_b", "-")
            c = cell.get("length_c", "-")

            # Symmetry
            symmetry = data.get("symmetry", {}).get("space_group_name_hm", "-")
            sg = data.get("symmetry", {}).get("int_tables_number", "-")

            # FASTA
            fasta = get_fasta(pdb_id)

            # Expression System
            expression_system = get_expression_system(pdb_id)

            return [
                pdb_id,
                expression_system,
                resolution,
                symmetry,
                sg,
                f"{alpha},{beta},{gamma}",
                f"{a},{b},{c}",
                xtal_details,
                xtal_temp,
                xtal_method,
                citation,
                ", ".join(author_list),
                fasta
            ]
        except Exception as e:
            print(f"Error processing PDB ID {pdb_id}: {e}")
            return None

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_pdb = {executor.submit(process_pdb, pdb_id): pdb_id for pdb_id in pdb_dict}
        for future in as_completed(future_to_pdb):
            pdb_id = future_to_pdb[future]
            result = future.result()
            if result:
                xtal_data.append(result)
            else:
                error_list.append(pdb_id)

    if error_list:
        print("Errors in the following IDs:", error_list)

    columns = [
        "PDB_ID",
        "EXPRESSION_SYSTEM",
        "RESOLUTION",
        "SYMMETRY",
        "SG",
        "ANGLE",
        "LENGTH",
        "XTAL_DETAILS",
        "XTAL_TEMP",
        "XTAL_METHOD",
        "CITATION",
        "AUTHOR_LIST",
        "FASTA"
    ]

    df = pd.DataFrame(xtal_data, columns=columns)
    df = df.sort_values(by='PDB_ID').reset_index(drop=True)
    df.to_csv(filename, index=False)
    print(f"CSV file '{filename}' created successfully with {len(xtal_data)} entries.")

def main(uniprot_id, filename):
    pdb_dict = get_list_of_pdbs(uniprot_id, result_type="entry")
    if not pdb_dict:
        print("No PDB IDs found for the given UniProt ID.")
        return
    make_xtal_csv(pdb_dict, filename)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch crystallographic data for PDB IDs based on a UniProt ID.")
    parser.add_argument("-uniprot", help="UniProt ID", type=str, required=True)
    parser.add_argument("-filename", help="Output CSV filename", type=str, required=True)
    args = parser.parse_args()
    main(args.uniprot, args.filename)