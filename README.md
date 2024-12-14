# misc_xtal_stuff
A small collection of scripts for protein crystallography that have made my life somewhat easier. All scripts are in their own folder.

# Prerequisits
Well, this depends on what script you wish to run. In general, if you have all of these on your system you should be fine:
- Pymol
- Phenix
- Bioconda

If you want to use the _table1_ script and get the ligand B-factor you need to somehow load this in your _.pymolrc_:
```python
@cmd.extend
def average_b_factor(selection="sele"):
    atom_count = cmd.select(selection)
    total_b_factor = 0
    for atom in cmd.get_model(selection).atom:
        total_b_factor += atom.b
    avg_b_factor = total_b_factor/atom_count
    print(avg_b_factor)
    return avg_b_factor
```

I have only tested this on a linux system, if you encounter any issues let me know.

# Brief Overview

- Table1
  - As the name implies this generates a basic table1 (data collection and refinement statistics) in _.csv_ format.
  - It reads the data collection statistics from the _CORRECT.LP_ file using XDS (Maybe I will add XSCALE some day). The refinement data is gathered from the provided _.pdb_ file it self or from running *phenix.model_statistics*. The ligand B-Factor is retrieved from pymol.
  - You can choose if you only want data collection statistics, refinement statistics or both.
  - Usage:
    ```bash
    $ ./table1.py -c CORRECT.LP -p model.pdb -f ligand1.cif ligand2.cif -l ligand_id
    ```

- Quick Validation
  - Well it is no more than that. Any other program is better, use them. Too often I forget to mutate residues in loops I built or misclick and add the wrong amino acid... It also checks if there are any polar interactions closer than 2.2 A (excluding inorganics) and provides a brief summary of molprobity (using *phenix.model_statistics*). On the plus side, it is very fast and checks for obvious mistakes that are frustrating once you have waited for the OneDep report.
  - Usage
    ```bash
    $ ./quick_validation.py reference_fasta model.pdb -c ligand.cif
    ```

- Xtal Conditions
  - You give it a Uniprot-ID it gives you a _.csv_ with all deposited structures and some information like: Resolution, space group, conditions, ...
  - Usage
    ```bash
    $ ./get_xtal_conditions.py -uniprot UNIPROT_ID -filename OUTPUT.CSV
    ```

- Structure Contacts
  - Quickly check which residues are in contacts with symmetry mates or with other monomers in the ASU. A default cutoff of 4 A is chosen, but can be changed.
  - For pure pymol(-python) scripts I call it directly with pymol, rather than dealing with installing pymol to an environment.
  - Usage
    ```bash
    $ pymol contact_residues.py --pdb PDBID [--cutoff 4]
    ```
# Bugs, Errors and Missing Functionality
If anything does not work, is wrong or is missing let me know. If I have time I will try to correct and implement.
