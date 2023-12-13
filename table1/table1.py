#!/usr/bin/env python
import subprocess
import argparse

class Table1:

    def __init__(self, xds_correct=None, pdb_file=None, cif_files=None, ligand_name=None):
        self.dataCollection = {"Wavelength": None,
                               "Space Group": None,
                               "Cell Dimensions": {"a": None, "b": None, "c": None,
                                                   "alpha": None, "beta": None, "gamma": None},
                                "Resolution Range": {"Total": {"High": None, "Low": None}, "High": {"High": None, "Low": None}},
                                "Redundancy": {"Total": None, "High": None},
                                "Completeness": {"Total": None, "High": None},
                                "Mean I/sigma(I)": {"Total": None, "High": None},
                                "Rmeas": {"Total": None, "High": None},
                                "CC half": {"Total": None, "High": None},
                                "Wilson B factor": None}
        self.refinement = {"Resolution Included": {"Low": None, "High": None},
                           "Rwork/Rfree": {"Work": None, "Free": None},
                           "Bond RMSD": None,
                            "Angle RMSD": None,
                            "Ramachandran": {"Favored": None, "Allowed": None, "Outliers": None},
                            "Rotamer Outliers": None,
                            "All-atom Clashscore": None,
                            "Average B factor": {"Overall": None, "Protein": None, "Ligand": None, "Water": None}}
        
        self.xds_correct = xds_correct
        self.pdb_file = pdb_file
        self.cif_files = cif_files
        self.ligand_name = ligand_name
    
    def get_data_collection(self):
        """Read all statistics from CORRECT.LP file"""
        final_index = None
        with open(self.xds_correct, "r") as correct:
            correct_content = correct.readlines()
            for index, line in enumerate(correct_content):
                if line.strip().startswith("X-RAY_WAVELENGTH"):
                    line_list = line.strip().split()
                    self.dataCollection["Wavelength"] = float(line_list[1])
            
                if line.strip().startswith("SPACE_GROUP_NUMBER"):
                    line_list = line.strip().split()
                    self.dataCollection["Space Group"] = int(line_list[1])
                
                if line.strip().startswith("UNIT_CELL_CONSTANTS"):
                    line_list = line.strip().split()
                    self.dataCollection["Cell Dimensions"]["a"] = float(line_list[1])
                    self.dataCollection["Cell Dimensions"]["b"] = float(line_list[2])
                    self.dataCollection["Cell Dimensions"]["c"] = float(line_list[3])
                    self.dataCollection["Cell Dimensions"]["alpha"] = float(line_list[4])
                    self.dataCollection["Cell Dimensions"]["beta"] = float(line_list[5])
                    self.dataCollection["Cell Dimensions"]["gamma"] = float(line_list[6])
                
                if line.strip().startswith("INCLUDE_RESOLUTION_RANGE"):
                    line_list = line.strip().split()
                    self.dataCollection["Resolution Range"]["Total"]["Low"] = float(line_list[1])
                    self.dataCollection["Resolution Range"]["Total"]["High"] = float(line_list[2])

                if line.strip().startswith('STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL"'):
                    final_index = index
                
                if line.strip().startswith("WILSON LINE (using all data)"):
                    line_list = line.strip().split()
                    self.dataCollection["Wilson B factor"] = float(line_list[9])
            
            for index, line in enumerate(correct_content[final_index:]):
                
                if line.strip().startswith("total"):
                    total_line = line.strip().split()
                    hres_line = correct_content[final_index + index - 1].strip().split()

                    self.dataCollection["Resolution Range"]["High"]["Low"] = float(correct_content[final_index + index - 2].strip().split()[0])
                    self.dataCollection["Resolution Range"]["High"]["High"] = float(hres_line[0])

                    self.dataCollection["Redundancy"]["Total"] = float(total_line[1]) / float(total_line[2])
                    self.dataCollection["Redundancy"]["High"] = float(hres_line[1]) / float(hres_line[2])

                    self.dataCollection["Completeness"]["Total"] = float(total_line[4][:-1])
                    self.dataCollection["Completeness"]["High"] = float(hres_line[4][:-1])

                    self.dataCollection["Mean I/sigma(I)"]["Total"] = float(total_line[8])
                    self.dataCollection["Mean I/sigma(I)"]["High"] = float(hres_line[8])

                    self.dataCollection["Rmeas"]["Total"] = float(total_line[9][:-1])
                    self.dataCollection["Rmeas"]["High"] = float(hres_line[9][:-1])

                    try:
                        self.dataCollection["CC half"]["Total"] = float(total_line[10])
                        self.dataCollection["CC half"]["High"] = float(hres_line[10])
                    except ValueError:
                        self.dataCollection["CC half"]["Total"] = float(total_line[10][:-1])
                        self.dataCollection["CC half"]["High"] = float(hres_line[10][:-1])        
        
    def get_refinement(self):
        """Read all statistics from refinement"""

        # Get resolution from pdb file
        
        with open(self.pdb_file, "r") as pdb:
            for line in pdb:
                if "RESOLUTION RANGE HIGH" in line:
                    self.refinement["Resolution Included"]["High"] = float(line.strip().split()[-1])
                if "RESOLUTION RANGE LOW" in line:
                    self.refinement["Resolution Included"]["Low"] = float(line.strip().split()[-1])
                if "R VALUE            (WORKING SET) :" in line.strip():
                    self.refinement["Rwork/Rfree"]["Work"] = float(line.strip().split()[-1])*100
                if "FREE R VALUE                     :" in line:
                    self.refinement["Rwork/Rfree"]["Free"] = float(line.strip().split()[-1])*100


        if self.cif_files:
            cif_str = " ".join(str(cif) for cif in self.cif_files)
            ref_stats_command = f"phenix.model_statistics {self.pdb_file} {cif_str}"
        else:
            ref_stats_command = f"phenix.model_statistics {self.pdb_file}"
        ref_stats = subprocess.run(ref_stats_command, shell=True, capture_output=True, text=True) # Run phenix.model_statistics

        for index, line in enumerate(ref_stats.stdout.splitlines()): # Get all parameters from output
            if line.strip().startswith("Bond"):
                line_list = line.strip().split()
                self.refinement["Bond RMSD"] = float(line_list[2])
            if line.strip().startswith("Angle"):
                line_list = line.strip().split()
                self.refinement["Angle RMSD"] = float(line_list[2])
            if line.strip().startswith("Ramachandran Plot"):
                outl_index = index + 1
                allowed_index = index + 2
                favored_index = index + 3

                outl_list = ref_stats.stdout.splitlines()[outl_index].strip().split()
                self.refinement["Ramachandran"]["Outliers"] = float(outl_list[2])
                    
                allowed_list = ref_stats.stdout.splitlines()[allowed_index].strip().split()
                self.refinement["Ramachandran"]["Allowed"] = float(allowed_list[2])

                favored_list = ref_stats.stdout.splitlines()[favored_index].strip().split()
                self.refinement["Ramachandran"]["Favored"] = float(favored_list[2])

            if line.strip().startswith("Rotamer Outliers"):
                line_list = line.strip().split()
                self.refinement["Rotamer Outliers"] = float(line_list[3])
            if line.strip().startswith("All-atom Clashscore"):
                line_list = line.strip().split()
                self.refinement["All-atom Clashscore"] = float(line_list[-1])
            
            if line.strip().startswith("Overall:"):
                line_list = line.strip().split()
                self.refinement["Average B factor"]["Overall"] = float(line_list[3])
            if line.strip().startswith("Protein:"):
                line_list = line.strip().split()
                self.refinement["Average B factor"]["Protein"] = float(line_list[3])
            if line.strip().startswith("Water:"):
                line_list = line.strip().split()
                self.refinement["Average B factor"]["Water"] = float(line_list[3])
    
    def get_b_factor(self):

        mini_script = f"cmd.load('{self.pdb_file}'); cmd.select('dummy', 'resn {self.ligand_name}'); cmd.do('average_b_factor dummy')"
        b_factor = subprocess.run(["pymol", "-cd", mini_script], check=True, capture_output=True, text=True)
        self.refinement["Average B factor"]["Ligand"] = float(b_factor.stdout.splitlines()[-1])

    def make_statistics(self):
        """Run all methods to get statistics"""
        if self.xds_correct:
            self.get_data_collection()
        if self.pdb_file:
            self.get_refinement()
        if self.ligand_name:
            self.get_b_factor()


    def make_output_formating(self):
        """Create output string for table1"""

        if self.xds_correct:

            data_collection_out = [ ("Wavelength", f'{self.dataCollection["Wavelength"]:.2f}'),
                                    ("Space Group", f'{self.dataCollection["Space Group"]}'),
                                    ("Cell Dimensions", " "),
                                    ("a / b / c", f'{self.dataCollection["Cell Dimensions"]["a"]:.2f} / {self.dataCollection["Cell Dimensions"]["b"]:.2f} / {self.dataCollection["Cell Dimensions"]["c"]:.2f}'),
                                    ("alpha / beta / gamma", f'{self.dataCollection["Cell Dimensions"]["alpha"]:.2f} / {self.dataCollection["Cell Dimensions"]["beta"]:.2f} / {self.dataCollection["Cell Dimensions"]["gamma"]:.2f}'),
                                    ("Resolution Range", f'{self.dataCollection["Resolution Range"]["Total"]["Low"]:.2f} - {self.dataCollection["Resolution Range"]["Total"]["High"]:.2f} ({self.dataCollection["Resolution Range"]["High"]["Low"]:.2f} - {self.dataCollection["Resolution Range"]["High"]["High"]:.2f})'),
                                    ("Redundancy", f'{self.dataCollection["Redundancy"]["Total"]:.2f} ({self.dataCollection["Redundancy"]["High"]:.2f})'),
                                    ("Completeness", f'{self.dataCollection["Completeness"]["Total"]:.2f} ({self.dataCollection["Completeness"]["High"]:.2f})'),
                                    ("Mean I/sigma(I)", f'{self.dataCollection["Mean I/sigma(I)"]["Total"]:.2f} ({self.dataCollection["Mean I/sigma(I)"]["High"]:.2f})'),
                                    ("Rmeas", f'{self.dataCollection["Rmeas"]["Total"]:.2f} ({self.dataCollection["Rmeas"]["High"]:.2f})'),
                                    ("CC half", f'{self.dataCollection["CC half"]["Total"]:.2f} ({self.dataCollection["CC half"]["High"]:.2f})'),
                                    ("Wilson B factor", f'{self.dataCollection["Wilson B factor"]:.2f}'),
                                    ("", "")]
        
        else:
            data_collection_out = []

        if self.pdb_file:
            refinement_out= [   ("Resolution Included", f'{self.refinement["Resolution Included"]["Low"]:.2f} - {self.refinement["Resolution Included"]["High"]:.2f}'),
                                ("Rwork/Rfree", f'{self.refinement["Rwork/Rfree"]["Work"]:.2f} / {self.refinement["Rwork/Rfree"]["Free"]:.2f}'),
                                ("Bond RMSD", f'{self.refinement["Bond RMSD"]:.3f}'),
                                ("Angle RMSD", f'{self.refinement["Angle RMSD"]:.3f}'),
                                ("Ramachandran", " "),
                                ("Favored / Allowed / Outliers", f'{self.refinement["Ramachandran"]["Favored"]:.2f} / {self.refinement["Ramachandran"]["Allowed"]:.2f} / {self.refinement["Ramachandran"]["Outliers"]:.2f}'),
                                ("Rotamer Outliers", f'{self.refinement["Rotamer Outliers"]:.2f}'),
                                ("All Atom Clashscore", f'{self.refinement["All-atom Clashscore"]:.2f}'),
                                ("Average B factor", " "),
                                ("Overall", f'{self.refinement["Average B factor"]["Overall"]:.2f}'),
                                ("Protein", f'{self.refinement["Average B factor"]["Protein"]:.2f}'),
                                ("Water", f'{self.refinement["Average B factor"]["Water"]:.2f}')]
        if self.ligand_name and self.pdb_file:
            refinement_out.append(("Ligand", f'{self.refinement["Average B factor"]["Ligand"]:.2f}'))
        else:
            refinement_out = []

        total_out = data_collection_out + refinement_out

        return total_out  


    def create_out_file(self):
        """Create table1.csv file"""
        output = self.make_output_formating()

        with open("table1.csv", "w") as out_file:
            for line in output:
                out_file.write(line[0] + "," + line[1] + "\n")
        



        



if __name__ == "__main__":
    
    # Create parser
    parser = argparse.ArgumentParser(description="Create table1.csv file from CORRECT.LP and refinement statistics")

    # Add arguments
    # If correct file is given, run data collection statistics
    parser.add_argument("-c", "--correct", help="CORRECT.LP file from XDS", required=False)
    # If pdb file is given, run refinement statistics
    parser.add_argument("-p", "--pdb", help="PDB file from refinement", required=False)
    # cif files are optional
    parser.add_argument("-f", "--cif", help="CIF files from refinement", required=False, nargs="+")
    parser.add_argument("-l", "--ligand", help="Ligand name", required=False)

    # Parse arguments
    args = parser.parse_args()

    # Create table1 object
    table1 = Table1(xds_correct=args.correct, pdb_file=args.pdb, cif_files=args.cif, ligand_name=args.ligand)
    table1.make_statistics()
    table1.create_out_file()


