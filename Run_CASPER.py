"""This is a python file where you can run the following CASPER functions: On-target search, Off-target analysis,
    Multitargeting, and population comparisons.  See the comments in each section for details on how to setup and run
    each script."""

# ==== IMPORT SETTINGS ==== #
from CASPERQuick import CasperQuick
from OffTarget import OffTargetAlgorithm
from multitargeting import Multitargeting
from Comparison import Compare_Orgs

# ---------- OPTIONS FOR RUN ----------- #
# __CASPERQuick__
# __OffTarget__
# __Multitargeting__
# __Populations__
# -------------------------------------- #

run = "__OffTarget__"

# ---------- GENERAL SETTINGS ------------------ #
# These settings are needed for all Analyses
CASPER_Seq_Finder_files_directory = "/Volumes/Seagate_Drive/CrisprDB/"
output_file_path = "/Users/brianmendoza/Desktop/Sequences/"
endonuclease = 'spCas9'  # This is for the output file name
base_organism_code = "sce"  # use the KEGG code that was used to name your CASPER_Seq_Finder file

# ---------- SPECIFIC SETTINGS ------------------ #
# Settings specific for CASPERQuick:
regions_or_kegg_codes = ["sce:YBR043C"]
off_target_all = False

# Settings specific for OffTarget:
casperOffList_file_path = "CASPEROfflist.txt"
other_orgs_off = []  # This list is if you want to check off targets against a population of organisms
threshold_score = 0.001  # Off target scores are between 0 and 1 with 1 being a full match.  Threshold > 0.2 is recommended

# Settings specific for Multitargeting:
# None additional parameters needed.

# Settings specific for Populations:
other_orgs = ["ctx", "cace"]  # These are the additional organisms. Their files must be in the same directory as your base organism

# ================ CODE EXECUTION. USER MAY IGNORE BELOW ================= #
# ======================================================================== #
csf_file = CASPER_Seq_Finder_files_directory + base_organism_code + endonuclease + ".cspr"
if run == "__CASPERQuick__":
    C = CasperQuick(csf_file,output_file_path,off_target_all)
    C.loadGenesandTargets(regions_or_kegg_codes)
elif run == "__OffTarget__":
    O = OffTargetAlgorithm(threshold_score, endonuclease, base_organism_code, csf_file,
                           other_orgs_off, casperOffList_file_path, output_file_path)
elif run == "__Multitargeting__":
    M = Multitargeting(csf_file, output_file_path)
elif run == "__Populations__":
    P = Compare_Orgs(output_file_path, csf_file, base_organism_code, endonuclease, other_orgs)
else:
    print("ERROR: typo in run string object!")
