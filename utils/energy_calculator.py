from Bio.Seq import Seq
import RNA
import json
import numpy as np
import math



#==========================================================
# This is CONFIGURATION LOADING UTILITY
def load_config(file_path='config.json'):
    with open(file_path, 'r') as config_file:
        config = json.load(config_file)
    return config
#============================================================




# ==============================================================================================
# Here are two (2) functions related to DNA --> mRNA input processings
def transcribe_DNA_to_RNA(DNA_sequence):
    print(f'Input DNA sequence = {DNA_sequence}')
    print(f"Transcribing DNA to RNA...")
    DNA_sequence = Seq(DNA_sequence.upper().replace(" ", "").replace("\n", ""))
    mRNA_sequence = DNA_sequence.transcribe()
    print("Transcription complete!")
    print(f"The resulting mRNA sequence is {mRNA_sequence}.")

    return mRNA_sequence

def parse_relevant_regions(mRNA_sequence):

    print("Looking for the start codon AUG...")
    start_codon_index=mRNA_sequence.find('AUG')

    if start_codon_index == -1:
        print("Start codon (AUG) not found in the provided mRNA sequence!")
    else: print(f"The start codon (AUG) is found in the following bp index: {start_codon_index}!")

    sd_start = max(0, start_codon_index - 20) # assumption. Can be parameterized if needed.
    sd_end = max(0, start_codon_index - 5) # assumption. Can be parameterized if needed.

    regions = {
        "mRNA_sequence": mRNA_sequence,
        "start_codon_index": start_codon_index,
        "sd_region": mRNA_sequence[sd_start:sd_end],
        "standby_region": mRNA_sequence[max(0, start_codon_index - 35):sd_start]
        }

    return regions
# =======================================================================================================





# =======================================================================================================
# Here are the five (5) functions related to ∆G computation
def predict_mRNA_structure(mRNA_sequence, region_start, region_end):
    """
    Computes ΔG_mRNA for folding a region near the RBS.
    """
    print("Computing the mRNA Structural free energy...")
    sub_seq = str(mRNA_sequence[region_start:region_end])
    print(f"The region of interest within the provided mRNA sequence is {sub_seq}.")
    structure, mfe = RNA.fold(sub_seq)
    print(f"The structure of the region of interest is {structure}.")
    print(f"The minimal folding energy (MFE) of the region of interest is {mfe:.4f} kcal/mol")
    return structure, mfe

def compute_SD_binding_energy(Shine_Dalgarno_sequence, rRNA_sequence):
    """
    Calculates ΔG of hybridization between SD and 16S rRNA 3' end.
    """
    duplex = RNA.duplexfold(Shine_Dalgarno_sequence, rRNA_sequence)
    return duplex.energy

def compute_start_codon_energy(start_codon_energy = -1.19):
    """
    ΔG_start is generally constant for AUG-tRNA interaction.
    """
    start_codon_energy = start_codon_energy
    return start_codon_energy

def compute_spacing_penalty(spacing_penalty):
    """
    Penalize suboptimal spacing between SD and start codon.
    Will have to look up literature data.
    """
    penalty = spacing_penalty

    return penalty

def compute_standby_site_energy(standby_sequence):
    """
    Predict ΔG of the standby site's folding (ribosome needs this open).
    """
    structure, mfe = RNA.fold(standby_sequence)
    return structure, mfe
#======================================================================================================



#======================================================================================================
# I have two (2) functions for you to use to transition total ∆G into translation rate 
def compute_total_deltaG(mRNA_G, SD_binding_G, start_codon_G, spacing_penalty_G, standby_site_G):
    """
    Combine individual ΔG values into a single ΔG_total.
    """

    total_G = np.sum([mRNA_G, SD_binding_G, start_codon_G, spacing_penalty_G, standby_site_G])

    return total_G

def calculate_initial_translation_rate(total_G):

    K = None # Proportionality constant. must be found from the literature. 
    beta = None # Inverse effective temperature constant. must be found from the literature for P. peutida. 

    initial_translation_rate = K * math.e ** (-1 * beta * total_G)

    return initial_translation_rate
#=================================================================================================================