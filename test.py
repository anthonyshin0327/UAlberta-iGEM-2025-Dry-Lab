from utils.energy_calculator import *

# CONFIGURATION =============================
config=load_config(file_path='config.json')

DNA=config['DNA_sequence']
# All the parameters below that's mentioned to be "None" must be defined in config.json, just like DNA above.
# ===========================================


def main_test():
    mRNA_sequence = transcribe_DNA_to_RNA(DNA)
    regions = parse_relevant_regions(mRNA_sequence)

    mRNA_G = predict_mRNA_structure(
        mRNA_sequence=None,
        region_start=None,
        region_end=None
    )
    SD_binding_G = compute_SD_binding_energy(
        Shine_Dalgarno_sequence=None,
        rRNA_sequence=None
    )
    start_codon_G = compute_start_codon_energy(
        start_codon_energy=None
    )
    spacing_penalty_G = compute_spacing_penalty(
        spacing_penalty=None
    )
    standby_site_G = compute_standby_site_energy(
        standby_sequence=None
    )


    total_delta_G = compute_total_deltaG(
        mRNA_G=mRNA_G,
        SD_binding_G=SD_binding_G,
        start_codon_G=start_codon_G,
        spacing_penalty_G=spacing_penalty_G,
        standby_site_G=standby_site_G
    )

    initial_translation_rate = calculate_initial_translation_rate(total_G=total_delta_G)

    return initial_translation_rate

if __name__ == '__main__':
    main_test()

