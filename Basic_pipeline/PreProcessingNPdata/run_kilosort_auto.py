from kilosort import run_kilosort
#from SGLXMetaToCoords2 import MetaToCoords
from read_chmap_ephys import read_meta, get_probe
from pathlib import Path
from kilosort.io import save_probe
from kilosort.io import load_probe
import urllib.request
from pathlib import Path
from tqdm import tqdm
import pandas as pd
import os


TABLE_OF_EXPERIMENT = '//sil3/data/Large_scale_mapping_NP/Experiment_Excel.xlsx'

def get_table_path(Chosen_Animal, Chosen_Insertion):
    exps_info = pd.read_excel(TABLE_OF_EXPERIMENT)
    data_info = exps_info[(exps_info.Animal_ID == Chosen_Animal) & (exps_info.Insertion == Chosen_Insertion)].iloc[0].to_dict()

    return data_info['Base_path'],  data_info['Exp_name'], data_info['Insertion']

def _get_meta_path(Base_path, Exp_name, Insertion):
    Base_path = Base_path.replace('\\', '/')
    file_path = (
        f"/{Base_path[1:]}/{Exp_name}/Insertion{Insertion}/"
        f"catgt_{Exp_name}_{Insertion}_g0/{Exp_name}_{Insertion}_g0_tcat.imec0.ap.meta"
    )
    return(file_path)

def _get_bin_path(Base_path, Exp_name, Insertion):
    Base_path = Base_path.replace('\\', '/')
    file_path = (
        f"{Exp_name}_{Insertion}_g0_tcat.imec0.ap.bin"
    )
    return(file_path)

def _get_path(Base_path, Exp_name, Insertion):
    Base_path = Base_path.replace('\\', '/')
    file_path = (
        f"/{Base_path[1:]}/{Exp_name}/Insertion{Insertion}/"
        f"catgt_{Exp_name}_{Insertion}_g0/"
    )
    return(file_path)



if __name__ == "__main__":
    animal = 'PV132'

    for i in  range(6, 7):
        meta_path = Path(_get_meta_path(*get_table_path(animal, i)))
        base_path, animal_name, insertion_name = get_table_path(animal, i)
        path_folder = _get_path(base_path, animal_name, insertion_name)
        path_bin = _get_bin_path(base_path, animal_name, insertion_name)
        print(path_folder)
        print(path_bin)
        meta_ephys = read_meta(path_folder+path_bin)
        probe_dict = get_probe(meta_ephys)
        os.chdir(_get_path(*get_table_path(animal, i)))

        SAVE_PATH = _get_bin_path(*get_table_path(animal, i))

        settings = {'filename': SAVE_PATH, 'n_chan_bin': 385, 'Th_universal': 8, 'Th_learned': 6, 'duplicate_spike_ms': 0.5, 'highpass_cutoff': 300}

        ops, st, clu, tF, Wall, similar_templates, is_ref, est_contam_rate, kept_spikes = \
            run_kilosort(
                settings=settings, probe=probe_dict)