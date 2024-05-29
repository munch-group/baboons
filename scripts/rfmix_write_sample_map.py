import sys
import yaml
import pandas as pd

_, analysis, config_file, sample_map_file = sys.argv
with open(config_file) as f:
    config = yaml.safe_load(f)

ref_name_list = config['rfmix_analyzes']

meta_data_samples = pd.read_csv(config['sample_meta_data'], sep=" ")
ref_samples = meta_data_samples.loc[meta_data_samples.C_origin.isin(ref_name_list[analysis])]
pop_df = pd.DataFrame({"PDGP_ID": ref_samples.PGDP_ID,
                    "population": ref_samples.C_origin})
pop_df.to_csv(sample_map_file, index=False, header=False, sep="\t")