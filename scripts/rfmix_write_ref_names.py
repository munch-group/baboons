import sys
import yaml
import pandas as pd

_, analysis, config_file, ref_name_file = sys.argv
with open(config_file) as f:
    config = yaml.safe_load(f)

ref_name_list = config['rfmix_analyzes']

meta_data_samples = pd.read_csv(config['sample_meta_data'], sep=" ")
ref_samples = meta_data_samples.loc[meta_data_samples.C_origin.isin(ref_name_list[analysis])]
pop_df = pd.DataFrame({"PDGP_ID": ref_samples.PGDP_ID,
                    "population": ref_samples.C_origin})
pop_df.to_csv(sample_map_file, index=False, header=False, sep="\t")


for n in ref_name_list:
    meta_data_samples_sub = meta_data_samples.loc[meta_data_samples.C_origin.isin(n[1])]
    query_samples = meta_data_samples.loc[~(meta_data_samples.C_origin.isin(n[1])) &
                                        (meta_data_samples.C_origin != "Gelada, Captive")]
    pop_df = pd.DataFrame({"PDGP_ID": meta_data_samples_sub.PGDP_ID,
                           "population": meta_data_samples_sub.C_origin})
    pop_df.to_csv(path_to_output+"/"+n[0]+"/ref_names.txt",
                  index=False, header=False, sep="\t")