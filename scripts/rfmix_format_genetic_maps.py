import sys
import yaml
import pandas as pd


_, config_file, aut_genetic_map_file, x_genetic_map_file = sys.argv

with open(config_file) as f:
    config = yaml.safe_load(f)

df_l = []
for a in range(1, 21):
    genetic_map = config['rec_maps'][f'chr{a}']
    recomb_df = pd.read_csv(genetic_map, sep=" ")
    df_l.append(recomb_df)
pd.concat(df_l)[["chromosome", "position", "Genetic_Map(cM)"]].to_csv(aut_genetic_map_file, sep=" ", index=False)

df_l = []
for a in ["X"]:
    genetic_map = config['rec_maps'][f'chr{a}']
    recomb_df = pd.read_csv(genetic_map.format(a), sep=" ")
    df_l.append(recomb_df)
pd.concat(df_l)[["chromosome", "position", "Genetic_Map(cM)"]].to_csv(x_genetic_map_file, sep=" ", index=False)

