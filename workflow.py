import os
import pandas as pd
import yaml
from gwf.workflow import collect
from gwf import Workflow, AnonymousTarget

# get all targets with a given key
def target_output_files(targets, key):
    return [out for target in targets[key] for out in target.outputs]

###############################################################################
## Top workflow for all the stuff required to run the individual pioe lines
###############################################################################

gwf = Workflow()

# read config file for workflow
with open('workflow_config.yml') as f:
    config = yaml.safe_load(f)

# rfmix workflow
from rfmix.workflow import rfmix_workflow

# controls wheather submodule workflows merged with main workflow
# or run in isolation (only use False for initial submoduile setup)
merge_workflows = True

###############################################################################
##  Generate input files and run RFmix submodule workflow
###############################################################################

# full_list = ['Cynocephalus, Central Tanzania', 'Anubis, Kenya', 'Kindae, Zambia',
#     'Hamadryas, Ethiopia', 'Anubis, Tanzania',
#     'Cynocephalus, Western Tanzania', 'Papio, Senegal', 'Ursinus, Zambia',
#     'Anubis, Ethiopia']

# rfmix analyses
rfmix_analyses = config['rfmix_analyzes']

# rfmix output dir
rfmix_output_dir = "steps/rfmix_gen100/"

# compile reference/query sample lists for rfmix
meta_data_samples = pd.read_csv(config['sample_meta_data'], sep=" ")
analyzes = []
for analysis in rfmix_analyses:
    d = {}
    d["analysis"] = analysis
    os.makedirs(rfmix_output_dir+"/"+analysis, exist_ok=True)
    ref_samples = meta_data_samples.loc[meta_data_samples.C_origin.isin(rfmix_analyses[analysis])]
    query_samples = meta_data_samples.loc[~(meta_data_samples.C_origin.isin(rfmix_analyses[analysis])) &
                                        (meta_data_samples.C_origin != "Gelada, Captive")]
    d["ref_samples"] =list(ref_samples.PGDP_ID)
    d["query_samples"] = list(query_samples.PGDP_ID)
    analyzes.append(d)

# write sample/population info
for analysis in rfmix_analyses:
    analysis_dir = rfmix_output_dir + "/" + analysis
    sample_map_file = analysis_dir + "/ref_names.txt"
    gwf.target(f'sample_map_{analysis}', inputs=['workflow_config.yml'], outputs=[rfmix_output_dir+"/"+analysis+"/ref_names.txt"]) << f'''

    mkdir -p analysis_dir
    python scripts/rfmix_write_sample_map.py {analysis} workflow_config.yml {sample_map_file}
    '''

# write recombination maps
autosome_rec_map = rfmix_output_dir + "aut_genetic_map.txt"
x_rec_map = rfmix_output_dir + "X_genetic_map.txt"
gwf.target('format_genetic_maps', 
           memory='16gb',
           inputs=['workflow_config.yml'], 
           outputs=[autosome_rec_map, x_rec_map]) << f'''

mkdir -p output_dir
python scripts/rfmix_format_genetic_maps.py workflow_config.yml {autosome_rec_map} {x_rec_map}
'''



# run the rfmix pipeline
_gwf, rfmix_targets  = rfmix_workflow(
                          gwf if merge_workflows else Workflow(working_dir=os.getcwd()),
#                          merge_workflows and gwf or Workflow(working_dir=os.getcwd()),
                          analyzes=analyzes, 
                          output_dir=rfmix_output_dir, 
                          vcf_files=config['vcf_files'], 
                          autosome_rec_map=autosome_rec_map,
                          x_rec_map=x_rec_map
                          )

globals()['rfmix'] = _gwf

###############################################################################
## Next workflow...
###############################################################################

# # get relevant outputs from A for input to B
# input_files =  target_output_files(rfmix_targets, 'work')

