# baboons
Master repo for baboon project

## Repo structure

Either use as library repo with with `<package_name>` and `docs` (and maybe `notebook`) dirs. Or as data analysis repo with `scripts`, `notebooks`, `data`, `steps`, `results`, `figures`dirs.

## Conda environment

Create empty environment:

```
conda create -n <name>
```

Activate environment
```
conda activate <name>
```

Add channels first to make the appear in the exported environment:
Add any channels first to make the appear in the exported environment (you can always add more later):
```
conda config --env --append channels conda-forge
conda config --env --append channels bioconda
conda config --env --append channels etetoolkit
conda config --env --append channels gwforg
conda config --env --append channels kaspermunch
```

Install a subset of packages:
```
conda install jupyterlab pandas seaborn ...
conda install gwf jupyterlab pandas seaborn ipympl statsmodels bioconda pyfaidx tskit ete3 scikit-allel samtools bamtools vcftools
```

Export it to the binder folder:
```
conda env export --from-history > binder/environment.yml
```

List revisions to envirionment:
```
conda list --revisions
```

Restore a past revision:
```
conda install --revision 2
```

Once you land on stable environment this can make it more stable:

```
conda env export --from-history > binder/environment.yml
conda deactivate <name>
conda env remove -n <name>
conda env create -n <name> -f binder/environment.yml
```

## Create PyPi and Conda packages

```
conda install conda-build twine anaconda-client
```

```
bash pypi.sh
```

```
bash conda.sh [-c channel -c channel ...]
```


## Using submodules in Git

Submodules are useful if you have analysis components that are used across many projects. It would good to have repos to use for subprojects for these components (between us, I think we have them)

- HapMix
- RFmix
- Mosaic
- Relate
- ARGweaver
- ARGweaver Clues
- Relate Clues
- SMC++
- LDhat/helmet
- pyro
- GATK mapping and base-calling
- ShapeIt phasing

The repos can each have their own adjustable workflow as decribed [here](https://gwf.app/guide/patterns/#dynamically-generating-a-workflow), that are then easily put together in a larger workflow as described [here](https://gwf.app/guide/patterns/#large-workflows).

### Preliminaries

Make sure you install the current version of Git in your environments. The one on the cluster is ancient.

    conda install -c conda-forge git

I also set these configs for each environment to get nicer/safer commands (commands below assume these are set):

    git config --global diff.submodule log
    git config status.submodulesummary 1
    git config push.recurseSubmodules check

Now, say you have a project repo called "umbrella" that will contain other projects and that you have cloned that:

    git clone git@github.com:kaspermunch/umbrella.git

### Add a submodule

clone repository as submodule:

    git submodule add git@github.com:kaspermunch/tester.git

and pull the current state of the submodule repo:

    git submodule init tester
    git submodule update tester

This also generates a `.gitmodules` configuration file that git uses to keep track of submodules. Commit that the addalong with the submodule:

    git add .gitmodules tester
    git commit -m 'Added tester repo as submodule'
    git push

If you want to work on/change submodule repo you need to check out a branch to work on (main or some other). Always do this. If you decide to make changes later and forgot you did not check out a branch you could loose those changes:

    cd tester
    git checkout main  # (or some other branch)

### Making changes to the submodule

now you can then do some work on the tester repo (E.g. change the README.md) and add, commit as usual:

    cd tester
    # change README.md
    git add README.md
    git commit

### Publishing submodule changes to GitHub

To publish your submodule commit to the tester repo on GitHub you run:

    cd tester
    git push

### Getting submodule changes from GitHub

if you run "git pull" in the umbrella repo, you pull upstream changes to the umbrella repo including the recorded state (commit) of the tester submodule:

    git pull

but it does not pull the tester submodule itself. To do that you run pull in the submodule:

    cd tester
    git pull

## Branches and multiple contributors

As you could see above submodules are really just repos inside other repos. The parent repo just treats the submodule as a file holding which state (commit) the submodule is in. So for most use, working alone on project it is quite simple. In some cases however, you need some of the `git submodule` commands. Those commands updates the relationship between the submodule state recorded in the parent repo and the state and of the submodule.

    git submodule update --checkout
    git submodule update --rebase
    git submodule update --merge
    
**--checkout** Checkout the commit recorded in the superproject on a detached HEAD in the submodule. This is the default behavior, the main use of this option is to override submodule.$name.update when set to a value other than checkout.

**--merge** Merge the commit recorded in the superproject into the current branch of the submodule. If this option is given, the submodule’s HEAD will not be detached. If a merge failure prevents this process, you will have to resolve the resulting conflicts within the submodule with the usual conflict resolution tools. 

**--rebase** Rebase the current branch onto the commit recorded in the superproject. If this option is given, the submodule’s HEAD will not be detached. If a merge failure prevents this process, you will have to resolve these failures with git-rebase[1].     
    
    git submodule update --remote --checkout
    git submodule update --remote --rebase
    git submodule update --remote --merge       

The last one fetches upsteam changes and pulls the submoduile branch recorded in the parent repo when the submodule was added. If this branch is the current branch in the local submodule, then the command is equivalent to `git pull` in the submodule.

**--remote** Instead of using the superproject’s recorded SHA-1 to update the submodule, use the status of the submodule’s remote-tracking branch. In order to ensure a current tracking branch state, update --remote fetches the submodule’s remote repository before calculating the SHA-1. If you don’t want to fetch, you should use submodule update --remote --no-fetch.

Use this option to integrate changes from the upstream subproject with your submodule’s current HEAD. Alternatively, you can run git pull from the submodule, which is equivalent except for the remote branch name: update --remote uses the default upstream repository and submodule.<name>.branch, while git pull uses the submodule’s branch.<name>.merge. Prefer submodule.<name>.branch if you want to distribute the default upstream branch with the superproject and branch.<name>.merge if you want a more native feel while working in the submodule itself.




(read the difference between merge and rebase [here](https://www.atlassian.com/git/tutorials/merging-vs-rebasing))


## Multiple submodules

You can have as many submodules as you want. With more submodules, each update command updates all submodules.

# Combining GWF workflows from multiple Git submodules

In each submodule the `workflow.py` could look like this:

```python
import os.path
import os
from collections import defaultdict
from gwf import Workflow

def submoduleA_workflow(working_dir=os.getcwd(), input_files=None, output_dir=None, summarize=True):

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # dict of targets as info for other workflows
    targets = defaultdict(list)

    gwf = Workflow(working_dir=working_dir)

    work_output = os.path.join(output_dir, 'A_output1.txt')
    target = gwf.target(
        name='A_work',
        inputs=input_files,
        outputs=[work_output],
    ) << f"""
    touch {work_output}
    """
    targets['work'].append(target)

    if summarize:
        summary_output = os.path.join(output_dir, 'A_output2.txt')
        target = gwf.target(
            name='A_summary',
            inputs=[work_output],
            outputs=[summary_output]
        ) << f"""
        touch {summary_output}
        """
        targets['summary'].append(target)

    return gwf, targets

# we need to assign the workflow to the gwf variable to allow the workflow to be
# run separetely with 'gwf run' in the submoduleA dir
gwf, targets = submoduleA_workflow(input_files=['./input.txt'], output_dir='A_outputs')
```

Thw workflow can be then be run the normal way:

```
gwf run
```

This way of writing workflows allows allows multiple submodule workflows to be combined in a master workflow file.

Try to put the `workflow.py` is in a `submoduleA` folder and the one below is in a `submoduleB` folder.

```python
import os.path
import os
from collections import defaultdict
from gwf import Workflow

def submoduleB_workflow(working_dir=os.getcwd(), input_files=None, output_dir=None, summarize=True):

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # dict of targets as info for other workflows
    targets = defaultdict(list)

    gwf = Workflow(working_dir=working_dir)

    work_output = os.path.join(output_dir, 'B_output1.txt')
    target = gwf.target(
        name='B_work',
        inputs=input_files,
        outputs=[work_output],
    ) << f"""
    touch {work_output}
    """
    targets['work'].append(target)

    if summarize:
        summary_output = os.path.join(output_dir, 'B_output2.txt')
        target = gwf.target(
            name='B_summary',
            inputs=[work_output],
            outputs=[summary_output]
        ) << f"""
        touch {summary_output}
        """
        targets['summary'].append(target)

    return gwf, targets

# we need to assign the workflow to the gwf variable to allow the workflow to be
# run separetely with 'gwf run' in the submoduleB dir
gwf, targets = submoduleB_workflow(input_files=['./input.txt'], output_dir='./B_outputs')
```

The two submodule workflows above can be combined into parts of a master workflow. The master `workflow.py` sits in the parent dir of `submoduleA` and `submoduleB`:

```
├── submoduleA
│   └── workflow.py
├── submoduleB
│   └── workflow.py
└── workflow.py
```

If you write it like this:

```python
import os
from gwf.workflow import collect

from submoduleA.workflow import submoduleA_workflow
from submoduleB.workflow import submoduleB_workflow

working_dir = os.getcwd()

def target_output_files(targets, key):
    return [out for target in targets[key] for out in target.outputs]

# submodule A workflow
gwf, A_targets  = submoduleA_workflow(working_dir=working_dir,
                          input_files=['./input.txt'],
                          output_dir='./A_outputs')
globals()['submoduleA'] = gwf

# add an extra target to glue workflows together
gwf.target('extra', inputs=['./input.txt'], outputs=['./A_outputs/extra.txt']) << """
    touch ./A_outputs/extra.txt
"""

# get relevant outputs from A for input to B
input_files =  target_output_files(A_targets, 'work')

# submodule B workflow
gwf, B_targets = submoduleB_workflow(working_dir=working_dir, 
                          input_files=input_files,
                          output_dir='./B_outputs' )
globals()['submoduleB'] = gwf
```

You can run each component workflow like this:

```
gwf -f workflow.py:submoduleA run
gwf -f workflow.py:submoduleB run
```