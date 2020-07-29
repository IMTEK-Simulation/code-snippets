# Code snippets for managing dtool datasets

## modify a README.yml using python

Assuming the workdir is in `data`
```python
from ruamel.yaml import YAML

with open("../README.yml",) as f:
    yaml = YAML()
    yaml.explicit_start = True
    yaml.width = 80
    yaml.indent(mapping=2, sequence=4, offset=2)
    yml = yaml.load(f)

 yml.update(.......
........
........
    )
with open("../README.yml", "w") as f:
    yaml.dump(yml, f)

```

TODO: despite the `yaml.width` is set to `80`, it seems that the linebreaks do not pass the yamllinter


## Report the software version

At runtime of my simulation script, I update the README file with version informations on the software, including the repository containing the simulation script (assuming the simulation script is not modified on a dataset basis). 

```python

def get_repository_info():
    import Adhesion, CrackFront
    import os
    import subprocess

    repopath = os.path.dirname(os.path.abspath(__file__))
    return \
    [
    dict(
    name="this project",
    version=subprocess.check_output(
        ["git", "-C", repopath, "show", "-s", "--format=%h"]
        ).decode().replace("\n", ""),
    repository=subprocess.check_output(
        ["git", "-C", repopath, "config", "--get", "remote.origin.url"],
        ).decode().replace("\n", "")
    ),
    dict(name="CrackFront",
         version=CrackFront.__version__,
         repository="git@github.com/ComputationalMechanics/CrackFront.git"),
    dict(name="Adhesion",
             version=Adhesion.__version__,
             repository="git@github.com/ComputationalMechanics/Adhesion.git")
    
```