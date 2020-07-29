# Code snippets for managing dtool datasets

## modify a README.yml using python

Assuming the workdir is in `data`
```python
from ruamel.yaml import YAML

with open("../README.yml",) as f:
    yaml = YAML()
    yaml.explicit_start = True # ensures the file begins with "---"
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

## Create README from template

Creating README files interactively from the command line is a bit cumbersome. This script generates READMES from a template. 
Creation and expiration date are determined automatically from the name of the dataset. The name must start with the date 
in ISO-format, e.g. `dtool create $(date -I)_some_dataset; python create_dtool_readme.py $(date -I)_some_dataset`. The 
UUID from which the dataset is derived can be specified using the optional `-u` argument.

```#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Create dtool README files based on a template"""
from string import Template
import argparse
import pathlib
from datetime import date

basic_readme_template = Template("""\
---
project: $project
description: Short description
owners:
  - name: Firstname Lastname
    email: firstname.lastname@imtek.uni-freiburg.de
    username: fr_fnXXXX
    orcid: XXXX-XXXX-XXXX-XXXX
funders:
  - organization: European Research Council (ERC)
    program: program
    code: code
creation_date: $creation_date
expiration_date: $expiration_date
"""
)

uuid_template = Template("""\
derived_from:
  - uuid: $uuid
"""
)

minimum_storage_years = 10

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", type=pathlib.Path, help="directory")
    parser.add_argument("-u", "--uuid", help="uuid of the dataset")
    parser.add_argument("-p", "--print_only", help="Print readme only")
    args = parser.parse_args()
    iso_string = args.directory.name[:10]
    year, month, day = [int(x) for x in iso_string.split("-")]
    creation_date = date(year, month, day)
    expiration_date = date(year + minimum_storage_years, month, day)
    readme = basic_readme_template.substitute(creation_date=creation_date, expiration_date=expiration_date, project=args.directory.name)
    if args.uuid is not None:
        readme += uuid_template.substitute(uuid=args.uuid)
    print(readme)
    if not args.print_only:
        readme_location = pathlib.Path(args.directory, "README.yml")
        with open(readme_location, "w") as file:
            file.write(readme)
 ````
