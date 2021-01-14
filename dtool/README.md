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

```python
#!/usr/bin/env python
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

## Accessing dataset readme from python.
You can access directly from isilon. So you only need to know th uuid.

```python
yaml = YAML()
dataset = DataSet.from_uri(f"smb://isilon/{uuid}")
yml = yaml.load(dataset.get_readme_content())
```

## Accessing dataset items: 


```python
def get_item_abspath(dataset, relpath):
    for identifier, properties in dataset.generate_manifest()["items"].items():
        if properties["relpath"] == relpath:
            return dataset.item_content_abspath(identifier)
    return False
``` 

## rapidly extract all uuids used in a file (assuming they are all hardcoded)

```python
import re

uuids = []
UUID_v4_REGEX = '[0-9a-fA-F]{8}-[0-9a-fA-F]{4}-[4][0-9a-fA-F]{3}-[89abAB][0-9a-fA-F]{3}-[0-9a-fA-F]{12}'
with open('my_postprocessing_file.py') as myfile:
    for line in myfile.readlines():
        result = re.search(UUID_v4_REGEX, line)
        if result is not None:
            print(line)
            uuids.append(result.group())

print('UUIDS: ', uuids)

out = "derived_from:\n"

for s in uuids:
    out+="  - uuid: " + s + "\n"

print(out)
```



## Integrating dtool in FireWorks workflows

See https://github.com/IMTEK-Simulation/imteksimfw/blob/master/imteksimfw/fireworks/user_objects/firetasks/dtool_tasks.py.

Integration merges workflow metadata into dtool's README.yml and uses same scheme as above for derived datasets, i.e. listing UUIDs of source datasets under field 'derived_from'.
