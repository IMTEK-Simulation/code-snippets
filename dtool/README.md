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

## Select datasets with a query and extract some information into a pandas dataframe: 

```python
async def collect_data_from_query(query, verbose=False):
    res = await dl.query(query)
    ds = []
    for uri in pd.DataFrame(res).uri.to_list():
        try: 
            yaml = YAML()
            print(f"processing {uri}")
            dataset = dtoolcore.DataSet.from_uri(uri.replace("smb://as1412", "smb://isilon"))
            yml = yaml.load(dataset.get_readme_content())
            # see https://dtoolcore.readthedocs.io/en/latest/descriptive.html#working-with-items-in-a-dataset
            d = dict(yml["parameters"])
            d.update(yml["results"])
            d.update(uuid=dataset.uuid,
                     uri=dataset.uri)
            ds.append(d)
        except Exception as err: 
            print(f"uri: {uri}, error: {err}")
    df = pd.DataFrame(ds)
    df.set_index("uuid")
    return df
```

## Manage datasets in python

Here is an example where I wanted to rerun a bunch of short simulations with slightly different parameters. I hence automatized the logistics using dtool. 

```python
import subprocess

from dtoolcore import DataSet
import dtoolcore

from dtool_create.dataset import _get_readme_template
from ruamel.yaml import YAML
from datetime import date
from simulation_base import directory_name
from dtoolcore import create_proto_dataset
import io
import os
from dtoolcore.storagebroker import _get_abspath_from_uri

uris="""ecs://simdata/lasdjfhadjksfgjadhsfas
ecs://simdata/kladsfhdfskhdfsdfsadfsh"""

uris = uris.split("\n")

def get_item_abspath(dataset, relpath):
    if isinstance(dataset, str):
        dataset =  DataSet.from_uri(f"ecs://simdata/{dataset}")
        
    for identifier, properties in dataset.generate_manifest()["items"].items():
        if properties["relpath"] == relpath:
            return dataset.item_content_abspath(identifier)
    return False

# Open the readme template. if fpath is not specified it will find the standard (or last) readme template
readme_template = _get_readme_template(fpath=os.path.dirname(__file__) + "/dtool_readme.yml")
yaml = YAML()
yaml.explicit_start = True
yaml.indent(mapping=2, sequence=4, offset=2)

for uri in uris:

    # load the parameters from an existing dataset
    print(f"base dataset: {uri}")
    descriptive_metadata = yaml.load(readme_template)

    dataset = DataSet.from_uri(uri)

    parameters = yaml.load(dataset.get_readme_content())["parameters"]
    
    # modify them slightly
    parameters.update(dict(param1=6, param2=12, dump_fields=False))
    descriptive_metadata.update(dict(parameters=parameters))
    readme_stream = io.StringIO()
    yaml.dump(descriptive_metadata, readme_stream)
    readme = readme_stream.getvalue()
    
    # create dataset. This creates the directory
    dataset = create_proto_dataset(
        directory_name(parameters),
        ".",
        readme_content=readme
    )

    # Do the simulation 
    
    abspath = _get_abspath_from_uri(dataset.uri)
    errcode = subprocess.call(
        f"python3 simulation_base.py --simulate_dataset {abspath} "
        f"> {abspath}/data/stdout 2> {abspath}/data/stderr", shell=True)

    assert errcode == 0
    
    # Freeze and upload the simulation
    dataset.freeze()
    
    new_uri = dtoolcore.copy(dataset.uri, "ecs://simdata/")
    print(f"pushed: {new_uri}")
    
    # Delete local copy
    if True:
        errcode = subprocess.call(f"rm -rf {abspath}", shell=True)
        if errcode == 0:
            print("deletion succeeded")
        else:
            print(f"could not delete {abspath}")
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

## Compare local dataset repository against registered entries on lookup server

See README and script in subfolder `dtool_lookup_compare`. 

## Integrating dtool in FireWorks workflows

See https://github.com/IMTEK-Simulation/imteksimfw/blob/master/imteksimfw/fireworks/user_objects/firetasks/dtool_tasks.py.

Integration merges workflow metadata into dtool's README.yml and uses same scheme as above for derived datasets, i.e. listing UUIDs of source datasets under field 'derived_from'.
