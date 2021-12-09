
from ruamel.yaml import YAML


def set_derived_from(uuids, readme_file="../README.yml", mode="w"):
    """

    Parameters
    ----------
    uuids: list or set of str
        uuids to be added to derived_from
    readme_file: str
        path ro readme file, default `../README.yml`
    mode: {"w", "a"}
        `w`: overwrites `derived_from` field
        `a`: appends uuids to the list existing in the file

    Returns
    -------

    """
    # Configure yaml I/O
    yaml = YAML()
    yaml.explicit_start = True
    yaml.width = 80
    yaml.indent(mapping=2, sequence=4, offset=2)

    # open df
    with open(readme_file, encoding='utf8') as f:
        yml = yaml.load(f)
    if yml is None:  # empty readme file
        yml = dict()
    if mode == "w":
        derived_from = []
    elif mode == "a":
        if "derived_from" in yml.keys():
            derived_from = yml["derived_from"]
        else:
            derived_from = []
    else:
        raise ValueError(f"mode {mode} not recognized")

    yml.update(
        dict(
            derived_from=derived_from + [dict(uuid=uuid) for uuid in uuids])
        )

    # Rewrite readmefile
    with open(readme_file, "w") as f:
        yaml.dump(yml, f)
