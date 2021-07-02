


def get_item_abspath(dataset, relpath):
    for identifier, properties in dataset.generate_manifest()["items"].items():
        if properties["relpath"] == relpath:
            return dataset.item_content_abspath(identifier)
    return False