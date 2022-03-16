


def get_item_abspath(dataset, relpath):
    for identifier in dataset.identifiers:
        if dataset._manifest['items'][identifier]["relpath"] == relpath:
            return dataset.item_content_abspath(identifier)
    return False