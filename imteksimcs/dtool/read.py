


def list_items(dataset, pattern=None):
    r'''
    Parameters:
    -----------
    dataset: dtoolcore.DataSet
    pattern: string
    show only items with relpath matching the pattern
    '''
    item_relpaths = []
    for key, item in dataset._manifest['items'].items():
        relpath = dataset._manifest['items'][key]['relpath']
        if pattern and not pattern in relpath:
            continue
        item_relpaths.append(relpath)
    return item_relpaths

def get_item_abspath(dataset, relpath):
    for identifier in dataset.identifiers:
        if dataset._manifest['items'][identifier]["relpath"] == relpath:
            return dataset.item_content_abspath(identifier)
    return False
