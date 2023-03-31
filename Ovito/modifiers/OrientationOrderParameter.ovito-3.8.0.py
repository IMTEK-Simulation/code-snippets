from ovito.data import *

import numpy as np


def modify(frame: int, data: DataCollection):
    """Compute orientation order parameters using x, y, and z axes as directors."""

    if "Selection" in data.particles:
        selection = data.particles["Selection"]
    else:
        selection = np.ones(data.particles.count)
    print(f"Selection shape: {selection.shape}")

    index_selection = np.nonzero(selection)[0]
    print(f"Index selection shape: {index_selection.shape}")

    # Extract data from previously GyrationTensor modifiers
    tilt_cosine = data.particles['Tilt Cosine'][index_selection]
    print(f"tilt_cosine shape: {tilt_cosine.shape}")

    # orientation order parameter:
    # http://gisaxs.com/index.php/Orientation_order_parameter
    # S = (3*<cos^2 X>-1)/2

    # Orientation order parameter along x, y, and z axis
    data.attributes["Orientation Order Parameter"] = (3. * np.mean(np.square(tilt_cosine), axis=0) - 1.) / 2.