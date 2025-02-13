import napari
import numpy as np

def plot(img, cells, puncta1, puncta2):
    # plot 3D images and segmentation in napari.

    viewer = napari.Viewer()

    labels_as_image_layer = viewer.add_image(
        img[:, :, :, 0],
        name = 'Blue'
    )
    labels_as_image_layer = viewer.add_image(
        img[:, :, :, 1],
        name = 'Green'
    )
    labels_as_image_layer = viewer.add_image(
        img[:, :, :, 2],
        name = 'Red'
    )

    cells = np.array(cells, dtype=int)
    labels_as_image_layer = viewer.add_image(
        cells,
        name = 'Cells'
        )

    ld = np.array(puncta1, dtype=int)
    labels_as_image_layer = viewer.add_image(
        ld,
        name = 'puncta1'
        )

    dye = np.array(puncta2, dtype=int)
    labels_as_image_layer = viewer.add_image(
        dye,
        name = 'puncta2'
        )

    napari.run()