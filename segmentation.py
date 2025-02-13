import numpy as np
import scipy as sp
from scipy import ndimage as ndi
import pandas as pd
import skimage
from skimage import filters, morphology, measure, segmentation
import time

class img:
    # thresholding, morphological operations, segmentation.

    def __init__(self, file, channel, sigma, width, spacing, pixels):
        self.img = skimage.io.imread(file)
        self.idx = channel
        self.sigma = sigma
        self.width = width
        self.original_spacing = np.array([spacing[0]/pixels[0], (spacing[1]/pixels[1]), (spacing[2]/pixels[2])])

    def __repr__(self):

        data = skimage.util.img_as_float(self.img[:, :, :, self.idx])

        info = f'''shape: {data.shape}
dtype: {data.dtype}
range: ({data.min()}, {data.max()})
microscope spacing: {self.original_spacing}'''

        if self.idx == 0:
            return ("blue channel\n" + info)
        elif self.idx == 1:
            return ("green channel\n" + info)
        elif self.idx == 2:
            return ("red channel\n" + info)

    def threshold_cells(self):
        # thresholding and morphological operations for cells.

        spacing = self.original_spacing / self.original_spacing[2]  # scale spacing
        img = self.img[:, :, :, self.idx]
        gaussian = ndi.gaussian_filter(                             # Gaussian filter
            img,
            sigma=self.sigma)
        edges = filters.sobel(gaussian)                             # Sobel filter for edge detection
        thresholded = edges > filters.threshold_li(edges)           # binary thresholding (Li)
        fill_holes = ndi.binary_fill_holes(thresholded)             # fill in cells
        eroded = morphology.isotropic_erosion(                      # binary morphological erosion
            fill_holes,
            radius = self.width/2)
        dilated = morphology.isotropic_dilation(                    # binary morphological dilation
            eroded,
            radius = self.width/2)

        print(f'distance transformation...')
        labels = measure.label(                                     # labeling
            dilated,
            connectivity = 3)
        transformed = ndi.distance_transform_edt(                   # distance transformation
            labels,
            sampling = spacing)

        return transformed, labels

    def threshold_puncta(self):
        # thresholding and morphological operations for puncta.

        print(f'thresholding...')
        spacing = self.original_spacing / self.original_spacing[2] # scale spacing
        img = self.img[:, :, :, self.idx]
        gaussian = ndi.gaussian_filter(                            # Gaussian filter
            img,
            sigma = self.sigma
        )
        denoised = ndi.median_filter(                              # median filter
            gaussian,
            size = 3
        )
        edges = filters.sobel(denoised)                         # Sobel filter for edge detection
        thresholded = edges > filters.threshold_li(edges)     # binary thresholding (Otsu)
        remove_holes = morphology.remove_small_holes(           # remove holes smaller than area_threshold (morphological closing)
            thresholded,
            area_threshold=self.width**3
        )
        remove_objects = morphology.remove_small_objects(       # remove objects smaller than area_threshold (morphological opening)
            remove_holes,
            min_size=self.width**3)

        print(f'distance transformation...')
        labels = measure.label(                                 # labeling
            remove_objects,
            connectivity=3
        )
        transformed = ndi.distance_transform_edt(               # distance transformation
            labels,
            sampling=spacing
        )   # distance transformation

        return transformed, labels

    def segment_cells(self, transformed, labels):
        # segmenting cells.

        coords = pd.DataFrame(measure.regionprops_table(      # get centroids of labels
                labels,
                properties=['centroid']
                ))
        mask = np.zeros(transformed.shape, dtype=bool)        # create markers for segmentation (array of basins)
        for i in np.round(np.array(coords, dtype=int)):
            j, k, l = i
            mask[j, k, l] = True
        markers, num_features = ndi.label(mask)
        segmented = segmentation.watershed(                   # segmentation
            transformed,
            markers=mask,
            mask=labels                                       # only points where mask==True will be labeled
        )
        speckles = measure.label(segmented)
        centroids = pd.DataFrame(measure.regionprops_table(   # get centroids of labels
                speckles,
                properties=['centroid']
                ))
        area = pd.DataFrame(measure.regionprops_table(        # get area of labels
                speckles,
                properties=['area']
                ))

        return centroids, area, segmented, speckles

    def segment_puncta(self, transformed, labels):
        # segmenting puncta.

        coords = pd.DataFrame(measure.regionprops_table(      # get centroids of labels
                labels,
                properties=['centroid']
                ))
        mask = np.zeros(transformed.shape, dtype=bool)        # create markers for segmentation
        for i in np.round(np.array(coords, dtype=int)):
            j, k, l = i
            mask[j, k, l] = True
        markers, num_features = ndi.label(mask)
        segmented = segmentation.watershed(                   # segmentation
            transformed,
            markers=markers,
            mask=labels                                       # only points where mask==True will be labeled
        )
        speckles = measure.label(segmented)
        centroids = pd.DataFrame(measure.regionprops_table(   # get centroids of labels
                speckles,

                properties=['centroid']
                ))
        area = pd.DataFrame(measure.regionprops_table(        # get area of labels
                speckles,
                properties=['area']
                ))

        return centroids, area, segmented, speckles

def segment(file,
         channel,
         spacing,                                 # size of image in microns (z-stack interval, x, y)
         sigma = 3,                               # standard deviation for Gaussian kernel
         width = 10,                              # threshold for morphological operations
         pixels = [1, 1, 1]):                     # pixel resolution of image (z, x, y)

    image = img(file,
                channel,
                sigma,
                width,
                spacing,
                pixels)

    print(image)

    if channel == 0:
        transformed, labels = image.threshold_cells()
    else:
        transformed, labels = image.threshold_puncta()

    if channel == 0:
        centroids, area, segmented, speckles = image.segment_cells(transformed, labels)
    else:
        centroids, area, segmented, speckles = image.segment_puncta(transformed, labels)

    return transformed, centroids, area, segmented, speckles, labels
