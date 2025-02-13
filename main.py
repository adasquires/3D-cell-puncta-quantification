import time
from segmentation import segment
from plotting import plot
from quantification import quantify
from writing import write_to_csv

import skimage as skimage
from skimage import filters, morphology, measure, segmentation

start = time.time()

filename = input('Enter filename: ')
z, x, y = 0.5, 0.0625, 0.0625
spacing = [z, x, y]     # microscope spacing


if __name__ == "__main__":
    print(f'segmenting...')
    transformed_cell, centroids_cell, area_cell, segmented_cell, speckles_cell, labels_cell = segment(filename,
                                                                 spacing = spacing,
                                                                 channel = 0,
                                                                 sigma = 20,
                                                                 width = 20)

    transformed_puncta1, centroids_puncta1, area_puncta1, segmented_puncta1, speckles_puncta1, labels_puncta1 = segment(filename,
                                                                 channel = 1,
                                                                 spacing = spacing,
                                                                 sigma = 1,
                                                                 width = 5)

    transformed_puncta2, centroids_puncta2, area_puncta2, segmented_puncta2, speckles_puncta2, labels_puncta2 = segment(filename,
                                                                                         channel = 2,
                                                                                         spacing = spacing,
                                                                                         sigma = 2,
                                                                                         width = 3)

    image = skimage.io.imread(filename)
    print(f'segmenting...done ' + str(time.time() - start))

    print(f'quantifying...')
    quantified = quantify(centroids_puncta1 = centroids_puncta1,
                  area_puncta1 = area_puncta1,
                  centroids_cells = centroids_cell,
                  area_cells = area_cell,
                  centroids_puncta2 = centroids_puncta2,
                  area_puncta2 = area_puncta2,
                  threshold = 15,
                  dist = 180,
                  puncta1_min = 90,
                  puncta1_max = 5000,
                  cell_min = 25000,
                  cell_max = 60000,
                  spacing = spacing)
    print(f'quantifying...done ' + str(time.time() - start))

    print(f'writing to csv...')
    write_to_csv(filename[:-4], quantified.count(), quantified.volume(), quantified.coloc())
    print(f'writing to csv...done ' + str(time.time() - start))

    print(f'plotting...')
    plot(image, segmented_cell, segmented_puncta1, segmented_puncta2)
    print(f'plotting...done ' + str(time.time() - start))