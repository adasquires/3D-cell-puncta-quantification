import numpy as np
import pandas as pd
import scipy as sp

class quantify:
    def __init__(self,
                 centroids_puncta1,
                 area_puncta1,
                 centroids_puncta2,
                 area_puncta2,
                 centroids_cells,
                 area_cells,
                 threshold,
                 puncta1_min = 1000,
                 puncta1_max = 50000,
                 cell_min = 200000,
                 cell_max = 5000000,
                 dist = None,
                 spacing = None):

        self.centroids_puncta1 = centroids_puncta1
        self.area_puncta1 = area_puncta1
        self.centroids_puncta2 = centroids_puncta2
        self.area_puncta2 = area_puncta2
        self.centroids_cells = centroids_cells
        self.area_cells = area_cells
        self.threshold = threshold
        self.puncta1_min = puncta1_min
        self.puncta1_max = puncta1_max
        self.cells_min = cell_min
        self.cells_max = cell_max
        self.dist = dist
        self.spacing = spacing


    def __repr__(self):
        _, average_count = quantify.count(self)
        _, average_volume = quantify.volume(self)
        average_coloc = quantify.coloc(self)
        return(f'''{len(average_count)} cells
average LD count per cell: {np.mean(average_count)}
average LD volume per cell: {np.mean(average_volume)}
colocalization: {np.mean(average_coloc)}''')

    def filter(self, centroids, area, min, max):
        # remove smallest and largest features

        filtered_area = area[area>min].replace("nan", np.nan).dropna(ignore_index=False)
        filtered_area = filtered_area[filtered_area<max].replace("nan", np.nan).dropna(ignore_index=False)

        filtered_centroids = centroids.loc[filtered_area.index]

        return filtered_centroids, filtered_area

    def count(self):

        # filter smallest and largest features.
        filtered_centroids_puncta1, filtered_area_puncta1 = quantify.filter(self,
                                                                  self.centroids_puncta1,
                                                                  self.area_puncta1,
                                                                  self.puncta1_min,
                                                                  self.puncta1_max)

        filtered_centroids_cells, filtered_area_cells = quantify.filter(self,
                                                                    self.centroids_cells,
                                                                    self.area_cells,
                                                                    self.cells_min,
                                                                    self.cells_max)

        # compute distance matrix between centroids.
        dist_matrix = pd.DataFrame(sp.spatial.distance_matrix(filtered_centroids_cells,
                                                              filtered_centroids_puncta1))

        # get values from distance matrix below dist.
        neighbors = np.array(dist_matrix < self.dist)

        # create dictionary: keys are cell IDs, values are lists of puncta1 IDs.
        count = {}

        # add puncta1's within specified distance of cell to value associated with cell key.
        for i in range(neighbors.shape[0]):
            count[i] = []
            id = []
            for j in range(neighbors.shape[1]):
                if neighbors[i, j]:
                    id.append(j)
                count[i] = id

        # create list of number of puncta1's per cell.
        average_count = []

        for i in range(len(count)):
            average_count.append(len(count[i]))

        average_count = np.array(average_count)[~np.isnan(np.array(average_count))]

        # count : dictionary. keys = cell ID, values = puncta1 IDs.
        # average_count : list of average puncta1 count per cell.
        return count, average_count

    def volume(self):

        # filter largest and smallest features
        filtered_centroids_puncta1, filtered_area_puncta1 = quantify.filter(self,
                                                                  self.centroids_puncta1,
                                                                  self.area_puncta1,
                                                                  self.puncta1_min,
                                                                  self.puncta1_max)

        filtered_centroids_cells, filtered_area_cells = quantify.filter(self,
                                                                    self.centroids_cells,
                                                                    self.area_cells,
                                                                    self.cells_min,
                                                                    self.cells_max)

        # scale area size from pixels to mcirons
        scale = self.spacing[0] * self.spacing[1] * self.spacing[2]
        filtered_area_puncta1 = filtered_area_puncta1 * scale

        # get puncta1's sorted by cell
        count, _ = quantify.count(self)

        volume = {}

        # add area of puncta1's to value associated with cell key
        for i in count.items():
            cell = []
            for j in i[1]:
                area = filtered_area_puncta1.iloc[j][0]
                cell.append(area)
            volume[i[0]] = cell

        # get list of average puncta1 volume per cell
        average_volume = []
        for i in volume.values():
            average_volume.append(np.mean(i))

        average_volume = np.array(average_volume)[~np.isnan(np.array(average_volume))]

        # volume : dictionary. keys = cell ID. values = puncta1 volumes.
        # average_volume : list of average puncta1 volume per cell.
        return volume, average_volume

    def coloc(self):

        # filter smallest and largest features.
        filtered_centroids_puncta1, filtered_area_puncta1 = quantify.filter(self,
                                                                  self.centroids_puncta1,
                                                                  self.area_puncta1,
                                                                  self.puncta1_min,
                                                                  self.puncta1_max)

        filtered_centroids_cells, filtered_area_cells = quantify.filter(self,
                                                                    self.centroids_cells,
                                                                    self.area_cells,
                                                                    self.cells_min,
                                                                    self.cells_max)

        filtered_centroids_puncta2, filtered_area_puncta2 = quantify.filter(self,
                                                                  self.centroids_puncta2,
                                                                  self.area_puncta2,
                                                                  self.puncta1_min,
                                                                  self.puncta1_max)

        # compute distance matrix between cells and puncta1 centroids.
        cells_dist_matrix = pd.DataFrame(sp.spatial.distance_matrix(filtered_centroids_puncta1,
                                                                  filtered_centroids_cells))

        # get values from distance matrix below dist.
        cells_neighbors = np.array(cells_dist_matrix < self.dist)

        # compute distance matrix between LD and lipid dye centroids.
        puncta2_dist_matrix = pd.DataFrame(sp.spatial.distance_matrix(filtered_centroids_puncta1,
                                                                 filtered_centroids_puncta2))

        puncta2_neighbors = np.array(puncta2_dist_matrix < self.dist)

        # get values from distance matrix below threshold.
        dye_neighbors = np.array(puncta2_dist_matrix < self.threshold)

        # create dictionaries: keys are cell IDs, values are puncta1/puncta2 IDs associated with each cell.
        puncta2_count_per_cell = {}
        puncta1_count_per_cell = {}

        for i in range(cells_neighbors.shape[1]):
            puncta2_count_per_cell[i] = []
            puncta1_count_per_cell[i] = []

        # create dictionary of puncta1 count per cell.
        for i in range(cells_neighbors.shape[1]):
            puncta1_count = []
            for j in range(cells_neighbors.shape[0]):
                if cells_neighbors[j, i]:
                    puncta1_count.append(j)
            puncta1_count_per_cell[i] = puncta1_count

        # create dictionary of overlapped puncta2 count per cell.
        for i in puncta1_count_per_cell.items():
            puncta2_count = []
            for j in i[1]:
                for k in range(puncta2_neighbors.shape[1]):
                    if puncta2_neighbors[j, k]:
                        puncta2_count.append(k)
                        break
            puncta2_count_per_cell[i[0]] = puncta2_count

        # divide the number of overlapped dye by the number of puncta1's per cell.
        average_coloc = []
        for i in range(len(puncta2_count_per_cell)):
            try:
                average_coloc.append(len(puncta2_count_per_cell[i]) / len(puncta2_count_per_cell[i]))
            except:
                continue

        # average_coloc : list of colocalization (# overlapped with x / total # x) per cell.
        return average_coloc
