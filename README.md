# Overview
A tool for quantifying 3-dimensional fluorescent images.

Input: three channels...

> 0: cell channel
> 
> 1: puncta1 channel
> 
> 2: puncta2 channel

Output: CSV file, including...

> - average volume of puncta1 per cell.
> - average count of puncta1 per cell.
> - average colocalization (ratio, 1.0=perfect colocalization) of puncta1 and puncta2 per cell.

# Running
In terminal, type `python 3D-cell-puncta-quantification/main.py`.

After running, input `filename`.

# Parameters
in main.py
 - `spacing` for microscope spacing `[z, y, x]`

in segmentation.py
 - `sigma` for Gaussian filter
 - `width` for morphological erosion/dilation

in quantification.py
 - `threshold` for colocalization of puncta
 - `dist` for identifying puncta within cells
 - `puncta1_min` for filtering puncta
 - `puncta1_max` for filtering puncta
 - `cells_min` for filtering cells
 - `cells_max` for filtering cells
