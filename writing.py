import csv
import pandas as pd

def write_to_csv(filename, count, volume, coloc):
    # write CSV file.

    _, average_count = count        # average puncta count per cell
    _, average_volume = volume      # average puncta volume per cell (µm^3)
    average_coloc = coloc           # average colocalization per cell (overlapped count / total puncta count)
    average_count = [i for i in average_count if i != 0]

    values = pd.DataFrame()
    values['average_count'] = pd.Series(average_count)
    values['average_volume'] = pd.Series(average_volume)
    values['average_coloc'] = pd.Series(average_coloc)

    values.to_csv(str(filename) + ".csv")
