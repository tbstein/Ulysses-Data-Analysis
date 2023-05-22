import numpy as np
import spiceypy as spice

#https://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
def find_nearest_idx(array: list, value: float) -> int:
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx