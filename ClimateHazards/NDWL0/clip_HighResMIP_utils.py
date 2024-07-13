import numpy as np
from osgeo import gdal
from netCDF4 import Dataset
import sys

def get_clip_extent(target_grid, reference_grid):

    # Debugging
    #print(f"Arg1 = {target_grid}")
    #print(f"Arg2 = {reference_grid}")
    #exit(0)

    # I need the extent of the reference grid and obtain the coordinates of grid.
    # Note that I intend to clip a NetCDF file. Only the NetCDF library can tell
    # the orientation of the grid.
    fnc = Dataset(reference_grid)

    ref_lats = fnc.variables['lat'][:]
    ref_lons = fnc.variables['lon'][:]

    fnc.close()

    # But keep this because it's an easy way to extract other useful parameters
    ds = gdal.Open(reference_grid)

    ref_xulcorner, xres, param1, ref_yulcorner, param2, yres = ds.GetGeoTransform()
    ref_nlats = ds.RasterYSize
    ref_nlons = ds.RasterXSize

    # From the NetCDF, coordinates are of the center of the pixel
    ref_lons -= xres/2
    ref_lats += yres/2

    ds = None
    del ds

    # I need to check the orientation
    northernmost = False
    print(ref_lats[0], ref_lats[-1])
    if ref_lats[0] > ref_lats[-1]:
        ref_lats = ref_lats[::-1] + yres
        northernmost = True

    # I need the the extent of the target grid and the coordinates of the grid
    fnc = gdal.Open(target_grid)

    tgt_xulcorner, xres, param1, tgt_yulcorner, param2, yres = fnc.GetGeoTransform()
    tgt_nlats = fnc.RasterYSize
    tgt_nlons = fnc.RasterXSize

    fnc = None
    del fnc

    tgt_xlrcorner = tgt_xulcorner + xres * tgt_nlons
    tgt_ylrcorner = tgt_yulcorner + yres * tgt_nlats

    # Find the cell in the reference grid
    ulxpos = np.searchsorted(ref_lons, tgt_xulcorner) - 1
    lrxpos = np.searchsorted(ref_lons, tgt_xlrcorner) - 1

    # Consider the orientation of the latitudes
    if northernmost:
        ulypos = (ref_nlats-1) - (np.searchsorted(ref_lats, tgt_ylrcorner) - 1) + 1
        lrypos = (ref_nlats-1) - (np.searchsorted(ref_lats, tgt_yulcorner) - 1) + 1

    else:
        ulypos = np.searchsorted(ref_lats, tgt_ylrcorner) - 1
        lrypos = np.searchsorted(ref_lats, tgt_yulcorner) - 1

    # Increase the extent by two cells
    ulxpos -= 2
    lrxpos += 2

    if northernmost:
        ulypos += 2
        lrypos -= 2
        if ulypos > lrypos:
            dummy = ulypos
            ulypos = lrypos
            lrypos = dummy
            del dummy
    else:
        ulypos -= 2
        lrypos += 2

    return (ulypos, lrypos, ulxpos, lrxpos)