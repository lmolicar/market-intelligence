import numpy as np
from osgeo import gdal, osr
import sys


year = 1981
country = "Ethiopia"
countryupp = country.upper()

srcdir = "//CATALOGUE.CGIARAD.ORG/OneSustainableFinance/ANALYZER/CLIMATE/CHIRPS-2.0/global_daily/netcdf/p25"
dstdir = f"//CATALOGUE.CGIARAD.ORG/OneSustainableFinance/ANALYZER/CLIMATE/CHIRPS-2.0/{countryupp}_daily"

# in general, this will be a GeoTIFF file (the ROI defined for NDWL0)
target_grid = sys.argv[1]

# This is the grid that I want to clip
reference_grid = sys.argv[2]

# Debugging
#print(f"Arg1 = {target_grid}")
#print(f"Arg2 = {reference_grid}")
#exit(0)

# I need the extent of the reference grid and obtain the coordinates of grid
ds = gdal.Open(reference_grid)

ref_xulcorner, xres, param1, ref_yulcorner, param2, yres = ds.GetGeoTransform()
ref_nlats = ds.RasterYSize
ref_nlons = ds.RasterXSize

ds = None
del ds

ref_lats = [ref_yulcorner + k*yres for k in range(ref_nlats)]
ref_lats = ref_lats[::-1]
ref_lons = [ref_xulcorner + k*xres for k in range(ref_nlons)]


# I need the the extent of the target grid and the coordinates of the grid
fnc = gdal.Open(target_grid)

tgt_xulcorner, xres, param1, tgt_yulcorner, param2, yres = fnc.GetGeoTransform()
tgt_nlats = fnc.RasterYSize
tgt_nlons = fnc.RasterXSize

fnc = None
del fnc


# Find the cell in the reference grid that contains the upper left corner of the target grid
ulxpos = np.searchsorted(ref_lons, tgt_xulcorner) - 1
#ulypos = ref_nlats - np.searchsorted(ref_lats, tgt_yulcorner) - 1
ulypos = 400 - np.searchsorted(ref_lats, tgt_yulcorner) - 1


# Find the cell in the reference grid that contains the lower right corner of the target grid
tgt_xlrcorner = tgt_xulcorner + xres * tgt_nlons
tgt_ylrcorner = tgt_yulcorner + yres * tgt_nlats
lrxpos = np.searchsorted(ref_lons, tgt_xlrcorner) - 1
#ulypos = ref_nlats - np.searchsorted(ref_lats, tgt_yulcorner) - 1
lrypos = 400 - np.searchsorted(ref_lats, tgt_ylrcorner) - 1


# Increase the extent by one cell
ulxpos -= 1
ulypos -= 1
lrxpos += 1
lrypos += 1

print(f"{ulypos} {lrypos} {ulxpos} {lrxpos}")