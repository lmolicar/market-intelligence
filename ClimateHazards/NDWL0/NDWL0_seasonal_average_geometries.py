import numpy as np                  # For math and matrix operations
from osgeo import gdal              # For raster manipulation 
import glob                         # List directory using patterns
import datetime                     # For date and time objects
import pandas as pd                 # For data manipulation
from rasterstats import zonal_stats # For zonal statistics on rasters
import geopandas as gpd             # For geospatial data manipulation
import rasterio                     # For reaster manipulation
import argparse                     # To implement a command-line interface

parser = argparse.ArgumentParser(
    prog = "NDWL0_seasonal_average_geometries.py",
    description = "Compute seasonal average of NDWL",
    epilog = "\n\n"
)

parser.add_argument("--source-ndwl-dir", help = "Full path of the directory that contains the NDWL0 data.", required = True)
parser.add_argument("--years", help = "The years over which the average will be computed; specified as year1,year2, such that year1 <= year2.", required = True)
parser.add_argument("--months", help = "The season over which the average will be computed; specified as month1,month2, where month1 and month2 are the number of the month (i.e., 1, 2, ..., 12). When the season starts in one year and ends the next year (e.g., month1 = October and month2 = January) the program will look for data from the preceding year.", required = True)
parser.add_argument("--geometries", help = "Names of the files containing the features for which the average will be computed.", required = True)
parser.add_argument("--output", help = "Full path to the file where results will be stored. Data will be stored in GEOJSON format.", required = True)

args = parser.parse_args()

srcdir = args.source_ndwl_dir

iyear, fyear = [int(x) for x in args.years.split(",")]

imonth, fmonth = [int(x) for x in args.months.split(",")]

lista = glob.glob(f"{srcdir}/NDWL0*.tif")

geometries_gjson = args.geometries
ndwl0 = lista[0]

dates = np.array([ datetime.datetime.strptime("%4.4d-%2.2d-%2.2d" % (year, month, 1), "%Y-%m-%d")
         for year in range(iyear, fyear + 1) for month in range(1,13) ])


outvrt = '/vsimem/stacked.vrt' #/vsimem is special in-memory virtual "directory"

df = pd.DataFrame()

for year in range(iyear, fyear + 1):
    print(f"************** {year} ********************")
    date_initial = datetime.datetime(year-1 if imonth >= fmonth else year, imonth, 1)
    date_final = datetime.datetime(year, fmonth, 1)
    season_dates = dates[[x and y for x,y in zip(dates <= date_final, dates >= date_initial)]]

    # File names follow convention NDWL0-yyyy-mm.tif
    thefiles = [ f"{srcdir}/NDWL0-%4.4d-%2.2d.tif" % (curdate.year, curdate.month) for curdate in season_dates ]

    vrtsrc = gdal.BuildVRT(outvrt, thefiles, separate=True)
    data = vrtsrc.ReadAsArray()

    # Obtain affine transform
    with rasterio.open(thefiles[0]) as src:
        affine = src.transform

    # Get the nodata value
    nodata_value = vrtsrc.GetRasterBand(1).GetNoDataValue()
    if nodata_value != None:
        if np.isnan(nodata_value):
            data = np.ma.masked_where(np.isnan(data), data)
        else:
            data = np.ma.masked_where(data == nodata_value, data)

    season_average = np.ma.average(data, axis = 0)

    # Set all_touched to True. Otherwise, pixels smaller than the pixel may end up with empty values,
    # even when completely within the polygon but do not contain a pixel center.
    # Given the rasterization strategy of the zonal_stats function. this solution may, in some cases over-
    # or underestimate the value.

    # prior to estimating the seasonal average we need to 'tag' the pixels without data (fill values), i.e.,
    # assign an actual value to those pixels without data because zonal_stats doesn't know how to treat
    # masked arrays.
    season_average = np.ma.where(season_average.mask == True, -9999, season_average)
    result = zonal_stats(geometries_gjson, season_average, stats="mean",
                         affine = affine, all_touched = True, nodata = -9999)

    result_series = pd.Series([x['mean'] for x in result])

    if date_initial.year == date_final.year:
        seasonID = f"{date_initial.year}"
    else:
        seasonID = f"{date_initial.year}_{date_final.year}"

    df[f"NDWL0_{seasonID}"] = result_series
    #df[f"season_{date_initial.year}_{date_final.year}"] = pd.DataFrame({f"season_{date_initial.year}_{date_final.year}" : result_series})

# Let's load the vector file
geom = gpd.read_file(geometries_gjson)

# Add the columns to the geometry
column_names = df.columns.to_list()

for k in range(len(column_names)):
    geom[column_names[k]] = df[[column_names[k]]]


geom.to_file(args.output, driver='GeoJSON')