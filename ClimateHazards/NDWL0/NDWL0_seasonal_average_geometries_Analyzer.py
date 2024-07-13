import numpy as np                  # For math and matrix operations
from osgeo import gdal              # For raster manipulation 
import glob                         # List directory using patterns
import datetime                     # For date and time objects
import pandas as pd                 # For data manipulation
from rasterstats import zonal_stats # For zonal statistics on rasters
import geopandas as gpd             # For geospatial data manipulation
import rasterio                     # For reaster manipulation
import argparse                     # To implement a command-line interface
import os
import calendar


parser = argparse.ArgumentParser(
    prog = "NDWL0_seasonal_average_geometries_Analyzer.py",
    description = "Compute seasonal average of NDWL",
    epilog = "\n\n"
)


parser.add_argument("--params-file", help = "Full path the parameters file.", required = True)
parser.add_argument("--output-dir", help = "Full path to the file where final data will be stored. Data will be stored in GEOJSON format.", required = True)

args = parser.parse_args(["--params-file", "\\\\CATALOGUE.CGIARAD.ORG\\OneSustainableFinance\\ANALYZER\\NDWL0\\SCRIPTS\\crop_conditions_time.csv",
                          "--output-dir", "\\\\CATALOGUE.CGIARAD.ORG\\OneSustainableFinance\\ANALYZER\\NDWL0\\DUMMY_ANALYZER"])

run_params = pd.read_csv(args.params_file)

for recnum in range(len(run_params)):
    print()
    print("*" * 80)
    print(f"Country: {run_params.COUNTRY[recnum]}, crop: {run_params.CROP_NAME[recnum]}")
    print("*" * 80)
    print()

    for curyear in range(run_params.YEAR_START[recnum], run_params.YEAR_END[recnum] + 1):
        if np.isnan(run_params.SEASON_END[recnum]):
            
            # This means the season lasts one month
            
            season_begin = f"{calendar.month_abbr[run_params.SEASON_START[recnum]]}-{curyear}"
            season_end   = f"{calendar.month_abbr[run_params.SEASON_START[recnum]]}-{curyear}"
        else:
            if run_params.SEASON_START[recnum] < run_params.SEASON_END[recnum]:
                
                # This means the season lasts more than one month; it starts and ends the same year
                
                season_begin = f"{calendar.month_abbr[run_params.SEASON_START[recnum]]}-{curyear}"
                season_end   = f"{calendar.month_abbr[run_params.SEASON_END[recnum]]}-{curyear}"
            elif run_params.SEASON_START[recnum] >= run_params.SEASON_END[recnum]:
                
                # This means the season starts one year and ends the next
                
                if curyear < run_params.YEAR_END[recnum]:
                    season_begin = f"{calendar.month_abbr[run_params.SEASON_START[recnum]]}-{curyear}"
                    season_end   = f"{calendar.month_abbr[run_params.SEASON_END[recnum]]}-{curyear + 1}"
                else:
                    continue

        # Debugging
        #print(f"Current year {curyear} {run_params.COUNTRY[recnum]} {run_params.CROP_NAME[recnum]} season {season_begin}_{season_end}")

        dir_time = f"{season_begin}_{season_end}"
        newpath = f"{args.output_dir}/{run_params.COUNTRY[recnum]}/{run_params.CROP_NAME[recnum]}/{dir_time}"
        if not os.path.exists(newpath):
            os.makedirs(newpath, exist_ok = True)

        srcdir = run_params.DATA[recnum]

        iyear = run_params.YEAR_START[recnum]
        fyear = run_params.YEAR_END[recnum]

        imonth = run_params.SEASON_START[recnum]
        fmonth = run_params.SEASON_END[recnum]

        lista = glob.glob(f"{srcdir}/NDWL0*.tif")

        geometries_gjson = run_params.GEOMETRY[recnum]
        ndwl0 = lista[0]

        dates = np.array([ datetime.datetime.strptime("%4.4d-%2.2d-%2.2d" % (year, month, 1), "%Y-%m-%d")
                for year in range(iyear, fyear + 1) for month in range(1,13) ])


        outvrt = '/vsimem/stacked.vrt' #/vsimem is special in-memory virtual "directory"

        df = pd.DataFrame()

        print(f" Year {curyear}")
        if imonth < fmonth:
            date_initial = datetime.datetime(curyear, imonth, 1)
            date_final = datetime.datetime(curyear, fmonth, 1)  
        elif imonth > fmonth:
            if curyear < fyear:
                date_initial = datetime.datetime(curyear, imonth, 1)
                date_final = datetime.datetime(curyear + 1, fmonth, 1)
            else:
                continue

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

        df["NDWL"] = result_series
        #df[f"season_{date_initial.year}_{date_final.year}"] = pd.DataFrame({f"season_{date_initial.year}_{date_final.year}" : result_series})

        # Reclassify the values
        def reclass_risk(x):
            retval = None
            
            if x < 2:
                retval = "No significant stress"
            elif x >= 2 and x < 5:
                retval = "Moderate stress"
            elif x >= 5 and x < 8:
                retval = "Severe stress"
            elif x > 8:
                retval = "Extreme stress"

            return retval
        
        vreclass_risk = np.vectorize(reclass_risk)

        risk_class = vreclass_risk(result_series.to_numpy())

        # Let's load the vector file
        geom = gpd.read_file(geometries_gjson)

        out_df = pd.DataFrame({"Map_ID" : geom.Map_ID, 
                               "NDWL" : np.squeeze(df.to_numpy()),
                               "NDWL_category" : risk_class})

        out_df.to_json(f"{newpath}/NDWL.json", mode = "w", orient = "records")
