import geopandas as gpd
import time
import numpy as np
import pandas as pd
from progress.bar import Bar
import os
import sys

key = sys.argv[3]
lut_type = sys.argv[4]

# Two possible values: simple or classes
lut_type = "simple"

# Load target polygons layer
src_country_grid = sys.argv[1]
country_grid = gpd.read_file(src_country_grid)

# Load the protected areas layer
src_pa = sys.argv[2]
pa = gpd.read_file(src_pa)

# Intersect the two geometries
start = time.time()
print("\n\n")
print("************************************************************************")
print("                        Intersecting polygons")
print("************************************************************************")
cell_pa_intersec = pa.overlay(country_grid, how = "intersection")
end = time.time()
print("\nIntersection - Elapsed time: %.4f seconds" % (end - start))
print("\n\n")

# If the areas did not intersect
if len(cell_pa_intersec) == 0:
    merged = gpd.GeoDataFrame({key : country_grid[key],
                              "IOU" : 0,
                              "CATEGORY" : "No risk",
                              "geometry": country_grid.geometry})

# Otherwise
else:
    start = time.time()
    print("\n\n")
    print("************************************************************************")
    print("                        Computing IoU")
    print("************************************************************************")

    # Make a list of the unique IDs of the cells that did intersect the protected areas
    cell_did_intersec = np.unique(cell_pa_intersec[key])

    intersec_cell_mapid = []
    intersec_cell_iou = []

    numItems = len(cell_did_intersec)

    bar = Bar("\n\nProcessing", max = numItems)

    for each_mapid in cell_did_intersec:

        intersec_cell_mapid.append(each_mapid)

        cell_intersections = cell_pa_intersec[cell_pa_intersec[key] == each_mapid]
        cell_intersections.reset_index(inplace=True)

        # Obtain a list of all the protected areas that intersected
        unique_pa = np.unique(cell_intersections[['WDPA_PID']])

        # Iterate through each PA intersected
        numPAs = len(unique_pa)
        iou = np.zeros(numPAs)

        # Intersection area
        area_intersec = gpd.GeoSeries(cell_pa_intersec[cell_pa_intersec[key] == each_mapid].geometry).area.values[0]

        for index_pa in range(numPAs):

            curpa = unique_pa[index_pa]

            # ... and carry out union operation
            curpa_df = pa[pa['WDPA_PID'] == curpa]
            curcell_df = country_grid[country_grid[key] == each_mapid]
            pa_cell_union = curpa_df.overlay(curcell_df, how = 'union')
            dissolved = pa_cell_union.dissolve()
            area_union = dissolved.area.values[0]

            iou[index_pa] = area_intersec / area_union

        final_iou = np.average(iou)
        intersec_cell_iou.append(final_iou)
        bar.next()

    end = time.time()
    print("\nElapsed time: %.4f seconds" % (end - start))
    print("\n\n")

    bar.finish()

    output_df = pd.DataFrame({key : intersec_cell_mapid, "IOU" : intersec_cell_iou})
    cell_id = country_grid[key]
    #merged = pd.merge(country_grid[[key, "geometry"]],
    merged = pd.merge(country_grid,
                      output_df,
                      on = key,
                      how = 'left')
    
    true_values = merged['IOU'].values[~np.isnan(merged['IOU'].values)]

    quantiles = np.quantile(true_values, np.linspace(0.25, 1, 4))

    if lut_type == "classes":
        def lut(x, p = None):
            if np.isnan(x):
                return "No risk"

            if x <= p[0]:
                return "No significant stress"
            elif x > p[0] and x <= p[1]:
                return "Moderate stress"
            elif x > p[1] and x <= p[2]:
                return "Severe stress"
            elif x > p[2] and x <= p[3]:
                return "Extreme stress"
            
    elif lut_type == "simple":
        def lut(x, p = None):
            if np.isnan(x):
                return "No risk"

            if 0 < x:
                return "Potential risk"
            else:
                return "No risk"


    vect_lut = np.vectorize(lut, otypes = [str], excluded = "p")

    if lut_type == "simple":
        category = vect_lut(merged['IOU'])
    elif lut_type == "classes":
        category = vect_lut(merged['IOU'], p = quantiles)

    merged['CATEGORY'] = category
    
    # I need to convert all NaN values in IOU to zero
    merged['IOU'] = np.where(np.isnan(merged['IOU'].values), 0, merged['IOU'].values)

outputDir = sys.argv[5]

# add a command line parameter to specify whether to save the results
# as plain JSON file or include the geospatial information (GeoJSON)
geom = sys.argv[6]

# If geom == True then save as GeoJSON
if geom == "True":
    filename = os.path.basename(src_country_grid).split(".")[0]
    filename = filename.replace("_Mollweide", "_biodiversity_Mollweide")
    outputfile = f"{outputDir}/{filename}"
    merged.to_file(outputfile + '.geojson', driver="GeoJSON")

# If geom == False then save as plain JSON
elif geom == "False":
    filename = "biodiversity.json"
    outputfile = f"{outputDir}/{filename}"
    output_simple = merged[["Map_ID", "IOU", "CATEGORY"]]
    output_simple = output_simple.rename(columns = {"IOU" : "natural_capital_iou", "CATEGORY" : "natural_capital_category"})
    output_simple.to_json(outputfile, orient = "records", indent = 4)