from clip_HighResMIP_utils import get_clip_extent
import os
import sys

strExtent = sys.argv[1]
aoiID = sys.argv[2]

# Region of interest
targetgrid = sys.argv[3]

src_rootdir = sys.argv[4]
dst_rootdir = sys.argv[5]

north, south, east, west = [float(x) for x in strExtent.split(",")]


#******************************************************************************
# Future climate
#******************************************************************************


# If the destination root directory does not exist, create it
if not os.path.isdir(dst_rootdir):
    os.makedirs(dst_rootdir)

for varcode in ["pr", "rlds", "rsds", "tas", "tasmax", "tasmin"]:

    # If the destination directory does not exist, create it
    if not os.path.isdir(dst_rootdir + "/" + varcode):
        os.makedirs(dst_rootdir + "/" + varcode)

    flag = 0

    for year in range(2015, 2051):
        filename = f"{varcode}_day_EC-Earth3P-HR_highres-future_r1i1p2f1_gr_{year}0101-{year}1231.nc"
        newfilename = filename.replace(".nc", f"_{aoiID}.nc")
        src = f"{src_rootdir}/{varcode}/{filename}"
        dst = f"{dst_rootdir}/{varcode}/{newfilename}"

        # Obtain the extent from the first file only (the rest should be the same)
        # TODO: provide functionality to deal with a situation when it's not the case.
        if flag != 1:
            latmin, latmax, lonmin, lonmax = get_clip_extent(targetgrid, src)
            flag=1

        print(f"Input = {src}")   
        nco_cmd = f"ncks --dmn lat,{latmin},{latmax} --dmn lon,{lonmin},{lonmax} {src} {dst}"
        os.system(nco_cmd)


#******************************************************************************
# Historical climate
#******************************************************************************


# If the destination root directory does not exist, create it
if not os.path.isdir(dst_rootdir):
    os.makedirs(dst_rootdir)

for varcode in ["pr", "rlds", "rsds", "tas", "tasmax", "tasmin"]:

    # If the destination directory does not exist, create it
    if not os.path.isdir(dst_rootdir + "/" + varcode):
        os.makedirs(dst_rootdir + "/" + varcode)

    flag = 0

    for year in range(2015, 2051):
        filename = f"{varcode}_day_EC-Earth3P-HR_highres-future_r1i1p2f1_gr_{year}0101-{year}1231.nc"
        newfilename = filename.replace(".nc", f"_{aoiID}.nc")
        src = f"{src_rootdir}/{varcode}/{filename}"
        dst = f"{dst_rootdir}/{varcode}/{newfilename}"

        # Obtain the extent from the first file only (the rest should be the same)
        # TODO: provide functionality to deal with a situation when it's not the case.
        if flag != 1:
            latmin, latmax, lonmin, lonmax = get_clip_extent(targetgrid, src)
            flag=1

        print(f"Input = {src}")   
        nco_cmd = f"ncks --dmn lat,{latmin},{latmax} --dmn lon,{lonmin},{lonmax} {src} {dst}"
        os.system(nco_cmd)