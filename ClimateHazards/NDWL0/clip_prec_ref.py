from utils import get_clip_extent
import os

country = "Ethiopia"
countryupp = country.upper()

# Region of interest
targetgrid="//CATALOGUE.CGIARAD.ORG/OneSustainableFinance/ANALYZER/CLIMATE/HighResMIP/ETHIOPIA/FUTURE/pr/pr_day_EC-Earth3P-HR_highres-future_r1i1p2f1_gr_20150101-20151231_rotated_Ethiopia.nc"

# Source directory
srcdir="//CATALOGUE.CGIARAD.ORG/OneSustainableFinance/ANALYZER/CLIMATE/CHIRPS-2.0/global_daily/netcdf/p25"

# Destination directory
dstdir = f"//CATALOGUE.CGIARAD.ORG/OneSustainableFinance/ANALYZER/CLIMATE/CHIRPS-2.0/{countryupp}_daily"

# If the destination directory does not exist, create it
if not os.path.isdir(dstdir):
    os.makedirs(dstdir)

flag = 0

for year in range(1981, 2011):
    filename = f"chirps-v2.0.{year}.days_p25.nc"
    src = f"{srcdir}/{filename}"
    dst = f"{dstdir}/chirps-v2.0.{year}.days_p25_{country}.nc"

    if flag != 1:
        latmin, latmax, lonmin, lonmax = get_clip_extent(targetgrid, src)
        flag=1

    print(f"Input = {src}")   
    nco_cmd = f"ncks --dmn latitude,{latmin},{latmax} --dmn longitude,{lonmin},{lonmax} {src} {dst}"
    os.system(nco_cmd)
    #print(nco_cmd)

