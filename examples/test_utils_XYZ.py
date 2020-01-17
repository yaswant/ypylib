# utils.XYZ(...).method(...)
from ypylib.utils import XYZ
from numpy.random import normal, randint
import cartopy.crs as ccrs

# Generate random data N samples
N = 200  # number of samples
x = randint(60, 100, N)  # horizontal (x) limit
y = randint(0, 40, N)  # vertical (y) limit
z = normal(3, 1, N)  # random values with mean=3, std=1

# Instantiate XYZ class object and store (x, y, z) triplets
xyz = XYZ(x, y, z)

# 1. Simple plots
xyz.plot().show()  # line plot ------------------------------------------------
xyz.scatter().show()  # scatter plot ------------------------------------------
xyz.hexbin().show()  # hexbin plot --------------------------------------------
xyz.contour(delta=(5, 5)).show()  # contour plot ------------------------------
xyz.pcolormesh().show()  # colormesh plot -------------------------------------

# 2. Plot on map (Cartopy or Basemap) coordinates -----------------------------
# 2a. Resample on a 2x2 degree grid and show on a map with Robinson
# projection
plt = xyz.mapdata(
    delta=(2, 2),  # grid data on 2x2 degree resolution
    mask_land=True,
    stat='mean',  # specify which stat to genrate on grid (mean, std, etc..)
    plt_type='pcolormesh',  # this is default anyway
    projection=ccrs.Robinson(),  # use cartopy Robinson projection
    show_datapoints=True,  # show original data point locations
    describe_data=True)  # add basic stats of original data
plt.show()

# 2b. Same as above, but using basemap. ---------------------------------------
# There is a slight discrepancy here that the map limit is not respected for
# non-rectangular projection
plt = xyz.mapdata(
    delta=(2, 2),
    mask_land=True,
    stat='mean',
    plt_type='pcolormesh',
    use_cartopy=False, projection='robin',  # turn off cartopy
    gspacing=(45, 45),  # gridline spacing
    show_datapoints=False,
    describe_data=True)
plt.show()
