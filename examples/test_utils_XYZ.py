""" Demonstrate use of XYZ(x, y, z).mapdata(...)"""

from ypylib.utils import XYZ
from numpy.random import normal, randint
import cartopy.crs as ccrs

# - Generate random data N samples
N = 500  # number of samples
x = randint(60, 100, N)  # horizontal (x) limit
y = randint(0, 40, N)  # vertical (y) limit
z = normal(3, 1, N)  # random values with mean=3, std=1

# - Instantiate XYZ class object and store (x, y, z) triplets
xyz = XYZ(x, y, z)

# 1. Simple plots  (no gridding is applied to original data)
xyz.plot().show()  # -- line plot
xyz.scatter().show()  # -- scatter plot
xyz.hexbin().show()  # -- hexbin plot
xyz.contour(delta=(5, 5)).show()  # -- contour plot
xyz.pcolormesh().show()  # -- colormesh plot

# 2. Mapped plots (grid and project data on a specific mapping coordinate)

# 2a. Grid data at 2x2 degree resolution and show on map with Robinson
# projection
plt = xyz.mapdata(
    delta=(2, 2),  # grid data on 2x2 degree resolution; can be a scalar
    fillcontinents=True,  # fill continents with solid colour (under data)
    # mask_land=True,  # mask land surface with solid colour (over data)
    filloceans=True,  # fill oceans with solid colour (under data)
    # mask_ocean=True,  # mask ocean regions with solid colour (over data)
    # stat='mean',  # 'count', 'min', 'max', 'std', 'sum'
    plt_type='pcolormesh',  # 'contour', 'hexbin', 'scatter'
    projection=ccrs.Robinson(),  # use cartopy Robinson projection
    # show_datapoints=True,  # show original data point locations
    # show_gridpoints=True,  # show grid points of gridded data
    describe_data=True  # show basic stats from original data
)
plt.show()


# 2b. Same as above, but using basemap. ---------------------------------------
# Note: There is a slight discrepancy here that the map limit is not respected
# pseudo cylindrical projections
plt = xyz.mapdata(
    delta=2,
    fillcontinents=True,
    # mask_land=True,
    filloceans=True,
    # mask_ocean=True,
    # stat='mean',
    plt_type='pcolormesh',
    use_cartopy=False,  # must use this flag to use Basemap
    projection='robin',  # use Basemap projection notations
    gspacing=(45, 45),  # gridline spacing
    # show_datapoints=True,
    # show_gridpoints=True,
    describe_data=True
)
plt.show()
