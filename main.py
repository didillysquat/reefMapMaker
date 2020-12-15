"""
This is a script for creating maps with reef locations
It is a modification from the CBASS84 script I used to make the figure for that paper
We will aim for it to take command inputs of the map bounding latitudes and longitudes
By default we will plot the world reefs on the map
We will also want to be able to add additional reefs to the map through some sort of csv import
"""

import os
import pandas as pd
import matplotlib as mpl
mpl.use('TKAgg')
import matplotlib.pyplot as plt
# NB the pip cartopy install seems to be broken as it doesn't install the required libararies.
# The solution was to install using conda. conda install cartopy.
# I then had to downgrade shapely to 1.5.17. pip install shapely==1.5.17
from cartopy.mpl.gridliner import Gridliner
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
import cartopy
import matplotlib.gridspec as gridspec
from matplotlib import collections, patches
import sys
import random
from matplotlib.patches import Polygon, Circle
from matplotlib.collections import PatchCollection
from matplotlib.colors import ListedColormap
import numpy as np
import pickle
from datetime import datetime
from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature
import shapely


class MapWthInsetFigure:
    def __init__(self, x1_bound=32, x2_bound=45, y1_bound=12, y2_bound=30, fig_output_path=None):
        self.root_dir = os.path.dirname(os.path.realpath(__file__))
        self.date_time = str(datetime.now()).split('.')[0].replace('-','').replace(' ','T').replace(':','')
        self.bounds = [x1_bound, x2_bound, y1_bound, y2_bound]
        self.gis_input_base_path = os.path.join(self.root_dir, 'reef_gis_input')
        self.user_input_path = os.path.join(self.gis_input_base_path, 'user_input.csv')
        self.user_circles = self._make_user_circles()
        self.reef_shape_file_path = os.path.join(
            self.gis_input_base_path, '14_001_WCMC008_CoralReefs2018_v4/01_Data/WCMC008_CoralReef2018_Py_v4.shp'
        )
        self.fig = plt.figure(figsize=(8,5))
        self.large_map_ax = plt.subplot(projection=ccrs.PlateCarree(), zorder=1)

        # figure output
        if not fig_output_path:
            self.fig_out_path = os.path.join(self.root_dir, f'map_out_{self.date_time}.png')
        else:
            if(fig_output_path.split('.')[-1] in ['svg', 'png', 'jpg', 'jpeg']):
                if not os.path.exists(os.path.abspath(fig_output_path)):
                    os.makedirs(os.path.abspath(fig_output_path))
                self.fig_out_path = fig_output_path
            else:
                raise RuntimeError(
                    "Please ensure that a full path is given for "
                    "--fig_output_path with extension .svg, .png, .jpg or .jpeg")


        self.large_map_ax.set_extent(extents=(x1_bound, x2_bound, y1_bound, y2_bound))

    def _make_user_circles(self):
        circles = []
        df = pd.read_csv(self.user_input_path)
        df.set_index(keys=['site'], inplace=True, drop=True)
        for site in df.index:
            circles.append(Circle(xy=(df.at[site, 'Longitude_(E)'], df.at[site, 'Latitude_(N)']), radius=0.05, color='red', zorder=3))
        return circles

    def draw_map(self):
        land_110m, ocean_110m, boundary_110m = self._get_naural_earth_features_big_map()
        print('drawing annotations on big map')
        self._draw_natural_earth_features_big_map(land_110m, ocean_110m, boundary_110m)
        print('big map annotations complete')
        self._put_gridlines_on_large_map_ax()
        # TODO needs to be made dynamic and linked to input
        # self._annotate_map_with_sites()
        # We are going to work wi th this data set for the reefs
        # https://data.unep-wcmc.org/datasets/1
        # We are modifying the code from the Restrepo et al paper that I wrote

        print('here')
        self._add_reference_reefs()
        self._add_user_reefs()

        print('saving .png')
        plt.savefig(self.fig_out_path, dpi=600)

    def _add_user_reefs(self):
        """
        Add the user supplied reefs.
        Currently as Circles, but we should enable polygons too.
        """
        for circle in self.user_circles:
            self.large_map_ax.add_patch(circle)
        foo = 'bar'

    def _add_reference_reefs(self):
        # We were working with shapely features and adding geometries but there were so many problems
        # I advise to stay away form trying to get them working again.
        # Instead we are now creating individual Polygon patches from the coordinates
        # and adding them to the plot
        print('We reach here before iter')
        reader = Reader(self.reef_shape_file_path)
        count = 0
        # records_list = [r for r in reader.records() if r.geometry is not None]
        for r in reader.records(): # reader.records() produces a generator
            if r.bounds[0] > self.bounds[0]:
                if r.bounds[1] > self.bounds[2]:
                    if r.bounds[2] < self.bounds[1]:
                        if r.bounds[3] < self.bounds[3]:
                            # if ('SAU' in r.attributes['ISO3']):
                            try:
                                if r.geometry.geom_type.lower() == 'multipolygon':
                                    # Create multiple matplotlib polygon objects from the multiple shape polygons
                                    for polygon in r.geometry:
                                        coords = polygon.exterior.coords.xy
                                    # each of the individual coords is a tup of tups
                                        reef_poly = Polygon(
                                            [(x, y) for x, y in zip(list(coords[0]), list(coords[1]))], closed=True,
                                            fill=True, edgecolor='#003366', facecolor='#003366',
                                            alpha=1, zorder=2)
                                        self.large_map_ax.add_patch(reef_poly)
                                elif r.geometry.geom_type.lower() == 'polygon':
                                    coords = r.geometry.exterior.coords.xy
                                    reef_poly = Polygon([(x,y) for x, y in zip(list(coords[0]), list(coords[1]))], closed=True, fill=True, edgecolor='#003366', facecolor='#003366',
                                                        alpha=1, zorder=2)
                                    self.large_map_ax.add_patch(reef_poly)
                                else:
                                    foo = 'bar'
                            except Exception as e:
                                # The common error that occurs is Unexpected Error: unable to find ring point
                                # We've given up trying to catch this is a more elegant way
                                # The class of exception raised is Exception in shapefile.py
                                count += 1
                                print(e)
                                continue

        print(f'{count} error producing records that were in bounds were counted')

    def _get_naural_earth_features_big_map(self):
        land_110m = cartopy.feature.NaturalEarthFeature(category='physical', name='land',
                                                        scale='50m')
        ocean_110m = cartopy.feature.NaturalEarthFeature(category='physical', name='ocean',
                                                         scale='50m')
        boundary_110m = cartopy.feature.NaturalEarthFeature(category='cultural',
                                                            name='admin_0_boundary_lines_land', scale='110m')
        return land_110m, ocean_110m, boundary_110m

    def _get_naural_earth_features_zoom_map(self):
        land_10m = cartopy.feature.NaturalEarthFeature(category='physical', name='land',
                                                       scale='50m')
        ocean_10m = cartopy.feature.NaturalEarthFeature(category='physical', name='ocean',
                                                        scale='50m')

        return land_10m, ocean_10m

    def _draw_natural_earth_features_big_map(self, land_110m, ocean_110m, boundary_110m):
        """NB the RGB must be a tuple in a list and the R, G, B must be given as a value between 0 and 1"""
        # self.large_map_ax.add_feature(land_110m, facecolor=[(238 / 255, 239 / 255, 219 / 255)],
        #                               edgecolor='black', linewidth=0.2)
        self.large_map_ax.add_feature(land_110m, facecolor='white',
                                      edgecolor='black', linewidth=0.2, zorder=1)
        self.large_map_ax.add_feature(ocean_110m, facecolor=[(136 / 255, 182 / 255, 224 / 255)],
                                      edgecolor='black', linewidth=0.2)
        self.large_map_ax.add_feature(boundary_110m, edgecolor='gray', linewidth=0.2, facecolor='None')

    def _put_gridlines_on_large_map_ax(self):
        """ Although there is a GeoAxis.gridlines() method, this method does not yet allow a lot of
        bespoke options. If we want to only put the labels on the top and left then we have to
        generate a Gridliner object (normally returned by GeoAxis.gridlines() ourselves. We then need
        to manually change the xlabels_bottom and ylabels_right attributes of this Gridliner object.
        We then draw it by adding it to the GeoAxis._gridliners list.
        TODO we need to make these locations dynamic with the input bounds"""
        xlocs = [34.0, 38.0, 42.0]
        ylocs = [13.0, 17.0, 21.0, 25.0, 29.0]

        if xlocs is not None and not isinstance(xlocs, mticker.Locator):
            xlocs = mticker.FixedLocator(xlocs)
        if ylocs is not None and not isinstance(ylocs, mticker.Locator):
            ylocs = mticker.FixedLocator(ylocs)
        g1 = Gridliner(
            axes=self.large_map_ax, crs=ccrs.PlateCarree(), draw_labels=True,
            xlocator=xlocs, ylocator=ylocs)
        g1.top_labels = False
        g1.right_labels = False
        self.large_map_ax._gridliners.append(g1)
        # self.large_map_ax.gridlines(draw_labels=True)

    # def _annotate_map_with_sites(self):
    #     # TODO tie this in to an input of some sort so that it can be user provided
    #     for site in ['ICN']:
    #         if site != 'PrT':
    #             self.large_map_ax.plot(self.sites_location_dict[site][0], self.sites_location_dict[site][1],
    #                                    self.site_marker_dict[site], markerfacecolor=self.site_color_dict[site], markeredgecolor='black', markersize=6, markeredgewidth=0.2)
    #         else:
    #             self.large_map_ax.plot(self.sites_location_dict[site][0], self.sites_location_dict[site][1],
    #                                    self.site_marker_dict[site], markerfacecolor='black', markeredgecolor='black',
    #                                    markersize=8)



# If full then the whole Red Sea length will be plotted
mwif = MapWthInsetFigure()
mwif.draw_map()

