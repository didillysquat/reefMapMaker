#!/usr/bin/env python3
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
import argparse
# NB the pip cartopy install seems to be broken as it doesn't install the required libararies.
# The solution was to install using conda. conda install cartopy.
# I then had to downgrade shapely to 1.5.17. pip install shapely==1.5.17
from cartopy.mpl.gridliner import Gridliner
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
import cartopy
from matplotlib.patches import Polygon, Circle
from datetime import datetime
from cartopy.io.shapereader import Reader
import sys

__version__ = "v0.1.0"

class MapWthInsetFigure:
    def __init__(self):
        print(f"""
        
                         __ __  __             __  __       _             
                        / _|  \/  |           |  \/  |     | |            
          _ __ ___  ___| |_| \  / | __ _ _ __ | \  / | __ _| | _____ _ __ 
         | '__/ _ \/ _ \  _| |\/| |/ _` | '_ \| |\/| |/ _` | |/ / _ \ '__|
         | | |  __/  __/ | | |  | | (_| | |_) | |  | | (_| |   <  __/ |   
         |_|  \___|\___|_| |_|  |_|\__,_| .__/|_|  |_|\__,_|_|\_\___|_|   
                                        | | {__version__}                              
                                        |_|                               

        """)
        print("""
            reefMapMaker uses the Global Distribution of Coral Reefs data set.
            
            From https://data.unep-wcmc.org/datasets/1:
            
            This dataset shows the global distribution of coral reefs in tropical and subtropical regions.
            It is the most comprehensive global dataset of warm-water coral reefs to date, acting as a
            foundation baseline map for future, more detailed, work.
            This dataset was compiled from a number of sources by UNEP World Conservation Monitoring Centre
            (UNEP-WCMC) and the WorldFish Centre, in collaboration with WRI (World Resources Institute) and
            TNC (The Nature Conservancy). Data sources include the Millennium Coral Reef Mapping Project
            (IMaRS-USF and IRD 2005, IMaRS-USF 2005) and the World Atlas of Coral Reefs (Spalding et al. 2001).
            
            When using maps containing the reference coral reefs, please enusre you include the relevant citations
            to this UNEP-WCMC resource.
            
            From https://data.unep-wcmc.org/datasets/1:
            Citations:
            UNEP-WCMC, WorldFish Centre, WRI, TNC (2018). Global distribution of warm-water coral reefs, compiled from 
            multiple sources including the Millennium Coral Reef Mapping Project. Version 4.0. Includes contributions 
            from IMaRS-USF and IRD (2005), IMaRS-USF (2005) and Spalding et al. (2001). Cambridge (UK): UN Environment 
            World Conservation Monitoring Centre. URL: http://data.unep-wcmc.org/datasets/1\n
    
            Citations for the separate entities:
            IMaRS-USF (Institute for Marine Remote Sensing-University of South Florida) (2005). Millennium Coral Reef 
            Mapping Project. Unvalidated maps. These maps are unendorsed by IRD, but were further interpreted by UNEP 
            World Conservation Monitoring Centre. Cambridge (UK): UNEP World Conservation Monitoring Centre
            
            IMaRS-USF, IRD (Institut de Recherche pour le Developpement) (2005). Millennium Coral Reef Mapping Project. 
            Validated maps. Cambridge (UK): UNEP World Conservation Monitoring Centre
            
            Spalding MD, Ravilious C, Green EP (2001). World Atlas of Coral Reefs. Berkeley (California, USA): 
            The University of California Press. 436 pp.\n\n\n\n
        """)
        parser = argparse.ArgumentParser(
            description='A script to make maps with annotated coral reef locations',
            epilog='For support email: didillysquat@gmail.com')
        self._define_additional_args(parser)
        self.args = parser.parse_args()
        self.root_dir = os.path.dirname(os.path.realpath(__file__))
        self.date_time = str(datetime.now()).split('.')[0].replace('-','').replace(' ','T').replace(':','')

        # User input
        user_input_path = os.path.join(self.args.user_map_config_sheet_path)
        df_user_map_input = pd.read_excel(user_input_path, sheet_name='user_map_input', index_col=0)
        df_user_point_input = pd.read_excel(user_input_path, sheet_name='user_point_input', index_col=0)
        self.config_dict = {k: v for k, v in zip(df_user_map_input.index.values, df_user_map_input['value'].values)}
        print('Checking input config file')
        self._check_config_input_is_valid()
        self.circles_dict = {
            site: Circle(
                xy=(df_user_point_input.at[site, 'longitude_deg_e'], df_user_point_input.at[site, 'latitude_deg_n']),
                radius=df_user_point_input.at[site, 'radius_in_deg_lat'],
                color=df_user_point_input.at[site, 'color'], zorder=3
            ) for site in df_user_point_input.index
        }

        self.bounds = [self.config_dict['x1_bound'], self.config_dict['x2_bound'], self.config_dict['y1_bound'], self.config_dict['y2_bound']]
        self.gis_input_base_path = os.path.join(self.root_dir, 'reef_gis_input')

        # reference reefs
        self.reference_reef_shape_file_path = os.path.join(
            self.gis_input_base_path, '14_001_WCMC008_CoralReefs2018_v4/01_Data/WCMC008_CoralReef2018_Py_v4.shp'
        )
        if not os.path.exists(self.reference_reef_shape_file_path):
            raise RuntimeError(f'unable to find reference reef shape file.\n'
                               f'Please ensure that you have downloaded the 14_001_WCMC008_CoralReefs2018_v4 '
                               f'dataset from https://data.unep-wcmc.org/datasets/1 and have placed decompressed '
                               f'directory in the reef_gis_input directory.')

        # setup the map figure
        #TODO make the figure at the same ratios as the bounds are given so that we don't end up
        # with a whole load of white space on the image.
        self.fig = plt.figure(figsize=(8,5))
        self.large_map_ax = plt.subplot(projection=ccrs.PlateCarree(), zorder=1)
        self.large_map_ax.set_extent(extents=(self.bounds[0], self.bounds[1], self.bounds[2], self.bounds[3]))

        # figure output paths
        if not self.args.figure_output_directory:
            self.fig_out_path_svg = os.path.join(self.root_dir, f'map_out_{self.date_time}.svg')
            self.fig_out_path_png = os.path.join(self.root_dir, f'map_out_{self.date_time}.png')
        else:
            if not os.path.exists(os.path.abspath(self.args.figure_output_directory)):
                os.makedirs(os.path.abspath(self.args.figure_output_directory))
            self.fig_out_path_svg = os.path.join(self.args.figure_output_directory, f'map_out_{self.date_time}.svg')
            self.fig_out_path_png = os.path.join(self.args.figure_output_directory, f'map_out_{self.date_time}.png')

    def _check_config_input_is_valid(self):
        # lat long
        self._check_bounds()
        self._check_bool_params()
        self._set_default_color_params()
        self._check_gridline_labels()
        self._check_grid_line_coords()

    def _check_grid_line_coords(self):
        if self.config_dict['plot_grid_lines']:
            try:
                assert (not pd.isnull(self.config_dict['grid_line_x_position']))
                assert (len(self.config_dict['grid_line_x_position'].split(',')) > 0)
            except AssertionError:
                raise RuntimeError("Please provide valid x coordinates for your gridlines")
            try:
                assert (not pd.isnull(self.config_dict['grid_line_y_position']))
                assert (len(self.config_dict['grid_line_y_position'].split(',')) > 0)
            except AssertionError:
                raise RuntimeError("Please provide valid y coordinates for your gridlines")

    def _check_gridline_labels(self):
        try:
            assert (self.config_dict['grid_line_x_label_position'] in ['bottom', 'top'])
        except AssertionError:
            raise RuntimeError("grid_line_x_label_position must be 'bottom' or 'top'")
        try:
            assert (self.config_dict['grid_line_y_label_position'] in ['left', 'right'])
        except AssertionError:
            raise RuntimeError("grid_line_y_label_position must be 'left' or 'right'")

    def _set_default_color_params(self):
        # set default color params
        default_color_dict = {'sea_color': "#88b5e0", 'land_color': 'white', 'reference_reefs_color': "#003366"}
        for c_param in default_color_dict.keys():
            if (pd.isnull(self.config_dict[c_param]) or self.config_dict[c_param] == 'AUTO'):
                self.config_dict[c_param] = default_color_dict[c_param]

    def _check_bool_params(self):
        for bool_param in ['plot_sea', 'plot_reference_reefs', 'plot_land', 'plot_grid_lines', 'plot_boundaries']:
            try:
                assert (type(self.config_dict[bool_param]) is bool)
            except AssertionError:
                raise RuntimeError(f"{bool_param} must be either TRUE or FALSE.")

    def _check_bounds(self):
        try:
            assert (self.config_dict['x1_bound'] < self.config_dict['x2_bound'])
        except AssertionError:
            raise RuntimeError('Check the format of your user_map_config_sheet.\n'
                               'x1_bound should be less than x2_bound')
        try:
            assert (self.config_dict['y1_bound'] < self.config_dict['y2_bound'])
        except AssertionError:
            raise RuntimeError('Check the format of your user_map_config_sheet.\n'
                               'y1_bound should be less than y2_bound')
        try:
            for xbound in ['x1_bound', 'x2_bound']:
                assert ((self.config_dict[xbound] < 180) and (self.config_dict[xbound] > -180))
        except AssertionError:
            raise RuntimeError("One of your xbounds appears to be an invalid value.")
        try:
            for ybound in ['y1_bound', 'y2_bound']:
                assert ((self.config_dict[ybound] < 90) and (self.config_dict[ybound] > -90))
        except AssertionError:
            raise RuntimeError("One of your ybounds appears to be an invalid value.")

    @staticmethod
    def _define_additional_args(parser):
        parser.add_argument(
            '--user_map_config_sheet_path',
            help='The full path to the .xlsx file that contains the configurations for the map',
            required=True
        )
        parser.add_argument(
            '--figure_output_directory',
            help='The full path to the directory where the map output figures will be saved. '
                 'A .svg and a .png file will be created with a time stamp.',
            default=None
        )

    def draw_map(self):
        land_110m, ocean_110m, boundary_110m = self._get_naural_earth_features_big_map()
        print('Drawing annotations on map\n')
        self._draw_natural_earth_features_big_map(land_110m, ocean_110m, boundary_110m)
        print('Annotations complete\n')
        if self.config_dict['plot_grid_lines']:
            self._put_gridlines_on_large_map_ax()
        # TODO needs to be made dynamic and linked to input
        # self._annotate_map_with_sites()
        # We are going to work wi th this data set for the reefs
        # https://data.unep-wcmc.org/datasets/1
        # We are modifying the code from the Restrepo et al paper that I wrote
        self._add_reference_reefs()
        self._add_user_reefs()

        print(f'saving to {self.fig_out_path_png}')
        plt.savefig(self.fig_out_path_png, dpi=600)
        print(f'saving to {self.fig_out_path_svg}')
        plt.savefig(self.fig_out_path_svg, dpi=600)

    def _add_user_reefs(self):
        """
        Add the user supplied reefs.
        Currently as Circles, but we should enable polygons too.
        """
        print('plotting user reefs\n')
        for site, circle in self.circles_dict.items():
            print(f'plotting user reef {site}')
            self.large_map_ax.add_patch(circle)

    def _add_reference_reefs(self):
        # We were working with shapely features and adding geometries but there were so many problems
        # I advise to stay away form trying to get them working again.
        # Instead we are now creating individual Polygon patches from the coordinates
        # and adding them to the plot
        print('Annotating reference reefs\n')
        reader = Reader(self.reference_reef_shape_file_path)
        count = 0
        reef_count = 0
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
                                            fill=True, color=self.config_dict['reference_reefs_color'],
                                            alpha=1, zorder=2)
                                        self.large_map_ax.add_patch(reef_poly)
                                        reef_count += 1
                                        if reef_count % 100 == 0:
                                            print(f'{reef_count} reference reefs plotted')
                                elif r.geometry.geom_type.lower() == 'polygon':
                                    coords = r.geometry.exterior.coords.xy
                                    reef_poly = Polygon([(x,y) for x, y in zip(list(coords[0]), list(coords[1]))], closed=True, fill=True, color=self.config_dict['reference_reefs_color'],
                                                        alpha=1, zorder=2)
                                    self.large_map_ax.add_patch(reef_poly)
                                    reef_count += 1
                                    if reef_count % 100 == 0:
                                        print(f'{reef_count} reference reefs plotted')
                            except Exception as e:
                                # The common error that occurs is Unexpected Error: unable to find ring point
                                # We've given up trying to catch this is a more elegant way
                                # The class of exception raised is Exception in shapefile.py
                                count += 1
                                continue


        print(f'\n{count} error producing records were discounted from the reference reefs')

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
        if self.config_dict['plot_land']:
            self.large_map_ax.add_feature(land_110m, facecolor=self.config_dict['land_color'],
                                          edgecolor='black', linewidth=0.2, zorder=1)
        if self.config_dict['plot_sea']:
            self.large_map_ax.add_feature(ocean_110m, facecolor=self.config_dict['sea_color'],
                                          edgecolor='black', linewidth=0.2)
        if self.config_dict['plot_boundaries']:
            self.large_map_ax.add_feature(boundary_110m, edgecolor='gray', linewidth=0.2, facecolor='None')

    def _put_gridlines_on_large_map_ax(self):
        """ Although there is a GeoAxis.gridlines() method, this method does not yet allow a lot of
        bespoke options. If we want to only put the labels on the top and left then we have to
        generate a Gridliner object (normally returned by GeoAxis.gridlines() ourselves. We then need
        to manually change the xlabels_bottom and ylabels_right attributes of this Gridliner object.
        We then draw it by adding it to the GeoAxis._gridliners list."""

        xlocs = mticker.FixedLocator([float(_) for _ in self.config_dict['grid_line_x_position'].split(',')])
        ylocs = mticker.FixedLocator([float(_) for _ in self.config_dict['grid_line_y_position'].split(',')])
        g1 = Gridliner(
            axes=self.large_map_ax, crs=ccrs.PlateCarree(), draw_labels=True,
            xlocator=xlocs, ylocator=ylocs)
        if self.config_dict['grid_line_x_label_position'] == 'bottom':
            g1.top_labels = False
        else:
            g1.bottom_labels = False
        if self.config_dict['grid_line_y_label_position'] == 'left':
            g1.right_labels = False
        else:
            g1.left_labels = False
        self.large_map_ax._gridliners.append(g1)


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

