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
from matplotlib.collections import PatchCollection
from datetime import datetime
from cartopy.io.shapereader import Reader
import sys
import numpy as np
import time

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
        self.root_dir = os.getcwd()
        self.date_time = str(datetime.now()).split('.')[0].replace('-','').replace(' ','T').replace(':','')

        # User input
        # We will allow user configurations through both the command line and a config file
        # We should prioritise command line input
        # Probably the clearest and most logical way to do this will be to go parameter by parameter
        self.confic_parameters = ['bounds', 'plot_sea', 'sea_color', 'plot_reference_reefs',
                'reference_reef_color', 'reference_reef_patch_type', 'plot_land', 'land_color',
                'plot_grid_lines', 'grid_line_x_position', 'grid_line_y_position',
                'grid_line_x_label_position', 'grid_line_y_label_position', 'plot_boundaries']
        if self.args.config_sheet:
            # Config sheet
            if self.args.config_sheet.endswith('.tsv'):
                config_df = pd.read_csv(self.args.config_sheet, index_col=0, sep='\t')
            elif self.args.config_sheet.endswith('.xlsx'):
                config_df = pd.read_excel(self.args.config_sheet, sheet_name='user_config', index_col=0)
            self.config_dict = {k: v for k, v in zip(config_df.index.values, config_df['value'].values)}

            self._setup_config()

            self.df_user_point_input = pd.read_excel(
                self.args.site_sheet, sheet_name='user_site', index_col=0
            ).sort_values('radius_in_deg_lat', axis=0, ascending=False)
            print('Checking input config file')
            self._setup_config()
        else:
            # no config sheet provided
            # Use default dict
            self.config_dict = {}
            if self.args.bounds:
                bounds = [float(_) for _ in self.args.bounds.split(',')]
                self.config_dict.update(
                    {'x1_bound': bounds[0], 'x2_bound': bounds[1], 'y1_bound': bounds[2], 'y2_bound': bounds[3]}
                )
            else:
                self.config_dict.update({'x1_bound': -180, 'x2_bound': 180, 'y1_bound': -90, 'y2_bound': 90})
            self._check_bounds()
            self.config_dict.update({
                'plot_sea': True, 'sea_color': "#88b5e0", 'plot_reference_reefs': True,
                'reference_reef_color': '#003366', 'plot_land': True, 'land_color': 'white',
                'plot_grid_lines': True, 'grid_line_x_position': None, 'grid_line_y_position': None,
                'grid_line_x_label_position': 'bottom', 'grid_line_y_label_position': 'left', 'plot_boundaries': True
            })
            self.user_site_dict = {}
            self.df_user_point_input = None

        self.bounds = [
            self.config_dict['x1_bound'], self.config_dict['x2_bound'],
            self.config_dict['y1_bound'], self.config_dict['y2_bound']
        ]

        self.reference_reef_shape_file_path = self._find_shape_path()

        self.fig = self._setup_map_figure()

        self.large_map_ax = plt.subplot(projection=ccrs.PlateCarree(), zorder=1)
        self.large_map_ax.set_extent(extents=(self.bounds[0], self.bounds[1], self.bounds[2], self.bounds[3]))
        # Scalar for converting the user inputted coordinate radii to point format for plotting
        self.coord_to_point_scaler = self._calc_scaler()

        # figure output paths
        if not self.args.fig_out_dir:
            self.fig_out_path_svg = os.path.join(self.root_dir, f'map_out_{self.date_time}.svg')
            self.fig_out_path_png = os.path.join(self.root_dir, f'map_out_{self.date_time}.png')
        else:
            if not os.path.exists(os.path.abspath(self.args.fig_out_dir)):
                print(f'Creating {self.args.fig_out_dir}')
                os.makedirs(os.path.abspath(self.args.fig_out_dir))
            self.fig_out_path_svg = os.path.join(self.args.fig_out_dir, f'map_out_{self.date_time}.svg')
            self.fig_out_path_png = os.path.join(self.args.fig_out_dir, f'map_out_{self.date_time}.png')

    def _calc_scaler(self):
        """
        Users input the size of their reefs to be plotted in coordinate system radii, i.e. 0.01 deg
        We need to convert these values to point values that can then be squared to use the scatter
        plot size unit for plotting.
        The transformations required are quite mind bending so I have included some helpful links here:
        https://matplotlib.org/3.1.1/tutorials/advanced/transforms_tutorial.html
        https://stackoverflow.com/a/47403507/5516420
        https://stackoverflow.com/a/33945420/5516420

        Our approach is to use a set of data points that are 1 deg lon apart
        We convert these to display units and therefore workout what 1 deg lon is in display units (pixels)
        We then convert this pixel length to a point length with:
        points = (pixels * 72) / dpi
        This length is then the scalar and will need to be squared when supplied to the size argument of scatter
        """
        # points in data coordinates
        p1_da = (5, 5)
        p2_da = (6, 5)

        # get the points in display coordinates
        di = self.large_map_ax.transData.transform([p1_da, p2_da])

        length_in_pixels = di[1, 0] - di[0, 0]
        return (length_in_pixels * 72) / self.fig.dpi

    def _setup_map_figure(self):
        """
        Set fig size ratios according to lat lon ratios
        """
        big_fig_size = 10
        lat = self.bounds[3] - self.bounds[2]
        lon = self.bounds[1] - self.bounds[0]
        # figsize is w x h
        if lat > lon:
            fig = plt.figure(figsize=((lon / lat) * big_fig_size, big_fig_size))
        else:
            fig = plt.figure(figsize=(big_fig_size, (lat / lon) * big_fig_size))
        return fig

    def _find_shape_path(self):
        """
        If path supplied by user, find the shapefile path and set self.reference_reef_shape_file_path.
        else, search for the shape file in the current working directory.
        """
        if self.args.ref_reef_dir:
            # The user has supplied a path
            # This could be the directory containing the parent directory of the dataset
            # It could be the parent directory of the dataset itself
            if 'WCMC008' in self.args.ref_reef_dir:
                potential_shape_path = os.path.join(self.args.ref_reef_dir,
                                                    '01_Data/WCMC008_CoralReef2018_Py_v4.shp')
                if os.path.exists(potential_shape_path):
                    return potential_shape_path
            else:
                # Search for the directory in the supplied path
                dir_to_search = self.args.ref_reef_dir
                return self._search_for_ref_reef_parent_data_dir(dir_to_search=dir_to_search)
        else:
            # No user supplied path, search for the directory in the current working directory
            dir_to_search = self.root_dir
            return self._search_for_ref_reef_parent_data_dir(dir_to_search=dir_to_search)

    def _search_for_ref_reef_parent_data_dir(self, dir_to_search):
        candidate_dirs = []
        for (dirpath, dirnames, filenames) in os.walk(dir_to_search):
            candidate_dirs.extend([os.path.join(dirpath, _) for _ in dirnames if 'WCMC008' in _])

        if len(candidate_dirs) == 1:
            # Then we have found the parent directory
            data_set_parent_dir = candidate_dirs[0]
        else:
            self._report_unable_to_find_ref_reef_path()
        # Verify that the shape file exists
        potential_shape_path = os.path.join(data_set_parent_dir, '01_Data/WCMC008_CoralReef2018_Py_v4.shp')
        if os.path.exists(potential_shape_path):
            return potential_shape_path
        else:
            self._report_unable_to_find_ref_reef_path()

    def _report_unable_to_find_ref_reef_path(self):
        raise RuntimeError(
            "Could not automatically find the reference reef dataset. "
            "Please specify the directory of the dataset on the command line using --ref-reef-dir")

    def _setup_config(self):
        """
        Check the values of the user config inputs
        We will assume that the user will supply either command line arguments
        OR a .tsv/.xlsx config_sheet.
        In the case where a config parameter is referenced in both of these,
        we will use the command line supplied argument but produce a warning for the user that their config sheet
        setting is being overwritten.
        Where no value has been provided for a given config parameter, we will inform the user of this and
        provide the value for the given param to be used.
        :return:
        """
        # lat long
        self._check_bounds()
        self._check_bool_params()
        self._set_default_color_params()
        self._check_gridline_labels()
        self._check_grid_line_coords()
        self._check_ref_reef_path_type()

    def _check_ref_reef_path_type(self):
        if pd.isnull(self.config_dict['reference_reef_patch_type']) and self.config_dict['plot_reference_reefs']:
            print('plotting reefs using default of "polygon".\n'
                  'If you would rather plot as points, please set the reference_reef_patch_type to "point"')
        else:
            raise RuntimeError("Invalid value provided for reference_reef_patch_type.\n"
                               "Valid values are 'polygon' or 'point'.")

    def _check_grid_line_coords(self):
        if self.config_dict['plot_grid_lines']:
            try:
                assert(not pd.isnull(self.config_dict['grid_line_x_position']))
                assert(len(self.config_dict['grid_line_x_position'].split(',')) > 0)
                assert (not pd.isnull(self.config_dict['grid_line_y_position']))
                assert (len(self.config_dict['grid_line_y_position'].split(',')) > 0)
            except AssertionError:
                print("WARNING: There is a problem with the formatting of your grid line coordinates.")
                print('Default coordinates will be used.')
                self.config_dict['grid_line_x_position'] = None
                self.config_dict['grid_line_y_position'] = None

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
        default_color_dict = {'sea_color': "#88b5e0", 'land_color': 'white', 'reference_reef_color': "#003366"}
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
            raise RuntimeError('Check the format of your user config sheet.\n'
                               'x1_bound should be less than x2_bound')
        try:
            assert (self.config_dict['y1_bound'] < self.config_dict['y2_bound'])
        except AssertionError:
            raise RuntimeError('Check the format of your user config sheet.\n'
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
            '--config-sheet',
            help='The full path to the .xlsx file that contains the configurations for the map',
            required=False
        )
        parser.add_argument(
            '--fig-out-dir',
            help='The full path to the directory where the map output figures will be saved. '
                 'A .svg and a .png file will be created with a time stamp.',
            default=None
        )
        parser.add_argument(
            '--ref-reef-dir',
            help='The full path to the directory containing the reference reef shapefile data.'
                 'Default is current working directory. The dataset can be downloaded from: '
                 'https://data.unep-wcmc.org/datasets/1',
            default=None
        )
        parser.add_argument(
            '--bounds',
            help='Comma delimited  coordinate boundaries in decimal degrees (N/E)\n'
                 'in the order of westernmost, easternmost, southermost, northernmost\n.'
                 'E.g. For bounds that encapsulate the Red Sea: 32,45,12,30',
            default=None
        )
        parser.add_argument(
            '--points',
            help='When passed, reference reef will be plotted as a series of cirlces rather than Polygons.\n'
                 'This can be helpful for maps that cover a larger area as'
                 ' stroking paths can create some graphical aftefacts.',
            action='store_true'
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
        self._add_user_reefs()
        self._add_reference_reefs()


        print(f'saving to {self.fig_out_path_png}')
        plt.savefig(self.fig_out_path_png, dpi=600)
        print(f'saving to {self.fig_out_path_svg}')
        plt.savefig(self.fig_out_path_svg, dpi=600)

    def _add_user_reefs(self):
        """
        Add the user supplied reefs.
        Currently as scatter, but we should enable polygons too.
        Line widths should be proportional to the marker size
        We should attempt to add the points in order of the largest first to minimise overlap
        """
        print('plotting user reefs\n')
        line_widths = ((self.df_user_point_input['radius_in_deg_lat'].astype(
                float) * 2) * self.coord_to_point_scaler) * 0.1
        self.large_map_ax.scatter(
            x=self.df_user_point_input['longitude_deg_e'],
            y=self.df_user_point_input['latitude_deg_n'],
            s=(((self.df_user_point_input['radius_in_deg_lat'].astype(
                float) * 2) * self.coord_to_point_scaler) ** 2),
            facecolors=self.df_user_point_input['facecolor'],
            edgecolors=self.df_user_point_input['edgecolor'], zorder=3,
            linewidths=line_widths
        )


    def _add_reference_reefs(self):
        # We were working with shapely features and adding geometries but there were so many problems
        # I advise to stay away form trying to get them working again.
        # Instead we are now creating individual Polygon patches from the coordinates
        # and adding them to the plot
        print('Annotating reference reefs\n')
        reader = Reader(self.reference_reef_shape_file_path)
        error_count = 0
        reef_count = 0
        checked_count = 0
        point_coords_x = []
        point_coords_y = []
        start_time = time.time()
        for r in reader.records():  # reader.records() produces a generator
            try:
                if r.geometry.geom_type.lower() == 'multipolygon':
                    # Create multiple matplotlib polygon objects from the multiple shape polygons
                    for polygon in r.geometry:
                        checked_count += 1
                        if checked_count % 1000 == 0:
                            new_time = time.time()
                            print(f'{checked_count} reference reefs checked in {new_time-start_time}s')
                            # if checked_count == 20000:
                            #     self._add_and_make_ref_reef_poly(coords=(point_coords_x, point_coords_y))
                        if self._if_within_bounds(polygon.bounds):
                            # each of the individual coords is a tup of tups
                            coords = polygon.exterior.coords.xy
                            if self.args.points:
                                # Then we want to collect all of the xy points to plot as a scatter
                                # and plot them as a single scatter at the end
                                point_coords_x.append(np.array(coords[0]))
                                point_coords_y.append(np.array(coords[1]))
                            else:
                                # we want to plot the polygons as we go
                                self._make_and_add_ref_reef_patches(coords)
                            reef_count += 1
                            if reef_count % 100 == 0:
                                print(f'{reef_count} reference reefs plotted')
                                # if checked_count == 20000:
                                #     self._add_and_make_ref_reef_poly(coords=(point_coords_x, point_coords_y))
                elif r.geometry.geom_type.lower() == 'polygon':
                    checked_count += 1
                    if checked_count % 1000 == 0:
                        new_time = time.time()
                        print(f'{checked_count} reference reefs checked in {new_time - start_time}s')
                    if self._if_within_bounds(r.bounds):
                        coords = r.geometry.exterior.coords.xy
                        if self.args.points:
                            # Then we want to collect all of the xy points to plot as a scatter
                            # and plot them as a single scatter at the end
                            point_coords_x.append(np.array(coords[0]))
                            point_coords_y.append(np.array(coords[1]))
                        else:
                            # we want to plot the polygons as we go
                            self._make_and_add_ref_reef_patches(coords)
                        reef_count += 1
                        if reef_count % 100 == 0:
                            print(f'{reef_count} reference reefs plotted')
            except Exception as e:
                # The common error that occurs is Unexpected Error: unable to find ring point
                # We've given up trying to catch this is a more elegant way
                # The class of exception raised is Exception in shapefile.py
                error_count += 1
                continue

        if self.args.points:
            # _add_and_make_ref_reef_poly will only return num_indi_points when self.args.points
            num_indi_points = self._make_and_add_ref_reef_patches(coords=(point_coords_x, point_coords_y))

        print(f'\n{error_count} error producing records were discounted from the reference reefs')
        if self.args.points:
            print(f'{reef_count} reference reefs were added to the plot as {num_indi_points} individual points')
        else:
            print(f'{reef_count} reference reefs were added to the plot')

    def _make_and_add_ref_reef_patches(self, coords, sizes=None):
        """
        When self.args.points we will plot the reference reefs as cirlces on the map.
        We were originally doing this using Circle patches and then either adding them to the plot
        as individual patches or as part of a PatchCollection. However, this is very slow when
        working with larger numbers of reefs. It is slow to add the patches to the ax, but it is
        also very slow to write out the figure as .png and .svg

        It is much faster to use scatter to plot the points. This also leads to a much faster write speed for the
        figure.
        """
        if self.args.points:
            # For the size of the scatter
            coords_x = np.concatenate(coords[0])
            coords_y = np.concatenate(coords[1])
            if sizes:
                # Then sizes have been provided and we should use these
                # This is when we are plotting user specified reefs
                points_size_array_sqr = [(size*self.coord_to_point_scaler)**2 for size in sizes]
                self.large_map_ax.scatter(x=coords_x, y=coords_y, s=points_size_array_sqr, zorder=2, facecolors=self.config_dict['reference_reef_color'], edgecolors='none')
            else:
                # Then we are plotting the user reefs and we should work with a standard
                # size. A sensible size is perhaps 1/500th of the shortest size
                lat = self.bounds[3] - self.bounds[2]
                lon = self.bounds[1] - self.bounds[0]
                if lat > lon:
                    # Then we should work with 1/100 of the lon range
                    deg_size = lon/500
                else:
                    deg_size = lat/500
                point_size = self.coord_to_point_scaler * deg_size
                self.large_map_ax.scatter(x=coords_x, y=coords_y, s=point_size**2, zorder=2, facecolors=self.config_dict['reference_reef_color'], edgecolors='none')
            return len(coords_x)
        else:
            reef_poly = Polygon([(x, y) for x, y in zip(list(coords[0]), list(coords[1]))],
                                closed=True, fill=True, edgecolor=self.config_dict['reference_reef_color'], linewidth=1,
                                facecolor=self.config_dict['reference_reef_color'],
                                alpha=1, zorder=2)
            self.large_map_ax.add_patch(reef_poly)

    def _if_within_bounds(self, bounds):
        """
        Check bounds for a given polygon lie within the bounds for the map
        """
        if bounds[0] > self.bounds[0]:
            if bounds[1] > self.bounds[2]:
                if bounds[2] < self.bounds[1]:
                    if bounds[3] < self.bounds[3]:
                        return True
        return False
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
        """
        A note on how to plot gridlines using cartopy:
        Although there is a GeoAxis.gridlines() method, this method does not yet allow a lot of
        bespoke options. If we want to only put the labels on the top and left then we have to
        generate a Gridliner object (normally returned by GeoAxis.gridlines() ourselves. We then need
        to manually change the xlabels_bottom and ylabels_right attributes of this Gridliner object.
        We then draw it by adding it to the GeoAxis._gridliners list.
        """
        if self.config_dict['grid_line_x_position'] and self.config_dict['grid_line_y_position']:
            xlocs = mticker.FixedLocator([float(_) for _ in self.config_dict['grid_line_x_position'].split(',')])
            ylocs = mticker.FixedLocator([float(_) for _ in self.config_dict['grid_line_y_position'].split(',')])
            g1 = Gridliner(
                axes=self.large_map_ax, crs=ccrs.PlateCarree(), draw_labels=True,
                xlocator=xlocs, ylocator=ylocs)
        else:
            g1 = Gridliner(
                axes=self.large_map_ax, crs=ccrs.PlateCarree(), draw_labels=True)
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

