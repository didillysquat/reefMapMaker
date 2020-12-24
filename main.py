#!/usr/bin/env python3
"""
This is a script for creating maps with reef locations
It is a modification from the CBASS84 script I used to make the figure for that paper
We will aim for it to take command inputs of the map bounding latitudes and longitudes
By default we will plot the world reefs on the map
We will also want to be able to add additional reefs to the map through some sort of csv import

# TODO finish coding up verification of the config parameters
# make args for the config parameters
# Finish documentation
# produce examples
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
from matplotlib.colors import is_color_like
import re

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
        self.config_params = ['bounds', 'plot_sea', 'sea_color', 'plot_reference_reefs',
                'reference_reef_color', 'reference_reef_patch_type', 'plot_land', 'land_color',
                'plot_grid_lines', 'x_grid_line_pos', 'y_grid_line_pos',
                'x_grid_lab_pos', 'y_grid_lab_pos', 'plot_boundaries']
        self.param_defaults_dict = {
            'bounds': [-180, 180, -90, 90],
            'plot_sea': True,
            'sea_color': '#88b5e0',
            'plot_reference_reefs': True,
            'reference_reef_color': '#003366',
            'plot_land': True, 'land_color': 'white',
            'plot_grid_lines': True, 'x_grid_line_pos': None,
            'y_grid_line_pos': None, 'x_grid_lab_pos': 'bottom',
            'y_grid_lab_pos': 'left', 'plot_boundaries': True
        }
        if self.args.config_sheet:
            # Config sheet
            print('Checking user config...')
            if self.args.config_sheet.endswith('.tsv'):
                config_df = pd.read_csv(self.args.config_sheet, index_col=0, sep='\t')
            elif self.args.config_sheet.endswith('.xlsx'):
                config_df = pd.read_excel(self.args.config_sheet, sheet_name='user_config', index_col=0)
            self.config_dict = {k: v for k, v in zip(config_df.index.values, config_df['value'].values)}
            self._setup_config()
            print('User config checks complete.')

            if self.args.site_sheet:
                print('Reading in site sheet')
                if self.args.site_sheet.endswith('.tsv'):
                    self.site_df = pd.read_csv(self.args.config_sheet, index_col=0, sep='\t').sort_values('radius_in_deg_lat', axis=0, ascending=False)
                elif self.args.site_sheet.endswith('.xlsx'):
                    self.site_df = pd.read_excel(
                        self.args.site_sheet, sheet_name='user_site', index_col=0
                    ).sort_values('radius_in_deg_lat', axis=0, ascending=False)
                print('Checking user sites...')
                self._check_site_sheet()
                print('User site checks complete.')
            else:
                print('No user site sheet provided. User sites will not be plotted.')
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
                'plot_grid_lines': True, 'x_grid_line_pos': None, 'y_grid_line_pos': None,
                'x_grid_lab_pos': 'bottom', 'y_grid_lab_pos': 'left', 'plot_boundaries': True
            })
            self.user_site_dict = {}
            self.site_df = None

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

    def _check_site_sheet(self):
        """
        Check the site_sheet input for valid inputs
        """
        for site in self.site_df.index:
            self._check_site_lat_lon(site)

    def _check_site_lat_lon(self, site):
        lat = self.site_df.at[site, 'latitude_deg_n']
        lon = self.site_df.at[site, 'longitude_deg_e']
        if pd.isnull(lat) or pd.isnull(lon):
            raise RuntimeError(f"invalid lat or lon format for site {site}")
        try:
            self.site_df.at[site, 'latitude_deg_n'] = float(lat)
            self.site_df.at[site, 'longitude_deg_e'] = float(lon)
        except ValueError:
            print("WARNING: Unable to convert lat or lon value for {site} to decimal degrees.")
            print("Attempting to fix the format")
            try:
                if 'N' in lat:
                    new_lat = float(lat.replace('N', '').replace(chr(176), ''))
                    # lat_float should be positive
                    if new_lat < 0:
                        new_lat = new_lat * -1
                elif 'S' in lat:
                    new_lat = float(lat.replace('S', '').replace(chr(176), ''))
                    # lat_float should be negative
                    if new_lat > 0:
                        new_lat = new_lat * -1
                else:
                    # There was not an N or S found in the lat so we should raise error
                    raise RuntimeError(f'Unable to fix lat lon format for {site}')
                if 'E' in lon:
                    new_lon = float(lon.replace('E', '').replace(chr(176), ''))
                    # lon_float should be positive
                    if new_lon < 0:
                        new_lon = new_lon * -1
                elif 'W' in lon:
                    new_lon = float(lon.replace('W', '').replace(chr(176), ''))
                    # lon_float should be negative
                    if new_lon > 0:
                        new_lon = new_lon * -1
                else:
                    # There was not an N or S found in the lat so we should raise error
                    raise RuntimeError
            except Exception:
                # see if they are in proper dms format
                try:
                    new_lat = self.dms2dec(lat)
                    new_lon = self.dms2dec(lon)
                # if all this fails, convert to 999
                except Exception:
                    raise RuntimeError(f'Unable to fix lat lon format for {site}')
            print(f"For {site}:\n"
                  f"\tlatitude was converted from {lat} --> {new_lat}\n"
                  f"\tlongitude was converted from {lon} --> {new_lon}")

    @staticmethod
    def dms2dec(dms_str):
        """Return decimal representation of DMS

            dms2dec(utf8(48°53'10.18"N))
            48.8866111111F

            dms2dec(utf8(2°20'35.09"E))
            2.34330555556F

            dms2dec(utf8(48°53'10.18"S))
            -48.8866111111F

            dms2dec(utf8(2°20'35.09"W))
            -2.34330555556F

            """

        dms_str = re.sub(r'\s', '', dms_str)

        sign = -1 if re.search('[swSW]', dms_str) else 1

        numbers = [*filter(len, re.split('\D+', dms_str, maxsplit=4))]

        degree = numbers[0]
        minute = numbers[1] if len(numbers) >= 2 else '0'
        second = numbers[2] if len(numbers) >= 3 else '0'
        frac_seconds = numbers[3] if len(numbers) >= 4 else '0'

        second += "." + frac_seconds
        return sign * (int(degree) + float(minute) / 60 + float(second) / 3600)

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

        Config params are:
        'bounds', 'plot_sea', 'sea_color', 'plot_reference_reefs',
                'reference_reef_color', 'reference_reef_patch_type', 'plot_land', 'land_color',
                'plot_grid_lines', 'x_grid_line_pos', 'y_grid_line_pos',
                'x_grid_lab_pos', 'y_grid_lab_pos', 'plot_boundaries'
        """
        self._check_bounds()
        self._check_bool_params()
        self._check_color_params()
        self._check_grid_lab_pos()
        self._check_grid_line_coords()
        self._check_ref_reef_path_type()

    def _check_ref_reef_path_type(self):
        """
        Check user provided valid value for reference_reef_patch_type or set default
        """
        if self.config_dict['plot_reference_reefs']:
            cl_param_set, config_param_set = self._param_set(param='reference_reef_patch_type')
            self._set_config_param(param='reference_reef_patch_type', cl_param=cl_param_set, config_param=config_param_set)
            if not self.config_dict['reference_reef_patch_type'] in ['polygon', 'point']:
                self._set_default_param(param='reference_reef_patch_type')

    def _check_grid_line_coords(self):
        for g_line_cords_param in ['x_grid_line_pos', 'y_grid_line_pos']:
            cl_param_set, config_param_set = self._param_set(param=g_line_cords_param)
            self._set_config_param(param=g_line_cords_param, cl_param=cl_param_set, config_param=config_param_set)
        self._check_valid_grid_line_coords()

    def _check_valid_grid_line_coords(self):
        if self.config_dict['plot_grid_lines']:
            try:
                assert (len(self.config_dict['x_grid_line_pos'].split(',')) > 0)
                [float(_) for _ in self.config_dict['x_grid_line_pos'].split(',')]
                assert (len(self.config_dict['y_grid_line_pos'].split(',')) > 0)
                [float(_) for _ in self.config_dict['y_grid_line_pos'].split(',')]
            except AssertionError:
                print("WARNING: There is a problem with the formatting of your grid line coordinates.")
                print('Default coordinates will be used.')
                self.config_dict['x_grid_line_pos'] = None
                self.config_dict['y_grid_line_pos'] = None

    def _check_grid_lab_pos(self):
        """
        Check for a user input for the grid_lab_pos params and set to default if not.
        """
        for g_lab_pos_param in ['x_grid_lab_pos', 'y_grid_lab_pos']:
            cl_param_set, config_param_set = self._param_set(param=g_lab_pos_param)
            self._set_config_param(param=g_lab_pos_param, cl_param=cl_param_set, config_param=config_param_set)
        self._check_valid_grid_lab_pos()

    def _check_valid_grid_lab_pos(self):
        try:
            assert (self.config_dict['x_grid_lab_pos'] in ['bottom', 'top'])
        except AssertionError:
            raise RuntimeError("x_grid_lab_pos must be 'bottom' or 'top'")
        try:
            assert (self.config_dict['y_grid_lab_pos'] in ['left', 'right'])
        except AssertionError:
            raise RuntimeError("y_grid_lab_pos must be 'left' or 'right'")

    def _check_color_params(self):
        """
        Check the config and user site color parameters provided
        and set the config dict to the default value if necessary.
        """
        # Map config colour params
        for c_param in ['sea_color', 'land_color', 'reference_reef_color']:
            cl_param_set, config_param_set = self._param_set(param=c_param)
            self._set_config_param(param=c_param, cl_param=cl_param_set, config_param=config_param_set)
            if not is_color_like(self.config_dict[c_param]):
                self._set_default_param(param=c_param)
        # User site colour params
        # TODO check the site params


    def _check_bool_params(self):
        for bool_param in ['plot_sea', 'plot_reference_reefs', 'plot_land', 'plot_grid_lines', 'plot_boundaries']:
            cl_param_set, config_param_set = self._param_set(param=bool_param)
            self._set_config_param(param=bool_param, cl_param=cl_param_set, config_param=config_param_set)
            self._check_valid_bool_param(bool_param)

    def _check_valid_bool_param(self, bool_param):
        try:
            assert (type(self.config_dict[bool_param]) is bool)
        except AssertionError:
            raise RuntimeError(f"{bool_param} must be either TRUE or FALSE.")

    def _check_bounds(self):
        # Check to see if the bounds are set by either the config_sheet or the command line
        cl_param_set, config_param_set = self._param_set(param='bounds')
        # Then there is user supplied value for bounds
        self._set_config_param(param='bounds', cl_param=cl_param_set, config_param=config_param_set)
        self._check_bounds_in_correct_order()
        self._check_bounds_valid_values()

    def _check_bounds_valid_values(self):
        try:
            assert ((self.config_dict['bounds'][0] < 180) and (self.config_dict['bounds'][0] > -180))
            assert ((self.config_dict['bounds'][1] < 180) and (self.config_dict['bounds'][1] > -180))
        except AssertionError:
            raise RuntimeError("One of your longitude bounds appears to be an invalid value.")
        try:
            assert ((self.config_dict['bounds'][2] < 90) and (self.config_dict['bounds'][2] > -90))
            assert ((self.config_dict['bounds'][3] < 90) and (self.config_dict['bounds'][3] > -90))
        except AssertionError:
            raise RuntimeError("One of your latitude bounds appears to be an invalid value.")

    def _check_bounds_in_correct_order(self):
        try:
            assert (self.config_dict['bounds'][0] < self.config_dict['bounds'][1])
        except AssertionError:
            raise RuntimeError('Check the format of your user config sheet.\n'
                               'latitude and longitude valus must be in valid decimal degree format and\n'
                               'westernmost bound should be less than easternmost bound')
        try:
            assert (self.config_dict['bounds'][2] < self.config_dict['bounds'][3])
        except AssertionError:
            raise RuntimeError('Check the format of your user config sheet.\n'
                                'latitude and longitude valus must be in valid decimal degree format and\n'
                               'southernmost bounds should be less than northernmost bound')

    def _set_config_param(self, param, cl_param, config_param):
        if config_param and cl_param:
            print(
                f"WARNING: the {param} parameter has been supplied both in the user "
                f"config sheet ({self.config_dict[param]}) and on the command line ({getattr(self.args, param)}).\n"
                f"The command line-supplied argument will be used."
            )
            if param == 'bounds':
                try:
                    self.config_dict[param] = [float(_) for _ in getattr(self.args, param).split(',')]
                except ValueError:
                    raise RuntimeError("Unable to convert one of the bounds to decimal degree format.\n"
                                       "Please check the formatting of your bounds.")
            else:
                self.config_dict[param] = getattr(self.args, param)
            self._notify_user_set_config_dict_param(param=param)
        elif config_param:
            # Then param only supplied by config_sheet
            # Conduct the check using the config_bounds
            if param == 'bounds':
                try:
                    self.config_dict[param] = [float(_) for _ in self.config_dict[param].split(',')]
                except ValueError:
                    raise RuntimeError("Unable to convert one of the bounds to decimal degree format.\n"
                                       "Please check the formatting of your bounds.")
            else:
                # No need to update the config_dict. Value already correctly set
                pass
            self._notify_user_set_config_dict_param(param=param)
        elif cl_param:
            # Then param only supplied by command line
            if param == 'bounds':
                try:
                    self.config_dict[param] = [float(_) for _ in getattr(self.args, param).split(',')]
                except ValueError:
                    raise RuntimeError("Unable to convert one of the bounds to decimal degree format.\n"
                                       "Please check the formatting of your bounds.")
            else:
                self.config_dict[param] = getattr(self.args, param)
            self._notify_user_set_config_dict_param(param=param)
        else:
            # Then no valid argument set for the param
            # Set param from default dict
            self._set_default_param(param)

    def _set_default_param(self, param):
        print(f'No valid value for {param} provided. Using default.')
        self.config_dict[param] = self.param_defaults_dict[param]
        self._notify_user_set_config_dict_param(param=param)

    def _param_set(self, param):
        """
        Check if a given parameter had a value set via the config sheet or via the command line
        """
        config_param = False
        cl_param = False
        if param in self.config_dict:
            if not pd.isnull(self.config_dict[param]):
                config_param = True
        if getattr(self.args, param):
            cl_param = True
        return cl_param, config_param

    def _notify_user_set_config_dict_param(self, param):
        print(f'{param} set to: {self.config_dict[param]}')

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
        line_widths = ((self.site_df['radius_in_deg_lat'].astype(
                float) * 2) * self.coord_to_point_scaler) * 0.1
        self.large_map_ax.scatter(
            x=self.site_df['longitude_deg_e'],
            y=self.site_df['latitude_deg_n'],
            s=(((self.site_df['radius_in_deg_lat'].astype(
                float) * 2) * self.coord_to_point_scaler) ** 2),
            facecolors=self.site_df['facecolor'],
            edgecolors=self.site_df['edgecolor'], zorder=3,
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
        if self.config_dict['x_grid_line_pos'] and self.config_dict['y_grid_line_pos']:
            xlocs = mticker.FixedLocator([float(_) for _ in self.config_dict['x_grid_line_pos'].split(',')])
            ylocs = mticker.FixedLocator([float(_) for _ in self.config_dict['y_grid_line_pos'].split(',')])
            g1 = Gridliner(
                axes=self.large_map_ax, crs=ccrs.PlateCarree(), draw_labels=True,
                xlocator=xlocs, ylocator=ylocs)
        else:
            g1 = Gridliner(
                axes=self.large_map_ax, crs=ccrs.PlateCarree(), draw_labels=True)
        if self.config_dict['x_grid_lab_pos'] == 'bottom':
            g1.top_labels = False
        else:
            g1.bottom_labels = False
        if self.config_dict['y_grid_lab_pos'] == 'left':
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

