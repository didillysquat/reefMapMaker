# reefMapMaker
## Introduction
reefMapMaker produces a map figure (.png and .svg format) of a user specified region
that is annotated with a reference set of coral reef 
locations as well as user-specified points of interest.
The locations of the reference reefs are imported from the [Global Distribution of Coral 
Reefs data set](https://data.unep-wcmc.org/datasets/1).

This dataset is compiled from a number of sources by UNEP World Conservation Monitoring Centre
(UNEP-WCMC) and the WorldFish Centre, in collaboration with WRI (World Resources Institute) and
TNC (The Nature Conservancy).

reefMapMaker is written in Python and leverages 
the [cartopy](https://scitools.org.uk/cartopy/docs/latest/) and 
[matplotlib](https://matplotlib.org/) libraries to produce its maps.

reefMapMaker is designed to simplify the task of producing publication-ready map figures that
detail the locations of user-specified points of interest. Making such maps is often undertaken
using R- or Python-based packages/scripts. However, doing so can often be time consuming and require
considerable experience/familiarity with the mapping libraries in question.
 
reefMapmaker is designed to enable users with 0 experience in scripting languages to produce publication-ready
maps of reef/site locations such as the one shown below. The output of reefMapMaker is a set of .png and .svg files.
The .svg file can be imported into a vector-based graphics editor of the user's choice for further manipulation.

An example output:

![example_map_red_sea](./map_out_20201215T110504.svg) 

## Installation

Installation of reefMapMaker can be considered in two parts
1. Installation of the reefMapMaker code and dependencies
2. Installation of the reference reef shapefile and meta information

### 1. Installation of the reefMapMaker code and dependencies
This is easiest done using [conda](https://docs.conda.io/projects/conda/en/latest/).

#### For those familiar with conda:
`conda install -c didillysquat -c conda-forge reefMapMaker`

#### For those unfamiliar with conda
##### Install conda
conda is a package and environment manager program.

It is most commonly installed as part of Anaconda or Miniconda.

Downloading and installing either of Anaconda or Miniconda will give you access to conda
and allow you to install reefMapMaker.

Anaconda comes with many scientific packages preinstalled and is therefore much larger than Miniconda
that comes with only conda and Python and the packages they depend on.

See the [docs](https://docs.anaconda.com/anaconda/install/) for installation.

##### Create a new environment and install reefMapMaker
It is advisable to install reefMapMaker in a new conda environment.

A new environment called reefMapMaker_evnv (change this name to whatever you like)
can be created with reefMapMaker installed using the following single command:

`conda create --name reefMapMaker_env -c didillysquat -c coda-forge reefMapMaker`

### 2. Installation of the reference reef shapefile and meta information
Due to [license/use restrictions](https://www.unep-wcmc.org/policies/general-data-license-excluding-wdpa#data_policy),
the shapefile that contains the data for the reference
reefs cannot be included as part of this package, nor can it be automatically downloaded by this package.

To install the reference reef dataset, download the dataset from here:
https://data.unep-wcmc.org/datasets/1

N.B. If the download does not start when you click on
the download button, right click, then 'open in new tab'.

Once downloaded, decompress the .zip file and place the entire directory
(currently '14_001_WCMC008_CoralReefs2018_v4') someplace safe and note the location.
Do not change the names of the downloaded files.

By default, reefMapMaker will look for the reference reef datafiles in your current working
directory. Alternatively you can supply the path to the directory using the
--ref-reef-path.

## Usage
### Basic usage
reefMapMaker can be run with no inputs:

`python3 reefMapMaker`

The output will be a global map with the reef locations plotted using default parameters:

