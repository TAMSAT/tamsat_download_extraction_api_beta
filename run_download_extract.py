"""
Download and extract TAMSAT rainfall estimates.

Author: TAMSAT (tamsat@reading.ac.uk)
"""

# Import the 'download' and 'extract' tools
from tamsat_download_extract_api import download, extract


# Specify local data directory
localdata_dir = '/gws/nopw/j04/tamsat/tamsat/scripts/api/data'

# Download TAMSAT files - copy and paste into terminal to run
download({
    "timestep": 'daily',
    "resolution": 0.25,
    "start_date": '2021-06-01',
    "end_date": '2021-06-15',
    "version": 3.1,
    "localdata_dir": localdata_dir
    })


# Extract TAMSAT rainfall estimates at a point - copy and paste into terminal to run
extract({
    "extract_type": 'point',
    "longitude": 22.73,
    "latitude": -3.51,
    "timestep": 'daily',
    "resolution": 0.25,
    "start_date": '2021-06-01',
    "end_date": '2021-06-15',
    "version": 3.1,
    "localdata_dir": localdata_dir
    })
    

# Extract area-average TAMSAT rainfall estimates - copy and paste into terminal to run
extract({
    "extract_type": 'area_average',
    "N": 22.73,
    "S": 20.73,
    "W": 4.0,
    "E": 5.0,
    "timestep": 'daily',
    "resolution": 0.25,
    "start_date": '2021-06-01',
    "end_date": '2021-06-15',
    "version": 3.1,
    "localdata_dir": localdata_dir
    })
    

# Subset TAMSAT rainfall estimates for a given domain - copy and paste into terminal to run
extract({
    "extract_type": 'domain',
    "N": 22.73,
    "S": 20.73,
    "W": 4.0,
    "E": 5.0,
    "timestep": 'daily',
    "resolution": 0.25,
    "start_date": '2021-06-01',
    "end_date": '2021-06-15',
    "version": '3.1',
    "localdata_dir": localdata_dir
    })
