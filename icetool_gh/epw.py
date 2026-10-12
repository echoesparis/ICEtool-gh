"""Minimal EPW reader, mirroring upstream Step4 (:240-248)."""

import csv

import pandas as pd

EPW_NAMES = ['Year', 'month', 'day', 'hour', 'Minute', 'Data Source and Uncertainty Flags',
             'Dry Bulb Temperature [DegC]', 'Dew Point Temperature', 'Relative Humidity',
             'Atmospheric Station Pressure', 'Extraterrestrial Horizontal Radiation',
             'Extraterrestrial Direct Normal Radiation', 'Horizontal Infrared Radiation Intensity',
             'Global Horizontal Radiation [Wh/m2]', 'Direct Normal Radiation',
             'Diffuse Horizontal Radiation', 'Global Horizontal Illuminance',
             'Direct Normal Illuminance', 'Diffuse Horizontal Illuminance', 'Zenith Luminance',
             'Wind Direction', 'Wind Speed', 'Total Sky Cover', 'Opaque Sky Cover', 'Visibility',
             'Ceiling Height', 'Present Weather Observation', 'Present Weather Codes',
             'Precipitable Water', 'Aerosol Optical Depth', 'Snow Depth',
             'Days Since Last Snowfall', 'Albedo', 'Liquid Precipitation Depth',
             'Liquid Precipitation Quantity']


def read_epw(path):
    """Return the 8760-row hourly DataFrame with upstream's column names."""
    with open(path, newline='') as f:
        for i, row in enumerate(csv.reader(f, delimiter=',', quotechar='"')):
            if row and row[0].isdigit():
                break
    return pd.read_csv(path, skiprows=i, header=None, names=EPW_NAMES)


def read_epw_location(path):
    """Site metadata from the EPW LOCATION header line.

    Returns (latitude_deg, longitude_deg, tz_meridian_longitude_deg, elevation_m).
    """
    with open(path, newline='') as f:
        head = f.readline().split(',')
    return float(head[6]), float(head[7]), float(head[8]) * 15.0, float(head[9])
