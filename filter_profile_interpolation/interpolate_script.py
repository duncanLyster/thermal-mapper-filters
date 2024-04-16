'''
This script loads spectral data at a lower resolution than required for Thermal_Mapper_Calculations.py from a CSV file, performs 
linear interpolation to adjust the wavelength range, and ensures all values are within specified limits. It then corrects any interpolated 
emissivity values that are negative to zero, and saves the new interpolated dataset to another CSV file. The key operations include:
- Loading data from some 'spectrum.csv' - change file name as required.
- Adjusting wavelength bounds to fit within the 5 to 250 range. (This can be changed as needed.)
- Interpolating emissivities over a finely spaced wavelength grid of 0.01 increments.
- Setting negative emissivities to zero after interpolation.
- Saving the resulting data to 'interpolated_spectral_data.csv'. (Change file name as needed.)

Author: Duncan Lyster
'''

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

# Load the data from CSV
data = pd.read_csv('NaCl_spectrum.csv', header=None, names=['Wavelength', 'Emissivity'])

# Determine the bounds for the new wavelength range based on your data
min_wavelength = data['Wavelength'].min()
max_wavelength = data['Wavelength'].max()

# Adjust the bounds to ensure they are within your specified limits
if min_wavelength < 5:
    min_wavelength = 5
if max_wavelength > 250:
    max_wavelength = 250

# Define the new wavelength range based on the data's min and max, maintaining the 0.01 step
new_wavelengths = np.arange(min_wavelength, max_wavelength + 0.01, 0.01)

# Prepare the interpolation function
interpolator = interp1d(data['Wavelength'], data['Emissivity'], kind='linear', bounds_error=False, fill_value="extrapolate")

# Interpolate data
new_emissivities = interpolator(new_wavelengths)

# Set any negative values to zero
new_emissivities[new_emissivities < 0] = 0

# Save the interpolated data to a new CSV file
interpolated_data = pd.DataFrame({
    'Wavelength': new_wavelengths,
    'Emissivity': new_emissivities
})
interpolated_data.to_csv('interpolated_spectral_data.csv', index=False)
