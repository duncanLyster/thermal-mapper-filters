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
interpolated_data.to_csv('interpolated_Na_Cl_data.csv', index=False)
