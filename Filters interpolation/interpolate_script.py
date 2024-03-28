import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

# Load the data from CSV
data = pd.read_csv('unsorted.csv', header=None, names=['Wavelength', 'Emissivity'])

# Prepare the interpolation function
interpolator = interp1d(data['Wavelength'], data['Emissivity'], kind='linear', bounds_error=False, fill_value="extrapolate")

# Define the new wavelength range
new_wavelengths = np.arange(5, 250.01, 0.01)

# Interpolate data
new_emissivities = interpolator(new_wavelengths)

# Set any negative values to zero
new_emissivities[new_emissivities < 0] = 0

# Save the interpolated data to a new CSV file
interpolated_data = pd.DataFrame({
    'Wavelength': new_wavelengths,
    'Emissivity': new_emissivities
})
interpolated_data.to_csv('interpolated_data.csv', index=False)
