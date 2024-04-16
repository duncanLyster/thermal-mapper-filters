import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize


def fwhm_calc(x, y, peak):
    """
    Calculate the Full Width at Half Maximum (FWHM) of a curve given by (x, y).
    """
    half_max = peak / 2
    left_idx = np.argwhere(y > half_max).min()
    right_idx = np.argwhere(y > half_max).max()
    fwhm = x[right_idx] - x[left_idx]
    return fwhm

def objective(sigma, x, mean, fwhm_target, peak, n):
    """
    Objective function for optimization: difference between calculated and target FWHM.
    """
    y = super_gaussian_params(x, mean, sigma, peak, n, optimize=False)
    fwhm_calculated = fwhm_calc(x, y, peak)
    return abs(fwhm_calculated - fwhm_target)

def super_gaussian_params(x, mean, sigma, peak, n=1, optimize=True):
    """
    Compute the super-Gaussian function with specified mean, FWHM, peak value, and n.
    Uses 'sigma' directly if 'optimize' is False, indicating that this is an optimization step.
    """
    if optimize:
        # Only calculate sigma based on FWHM if optimizing
        initial_sigma = sigma  # Now 'sigma' is directly passed as the initial guess for optimization
        # Minimize the objective function to find the optimized sigma
        result = minimize(objective, initial_sigma, args=(x, mean, fwhm_input, peak, n), method='Nelder-Mead')
        if result.success:
            optimized_sigma = result.x[0]
        else:
            optimized_sigma = initial_sigma  # Fallback to initial guess if optimization fails
        # Recalculate y using the optimized sigma
        y = peak * np.exp(-((np.abs(x - mean) / optimized_sigma) ** n))
    else:
        # Use sigma directly for calculation if not in optimization step
        y = peak * np.exp(-((np.abs(x - mean) / sigma) ** n))
    
    return y


# Inputs
mean_input = 250  # µm
fwhm_input = 100  # µm
peak_input = 0.9 # Peak value of the super-Gaussian
n_input = 8 # Higher number for more square-like shape

# Define x-axis range
x_range = np.arange(5, 250, 0.01)

# Calculate super-Gaussian with the new parameters
y_range = super_gaussian_params(x_range, mean_input, fwhm_input, peak_input, n_input)

# Save the data to CSV
data_to_save = pd.DataFrame({'Wavelength (µm)': x_range, 'Transmission': y_range})
csv_path = 'super_gaussian_curve.csv'
data_to_save.to_csv(csv_path, index=False)

# Plot
plt.figure(figsize=(10, 5))
plt.plot(x_range, y_range, label=f'Super-Gaussian, Peak: {peak_input}', color='green')
plt.title('Super-Gaussian Curve with Optimized Parameters')
plt.xlabel('Wavelength (µm)')
plt.ylabel('Transmission')
plt.grid(True)
plt.legend()
plt.show()

csv_path