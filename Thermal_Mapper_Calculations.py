'''
This script conducts optical filter analysis for thermal mapping instruments. It was built for the Enceladus Thermal Mapper proposal but can 
be used for any multi-spectral thermal instrument. 

It executes several key functions:

1. Reads mission configuration from 'instrument_parameters.json', detailing properties like detector characteristics
   and telescope specifications.
2. Loads filter assembly configurations from 'assemblies.json' and transmission data from 'ETM_filter_transmission_profile_data.csv'.
3. Computes the transmission efficiency for each filter in different assemblies across a specified range of wavelengths.
4. Calculates the signal-to-noise ratio (SNR) for each filter under multiple thermal scenarios, using the Planck radiation formula
   to simulate scene radiance at various temperatures.
5. Visualises the filter transmissions and SNR values across different setups, aiding in the assessment of filter performance
   in potential observational conditions.

Outputs include graphical representations of filter characteristics and a detailed SNR analysis for decision-making regarding
instrument design and mission planning.

Author: Duncan Lyster
'''

import numpy as np
import json
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, FuncFormatter
    

# Initalise global variables (to be assigned later by the model parameters file)
# Detector properties
wavenumber_min = None
wavenumber_max = None
wavenumber_unit = None
spectral_resolution = None
Dstar = None
detector_absorption = None
detector_fnumber = None
detector_side_length = None
telescope_diameter = None
fov_horizontal_deg = None
fov_vertical_deg = None
mirror_reflectivity = None
number_of_mirrors = None
solid_angle = None

# Orbit properties
latitude_deg = None
longitude_deg = None
altitude_m = None
pointing_azimuth_deg = None
pointing_elevation_deg = None
roll_deg = None
target_mass_kg = None
target_radius_m = None

class Filter:
    def __init__(self, number, transmission_data):
        self.number = number
        self.transmission_data = transmission_data
        self.transmission = None
        self.signal_to_noise = None
        self.tdi_pixels = None
        self.name = None

def tick_formatter(val, pos):
    return f"{val:g}"

def load_and_assign_model_parameters(json_filepath):
    with open(json_filepath, "r") as file:
        parameters = json.load(file)

    # Assign each parameter as a global variable
    for key, value in parameters.items():
        if isinstance(value, dict):
            for sub_key, sub_value in value.items():
                globals()[sub_key] = sub_value
        else:
            globals()[key] = value

def load_filter_assemblies(json_filepath):
    with open(json_filepath, "r") as file:
        assemblies = json.load(file)
    return assemblies

def calculate_scene_radiance(temperature, emissivity, wavelength_array):
    # Calculate the scene radiance for a given temperature, emissivity, and wavelength
    # Planck's constant
    h = 6.626e-34
    # Speed of light
    c = 3e8
    # Boltzmann's constant
    k = 1.38e-23

    # Calculate the scene radiance
    scene_radiance = []
    for i, wavelength in enumerate(wavelength_array):
    
        # Convert wavelength from um to m
        wavelength_m = wavelength * 1e-6
        
        # Calculate the wavenumber
        wavenumber = 1 / wavelength_m

        # Calculate the radiance
        radiance = emissivity * (2 * h * c**2 * wavenumber**3) / (np.exp(h * c * wavenumber / (k * temperature)) - 1)

        scene_radiance.append(radiance)

    return scene_radiance

def format_snr(value):
    # Check if the value is essentially an integer
    if value < 1.0:
        return f"{value:<10.3g}"
    else:
        return f"{value:<10.2f}"

def main():
    # Load the model parameters
    load_and_assign_model_parameters("instrument_parameters.json")

    # Load filter assemblies
    assemblies = load_filter_assemblies("assemblies.json")

    #temperatures = [30, 45, 60, 100, 140, 180, 200]
    temperatures = [30, 60, 180, 200] # Stripped back for graphs
    scene_emissivity = 0.9

    # Read the first two rows separately for filter names and numbers
    with open("ETM_filter_transmission_profile_data.csv", 'r') as f:
        filter_names = next(f).strip().split(',')[1:]  # Skip the first entry (wavelength)
        filter_numbers = next(f).strip().split(',')[1:]  # Skip the first entry (filter number)

    # Load the rest of the CSV into a pandas DataFrame
    df = pd.read_csv("ETM_filter_transmission_profile_data.csv", skiprows=2)
    df.fillna(0, inplace=True)  # Handle NaNs
    csv_data = df.to_numpy()
    wavelength_array = csv_data[:, 0]  # Wavelengths
    transmission_data = csv_data[:, 1:]  # Transmission data

    filter_numbers = [int(num) for num in filter_numbers]

    # Determine subplot grid size and setup subplots
    num_assemblies = len(assemblies)
    cols = 2
    rows = (num_assemblies + cols - 1) // cols  # Ensure enough rows
    fig, axs = plt.subplots(rows, cols, figsize=(5 * cols, 4 * rows))
    axs = axs.flatten()  # Flatten to simplify indexing

    # Define colors for each filter for consistency across plots
    colors = plt.cm.tab10(np.linspace(0, 1, len(filter_names)))  # Adjust as needed based on the number of filters

    # Dictionary to accumulate unique legend handles and labels across all subplots
    unique_legend_entries = {}

    for ax, assembly in zip(axs, assemblies):
        print(f"\n{assembly['assembly_name']}\n{'-'*50}")

        filter_assembly = []

        # Create Filter objects for each filter in the current assembly
        for filter_data in assembly["filters"]:
            filter_number = filter_data["filter_number"]
            filter_obj = Filter(filter_data["filter_number"], filter_data)
            filter_obj.tdi_pixels = filter_data["tdi_pixels"]
            filter_obj.name = filter_names[filter_number - 1] if (filter_number - 1) < len(filter_names) else "Unknown"
            filter_obj.transmission = transmission_data[:, filter_number - 1] / 100
            filter_assembly.append(filter_obj)

        pixel_integration_time = 0.6  # Calculated by hand for Nightingale
        overall_system_transmission = mirror_reflectivity**number_of_mirrors * detector_absorption
        SNR_factor = Dstar * detector_side_length * solid_angle * overall_system_transmission

        snr_table = []


        for filter_obj in filter_assembly:
            if filter_obj.tdi_pixels == 1: # Skip adding SNR information where tdi_pixels is 1 (e.g emission spectra)
                continue
            snr_row = [f"Filter {filter_obj.number:<4} {filter_obj.name:<35} Width (px): {filter_obj.tdi_pixels:<5}"]
            for temp in temperatures:
                scene_radiance = calculate_scene_radiance(temp, scene_emissivity, wavelength_array=wavelength_array)
                filter_obj.integrated_transmission = np.trapz(filter_obj.transmission * scene_radiance, wavelength_array)
                filter_obj.signal_to_noise = SNR_factor * (pixel_integration_time * filter_obj.tdi_pixels)**0.5 * filter_obj.integrated_transmission
                formatted_snr = format_snr(filter_obj.signal_to_noise)
                snr_row.append(f"{temp}K: {formatted_snr}")
            snr_table.append(snr_row)

        # Print SNR table
        for row in snr_table:
            print("{:<55} {}".format(row[0], " ".join(row[1:])))

        for filter_obj in filter_assembly:
            color = colors[filter_obj.number - 1] # Assign a color based on the filter number
            linestyle = '--' if filter_obj.tdi_pixels == 1 else '-'
            label = f"{filter_obj.name}"
            mask = filter_obj.transmission > 0 if filter_obj.tdi_pixels == 1 else slice(None)
            line, = ax.plot(wavelength_array[mask], filter_obj.transmission[mask], color=color, label=label, linestyle=linestyle)

            if label not in unique_legend_entries:  # Avoid duplicate legend entries
                unique_legend_entries[label] = line

        # Generate the temperature label string
        temperature_labels = ', '.join([f"{temp}K" for temp in temperatures]) + "\nnormalised blackbody curves"

        # Plot the temperatures with a single label
        for index, temp in enumerate(temperatures):
            scene_radiance = calculate_scene_radiance(temp, scene_emissivity, wavelength_array=wavelength_array)
            normalised_radiance = scene_radiance / np.max(scene_radiance)
            
            if index == 0:
                # Only the first plot gets the consolidated label
                ax.plot(wavelength_array, normalised_radiance, 'lightgrey', linestyle=':', label=temperature_labels)
            else:
                # Other plots do not get a label
                ax.plot(wavelength_array, normalised_radiance, 'lightgrey', linestyle=':')

        ax.set_xlabel("Wavelength (µm)", fontsize=14)
        ax.set_ylabel("Transmission", fontsize=14)
        ax.set_xscale('log')
        ax.tick_params(axis='both', which='major', labelsize=13)
        ax.xaxis.set_major_locator(LogLocator(base=10.0, subs=(1.0, 2.0, 5.0)))
        ax.xaxis.set_major_formatter(FuncFormatter(tick_formatter))
        ax.set_title(f"{assembly['assembly_name']} Filter Assembly", fontsize=15)
    
    fig.subplots_adjust(bottom=0.1)  # Adjust this value as needed to make space for the legend
    fig.legend(unique_legend_entries.values(), unique_legend_entries.keys(), loc='lower center', ncol=3, bbox_to_anchor=(0.5, 0.01), fontsize=12)
    plt.tight_layout(rect=[0, 0.1, 1, 0.95])  # Adjust the second value as needed based on the space required for the legend
    plt.show()

# Call the main program to start execution
if __name__ == "__main__":
    main()