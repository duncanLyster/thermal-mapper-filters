''' 
This script is built for performing filter calculations for the Enceladus Thermal Mapper (ETM) instrument.
It reads in mission configuration info from a .json file, filter assembly configuration info from a .json file, and a .csv file containing the filter transmission data.
It then calculates the filter transmission for each filter assembly, and signal-to-noise ratio for each filter assembly in a range of wavelengths/scenarios. 

To Do: 
Read names from .csv file instead of .json. 
Make .json more readable.
'''

import numpy as np
import json
import pandas as pd
import matplotlib.pyplot as plt
    

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
    load_and_assign_model_parameters("instrument.json")

    # Load filter assemblies
    assemblies = load_filter_assemblies("assemblies.json")

    temperatures = [30, 45 ,60, 100, 140, 180, 200]
    scene_emissivity = 0.9

    # Read the first two rows separately for filter names and numbers
    with open("filter_transmission.csv", 'r') as f:
        filter_names = next(f).strip().split(',')[1:]  # Skip the first entry (wavelength)
        filter_numbers = next(f).strip().split(',')[1:]  # Skip the first entry (filter number)

    # Load the rest of the CSV directly into a pandas DataFrame for easier NaN handling
    df = pd.read_csv("filter_transmission.csv", skiprows=2)
    
    # Check for NaN values and fill them with 0
    df.fillna(0, inplace=True) 
    
    # Convert the DataFrame back into a NumPy array if necessary
    csv_data = df.to_numpy()
    wavelength_array = csv_data[:, 0]  # Wavelengths
    transmission_data = csv_data[:, 1:]  # Transmission data
    

    # Convert filter numbers to integers
    filter_numbers = [int(num) for num in filter_numbers]

    for assembly in assemblies:
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

        snr_messages = []

        plt.figure(figsize=(10, 8))

        snr_table = []
        for filter_obj in filter_assembly:
            # Skip adding SNR information for filters where tdi_pixels is 1
            if filter_obj.tdi_pixels == 1:
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
            linestyle = '--' if filter_obj.tdi_pixels == 1 else '-'
            mask = filter_obj.transmission > 0 if filter_obj.tdi_pixels == 1 else slice(None)
            plt.plot(wavelength_array[mask], filter_obj.transmission[mask], label=f"{filter_obj.name}", linestyle=linestyle)

        # Generate the temperature label string
        temperature_labels = ', '.join([f"{temp}K" for temp in temperatures]) + "\nnormalised blackbody curves"

        # Plot the temperatures with a single label
        for index, temp in enumerate(temperatures):
            scene_radiance = calculate_scene_radiance(temp, scene_emissivity, wavelength_array=wavelength_array)
            normalised_radiance = scene_radiance / np.max(scene_radiance)
            
            if index == 0:
                # Only the first plot gets the consolidated label
                plt.plot(wavelength_array, normalised_radiance, 'lightgrey', linestyle=':', label=temperature_labels)
            else:
                # Other plots do not get a label
                plt.plot(wavelength_array, normalised_radiance, 'lightgrey', linestyle=':')

        plt.xlabel("Wavelength (um)")
        plt.ylabel("Radiance / Transmission")
        plt.xscale('log')
        plt.title(f"{assembly['assembly_name']} Filter Assembly")
        plt.subplots_adjust(bottom=0.2)  # Adjust the bottom margin
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3)
        plt.show()

# Call the main program to start execution
if __name__ == "__main__":
    main()