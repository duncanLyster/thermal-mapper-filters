''' 
This script is built for performing filter calculations for the Enceladus Thermal Mapper (ETM) instrument.
It reads in mission configuration info from a .json file, filter assembly configuration info from a .json file, and a .csv file containing the filter transmission data.
It then calculates the filter transmission for each filter assembly, and signal-to-noise ratio for each filter assembly in a range of wavelengths/scenarios. 
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

def main():
    # Load the model parameters
    load_and_assign_model_parameters("instrument.json")

    # Load filter assemblies
    assemblies = load_filter_assemblies("assemblies.json")

    scene_temperature = 180
    scene_emissivity = 0.7

    # Preload the filter transmission data once
    df = pd.read_csv("filter_transmission.csv", skiprows=1)
    wavelength_array = df.iloc[:, 0].values


    # Load the filter assembly configuration from filter_assembly.json, assign to Filter objects
    filter_assembly = []

    for assembly in assemblies:
        filter_assembly = []

        print(f"\n{assembly['assembly_name']} at {scene_temperature}K\n{'-'*50}")

        # Create Filter objects for each filter in the current assembly
        for filter_data in assembly["filters"]:
            filter_obj = Filter(filter_data["filter_number"], filter_data)
            filter_obj.tdi_pixels = filter_data["tdi_pixels"]
            filter_obj.name = filter_data["filter_name"]
            filter_obj.transmission = df.iloc[:, filter_obj.number].values / 100
            filter_assembly.append(filter_obj)

        # Calculate SNRs for the current assembly
        scene_radiance = calculate_scene_radiance(scene_temperature, scene_emissivity, wavelength_array=wavelength_array)
        
        pixel_integration_time = 0.6  # Example value
        overall_system_transmission = mirror_reflectivity**number_of_mirrors * detector_absorption
        SNR_factor = Dstar * detector_side_length * solid_angle * overall_system_transmission

        for filter_obj in filter_assembly:
            filter_obj.integrated_transmission = np.trapz(filter_obj.transmission * scene_radiance, wavelength_array)
            filter_obj.signal_to_noise = SNR_factor * (pixel_integration_time * filter_obj.tdi_pixels)**0.5 * filter_obj.integrated_transmission
            print(f"Filter {filter_obj.number:<4} {filter_obj.name:<20} {filter_obj.tdi_pixels:<3} px wide   SNR: {filter_obj.signal_to_noise:<.2f}")

        # Plot the filter transmission for the current assembly with the normalised blackbody spectrum at the scene temperature overlaid
        # Normalise the blackbody spectrum to the maximum transmission value
        scene_radiance = scene_radiance / np.max(scene_radiance)
        plt.figure()
        plt.plot(wavelength_array, scene_radiance, label="Blackbody Radiance at 180K (normalised)")
        for filter_obj in filter_assembly:
            plt.plot(wavelength_array, filter_obj.transmission, label=f"Filter {filter_obj.number} ({filter_obj.name})")
        plt.xlabel("Wavelength (um)")
        plt.ylabel("Radiance / Transmission")
        plt.title(f"{assembly['assembly_name']} Filter Transmission")
        plt.legend()
        plt.show()

# Call the main program to start execution
if __name__ == "__main__":
    main()