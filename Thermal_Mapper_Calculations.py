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
latitude_deg = None
longitude_deg = None
altitude_m = None
pointing_azimuth_deg = None
pointing_elevation_deg = None
roll_deg = None
target_mass_kg = None
target_radius_m = None

class Filter:
    def __init__(self, index, transmission_data):
        self.index = index
        self.transmission_data = transmission_data
        self.transmission = None
        self.signal_to_noise = None

def load_and_assign_model_parameters(json_filepath):
    with open(json_filepath, "r") as file:
        parameters = json.load(file)

    # Assign each parameter as a global variable
    for key, value in parameters.items():
        globals()[key] = value

def main():
    # Load the model parameters
    load_and_assign_model_parameters("instrument.json")

    # Load the filter assembly configuration from filter_assembly.json, assign to Filter objects
    filter_assembly = []
    with open("test_assembly.json", "r") as file:
        filter_assembly_config = json.load(file)
        for i, filter_data in enumerate(filter_assembly_config["filters"]):
            filter_assembly.append(Filter(i + 1, filter_data))

    # Print out the number of filters in the filter assembly
    print(f"Number of filters in filter assembly: {len(filter_assembly)}")

    # Preload the filter transmission data once
    df = pd.read_csv("filter_transmission.csv", skiprows=1)

    # Load the filter transmission data from the file filter_transmission.csv and assign to Filter objects
    for filter in filter_assembly:
        filter.transmission = df.iloc[:, filter.index].values

    # Plot the filter transmission data
    for filter in filter_assembly:
        plt.plot(df.iloc[:, 0], filter.transmission, label=f"Filter {filter.index}")

    # plt.xlabel("Wavelength (nm)")
    # plt.ylabel("Transmission")
    # plt.legend()
    # plt.show()


# Call the main program to start execution
if __name__ == "__main__":
    main()