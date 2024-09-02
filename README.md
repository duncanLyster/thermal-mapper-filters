# Multispectral Thermal Mapper Filters Toolkit

Welcome to the Multispectral Thermal Mapper Filters Toolkit! This repository houses tools and scripts designed to analyze optical filters for various multispectral thermal mapping instruments. Originally developed for the Enceladus Thermal Mapper proposal, the toolkit is adaptable for use with any multispectral thermal instrument.

## Contents
- `assemblies.json` - Configuration for filter assemblies.
- `ETM_filter_transmission_profile_data.csv` - Filter transmission data used in calculations.
- `Filters Profile Interpolation/` - Contains scripts and example files for filter profile interpolation including `interpolate_script.py`.
- `instrument_parameters.json` - Configuration file detailing the instrument parameters like detector characteristics and telescope specifications.
- `LICENSE` - The licensing file which details the terms under which this software can be used.
- `model_filter_profile_generator.py` - Main script for generating filter profiles.
- `README.md` - The file you are currently reading.
- `SNRs table.xlsx` - Excel spreadsheet containing calculated signal-to-noise ratios for different filter setups.
- `Thermal_Mapper_Calculations.py` - Script that performs detailed SNR calculations and generates visualizations of filter performances.
- `requirements.txt` - Required Python dependencies for the project.

## Installation
To get started with this toolkit, clone this repository to your local machine using:
git clone [repository URL]

Ensure you have Python installed, and then install required dependencies:
pip install -r requirements.txt

## Usage
To run the main filter profile generator:
python model_filter_profile_generator.py

For interpolation of filter profiles:
cd filter_profile_interpolation
python interpolate_script.py

To calculate SNRs and generate visual outputs:
python Thermal_Mapper_Calculations.py

## Contributing
Contributions to this toolkit are welcome. To contribute, please fork the repository, make your changes, and submit a pull request. For more details, see `CONTRIBUTING.md`.

## License
This project is licensed under the MIT License - see the `LICENSE` file for details.

## Authors
- Duncan Lyster - Initial work and ongoing maintenance.

For any additional information or support, please open an issue in this repository.