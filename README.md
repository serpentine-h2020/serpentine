# SERPENTINE
[![DOI](https://zenodo.org/badge/498388422.svg)](https://zenodo.org/badge/latestdoi/498388422)


This repository contains the Python tools and software developed by the [solar energetic particle analysis platform for the inner heliosphere (SERPENTINE)](https://serpentine-h2020.eu/) project for the downloading of data and the performing of analysis and visualisation.

## Contents

Jupyter Notebooks:
- [Multi-Spacecraft Constellation Plotter Solar-MACH](https://github.com/serpentine-h2020/serpentine/tree/main/notebooks/solarmach)
- [Solar Energetic Particle Analysis Tools](https://github.com/serpentine-h2020/serpentine/tree/main/notebooks/sep_analysis_tools)

## Installation

1. Make sure you have a recent version of `conda` installed (we recommend [miniforge](https://github.com/conda-forge/miniforge)). To test this, open your terminal/command line/conda prompt and try to run the command `conda`.
2. [Download this file](https://github.com/serpentine-h2020/serpentine/archive/refs/heads/main.zip) and extract to a folder of your choice (or clone the repository https://github.com/serpentine-h2020/serpentine if you know how to use `git`).
3. Open your terminal/command line/conda prompt, navigate to the downloaded/extracted folder (which contains the file `requirements.txt`), and run the following:

    ``` bash
    $ conda create --name serpentine python=3.12
    $ conda activate serpentine
    $ pip install -r requirements.txt
    ```

## Usage

1. Open your terminal/command line/conda prompt.
2. In the terminal, navigate to the downloaded/extracted folder.
3. Make sure the corresponding conda environment is activated by running `conda activate serpentine` in the terminal.
4. Run `jupyter-lab'`, your standard web-browser should now open the JupyterLab interface.
5. In the *File Browser* (click *View* -> *File Browser* if it's not shown) double-click on the `notebooks` folder, then `sep_analysis_tools` or `solarmach`, and finally the corresponding `.ipynb` file for a specific tool.

## Other SERPENTINE Software

- [Solar Magnetic Connection Haus tool (Solar-MACH)](https://github.com/jgieseler/solarmach)
- [SEPpy - A compendium of Python data loaders and analysis tools for in-situ measurements of Solar Energetic Particles (SEP) ](https://github.com/serpentine-h2020/seppy)

## Citation

- **If you use the [Multi-Spacecraft Constellation Plotter Solar-MACH](https://github.com/serpentine-h2020/serpentine/tree/main/notebooks/solarmach) in your publication, please cite the following paper:**
 
  Gieseler, J., Dresing, N., Palmroos, C., von Forstner, J. L. F., Price, D. J., Vainio, R., Kouloumvakos A., Rodríguez-García L., Trotta D., Génot V., Masson A., Roth M., Veronig A. (2023).
Solar-MACH: An open-source tool to analyze solar magnetic connection configurations. *Front. Astronomy Space Sci.* 9. [doi:10.3389/fspas.2022.1058810](https://doi.org/10.3389/fspas.2022.1058810)
- **If you use the [Solar Energetic Particle Analysis Tools](https://github.com/serpentine-h2020/serpentine/tree/main/notebooks/sep_analysis_tools) in your publication, please cite the following paper:**

  Palmroos, C., Gieseler, J., Dresing N., Morosan D. E., Asvestari E., Yli-Laurila A., Price D. J., Valkila S., Vainio R. (2022). Solar energetic particle time series analysis with Python. *Front. Astronomy Space Sci.* 9. [doi:10.3389/fspas.2022.1073578](https://doi.org/10.3389/fspas.2022.1073578)

## Acknowledgements
<img align="left" height="45px" src="eu_logo.png"> 
This project has received funding from the European Union’s Horizon 2020 research and innovation programme under grant agreement No 101004159.
