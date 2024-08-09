# Installation 
1. Make sure you have a recent version of `conda` installed (we recommend [miniforge](https://github.com/conda-forge/miniforge)). To test this, open your terminal/command line/conda prompt and try to run the command `conda`.
2. [Download this file](https://github.com/serpentine-h2020/serpentine/archive/refs/heads/main.zip) and extract to a folder of your choice (or clone the repository https://github.com/serpentine-h2020/serpentine if you know how to use `git`).
3. Open your terminal/command line/conda prompt, navigate to the downloaded (extracted) folder `notebooks/sep_analysis_tools` that contains the file `requirements.txt`, and run the following:

    ``` bash
    $ conda create --name serpentine python=3.11
    $ conda activate serpentine
    $ pip install -r requirements.txt
    ```


# Run 
1. Open your terminal/command line/conda prompt.
2. In the terminal, navigate to the downloaded (extracted) folder `notebooks/sep_analysis_tools` that contains some `.ipynb` files.
3. Make sure the corresponding conda environment is activated by running `conda activate serpentine` in the terminal.
4. Run `jupyter notebook`
5. Your standard web-browser should now open the Jupyter interface, where you can double click on the corresponding `.ipynb` files to launch them.
