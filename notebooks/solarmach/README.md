# Installation 
1. Make sure you have a recent version of [Anaconda](https://www.anaconda.com/products/distribution) installed.
2. [Download this file](https://github.com/serpentine-h2020/serpentine/archive/refs/heads/main.zip) and extract to a folder of your choice (or clone the repository https://github.com/serpentine-h2020/serpentine if you know how to use `git`).
3. Open your terminal/command line/Anaconda prompt, navigate to the downloaded (extracted) folder `notebooks/solarmach` that contains the file `requirements.txt`, and run the following:

``` bash
$ conda create --name serpentine python=3.9
$ conda activate serpentine
$ pip install -r requirements.txt
```


# Run 
1. Open your terminal/command line/Anaconda prompt.
2. In the terminal, navigate to the downloaded (extracted) folder `notebooks/solarmach` that contains the file `solarmach.ipynb`
3. Make sure the corresponding conda environment is activated by running `conda activate serpentine` in the terminal.
4. Run `jupyter notebook solarmach.ipynb`
5. Your standard web-browser should now open the Jupyter Notebook.
