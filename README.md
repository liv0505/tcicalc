# tcicalc
Calculate useful indexes for tropical cyclones:
* Potential intensity: equation(3) in [Bister and Emanuel (2002)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2001JD000776) under pseudoadia-batic assumptions
* Entropy deficit: equation(2) in [Tang and Emanuel(2012b).](https://journals.ametsoc.org/doi/abs/10.1175/BAMS-D-11-00165.1)
* Ventilation index: equation(1) in [Tang and Emanuel(2012b).](https://journals.ametsoc.org/doi/abs/10.1175/BAMS-D-11-00165.1)
* Genesis potential index 2004: [Emanuel and Nolan 2004](http://texmex.mit.edu/pub/emanuel/PAPERS/em_nolan_extended_2004.pdf)
* Genesis potential index 2010: equation(1) in [Emanuel 2010](https://agupubs.onlinelibrary.wiley.com/doi/10.3894/JAMES.2010.2.1)

## File descriptions:
**1. pcmin.f90** is a collection of five subroutines written in FORTRAN 90 by Prof. Emanuel. It computes the maximum wind speed and minimum central pressure achievable in tropical cyclones. See Prof. Emanuel's [original fortran version](ftp://texmex.mit.edu/pub/emanuel/TCMAX/).
* The first four subroutines (`pcmin3, pcmin3_kflag, pcmin2 and pcmin`) could all be used to calculate the potential intensity of a tropical cyclone. The detailed description of PI could be found [here](https://emanuel.mit.edu/limits-hurricane-intensity). 
* The fifth subroutine calculates `CAPE`, which will be used in other four subroutines in this f90 file.
* To call pcmin.f90 subroutines in python, we need to build a python module by `numpy.f2py`. See [F2PY Users Guide](https://numpy.org/doc/stable/f2py/) for details.

**2. calculator.py** is the main python script that calls necessary modules from other `*.py` files to calculate TC indexes using data from **ERA5_Data** folder and save the results to **TC_example** folder. Detailed description could be found in this code.

**3. get_figure.py** is the python script that plot the results and save the figures to **TC_example** folder.

## Quick Start:
Before start, you should modify the `datadir` and `outdir` (directories where you store input data and output results.) in **calculator.py**

**Step 1. Compile FORTRAN 90 file:**<br>
Enter the following code into your terminal:<br>
>`f2py -c pcmin.f90 -m tcpi`<br>

And you will get a `tcpi*.so file`, (e.g. `tcpi.cpython-37m-darwin.so`), which could be called in python.<br>

**Step 2. Run calculator.py:**<br>
First, modify the year and directory names in **calculator.py**, i.e. `year`, `datadir` and `outdir`  in main().<br>
Enter the following code into your terminal:<br>
>`python /...(rootpath).../calculator.py`<br>

You will find the results in the `outdir` you just set (my sample results are in **TC_example** folder).<br>

**Step 3. Run get_figure.py:**<br>
First, modify the year and directory name in **get_figure.py**, i.e. `year` and `outdir`  in main().<br>
Enter the following code into your terminal:<br>
>`python /...(rootpath).../get_figure.py`<br>

The figures could be found in the `outdir` you just set (my sample results are in **TC_example** folder).<br>
