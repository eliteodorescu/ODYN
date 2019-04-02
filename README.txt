*******************************************************************************************************
					Welcome to ODYN! 
- Open-source software analysis tool to investigate space plasma turbulence and nonlinear DYNamics -
	****************************************************************************
	*   The latest version of ODYN can be found here:                          * 
	*   http://www.spacescience.ro/projects/odyn                               *   
	*                                                                          *
	*   You may use, redistribute, modify and distribute modified              *
	*   versions of ODYN according to the COPYRIGHT included in the package.   *
	****************************************************************************

=======================================================================================================
INSTRUCTIONS TO INSTALL
-------------------------------------------------------------------------------------------------------

ODYN is developed under Python open-source programming language 
Under all platforms (Linux/Unix, Max OS X, Windows) several packages need to be 
pre-installed in order for ODYN to run without errors:
- Python
- Numpy 
- Scipy 
- Matplotlib 
- Jupyter 
- h5py
- Networkx
- ffnet
- Spacepy
- CDF NASA library

The easiest way to install most of the packages is through Anaconda package 
(available for Linux/UNIX, Windows and Mac OS X): it includes 
the first 7 packages listed above. 
The 'ffnet' package requires fortran and c compilers and once 
these are correctly installed the installation of the last 3 packages is straightforward
=======================================================================================================

=======================================================================================================
STRUCTURE OF ODYN
-------------------------------------------------------------------------------------------------------

There are 5 sub-directories in the folder named ODYN
	-'AnalysisMethods' - all methods are collected in the Analysis.py file
	-'Config' - ODYN configurator file can be edited from this directory: CONFIG_FILE.txt
	-'Data' - default folder where test data are stored
	-'Results' - default folder where results of the analyses are saved   
	-'Notebooks' - contains the editable Jupyter notebook - analysis parameters and 
other usefull key parameters can be easily modified in this file
=======================================================================================================

=======================================================================================================
RUN INSTRUCTIONS
-------------------------------------------------------------------------------------------------------

In a terminal type: 'jupyter notebook'
A web-type page will open in the computer's default web browser.
Navigate to the folder ODYN -> Notebooks and click on the file with the .ipynb extention
The notebook is structured in cells that serve specific purposes, indicated in a header of each cell.
For example, the cell with the header:

######################################################################################
'''DATA PRE-PROCESSING: data-gaps management - interpolate, flagg erroneous data'''
######################################################################################

executes commands and function calls that deal with data-gaps, flagged data, etc.
and return time series of time-stamps and data with constant resolution (the missing data will 
be labeled with the chosen user flag, e.g. the default is set to 9999.)
To execute the commands of a cell, go to the notebook's toolbar to 'Cells->Run Cells'.
If the execution is done, in the upper left corned outside of the cell it will show: 'In [1]',
while if the execution is still in progress it will show: 'In [*]'.
It is possible to run all cells by choosing 'Cells->Run All'. Also, if something goes wrong,
restart the kernel and try again by going to 'Kernel->Restart & Run All'.
=======================================================================================================

=======================================================================================================
DEFAULTS IN ODYN
-------------------------------------------------------------------------------------------------------

Throughout the notebook there are several key parameters that can be modified - such parameters 
are preceded by a comment that starts with '''!!! REMOVE '#' AND CHOOSE ... !!!'''. 
There are also several boolean variables through which the user can choose to execute or not
a command or a set of commands by simply switching 'True' with 'False'. Such variables are preceded by:

'''!!!REMOVE '#' AND SET True OR False!!!'''

The values of all the above mentioned configurable parameters can also be modified through 
the ODYN configurator: "CONFIG_FILE.txt". All parameters have default values that can be easily 
changed in the preferred text editor.
 

BY DEFAULT, the follwing is happening in ODYN (if choosing to 'Run All Cells' as a first step):
    - Python3 is the default version (change to Python2 in the first cell, 
    by changing pyv='Python2', if this is what you have installed)
    - VEX test data are analyzed (there are 4 possible spacraft data 
    that can be read, Venus Express, Ulysses and Cluster, and USER_DATA for which a template 
    file: USER_DATA_TEMPLATE.txt and an example: USER_DATA_SAMPLE.txt have been provided, 
    such that the user can analyse any data, experimental or virtual, if written in the provided format)
    - computation of PSD is set to 'True' - PSD is computed for all variables
    - computation of PDF, SF, Flatness, ROMA and MI is set to 'False' - if switched to 'True',
    PDF, SF, Flatness, ROMA and MI are computed only for one chosen variable (the default 
    is the first data column, 0)
    - plotting of results is generally set to 'True' - the graphs are shown, on screen, under each cell
    - saving of results is generally set to 'False' - except for saving the masked time-series 
    representation (in the default 'Results' folder)
    
    ''' WARNING '''
    - in order for ODYN to work under Windows, the Spacepy package is not imported by default. To
    enable this feature, edit the 'ANALYSIS.py' file: 
	1) go to the 'AnalysisMethods' folder and open 'ANALYSIS.py' 
	2) uncomment the second line by removing the '#' sign: 
	the line:
		#from spacepy import pycdf
	will become:
		from spacepy import pycdf 	
=======================================================================================================

========
Authors
--------

Package Author
--------------
* Eliza Teodorescu - Researcher at the Institute of Space Science - ISS, Magurele, ROMANIA

Code Contribution
-----------------
* Marius Echim - Researcher at the Institute of Space Science - ISS, Magurele, ROMANIA
* Gabriel Voitcu - Researcher at the Institute of Space Science - ISS, Magurele, ROMANIA
* Costel Munteanu - Researcher at the Institute of Space Science - ISS, Magurele, ROMANIA
* Catalin Negrea - Researcher at the Institute of Space Science - ISS, Magurele, ROMANIA



.. _Python: http://www.python.org
.. _Numpy: http://www.numpy.org
.. _Scipy: http://www.scipy.org
.. _Matplotlib: http://matplotlib.org
.. _Jupyter: http://jupyter.org
.. _Spacepy: https://pythonhosted.org/SpacePy
=======================================================================================================
