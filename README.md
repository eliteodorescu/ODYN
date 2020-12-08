# ODYN
Open-source software analysis tool to investigate space plasma turbulence and nonlinear DYNamics

ODYN – is a versatile modularized software library that wraps a comprehensive set of advanced data analysis methods meant to facilitate the study of turbulence, nonlinear dynamics and intermittency in space plasmas. Python programming language is used for the algorithmic implementation of models and methods devised to understand fundamental phenomena of space plasma physics e.g. elements of spectral analysis, probability distribution functions and their moments, multifractal analysis or information theory.
The software includes: 
-	Satellite data reading routines for the magnetic field of three space missions: Cluster, Venus Express and Ulysses. There is a template provided that enables any user data to be analyzed with ODYN
-	Satellite data pre-processing routines that include interpolation or management of data-gaps (in the sense that data-gaps could be ignored), or the possibility to flag the erroneous data
-	Visualization of experimental or simulated data
-	Software algorithms dedicated to classical and advanced analysis methods tailored, in principle, for the study of turbulence and non-linear dynamics of space plasma data:
o	PSD – computation of power spectral density with standard scipy tools
o	PDF – computation of probability density distributions
o	SF – computation of the moments of PDFs, namely, the structure functions and the estimation of Flatness/Kurtosis
o	PF – computation of partition functions and the estimation of the multifractal spectrum f(a)
o	ROMA – computation of the Rank-ordered Multifractal spectrum with the novel ROMA method
o	MI – computation of mutual information between two variables and the possibility to generate and extract the baseline thorough the randomization of the original data
-	Automatic data analysis feature. For example, the PDFs (or any of the implemented methods) can be computed automatically for an initial dataset.
-	User customization of the analysis parameters. A simple USER Configurator is provided that allows the user to adjust analysis parameters as needed, e.g., the segment length for PSD estimation, the binning of PDF histograms, the range of scales for PDF computation, the ranks of SF, and many more.

