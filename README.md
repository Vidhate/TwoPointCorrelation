# TwoPointCorrelation
This repository contains multiple codes which calculate the 2 point correlation (DD/RR)-1<br>
<ul>
  <li><b>"BruteForce"</b> directory contains a C++ code developed by me for datasets with small number of points.<br></li>
  
  <li><b>"Comparison Plots"</b> holds 2 types of files. Files without the prefix "set2" are comparisons of Brute Force and KDTree. Files with the prefix "set2" are comparisons of all the 3 methods currently available in the repo.<br></li>
  
  <li><b>"Data"</b> directory contains the data input files used to calculate the 2 point correlations</li>
  
  <li><b>"KDTree"</b> directory uses in-built functions in Python imported from Sklearn to do the computations using a KDTree. <br></li>
  
  <li><b>"NFW Profile"</b> directory holds code to compute the NFW Dark Matter Halo Profile as described in their <a href="https://arxiv.org/abs/astro-ph/9508025">paper</a>, in Fourier space.</li><br>
  
  <li><b>"Normalize PS"</b> directory code computes the proper normalization to the Linear Theory Power Spectrum at redshift=0, from an input of a Transfer Function. Theory is described in Scott Dodelson's Modern Cosmology - Chapter on <i>Inhomogeinities</i>. It outputs the normalised Power Spectrum in the same directory, for later computations of Mass Function, Bias and Halo Profile. </li><br>
  
<li><b>"PKDTree"</b> directory uses a modification of the in-built scipy function. Modification done by Patrick Varilly. A few minute commits have been made by me to make the algorithm Python3 compatible (download latest from https://github.com/Vidhate/periodic_kdtree).</li><br>

<li><b>"Ratios"</b> directory hosts plots made by finding the ratios of correlation values from different input files. It also contains the ipynb code for plotting the ratios.</li><br>

<li><b>"Shadab"</b> directory contains code developed by Shadab Alam to calculate 2 point correlations.<br>

<li><b>"Tinker MF & Bias"</b> houses code to compute the Mass Function and large scale Bias as described by Tinker et.al. in their <a href="https://arxiv.org/abs/0803.2706">Mass Function paper[2008]</a> and their <a href="https://arxiv.org/abs/1001.3162">Bias paper[2010]</a>.</li>
  
<li><b>"pickleMocks.py"</b> helps to extract positions from the data files and pickle dump them as numpy arrays.<br>
</ul>

## Prerequisites
<ul>
  <li>C++ compiler
  <li>Python (package <i>sklearn</i> should be installed)
  <li>Jupyter Notebook
  <li>Python3 compatible periodic_kdtree should be available (download: https://github.com/Vidhate/periodic_kdtree and setup)
</ul>

## Usage
To run for your own data sets please note the following things:
<ul>
  <li>Pickle your input data by running <i>pickleMocks.py</i> with the correct paths to input dataset</li>
  <li>All the codes have a variable <i>path</i> or <i>paths</i> that specifies where to derive the pickle dumps of input data. The <i>path</i> variable should point to the directory in which your input file is. The <i>names</i> or <i>fileNames</i> vector/list variable should have the name of the input files you will be using.</li>
  <li>To run the C++ code kindly do <br>
    <i> >> g++ BruteForce.cpp<br>
      >> ./a.out</i>
  </li>
  <li>Rest of the codes can be run in the Jupyter Notebook.</li>
  <li>Please make sure that the <i>periodic_kdtree</i> repo is available on your system and the path is set properly before running the Jupyter Notebook in the PKDTree directory.</li>
</ul>
  
## Output
The output of correlation values found from each technique will be stored in the same directory which holds the code in the <i>.txt</i> format, with the same name as the input file name prefixed with an appropriate abbreviation of the technique.<br>
The output of comparison plots is stored in the folder "Comparison Plots". The prefix "setX" can be changed in the code.
