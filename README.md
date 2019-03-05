# TwoPointCorrelation
This repository contains multiple codes which calculate the 2 point correlation (DD/RR)-1<br>
<ul>
<li><b>"BruteForce"</b> directory contains a C++ code developed by me for datasets with small number of points.<br>
<li><b>"KDTree"</b> directory uses in-built functions in Python imported from Sklearn to do the computations using a KDTree.<br>
<li><b>"PKDTree"</b> directory uses a modification of the in-built scipy function. Modification done by Patrick Varilly. A few minute commits have been made by me to make the algorithm Python3 compatible (download latest from https://github.com/Vidhate/periodic_kdtree).
<li><b>"Shadab"</b> directory contains code developed by Shadab Alam to calculate 2 point correlations.<br>
<li><b>"Comparison Plots"</b> holds 2 types of files. Files without the prefix "set2" are comparisons of Brute Force and KDTree. Files with the prefix "set2" are comparisons of all the 3 methods currently available in the repo.<br>
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
