# TwoPointCorrelation
This repository contains numerical and analytical calculations to compute 2 point correlations as Xi(r)=n(DD)/n(RR)-1 in the context of Dark Matter Halos and Large Scale Structure in the Universe.<br>
<ul>
  <li><b>"BruteForce"</b> directory contains a C++ code developed by me for computing 2PC on datasets with small number of points.<br></li>
  
  <li><b>"Comparison Plots"</b> The Comparison code housed in this directory functions to plot multiple outputs of 2PC computations given their paths. Helpful in comparing Analytic and NUmerical computations. Recently modified to check ratios wrt a baseline file for better comparisons.<br></li>
  
  <li><b>"Data"</b> directory contains the data input files used to calculate the 2 point correlations.</li>
  
  <li><b>"NFW Profile"</b> directory holds code to compute the NFW Dark Matter Halo Profile as described in their <a href="https://arxiv.org/abs/astro-ph/9508025">paper</a>, in Fourier space.</li><br>
  
  <li><b>"Normalize PS"</b> directory code computes the proper normalization to the Linear Theory Power Spectrum at redshift=0, from an input of a Transfer Function. Theory is described in Scott Dodelson's Modern Cosmology - Chapter on <i>Inhomogeinities</i>. It outputs the normalised Power Spectrum in the same directory, for later computations of Mass Function, Bias and Halo Profile. Note that the two modes in which the code can output - 1. Finessed (interpolated at more values) 2. Same as k (no interpolations)</li><br>
  
<li><b>"PKDTree"</b> Computes 2PC considering Periodicity in the given dataset. This is brought about by wrapping the original dataset in its own images to make a bigger volume, while computing 2PC on the original volume(dataset) at the center.</li><br>

<li><b>"Shadab"</b> directory contains code developed by Shadab Alam to calculate 2 point correlations.<br>

<li><b>"Tinker MF & Bias"</b> houses code to compute the Mass Function and large scale Bias as described by Tinker et.al. in their <a href="https://arxiv.org/abs/0803.2706">Mass Function paper[2008]</a> and their <a href="https://arxiv.org/abs/1001.3162">Bias paper[2010]</a>.</li>
  
<li><b>"pickleMocks.py"</b> helps to extract positions from the data files and pickle dump them as numpy arrays.<br></li>
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
</ul>
  
## Output
The output of all the computations will be stored in the same directory which holds the code in the <i>.txt</i> format, with the same name as the input file name prefixed with an appropriate abbreviation of the technique.<br>
The output of comparison plots is stored in the folder "Comparison Plots". The prefix "setX" can be changed in the code.
