# Main function #
## Spin2Win ##
GUI (split in Spin2Win.m and Spin2Win.fig)
Requires:
  * Main\_calc
  * Plotter
  * ReadFile
  * import\_n\_format
  * load\_default\_iso
  * load\_default\_ort
  * merge\_indata
  * split\_indata

# Core calculation functionality #
## Main\_calc ##
Preassure and stress calculator for composite and isotropic materials
Requires:
  * d\_calc.m - Interface pressure calculator
  * en\_cr\_calc.m - Stress distribution calculator
  * _M\_calc - Energy, Mass, Cost calculator_

## Pd\_calc.m ##
Interface pressure calculator

## Ten\_cr\_calc.m ##
Stress distribution calculator

## E\_M\_calc.m ##
Energy, Mass, Cost calculator

# Plotting functionality and GUI #
## Plotter.m ##
Plots data calcualted by Main\_calc

# Data handling functions #
## Figure\_data\_extract.m ##
Figure data extractor for subplotted figures
Made to extract Spin2Win figure data.
Extracts all data from currently active figure

Note: Currently standalone, to be implemented into a data extraction tool.

## ReadFile.m ##
Read old run files and finds line numbers of diffrent runs

## split\_indata.m ##
Splits raw indata to geometrical and settings data

## merge\_indata.m ##
Merges geometrical and settings data into a single matrix for easy transport between sub functions
_Suggestion: Readability improvement - use structs instead?_

## load\_default\_ort.m ##
Creates default orthotropic material-dataset
and returns them as a matrix.

## load\_default\_iso.m ##
Creates default isotropic material-dataset
and returns them as a matrix.

## import\_n\_format.m ##
Format input data from ReadFile.m raw data and creates idat-matrix

_Suggestion: Readability improvement - use structs instead?_