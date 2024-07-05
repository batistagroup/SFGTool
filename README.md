# SFGTool
Python script for calculating vSFG, double resonant SFG, and PM-IRRAS used in a Tutorial paper (in prep.). <br>
For now it can only read in Gaussian output.
## Prerequisites
- Python>=3
- numpy
- matplotlib
## Basic usage
Run the script with <br>
<code>python sfg.py sfg.inp</code>
## SFG Examples
Example SFG inputs (Gaussian outputs and input parameter files) are included in the <code>exmaple</code> folder. <br>
To run the SFG exmaple, go to the <code>exmaple</code> folder, and run <br>
<code>bash generate_results.sh</code>
## PM-IRRAS Example
PM-IRRAS exmaple is provided in the <code>pmirras_example</code> folder <br>
To run the SFG exmaple, go to the <code>pmirras_exmaple</code> folder, and run <br>
<code>python ../sfg.py pmirras.inp</code>
