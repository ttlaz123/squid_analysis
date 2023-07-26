# Installation:
Go to folder you want to save the code, then type
`git clone -b tuning_clean --single-branch https://github.com/ttlaz123/squid_analysis/`
# squid_analysis
NOTE: Python 3.6+ required
For working example, run 

`cd squid_analysis`

windows systems: ` python .\code\squid_analysis.py -i -t .\input_data\rowsel_on_small\ -u .\input_data\rowsel_off_small\`

unix systems: ` python ./code/squid_analysis.py -i -t ./input_data/rowsel_on_small/ -u ./input_data/rowsel_off_small/`

* `-f` reads the `sq1.bias` files faster, but is more limited in what it can read (must be space separated and properly formatted.
* `-s` flips the sign of the SAFB data
* `-p` provides chosen bias points
* `-l` chooses the maximum number of columns to check
* `-h` for more detailed help info

### To choose optimal bias points, run the following:
`python .\code\squid_analysis.py -f -t .\input_data\rowsel_on_small\ -u .\input_data\rowsel_off_small\ -c .\tune_cfg\slac_cd19.cfg -i -p .\output_data\rowsel_on_small\units_dac\rowsel_on_small_col_biases.csv` 
(Yes the `csv` file should be in input data; I'll fix it later, but you can put any file anywhere as long as you point to it correctly while calling the script).

There are a few ways of choosing the optimal biases automatically, `naive`, `bias_current`, and `device_current`. Will add options for these soon.

## RS Servo analysis
` python ./code/squid_analysis.py -r -t ./input_data/rowsel_on_small/ `

Same as standard crosstalk analysis, but `-r` instead of `-i` indicaes that you wish to run the rs_servo analysis instead

### See output plots in `output_data/`

Ultimately, should cleanly make the plots seen in this posting: `http://bicep.rc.fas.harvard.edu/bicep_array/analysis_logbook/20230525_squid_tuning_devcur_analysis/squid_analysis/`



