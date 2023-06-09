# squid_analysis
For working example, run 

`cd squid_analysis`

` python .\code\squid_analysis.py -i -t .\input_data\rowsel_on_small\ -u .\input_data\rowsel_off_small\`

`-f` reads the `sq1.bias` files faster, but is more limited in what it can read (must be space separated and properly formatted.
`-s` flips the sign of the SAFB data

See output plots in `output_data/`

Ultimately, should cleanly make the plots seen in this posting: `http://bicep.rc.fas.harvard.edu/bicep_array/analysis_logbook/20230525_squid_tuning_devcur_analysis/squid_analysis/`
