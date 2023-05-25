# squid_analysis

For analyzing BA squid tuning data

By Tom Liu, adapting script from David Goldfinger and Shawn Henderson

Includes sample input and output data from L9 (cols 0-15) and L1 (cols 16-31). Excludes `*_RCs_sq1servo_sa.bias` files due to size.

Typical run involves `python squid_analysis.py -c /path/to/config/file -t /path/to/sq1/bias/file -u /path/to/sq1.bias/file/with/rsoff -i`. Use `python squid_analysis.py -h` to see more options.

Requirements: `Python version > 3.6`, 
Standard packages: `scipy`, `glob`, `pandas`, etc,
Input data: Standard data output from MAS with folder structure of `ctime/ctime_*files*`
