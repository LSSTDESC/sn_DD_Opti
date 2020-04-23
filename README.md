# sn_DD_Opti

A framework to optimise DD mini-surveys

```
This software was developed within the LSST DESC using LSST DESC resources, and so meets the criteria
 
given in, and is bound by, the LSST DESC Publication Policy for being a "DESC product". 
We welcome requests to access code for non-DESC use; if you wish to use the code outside DESC 
please contact the developers.

```
## Release Status

This code is under development and has not yet been released.



## Feedback, License etc

If you have comments, suggestions or questions, please [write us an issue](https://github.com/LSSTDES
C/sn_DD_opti/issues).

This is open source software, available for re-use under the modified BSD license.

```
Copyright (c) 2020, the sn_DD_opti contributors on GitHub, https://github.com/LSSTDESC/sn_DD_opti/graphs/co
ntributors.
All rights reserved.
```

## ** How to run this package **

### Modify your PYTHONPATH
- export PYTHONPATH=sn_DD_opti:$PYTHONPATH

### Two GUIs can be displayed
#### Nvisits vs zlim
- python scripts/showGui.py --show Visits --cadence cad (cad is the cadence: default: 3)

#### Budget vs zlim
- python scripts/showGui.py --show Budget
