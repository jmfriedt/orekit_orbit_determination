## Description

This repository provides examples for orbit determination using the open-source
orbital dynamics library Orekit (https://www.orekit.org/) and in particular the
Python wrapper developed by @petrushy: 
https://gitlab.orekit.org/orekit-labs/python-wrapper. 

Copied from https://github.com/GorgiAstro/orbit-determination-examples/ and 
adapted to get rid of the Jupyter Lab notebook dependency and to work with 
current (2024) version of Orekit (v12).

Test with Debian/GNU Linux sid:

1. Install Anaconda3

https://docs.anaconda.com/free/anaconda/install/linux/

2. Install Orekit for Conda

See https://www.orekit.org/download.html

```
conda install -c conda-forge orekit
``` 

3. Check Conda is working

Launch conda installed in ``/yourpath``:
```
eval "$(/yourpath/anaconda3/bin/conda shell.bash hook)" 

```

python3:
```
import orekit
```
must not complain of a missing library
