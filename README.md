## Description

This repository provides examples for orbit determination using the open-source
orbital dynamics library Orekit (https://www.orekit.org/) and in particular the
Python wrapper developed by @petrushy: 
https://gitlab.orekit.org/orekit-labs/python-wrapper. 

Copied from https://github.com/GorgiAstro/orbit-determination-examples/ and 
adapted to get rid of the Jupyter Lab notebook dependency and to work with 
current (2024) version of Orekit (v12).

Tested with Debian/GNU Linux sid:

1. Install Anaconda3

See https://docs.anaconda.com/free/anaconda/install/linux/

2. Install Orekit (Python wrapper) for Conda

See https://www.orekit.org/download.html

```
conda install -c conda-forge orekit
``` 
will also install (Java) dependences.

3. Check Conda is working

Launch conda installed in ``/yourpath``:
```
eval "$(/yourpath/anaconda3/bin/conda shell.bash hook)" 
python3
import orekit
```
must not complain of a missing library

4. When first executing, uncomment (e.g. line 28 of ``00-generate-position-data.py``)
```
download_orekit_data_curdir()
```
to download https://gitlab.orekit.org/orekit/orekit-data/-/archive/master/orekit-data-master.zip, and
unzip to acceleretate the next examples. This line can be commented back once the dataset has been
downloaded.

Examples taken from https://github.com/GorgiAstro/laser-orbit-determination/blob/master/02-orbit-determination-example.ipynb
and https://github.com/GorgiAstro/tle-fitting/blob/master/01-tle-fitting-example.ipynb, converted to Python and updated
to Orekit v.12 API.
