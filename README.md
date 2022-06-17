## Bachelor Research Project ##
### Automating the analysis of first binary formation in core collapse ###
This repository contains all the code and scripts for my bachelor research project. The Python code for simulating core collapse in a "globular" cluster of 1000 stars is based on an example script (see initial commit) my supervisor gave me.

To run it, you first need to install [AMUSE](https://amusecode.org/) and its prerequisites, and make. Then make sure you are in the AMUSE virtual environment (if you created it) and execute `make run` in the root directory of the repo.

The parameters for both the simulation and the analysis Python scripts are laid out when running `python3 main.py --help` in the respective folder.

The integration tests in the testing directory are mostly out of sync with the state of the code. If you want to work with the code, I highly suggest you make sure they pass.
