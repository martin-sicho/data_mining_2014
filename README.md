Data mining project - 2014
================

This repo was created for a school project. It is an implementation of a very simple classification and QSAR modelling tool.

WARNING: The code is mostly ugly :)

### How to work with it
Each module has `params.py` file where the parameters of the module can be specified. 
You will probably want to configure at least the `datageneration` module, which contains definitions of file paths and some basic parameters such as the protein target, clustering threshold...

Please note that you will have to download the decoys as an `.sdf` file and point to the proper loacation via `datamanipulation.params.DECOYS_SDF_FILE_PATH`. Actives should be downloaded automatically.

The modules are loaded from the `main.py` file in the root of the repo.
The `main.py` file includes a method that builds models, a prediction method for pIC<sub>50</sub> predicitons and the `main` method. The `main` method includes a simple usage example.
