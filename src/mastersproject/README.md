# Overview

Information about the In-Situ Stimulation and Circulation experiment at the Grimsel Test Site can be found here: http://www.grimsel.com/gts-phase-vi/isc-in-situ-stimulation-circulation-experiment

## Data extraction
The class which extracts raw data from the In-Situ Stimulation and Circulation (ISC) 
experiment at the Grimsel Test Site (GTS) is found in `GTS/ISC_data/isc.py`.
The raw data is loacted in `GTS/01BasicInputData`.

Various extractions of this data can be found in `GTS/Method prototypes` (which exist for the purpose of quick testing).

# Model setup

The most advanced model currently being tested is Biot's equations with contact mechanics for fractured porous media.

The setup, which bridges the ISC data with the porepy project (our group's modelling tool, see: https://github.com/pmgbergen/porepy)
is located in `GTS/isc_modelling/contact_mechanics_biot.py`.

Results (from various stages of the testing process) can be found in `GTS/isc_modelling/results`. The files can be opened in e.g. ParaView.
