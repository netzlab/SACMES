# SACMES
This repository contains the master script for the SACMES program and all updates. 

SACMES stands for Simultaneous Analysis and Control of Electrochemical Systems.

Each script analyzes data from specific electrochemical techniques in real time and exports the data into space delimited txt files. Each script runs off of an multithreaded animation module known as ElectrochemicalAnimation and can analyze any number of electrodes simultaneously.

# SACMES_SWV.py #
SACMES_SWV analysis data from Square-Wave Voltammograms and offers real time control of data analysis.

11/12/2019
SACMES_SWV_FrequencyMap.py has been consolidated into SACMES_SWV.py in a single script.

11/11/2019
New settings toolbar includes options for changing parameters used to read the user's data file including different file extensions (.txt, .csv, and .DTA), the delimiter between data columns, the spacing between current columns, the column number for potentials, and the column number for currents. 

# SACMES_CV.py #
SACMES_CV analyzes data from Cyclic Voltammgrams

# SACMES_CA.py
SACMES_CA analyzes data from Chronoamperograms

