# rtip_veros

This repository contains files and data that underlie the publication "Risk of tipping the overturning circulation due to increasing rates of ice melt" by Johannes Lohmann and Peter D. Ditlevsen.

The control simulations, as well as sensitivity with respect to salinity restoring time and wind forcing, have been analyzed using the files "climatology.py" and "control_40l.py".

The Veros setup file used for the freshwater hysteresis simulations is "global_flexible_res_freshw_hyst.py".

Analysis of the data concerning the ramped forcing simulations using different rates of change, as well as initial conditions has been done with the script "rtip_40l.py". 

An effort to cluster or linearly separate different ensemble realizations and their tipping outcome has been performed with the algorithms in the scripts "svm_bruteforce.py" and "svm_bruteforce_3d.py" (see also supplemental material).

Also provided are .txt files containing data used to plot the figures of the main text.
