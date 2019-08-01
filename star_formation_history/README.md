Star formation history plot
===========================

This plots the star formation history of your simulation (from SFR.txt),
possibly against other simulations, with observational data overplotted.
For references please see the `loadObservationalData.py` script.

To plot simulations, add them to the dictionary in plotSFH.py, called
'simulations'. They are key-value pairs, with the key being the path
to the simulation directory (i.e. the directory where SFH.txt lives),
and the value being the 'fancy' name you would like in the legend.
