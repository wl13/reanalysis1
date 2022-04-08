### Description ###

'estimate-dinuc-equilibria.py' is a python script written by Alexander T. Ho to estimate dinucleotide equilibrium frequencies given a mutational matrix. Information from the mutational matrix is extracted to define a series of simultaneous equations where mutational equilibrium is equal to the rate of mutational gain minus the rate of mutational loss. These simultaneous equations are resolved using NumPy. Equilibrium of the final dinucleotide (TT) is estimated as 1 minus the sum of all the other dinucleotide equilibria so that all dinucleotide equilibrium estimates sum to 1. Bootstraps are created by resampling the input data with replacement and repeating the described methodology. The input file required and output file created are described below.

### Dependencies/Imported packages ###

- import os
- import numpy as np
- import csv
- import random
- import scipy.stats

### Input file (edit path within script) ###

Complete mutational matrix TSV file with headers:

  - ref (the reference dinucleotide e.g. 'AA')
  - alt (the alternate dinucleotide e.g. 'AC')
  - count (the number of ref to alt events e.g. AA to AC events)
  - total (the total occurence of the ref dinucleotide in the sequence of interest e.g. whole genome or smaller region)

Note that all dinucleotide ref/alt combinations must be included even if the number of events is 0.

Example input TSV (first 5 lines):

		ref	alt	count	total
		aa	ac	18	13775613
		aa	ag	39	13775613
		aa	at	35	13775613
		aa	ca	19	13775613

### Output file (edit path within script) ###

TSV file with headers:

- dinucleotide (the dinucleotide for which mutational equilibrium has been estimated)
- estimate (the mutational equilibrium estimate)
- bootstrapped_mean (the mean estimate of X bootstraps, note the number of bootstraps may be adjusted in the script)
- bootstrapped_lower_95_ci (the lower 95% CI limit of X bootstraps)
- bootstrapped_upper_95_ci (the upper 95% CI limit of X bootstraps)

Example output TSV (first 5 lines)

		dinucleotide	estimate	bootstrapped_mean	bootstrapped_upper_95_ci	bootstrapped_lower_95_ci
		AA	0.22440274995551400	0.22623445068256500	0.2222291691149060	0.23023973225022400
		AA	0.045473718539784600	0.044816864175406	0.043642893702403300	0.045990834648408700
		AG	0.03334210399338560	0.03384595464423510	0.0330570288459714	0.034634880442498900
		AT	0.17358395823923800	0.1737686329450380	0.17100939880247200	0.1765278670876050
