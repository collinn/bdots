# bdots 03-02-2021
- Ready for launch!

# bdots 03-24-2021
- Temporary inclusion of `doubleGauss2` fitting function. A bit slower than `doubleGauss`, though with generally better fits
- added `paramDT` function to create `data.table` of subject/groups and fitted coefs which can be input to `bdotsRefit`
- Modified curve fitting to be less critical of marginal fits, decreasing number of `NULL` fits
- updated vignettes
- added temporary `writeCSV` function to write out `bdotsBoot` curve information for plotting

# bdots 05-03-2021
- Fixes to choosing AR1 status during refit steps

# bdots 07-29-2021
- Major updates to fit, parameter, and bootstrap plots with ggplot2
- Handle case in which observation has zero variance in outcome
- Relaxed fit criteria in refitting step
- Fixed bugs in `bdotsRefit` walkthrough

# bdots 05-26-2022
- Major updates to plotting functions
- Added correlation function, `bdotsCorr`
- Fixed bugs in `bdotsRefit`, added option to save/restore refit 
- Added new curve functions

# bdots 07-25-2022
- Time values no longer need to be identical between subjects for correct plotting
- Corrected error writing out bootstrapped values when there are no significant times
