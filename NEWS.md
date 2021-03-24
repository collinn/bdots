# bdots 03-02-2021
- Ready for launch!

# bdots 03-24-2021
- Temporary inclusion of `doubleGauss2` fitting function. A bit slower than `doubleGauss`, though with generally better fits
- added `paramDT` function to create `data.table` of subject/groups and fitted coefs which can be input to `bdotsRefit`
- Modified curve fitting to be less critical of marginal fits, decreasing number of `NULL` fits
- updated vignettes
- added temporary `writeCSV` function to write out `bdotsBoot` curve information for plotting

