# Fit runs

One directory per fit.  Inside each:

    fitpar.npy      the 37-slot parameter vector
    summary.json    chi2 of EVERY data set, fitted or not, plus slopes and ratios
    run.log         the console output of the fit and the plotting
    plots/          the figures

## Reading a plot

The file name says it:

    07_clas12_xs_out_U.png
       |     |      |   |
       |     |      |   +--- observable: U, T, LT, TT, LTp, or an asymmetry
       |     |      +------- in  = the set was fitted
       |     |               out = it was NOT, so the curve is a blind prediction
       |     +-------------- which data set
       +-------------------- order, so the directory sorts the way the paper reads

`00_chi2_summary.png` is the overview: one bar per set, dark if fitted,
brown if left out.

## The runs

    00_baseline       what the fit contained before tonight
    01_all            every data set we have
    02_no_clas12xs    everything except the preliminary CLAS12 cross sections
    03..13_drop_X     as 02, with set X removed -- so the X plots in that run
                      are its blind prediction

`COMPARISON.txt` puts all runs in one table.
