# Sahul-wide human population growth
<img align="right" src="www/popgr.png" alt="population growth icon" width="200" style="margin-top: 20px">

Sahul-wide human population growth estimates from temporal distribution of archaeological dates between 40 ka and 5 ka

This R code recalculates Sahul-wide human population growth using a resampling approach (values originally presented in Williams, A.N., 2013. <a href="https://royalsocietypublishing.org/doi/full/10.1098/rspb.2013.0486">A new population curve for prehistoric Australia</a>. <em>Proc. R. Soc. Lond. B</em> 280, 20130486)

## General approach

1. calibrate raw radiocarbon dates to calendar year using the <code>rcarbon</code> package (using the <a href="https://c14.arch.ox.ac.uk/oxcalhelp/hlp_curves.html"><em>shcal20</em></a> calibration curve for non-marine dates, and the <a href="https://c14.arch.ox.ac.uk/oxcalhelp/hlp_curves.html"><em>marine20</em></a> calibration curve for marine dates (e.g., marine shells)
2. bin calibrated dates into 200-year interval bins (iteratively, sampling from calibrated range per date)
3. apply smoothing spline (df = 25) to number of dates per bin (iteratively)
4. apply taphonomic correction (i.e., number of dates necessarily declines with age because of taphonomic loss - mean correction data from: Williams, A.N., 2012. <a href="https://www.sciencedirect.com/science/article/abs/pii/S0305440311002482">The use of summed radiocarbon probability distributions in archaeology: a review of methods</a>. <em>J. Archaeol. Sci.</em> 39, 578-589)
5. calculate Williams' (2013) mean annual population growth (GRann) and instantaneous rate of exponential growth (<em>r</em>) for each resampled, corrected series: 

GRann = 0.5(<em>n</em><sub><em>i</em>-1</sub>-<em>n</em><sub><em>i</em></sub>)/<em>n</em><sub><em>i</em></sub>

<em>r</em> = log(<em>n</em><sub><em>i</em></sub>/<em>n</em><sub><em>i</em>-1</sub>)

where <em>n</em><sub><em>i</em></sub> = number of dates in the older bin, <em>n</em><sub><em>i</em>-1</sub> = number of dates in the younger bin


6. correlate GRann & <em>r</em> for each iteration
7. calculate 95% confidence interval of correlation
8. plot reconstructed population growth curves per temporal bin (± 95% confidence limits)
9. plot taphonomically corrected number of dates per temporal bin ((± 95% confidence limits)

## Required R packages

- <code>ggplot2</code>
- <code>stratigraphr</code>
- <code>DescTools</code>
- <code>rcarbon</code>
- <code>stringr</code>
