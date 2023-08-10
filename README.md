# Supplementary materials: Fiddler crabs detect second-order motion in both intensity and polarization (submitted/under review)
### Supplementary information for: *Fiddler crabs detect second-order motion in both intensity and polarization.* Samuel P. Smithers, Maisie F. Brett, Martin J. How, Nicholas E. Scott-Samuel, and Nicholas W. Roberts.

## Repository contents
- Raw date from intensity and polarization experiments (.csv files)
- ```Stats and graphs for fiddler crabs see SO motion.r```: R code used to generate the figures for, and perform the statistical analysis on, the data from the intensity and polarization experiments as reported in the manuscript.
- Movie 1: Original RGB video of mudflat showing fiddler crabs and flying tern.
-	Movie 2: Original degree of polarization video of mudflat showing fiddler crabs and flying tern.
-	Movie 3: Original receptor contrast video of mudflat showing fiddler crabs and flying tern. 
-	Movie 4: *Gelasimus vomeris* resolution RGB video of mudflat showing fiddler crabs and flying tern.
-	Movie 5: *Gelasimus vomeris* resolution degree of polarization video of mudflat showing fiddler crabs and flying tern.
-	Movie 6: *Gelasimus vomeris* resolution receptor contrast video of mudflat showing fiddler crabs and flying tern.
-	Movie 7: Expanding negative control stimulus example. 
-	Movie 8: Expanding first-order stimulus example. 
-	Movie 9: Expanding mean grey stimulus example. 
-	Movie 10: Expanding second-order flicker stimulus example. 
-	Movie 11: Flicker control example.
-	Movies 12-16: Example recordings and motion traces for a response
-	Movies 17-21: Example recordings and motion traces for no response
-	Movies 22-26: Example recordings and motion traces for a rejection

## Viewing the supplementray movies
To view the supplementary movies click on the movie and on the next page click "View raw" to download the movie to your device. The movies can be viewed using any device/program capable of playing .mp4 files.  

## Statistical analysis & plotting (R)
### Requirements & dependencies
- R version >= 4
- The following packages must be installed:
  - lme4
  - dplyr
  - DHARMa
  - ggplot2
  - Hmisc
  - Matrix
  - bannerCommenter (optional)

#### Running instructions
1. Download ```Stats and graphs for fiddler crabs see SO motion.r``` from this repository. You will also need to download the "*Data from intensity experiment.csv*" and "*Data from polarization experiment.csv*" data files. Alternatively, you can clone this repository ([instructions for Rstudio here](https://datacarpentry.org/rr-version-control/03-git-in-rstudio/index.html)). 
2. Install any uninstalled dependencies listed above.
3. Open ```Stats and graphs for fiddler crabs see SO motion.r``` and follow instructions within the script.

