##Resubmission
This is a resubmission. In this version I have:

* removed the print command in the function updt_winner_weights (it is now a warning)
* added \value to the run_exp.Rd file and explained the function's results (note that in the first submission several helper functions were exported which are not needed by the end-user; these are now internal functions and thus the documentation was also removed for them)
* added more details about the package and the background of the PASS-T model, so that it is more clear why the package is useful
* added references for the PASS-T model and empirical studies related to it

## Test environments
* local Arch GNU/Linux install, R 3.6.1
* ubuntu 16.04 (on travis-ci), R 3.6.1
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs.

There was 1 NOTE on win-builder devel and release:

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Johannes Titz <johannes.titz@gmail.com>'

  As I understand it this is just a reminder for CRAN maintainers.
  
## Downstream dependencies
There are no downstream dependencies
