# Patch for R devel

This is a submission to fix a problem that occurs in my package due to how class
now works in the development version of R. I received an e-mail by Kurt Hornik
to fix the problem in my package. I now use methods::is instead of class.

## Test environments
* local Arch GNU/Linux install, R 3.6.1
* ubuntu 16.04 (on travis-ci), R 3.6.1, devel
* win-builder (release, devel)

## R CMD check results
There were no ERRORs, WARNINGs or NOTES.
  
## Downstream dependencies
There are no downstream dependencies
