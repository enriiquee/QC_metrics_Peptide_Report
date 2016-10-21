#!/bin/bash

# To use this script, pandoc must be installed
# sudo apt install pandoc

R -e "rmarkdown::render('QC_metrics_for_Peptide_Reports.Rmd')"

