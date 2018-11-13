#!/bin/bash

R -e "devtools::load_all('pkg');devtools::document('pkg')"
R CMD Rd2pdf --force --no-preview -o manual.pdf ./pkg


