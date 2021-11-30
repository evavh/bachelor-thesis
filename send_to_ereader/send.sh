#!/bin/bash

date=$(date +"%Y-%m-%d")

pdflatex --output-directory=/tmp --jobname="simulation-main.py_${date}" send_to_ereader/simulation-main.py.tex && rmapi put /tmp/simulation-main.py_${date}.pdf Uni/Code

pdflatex --output-directory=/tmp --jobname="analysis-main.py_${date}" send_to_ereader/analysis-main.py.tex && rmapi put /tmp/analysis-main.py_${date}.pdf Uni/Code
