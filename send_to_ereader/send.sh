#!/bin/bash

date=$(date +"%Y-%m-%d")

pdflatex --output-directory=/tmp --jobname="main.py_${date}" send_to_ereader/main.py.tex && rmapi put /tmp/main.py_${date}.pdf Uni/Code
