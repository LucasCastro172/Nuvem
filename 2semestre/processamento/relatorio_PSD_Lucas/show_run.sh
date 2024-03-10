#!/bin/bash
output=relatorio_PSD
file=relatorio_PSD
pdflatex ${file}.tex #compilar o arquivo e gerar o dvi
bibtex ${file}    #gerar a bibliografia
pdflatex ${file}.tex #compilar o arquivo e gerar o dvi
pdflatex ${file}.tex #compilar o arquivo e gerar o dvi


#dvipdf ${file}.dvi #converter o arquivo dvi para pdf

mv ${file}.pdf ${output}.pdf

evince ${output}.pdf & #visualizar o arquivo pdf

