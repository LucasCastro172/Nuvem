#! /bin/sh
#File: model.sh

# Set messages on
set -x

#Experiment Number
#num=1

datafile=camadas_planas.dat


trimodel xmin=0 zmin=0 xmax=6 zmax=2 \
	1 xedge=0,6 \
          zedge=0,0  \
          sedge=0,0  \
	2 xedge=0,6 \
	  zedge=0.8,0.8 \
	  sedge=0,0 \
	3 xedge=0,6 \
	  zedge=1.6,1.6 \
	  sedge=0,0 \
     4 xedge=0,6 \
        zedge=2,2 \
        sedge=0,0 \
	  sfill=1,0.3,0,0,0.44,0,0 \
	  sfill=1,1,0,0,0.34,0,0 \
	  sfill=1,1.8,0,0,0.27,0,0 \
	  kedge=1,2,3,4 \
        >$datafile

# Cria um Arquivo *.PS do Modelo
 spsplot <$datafile >camadas_planas.eps \
         title=" Modelo Plano" \
 	titlecolor=black axescolor=black gridcolor=white \
         labelz="Profundidade (km)" labelx="Distancia" legend=1\
 	xbeg=0.0 xend=6.0  \
 	zbeg=0  zend=2.0 \
 	dxnum=1.0 dznum=0.5 \
 	gridx= solid gridz= dash\
 	gedge=0.1 gtri=2.0 gmin=0.2 gmax=0.8 \
        wbox=10.0 hbox=5.0 &

sxplot < $datafile -geom 1350x450+0+0 edgecolor=black tricolor=none title="Modelo em profundidade" \
		label1="Profundidade (km)" label2="Distancia (km)" &

#evince modelo.eps &

exit

