#!/bin/bash
#
# prototipo de modelagem por diferencas finitas
#
# script para modelagem:
#
#arquivo com modelo de velocidade:
#vmodel=marm_vel.hdr
vmodel=homog.hdr
#posicao da fonte em metros:
zfonte=250.0
xfonte=4500.0
#posicao do geofone em metros:
zgeof=1250.0
xgeof=8500.0
#frequencia pico do pulso fonte (Hz):
freq=15.0
#tempo total de modelagem (s):
ttime=3.5
#arquivo de saida com evolucao do campo: 
sumovie=afd_movie.su
#
export OMP_NUM_THREADS=4
#

for i in {1..5}
do
#
afd_2d.x<<EOF
${vmodel}
${zfonte} ${xfonte}
${zgeof} ${xgeof}
${freq}
${ttime}
${sumovie}
EOF
done
#
# para visualizacao usando Seismic Un*x
#
suxmovie < ${sumovie} \
	   loop=2 title='Frame '%g clip=0.01 \
	   label1='Profundidade (m)' label2='Distancia (m)' 
