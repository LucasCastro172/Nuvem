#! /bin/sh
#-------------------------------------------------------------------------------
# Adiciona ruido a secao sismica
#-------------------------------------------------------------------------------
indata=camadas_planas_cdp.su    # dado de entrada
outdata=ruido_cdp.su           # dado com ruido
sn=500                           # razao sinal/ruido
noisetype=gauss                 # tipo de ruido
seed=from_clock                 # semente de partida

suaddnoise <$indata >$outdata  sn=$sn  noise=$noisetype  seed=$seed	


exit

	

