#!/bin/sh

#Vai plotar imagens dos shots para verificar se está ok a aquisição.


infile=camadas_planas.su
#temp1=camada_plana_cdp1.su
outfile=camadas_planas_cdp.su

#Inserir a informações referentes aos cdps
suchw<$infile key1=cdp key2=gx key3=sx b=1 c=1 d=2 > temp1

#Renumeração dos cdps
suchw <temp1 key1=cdp key2=cdp key3=cdp a=-1237 b=1 d=25 | susort > $outfile cdp offset

#rm temp*

sugethw <$outfile key=tracf,fldr,sx,gx,offset,cdp > info.txt


sukeycount <$outfile key=cdp> fold.txt


exit

