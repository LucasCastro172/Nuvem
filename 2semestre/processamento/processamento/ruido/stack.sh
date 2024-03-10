#! /bin/sh

#------------------------------------------------
# Empilhamento de dados
#------------------------------------------------

sustack <nmo_rui.data >stack_rui.data
supsimage < stack_rui.data > stack_rui.eps
supsimage <stack_rui.data >empilhamento.ps \
title="Stacked Data" label1="Two Way Traveltime [s]" \
wbox=6 hbox=2 perc=95 x2beg=-0.735 d2=0.025 f2=-0.735 &
exit
