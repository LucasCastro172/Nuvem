#! /bin/sh

#------------------------------------------------
# Stack Data...
#------------------------------------------------

sustack <astur_nmo.data >astur_stack.data
supsimage < astur_stack.data > astur_stack.eps
supsimage <astur_stack.data >empilhamento.ps \
title="Stacked Data" label1="Two Way Traveltime [s]" \
wbox=6 hbox=2 perc=95 x2beg=-0.735 d2=0.025 f2=-0.735 &
exit
