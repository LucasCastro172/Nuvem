#! /bin/sh
#­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­
# Aplica a correcao NMO nos dados CMP ...
#­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­
sunmo <asturseis_cmp.su astur_velan.data smute=1000 >astur_nmo.data
supsimage < astur_nmo.data > correção-NMO-nos-dados-CMP.eps
exit
