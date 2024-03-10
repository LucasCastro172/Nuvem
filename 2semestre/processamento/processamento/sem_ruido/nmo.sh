#! /bin/sh
#­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­
# Aplica a correcao NMO nos dados CMP ...
#­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­
sunmo <camadas_planas_cdp.su astur_velan.data smute=1000 >astur_nmo.data
supsimage < astur_nmo.data > correção-NMO-nos-dados-CMP.eps
exit
