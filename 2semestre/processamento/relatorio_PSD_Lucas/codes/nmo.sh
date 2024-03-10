#! /bin/sh
#­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­
# Aplica a correcao NMO nos dados CMP ...
#­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­­
sunmo <ruido_cdp.su velan_ruido.data smute=1000 >nmo_rui.data
supsimage < nmo_rui.data > correção-NMOruido.eps
exit
