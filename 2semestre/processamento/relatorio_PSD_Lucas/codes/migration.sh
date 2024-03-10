#! /bin/sh

#------------------------------------------------
# Migrate Stacked Data...
#------------------------------------------------

sumigps <stack_rui.data >migration_rui.data \
tmig=0.0 \
vmig=2000 dx=50
supsimage < migration_rui.data > Migracao_dos_dados.eps
exit
