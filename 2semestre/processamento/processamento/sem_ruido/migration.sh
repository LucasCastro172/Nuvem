#! /bin/sh

#------------------------------------------------
# Migrate Stacked Data...
#------------------------------------------------

sumigps <astur_stack.data >astur_migration.data \
tmig=0.0 \
vmig=2000 dx=50
supsimage < astur_migration.data > Migração_dos_dados.eps
exit
