#! /bin/sh
#________________________________________________________________________________________
# O objetivo deste programa é de fazer análise de velocidades.
# Para determinar os demais cmp usamos uma progressão aritmética: a_{n}=a_{1} + (n-1)*r 
#________________________________________________________________________________________
# set -x

cmps=ruido_cdp.su 

if [ ! -f $cmps ]
then    echo "Sort to CMP first!"
        pause EXIT
        exit
fi   # encerra "if"

echo "Analise de Velocidades"

rm -f panel.* picks.* par.* tmp*

#------------------------------------------------
# Definicao de variaveis:
#------------------------------------------------

indata=$cmps
outdata=velan_ruido.data #vpick

nt=751
dt=0.004

nv=420      # Numero de velocidades
dv=10       # Intervalo
fv=1500     # Primeira velocidade

>$outdata   # Arquivo vazio
>par.cmp    # Arquivo vazio

#------------------------------------------------
# Analise de velocidade Interativo
#------------------------------------------------

echo "Entre o numero de CMPs para o Picking:" >/dev/tty
read nrpicks

i=1
while [ $i -le $nrpicks ]
do
    echo "Especifique o numero de CMP $i" >/dev/tty
    read picknow
    echo "Preparando a posicao $i de $nrpicks para Picking "
    echo "A posicao e CMP $picknow "

#------------------------------------------------
# Plotando o CDP:
#------------------------------------------------

    suwind <$indata key=cdp min=$picknow \
            max=$picknow >panel.$picknow
    suxwigb <panel.$picknow xbox=312 ybox=10 \
             wbox=300 hbox=600 \
             title="Secao CMP $picknow" \
             perc=94 key=offset verbose=0 &

#------------------------------------------------
# Preparando o CVS (Constant Velocity Stack)
#------------------------------------------------

    >tmp1
    j=1
    k=`expr $picknow + 10`
    l=`bc -l <<-END
    $dv * $nv / 120
    END`

    suwind <$indata key=cdp min=$picknow \
            max=$k >tmp0

    while [ $j -le 10 ]
    do
	vel=`bc -l <<-END
    	$fv + $dv * $j * $nv / 10
    	END`

	sunmo <tmp0 vnmo=$vel |

	sustack >>tmp1
        sunull ntr=2 nt=$nt dt=$dt >>tmp1
        j=`expr $j + 1`
    done
    suximage <tmp1 xbox=834 ybox=10 wbox=400 hbox=600 \
	     title="Constant Velocity Stack CMP $picknow" \
	     label1="Tempo [s]" label2="Velocidade [m/s]" \
	     f2=$fv d2=$l verbose=0 mpicks=picks.$picknow \
	     perc=99 n2tic=5 cmap=rgb0 &

#------------------------------------------------
# Plotando o mapa de Semblance:
#------------------------------------------------

    echo "  Coloque o ponteiro do mouse sobre o mapa de Semblance ou sobre"
    echo "  os paineis CVS e pressione 's' para salvar a velocidade marcada." 
    echo "  Pressione 'q' no mapa de semblance depois de marcar todas as velocidades"

    suvelan <panel.$picknow nv=$nv dv=$dv fv=$fv |
    suximage xbox=1 ybox=10 wbox=300 hbox=600 \
	     units="Semblance" f2=$fv d2=$dv \
	     label1="Tempo [s]" label2="Velocidades [m/s]" \
	     title="Mapa de Semblance - CMP $picknow" cmap=hsv2 \
	     legend=1 units=Semblance verbose=0 gridcolor=black \
	     grid1=solid grid2=solid mpicks=picks.$picknow
    
    sort <picks.$picknow -n |
    mkparfile string1=tnmo string2=vnmo >par.$i

#------------------------------------------------------------
# Plotando a seção corrigida de NMO e o perfil de velocidades
#------------------------------------------------------------

    >tmp2
    echo "cdp=$picknow" >>tmp2
    cat par.$i >>tmp2

    sunmo <panel.$picknow par=tmp2 |
    suxwigb title="Secao CMP depois do NMO" xbox=1 ybox=10 \
	     wbox=300 hbox=600 verbose=0 key=offset perc=94 &

        sed <par.$i '
	s/tnmo/xin/
	s/vnmo/yin/
	        ' >par.uni.$i
    unisam nout=$nt fxout=0.0 dxout=$dt \
	   par=par.uni.$i method=mono |
    xgraph n=$nt nplot=1 d1=$dt f1=0.0 \
           label1="Tempo [s]" label2="Velocidade [m/s]" \
	   title="Funcao Velocidade CMP $picknow" \
           -geometry 200x600+934+10 style=seismic &

    echo "Deseja salvar as velocidades marcadas? (y/n) " >/dev/tty
    read response

    rm tmp*

    case $response in
	n) i=$i echo "Removendo as velocidades marcadas" ;;
        y) i=`expr $i + 1`
           echo "$picknow  $i" >>par.cmp ;;
        *) echo "Digite y para yes e n para nao";;
    esac

done

#------------------------------------------------
# Criando o arquivo de saida para as velocidades
#------------------------------------------------

mkparfile <par.cmp string1=cdp string2=# >par.0

i=0
while [ $i -le $nrpicks ]
do
	cat par.$i >>$outdata
	i=`expr $i + 1`
done

rm -f panel.* picks.* par.* tmp*

exit

