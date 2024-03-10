#! /bin/sh
# File: acq.sh
# Set messages on
##set -x

# Assign values to variables
nangle=131 
fangle=-65 
langle=65 
nt=751 
dt=0.004 

# Model
#num=4
#echo " --Model number = $num"

# Name input model file
inmodel=camadas_planas.dat

# Name output seismic file
outseis=camadas_planas.su

# Remove survey file
rm -f survey.txt
# Name survey file
survey=survey.txt

#=================================================
# Create the seismic traces with "triseis"
#   i-loop = 40 source positions
#   j-loop = 60 geophone positions (split-spread)
#            per shot position
#   k-loop = layers 2 through 24
#            (do not shoot layers 1 and 25)

echo " --Begin looping over triseis."

i=0
while [ "$i" -ne "40" ]
do
	fs=`bc -l <<-END
	$i * 0.05 
	END`
	sx=`bc -l <<-END
	$i * 50 + 2000
	END`
	fldr=`bc -l <<-END
	$i + 1
	END`
		j=0
		while [ "$j" -ne "60" ]
 		do
   		fg=`bc -l <<-END
    		$i * 0.05 + $j *0.05 
		END`
		gx=`bc -l <<-END
		$i * 50 + $j * 50 + 525
		END`
		offset=`bc -l <<-END
		$j * 50 - 1475
		END`
		tracl=`bc -l <<-END
		$i * 60 + $j + 1
		END`
		tracf=`bc -l <<-END
		$j + 1
		END`
		echo " Sx=$sx   Gx=$gx   fldr=$fldr   Offset=$offset   tracl=$tracl\
		fs=$fs   fg=$fg"
		echo " Sx=$sx   Gx=$gx   fldr=$fldr   Offset=$offset   tracl=$tracl\
		fs=$fs   fg=$fg" >> $survey
			k=2
			while [ "$k" -ne "4" ]
			do
			triseis < $inmodel  xs=2,3.95 xg=0.525,5.425  zs=0,0 zg=0,0 \
			nangle=$nangle fangle=$fangle langle=$langle \
			kreflect=$k krecord=1 fpeak=40 lscale=0.5 \
			ns=1 fs=$fs ng=1 fg=$fg nt=$nt dt=$dt |
			suaddhead nt=$nt | sushw key=dt,tracl,tracr,fldr,tracf,trid,offset,sx,gx \
			a=4000,$tracl,$tracl,$fldr,$tracf,1,$offset,$sx,$gx >> temp$k
			k=`expr $k + 1`
			done
		j=`expr $j + 1`
		done
i=`expr $i + 1`
done

echo " --End looping over triseis."

#=================================================

# Sum contents of the temp files
echo " --Sum files."
susum temp2 temp3 > $outseis 

# Remove temp files
echo " --Remove temp files."
rm -f temp*

# Report output file
echo " --Output file   ** $outseis **"

# Exit politely from shell script
echo " --Finished!"

sugethw <camadas_planas.su key=tracl,fldr,tracf,sx,gx,offset> teste.txt

exit
