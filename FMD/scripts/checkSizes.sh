#!/bin/bash

extra=" DPMJET		\
        TDPMjet		\
        EPEMGEN		\
        TEPEMGEN	\
        HBTP		\
        THbtp		\
        HERWIG		\
        THerwig		\
        HIJING		\
        THijing		\
        ISAJET		\
        TIsajet		\
        LHAPDF		\
        MEVSIM		\
        TMEVSIM		\
        MICROCERN	\
        PDF		\
        PYTHIA6		\
        TPHIC"
base="  ALIFAST		\
	ALIROOT		\
	ANALYSIS	\
	CONTAINERS	\
        CRT		\
	DISPLAY		\
	EMCAL		\
	EVE		\
        EVGEN		\
	FASTSIM		\
	FLOW		\
	FMD		\
        HBTAN		\
	HLT		\
	ITS		\
	JETAN		\
        LHC		\
	MONITOR		\
	MUON		\
	PHOS		\
        PMD		\
	PWG0		\
	PWG2		\
	PWG3		\
        RALICE		\
	RAW		\
	RICH		\
	SHUTTLE		\
        START		\
	STEER		\
	STRUCT		\
	TOF		\
        TPC		\
	TRD		\
	VZERO		\
	ZDC"

cat <<EOF > exclude
*/tgt_*/*
*/html/*
.#*
*/CVS*
*~ 
*.root
*.so
*.o
EOF

get_size()
{
    s=`du -X exclude -kc $1 | tail -n 1 | awk 'BEGIN {FS=" "}{print $1}'`
    printf "\t%-30s\t%10d kB\n" $1 $s 
    total=`echo ${total} + ${s} | bc`
}
    
echo "Extras:"
total=0
for e in $extra ; do 
    get_size $e
done 
for i in `seq 1 56` ; do echo -n "-" ; done
mb=`echo $total / 1024 | bc` 
printf "\n\t%-30s\t%10d kB = %10d MB\n" "Total" $total $mb

echo "Base:"
total=0
for b in $base ; do 
    get_size $b 
done 
for i in `seq 1 56` ; do echo -n "-" ; done
mb=`echo $total / 1024 | bc` 
printf "\n\t%-30s\t%10d kB = %10d MB\n" "Total" $total $mb

        
rm -f exclude 