
cent=10
hist=0
#suffix=K0MCb
#suffix=K0MCa
suffix=K0Real
#suffix=LMCb
#suffix=LMCa
#suffix=LReal
#suffix=AntiLReal
#suffix=Test
partID=0
#data=LHC10h_000139507_p2
#data=LHC10h_000139310_p2
#data=LHC10h_000139314_p2
#data=LHC10h_000139328_p2
#data=LHC10h_000139466_p2
#data=old/LHC10h_000139507_p2
#data=LHC10hAll/final/
#data=linkData
#data=NoPIdoldAll/
#data=LHC10h_000139507_p2
#data=noPid/LHC10h_000139310_p2
#data=LHC10h_000139507_p2
#data=LHC11a10bAll/final
#data=LHC11a10a_000139510
data=LHC10hAll/final/
#data=LHC11a10bAll/final/
#data=LHC11a10aAll/final/
#data=LHC11a10b_plus_000139466
#data=LHC10h_000137718_p2/
out=results/
#out=resultsPt


#./run.sh -f LHC10h_000139172_p2 -b 00 -s 4 -x K0cent00MC -m -p 0
./run.sh -f ${data} -b ${cent} -s ${hist} -x ${suffix} -m -p ${partID}

mkdir ./${out}
mkdir ./${out}/${suffix}/
mkdir ./${out}/${suffix}/${cent}
outfold=./${out}/${suffix}/${cent}
mv Yield* ${outfold}
mv Diag* ${outfold}
mv Masses* ${outfold}

echo "//////////////////////////////////////"
echo "output directory is: " ${outfold}
