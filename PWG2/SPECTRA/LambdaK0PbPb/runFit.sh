
cent=04
hist=2
suffix=K0MC
#suffix=K0Real
#suffix=LMC
#suffix=LReal
partID=0
#data=LHC10h_000139172_p2
#data=LHC10h_pass2
data=LHC11a10a_final
out=results
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
echo "output directory is: " ${out fold}
