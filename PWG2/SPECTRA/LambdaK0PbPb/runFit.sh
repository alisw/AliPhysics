
cent=00
# Set hist equal to 0 for standard pt analysis, other values for QA histograms
hist=0
#suffix=K0MC
#suffix=K0Real
#suffix=LMC
suffix=LReal
# Set partID = 0 for K0s, 1 for Lambda, 2 for anti-Lambda
# or 3 for special case of summed anti-Lambda + Lambda 
partID=1
# Set data equal to the subdirectory within output[or output10binsNew - check in run.sh] which contains root files from task
#data=LHC10h_000139172_p2
#data=LHC10h_pass2
#data=LHC11a10a_final
data=Twiki23Aug
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
echo "output directory is: " ${outfold}
