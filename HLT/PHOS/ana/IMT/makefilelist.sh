files=/tmp/data2/digits*/*root
filecnt=0;

rm  filelist.txt

#for file in $files; do
#    echo $file
#    $((filecnt++))
#    echo   $((filecnt))
# done

for file in $files; do
    $((filecnt++))
done

echo $((filecnt)) > filelist.txt

for file in $files; do
    echo $file >> filelist.txt
done