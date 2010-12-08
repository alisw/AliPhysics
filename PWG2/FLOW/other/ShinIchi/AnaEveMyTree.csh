#! /bin/csh -f 

source root.setup
@ ii = 0
# foreach run (137549)
# foreach file (1 4 5 6 7 8 9 10 11 12 14 15 16 18 19 23 24 25 26 27 28 29 30 31 32 33 35 37 38 39 40 42 43 44 45 46 47 49 51 52 53 57 58 60 62 64 65 66 67 68 69 70 71 72 73 74 75)
foreach run (137431)
foreach file (1 2 9 11 27 28 39 40 41 99)
foreach nn (1 2 3 4 5)
  @ ii = $ii + 1
  echo "++++++++++++++++++++++++++++++++++++++++++++++"
  echo "+++++" num $ii run $run file $file step $nn "+++++"
  echo "++++++++++++++++++++++++++++++++++++++++++++++"
  date
  echo "++++++++++++++++++++++++++++++++++++++++++++++"

root.exe -b > AnaEveMyTree_${nn}.log <<EOF
  .L AnaEveMyTree.C
  AnaEveMyTree a(${run},${file})
  a.Loop()
  .q
EOF

end
mv AnaEveMyTree.root   root${run}/AnaEveMyTree_${file}.root
mv AnaEveMyTree.cal     cal${run}/AnaEveMyTree_${file}.cal
mv AnaEveMyTree_1.log   log${run}/AnaEveMyTree_${file}_1.log
mv AnaEveMyTree_2.log   log${run}/AnaEveMyTree_${file}_2.log
mv AnaEveMyTree_3.log   log${run}/AnaEveMyTree_${file}_3.log
mv AnaEveMyTree_4.log   log${run}/AnaEveMyTree_${file}_4.log
mv AnaEveMyTree_5.log   log${run}/AnaEveMyTree_${file}_5.log
end
end

echo "++++++++++++++++++++++++++++++++++++++++++++++"
date
echo "++++++++++++++++++++++++++++++++++++++++++++++"
echo finished
echo "++++++++++++++++++++++++++++++++++++++++++++++"
