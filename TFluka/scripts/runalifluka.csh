rm -rf tmp
mkdir tmp
cd tmp
ln -s $FLUPRO/xnloan.dat .
ln -s $FLUPRO/sigmapi.bin .
ln -s $FLUPRO/nuclear.bin .
ln -s $FLUPRO/neuxsc_72.bin neuxsc.bin
#cp ../fort.16 .
#cp $FLUPRO/random.dat ranmu001
#ln -s ranmu001 fort.1
#ln -s ranmu002 fort.2
cp $FLUPRO/random.dat old.seed
ln -s mu001.out fort.11
#ln -s $FLUPRO/libec_thihecufealw_10t.pemf fort.12
#ln -s ../alice.pemf fort.12
ln -s ../alice.pemf .
#ln -s mu001.err fort.15
ln -s $FLUPRO/fluodt.dat .
ln -s $FLUPRO/elasct.bin .
ln -s ../muon.inp .
aliroot
cd ..
