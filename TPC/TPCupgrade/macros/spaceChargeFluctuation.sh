
#
# shell scipt to 
#
# argument 1 -  path to the aliroot iinitialization script
# argument 2 -  action argument
#               see   $ALICE_ROOT/TPC/Upgrade/macros/spaceChargeFluctuation.C
#                if (mode==0) GenerateMapRawIons(arg0);  
#                if (mode==1) DoMerge();  
#                if (mode==2) spaceChargeFluctuationToyMC(arg0,arg1);
# argument 3 -  according the C++ code
#export flucPath=$HOME/AliRoot/TPCdev/TPC/Upgrade/macros
              
source $1 
aliroot -b -q $flucPath/NimStyle.C $ALICE_ROOT/TPC/Upgrade/macros/spaceChargeFluctuation.C+\($2,$3,$4,$5\)
exit;


#
#Example usage hera (here we have 1000 nodes, but no acces to user disk)
# (to install aliroot on hera, the best is to use the rsync)
# jobs to be submitted form the lxsub0x
export baliceTPC=/hera/alice/miranov/.baliceHera
export flucPath=$HOME/AliRoot/TPCdev/TPC/Upgrade/macros/
export batchCommand="qsub -cwd -V -l h_rt=24:0:0,h_rss=4G  "
#
#
# Example usage local 
# jobs to be submitted form the lxb1001 or lxb1002
#(here we have 80 nodes and user disk)
#
export baliceTPC=/u/miranov/.baliceTPC
export flucPath=$HOME/AliRoot/TPCdev/TPC/Upgrade/macros/
export batchCommand="qsub -cwd  -V "
#
# 0.) sumbmit jobs to process raw data and accumlate space charge
#
#
prefix="/hera/alice/local/filtered/alice/data/"
wdir=`pwd`
for a in `tail -n 1900  rawAll.list`; do
    dname=`echo $a| sed s_"$prefix"__g| sed s_"/"_"\_"_g  | sed s_".root"__`
    echo $a $dname
    mkdir $wdir/$dname
    cd $wdir/$dname
    ln -sf $a raw.root
    $batchCommand    -o  filter.log  $flucPath/spaceChargeFluctuation.sh $baliceTPC  0  0 0
    cd $wdir
done;
#
# 1.) submit merging jobs example
#
ls `pwd`/*/histo.root | grep -v dirmerge  > histo.list
wdir=`pwd`
split -l 50 -d  histo.list merge
for a in `ls merge*`; do
    mkdir dir$a
    mv $a dir$a/histo.list
done;

wdir=`pwd`
for a in `ls -d dirmerge*`; do
    cd $wdir/$a
    $batchCommand    -o  filter.log  $flucPath/spaceChargeFluctuation.sh $baliceTPC  1  1 0
    cd $wdir
done;


#
# 3.)  submit fluctuation code analysis
#


wdir=`pwd`
for a in `ls -d dirmerge*`; do
    cd $wdir/$a
    $batchCommand    -o  filter.log  $flucPath/spaceChargeFluctuation.sh $baliceTPC  3  1000  10
    cd $wdir
done;


#
# 4.)  submit fluctuation code dist scan
#
rm dir*/SpaceCharg*root 
rm dir*/filter*.log
rm dir*/*.sh.e*
 
wdir=`pwd`
for a in `ls -d dirmerge* | grep -v dirmergeAll`; do
    cd $wdir/$a
    $batchCommand    -o  filterFluc0.log  $flucPath/spaceChargeFluctuation.sh $baliceTPC  4  1 0 0
    $batchCommand    -o  filterFlucP.log  $flucPath/spaceChargeFluctuation.sh $baliceTPC  4  1 0 1
    $batchCommand    -o  filterFlucM.log  $flucPath/spaceChargeFluctuation.sh $baliceTPC  4  1 0 2
    cd $wdir
done;

wdir=`pwd`
ls $wdir/dirmerge*/fluct*.root| grep -v mergeAll >  $wdir/dirmergeAll/fluctuation.list 
cd $wdir/dirmergeAll
rm dir*/SpaceCharg*root 
rm dir*/filter*.log
rm dir*/*.sh.e*

for idir in {0..100}; do
mkdir $wdir/dirmergeAll/dir$idir
cd  $wdir/dirmergeAll/dir$idir
cp $wdir/dirmergeAll/fluctuation.list $wdir/dirmergeAll/dir$idir/
for i in {0..14..2} ; do
   $batchCommand    -o  filterFluc0$i.log  $flucPath/spaceChargeFluctuation.sh $baliceTPC  4  1 $i 0
   $batchCommand    -o  filterFlucP$i.log  $flucPath/spaceChargeFluctuation.sh $baliceTPC  4  1 $i 1
   $batchCommand    -o  filterFlucM$i.log  $flucPath/spaceChargeFluctuation.sh $baliceTPC  4  1 $i 2
done;
cd $wdir
done

#
# 4.b) submit the code for the epsilon scan
#
#example directory
cd /hera/alice/miranov/SpaceCharge/Fluctuations/PbPbWithGain
wdir=`pwd`

ls $wdir/dirmerge*/fluct*.root| grep -v mergeAll >  $wdir/dirmergeAll/fluctuation.list 
cd $wdir/dirmergeAll
for epsilon in {10,20}; do
   #create and clean  directories for epsilon
   mkdirhier  $wdir/dirmergeAll/dEpsilon$epsilon
   rm -rf $wdir/dirmergeAll/dEpsilon$epsilon/*
   cp $wdir/dirmergeAll/fluctuation.list  $wdir/dirmergeAll/dEpsilon$epsilon/
done;
#submit epsilon scan jobs
for epsilon in {10,20}; do         # loop over epsilons
    for idir in {0..40}; do         # loop  create different random  ion pileup frames
       mkdir $wdir/dirmergeAll/dEpsilon$epsilon/dir$idir
       cd  $wdir/dirmergeAll/dEpsilon$epsilon/dir$idir
       cp $wdir/dirmergeAll/fluctuation.list $wdir/dirmergeAll/dEpsilon$epsilon/dir$idir/
       for i in {0..14..2} ; do         # specify differnt mulitpliicty of ions in pilepy frame run B0, B+, B-
          scaling=$(($epsilon/5))
          $batchCommand    -o  filterFluc0$i.log  $flucPath/spaceChargeFluctuation.sh $baliceTPC  4  $scaling $i 0
          $batchCommand    -o  filterFlucP$i.log  $flucPath/spaceChargeFluctuation.sh $baliceTPC  4  $scaling $i 1
          $batchCommand    -o  filterFlucM$i.log  $flucPath/spaceChargeFluctuation.sh $baliceTPC  4  $scaling $i 2
       done;
       cd $wdir
    done;
done; 




#
# 5.)  submit drawing jobs
#
wdir=`pwd`
for a in `ls -d dirmerge*`; do
    cd $wdir/$a
    $batchCommand    -o  filterFluc.log  $flucPath/spaceChargeFluctuation.sh $baliceTPC  5  0 10000 
    cd $wdir
done;


#
# 6.)  submit drawing fluctuation jobs
#
wdir=`pwd`
for a in `ls -d dir*`; do
    cd $wdir/$a
    rm localFit.root
    $batchCommand    -o  drawFlucFit.log  $flucPath/spaceChargeFluctuation.sh $baliceTPC  6  -1.5 1.5 0 0
    cd $wdir
done;

#
# 7.)  submit drawing fluctuation jobs
#
wdir=`pwd`
for a in `ls -d dir*`; do
    cd $wdir/$a
    rm localBins.root
    $batchCommand    -o  drawFlucBin.log  $flucPath/spaceChargeFluctuation.sh $baliceTPC  7  100000 10000 0 0
    cd $wdir
done;


#
# 6.) make chain files
#
for i in 0  2 4 6 8 10 12 14 ; do
   ls `pwd`/dir*/SpaceChargeTrackFluc$i\_1.root > track$i\_1.list
   ls `pwd`/dir*/SpaceChargeTrackFluc$i\_0.root > track$i\_0.list
done;
