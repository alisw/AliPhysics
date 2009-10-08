#
# marian.ivanov@cern.ch
#
# Make workspace structure
# Create a list for each run 
# and make directory structure
# This is fast
mydir=`pwd`
for adir in `cat run.list`; do
    mkdir $adir;
    rm raw${adir}.txt
    rm esd${adir}.txt
    cat  $mydir/raw.list | grep $adir >raw${adir}.txt
    cat  $mydir/esd.list | grep $adir >esd${adir}.txt
    cp raw${adir}.txt   ${adir}/raw.txt
    cp esd${adir}.txt   ${adir}/esd.txt
done 


