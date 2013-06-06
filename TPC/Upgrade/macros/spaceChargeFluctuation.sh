#
# shell scipt to 
#
# arument 1 -  path to the aliroot iinitialization script
#
source $1 
aliroot -b -q $ALICE_ROOT/TPC/Upgrade/macros/spaceChargeFluctuation.C+
exit;



prefix="/hera/alice/local/filtered/alice/data/"
wdir=`pwd`
for a in `cat rawAll.list`; do
    dname=`echo $a| sed s_"$prefix"__g| sed s_"/"_"\_"_g  | sed s_".root"__`
    echo $a $dname
    mkdir dname
done;

