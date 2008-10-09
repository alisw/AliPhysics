#/bin/bash

if [ $# -ne 1 ]; then
        echo 1>&2 Usage: $0 RCU 
        echo 1>&2 Example: $0 1
        exit 127
fi

#what has been asked for?
RCU=$1


BRANCH="A"
date 
echo Branch $BRANCH
FEC=1
while [ $FEC -le 9 ]; do
 FILE=APD/set_rcu_${RCU}_bias_branch_${BRANCH}_FEC_${FEC}.scr
 echo
 echo File $FILE will be used
 echo
 rcu-sh b $FILE

 let FEC=FEC+1
done

BRANCH="B"
date 
echo Branch $BRANCH
FEC=1
while [ $FEC -le 9 ]; do
 FILE=APD/set_rcu_${RCU}_bias_branch_${BRANCH}_FEC_${FEC}.scr
 echo
 echo File $FILE will be used
 echo
 rcu-sh b $FILE

 let FEC=FEC+1
done

#fini
date
echo
echo That is all. Fin.
echo