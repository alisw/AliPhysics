#!/bin/bash
# List all CDB paths present in the current OCDB folder
# The script is meant to check the different CDB path in the reference OCDB
# ($ALICE_ROOT) to make sure all entries are copied to the destination OCDB
# folder. In particular we are interested in the second of the three levels
# composing the CDB-path since usually the copy from one OCDB to the other
# is done by lines of the type
# TList* l = cdb->GetAll("*/Calib/*");
# TList* l = cdb->GetAll("*/Align/*");
# (see macro CDBToGrid.C)
#

stdfname="Run0_999999999"
msg="List of second level strings for CDB entries:"
cd $1
for i in `find ./ -name 'Run*' | grep -v svn`
do
	fname=`echo "$i" | cut -d/ -f5`
	if(expr match "$fname" "Run0." > /dev/null)
	then
		if [[ $fname == Run0_99999999* ]]
		then
			second=`echo "$i" | cut -d/ -f3`
			if [[ $second != Align ]] && [[ $second != Calib ]] && [[ $second != Config ]]
			then
				#printf "$second in $i\n"
				msg="$msg
				$second for entry $i" 
			fi
		else
			echo "Warning: $i is not a good file name for OCDB"
		fi
	fi
done
echo "$msg"

