#!/bin/bash 

for i in $* ; do 
    t=`file $i` 
    case $t in 
	*ROOT*) continue ;; 
	*Zip*) mv $i ${i}.zip ;; # fall - through 
	*) mv $i $i.broken ; continue ;; 
    esac 
    
    d=`mktemp -d` 
    unzip ${i}.zip -d ${d}
    f=`echo ${i} | sed 's/_[0-9][_0-9]*\.root/.root/'` 
    if test -f ${d}/${f} ; then 
	mv ${d}/${f} ${i}
	echo "Extracted ${f} to ${i}"
    else 
	echo "${f} not found in ${i}.zip" 
	unzip -l ${i}.zip 
    fi 
    rm -rf ${d}
done 

	
	    
