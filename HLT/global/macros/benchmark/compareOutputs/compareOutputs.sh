#! /bin/bash
awk '
	BEGIN{
		counted=0
		current=""
	} 
	{
		if( $1 ~ "^_" ){
			current = substr($1,2)
			classes[current]++
			counted = 0
		}
		else{
			if( $NF!=0 ){
				if( counted == 0 ){
					conflicts[current][0] ++
					counted = 1
				}
				conflicts[current][$1] ++
				print $0 >> "conflicts_" current ".txt"
			}
		}
	}
	END{
		for(c in classes){
			if( conflicts[c][0] != 0 ){
				print c ":  " conflicts[c][0] " conflicts:"
				for( f in conflicts[c] ){
					if( f != 0 )
						print "    " f ": " conflicts[c][f]
				}
			}
			else{
				print c ": no conflicts"
			}
			print ""
		}
	} 
' < ${1:-"comparison.txt"}