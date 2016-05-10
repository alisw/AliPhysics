#!/bin/bash
#
# Script to download generated AODs in parallel using GNU parallel [1]. 
# GNU parallel must be installed on the system.
#
# [1] http://www.gnu.org/software/parallel/
#
# --- Some global variables ------------------------------------------
#
# Delay between downloads in seconds 
delay=0.1
# Storage elements to skip 
skip=

# --- Download function ----------------------------------------------
#
# Function to download a single file and store it with a unique name
# based on the parent and grand-parent directory names on the Grid.
#
copyit()
{
    local src=$1
    local dir=`dirname $1`
    local base=`basename $1 .root`
    local par=`basename $dir`
    local grnd=`basename $(dirname $dir)` 
    local dest=${base}_${grnd}_${par}.root
    # mkdir -p ${par}
    if test -f $dest ; then
	# echo "$src -> $dest ... exists"
	return
    fi
    # Below can be used to filter on specific locations
    local at=
    if test "x$skip" != "x" ; then 
	local where=($(alien_whereis $dir/root_archive.zip 2>/dev/null | \
		      sed -n 's/[[:space:]]*SE => \([^[:space:]]*\) .*/\1/p'))
	for s in $skip ; do  
	    if test "x${where[0]}" = "x${s}" ; then
		at="@${where[1]}"
		break
	    fi
	done
    fi 
    alien_cp -s alien:${src}${at} file:${dest} >> download.log 2>&1 
    sleep $delay
}
# Export the function
export -f copyit

# --- Help function --------------------------------------------------
usage()
{
    cat <<-EOF
	Usage: $0 -s START -e END -f FILENAME -d TOP_DIRECTORY
	
	Options:
	  -s,--start	NUMBER		Offset in list [$start]
	  -e,--end	NUMBER		Last file number [$stop]
	  -f,--file	FILENAME	File to get [$file]
	  -d,--top	DIRECTORY	Grid top directory [$dir]
	  -j,--jobs	JOBS		Number of parallel jobs [$jobs]
	  -D,--delay    SECONDS		Delay between downloads [$delay]
	  -S,--skip	SE-NAME		Storage elements to skip [$skip]
	  -h,--help			This help
	  -z,--zenity			Use zenity (GUI) to show progress
	EOF
}

# --------------------------------------------------------------------
#
# Process the command line
#
start=0
stop=10000
file=
dir=
zen=
jobs=4
while test $# -gt 0 ; do
    case $1 in
	-h|--help)	usage 		; exit 0 ;;
	-s|--start)	start=$2 	; shift ;;
	-e|--end)	stop=$2  	; shift ;;
	-f|--file)	file=$2  	; shift ;;
	-d|--top)	dir=$2   	; shift ;;
	-j|--jobs)	jobs=$2		; shift ;; 
	-D|--delay)	delay=$2 	; shift ;;
	-z|--zenity)	zen=1		;;
	*) echo "Unknown option: $1" > /dev/stderr ; exit 1 ;;
    esac
    shift
done

# --------------------------------------------------------------------
#
# Sanity checks on arguments
#
if test $start -gt $stop ; then
    echo "Invalid range $start < $stop" > /dev/stderr 
    exit 1
fi
if test "x$file" = "x" ; then
    echo "No pattern given" > /dev/stderr
    exit 1
fi
if test "x$dir" = "x" ; then
    echo "No directory given" > /dev/stderr
    exit 1    
fi
if [[ $dir =~ ^alien: ]] ; then
    echo "No need to specify protocol"
    dir=`echo $dir | sed -e 's,^alien:,,' -e 's,^/*\(/\),\1,'`
fi
if ! [[ $dir =~ ^/ ]] ; then
    # if not absolute, assume it's in users home directory
    ausr=`alien_whoami | tr -d '[[:space:]]'`
    alet=`echo $ausr | sed 's/^\(.\).*/\1/'` 
    ahme="/alice/cern.ch/user/${alet}/${ausr}"
    dir=${ahme}/${dir}
    echo "Not absolute, path assuming in users home: $ahme"
fi
if ! alien_ls $dir > /dev/null 2>&1 ; then
    echo "Directory $dir does not exist" > /dev/stderr
    exit 1
fi

# --- Derived variables ----------------------------------------------
suf=`basename $file | sed 's/.*\.//'`
n=$(($stop-$start))
total=`alien_find $dir $file | tail -n 1 | sed 's/ *files found//'`

# --- Show settings --------------------------------------------------
cat <<EOF
Directory:	${dir}
File name:	${file}
Available:	${total}
First:		${start}
Second:		${stop}
To get:		${n}
SEs to skip:	${skip}
Zenity GUI:	${zen}
Jobs:		${jobs}
EOF

# --- Now run it -----------------------------------------------------
if test $zen -gt 0 ; then 
    alien_find $dir $file | \
	grep .${suf} | \
	tr -d ' ' | \
	head -n $stop | \
	tail -n $n | \
	parallel --eta --bar -j ${jobs} copyit {} \
		 2> >(zenity --width=700 --height=100 --progress --auto-close --time-remaining --no-cancel)
else 
    alien_find $dir $file | \
	grep .${suf} | \
	tr -d ' ' | \
	head -n $stop | \
	tail -n $n | \
	parallel --eta --bar -j ${jobs} copyit {}
fi 

# --- End of job -----------------------------------------------------
#
# EOF
#



