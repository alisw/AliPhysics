#!/bin/bash

# --- Make a message -------------------------------------------------
msg() { 
    local lvl=$1 ; shift 
    if test $lvl -le $verb ; then 
	echo -e "$@" >/dev/stderr
    fi 
}

# --- Extract title --------------------------------------------------
extractTitle()
{
    local dir=$1 
    if test -f ${dir}/.doc ; then 
	grep -v -i '^long:' ${dir}/.doc | sed 's/^[sS]hort://'
    else
	echo $dir
    fi
}
# --- Extract description --------------------------------------------------
extractDescription()
{
    local dir=$1 
    if test -f ${dir}/.doc ; then 
	grep -v -i '^short:' ${dir}/.doc | sed 's/^[lL]ong://'
    else
	echo $dir
    fi
}


# --- Loop dir -------------------------------------------------------
row=1
loopDir()
{
    
    # --- Our inputs -------------------------------------------------
    local level=$1 ; shift 
    local here=$1 ; shift 
    local parent=$1 ; shift

    # --- Check that we want to process this dir ----------------------
    if test $level -ge $maxCol ; then return ; fi
    local basehere=`basename $here`
    case x$basehere in 
	xqa|xjsRoot*) return ;; 
    esac

    if test -f ${here}/.nosub ; then return ; fi
    
    msg 0 "Processing @ level $level: $here ($basehere)" 

    # --- Loop 1st pass ----------------------------------------------
    local subid=1
    for sub in ${here}/* ; do 
	msg 1 "Processing $sub @ level $level" 
	case $sub in 
	     *jsRoot*) continue ;; 
	esac
	# --- If this is not a sub-directory, go on ------------------
	if test ! -d $sub ; then continue ; fi 
	

	# --- Get some stuff for sub-dir -----------------------------
	subshort=`extractDescription $sub` 
	if test "x$subShort"  = "x" ; then 
	    subshort=`extractTitle $sub` 
	fi
	subbase=`basename $sub` 
	subtitle=`echo $subbase | sed 's/^0*//'`
	local subspan=$maxCol
	local subpar=${parent}${subid}
	local subvis="table-row"
	if test $level -gt 0 ; then 
	    subvis="none"
	fi
	let subspan=$subspan-$level

	# --- Write out our row --------------------------------------
	cat <<EOF
        <tr style='display:${subvis}; cursor: move' 
            id='r$subpar' 
            onClick='toggle("r${subpar}",false);';>
EOF
	
	# --- Put in fillers -----------------------------------------
	for i in `seq 1 $level` ; do 
	    echo "        <td class='filler'>&nbsp;</td>"
	done

	# --- Now check what we should put in front ------------------
	local subpre="&nbsp;"
	local sublvl=$level
	let   sublvl=$sublvl+1
	if test $sublvl -lt $maxCol ; then 
	    hasSubSub=0
	    for subsub in ${sub}/* ; do 
		if test -d $subsub ; then hasSubSub=1 ; break ; fi
	    done
	    if test $hasSubSub -gt 0 ; then 
		subpre="&blacktriangleright;"
	    fi
	fi
	
	# --- Put in our cells ---------------------------------------
	cat <<EOF
        <td colspan='$subspan' id='n${subpar}' class='name'>
          <span style='width: 4em' id='a${subpar}'>${subpre}</span> 
EOF
	if test $link -gt 0 && test -f $sub/index.html 
	    # && test $level -ge $maxCol; 
	then 
	    echo "        <span class='lnk' onclick='showSub(\"$sub\")'>$subtitle</span>"
	    # echo "        <a href='$sub/index.html'>$subbase</a>"
	else
            # echo "        $subbase"
	    echo "        $subtitle"
	fi
	echo "</td>"
	if test $size -gt 0 ; then
	    subsize=`du ${duopt} -s $sub | cut -f1` 
	    echo "<td class='subsize'>$subsize</td>"
	fi 
	echo -e "<td id='desc$row' class='desc'>$subshort</td>\n     </tr>"

	let   row=$row+1
	let   subid=$subid+1
	loopDir "$sublvl" "$sub" "${subpar}."
    done
}

# --- Help -----------------------------------------------------------
usage()
{
    cat <<EOF
Usage: $0 [OPTIONS]

Options:
	-h,--help		 Show this help
	-m,--max-depth	 NUM	 Do not scan deeper than NUM ($maxCol)
	-t,--title       STRING  Title of page ($title)
	-d,--description STRING  Description on page ($desc)
	-o,--output      FILE    Output file ($out)
	-i,--input       DIR     Starting diretory ($inp)
	-l,--link                Link to index.html in subdirs ($link)
	-s,--size                Show size of directories
EOF
}

# === Executable code ================================================
maxCol=5
title="Local data"
desc="Collection of downloaded and locally generated ALICE data"
out=index.html
inp=.
link=0
unit=m
frame=0
size=0
base=$QA_FWD
if test "X$base" = "X" ; then 
    base=$ALICE_PHYSICS/PWGLF/FORWARD/analysis2/qa
fi 
verb=0

while test $# -gt 0 ; do
    case $1 in 
	-h|--help)        usage $0 ; exit 0 ;;
	-b|--base)        base=$2; shift ;;
	-m|--max-depth)   maxCol=$2 ; shift ;; 
	-t|--title)       title="$2" ; shift ;; 
	-d|--description) desc="$2" ; shift ;;
	-o|--output)      out="$2" ; shift ;; 
	-i|--input)       inp="$2" ; shift ;; 
	-u|--unit)        unit=`echo $2 | tr '[a-z]' '[A-Z]'` ; shift ;;
	-l|--link)        link=1 ;;
	-f|--frame)       frame=1 ;;
	-s|--size)        size=1 ;;
	--no-size)        size=0 ;;
	-v|--verbose)     let verb=$verb+1 ;; 
	*) echo "$0: Unknown option '$1'" > /dev/stderr ; exit 1;; 
    esac
    shift
done

case $unit in 
    M) ut="(MB)"; duopt="-BM" ;; 
    G) ut="(GB)"; duopt="-BG" ;;
    h) ut=""    ; duopt="-h"  ;;
    K) ut="(kB)"; duopt="-BK" ;;
    *) ut="(?)" ; duopt="-h"  ;;
esac
cat <<EOF > ${out}
<!DOCTYPE html>
<html>
  <head>
    <title>$title</title>
    <link rel='stylesheet' href='style.css'>
    <link rel='shortcut icon' href='fmd_favicon.png' 
          type='image/x-png'>
    <script type="text/javascript" src="script.js"></script>
  </head>
  <body>
    <h1 style='vertical-align: middle'>
      <img style='width: 100px;'
           src='fmd_logo.png'>
      $title
    </h1>
    <p>
      $desc
    </p>
    <div id="nav">
      <table>
        <tr>
          <th colspan="$maxCol" style='min-width:300px'>Directory</th>
EOF
if test $size -gt 0 ; then 
    echo "          <th>Size $ut</th>" >> ${out}
fi
cat <<EOF >> ${out}
          <th>Description</th>
        </tr>
EOF
loopDir 0 "${inp}" "" >> ${out}
cat <<EOF >> ${out}
        <tr style='border-top:thin solid gray'>
          <td colspan="$maxCol"></td>
EOF
if test $size -gt 0 ; then 
    totalSize=`du ${duopt} -s ${inp} | cut -f1`
    echo "           <td class='subsize'>$totalSize</td>"
fi
date=`date`
let md=$maxCol-1
cat <<EOF >> ${out}
         <td></td>
        </tr>
      </table>
      <p>
        <button onClick='hideAll();'>Collapse all</button>
        <button onClick='expandAll();'>Expand levels:</button>
        <output id="currentMax" for="maxExpand">2</output>  
        <input type="range" id="maxExpand" min="1" max="$md" value="2" list="maxList"
               onchange="currentMax.value=this.value">
        <datalist id="maxList">
EOF
for lvl in `seq 1 $md` ; do 
cat <<EOF >> $out
          <option>$lvl</option>
EOF
done
cat <<EOF >> ${out}
        </datalist>
      </p>
      <div class='change'>Last update: ${date}</div>
    </div>
EOF
if test $frame -gt 0 ; then 
    cat <<EOF >> ${out}
    <div id="frame">
       <div id="close" onclick="closeDisplay()">Close</div>
       <!-- <div id="dframe"> -->
         <iframe name="display" id="iframe"></iframe>
       <!-- </div> -->
    </div>
EOF
fi
cat <<EOF >> ${out}
  </body>
</html>
EOF
if test ! -d $base ; then 
    exit 0
fi
cp $base/style.css .
cp $base/script.js . 
cp $base/fmd_favicon.png .
cp $base/fmd_logo.png . 
#
# EOF
#

