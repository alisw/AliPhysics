#!/bin/bash


startDoc() {
    cat <<EOF
%% --- Class, packages, and theme ------------------------------------
\\documentclass[compress]{beamer}
\\usetheme{alice}
\\usepackage[english,british]{babel}
\\usepackage{pdfpages}
\\mode<presentation>
\\title{Survey of Cuts in Good and Bad run from LHC10h}
\\author{Christian Holm Christensen\\inst{1}}
\\institute{\\inst{1}Niels Bohr Institute}
\\date{\\today}
\\place{HEHI Group Meeting}
\\subject{High Energy Heavy Ion Physics}
\\setbeamertemplate{navigation symbols}{}
\\newcommand\\Lw{\\xi}
%% --- The document --------------------------------------------------
\\begin{document}

\\aliceTitlePage{}
EOF
    if test -f `dirname $out`/intro.tex ; then 
	echo "%% --- Introduction ----------- "
	echo "\\include{intro}"
    fi
    echo "%% --- The different runs ------------ "
}

endDoc()
{
    cat <<EOF
\\end{document}
EOF
}

getCutText()
{
    local w=$1  ; shift
    local n=$1  ; shift 
    local d=$1  ; shift

    local m=`echo $d | sed "s/.*_${w}\([a-z][a-z]*\).*/\1/"`
    local v=`echo $d | sed "s/.*_${w}${m}\([0-9d][0-9d]*\).*/\1/"`
    local c=`echo $v | tr 'd' '.'` 
    c=`echo $c | sed -e 's/^0*\([0-9]\)/\1/' -e 's/0*$//' -e 's/\.$//'` 
    # echo "m=$m v=$v c=$c" > /dev/stderr 
    local t=
    case $m in 
	mpv)  t=`printf "%s&=%4.2f\\Delta_{p}"                "$n" "$c"` ;; 
	xi)   t=`printf "%s&=\\Delta_{p}-%4.1f\\Lw"           "$n" "$c"` ;; 
        sig)  t=`printf "%s&=\\Delta_{p}-%4.1f(\\Lw+\\sigma)" "$n" "$c"` ;;
	prob) t=`printf "%s&: P(\\Delta<c)=10^{-%d}"          "$n" "$c"` ;; 
	fix)  t=`printf "%s&=%4.2f"                           "$n" "$c"` ;;
	*)    echo "Cut type $m unknown" >/dev/stderr ; t="?";; 
    esac
    
    echo "$t"
}


lastThree=
lastRun=
lastMethod=
processOne()
{
    local dir=`dirname $1` 
    local file=`basename $1` 
    local tgt=${dir}/${file}
    local base=`echo $tgt | sed -e 's/\.png//'`
    local run=`echo $dir | sed -e 's/_.*//'`
    local nstrip=`echo $dir | sed -e 's/.*\([23]\)strip.*/\1/'` 
    local dcCut=`getCutText "dc" "t_{\\text{hit}}" "$dir"`
    local slCut=`getCutText "sl" "m_{\\text{low}}" "$dir"`
    local shCut=`getCutText "sh" "m_{\\text{high}}" "$dir"`
    local strCut="N_{\\text{merged}}&\\leq${nstrip}"
    local title=`echo "${slCut} ${shCut} ${dcCut}" | tr -d '&'`
    title=`echo "$title" | tr -d '\\\\' | sed 's/_{text{\([^}]*\)}}/\1/g'`
    title=`echo "$title" | sed 's/Lw/xi/g'` 
    case ${run} in 
	137848) runText="Bad run" ;; 
	138190) runText="Good run" ;;
	last) return
    esac
    if test ! "X$lastRun" = "X$run"; then 
        echo "%% --- New Run ---"
	echo "\\section{${run}}" 
    fi
    lastRun=$run
    # if test ! "X$lastMethod" = "X$method"; then 
    # 	echo "\\subsubsection{${method}}" 
    # fi
    lastMethod=$method
    # cp $i doc/$tgt

    cat <<EOF > ${base}.tex
%% --- $1 ------------------------------
\\subsection{${title}}
\\begin{frame}{}
  \\begin{columns}
    \\begin{column}{.65\\linewidth}
      \\includegraphics[keepaspectratio,width=\\linewidth]{%
        $i}
    \\end{column}
    \\begin{column}{.33\\linewidth}
      \\footnotesize
      ${runText}

      \\begin{align*}
        ${strCut}\\\\
        ${slCut}\\\\
        ${shCut}\\\\
        ${dcCut}\\\\
      \\end{align*}
    \\end{column}
  \\end{columns}
\\end{frame}
EOF
    echo "\\input{${base}}"
}
# --- Process a set -------------------------------------------------
processSet()
{
  local glob="$1"
  # echo "glob=$glob"
  # ls ${glob}/dNdeta_CENT.png
  local l=`ls ${glob}/dNdeta_CENT.png 2> /dev/null`
  local n=`echo $l | wc -w` 
  local j=1
  for i in $l ; do 
      local d=`dirname $i` 
      printf "\r%3d/%3d %-60s ..." "$j" "$n" "$d" > /dev/stderr 
      processOne $i 
      printf " done" > /dev/stderr 
      let j=$j+1
  done
  echo "" > /dev/stderr 
}

# --- Make it --------------------------------------------------------
runLatex()
{
  local file=`basename $1` 
  local dir=`dirname $1`
  echo "Running PDFLaTeX on $1" > /dev/stderr 
  (cd $dir && \
      TEXINPUTS=".:../:" \
      pdflatex -interaction scrollmode -halt-on-error $file)
}

# --- Friendly message -----------------------------------------------
usage()
{
    cat <<EOF
Usage: $0 [OPTIONS] 

	-h,--help		This help 
	-p,--pattern GLOB	Input glob pattern
	-o,--out     FILE       Output LaTeX file
EOF
}

out="doc/collecteddNdeta.tex"
pattern="*dndeta_*"
while test $# -gt 0 ; do
    case $1 in 
	-h|--help) usage; exit 0 ;; 
	-o|--out)  out=$2 ; shift ;; 
	-p|--pattern) pattern="$2" ; shift ;; 
    esac
    shift 
done



mkdir -p `basename $out`
startDoc                    >  ${out}
processSet "${pattern}"     >> ${out}
endDoc                      >> ${out}

runLatex $out > /dev/null
# runLatex $out > /dev/null

#
# EOF
#
