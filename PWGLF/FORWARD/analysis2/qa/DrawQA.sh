#!/bin/sh 
#
# @file   DrawQA.sh
# @author Christian Holm Christensen <cholm@nbi.dk>
# @date   Thu Nov 17 11:31:03 2011
# 
# @brief  Draw most QA stuff
# 
# @deprecated Use QATrender instead
# @ingroup pwg2_forward_qa_scripts
file=trending.root 

# --- Help message ---------------------------------------------------
usage ()
{
    cat <<EOF
Usage: $0 [OPTIONS] 

Options:
	-f,--file    FILE	Input file 
	-t,--title   STRING	Title string 
	-o,--output  FILE	Write output on this file 
	-h,--help		Show this help
EOF
}

# --- Command line parsing -------------------------------------------
while test $# -gt 0 ; do 
    case $1 in  
	-f|--file) file=$2 ; shift ;;
	-t|--title) title=$2 ; shift ;; 
	-h|--help) usage ; exit 0 ;; 
	*)  echo "Uknown option $1" > /dev/stderr ; exit 1 ;; 
    esac
    shift 
done 

# --- Run aliroot ----------------------------------------------------
scr=$ALICE_ROOT/PWG2/FORWARD/analysis2/qa/DrawQA.C 

root -l -b -q ${scr}\(\"$file\"\)

# --- Make LaTeX code ------------------------------------------------
if test "x$title" = "x" ; then 
    title=`echo "QA from $file" | sed -e 's/_/\\_/g'`
fi
if test "x$output" = "x" ; then 
    output=`basename $file .root` 
else 
    output=`basename $file .pdf` 
fi
doc=${output}.tex

echo $title
cat <<EOF > $doc
\documentclass[landscape,12pt,a4paper]{article}
\usepackage[a4paper,margin=2cm]{geometry}
\usepackage{graphicx}
\title{$title}
\author{FMD Team}
\date{\today}
\begin{document}
\maketitle
\clearpage

EOF
pngs="fit_results neighbors beforeAfter 123 recAnaELoss occupancy elossVsPoisson"
for i in $pngs ; do
    case $i in  
	fit_results) t="Energy loss fits" ;; 
	neighbors)   t="Correlation of neighbors" ;; 
	beforeAfter) t="Effect of sharing correction" ;; 
	123)         t="Energy loss from single, double, and tripple hits" ;; 
	recAnaELoss) t="Energy loss from reconstruction and used in analysis" ;;
	occupancy)   t="Calculated occupancy" ;;
	elossVsPoisson) t="Correlation of Poisson and Energy loss methods" ;;
	*)           t="Unkknown" ;;
    esac
    cat <<EOF >> $doc
\section*{$t}
\begin{center}
\includegraphics[keepaspectratio,height=.9\textheight]{$i}
\end{center}
\clearpage
EOF
done

cat <<EOF >> $doc
\end{document}
EOF

pdflatex $doc    

rm -f $doc $output.aux $output.log FitResults.pdf 
for i in $pngs ; do rm -f $i.png ; done 

# 
# EOF
# 

