The FLOW.tex (FLOW.pdf) is the main documentation of the ALICE package.
otherdocs contains publications which explain the various flow methods used and
misc documentation and presentations 

To re-create the flow package manual from the FlowPackageManual.tex source file, 
including the html files and index, put the following in a shell-script and run under linux / unix:

file=FlowPackageManual
echo "Processing $file"
pdflatex --interaction=batchmode $file
bibtex $file
makeindex FlowPackageManual.idx
pdflatex --interaction=batchmode $file
pdflatex --interaction=batchmode $file
#parse as self-container html
pandoc -r latex -w html -S -s -m -N --toc --highlight-style tango --indented-code-classes numberLines --self-contained -o FlowPackageManual.html FlowPackageManual.tex
#cleanup for git
rm *.aux
rm *.bbl
rm *.blg
rm *.idx
rm *.ilg
rm *.ind
rm *.log
rm *.toc
rm *.out
