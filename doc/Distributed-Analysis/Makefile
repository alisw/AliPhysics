all:	DistributedAnalysis.pdf

DistributedAnalysis.pdf:	DistributedAnalysis.dvi DistributedAnalysis.tex
	dvipdf DistributedAnalysis.dvi DistributedAnalysis.pdf
	rm *.dvi

DistributedAnalysis.dvi:		DistributedAnalysis.tex
	latex DistributedAnalysis.tex

clean:
	rm -rf *~ *.pdf *.dvi *.log *.aux *.tmp

