

FILENAME=Photos_interface_design


all:	$(FILENAME).dvi
#	dvips $(FILENAME).dvi
	dvipdf $(FILENAME).dvi
#	gv $(FILENAME).ps

$(FILENAME).dvi: $(FILENAME).tex $(FILENAME).cls $(FILENAME).bib
	latex $(FILENAME).tex
	bibtex $(FILENAME)
	latex $(FILENAME).tex
	latex $(FILENAME).tex


clean:
	rm -f $(FILENAME).dvi
	rm -f *~ *.aux *.log *.toc *.bak *.out *.blg *.bbl

clobber: clean
	rm -f $(FILENAME).ps
