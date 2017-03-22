FILE		= primer
EPSS		= picts/Parab.eps			\
		  picts/ReconstructionFramework.eps	\
		  picts/SimulationFramework.eps	        \
		  picts/alienfig7.eps			\
		  picts/aligen.eps			\
		  picts/cocktail.eps			\
		  picts/evglib.eps			\
		  picts/folderstruct.eps		\
		  picts/kaon.eps			\
		  picts/loaderdiagram.eps		\
		  picts/partcuts.eps			\
		  picts/phipp.eps			\
		  picts/vmc.eps			        \
		  picts/tagana.eps		        \
		  picts/tagsources.eps	                \
		  picts/tagpost.eps
PDFS		= $(EPSS:%.eps=%.pdf)
LATEX		= latex
PDFLATEX	= pdflatex
DVIPS		= dvips 
PSTOPDF		= epstopdf 

%.aux:	%.tex
	$(LATEX) $< 

%.dvi:	%.aux
	$(LATEX) $* 

%.ps:	%.dvi 
	$(DVIPS) $< -o $@ 

%.pdf:	%.tex
	$(PDFLATEX) $*
	$(PDFLATEX) $*

%.pdf:	%.eps
	$(PSTOPDF) $< 

default: $(FILE).pdf

clean:
	rm -f 	$(FILE).dvi 	\
		$(FILE).aux 	\
		$(FILE).toc 	\
		$(FILE).log 	\
		$(FILE).ps 	\
		$(FILE).out	\
		$(FILE).lol	\
		$(FILE).pdf $(PDFS)


dist:
	mkdir -p aliroot-primer 
	cp -a $(FILE).tex aliroot-primer/
	cp -a $(EPSS)     aliroot-primer/
	cp -a Makefile    aliroot-primer/
	cp -a scripts     aliroot-primer/
	tar -cvf - aliroot-primer | gzip - > aliroot-primer.tar.gz 
	rm -rf aliroot-primer 

$(FILE).dvi:$(FILE).aux	$(EPSS)
$(FILE).pdf:$(FILE).tex	$(PDFS)
