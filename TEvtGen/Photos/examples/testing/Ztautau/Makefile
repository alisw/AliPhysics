############## Make tests ###########################
include ../test.inc

NAME_BENCHMARK=mc-tester.root-pispectr-RCall
# Or use other benchmark reference files:
# NAME_BENCHMARK=mc-tester.root-pispectr-RCtau
# NAME_BENCHMARK=mc-tester.root-pispectr-RCno

CONF=pythia_Z.conf
NAME=Ztautau

PROGRAM=$(PHOTOS_TAUOLA)
PYTHIA_MODE=$(Z_MU_MODE)
TAUOLA_MODE=$(PIPI_MODE)

tester.$(NAME).pdf: mc-tester.$(NAME).root
	$(COMPARE) $(NAME_BENCHMARK) mc-tester.$(NAME).root
	mv tester.pdf tester.$(NAME).pdf
	$(VIEWER) tester.$(NAME).pdf &

mc-tester.$(NAME).root: $(PROGRAM) $(CONF)
	$(PROGRAM) $(CONF) $(PYTHIA_MODE) $(EVENTS) $(TAUOLA_MODE) 2>/dev/null
	mv mc-tester.root mc-tester.$(NAME).root

clean:
	rm -rf *~ booklet* tester* mc-results* 
	rm -rf MC-TESTER.DAT *_C.so *_C.d

clobber: clean
	rm -f mc-tester.$(NAME).root 
