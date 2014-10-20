############## Make tests ###########################
include ../test.inc

NAME_BENCHMARK=mc-tester.root-1000M-NLO-EXP
# Or use other benchmark reference files:
#NAME_BENCHMARK=mc-tester.root-1000M-LO-EXP
#NAME_BENCHMARK=mc-tester.root-1000M-LO-1-phot
#NAME_BENCHMARK=mc-tester.root-1000M-NLO-1-phot

CONF=pythia_Z.conf
NAME=ScalarNLO

PROGRAM=$(PHOTOS_TAUOLA)
PYTHIA_MODE=$(Z_MU_MODE)
TAUOLA_MODE=$(K-STAR_MODE)

ONE_PHOTON_MODE= 0
USE_NLO        = 1

tester.$(NAME).pdf: mc-tester.$(NAME).root
	$(COMPARE) $(NAME_BENCHMARK) mc-tester.$(NAME).root
	mv tester.pdf tester.$(NAME).pdf
	$(VIEWER) tester.$(NAME).pdf &

mc-tester.$(NAME).root: $(PROGRAM) $(CONF)
	$(PROGRAM) $(CONF) $(PYTHIA_MODE) $(EVENTS) $(TAUOLA_MODE) $(ONE_PHOTON_MODE) $(USE_NLO)
	mv mc-tester.root mc-tester.$(NAME).root

clean:
	rm -rf *~ booklet* tester* mc-results* 
	rm -rf MC-TESTER.DAT *_C.so *_C.d

clobber: clean
	rm -f mc-tester.$(NAME).root 
