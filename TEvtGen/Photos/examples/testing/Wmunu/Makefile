############## Make tests ###########################
include ../test.inc

NAME_BENCHMARK=mc-tester.root-Wmunu-photos-100M-sym-10MeV
# Or use other benchmark reference files:

#WARNING! Comparison with this benchmark file requires:
#         - uncommenting few lines at the end of SETUP.C file
#         - usage of pythia_W_resonance_peak.conf -
#           comment out line 'CONF' below.
#NAME_BENCHMARK=mc-tester.root-Wmunu-photos-F77-100M-1-phot-5MeV

CONF=pythia_W.conf
#CONF=pythia_W_resonance_peak.conf
NAME=Wmunu

PROGRAM=$(PHOTOS_ONLY)
PYTHIA_MODE=$(DEFAULT_MODE)
SPECIAL_MODE=0
ONE_PHOTON_MODE=0

tester.$(NAME).pdf: mc-tester.$(NAME).root
	$(COMPARE) $(NAME_BENCHMARK) mc-tester.$(NAME).root
	mv tester.pdf tester.$(NAME).pdf
	$(VIEWER) tester.$(NAME).pdf &

mc-tester.$(NAME).root: $(PROGRAM) $(CONF)
	$(PROGRAM) $(CONF) $(PYTHIA_MODE) $(EVENTS) $(SPECIAL_MODE) $(ONE_PHOTON_MODE)
	mv mc-tester.root mc-tester.$(NAME).root

clean:
	rm -rf *~ booklet* tester* mc-results* 
	rm -rf MC-TESTER.DAT *_C.so *_C.d

clobber: clean
	rm -f mc-tester.$(NAME).root 
