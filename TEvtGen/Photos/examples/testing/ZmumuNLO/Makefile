############## Make tests ###########################
include ../test.inc

NAME_BENCHMARK=mc-tester.root-Zmumu-NLO-100M
# Or use other benchmark reference files:
#NAME_BENCHMARK=mc-tester.root-Zmumu-LO-100M
#NAME_BENCHMARK=mc-tester.root-Zmumu-1-phot-NLO-100M
#NAME_BENCHMARK=mc-tester.root-Zmumu-1-phot-LO-100M
#NAME_BENCHMARK=mc-tester.root-Zmumu-CEEX1-100M-nosym-1GeV
#NAME_BENCHMARK=mc-tester.root-Zmumu-CEEX2-100M-nosym-1GeV
#NAME_BENCHMARK=mc-tester.root-Zmumu-kkmc-F77-10M-sym-1GeV
#NAME_BENCHMARK=mc-tester.root-Zmumu-photos-F77-100M-sym-10MeV

#WARNING! Comparison with this benchmark file requires:
#         - uncommenting few lines at the end of SETUP.C file
#         - usage of pythia_Z_resonance_peak.conf -
#           comment out line 'CONF' below.
#NAME_BENCHMARK=mc-tester.root-Zmumu-photos-F77-100M-1-phot-1GeV-NLO
#NAME_BENCHMARK=mc-tester.root-Zmumu-koralz-100M-1-phot-1GeV-NLO
#NAME_BENCHMARK=mc-tester.root-Zmumu-kkmc-F77-100M-nosym-1GeV-NLO

CONF=pythia_Z.conf
#CONF=pythia_Z_resonance_peak.conf
NAME=ZmumuNLO

PROGRAM=$(PHOTOS_ONLY)
PYTHIA_MODE=$(Z_MU_MODE)
SPECIAL_MODE=$(NLO_MODE)
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
