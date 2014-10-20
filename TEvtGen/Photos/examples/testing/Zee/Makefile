############## Make tests ###########################
include ../test.inc

NAME_BENCHMARK=mc-tester.root-Zee-10M-nosym-10MeV
# Or use other benchmarks reference files:
# NAME_BENCHMARK=mc-tester.root-Zee-kkmc-sym-1GeV
# NAME_BENCHMARK=mc-tester.root-Zee-photos-100M-sym-10MeV

CONF=pythia_Z.conf
NAME=Zee

PROGRAM=$(PHOTOS_ONLY)
PYTHIA_MODE=$(PP_MODE)

tester.$(NAME).pdf: mc-tester.$(NAME).root
	$(COMPARE) $(NAME_BENCHMARK) mc-tester.$(NAME).root
	mv tester.pdf tester.$(NAME).pdf
	$(VIEWER) tester.$(NAME).pdf &
	@echo
	@echo "Error messages redirected to /dev/null due to large output from Pythia"
	@echo

mc-tester.$(NAME).root: $(PROGRAM) $(CONF)
	$(PROGRAM) $(CONF) $(PYTHIA_MODE) $(EVENTS) 2>/dev/null
	mv mc-tester.root mc-tester.$(NAME).root

clean:
	rm -rf *~ booklet* tester* mc-results* 
	rm -rf MC-TESTER.DAT *_C.so *_C.d

clobber: clean
	rm -f mc-tester.$(NAME).root 
