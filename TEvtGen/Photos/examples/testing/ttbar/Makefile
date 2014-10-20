############## Make tests ###########################
include ../test.inc

NAME_BENCHMARK=mc-tester.root-ttbar-10M


CONF=pythia_ttbar.conf
NAME=ttbar

PROGRAM=$(PHOTOS_ONLY)
PYTHIA_MODE =$(PP_MODE)
SPECIAL_MODE=$(TOP_MODE)

tester.$(NAME).pdf: mc-tester.$(NAME).root
	$(COMPARE) $(NAME_BENCHMARK) mc-tester.$(NAME).root
	mv tester.pdf tester.$(NAME).pdf
	$(VIEWER) tester.$(NAME).pdf &
	@echo
	@echo "IMPORTANT:"
	@echo "Following example is missing hadronization of b b_bar pair."
	@echo
	@echo "Error messages redirected to /dev/null due to large output from Pythia"
	@echo

mc-tester.$(NAME).root: $(PROGRAM) $(CONF)
	$(PROGRAM) $(CONF) $(PYTHIA_MODE) $(EVENTS) $(SPECIAL_MODE) 2>/dev/null
	mv mc-tester.root mc-tester.$(NAME).root

clean:
	rm -rf *~ booklet* tester* mc-results* 
	rm -rf MC-TESTER.DAT *_C.so *_C.d

clobber: clean
	rm -f mc-tester.$(NAME).root
