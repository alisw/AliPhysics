#-*- Mode: Makefile -*-


reveng-@MODULE@:		@MODULE@/check/classDiagram.dot

@MODULE@/check/classDiagram.dot:	$(PACKREVENG)
	@$(REV_ENG) $^
	@-mv classDiagram.dot $@

revdisp-@MODULE@:	reveng-@MODULE@
	@echo revdisp for @MODULE@
	@cd @MODULE@/check ; \
      $(IRST_INSTALLDIR)/webreveng/create-class-diagram-pages.sh
	@sed -e "s/STEER/@MODULE@/g" < $(IRST_INSTALLDIR)/webreveng/WWW/STEER/HomePage.html > @MODULE@/check/HomePage.html

PACKREVENG =

#  Clean option for the whole module
clean-@MODULE@:
ifndef ALIQUIET
		@echo "***** Cleaning @MODULE@ *****"
endif
		$(MUTE)rm @MODULE@/module.mk
		$(MUTE)rm -rf @MODULE@/tgt_$(ALICE_TARGET) 
		$(MUTE)rm -f $(@MODULE@LIBS)
		$(MUTE)rm -f $(@MODULE@BINS)

clean-check-@MODULE@:
ifndef ALIQUIET
		@echo "***** Cleaning code check for @MODULE@ *****"
endif
		$(MUTE)rm -f `find @MODULE@/check -name '*.i'` `find @MODULE@/check -name '*.ii'` `find @MODULE@/check -name '*.viol'`

clean-smell-@MODULE@:
ifndef ALIQUIET
		@echo "***** Cleaning code smell for @MODULE@ *****"
endif
		$(MUTE)rm -f `find @MODULE@/check -name '*.ml'` `find @MODULE@/check -name '*.smell'`

clean-reveng-@MODULE@:
ifndef ALIQUIET
		@echo "***** Cleaning reverse engineering files for @MODULE@ *****"
endif
		$(MUTE)rm -f `find @MODULE@/check -name '*.dot'`

