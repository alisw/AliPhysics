
#  Clean option for the whole module
clean-@MODULE@:
ifndef ALIQUIET
		@echo "***** Cleaning @MODULE@ *****"
endif
		$(MUTE)rm @MODULE@/module.mk
		$(MUTE)rm -rf @MODULE@/tgt_$(ALICE_TARGET) 
		$(MUTE)rm -f $(@MODULE@LIBS)
		$(MUTE)rm -f $(@MODULE@BINS)
