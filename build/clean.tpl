
#  Clean option for the whole module
clean-@MODULE@:
		@echo "***** Cleaning @MODULE@ *****"
		rm @MODULE@/module.mk
		rm -rf @MODULE@/tgt_$(ALICE_TARGET) 
		rm -f $(@MODULE@LIBS)
		rm -f $(@MODULE@BINS)
