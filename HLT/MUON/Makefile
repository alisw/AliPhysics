#/**************************************************************************
# * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
# *                                                                        *
# * Author: Artur Szostak                                                  *
# * Email:  artur@alice.phy.uct.ac.za                                      *
# *                                                                        *
# * Permission to use, copy, modify and distribute this software and its   *
# * documentation strictly for non-commercial purposes is hereby granted   *
# * without fee, provided that the above copyright notice appears in all   *
# * copies and that both the copyright notice and this permission notice   *
# * appear in the supporting documentation. The authors make no claims     *
# * about the suitability of this software for any purpose. It is          *
# * provided "as is" without express or implied warranty.                  *
# **************************************************************************/

include build/config.mk

.PHONY : all debug verbose depend doc clean clean-all incmajor incminor incbuild dist debug-dist

MODULES = $(wildcard $(BUILD_DIR)/*.module.mk)

############################################################################

all : $(MODULES)
	$(VERBOSE) for MODULE in $(MODULES) ; do $(MAKE) -f $$MODULE ; done

############################################################################

debug :
	@echo "---------- Entering debug build mode ----------"
	@$(MAKE) DEBUGGING=ON all

############################################################################

verbose :
	@echo "---------- Entering verbose build mode ----------"
	@$(MAKE) VERBOSE_MAKE=ON all

############################################################################

depend : $(MODULES)
	$(VERBOSE) for MODULE in $(MODULES) ; do $(MAKE) -f $$MODULE depend ; done

############################################################################

doc :
ifndef VERBOSE_MAKE
	@echo "making source code documentation"
endif
	$(VERBOSE) mkdir -p srcdoc/html/
	$(VERBOSE) cp doc/doxygen/footer.html doc/doxygen/header.html doc/doxygen/style.css srcdoc/html/
	$(VERBOSE) doxygen doc/doxygen/dHLTconfig

############################################################################

clean :
ifndef VERBOSE_MAKE
	@echo "cleaning up $(OUTPUT_DIR)"
endif
	$(VERBOSE) rm -f -r $(OUTPUT_DIR)

############################################################################

clean-all : clean
ifndef VERBOSE_MAKE
	@echo "cleaning up $(BIN_DIR)"
endif
	$(VERBOSE) rm -f -r $(BIN_DIR)
ifndef VERBOSE_MAKE
	@echo "cleaning up $(LIB_DIR)"
endif
	$(VERBOSE) rm -f -r $(LIB_DIR)
ifndef VERBOSE_MAKE
	@echo "cleaning up $(INCLUDE_DIR)"
endif
	$(VERBOSE) rm -f -r $(INCLUDE_DIR)
	$(VERBOSE) rm -f $(TOP_DIR)/core $(BIN_DIR)/core $(MACRO_DIR)/core

############################################################################

incmajor :
ifndef VERBOSE_MAKE
	@echo Incrementing major version number.
endif
	$(VERBOSE) let NUM=`cat $(SRC_DIR)/Version/major.number`+1 && echo $$NUM > $(SRC_DIR)/Version/major.number

############################################################################

incminor :
ifndef VERBOSE_MAKE
	@echo Incrementing minor version number.
endif
	$(VERBOSE) let NUM=`cat $(SRC_DIR)/Version/minor.number`+1 && echo $$NUM > $(SRC_DIR)/Version/minor.number

############################################################################

incbuild : $(MODULES)
	$(VERBOSE) for MODULE in $(MODULES) ; do \
		$(MAKE) -q -f $$MODULE ; \
		if test $$? -eq 1 ; then \
			$(MAKE) -f $$MODULE incbuild ; \
		fi; \
	done

############################################################################

dist : incbuild
	@echo "Making distribution $(DISTRIBUTION_FILE)"
	$(VERBOSE) mkdir -p $(DISTRIBUTION_SUBDIR)
	$(VERBOSE) for MODULE in $(MODULES) ; do $(MAKE) -f $$MODULE create-distrib ; done
	$(VERBOSE) cp -r $(INCLUDE_DIR) $(DISTRIBUTION_SUBDIR)/
	$(VERBOSE) cd $(DISTRIBUTION_DIR) && tar -czf $(TOP_DIR)/$(DISTRIBUTION_FILE) $(DISTRIBUTION_DIRNAME)
	
############################################################################

debug-dist : 
	@$(MAKE) DEBUGGING=ON dist

