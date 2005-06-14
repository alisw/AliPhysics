############################################################################
#
# Author: Artur Szostak
# Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
#
############################################################################

# Implementation of make rules that are common to all modules.

.PHONY : all depend create-distrib symlink

all : $(TARGET_NAME) symlink $(addprefix $(INCLUDE_DIR)/,$(INCLUDE_HEADERS))

create-distrib : all
	@mkdir -p $(DISTRIBUTION_TARGET_DIR)
	$(VERBOSE) cp $(TARGET_NAME) $(DISTRIBUTION_TARGET_DIR)/

depend : $(DEPENDS)

symlink : $(SYMLINK_NAME)

############################################################################

$(TARGET_NAME) : $(OBJECTS)
#	Create the output directory if it is missing
	@mkdir -p $(@D)
ifdef VERBOSE_MAKE
	@echo -------------------------------------------------------------------------------
else
	@echo "linking $(^F) => $(@F)"
endif
	$(VERBOSE) $(CXX) $(LINKFLAGS) $^ -o $@
	@echo "Made $@"

############################################################################

ifdef DICTIONARY

# The following file must be generated for rootcint. Refer to the comment
# where the command is used.
$(OUTPUT_TARGET_DIR)/dummy.h :
	@mkdir -p $(@D)
	@touch $@

# Note the -p, +P, -P options makes rootcint use the real C++ preprocessor
# to preprocess the headers.
$(OUTPUT_TARGET_DIR)/$(DICTIONARY) : $(DICTIONARY_HEADERS) $(STUB_HEADERS) $(LINKDEF) $(OUTPUT_TARGET_DIR)/dummy.h
#	Create the output directory if it is missing
	@mkdir -p $(@D)
ifndef VERBOSE_MAKE
	@echo "generating $(@F)"
endif
#	$(VERBOSE) rootcint -f $@ -c $(INCLUDES) $(MACROS) -p +STUB $(STUB_HEADERS) -STUB $(DICTIONARY_HEADERS) $(LINKDEF)
# Since rootcint contains a bug and does not allow one to use -p with the
# (+/-)STUB option but does allow (+/-)P combinations we need to run the
# command as given bellow. Note the dummy file is there since rootcint
# requires at least one file before the +P or +STUB options.
	$(VERBOSE) rootcint -f $@ -c $(INCLUDES) $(MACROS) $(OUTPUT_TARGET_DIR)/dummy.h +P +STUB $(STUB_HEADERS) -STUB $(DICTIONARY_HEADERS) -P $(LINKDEF)


endif

############################################################################

$(OUTPUT_TARGET_DIR)/%.d : %.cxx
#	Create the output directory if it is missing
	@mkdir -p $(@D)
ifndef VERBOSE_MAKE
	@echo "making dependancies $(<F) => $(@F)"
endif
#	Note the sed command replaces the target xx.o with xx.o xx.d
	$(VERBOSE) $(CXX) $(CXXFLAGS) -MM $< | \
		sed "s/\(.*\)\.o/$(subst /,\/,$(@D)/)\1\.o $(subst /,\/,$(@D)/)\1\.d/g" > $@

############################################################################

$(OUTPUT_TARGET_DIR)/%.o : %.cxx
#	Create the output directory if it is missing
	@mkdir -p $(@D)
ifndef VERBOSE_MAKE
	@echo "compiling $(<F) => $(@F)"
endif
	$(VERBOSE) $(CXX) $(CXXFLAGS) -c $< -o $@

############################################################################

$(INCLUDE_DIR)/%.hpp : %.hpp
#	Create the output directory if it is missing
	@mkdir -p $(@D)
ifndef VERBOSE_MAKE
	@echo "copying header $(@F) to include directory."
endif
	$(VERBOSE) cp $< $@

############################################################################

$(SYMLINK_NAME) : $(TARGET_NAME)
#	Create the output directory if it is missing
	@mkdir -p $(@D)
ifndef VERBOSE_MAKE
	@echo "creating symbolic link to $<"
endif
	$(VERBOSE) rm -f $@
	$(VERBOSE) ln -s $< $@

############################################################################

incbuild : $(OUTPUT_TARGET_DIR)/build.number
ifndef VERBOSE_MAKE
	@echo Incrementing build version number for $(TARGET_NAME).
endif
	$(VERBOSE) let NUM=`cat $<`+1 && echo $$NUM > $<

$(OUTPUT_TARGET_DIR)/build.number :
	@mkdir -p $(@D)
	@echo 0 > $@

############################################################################

$(SRC_DIR)/Version/Version.cxx : $(OUTPUT_TARGET_DIR)/VersionNumbers.hpp

$(OUTPUT_TARGET_DIR)/VersionNumbers.hpp : major.number minor.number $(OUTPUT_TARGET_DIR)/build.number
ifndef VERBOSE_MAKE
	@echo Generating $(@F)
endif
	$(VERBOSE) echo // This file was generated with \"rules.mk\". Do not modify. > $@
	$(VERBOSE) echo \#define MAJOR_VERSION `cat $(word 1,$^)` >> $@
	$(VERBOSE) echo \#define MINOR_VERSION `cat $(word 2,$^)` >> $@
	$(VERBOSE) echo \#define BUILD_NUMBER  `cat $(word 3,$^)` >> $@

############################################################################

ifneq ($(findstring clean, $(MAKECMDGOALS)), clean)
-include $(DEPENDS)
endif

