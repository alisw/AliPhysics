# $Id$
# Makefile template for the Alice HLT Alice HLT Online Monitoring GUI


REAL_SRC		= AliHLTGUI.cxx \
			  AliHLTGUI.h \
			  main.cxx \
			  AliHLTGUIMainForm.ui \
			  AliceLogo.gif \
			  KIPLogo.gif

REAL_SRC_DIR		= $(top_srcdir)/src


QT_BUILD_DIR		= qtbuild
BUILT_SOURCES		= qtbuild/flags

bin_SCRIPTS		= src/startAliHLTGUI \
			  qtbuild/AliHLTGUI

SUBDIRS 		= . doc
DIST_SUBDIRS		= $(SUBDIRS)

EXTRA_DIST		= $(REAL_SRC:%=src/%) \
			  src/AliHLTGUI.pro \
			  src/flags.in \
			  src/startAliHLTGUI.in

CLEANFILES		= 

LFLAGS                  = $(LDFLAGS)

$(QT_BUILD_DIR):
	test -d $(QT_BUILD_DIR) || mkdir $(QT_BUILD_DIR)
#	$(foreach i, $(REAL_SRC:%=$(REAL_SRC_DIR)/%), (cd $(QT_BUILD_DIR) && ln -sf ../$i );)
#	ln -sf ../src/flags $(QT_BUILD_DIR)/flags
#	ln -sf ../src/AliHLTGUI.pro $(QT_BUILD_DIR)/AliHLTGUI.pro

$(REAL_SRC:%=$(QT_BUILD_DIR)/%): $(REAL_SRC:%=$(REAL_SRC_DIR)/%)
	test -h $@ || ln -s $(patsubst $(QT_BUILD_DIR)/%, ../$(REAL_SRC_DIR)/%, $@) $@

$(QT_BUILD_DIR)/flags: src/flags
	test -d $(QT_BUILD_DIR) || mkdir $(QT_BUILD_DIR)
	test -h $@ || ln -s ../$< $@

$(QT_BUILD_DIR)/AliHLTGUI.pro: $(top_srcdir)/src/AliHLTGUI.pro
	test -d $(QT_BUILD_DIR) || mkdir $(QT_BUILD_DIR)
	test -h $@ || ln -s ../$< $@

$(QT_BUILD_DIR)/Makefile: $(QT_BUILD_DIR)/AliHLTGUI.pro
	(cd $(QT_BUILD_DIR) \
	&& export PATH=@QT_DIR@/bin:$(PATH) \
	&& export QTDIR=@QT_DIR@ \
	&& qmake  AliHLTGUI.pro)
	cat $@ | sed -e 's|^LFLAGS\([ \t]*\)=|LFLAGS\1= $$(LDFLAGS)|' > $@.tmp
	mv $@.tmp $@
	cat $@ | sed -e 's|^CXXFLAGS\([ \t]*\)=|CXXFLAGS\1+=|' > $@.tmp
	mv $@.tmp $@

$(QT_BUILD_DIR)/AliHLTGUI: $(QT_BUILD_DIR) $(REAL_SRC:%=$(QT_BUILD_DIR)/%) $(QT_BUILD_DIR)/Makefile
	(cd $(QT_BUILD_DIR)  \
	&& export PATH=@QT_DIR@/bin:$(PATH) \
	&& export QTDIR=@QT_DIR@ \
	&& $(MAKE) CXXFLAGS="$(CXXFLAGS) -DALIHLTTPCDISPLAY_VERSION=@ALIHLTTPCDISPLAY_VERSION@ -DPACKAGE_VERSION=\"@PACKAGE_VERSION@\"")

clean-local:
	rm -rf $(QT_BUILD_DIR)

#
# EOF
#
