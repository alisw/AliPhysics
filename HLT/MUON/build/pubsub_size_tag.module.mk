############################################################################
#
# Author: Gareth de Vaux
# Email:  devaux@lhc.phy.uct.ac.za | dhlt@lordcow.org
#
############################################################################

BINARY  = pubsub_size_tag

SOURCES = pubsub_size_tag.cxx

MACROS =

############################################################################

include build/config.mk
include $(BUILD_DIR)/rules.mk
