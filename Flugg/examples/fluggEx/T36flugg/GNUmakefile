# $Id$
# Flugg tag $Name$

# Id: GNUmakefile,v 2.1 1998/07/02 15:26:52 gunter Exp 
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := mainT36
G4TARGET := $(name)
G4EXLIB := true

ifndef FLUGGINSTALL
  FLUGGINSTALL = ../../..
endif

#the target of the first rule of the first makefile is the default goal
#so "make" will do lib and bin targets.
.PHONY: all
all: lib bin

include $(FLUGGINSTALL)/config/binmake.gmk
