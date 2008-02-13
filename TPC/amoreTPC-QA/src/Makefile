include Makefile.inc

INCLUDES = $(all_includes)
METASOURCES = AUTO

SUBDIRS = common publisher subscriber ui

default: all

all:
	@cd common; make;
	@cd publisher; make;
	@cd subscriber; make;
	@cd ui; make;

install: all
	@cd common; make install;
	@cd publisher; make install;
	@cd subscriber; make install;
	@cd ui; make install;

uninstall:
	@cd common; make uninstall;
	@cd publisher; make uninstall;
	@cd subscriber; make uninstall;
	@cd ui; make uninstall;

clean:
	@cd common; make clean;
	@cd publisher; make clean;
	@cd subscriber; make clean;
	@cd ui; make clean;
