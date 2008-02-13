INCLUDES = $(all_includes)
INCLUDES += -I/usr/include/amore/

METASOURCES = AUTO

SUBDIRS = src

timestamp=$(shell date +%Y%m%d)
snapname=amoreTPC-QA$(timestamp)
snaptgz=$(snapname).tgz

default: all

all:
	@cd src; make;

install: all
	@cd src; make install;

uninstall:
	@cd src; make uninstall;

clean:
	@cd src; make clean;

dist:   clean
	@mkdir amoreTPC-QA
	@cp -r src AUTHORS ChangeLog README TODO VERSION Makefile amoreTPC-QA
	@tar cfz $(snaptgz) amoreTPC-QA;
	@rm -rf amoreTPC-QA

rpm:
