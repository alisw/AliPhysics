//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 1998      Caltech, UCSB
//
// Module: EvtSymTable.cc
//
// Description: Class to hold the symbols that are defined
//              in the DECAY files.
// Modification history:
//
//    RYD     May 8, 1997         Module created
//
//------------------------------------------------------------------------
//
#include "EvtGenBase/EvtPatches.hh"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <ctype.h>
#include "EvtGenBase/EvtSymTable.hh"
#include "EvtGenBase/EvtReport.hh"
#include <string>
using std::endl;
using std::fstream;

std::map<std::string,std::string> EvtSymTable::_symMap;


EvtSymTable::~EvtSymTable(){}

EvtSymTable::EvtSymTable() {

}

void EvtSymTable::define(const std::string& symname,std::string d) {

  if ( _symMap.find(symname)!=_symMap.end() ) {
    report(Severity::Info,"EvtGen") << "Symbol:"<<symname.c_str()<<
      " redefined, old value:"<<_symMap[symname].c_str()<<" new value:"<<d.c_str()<<endl;
    _symMap[symname]=d;
    return;
  }

  _symMap[symname]=d;
  return;
}

std::string EvtSymTable::get(const std::string& symname,int& ierr) {

  ierr=0;

  if ( _symMap.find(symname)!=_symMap.end() ) return _symMap[symname];

  // If no matching symbol found just return the string

  return symname;
}

