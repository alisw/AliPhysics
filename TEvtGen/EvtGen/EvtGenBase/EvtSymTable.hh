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
// Module: EvtGen/EvtSymTable.hh
//
// Description: Class to hold the symbols that are defined
//              in the DECAY files.
//
// Modification history:
//
//    RYD     May 8, 1998         Module created
//
//------------------------------------------------------------------------

#ifndef EVTSYMTABLE_HH
#define EVTSYMTABLE_HH

#include <string>
#include <map>

class EvtSymTable {

public:

  EvtSymTable(); 
  ~EvtSymTable();

  static void define(const std::string& name,std::string d);

  static std::string get(const std::string& name,int &ierr);

private:

  static std::map<std::string,std::string> _symMap;
  
}; 

#endif






