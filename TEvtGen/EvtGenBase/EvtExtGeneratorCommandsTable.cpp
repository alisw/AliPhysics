//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package. If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 2012      University of Warwick, UK
//
// Module: EvtExtGeneratorCommandsTable
//
// Description:  Table of commands to pass to external generators
//
// Modification history:
//
//    Daniel Craik       March 2012            Module created
//
//------------------------------------------------------------------------

#include "EvtGenBase/EvtExtGeneratorCommandsTable.hh"

EvtExtGeneratorCommandsTable::EvtExtGeneratorCommandsTable() {
  _commandMap.clear();
}

EvtExtGeneratorCommandsTable::~EvtExtGeneratorCommandsTable() {
  _commandMap.clear();
}

EvtExtGeneratorCommandsTable* EvtExtGeneratorCommandsTable::getInstance() {

  static EvtExtGeneratorCommandsTable* theCommandMap = 0;

  if (theCommandMap == 0) {
    theCommandMap = new EvtExtGeneratorCommandsTable();
  }

  return theCommandMap;

}
