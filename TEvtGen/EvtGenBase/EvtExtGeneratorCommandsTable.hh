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

#include <map>
#include <vector>
#include <string>

#ifndef EVTEXTGENERATORCOMMANDSTABLE_HH
#define EVTEXTGENERATORCOMMANDSTABLE_HH

typedef std::map<std::string, std::string> Command;
typedef std::vector<Command> GeneratorCommands;
typedef std::map<std::string, GeneratorCommands > GlobalCommandMap;

class EvtExtGeneratorCommandsTable {

public:

  static EvtExtGeneratorCommandsTable* getInstance();

  void addCommand(std::string extGenerator, Command command) { _commandMap[extGenerator].push_back(command); }
  const GeneratorCommands& getCommands(std::string extGenerator) { return _commandMap[extGenerator]; }

protected:

  EvtExtGeneratorCommandsTable();
  ~EvtExtGeneratorCommandsTable();

private:

  GlobalCommandMap _commandMap;

  EvtExtGeneratorCommandsTable(const EvtExtGeneratorCommandsTable&) {};

};

#endif
