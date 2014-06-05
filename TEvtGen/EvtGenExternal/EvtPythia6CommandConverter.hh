//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package. If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 2012      University of Warwick, UK
//
// Module: EvtPythia6CommandConverter
//
// Description:  Function to replace Pythia 6 commands with the
//               corresponding Pythia 8 commands.
//
// Modification history:
//
//    Daniel Craik       March 2012            Module created
//
//------------------------------------------------------------------------

#include "EvtGenBase/EvtExtGeneratorCommandsTable.hh"

#include <string>
#include <vector>

#ifndef EVTPYTHIA6COMMANDCONVERTER_HH
#define EVTPYTHIA6COMMANDCONVERTER_HH

std::vector<std::string> convertPythia6Command(Command command);

#endif
