//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 2002      Caltech, LLNL
//
// Module: EvtGen/EvtModelAlias.hh
//
// Description:Class to keep track of model aliases 
//             read in from the decay table
//
// Modification history:
//
//    Lange     January 19, 2002         Module created
//
//------------------------------------------------------------------------

#ifndef EVTMODELALIAS_HH
#define EVTMODELALIAS_HH

#include <vector>
#include <string>

class EvtModelAlias{

public:

  EvtModelAlias() {}; 
  EvtModelAlias(std::string alias, std::string model, std::vector<std::string> args); 
  ~EvtModelAlias() {};
  EvtModelAlias(const EvtModelAlias &copyMe);
  EvtModelAlias operator=(const EvtModelAlias &copyMe);
  bool matchAlias(const std::string &cand) {if (cand==_aliasName) return true; 
                                                          return false;}
  std::string getName() { return _model;} 
  std::vector<std::string> getArgList();
private:

  std::string _aliasName;
  std::string _model;
  std::vector<std::string> _modelArgs;

};
#endif
