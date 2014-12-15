//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package. If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//
// Module: EvtGenModels/EvtBcBsNPi.hh
//
// Description: Decay model for Bc -> Bs + npi
//
// Modification history:
//
//    A.Berezhnoy, A.Likhoded, A.Luchinsky  July 2011   Module created
//
//------------------------------------------------------------------------

#ifndef EvtBcBsNpi_HH
#define EvtBcBsNpi_HH

#include "EvtGenModels/EvtBcToNPi.hh"
#include "EvtGenBase/EvtDecayBase.hh"
#include <string>

class EvtBcBsNPi : public EvtBcToNPi {
public:

  EvtBcBsNPi();
  virtual ~EvtBcBsNPi();

  virtual void init();	
  virtual void initProbMax();
  
  virtual std::string getName();
  virtual EvtDecayBase* clone();

};

#endif
