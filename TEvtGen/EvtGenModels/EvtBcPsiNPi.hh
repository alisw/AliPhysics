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
// Description: Decay model for Bc -> J/psi + npi
//
// Modification history:
//
//    A.Berezhnoy, A.Likhoded, A.Luchinsky  April 2011   Module created
//
//------------------------------------------------------------------------

#ifndef EvtBcPsiNpi_HH
#define EvtBcPsiNpi_HH

#include "EvtGenModels/EvtBcToNPi.hh"
#include "EvtGenBase/EvtDecayBase.hh"

#include <string>

class EvtBcPsiNPi : public EvtBcToNPi {

public:

  EvtBcPsiNPi();

  virtual ~EvtBcPsiNPi();

  virtual void init();
  virtual void initProbMax();
	
  virtual std::string getName();
  virtual EvtDecayBase* clone();

};

#endif
