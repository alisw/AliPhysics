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
// Module: EvtGen/EvtModelReg.hh
//
// Description: 
//
// Modification history:
//
//    Lange     April 18, 2002         Module created
//
//------------------------------------------------------------------------

#ifndef EVTMODELREG_HH
#define EVTMODELREG_HH

#include <list>

class EvtDecayBase ;

class EvtModelReg{

public:

  EvtModelReg(const std::list<EvtDecayBase*>* extraModels=0);

private:


};


#endif



