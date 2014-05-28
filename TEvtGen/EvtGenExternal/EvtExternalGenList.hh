//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package. If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 2012      University of Warwick, UK
//
// Module: EvtExternalGenFactory
//
// Description: A factory type method to create engines for external physics
// generators like Pythia.
//
// Modification history:
//
//    John Back       Sept 2012           Module created
//
//------------------------------------------------------------------------------
//

#ifndef EVTEXTERNALGENLIST_HH
#define EVTEXTERNALGENLIST_HH

#include "EvtGenBase/EvtAbsRadCorr.hh"
#include "EvtGenBase/EvtDecayBase.hh"

#include <list>

class EvtExternalGenList {

public:

  EvtExternalGenList(bool convertPythiaCodes = true, std::string pythiaXmlDir = "", 
		     std::string photonType = "gamma", bool useEvtGenRandom = true);

  virtual ~EvtExternalGenList();

  std::list<EvtDecayBase*> getListOfModels();

  EvtAbsRadCorr* getPhotosModel();

protected:

private:

};

#endif

