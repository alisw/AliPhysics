//--------------------------------------------------------------------------
//
//
// Copyright Information: See EvtGen/COPYRIGHT
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
//
// Module: EvtGen/EvtBtoXsgammaAbsModel.hh
//
// Description:
//      B->Xs gamma model base class.
//
// Modification history:
//
//    Jane Tinslay            March 21, 2000      Module Created
//                        
//------------------------------------------------------------------------

#ifndef EVTBTOXSGAMMAABSMODEL_HH
#define EVTBTOXSGAMMAABSMODEL_HH

class EvtBtoXsgammaAbsModel {

public:
  
  EvtBtoXsgammaAbsModel() {}

  virtual ~EvtBtoXsgammaAbsModel();

  virtual void init(int, double*);

  virtual double GetMass(int code)=0;

private:

};

#endif

