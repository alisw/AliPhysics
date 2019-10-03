#ifndef ALIETACLUSTERMCSTATSTRUCT_H
#define ALIETACLUSTERMCSTATSTRUCT_H

#include <TObject.h>

/// cluster information structure
class ClusterMCStatStruct: public TObject
{

  public:

  /// object name (re-implemented)
  virtual const char*	GetName() const
  { return "ClusterMCStatStruct"; }

  /// default contructor
  ClusterMCStatStruct( void ):
  daughter1(0),
  daughter2(0),
  EnergyMC(0),
  pTMC(0),
  pidMC(0),
  parentID(0),
  etaMC(0),
  phiMC(0),
  thetaMC(0),
  vXmc(0), vYmc(0), vZmc(0),
  EventsMC(0)
{}

 //characteristics pion
Int_t    daughter1;
Int_t    daughter2;
 Float_t  EnergyMC;
 Float_t  pTMC;
 Int_t    pidMC;
 Int_t    parentID;
 Float_t  etaMC;
 Float_t  phiMC;
 Float_t  thetaMC;
 Float_t  vXmc;
 Float_t  vYmc;
 Float_t  vZmc;
 Int_t    EventsMC;
  //
  ClassDef(ClusterMCStatStruct,1)

};

#endif
