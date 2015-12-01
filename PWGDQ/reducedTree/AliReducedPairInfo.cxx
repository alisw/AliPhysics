/*
***********************************************************
  Implementation of AliReducedPairInfo class.
  Contact: iarsene@cern.ch
  2015/04/08
  *********************************************************
*/

#ifndef ALIREDUCEDPAIRINFO_H
#include "AliReducedPairInfo.h"
#endif

ClassImp(AliReducedPairInfo)


//_______________________________________________________________________________
AliReducedPairInfo::AliReducedPairInfo() :
  AliReducedBaseTrack(),
  fCandidateId(-1),
  fPairType(-1), 
  fLegIds(),
  fMass(),
  fLxy(0.0),
  fPointingAngle(0.0),
  fChisquare(0.0)
{
  //
  // Constructor
  //
  fLegIds[0] = 0; fLegIds[1] = 0;
  fMass[0]=-999.; fMass[1]=-999.; fMass[2]=-999.; fMass[3]=-999.;
}


//_______________________________________________________________________________
AliReducedPairInfo::AliReducedPairInfo(const AliReducedPairInfo &c) :
  AliReducedBaseTrack(c),
  fCandidateId(c.CandidateId()),
  fPairType(c.PairType()),
  fLegIds(),
  fMass(),
  fLxy(c.Lxy()),
  fPointingAngle(c.PointingAngle()),
  fChisquare(c.Chi2())
{
  //
  // copy constructor
  //
  fLegIds[0] = c.LegId(0);
  fLegIds[1] = c.LegId(1);
  fMass[0] = c.Mass(0); fMass[1] = c.Mass(1); fMass[2] = c.Mass(2); fMass[3] = c.Mass(3);
}


//_______________________________________________________________________________
AliReducedPairInfo::~AliReducedPairInfo() {
  //
  // destructor
  //
}
