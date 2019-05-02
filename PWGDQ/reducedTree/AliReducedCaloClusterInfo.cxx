/*
***********************************************************
  Implementation of AliReducedBaseTrack class.
  Contact: iarsene@cern.ch
  2015/04/08
  *********************************************************
*/

#ifndef ALIREDUCEDCALOCLUSTERINFO_H
#include "AliReducedCaloClusterInfo.h"
#endif

ClassImp(AliReducedCaloClusterInfo)

//_____________________________________________________________________________
AliReducedCaloClusterInfo::AliReducedCaloClusterInfo() :
 fFlags(0),
 fClusterID(-999),
 fType(kUndefined),
 fEnergy(-999.),
 fTrackDx(-999.),
 fTrackDz(-999.),
 fM20(-999.),
 fM02(-999.),
 fDispersion(-999.),
 fPosition(),
 fTOF(-999.),
 fNCells(0),
 fNMatchedTracks(0)
{
  //
  // default constructor
  //
  for(Int_t i=0;i<3;++i) fPosition[i] = -999.;
}


//_____________________________________________________________________________
AliReducedCaloClusterInfo::~AliReducedCaloClusterInfo()
{
  //
  // destructor
  //
}
