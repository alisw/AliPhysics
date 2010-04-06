//-------------------------------------------------------------------------
//
// A calss for keeping the MC track information used in  
// the comparison tasks by:  Andrei.Zalite@cern.ch
//
//-------------------------------------------------------------------------

#include "AliMCComparisonTrack.h"

AliMCComparisonTrack::AliMCComparisonTrack() :
  fMCLabel(0), fMCPdg(0), fPz(0), fPt(0), fPhi(0), 
  fLocalX(0), fLocalY(0), fZ(0)
  {}
