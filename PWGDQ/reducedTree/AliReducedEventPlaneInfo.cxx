/*
***********************************************************
  Implementation of AliReducedEventPlaneInfo class.
  Contact: iarsene@cern.ch
  2015/04/08
  *********************************************************
*/

#ifndef ALIREDUCEDEVENTPLANEINFO_H
#include "AliReducedEventPlaneInfo.h"
#endif

ClassImp(AliReducedEventPlaneInfo)


//_______________________________________________________________________________
AliReducedEventPlaneInfo::AliReducedEventPlaneInfo() :
 fQvector(),
 fEventPlaneStatus()
{
  //
  // default constructor
  //
  for(Int_t idet=0; idet<kNdetectors; ++idet) {
    for(Int_t ih=0; ih<fgkNMaxHarmonics; ++ih) {
      fEventPlaneStatus[idet][ih] = kUnset;
      for(Int_t ic=0; ic<2; ++ic)
	fQvector[idet][ih][ic] = 0.0;
    }
  }
}


//____________________________________________________________________________
AliReducedEventPlaneInfo::~AliReducedEventPlaneInfo()
{
  //
  // De-Constructor
  //
  ClearEvent();
}


//_____________________________________________________________________________
void AliReducedEventPlaneInfo::ClearEvent() {
  //
  // clear the event
  //
  for(Int_t idet=0; idet<kNdetectors; ++idet) {
    for(Int_t ih=0; ih<fgkNMaxHarmonics; ++ih) {
      fEventPlaneStatus[idet][ih] = kUnset;
      for(Int_t ic=0; ic<2; ++ic)
	fQvector[idet][ih][ic] = 0.0;
    }
  }
}


//____________________________________________________________________________
void AliReducedEventPlaneInfo::CopyEvent(const AliReducedEventPlaneInfo* event) {
  //
  // copy information from another event to this one
  //
  for(Int_t idet=0; idet<kNdetectors; ++idet) {
    for(Int_t ih=0; ih<fgkNMaxHarmonics; ++ih) {
      fQvector[idet][ih][0] = event->Qx(idet, ih+1);
      fQvector[idet][ih][1] = event->Qy(idet, ih+1);
      fEventPlaneStatus[idet][ih] = event->GetEventPlaneStatus(idet, ih+1);
    }
  }
}
