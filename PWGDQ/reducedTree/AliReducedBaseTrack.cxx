/*
***********************************************************
  Implementation of AliReducedBaseTrack class.
  Contact: iarsene@cern.ch
  2015/04/08
  *********************************************************
*/

#ifndef ALIREDUCEDBASETRACK_H
#include "AliReducedBaseTrack.h"
#endif

ClassImp(AliReducedBaseTrack)

//_______________________________________________________________________________
AliReducedBaseTrack::AliReducedBaseTrack() :
  fTrackId(0),
  fIsCartesian(kFALSE),
  fCharge(0),
  fFlags(0),
  fQualityFlags(0),
  fMCFlags(0),
  fIsMCTruth(kFALSE)
{
  //
  // Constructor
  //
  fP[0]=0.0; fP[1]=0.0; fP[2]=0.0;
}

//_______________________________________________________________________________
AliReducedBaseTrack::AliReducedBaseTrack(const AliReducedBaseTrack &c) :
  TObject(c),
  fTrackId(c.fTrackId),
  fIsCartesian(c.IsCartesian()),
  fCharge(c.Charge()),
  fFlags(c.GetFlags()),
  fQualityFlags(c.GetQualityFlags()),
  fMCFlags(c.GetMCFlags()),
  fIsMCTruth(c.IsMCTruth())
{
  //
  // Copy constructor
  //
  if(c.IsCartesian()) {fP[0]=c.Px();fP[1]=c.Py();fP[2]=c.Pz();}
  else {fP[0]=c.Pt();fP[1]=c.Phi();fP[2]=c.Eta();}
}

//_______________________________________________________________________________
AliReducedBaseTrack::~AliReducedBaseTrack()
{
  //
  // De-Constructor
  //
}

//_______________________________________________________________________________
UInt_t AliReducedBaseTrack::GetFirstHalfOfQualityFlags() const {
   //
   // Get the first half of the fQualityFlags (bits 0->31)
   //
   UInt_t map = 0;
   for(UShort_t i=0; i<32; ++i) 
      if(TestQualityFlag(i)) 
         map |= (UInt_t(1)<<i);
   return map;   
}

//_______________________________________________________________________________
UInt_t AliReducedBaseTrack::GetSecondHalfOfQualityFlags() const {
   //
   // Get the second half of the fQualityFlags (bits 32->63)
   //
   UInt_t map = 0;
   for(UShort_t i=0; i<32; ++i) 
      if(TestQualityFlag(i+32)) 
         map |= (UInt_t(1)<<i);
      return map;   
}
