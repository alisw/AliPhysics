#ifndef ALITOFT0V2_H
#define ALITOFT0V2_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//----------------------------------------------------------------------------//
//                                                                            //
//   Description: class to performe an event time measurment with TOF.        //
//                allowing for multiple-event-time --> PileUp                 //
//                                                                            //
//----------------------------------------------------------------------------//

#include "AliTOFT0v1.h"

class AliTOFT0v2: public AliTOFT0v1 {
 public:
  AliTOFT0v2(AliESDpid *extPID=NULL); // default constructor
  AliTOFT0v2(AliESDEvent *event,AliESDpid *extPID=NULL); // overloaded constructor
  virtual ~AliTOFT0v2() ; // dtor


  virtual void Reset();
  virtual Double_t* DefineT0(Option_t *option,Float_t pMinCut=3,Float_t pMaxCut=5,Int_t isCalibrationMode=0); 

 private:
  static const Int_t fgkMaxTracks = 2000;
  Int_t fStartTimeIndex[fgkMaxTracks]; //! index of the associated start time
  Int_t fTrackLabel[fgkMaxTracks];     //! track label map TOF -> ESD
  Int_t fCurrentStartTime;             //! number of calls to DefineT0

  ClassDef(AliTOFT0v2,1);  // Calculate the time zero using TOF detector */

};
#endif 
