#ifndef ALIT0CALIBRATOR_H
#define ALIT0CALIBRATOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
/* $Id$ */
////////////////////////////////////////////////
//  class for T0 calibration during reconstruction                //
////////////////////////////////////////////////


#include "TNamed.h"
//#include "AliT0RecoParam.h"

class AliT0Calibrator: public TNamed
 {
 public:

   AliT0Calibrator();
  AliT0Calibrator( const AliT0Calibrator&r );
  AliT0Calibrator& operator=(const AliT0Calibrator&r); 
  virtual ~AliT0Calibrator() {};
  //  const AliT0RecoParam* fRecoParam;     // Pointer to T0 Recon. Pars
  //   const AliT0RecoParam* GetRecoParam() const { return fRecoParam; }


  Int_t WalkCorrection(Int_t refAmp, Int_t ipmt, Int_t qt, Int_t time) ;
  void SetEq(Int_t eq) { fEqualized= eq; };
 protected:

  Int_t           fTimeDelayCFD[24];  //CFD[i]-CFD[0]
  Float_t           fMaxValue[24];  //CFD[i]-CFD[0]
  Float_t         fChannelWidth  ;   //channel width
  TObjArray       fWalk;             //walk correction function
  Int_t           fEqualized;        //if != 0 time centered around 0
      
 // const AliT0RecoParam* fRecoParam; ///< reference to reco parameters
   
  ClassDef(AliT0Calibrator, 2)   // class for the T0 reconstruction

};


#endif
