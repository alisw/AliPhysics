#ifndef ALIT0CALIBRATOR_H
#define ALIT0CALIBRATOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliT0CalibData.h"
#include "TNamed.h"

class AliT0Calibrator: public TNamed
 {
 public:

  AliT0Calibrator();
  AliT0Calibrator( const AliT0Calibrator& );
  AliT0Calibrator& operator=(const AliT0Calibrator&); 
  virtual ~AliT0Calibrator() {};
  

  Int_t WalkCorrection(Int_t ipmt, Int_t qt, Int_t time) ;
  //  Int_t EquivalizeChannel(Int_t ipmt)  ;
 protected:

  Int_t           fTimeDelayCFD[24]; 
  Float_t           fChannelWidth  ;  
  TObjArray       fWalk;
  
  ClassDef(AliT0Calibrator, 1)   // class for the T0 reconstruction

};


#endif
