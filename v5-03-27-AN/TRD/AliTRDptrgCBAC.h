#ifndef AlITRDPTRGCBAC_H
#define AliTRDPTRGCBAC_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

// --------------------------------------------------------
// 
// Pre-Trigger Control-Box A or C class
//
// --------------------------------------------------------

#include "AliTRDptrgParam.h"
#include <TObjArray.h>
#include <TObject.h>

class AliRunLoader;


class AliTRDptrgCBAC : public TObject {
 public:
  AliTRDptrgCBAC(AliRunLoader *rl = 0x0);
  AliTRDptrgCBAC(AliRunLoader *rl, AliTRDptrgFEBPosition_t position,
                 AliTRDptrgOperatingMode_t operatingMode, 
                 AliTRDptrgParam *param);
  ~AliTRDptrgCBAC();
  
  Int_t* Simulate();

 protected:
  Bool_t LoadParams(); // load AliTRDprtgParam content

  AliRunLoader *fRunLoader;  //!
  TObjArray fLUTArray; // Array with Look-Up-Tables
  TObjArray fFEBArray; // front end boxes connected to T0 (fFEB[0]) and VO (4x)
  AliTRDptrgFEBPosition_t fPosition; // Control box position (A or C side)
  AliTRDptrgOperatingMode_t fOperatingMode; // working on Digits or Hits?
  AliTRDptrgParam* fParam; // parameters
 private:
  AliTRDptrgCBAC& operator=(const AliTRDptrgCBAC &rhs); // not implemented
  AliTRDptrgCBAC(const AliTRDptrgCBAC &rhs); // not implemented		 

  ClassDef(AliTRDptrgCBAC, 1);
};

#endif
