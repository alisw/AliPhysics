#ifndef ALIMAGFC_H
#define ALIMAGFC_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//     Constant magnetic field class
//     Used by AliRun class
//     Author:
//-------------------------------------------------------------------------

#include "AliMagF.h"

class AliMagFC  : public AliMagF
{
  //Alice Constant Magnetic Field

public:
  AliMagFC(){}
  AliMagFC(const char *name, const char *title, Int_t integ, 
	   Float_t factor, Float_t fmax);
  virtual ~AliMagFC() {}
  virtual void Field(Float_t *x, Float_t *b);
  virtual void ReadField() {}
  virtual void ZDCField(Float_t *x, Float_t *b);
  ClassDef(AliMagFC,1)  //Class for all Alice Constant MagField 
};


//ZDC part -------------------------------------------------------------------

// ************************ LHC optics v6.4 *****************************
static const Float_t kG1=20.443;
static const Float_t kFDIP=-37.85;
static const Float_t kFCORN2=-9.6979; 
//
// ZBEG       Beginning of the inner triplet
// D1BEG      Beginning of separator dipole 1
// D2BEG      Beginning of separator dipole 2
// CORBEG     Corrector dipole beginning (because of dimuon arm)
//
static const Float_t kCORBEG2 = -1972.5,kCOREND2 = kCORBEG2 - 153., kCOR2RA2 = 4.5 * 4.5;
//
static const Float_t kZBEG  = -2296.5;
static const Float_t kZ1BEG = kZBEG +   0.,   kZ1END = kZ1BEG - 637.,kZ1RA2 = 3.5 * 3.5;
static const Float_t kZ2BEG = kZBEG - 908.5,  kZ2END = kZ2BEG - 550.,kZ2RA2 = 3.5 * 3.5;
static const Float_t kZ3BEG = kZBEG - 1558.5, kZ3END = kZ3BEG - 550.,kZ3RA2 = 3.5 * 3.5;
static const Float_t kZ4BEG = kZBEG - 2430.,  kZ4END = kZ4BEG - 637.,kZ4RA2 = 3.5 * 3.5;
static const Float_t kD1BEG = - 5838.3    ,kD1END = kD1BEG - 945., kD1RA2 = 4.5 * 4.5;
static const Float_t kD2BEG = - 12167.8   ,kD2END = kD2BEG - 945., kD2RA2 = 4.5 * 4.5;
//
static const Float_t kXCEN1D2 = -9.7     ,kYCEN1D2 = 0.;
static const Float_t kXCEN2D2 =  9.7     ,kYCEN2D2 = 0.;

//ZDC part -------------------------------------------------------------------

#endif
