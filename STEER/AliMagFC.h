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
  AliMagFC():AliMagF(),fCompensator(kFALSE){}
  AliMagFC(const char *name, const char *title, Int_t integ, 
	   Float_t factor, Float_t fmax);
  virtual ~AliMagFC(){}
  virtual void Field(Float_t *x, Float_t *b) const;
  virtual void ReadField() {}
  virtual void ZDCField(Float_t *x, Float_t *b) const;
  virtual void SetCompensatorMagnet(Bool_t flag) {fCompensator = flag;}
 private:
  Bool_t  fCompensator; // Flag for compensator magnetic field (kTrue -> ON)
  ClassDef(AliMagFC,2)  //Class for all Alice Constant MagField 
};


//ZDC part -------------------------------------------------------------------

// ************************ LHC optics v6.5 *****************************
static const Float_t kG1=20.707;
static const Float_t kFDIP=-37.71;
static const Float_t kFCORN2=-9.667; 
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


// ************************ Left line *****************************************
// ************************ LHC optics v6.5 (official version for ions)****************
static const Float_t kG1l=20.707;
static const Float_t kFDIPl=-37.71;
static const Float_t kFCORN2l=-11.72; 
//
//
// ZBEG       Beginning of the inner triplet
// D1BEG      Beginning of separator dipole 1
// D2BEG      Beginning of separator dipole 2
// CORBEG     Corrector dipole beginning (because of dimuon arm)
//
static const Float_t kCORBEG2l = 1972.5,kCOREND2l = kCORBEG2l + 153., kCOR2RA2l = 4.5 * 4.5;// second corrector
static const Float_t kZBEGl  = 2296.5;// inner triplet beginning
static const Float_t kZ1BEGl = kZBEGl +   0.,   kZ1ENDl = kZ1BEGl + 637.,kZ1RA2l = 3.5 * 3.5;// Q1
static const Float_t kZ2BEGl = kZBEGl + 908.5,  kZ2ENDl = kZ2BEGl + 550.,kZ2RA2l = 3.5 * 3.5;// Q2A
static const Float_t kZ3BEGl = kZBEGl + 1558.5, kZ3ENDl = kZ3BEGl + 550.,kZ3RA2l = 3.5 * 3.5;// Q2B
static const Float_t kZ4BEGl = kZBEGl + 2400.,  kZ4ENDl = kZ4BEGl + 637.,kZ4RA2l = 3.5 * 3.5;// Q3
static const Float_t kD1BEGl = 5838.3	 ,kD1ENDl = kD1BEGl + 945., kD1RA2l = 3.375 * 3.375;// D1
static const Float_t kD2BEGl = 12167.8   ,kD2ENDl = kD2BEGl + 945., kD2RA2l = 3.75 * 3.75;// D2
static const Float_t kXCEN1D2l = -9.4	  ,kYCEN1D2l = 0.;// D2
static const Float_t kXCEN2D2l =  9.4	  ,kYCEN2D2l = 0.;// D2


/*
// ************************ LHC optics v6.500 (official version for pp 7TeV) **********
static const Float_t kG1l=22.000;
static const Float_t kFDIPl=-37.804;
//static const Float_t kFCORN1l=13.201; 
static const Float_t kFCORN2l=-11.751;
//
static const Float_t kCORBEG1l = 945.,kCOREND1l = kCORBEG1l + 260., kCOR1RA2l = 4.5 * 4.5;// first corrector 
static const Float_t kCORBEG2l = 1972.5,kCOREND2l = kCORBEG2l + 153., kCOR2RA2l = 4.5 * 4.5;// second corrector
static const Float_t kZBEGl  = 2296.5;// inner triplet beginning
static const Float_t kZ1BEGl = kZBEGl +   0.,	kZ1ENDl = kZ1BEGl + 637.,kZ1RA2l = 3.5 * 3.5;// Q1
static const Float_t kZ2BEGl = kZBEGl + 908.5,  kZ2ENDl = kZ2BEGl + 550.,kZ2RA2l = 3.5 * 3.5;// Q2A
static const Float_t kZ3BEGl = kZBEGl + 1558.5, kZ3ENDl = kZ3BEGl + 550.,kZ3RA2l = 3.5 * 3.5;// Q2B
static const Float_t kZ4BEGl = kZBEGl + 2400.,  kZ4ENDl = kZ4BEGl + 637.,kZ4RA2l = 3.5 * 3.5;// Q3
static const Float_t kD1BEGl = 5838.3	 ,kD1ENDl = kD1BEGl + 945., kD1RA2l = 3.375 * 3.375;// D1
static const Float_t kD2BEGl = 12167.8   ,kD2ENDl = kD2BEGl + 945., kD2RA2l = 3.75 * 3.75;// D2
static const Float_t kXCEN1D2l = -9.4	  ,kYCEN1D2l = 0.;// D2
static const Float_t kXCEN2D2l =  9.4	  ,kYCEN2D2l = 0.;// D2
*/

//ZDC part -------------------------------------------------------------------

#endif
