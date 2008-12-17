#ifndef ALIMAGFC_H
#define ALIMAGFC_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//     Constant magnetic field class
//     Used by AliRun class
//     Author:
//-------------------------------------------------------------------------

#include "AliMagF.h"

enum BeamType_t {kBeamTypeAA, kBeamTypepp};
class AliMagFC  : public AliMagF
{
  //Alice Constant Magnetic Field

public:
  AliMagFC();
  AliMagFC(const char *name, const char *title, Int_t integ, 
	   Float_t factor, Float_t fmax);
  virtual ~AliMagFC(){}
  virtual void Field(const float *x, float *b)      const;
  virtual void Field(const double *x, double *b)    const;
  virtual void ReadField() {}
  virtual void ZDCField(const float *x, float *b)   const;
  virtual void ZDCField(const double *x, double *b) const;
  virtual void SetBeamType(BeamType_t type)      {fBeamType    = type;}
  virtual void SetBeamEnergy(Float_t energy)     {fBeamEnergy  = energy;}
  virtual void SetCompensatorMagnet(Bool_t flag) {fCompensator = flag;}

private:
  Bool_t     fCompensator; // Flag for compensator magnetic field (kTrue -> ON)
  BeamType_t fBeamType;    // Beam type: A-A (fBeamType=0) or p-p (fBeamType=1)
  Float_t    fBeamEnergy;  // Beam energy in GeV
  mutable Float_t    fQuadGradient;// Gradient field for inner triplet quadrupoles
  mutable Float_t    fDipoleField; // Field value for D1 and D2 dipoles
  mutable Float_t    fCCorrField;  // Side C 2nd compensator field
  mutable Float_t    fACorr1Field; // Side A 1st compensator field 
  mutable Float_t    fACorr2Field; // Side A 2nd compensator field
  
  ClassDef(AliMagFC,3)  //Class for all Alice Constant MagField 
};


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  ZDC part  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// ************************ LHC optics v6.5 *****************************
// ---- Position of the magnetic elements of LHC beam optics ----
// -> SIDE C
static const Float_t kCCorrBegin = -1972.5, kCCorrEnd = kCCorrBegin - 153., kCCorrSqRadius = 4.5*4.5;
//
static const Float_t kCTripletBegin  = -2296.5;
static const Float_t kCQ1Begin = kCTripletBegin,        kCQ1End = kCQ1Begin-637., kCQ1SqRadius = 3.5*3.5;
static const Float_t kCQ2Begin = kCTripletBegin-908.5,  kCQ2End = kCQ2Begin-550., kCQ2SqRadius = 3.5*3.5;
static const Float_t kCQ3Begin = kCTripletBegin-1558.5, kCQ3End = kCQ3Begin-550., kCQ3SqRadius = 3.5*3.5;
static const Float_t kCQ4Begin = kCTripletBegin-2400.,  kCQ4End = kCQ4Begin-637., kCQ4SqRadius = 3.5*3.5;
//
static const Float_t kCD1Begin = -5838.3,  kCD1End = kCD1Begin-945., kCD1SqRadius = 4.5*4.5;
static const Float_t kCD2Begin = -12167.8, kCD2End = kCD2Begin-945., kCD2SqRadius = 4.5*4.5;
static const Float_t kCD2XCentre1 = -9.7;
static const Float_t kCD2XCentre2 =  9.7;
//
// -> SIDE A
// NB -> kACorr1Begin = 919. to be checked
static const Float_t kACorr1Begin = 919., kACorr1End = kACorr1Begin+260., kCCorr1SqRadius = 4.*4.;
static const Float_t kACorr2Begin = 1972.5, kACorr2End = kACorr2Begin+153., kCCorr2SqRadius = 4.5*4.5;
static const Float_t kATripletBegin  = 2296.5;
static const Float_t kAQ1Begin = kATripletBegin,	kAQ1End = kAQ1Begin+637., kAQ1SqRadius = 3.5*3.5;
static const Float_t kAQ2Begin = kATripletBegin+908.5,  kAQ2End = kAQ2Begin+550., kAQ2SqRadius = 3.5*3.5;
static const Float_t kAQ3Begin = kATripletBegin+1558.5, kAQ3End = kAQ3Begin+550., kAQ3SqRadius = 3.5*3.5;
static const Float_t kAQ4Begin = kATripletBegin+2400.,  kAQ4End = kAQ4Begin+637., kAQ4SqRadius = 3.5*3.5;
//
static const Float_t kAD1Begin = 5838.3,  kAD1End = kAD1Begin+945., kAD1SqRadius = 3.375*3.375;
static const Float_t kAD2Begin = 12167.8, kAD2End = kAD2Begin+945., kAD2SqRadius = 3.75*3.75;
static const Float_t kAD2XCentre1 = -9.4;
static const Float_t kAD2XCentre2 =  9.4;

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  ZDC part  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#endif

