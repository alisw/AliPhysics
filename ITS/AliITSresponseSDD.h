#ifndef ALIITSRESPONSESDD_H
#define ALIITSRESPONSESDD_H
 
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 

#include <TObject.h>
#include <AliLog.h>

/* $Id$ */

/////////////////////////////////////////////////////////////
//  Base settings for the ITS response classes.            //  
//  The data member of this class are static and set once  //
//  for all the modules.                                   //    
///////////////////////////////////////////////////////////// 

class AliITSresponseSDD : public TObject {
 public:
  enum {kVDCorr2Side = BIT(14),kVDCorrMult = BIT(15)};   // if bit set, the object contains separate corrections for 2 sides
  //
  AliITSresponseSDD();
  virtual ~AliITSresponseSDD(){};

  virtual void SetSideATimeZero(Float_t tzero){
    SetLayer3ATimeZero(tzero);
    SetLayer4ATimeZero(tzero);
  }
  virtual void SetSideCTimeZero(Float_t tzero){
    SetLayer3CTimeZero(tzero);
    SetLayer4CTimeZero(tzero);
  }
  virtual void SetLayer3ATimeZero(Float_t tzero){
    for(Int_t iLad=1; iLad<=kNLaddersLay3; iLad++) SetHalfLadderATimeZero(3,iLad,tzero);      
  }
  virtual void SetLayer3CTimeZero(Float_t tzero){
    for(Int_t iLad=1; iLad<=kNLaddersLay3; iLad++) SetHalfLadderCTimeZero(3,iLad,tzero);
  }
  virtual void SetLayer4ATimeZero(Float_t tzero){
    for(Int_t iLad=1; iLad<=kNLaddersLay4; iLad++) SetHalfLadderATimeZero(4,iLad,tzero);      
  }
  virtual void SetLayer4CTimeZero(Float_t tzero){
    for(Int_t iLad=1; iLad<=kNLaddersLay4; iLad++) SetHalfLadderCTimeZero(4,iLad,tzero);
  }
  virtual void SetHalfLadderATimeZero(Int_t lay, Int_t lad, Float_t tzero);
  virtual void SetHalfLadderCTimeZero(Int_t lay, Int_t lad, Float_t tzero);
  virtual void SetModuleTimeZero(Int_t modIndex, Float_t tzero){
    if(modIndex<kNSPDmods || modIndex>kNSPDmods+kNSDDmods) AliError(Form("SDD module number %d out of range",modIndex));
    fTimeZero[modIndex-kNSPDmods]=tzero;
  }

  virtual void SetDeltaVDrift(Int_t modIndex, Float_t dv, Bool_t rightSide=kFALSE) {
    int ind = GetVDIndex(modIndex,rightSide);
    if (ind>=0) fDeltaVDrift[ind] = dv;
  }

  virtual Float_t GetDeltaVDrift(Int_t modIndex,Bool_t rightSide=kFALSE) const {
    int ind = GetVDIndex(modIndex,rightSide);
    return ind<0 ? 0.:fDeltaVDrift[ind];
  }
  // 
  Bool_t IsVDCorr2Side()                       const {return TestBit(kVDCorr2Side);}
  Bool_t IsVDCorrMult()                        const {return TestBit(kVDCorrMult);}
  void   SetVDCorr2Side(Bool_t v=kTRUE)              {SetBit(kVDCorr2Side,v);}
  void   SetVDCorrMult(Bool_t v=kTRUE)               {SetBit(kVDCorrMult,v);}
  //
  static Float_t DefaultTimeOffset() {return fgkTimeOffsetDefault;}
  virtual void SetTimeOffset(Float_t to){fTimeOffset = to;}
  virtual Float_t GetTimeOffset()const {return fTimeOffset;}
  virtual Float_t GetTimeZero(Int_t modIndex) const {
    if(modIndex<kNSPDmods || modIndex>kNSPDmods+kNSDDmods){
      AliError(Form("SDD module number %d out of range",modIndex));
      return 0.;
    }
    return fTimeZero[modIndex-kNSPDmods];
  }

  virtual void SetADC2keV(Float_t conv){fADC2keV=conv;}
  virtual Float_t GetADC2keV()const {return fADC2keV;}
  virtual void SetADCtokeV(Int_t modIndex, Float_t conv){
    if(modIndex<kNSPDmods || modIndex>kNSPDmods+kNSDDmods) AliError(Form("SDD module number %d out of range",modIndex));
    fADCtokeV[modIndex-kNSPDmods]=conv;
  }
  virtual Float_t GetADCtokeV(Int_t modIndex) const {
    if(modIndex<kNSPDmods || modIndex>kNSPDmods+kNSDDmods){
      AliError(Form("SDD module number %d out of range",modIndex));
      return 0.;
    }
    return fADCtokeV[modIndex-kNSPDmods];
  }

  virtual void SetChargevsTime(Float_t slope){fChargevsTime=slope;}
  virtual Float_t GetChargevsTime()const {return fChargevsTime;}

  static Float_t DefaultADC2keV() {return fgkADC2keVDefault;}
  static Float_t DefaultChargevsTime() {return fgkChargevsTimeDefault;}

  static Float_t GetCarlosRXClockPeriod() {return fgkCarlosRXClockPeriod;}
  void PrintChargeCalibrationParams() const;
  void PrintTimeZeroes() const;
  void PrintVdriftCorerctions() const;


 protected:
    //
  virtual Int_t GetVDIndex(Int_t modIndex, Bool_t rightSide=kFALSE) const {
    int ind = modIndex - kNSPDmods;
    if(ind<0 || ind>=kNSDDmods) {AliError(Form("SDD module number %d out of range",modIndex)); return -1;}
    return (rightSide && IsVDCorr2Side()) ? ind + kNSDDmods : ind;
  }


 protected:

  enum {kNSPDmods = 240};
  enum {kNSDDmods = 260};
  enum {kNLaddersLay3 = 14};
  enum {kNLaddersLay4 = 22};


  static const Float_t fgkTimeOffsetDefault;   // default for fTimeOffset
  static const Float_t fgkADC2keVDefault;      // default for fADC2keV
  static const Float_t fgkChargevsTimeDefault; // default for fChargevsTime
  static const Float_t fgkCarlosRXClockPeriod; // clock period for CarlosRX

  Float_t  fTimeOffset;             // Time offset due to electronic delays 
                                    // --> obsolete, kept for backw. comp. 
  Float_t  fTimeZero[kNSDDmods];    // Time Zero for each module
  Float_t  fDeltaVDrift[2*kNSDDmods];  // Vdrift correction (um/ns) for each module left (<kNSDDmods) and right (>=kNSDDmods) sides
  Float_t  fADC2keV;                // Conversion factor from ADC to keV
                                    // --> obsolete, kept for backw. comp. 
  Float_t  fChargevsTime;           // Correction for zero suppression effect
  Float_t  fADCtokeV[kNSDDmods]; // ADC to keV conversion for each module

 private:

  AliITSresponseSDD(const AliITSresponseSDD &ob); // copy constructor
  AliITSresponseSDD& operator=(const AliITSresponseSDD & /* source */); // ass. op.

  ClassDef(AliITSresponseSDD,20) 
     
    };
#endif
