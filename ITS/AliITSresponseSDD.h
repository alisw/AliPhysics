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

    virtual void SetTimeOffset(Float_t to){fTimeOffset = to;}
    virtual Float_t GetTimeOffset()const {return fTimeOffset;}
    virtual Float_t GetTimeZero(Int_t modIndex){
      if(modIndex<kNSPDmods || modIndex>kNSPDmods+kNSDDmods){
	AliError(Form("SDD module number %d out of range",modIndex));
	return 0.;
      }
      return fTimeZero[modIndex-kNSPDmods];
    }
    static Float_t DefaultTimeOffset() {return fgkTimeOffsetDefault;}

    virtual void SetADC2keV(Float_t conv){fADC2keV=conv;}
    virtual Float_t GetADC2keV()const {return fADC2keV;}
    static Float_t DefaulttADC2keV() {return fgkADC2keVDefault;}

    static Float_t GetCarlosRXClockPeriod() {return fgkCarlosRXClockPeriod;}
 

 protected:

    enum {kNSPDmods = 240};
    enum {kNSDDmods = 260};
    enum {kNLaddersLay3 = 14};
    enum {kNLaddersLay4 = 22};


    static const Float_t fgkTimeOffsetDefault; // default for fTimeOffset
    static const Float_t fgkADC2keVDefault; // default for fADC2keV
    static const Float_t fgkCarlosRXClockPeriod;  // clock period for CarlosRX

    Float_t  fTimeOffset;          // Time offset due to electronic delays 
                                   // --> obsolete, kept for backw. comp. 
    Float_t  fTimeZero[kNSDDmods]; // Time Zero for each module
    Float_t  fADC2keV;             // Conversion factor from ADC to keV

 private:

   AliITSresponseSDD(const AliITSresponseSDD &ob); // copy constructor
   AliITSresponseSDD& operator=(const AliITSresponseSDD & /* source */); // ass. op.

    ClassDef(AliITSresponseSDD,16) 
    
    };
#endif
