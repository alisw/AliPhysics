#ifndef AliTRDCOMMONPARAM_H
#define AliTRDCOMMONPARAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class containing constant common parameters                               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

/* $Id$ */

#include "TObject.h"

class AliTRDpadPlane;

class AliTRDCommonParam : public TObject
{

  public:
  
    enum { kNplan = 6, kNcham = 5, kNsect = 18, kNdet = 540 };
    
    AliTRDCommonParam(const AliTRDCommonParam &p);   
    AliTRDCommonParam &operator=(const AliTRDCommonParam &p); 
    virtual        ~AliTRDCommonParam();

    static AliTRDCommonParam *Instance();
    static  void    Terminate();

    virtual void    Copy(TObject &p) const;
    
    void            SetExB(Int_t exbOn = 1)                        { fExBOn             = exbOn; }
    void            SetSamplingFrequency(Float_t freq)             { fSamplingFrequency = freq;  }
    
    Bool_t          ExBOn() const                                  { return fExBOn;              }
    
    AliTRDpadPlane *GetPadPlane(Int_t p, Int_t c) const;
    Int_t           GetRowMax(Int_t p, Int_t c, Int_t /*s*/) const;
    Int_t           GetColMax(Int_t p) const;
    Double_t        GetRow0(Int_t p, Int_t c, Int_t /*s*/) const;
    Double_t        GetCol0(Int_t p) const;
    Float_t         GetSamplingFrequency() const                   { return fSamplingFrequency;  }

  protected:

    static AliTRDCommonParam *fgInstance;         //  Instance of this class (singleton implementation)
    static Bool_t             fgTerminated;       //  Defines if this class has already been terminated
    
    void Init();
    
    Int_t                     fExBOn;             //  Switch for the ExB effects

    Float_t                   fSamplingFrequency; //  Sampling Frequency in MHz
  
    TObjArray                *fPadPlaneArray;     //! Array of pad plane objects
  
  private:

    // This is a singleton, constructor is private!  
    AliTRDCommonParam();
  
    ClassDef(AliTRDCommonParam,3)                  // The constant parameters common to simulation and reconstruction
};

#endif
