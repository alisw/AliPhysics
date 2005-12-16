#ifndef AliTRDCOMMONPARAM_H
#define AliTRDCOMMONPARAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class containing constant common parameters                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

/* $Id$ */

#include "TObject.h"

class AliTRDpadPlane;

class AliTRDCommonParam : public TObject
{
  public:
    static AliTRDCommonParam* Instance();
    static void Terminate();
  
    enum { kNplan = 6, kNcham = 5, kNsect = 18, kNdet = 540 };
    
    AliTRDCommonParam(const AliTRDCommonParam &p);   
    AliTRDCommonParam &operator=(const AliTRDCommonParam &p); 
    virtual void Copy(TObject &p) const;
    
    void SetField(Float_t field)                        { fField          = field;    };
    
    virtual void     SetExB(Int_t exbOn = 1)                        { fExBOn          = exbOn;    };
    
    virtual void     SetPadResponse(Int_t prfOn = 1)                { fPRFOn          = prfOn;    };
    
    Float_t  GetField()                               const { return fField; };
    
    Bool_t   ExBOn()                                  const { return fExBOn;         };
    
    Bool_t   PRFOn()                                  const { return fPRFOn;         };
    virtual Int_t     PadResponse(Double_t signal, Double_t dist, Int_t plane, Double_t *pad) const;
    
    virtual AliTRDpadPlane *GetPadPlane(Int_t p, Int_t c) const;
    Int_t    GetRowMax(Int_t p, Int_t c, Int_t /*s*/) const;
    Int_t    GetColMax(Int_t p) const;
    Double_t GetRow0(Int_t p, Int_t c, Int_t /*s*/) const;
    Double_t GetCol0(Int_t p) const;
  
  protected:
    static AliTRDCommonParam* fgInstance;     // Instance of this class (singleton implementation)
    static Bool_t fgTerminated;               // Defines if this class has already been terminated and therefore does not return instances in GetInstance anymore
    
    void Init();
    void SamplePRF();
    
    Float_t              fField;                              //  Magnetic field
    
    Int_t                fExBOn;                              //  Switch for the ExB effects
    
    Int_t                fPRFOn;                              //  Switch for the pad response
    Float_t             *fPRFsmp;                             //! Sampled pad response
    Int_t                fPRFbin;                             //  Number of bins for the PRF
    Float_t              fPRFlo;                              //  Lower boundary of the PRF
    Float_t              fPRFhi;                              //  Higher boundary of the PRF
    Float_t              fPRFwid;                             //  Bin width of the sampled PRF
    Int_t                fPRFpad;                             //  Distance to next pad in PRF
  
    TObjArray  *fPadPlaneArray;                      //  Array of pad plane objects
  
  private:
    // this is a singleton, constructor is private!  
    AliTRDCommonParam();
    ~AliTRDCommonParam();
  
    ClassDef(AliTRDCommonParam, 0)
};

#endif
