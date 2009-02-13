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

#include "AliTRDSimParam.h"

class TRootIoCtor;

class AliTRDpadPlane;

class AliTRDCommonParam : public TObject
{

  public:
  
    enum { kNlayer  = 6
         , kNstack  = 5
         , kNsector = 18
         , kNdet    = 540 };

    enum { kXenon =   0
	 , kArgon =   1   };
    
    AliTRDCommonParam(TRootIoCtor *);
    AliTRDCommonParam(const AliTRDCommonParam &p);   
    AliTRDCommonParam &operator=(const AliTRDCommonParam &p); 
    virtual        ~AliTRDCommonParam();

    static AliTRDCommonParam *Instance();
    static  void    Terminate();

    virtual void    Copy(TObject &p) const;
    
    void            SetExB(Int_t exbOn = 1)                        { fExBOn             = exbOn;    }
    void            SetSamplingFrequency(Float_t freq)             { fSamplingFrequency = freq;     }
    void            SetXenon()                                     { fGasMixture        = kXenon; 
                                                                     AliTRDSimParam::Instance()->ReInit(); }
    void            SetArgon()                                     { fGasMixture        = kArgon; 
                                                                     AliTRDSimParam::Instance()->ReInit(); }

    Bool_t          ExBOn() const                                  { return fExBOn;                 }
    Bool_t          IsXenon() const                                { return (fGasMixture == kXenon) 
                                                                     ? kTRUE : kFALSE;              }
    Bool_t          IsArgon() const                                { return (fGasMixture == kArgon) 
                                                                     ? kTRUE : kFALSE;              }

    Int_t           GetGasMixture() const                          { return fGasMixture;            }
    Float_t         GetSamplingFrequency() const                   { return fSamplingFrequency;     }

    Float_t         GetOmegaTau(Float_t vdrift);
    Bool_t          GetDiffCoeff(Float_t &dl, Float_t &dt, Float_t vdrift);

    Double_t        TimeStruct(Float_t vdrift, Double_t xd, Double_t z);

  protected:

    void            SampleTimeStruct(Float_t vdrift);

    static AliTRDCommonParam *fgInstance;          //  Instance of this class (singleton implementation)
    static Bool_t             fgTerminated;        //  Defines if this class has already been terminated    

    Int_t                     fExBOn;              //  Switch for the ExB effects

    Float_t                   fDiffusionT;         //  Transverse drift coefficient
    Float_t                   fDiffusionL;         //  Longitudinal drift coefficient
    Float_t                   fDiffLastVdrift;     //  The structures are valid for fLastVdrift (caching)

    Float_t                  *fTimeStruct1;        //! Time Structure of Drift Cells
    Float_t                  *fTimeStruct2;        //! Time Structure of Drift Cells
    Float_t                   fVDlo;               //  Lower drift velocity, for interpolation
    Float_t                   fVDhi;               //  Higher drift velocity, for interpolation
    Float_t                   fTimeLastVdrift;     //  The structures are valid for fLastVdrift (caching)

    Float_t                   fSamplingFrequency;  //  Sampling Frequency in MHz

    Int_t                     fGasMixture;         //  Gas mixture: 0-Xe/C02 1-Ar/CO2. 
  
  private:

    // This is a singleton, constructor is private!  
    AliTRDCommonParam();
  
    ClassDef(AliTRDCommonParam,7)                  // The constant parameters common to simulation and reconstruction

};
#endif
