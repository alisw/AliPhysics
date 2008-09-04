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

class TRootIoCtor;

class AliTRDpadPlane;

class AliTRDCommonParam : public TObject
{

  public:
  
    enum { kNlayer = 6, kNstack = 5, kNsector = 18, kNdet = 540 };
    
    AliTRDCommonParam(TRootIoCtor *);
    AliTRDCommonParam(const AliTRDCommonParam &p);   
    AliTRDCommonParam &operator=(const AliTRDCommonParam &p); 
    virtual        ~AliTRDCommonParam();

    static AliTRDCommonParam *Instance();
    static  void    Terminate();

    virtual void    Copy(TObject &p) const;
    
    void            SetExB(Int_t exbOn = 1)                        { fExBOn             = exbOn;  }
    void            SetSamplingFrequency(Float_t freq)             { fSamplingFrequency = freq;   }

    Bool_t          ExBOn() const                                  { return fExBOn;               }
    Float_t         GetSamplingFrequency() const                   { return fSamplingFrequency;   }

  protected:

    void            Init();

    static AliTRDCommonParam *fgInstance;          //  Instance of this class (singleton implementation)
    static Bool_t             fgTerminated;        //  Defines if this class has already been terminated    
    Int_t                     fExBOn;              //  Switch for the ExB effects
    Float_t                   fSamplingFrequency;  //  Sampling Frequency in MHz
  
  private:

    // This is a singleton, constructor is private!  
    AliTRDCommonParam();
  
    ClassDef(AliTRDCommonParam,6)                  // The constant parameters common to simulation and reconstruction

};

#endif
