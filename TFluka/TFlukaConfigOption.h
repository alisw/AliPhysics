#ifndef TFLUKACONFIGOPTION
#define TFLUKACONFIGOPTION

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                                                                           //
// Class to store FLUKA and VMC configuration options:                       // 
// Cuts, Processes, User Scoring                                             // 
//                                                                           //
//                                                                           //
// Author: andreas.morsch@cern.ch                                            // 
//                                                                           //  
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
//
//
//


#include <TNamed.h>


class TFlukaConfigOption : public TNamed
{
public:
    // Constructors
    TFlukaConfigOption();
    TFlukaConfigOption(const char* cutName, Double_t cut);
    TFlukaConfigOption(const char* cutName, Double_t cut, Int_t imed);
    TFlukaConfigOption(const char* procName, Int_t flag);
    TFlukaConfigOption(const char* procName, Int_t flag, Int_t imed);
    // Getters
    Double_t Cut()  const                {return fCutValue;}
    Int_t    Flag() const                {return fProcessFlag;}
    Int_t    Medium()                    {return fMedium;}    
    Bool_t   HasMediumAssigned()         {return  (fMedium > -1);}
    // Setters
    void     SetCut(Double_t val)        {fCutValue     =   val;}
    void     SetFlag(Int_t val)          {fProcessFlag  =   val;}
    void     SetMedium(Int_t imed)       {fMedium       =   imed;}    

 protected:
    Double_t fCutValue;                // User cut
    Int_t    fProcessFlag;             // User flag assigned to processes
    Int_t    fMedium;                  // Materials assigned to user settings
    ClassDef(TFlukaConfigOption, 1)    // Fluka Configuration Option
};
	
#endif
	
