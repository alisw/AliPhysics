#ifndef ALIG3Medium_H
#define ALIG3Medium_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "TNamed.h"

const static Int_t kNPars=33;

class AliG3Medium : public TNamed
{
public:
    AliG3Medium();
    AliG3Medium(Int_t imed, Int_t imat, const char* name, Int_t isvol,
		 Int_t ifield, Float_t fieldm, Float_t tmaxfd,
		 Float_t stemax, Float_t deemax,
		 Float_t epsil, Float_t stmin);
    
    virtual ~AliG3Medium(){;}
    // Dump medium parameters
    virtual void    Dump();
    // Get id
    virtual Int_t   Id();
    // Get parameters
    virtual Int_t   IdMat()   {return fIdMat;}
    virtual Int_t   Isvol()   {return fIsvol;}
    virtual Int_t   Ifield()  {return fIfield;}
    virtual Float_t Fieldm()  {return fFieldm;}    
    virtual Float_t Tmaxfd()  {return fTmaxfd;}        
    virtual Float_t Stemax()  {return fStemax;}    
    virtual Float_t Deemax()  {return fDeemax;}        
    virtual Float_t Epsil()   {return fEpsil;}
    virtual Float_t Stmin()   {return fStmin;}
    virtual void    SetPar(Int_t ipar, Float_t par) {fPars[ipar-1]=par;}
    virtual Float_t GetPar(Int_t ipar);
    // Set and get link to widget entry
    virtual Int_t ItemId() {return fItem;}
    virtual void  SetItemId(Int_t id) {fItem=id;}
    
 private:
    Float_t fPars[kNPars];   // special medium parameters
    Int_t   fId;             // Id number of the Medium
    Int_t   fIdMat;          // Associated material
    Int_t   fIsvol;          // Sensitivity flag 
    Int_t   fIfield;         // Magnetic Field Flag
    Float_t fFieldm;         // Maximum Field Strength
    Float_t fTmaxfd;         // Max. Ang. Deviation
    Float_t fStemax;         // Maximum Step   
    Float_t fDeemax;         // Max. Frac. Energy Loss",
    Float_t fEpsil;          // Crossing Precission 
    Float_t fStmin;          // Minimum Step Size
    //
    Int_t   fItem;           // Link to Widget Entry

  AliG3Medium(const AliG3Medium&) {}
  AliG3Medium & operator=(const AliG3Medium&) {return *this;}

    ClassDef(AliG3Medium,1) // Tracking Medium Object for GUI 
};

#endif
