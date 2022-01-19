/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */
#ifndef AliMCGenParticleContainer_H
#define AliMCGenParticleContainer_H

#include "AliMCParticle.h"
#include "AliPID.h"

class AliMCGenParticleContainer : public TObject
{
public:
    AliMCGenParticleContainer(); // default constructor
    virtual ~AliMCGenParticleContainer();//destructor.
    
    void Fill(AliMCParticle* particle, Int_t label, Int_t lMotherPDG, Bool_t IsPrimary);
    void Reset();
    
    //setter 
    
    
    //getter
    //-----------BASIC-INFO---------------------------
    
    Int_t       GetPID() const; // get MC PID
    Double_t    GetMass() const; // get mass
    
    Bool_t      GetIsPhysicalPrimary()  const { return  fTreeMCVarIsPhysicalPrimary  ;}
    
    Int_t       GetPdgCode()            const { return  fTreeMCVarPdgCode  ;}
    Int_t       GetMotherPdgCode()      const { return  fTreeMCVarMotherPdgCode  ;}
    Int_t       GetLabel()              const { return  fTreeMCVarLabel  ;}
    
    Int_t       GetCharge()             const { return  fTreeMCVarCharge  ;}
    Double_t    GetRap()                const { return  fTreeMCVarRap  ;}
    Double_t    GetEta()                const { return  fTreeMCVarEta  ;}
    Double_t    GetPtot()               const { return  fTreeMCVarPtot  ;}
    Double_t    GetPt()                 const { return  fTreeMCVarPt  ;}
    Double_t    GetPx()                 const { return  fTreeMCVarPx  ;}
    Double_t    GetPy()                 const { return  fTreeMCVarPy  ;}
    Double_t    GetPz()                 const { return  fTreeMCVarPz  ;}
    Double_t    GetPhi()                const { return  fTreeMCVarPhi  ;}
    
private:
    Bool_t      fTreeMCVarIsPhysicalPrimary;//
    
    Int_t       fTreeMCVarPdgCode;//
    Int_t       fTreeMCVarMotherPdgCode;//
    Int_t       fTreeMCVarLabel;//
    
    Int_t       fTreeMCVarCharge;//
    Double_t    fTreeMCVarRap;//
    Double_t    fTreeMCVarEta;//
    Double_t    fTreeMCVarPtot;//
    Double_t    fTreeMCVarPt;//
    Double_t    fTreeMCVarPx;//
    Double_t    fTreeMCVarPy;//
    Double_t    fTreeMCVarPz;//
    Double_t    fTreeMCVarPhi;//
    
    //------------------------------------------------
    ClassDef(AliMCGenParticleContainer, 1);
};
#endif
