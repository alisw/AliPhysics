#ifndef ALIDECAYERPYTHIA_H
#define ALIDECAYERPYTHIA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

// Implementation of AliDecayer using Pythia
// Method forwarding to the AliPythia instance.
// Author: andreas.morsch@cern.ch

#include "AliDecayer.h"
class AliPythia;
class TClonesArrray;

class AliDecayerPythia :
public AliDecayer
{
 public:
    AliDecayerPythia();
    AliDecayerPythia(const AliDecayerPythia &decayer);
    //    
    virtual         ~AliDecayerPythia(){;}
    virtual void    Init();
    virtual void    Decay(Int_t idpart, TLorentzVector *p);
    virtual Int_t   ImportParticles(TClonesArray *particles);
    virtual void    SetForceDecay(Decay_t decay) {fDecay=decay;}
    virtual void    SetForceDecay(Int_t decay)
      {SetForceDecay((Decay_t) decay);}
    virtual void    ForceDecay();
    virtual void    SetPatchOmegaDalitz() {fPatchOmegaDalitz = 1;}
    virtual void    SetDecayerExodus()    {fDecayerExodus = 1;}
    virtual void    HeavyFlavourOff() {fHeavyFlavour = kFALSE;}
    virtual void    DecayLongLivedParticles()  {fLongLived    = kTRUE;}
    virtual Float_t GetPartialBranchingRatio(Int_t ipart);
    virtual Float_t GetLifetime(Int_t kf);
    virtual void    SwitchOffBDecay();
    virtual void    SwitchOffPi0() {fPi0 = 0;}
    virtual void    SwitchOffParticle(Int_t kf);
    virtual void    WriteDecayTable();
    virtual void    ReadDecayTable();
    
 private:
    void     Lu1Ent(Int_t flag, Int_t idpart, 
		    Double_t mom, Double_t theta, Double_t phi);
    Int_t    CountProducts(Int_t channel, Int_t particle);
    void     ForceParticleDecay(Int_t particle, Int_t product, Int_t mult);
    void     ForceParticleDecay(Int_t particle, const Int_t* products, Int_t* mult, Int_t npart, Bool_t flag = 0);
    void     ForceHadronicD(Int_t optUse4Bodies=1, Int_t optUseDtoV0=0);
    void     ForceOmega();
    void     ForceLambda();
    void     SwitchOffHeavyFlavour();
    void     ForceBeautyUpgrade();
    Float_t  GetBraPart(Int_t kf);
    void     Copy(TObject &decayer) const;

    AliDecayerPythia &operator=(const AliDecayerPythia &decayer) 
    {decayer.Copy(*this);return(*this);}
    
    
 private:
    AliPythia*  fPythia;          //! Pointer to AliPythia
    Decay_t     fDecay;           //  Forced decay mode
    Float_t     fBraPart[501];    //! Branching ratios
    Bool_t      fHeavyFlavour;    //! Flag for heavy flavors
    Bool_t      fLongLived;       //! Flag for long lived particle decay
    Bool_t      fPatchOmegaDalitz;//! Flag to patch the omega Dalitz decays 
    Bool_t      fDecayerExodus;  //! Flag for EXODUS decayer
    Bool_t      fPi0;             //! Flag for pi0 decay 
    static Bool_t fgInit;         //! initialization flag 
    
    ClassDef(AliDecayerPythia, 5) // AliDecayer implementation using Pythia  
};
#endif







