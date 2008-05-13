/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               *
 **************************************************************************/

//-------------------------------------------------------------------------
//                      Class AliRsnParticle
//
//           A simple object which describes a reconstructed track
//           with some references to its related generated particle
//           and some facilities which could help in composing its
//           4-momentum, for resonance study.
//
// author: A. Pulvirenti             (email: alberto.pulvirenti@ct.infn.it)
//-------------------------------------------------------------------------

#ifndef ALIRSNPARTICLE_H
#define ALIRSNPARTICLE_H

class TParticle;

class AliRsnParticle : public TObject
{
public:

    AliRsnParticle();
    virtual ~AliRsnParticle() { }
    
    void     Adopt(TParticle *part);
    
    Int_t    PDG() const {return fPDG;}
    Int_t    Mother() const {return fMother;}
    Short_t  MotherPDG() const {return fMotherPDG;}
    void     SetPDG(Int_t pdg) {fPDG = pdg;}
    void     SetMother(Int_t mlabel) {fMother = mlabel;}
    void     SetMotherPDG(Int_t pdg) {fMotherPDG = (Short_t)pdg;}

private:
	
    Int_t    fPDG;          // PDG code
	Int_t    fMother;       // GEANT label of mother particle
	Short_t  fMotherPDG;    // PDG code of mother particle
	
	ClassDef(AliRsnParticle,1)
};

#endif
