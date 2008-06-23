/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               *
 **************************************************************************/

//-------------------------------------------------------------------------
//                      Class AliRsnAnalysis
//             Reconstruction and analysis of K* Rsn
// ........................................
// ........................................
// ........................................
// ........................................
//
// author: A. Pulvirenti             (email: alberto.pulvirenti@ct.infn.it)
//-------------------------------------------------------------------------

#ifndef ALIRSNPAIRDEF_H
#define ALIRSNPAIRDEF_H

class AliRsnDaughter;
#include "AliRsnPID.h"

class AliRsnPairDef : public TObject
{
public:

    AliRsnPairDef();
    AliRsnPairDef(Char_t sign1, AliRsnPID::EType type1, Char_t sign2, AliRsnPID::EType type2,
                  Int_t nbins, Double_t min, Double_t max);
    AliRsnPairDef(const AliRsnPairDef &copy);
    const AliRsnPairDef& operator=(const AliRsnPairDef &copy);
    virtual ~AliRsnPairDef() { }

    // getters
    Char_t           GetCharge(Int_t i) const {if (i>=0&&i<2) return fCharge[i]; else return 0;}
    AliRsnPID::EType GetType(Int_t i) const {if (i>=0&&i<2) return fType[i]; else return AliRsnPID::kUnknown;}
    Double_t         GetMass(Int_t i) const {if (i>=0&&i<2) return fMass[i]; else return 0.0;}
    Int_t            GetNBins() const {return fNBins;}
    Double_t         GetMin() const {return fMin;}
    Double_t         GetMax() const {return fMax;}
    Int_t            GetMotherPDG() const {return fMotherPDG;}

    // setters
    Bool_t SetPairElement(Int_t i, Char_t charge, AliRsnPID::EType pid);
    Bool_t SetPair(Char_t charge1, AliRsnPID::EType pid1, Char_t charge2, AliRsnPID::EType pid2);
    void   SetMotherPDG(Int_t pdg) {fMotherPDG = pdg;}
    void   SetNBins(Int_t value) {fNBins = value;}
    void   SetMin(Double_t value) {fMin = value; CheckEdges();}
    void   SetMax(Double_t value) {fMax = value; CheckEdges();}
    void   SetBins(Int_t nbins, Double_t min, Double_t max) {fNBins=nbins; fMin=min; fMax=max; CheckEdges();}

    // pair information methods
    Bool_t IsLikeSign() {return (fCharge[0] == fCharge[1]);}
    Bool_t HasEqualTypes() {return (fType[0] == fType[1]);}

    // working routines
    void     CheckEdges();
    Double_t ComputeWeight(AliRsnDaughter *d0, AliRsnDaughter *d1);

private:

    // pair parameters
    Int_t             fMotherPDG;  // PDG code of true mother (if known)
    Double_t          fMass[2];    // mass of particles
    Char_t            fCharge[2];  // charge of particles
    AliRsnPID::EType  fType[2];    // PID of particles

    // histogram parameters
    Int_t             fNBins;      // number of histogram bins
    Double_t          fMin;        // lower edge
    Double_t          fMax;        // upper edge

    // ROOT dictionary
    ClassDef(AliRsnPairDef, 1)
};

#endif
