#ifndef ALIGENMUONCOCKTAIL_H
#define ALIGENMUONCOCKTAIL_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id$ */
//
// Classe to create the MUON coktail for physics in the Alice muon spectrometer
// The followoing muons sources are included in this cocktail: 
//     jpsi, upsilon, non-correlated open and beauty, and muons from pion and kaons.
// The free parameeters are :
//      pp reaction cross-section
//      production cross-sections in pp collisions and 
//      branching ratios in the muon channel
// Hard probes are supposed to scale with Ncoll and hadronic production with (0.8Ncoll+0.2*Npart)
// There is a primordial trigger wiche requires :
//      a minimum muon multiplicity above a pT cut in a theta acceptance cone
//
// Gines Martinez, jan 2004, Nantes  martinez@in2p3.fr




#include "AliGenCocktail.h"

class AliFastGlauber;
class AliGenCocktailEntry;


class AliGenMUONCocktail : public AliGenCocktail
{
 public:
    AliGenMUONCocktail();
     virtual ~AliGenMUONCocktail();
    virtual void Init();
    virtual void Generate();
    Int_t   GetMuonMultiplicity()  const {return fMuonMultiplicity;}
    Int_t   GetNSucceded()         const {return fNSucceded;}
    Int_t   GetNGenerated()        const {return fNGenerated;}
    Float_t GetNumberOfCollisions()const {return fNumberOfCollisions;} 
    Float_t GetNumberOfParticipants() const {return fNumberOfParticipants;}
    Float_t GetMuonPtCut()         const { return fMuonPtCut;}

    void    SetMuonMultiplicity(Int_t MuonMultiplicity) { fMuonMultiplicity= MuonMultiplicity;}
    void    SetNumberOfCollisions(Float_t NumberOfCollisions) { fNumberOfCollisions= NumberOfCollisions;}
    void    SetNumberOfParticipants(Float_t NumberOfParticipants) { fNumberOfParticipants= NumberOfParticipants;}
    void    SetImpactParameterRange(Float_t bmin=0., Float_t bmax=5.) { fLowImpactParameter = bmin; fHighImpactParameter=bmax;}
    void    SetMuonPtCut(Float_t PtCut) { fMuonPtCut = PtCut;}
    void    SetMuonThetaCut(Float_t ThetaMin, Float_t ThetaMax) 
      { fMuonThetaMinCut=ThetaMin; 
        fMuonThetaMaxCut=ThetaMax; }
    void    SetHadronicMuons(Bool_t HadronicMuons) { fHadronicMuons = HadronicMuons;}
    void    SetInvMassRange(Float_t MassMin, Float_t MassMax) 
      { fInvMassMinCut=MassMin; 
        fInvMassMaxCut=MassMax;
        fInvMassCut = kTRUE; }
 private:
    AliGenMUONCocktail(const AliGenMUONCocktail &cocktail); 
    AliGenMUONCocktail& operator=(const AliGenMUONCocktail & rhs);

    //
 private:
    AliFastGlauber *  fFastGlauber; //! Fast glauber calculations
    Float_t fTotalRate;             // Total rate of the full cocktail processes
    Int_t   fMuonMultiplicity;      // Muon multiplicity for the primordial trigger
    Float_t fMuonPtCut;             // Transverse momentum cut for muons
    Float_t fMuonThetaMinCut;       // Minimum theta cut for muons
    Float_t fMuonThetaMaxCut;       // Maximum theta cut for muons
    Int_t   fNSucceded;             // Number of Succes in the dimuon pair generation in the acceptance
    Int_t   fNGenerated;            // Number of generated cocktails
    Float_t fLowImpactParameter;    // Lowest simulated impact parameter
    Float_t fHighImpactParameter;   // Highest impact parameter
    Float_t fAverageImpactParameter;// AVergae Impact parameter in the impact parameter range
    Float_t fNumberOfCollisions;    // Average number of collisions in the impact parameter range
    Float_t fNumberOfParticipants;  // Average number of participants in the impact parameter range
    Bool_t  fHadronicMuons;         // If kTRUE hadronic muons are included in the cocktail. Default is true.
    Bool_t  fInvMassCut;            // If kTRUE cut on the Invariant mass is required. Default is false
    Float_t fInvMassMinCut;	    // Minimum invariant mass cut
    Float_t fInvMassMaxCut;	    // Maximum invariant mass cut
   
    ClassDef(AliGenMUONCocktail,1)  //  MUON cocktail for physics in the Alice muon spectrometer
};

#endif





