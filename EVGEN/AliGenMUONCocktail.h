#ifndef ALIGENMUONCOCKTAIL_H
#define ALIGENMUONCOCKTAIL_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

// Classe to create the MUON coktail for physics in the Alice muon spectrometer
// Gines Martinez, jan 2004, Nantes  martinez@in2p3.fr

#include "AliGenCocktail.h"

class AliGenCocktailEntry;


class AliGenMUONCocktail : public AliGenCocktail
{
 public:
    AliGenMUONCocktail();
    AliGenMUONCocktail(const AliGenMUONCocktail &cocktail); 
    virtual ~AliGenMUONCocktail();
    virtual void Init();
    virtual void Generate();
    Int_t   GetMuonMultiplicity() {return fMuonMultiplicity;}
    Int_t   GetNSucceded()  {return fNSucceded;}
    Int_t   GetNGenerated() {return fNGenerated;}
    Float_t GetNumberOfCollisions()   {return fNumberOfCollisions;}
    Float_t GetNumberOfParticipants() {return fNumberOfParticipants;}
    Float_t GetMuonPtCut()  { return fMuonPtCut;}

    void    SetMuonMultiplicity(Int_t MuonMultiplicity) { fMuonMultiplicity= MuonMultiplicity;}
    void    SetNumberOfCollisions(Float_t NumberOfCollisions) { fNumberOfCollisions= NumberOfCollisions;}
    void    SetNumberOfParticipants(Float_t NumberOfParticipants) { fNumberOfParticipants= NumberOfParticipants;}
    void    SetMuonPtCut(Float_t PtCut) { fMuonPtCut = PtCut;}
    void    SetMuonThetaCut(Float_t ThetaMin, Float_t ThetaMax) 
      { fMuonThetaMinCut=ThetaMin; 
        fMuonThetaMaxCut=ThetaMax; } 

 protected:
 
    //
 private:
    Float_t fTotalRate;  // Total rate of the full cocktail processes
    Int_t   fMuonMultiplicity; // Muon multiplicity for the primordial trigger
    Float_t fMuonPtCut;       // Transverse momentum cut for muons
    Float_t fMuonThetaMinCut; // Minimum theta cut for muons
    Float_t fMuonThetaMaxCut; // Maximum theta cut for muons
    Int_t   fNSucceded;  //  Number of Succes in the dimuon pair generation in the acceptance
    Int_t   fNGenerated; // Number of generated cocktails
    Float_t fNumberOfCollisions; // Average Number of collisions in the centrality class 
    Float_t fNumberOfParticipants; // Average Number of participants in the centrality class 
    
    ClassDef(AliGenMUONCocktail,1) //  MUON cocktail for physics in the Alice muon spectrometer
};

#endif





