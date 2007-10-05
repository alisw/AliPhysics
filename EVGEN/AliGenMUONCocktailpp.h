#ifndef ALIGENMUONCOCKTAILPP_H
#define ALIGENMUONCOCKTAILPP_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliGenCocktail.h"
#include "AliDecayer.h"

class AliGenCocktailEntry;

class AliGenMUONCocktailpp : public AliGenCocktail
{
 public:

    AliGenMUONCocktailpp();
    virtual ~AliGenMUONCocktailpp();    
    virtual void Init();
    virtual void Generate();    
    Int_t   GetNSucceded()         const {return fNSucceded;}    
    Int_t   GetNGenerated()        const {return fNGenerated;}
    Int_t   GetMuonMultiplicity()  const {return fMuonMultiplicity;}
    Float_t GetMuonPtCut()         const {return fMuonPtCut;}
    Float_t GetMuonPCut()          const {return fMuonPCut;}    
    Float_t GetMuonThetaMin()      const {return fMuonThetaMinCut;}
    Float_t GetMuonThetaMax()      const {return fMuonThetaMaxCut;}	    
    Float_t GetMuonOriginCut()     const {return fMuonOriginCut;}	    
    Float_t GetDecayModeResonance()const {return fDecayModeResonance;}
    Float_t GetDecayModePythia()   const {return fDecayModePythia;}
    
    void    SetMuonMultiplicity(Int_t MuonMultiplicity) { fMuonMultiplicity = MuonMultiplicity;}
    void    SetMuonPtCut(Float_t PtCut) { fMuonPtCut = PtCut;}
    void    SetMuonPCut(Float_t PCut) { fMuonPCut = PCut;}    
    void    SetMuonOriginCut(Float_t originCut) { fMuonOriginCut = originCut;}
    void    SetMuonThetaRange(Float_t ThetaMin, Float_t ThetaMax){
	fMuonThetaMinCut=ThetaMin;
	fMuonThetaMaxCut=ThetaMax; }    
    void    SetDecayer(AliDecayer* decayer){fDecayer = decayer;}
    void    SetDecayModeResonance(Decay_t decay){ fDecayModeResonance = decay;}
    void    SetDecayModePythia(Decay_t decay){ fDecayModePythia = decay;}
 protected:

    //
 private:
    AliGenMUONCocktailpp(const AliGenMUONCocktailpp &cocktail); 
    AliGenMUONCocktailpp & operator=(const AliGenMUONCocktailpp &cocktail); 
    AliDecayer* fDecayer;
    Decay_t fDecayModeResonance; //decay mode in which resonances are forced to decay, default: kAll
    Decay_t fDecayModePythia; //decay mode in which particles in Pythia are forced to decay, default: kAll
    Float_t fTotalRate;// Total rate of the full cocktail processes
    Int_t   fMuonMultiplicity; // Muon multiplicity for the primordial trigger
    Float_t fMuonPtCut;// Transverse momentum cut for muons
    Float_t fMuonPCut;// Momentum cut for muons    
    Float_t fMuonThetaMinCut;// Minimum theta cut for muons
    Float_t fMuonThetaMaxCut; // Maximum theta cut for muons
    Float_t fMuonOriginCut; //use only muons whose "part->Vz()" value is larger than fMuonOrigin
    Int_t   fNSucceded;// Number of Succes in the (di)-muon generation in the acceptance
    Int_t   fNGenerated;// Number of generated cocktails
    

    ClassDef(AliGenMUONCocktailpp,2)  //  cocktail for physics in the Alice
};

#endif





