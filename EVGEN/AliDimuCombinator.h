#ifndef ALIDIMUCOMBINATOR_H
#define ALIDIMUCOMBINATOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
#include <TParticle.h>
#include <TObject.h>

class TClonesArray;
class TParticle;


class AliDimuCombinator:
    public TObject 
{
public:
    AliDimuCombinator(TClonesArray* Partarray);
    AliDimuCombinator(const AliDimuCombinator &combinator);
//    
//  Iterators
//  Single muons
    TParticle* FirstMuon();
    TParticle* NextMuon();
//  Single muons selected
    TParticle* FirstMuonSelected();
    TParticle* NextMuonSelected();
//  Dimuons    
    void FirstMuonPair(TParticle* & muon1, TParticle* & muon2);
    void NextMuonPair(TParticle* & muon1, TParticle* & muon2);
//  Dimuons selected    
    void FirstMuonPairSelected(TParticle* & muon1, TParticle* & muon2);
    void NextMuonPairSelected(TParticle* & muon1, TParticle* & muon2);
//  Loop over all prticles    
    void ResetRange();
//  Set two ranges for dimuon loop    
    void SetFirstRange (Int_t from, Int_t to);
    void SetSecondRange(Int_t from, Int_t to);    
//  Cuts
    void SetPtMin(Float_t ptmin) {fPtMin=ptmin;}
    void SetEtaCut(Float_t etamin, Float_t etamax)
	{fEtaMin=etamin; fEtaMax=etamax;}
    Bool_t Selected(TParticle* part);
    Bool_t Selected(TParticle* part1, TParticle* part2);
// Kinematics
    Float_t Mass(TParticle* part1, TParticle* part);
    Float_t PT(TParticle* part1, TParticle* part);
    Float_t Pz(TParticle* part1, TParticle* part);
    Float_t Y(TParticle* part1, TParticle* part);
// Response
    void SmearGauss(Float_t width, Float_t & value);
// Weight
    Bool_t  Correlated(TParticle* part1, TParticle* part2);
    void    SetRate(Float_t rate){fRate1=rate;}
    void    SetRate(Float_t rate1, Float_t rate2 ){fRate1=rate1; fRate2=rate2;}
    Float_t Weight(TParticle* part);
    Float_t Weight(TParticle* part1, TParticle* part);
    Float_t DecayProbability(TParticle* part);
    AliDimuCombinator & operator=(const AliDimuCombinator & rhs);
    
 private:
    void FirstPartner();
    void NextPartner();
    void FirstPartnerSelected();
    void NextPartnerSelected();
    Int_t Origin(TParticle* part);
    TParticle* Parent(TParticle* part);
    TParticle* Partner();
    Int_t Type(TParticle *part) {return part->GetPdgCode();}
private:
    TClonesArray *fPartArray;      // MC Event
    Int_t fNParticle;              // Number of particles
    Int_t fimuon1;                 // Index of first muon
    Int_t fimuon2;                 // Index of second muon
    Int_t fimin1;                  // Lowest index for first   muon  
    Int_t fimin2;                  // Lowest index for second  muon 
    Int_t fimax1;                  // Highest index for first  muon  
    Int_t fimax2;                  // Highest index for second muon 
    Float_t fRate1;                // weight factor  
    Float_t fRate2;                // weight factor
    TParticle *fmuon1;             // First muon
    TParticle *fmuon2;             // Second muon
    Float_t fPtMin;                // pT-cut 
    Float_t fEtaMin;               // Minimum pseudorapidity cut
    Float_t fEtaMax;               // Maximum pseudorapidity cut
    
    ClassDef(AliDimuCombinator,1)  // Tools for dimuon combinatoric studies
};
#endif





