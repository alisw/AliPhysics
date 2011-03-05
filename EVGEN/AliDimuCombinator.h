#ifndef ALIDIMUCOMBINATOR_H
#define ALIDIMUCOMBINATOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


//  Class for dimuon analysis and fast dimuon simulation.
//  It uses the AliRun particle tree.
//  Comments and suggestions to andreas.morsch@cern.ch


#include <TObject.h>

class TClonesArray;
class TParticle;


class AliDimuCombinator:
    public TObject 
{
public:
    AliDimuCombinator();
    void  Copy(TObject &combi) const;
//    
//  Iterators
//  Access to particle stack
    TParticle* Particle(Int_t i) const;
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
    Bool_t Selected(const TParticle* part) const;
    Bool_t Selected(const TParticle* part1, const TParticle* part2) const;
// Kinematics
    Float_t Mass(const TParticle* part1, const TParticle* part) const;
    Float_t PT(const TParticle* part1, const TParticle* part) const;
    Float_t Pz(const TParticle* part1, const TParticle* part) const;
    Float_t Y(const TParticle* part1, const TParticle* part) const;
// Response
    void SmearGauss(Float_t width, Float_t & value) const;
// Weight
    Bool_t  Correlated(const TParticle* part1, const TParticle* part2) const;
    void    SetRate(Float_t rate) {fRate1=rate;}
    void    SetRate(Float_t rate1, Float_t rate2 ) {fRate1=rate1; fRate2=rate2;}
    Float_t Weight(const TParticle* part) const;
    Float_t Weight(const TParticle* part1, const TParticle* part) const;
    Float_t DecayProbability(const TParticle* part) const;
    
 private:
    void FirstPartner();
    void NextPartner();
    void FirstPartnerSelected();
    void NextPartnerSelected();
    Int_t Origin(const TParticle* part) const;
    TParticle* Parent(const TParticle* part) const;
    TParticle* Partner() const;
    Int_t Type(const TParticle *part) const;
    AliDimuCombinator(const AliDimuCombinator &combinator);
    AliDimuCombinator & operator=(const AliDimuCombinator & rhs);

 private:
    Int_t fNParticle;              // Number of particles
    Int_t fImuon1;                 // Index of first muon
    Int_t fImuon2;                 // Index of second muon
    Int_t fImin1;                  // Lowest index for first   muon  
    Int_t fImin2;                  // Lowest index for second  muon 
    Int_t fImax1;                  // Highest index for first  muon  
    Int_t fImax2;                  // Highest index for second muon 
    Float_t fRate1;                // weight factor  
    Float_t fRate2;                // weight factor
    TParticle *fMuon1;             // First muon
    TParticle *fMuon2;             // Second muon
    Float_t fPtMin;                // pT-cut 
    Float_t fEtaMin;               // Minimum pseudorapidity cut
    Float_t fEtaMax;               // Maximum pseudorapidity cut
    
    ClassDef(AliDimuCombinator,1)  // Tools for dimuon combinatoric studies
};
#endif








