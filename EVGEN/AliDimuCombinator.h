#ifndef _AliDimuCombinator_H
#define _AliDimuCombinator_H
#include "TParticle.h"
#include <TBrowser.h>
#include <TList.h>
#include <TTree.h>
#include <TROOT.h>


class AliDimuCombinator:
    public TObject 
{
public:
    AliDimuCombinator(TClonesArray* Partarray){
	fPartArray=Partarray;
	fNParticle=fPartArray->GetEntriesFast();
	
	fimuon1 =0;
	fimuon2 =0;
	fmuon1  =0;
	fmuon2  =0;
	fimin1  = 0;
	fimin2  = 0;
	fimax1  = fNParticle;
	fimax2  = fNParticle;
	fPtMin  =0;
	fEtaMin =-10;
	fEtaMax =-10;
	fRate1=1.;
	fRate2=1.;
    }
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
    void SetEtaCut(Float_t etamin, Float_t etamax){fEtaMin=etamin; fEtaMax=etamax;}      Bool_t Selected(TParticle* part);
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
    Float_t Decay_Prob(TParticle* part);
    
 private:
    void FirstPartner();
    void NextPartner();
    void FirstPartnerSelected();
    void NextPartnerSelected();
    Int_t Origin(TParticle* part);
    TParticle* Parent(TParticle* part);
    TParticle* Partner();
    Int_t Type(TParticle *part){return part->GetPdgCode();}
private:
    TClonesArray *fPartArray;
    Int_t fNParticle;
    Int_t fimuon1;
    Int_t fimuon2;
    Int_t fimin1;
    Int_t fimin2;
    Int_t fimax1;
    Int_t fimax2;
    Float_t fRate1;
    Float_t fRate2;
    TParticle *fmuon1;
    TParticle *fmuon2;
    Float_t fPtMin;
    Float_t fEtaMin;
    Float_t fEtaMax;
  ClassDef(AliDimuCombinator,1) // Dimuon Combinator
};
#endif




