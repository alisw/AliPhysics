#ifndef ALIMUONANALYSIS_H
#define ALIMUONANALYSIS_H
//________________________________
///////////////////////////////////////////////////////////
//
// class AliMuonAnalysis
//
// Flow Analysis
//
//
// S.Radomski@gsi.de
// Piotr.Skowronski@cern.ch
//
///////////////////////////////////////////////////////////

#include "AliAnalysis.h"

class AliESD;
class AliAOD;
class AliStack;
class AliAODParticleCut;
class TFile;
class TH1F;
class TH2F;

class AliMuonAnalysis: public AliAnalysis
{ 
  public: 
     AliMuonAnalysis();
     virtual ~AliMuonAnalysis();

    Int_t Init();
    Int_t ProcessEvent(AliAOD* aodrec, AliAOD* aodsim);
    Int_t Finish();
   
    void SetParticleCut(AliAODParticleCut* pcut){fPartCut = pcut;}

    void GetInvMass(AliAOD* aod);

  protected:
    
  private:

    TFile *fHistoFile;         // histogramm file pointer
    TH1F *fHPtMuon;            // Muon Pt distribution
    TH1F *fHPtMuonPlus;        // Muon Plus Pt distribution
    TH1F *fHPtMuonMinus;       // Muon Minus Pt distribution
    TH1F *fHPMuon;             // Muon momentum distribution
    TH1F *fHInvMassAll;        // Invariant mass distribution
    TH1F *fHRapMuon;           // Muon rapidity distribution
    TH1F *fHRapResonance;      // Muon rapidity distribution around resonance
    TH1F *fHPtResonance;       // Muon Pt distribution around resonance
    TH2F *fHInvMassAllvsPt;    // Invariant mass vs Pt distribution

    AliAODParticleCut* fPartCut;//Particle Cut
    ClassDef(AliMuonAnalysis,1)
};

#endif
