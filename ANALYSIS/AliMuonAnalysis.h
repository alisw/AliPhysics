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

    TFile *fHistoFile;
    TH1F *fHPtMuon;
    TH1F *fHPtMuonPlus;
    TH1F *fHPtMuonMinus;
    TH1F *fHPMuon;
    TH1F *fHInvMassAll;
    TH1F *fHRapMuon;
    TH1F *fHRapResonance;
    TH1F *fHPtResonance;
    TH2F *fHInvMassAll_vs_Pt;

    AliAODParticleCut* fPartCut;//Particle Cut
    ClassDef(AliMuonAnalysis,1)
};

#endif
