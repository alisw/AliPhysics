#ifndef ALIANALYSISTASKSEMUONSHF_H
#define ALIANALYSISTASKSEMUONSHF_H

#include <TString.h>
#include <TList.h>
#include <TClonesArray.h>

#include "AliMuonsHFHeader.h"
#include "AliAODMuonTrack.h"
#include "AliAODMuonPair.h"
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskSEMuonsHF : public AliAnalysisTaskSE {
 public:

  AliAnalysisTaskSEMuonsHF();
  AliAnalysisTaskSEMuonsHF(const char *name);
  virtual ~AliAnalysisTaskSEMuonsHF();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *opt);
  virtual void Terminate(Option_t *opt);

  void SetAnaMode(Int_t mode) { fAnaMode = (mode<3 ? mode: 0); }
  void SetIsOutputTree(Bool_t ist) { fIsOutputTree = ist; }
  void SetIsUseMC(Bool_t isMC) { fIsUseMC = isMC; }

  void SetSingleMuonCuts(Double_t cuts[10]) {
    for (Int_t i=0; i<10; i++) fSingleMuonCuts[i]=cuts[i];
  }

 private:

  void CreateOutoutHistosAtEvnetLevel();
  void CreateOutputHistosSingleMuon();
  void CreateOutputHistosDimuon();
  void CreateOutputHistosSingleMuonMC(Int_t *nbins, Double_t *xlow, Double_t *xup,
                                      TString *name, TString *axis, TString *unit);
  void CreateOutputHistosDimuonMC(Int_t *nbins, Double_t *xlow, Double_t *xup,
                                  TString *name, TString *axis, TString *unit,
                                  TString *dimuName, TString *dimuTitle);

  void FillDistributionsAtEventLeavel();
  void FillDistributionsSingleMuon(AliAODMuonTrack *track, Int_t src=-1);
  void FillDistributionsDimuon(AliAODMuonPair *pair, Int_t src=-1);
  void FillDistributionsSingleMuonMC(Double_t *disMu, Int_t src);
  void FillDistributionsDimuonMC(Double_t *disDimu, Int_t dimuK, Int_t src);

  enum {  // histos at event level
    kHistVx,    // x position of vtx
    kHistVy,    // y position of vtx
    kHistVz,    // z position of vtx
    kHistMult,  // multiplicity of single muon track
    kNHistEv    // number of histos at event level
  };
  enum {  // distributions of single mu
    kHistP,    // histo of single mu p
    kHistPt,   // histo of single mu pt
    kHistEta,  // histo of single mu eta
    kHistDca,  // histo of single mu DCA
    kHistTrg,  // histo of single mu pass trigger
    kNHistMu   // number of histos of single mu
  };
  enum {  // distribution of dimuon
    kInvM,      // invariance mass if dimuon
    kPtPair,    // pt of dimu pair
    kNHistDimu  // number of histos of dimuon
  };

  enum {  // MC sources of single muon
    kBeautyMu,      // mu<-b
    kCharmMu,       // mu<-c
    kPrimaryMu,     // primary muon
    kSecondaryMu,   // mu<-secondary paritcle
    kNotMu,         // not muon
    kNSingleMuSrcs  // number of single muon sources
  };
  enum {  // MC soureces of dimuon
    kBBdiff,     // dimuon<-BB_diff
    kBchain,     // dimuon<-B-chain
    kDDdiff,     // dimuon<-DD_diff
    kDchain,     // dimuon<-D-chain
    kResonance,  // dimuon<-resonances
    kUncorr,     // uncorr dimuon
    kNDimuSrcs   // number of dimuon sources
  };
  enum {  // different kinds of correlated dimu
    kMuNMuP,  // mu-mu+
    kMuNMuN,  // mu-mu-
    kMuPMuP,  // mu+mu+
    kNMuMus   // number of corr dimu kinds
  };

  Int_t fAnaMode;  // = 0, ana both single muon and dimuon
                   // = 1, ana single muon
                   // = 2, ana dimuon
  Bool_t fIsOutputTree;  // flag used to switch on/off tree output
  Bool_t fIsUseMC;       // flag used to switch on/off MC ana

  Double_t fSingleMuonCuts[10];  // 0, max of 3-momentum
                                 // 1, min of 3-momentum
                                 // 2, pt_Min
                                 // 3, pt_Max
                                 // 4, eta_Min
                                 // 5, eta_Max
                                 // 6, dca_Min
                                 // 7, dca_Max
                                 // 8, about trigger
                                 // 9, about trigger

  AliMuonsHFHeader *fHeader;   // output clones array for info at ev level
  TClonesArray *fMuTrkClArr;   // output clones array for single mu
  TClonesArray *fMuPairClArr;  // output clones array for dimu

  TList *fListHisAtEvLevel;   // output list of histos at event level
  TList *fListHisSingleMuon;  // output list of histos for single mu
  TList *fListHisDimuon;      // output list of histos for dimuon

  ClassDef(AliAnalysisTaskSEMuonsHF, 5);
};

#endif
