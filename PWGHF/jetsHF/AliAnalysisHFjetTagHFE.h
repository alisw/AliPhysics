#ifndef ALIANALYSISTASKHFJETTAGHFE_H
#define ALIANALYSISTASKHFJETTAGHFE_H

// $Id$

class TH1;
class TH2;
class TH3;
class THnSparse;
class TClonesArray;
class TArrayF;
class AliAODEvent; // sample
class AliJetContainer;
class AliParticleContainer;
class AliClusterContainer;
class AliPIDResponse; 
class AliAODMCHeader;
class AliAODMCParticle; // sample
class AliMultSelection;
class TRandom;

#include "TObject.h"
#include "TObjArray.h"
#include "TClonesArray.h"
//#include "AliAODMCParticle"
#include "AliAnalysisTaskEmcalJet.h"

class AliAnalysisHFjetTagHFE : public AliAnalysisTaskEmcalJet {
 public:

  AliAnalysisHFjetTagHFE();
  AliAnalysisHFjetTagHFE(const char *name);
  virtual ~AliAnalysisHFjetTagHFE();

  void                        UserCreateOutputObjects();
  void                        Terminate(Option_t *option);

  void SetCentralityMimHFEjet(Int_t centMim) {fcentMim = centMim;};
  void SetCentralityMaxHFEjet(Int_t centMax) {fcentMax = centMax;};
  void SetDebugHFEjet(Bool_t dbHFEj) {idbHFEj = dbHFEj;};
  void SetHybridTrack(Bool_t Hybrid){iHybrid = Hybrid;};
  void SetMinSig(Double_t mimSig){fmimSig = mimSig;};
  void SetMinEop(Double_t mimEop){fmimEop = mimEop;};
  void SetMCdata(Bool_t mcData) {fmcData = mcData;};
  void SetInvMassCut0(Double_t InvmassCut) {fInvmassCut = InvmassCut;};
  void SetInvMassCut1(Double_t ptAssocut) {fptAssocut = ptAssocut;};

 protected:
  void                        ExecOnce();
  Bool_t                      FillHistograms()   ;
  Bool_t                      Run()              ;
  void                        CheckClusTrackMatching();


  AliVEvent   *fVevent;  //!event object
  AliMultSelection *fMultSelection;
  TClonesArray  *ftrack;
  TClonesArray  *fCaloClusters;
  AliAODMCHeader *fMCheader; 
  AliPIDResponse *fpidResponse; //!pid response

    Float_t fcentMim; // mim. centrality
    Float_t fcentMax; // max. centrality
    Bool_t idbHFEj;
    Bool_t iHybrid;
    Double_t fmimSig; // max. centrality
    Double_t fmimEop; // max. centrality
    Double_t fInvmassCut;  
    Double_t fptAssocut;  
    Bool_t fmcData;
    Int_t NembMCpi0;
    Int_t NembMCeta;
    Int_t NpureMCproc;

  // General histograms
  TH1                       **fHistTracksPt;            //!Track pt spectrum
  TH1                       **fHistClustersPt;          //!Cluster pt spectrum
  TH1                       **fHistLeadingJetPt;        //!Leading jet pt spectrum
  TH2                       **fHistJetsPhiEta;          //!Phi-Eta distribution of jets
  TH2                       **fHistJetsPtArea;          //!Jet pt vs. area
  TH2                       **fHistJetsPtLeadHad;       //!Jet pt vs. leading hadron
  TH2                       **fHistJetsCorrPtArea;      //!Jet pt - bkg vs. area
  TH3                        *fHistPtDEtaDPhiTrackClus; //!track pt, delta eta, delta phi to matched cluster
  TH3                        *fHistPtDEtaDPhiClusTrack; //!cluster pt, delta eta, delta phi to matched track

  TH1F                        *fHistClustDx; //!
  TH1F                        *fHistClustDz; //!

  TH1F                        *fHistMultCent; //!
  TH2F                        *fHistZcorr; //!
  TH1F                        *fHistCent; //!
  TH2F                        *fHistTPCnSigma;
  TH2F                        *fHistEopNsig;
  TH2F                        *fHistEop;
  TH2F                        *fHistEopHad;
  TH1F                        *fHistJetOrg;
  TH2F                        *fHistJetOrgArea;
  TH1F                        *fHistJetBG;
  TH1F                        *fHistJetSub;
  TH1F                        *fHisteJetOrg;
  TH1F                        *fHisteJetBG;
  TH1F                        *fHisteJetSub;
  TH1F                        *fHistIncEle;
  TH1F                        *fHistIncEleInJet0;
  TH1F                        *fHistIncEleInJet1;
  TH1F                        *fHistHfEleMC;
  TH1F                        *fHistHfEleMCreco;
  TH1F                        *fHistPhoEleMC;
  TH1F                        *fHistPhoEleMCpi0;
  TH1F                        *fHistPhoEleMCeta;
  TH1F                        *fHistPhoEleMCreco;
  TH1F                        *fHistPhoEleMCrecopi0;
  TH1F                        *fHistPhoEleMCrecoeta;
  TH1F                        *fHistMCorgPi0;
  TH1F                        *fHistMCorgEta;
  TH2F                        *fHistIncjet;
  TH2F                        *fHistIncjetFrac;
  TH2F                        *fHistIncjetOrg;
  TH2F                        *fHistIncjetBG;
  TH2F                        *fHistHFjet;
  TH1F                        *fHistHFdijet;
  TH2F                        *fHistULSjet;
  TH2F                        *fHistLSjet;
  TH2F                        *fHistHFjetOrder;
  TH2F                        *fHistDiJetPhi; 
  TH2F                        *fHistDiJetMomBalance; 
  TH2F                        *fHistDiJetMomBalance_All; 
  TH2F                        *fHistDiJetPhi_MC; 
  TH2F                        *fHistDiJetMomBalance_MC; 
  TH2F                        *fInvmassULS;
  TH2F                        *fInvmassLS;
  TH2F                        *fInvmassHFuls;
  TH2F                        *fInvmassHFls;
  TH1F                        *fLxy_uls;
  TH1F                        *fLxy_ls;
  THnSparse                   *HFjetCorr0;
  THnSparse                   *HFjetCorr1;
  THnSparse                   *HFjetParticle;
  TH1F                        *fQAHistJetPhi;
  TH1F                        *fQAHistTrPhiJet;
  TH1F                        *fQAHistTrPhi;
  TH1F                        *fQAHistNits;
  TH2F                        *fQAHistEleDCAxy;
  TH2F                        *fQAHistEleDCAz;
  TH1F                        *fHistClustE;
  TH1F                        *fHistClustEtime;
  TH2F                        *fEMCClsEtaPhi;
  TH1F                        *fHistBGfrac;
  TH1F                        *fHistBGfracHFEev;
  TF1                         *fPi0Weight;
  TF1                         *fEtaWeight;
  TRandom                     *generator;

  AliJetContainer            *fJetsCont;                   //!Jets
  AliJetContainer            *fJetsContPart;                   //!Jets particle
  AliParticleContainer       *fTracksCont;                 //!Tracks
  AliClusterContainer        *fCaloClustersCont;           //!Clusters  
  Bool_t tagHFjet(AliEmcalJet* jet, double *epT, int MCpid, double &maxpT_e);
  //void SelectPhotonicElectron(Int_t itrack, AliVTrack *track, Bool_t &fFlagPhotonicElec);
  void SelectPhotonicElectron(Int_t itrack, AliVTrack *track, Bool_t &fFlagPhotonicElec, Bool_t &fFlagConvinatElec);
  Double_t CalRandomCone(Double_t HFjetPhi[], Double_t HFjetEta[], Double_t HFjetArea);
  Bool_t isHeavyFlavour(int Mompdg);
  Bool_t isPhotonic(int Mompdg);
  //void MakeParticleLevelJet(THnSparse *pJet);
  void MakeParticleLevelJet();
  //void SetCentralityMim(Int_t centMim) {fcentMim = centMim;};
  //void SetCentralityMax(Int_t centMax) {fcentMax = centMax;};
  void FindMother(AliAODMCParticle* part, int &label, int &pid, double &ptmom);

 private:


  //TClonesArray  *ftrack;

  AliAODEvent 		*fAOD;			
  TClonesArray 		*fMCarray;
  AliAODMCParticle 	*fMCparticle;
  AliAODMCParticle 	*fMCparticleMother;

  //Bool_t fmcData;

  AliAnalysisHFjetTagHFE(const AliAnalysisHFjetTagHFE&);            // not implemented
  AliAnalysisHFjetTagHFE &operator=(const AliAnalysisHFjetTagHFE&); // not implemented

  ClassDef(AliAnalysisHFjetTagHFE, 6) // jet sample analysis task
};
#endif
