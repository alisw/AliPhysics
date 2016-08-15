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
  void SetDebugHFEjet(Bool_t dbHFEj) {fdbHFEj = dbHFEj;};

 protected:
  void                        ExecOnce();
  Bool_t                      FillHistograms()   ;
  Bool_t                      Run()              ;
  void                        CheckClusTrackMatching();


  AliVEvent   *fVevent;  //!event object
  AliMultSelection *fMultSelection;
  TClonesArray  *ftrack;
  TClonesArray  *fCaloClusters;
  AliPIDResponse *fpidResponse; //!pid response

    Float_t fcentMim; // mim. centrality
    Float_t fcentMax; // max. centrality
    Bool_t fdbHFEj;

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
  TH1F                        *fHistCent; //!
  TH2F                        *fHistTPCnSigma;
  TH2F                        *fHistEop;
  TH1F                        *fHistJetOrg;
  TH1F                        *fHistJetBG;
  TH1F                        *fHistJetSub;
  TH1F                       *fHisteJetOrg;
  TH1F                        *fHisteJetBG;
  TH1F                        *fHisteJetSub;
  TH1F                        *fHistIncEle;
  TH1F                        *fHistHfEleMC;
  TH1F                        *fHistHfEleMCreco;
  TH1F                        *fHistPhoEleMC;
  TH1F                        *fHistPhoEleMCreco;
  TH2F                        *fHistIncjet;
  TH2F                        *fHistIncjetFrac;
  TH2F                        *fHistIncjetOrg;
  TH2F                        *fHistIncjetBG;
  TH2F                        *fHistHFjet;
  TH2F                        *fHistULSjet;
  TH2F                        *fHistLSjet;
  TH2F                        *fInvmassULS;
  TH2F                        *fInvmassLS;
  THnSparse                  *HFjetCorr0;
  THnSparse                  *HFjetCorr1;
  THnSparse                  *HFjetParticle;
  TH1F                        *fQAHistJetPhi;
  TH1F                        *fQAHistTrPhiJet;
  TH1F                        *fQAHistTrPhi;


  AliJetContainer            *fJetsCont;                   //!Jets
  AliParticleContainer       *fTracksCont;                 //!Tracks
  AliClusterContainer        *fCaloClustersCont;           //!Clusters  
  AliJetContainer            *fJetsContPart;                   //!Jets particle
  Bool_t tagHFjet(AliEmcalJet* jet, double *epT, int MCpid, double &maxpT_e);
  //void SelectPhotonicElectron(Int_t itrack, AliVTrack *track, Bool_t &fFlagPhotonicElec);
  void SelectPhotonicElectron(Int_t itrack, AliVTrack *track, Bool_t &fFlagPhotonicElec, Bool_t &fFlagConvinatElec);
  Bool_t isHeavyFlavour(int Mompdg);
  Bool_t isPhotonic(int Mompdg);
  void MakeParticleLevelJet(THnSparse *pJet);
  //void SetCentralityMim(Int_t centMim) {fcentMim = centMim;};
  //void SetCentralityMax(Int_t centMax) {fcentMax = centMax;};

 private:


  //TClonesArray  *ftrack;

  AliAODEvent 		*fAOD;			
  TClonesArray 		*fMCarray;
  AliAODMCParticle 	*fMCparticle;
  AliAODMCParticle 	*fMCparticleMother;

  Bool_t fmcData;

  AliAnalysisHFjetTagHFE(const AliAnalysisHFjetTagHFE&);            // not implemented
  AliAnalysisHFjetTagHFE &operator=(const AliAnalysisHFjetTagHFE&); // not implemented

  ClassDef(AliAnalysisHFjetTagHFE, 6) // jet sample analysis task
};
#endif
