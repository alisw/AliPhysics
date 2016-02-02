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

  TH1                        *fHistClustDx; //!
  TH1                        *fHistClustDz; //!

  TH2                        *fHistTPCnSigma;
  TH2                        *fHistEop;
  TH1                        *fHistJetOrg;
  TH1                        *fHistJetBG;
  TH1                        *fHistJetSub;
  TH1                        *fHisteJetOrg;
  TH1                        *fHisteJetBG;
  TH1                        *fHisteJetSub;
  TH1                        *fHistIncEle;
  TH1                        *fHistHfEleMC;
  TH1                        *fHistHfEleMCreco;
  TH1                        *fHistPhoEleMC;
  TH1                        *fHistPhoEleMCreco;
  TH2                        *fHistIncjet;
  TH2                        *fHistIncjetFrac;
  TH2                        *fHistIncjetOrg;
  TH2                        *fHistIncjetBG;
  TH2                        *fHistHFjet;
  TH2                        *fHistULSjet;
  TH2                        *fHistLSjet;
  TH2                        *fInvmassULS;
  TH2                        *fInvmassLS;
  THnSparse                  *HFjetCorr0;
  THnSparse                  *HFjetCorr1;
  THnSparse                  *HFjetParticle;
  TH1                        *fQAHistJetPhi;
  TH1                        *fQAHistTrPhiJet;
  TH1                        *fQAHistTrPhi;


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
  void SetCentralityMim(Int_t centMim) {fcentMim = centMim;};
  void SetCentralityMax(Int_t centMax) {fcentMax = centMax;};

 private:

  //TClonesArray  *ftrack;

  AliAODEvent 		*fAOD;			
  TClonesArray 		*fMCarray;
  AliAODMCParticle 	*fMCparticle;
  AliAODMCParticle 	*fMCparticleMother;

  Bool_t fmcData;

  AliAnalysisHFjetTagHFE(const AliAnalysisHFjetTagHFE&);            // not implemented
  AliAnalysisHFjetTagHFE &operator=(const AliAnalysisHFjetTagHFE&); // not implemented

  ClassDef(AliAnalysisHFjetTagHFE, 5) // jet sample analysis task
};
#endif
