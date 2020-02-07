#ifndef ALIANALYSISTASKCHECKESDTRACKS
#define ALIANALYSISTASKCHECKESDTRACKS

/* Copyright(c) 1998-2012, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// Class AliAnalysisTaskCheckESDTracks
// AliAnalysisTaskSE to extract QA and performance histos for tracks
// 
//
// Authors:
//          F. Prino, C, Zampolli
//          
//*************************************************************************

class TList;
class TTree;
class TH1F;
class TH2F;
class TH3F;
class TString;
class AliESDEvent;

#include "AliESDtrackCuts.h"
#include "AliAnalysisTaskSE.h"
#include "AliPID.h"

class AliAnalysisTaskCheckESDTracks : public AliAnalysisTaskSE {

 public:
  
  AliAnalysisTaskCheckESDTracks();
  virtual ~AliAnalysisTaskCheckESDTracks();

  virtual void   UserExec(Option_t *option);
  virtual void   UserCreateOutputObjects();
  virtual void   Terminate(Option_t *option);

  void SetFillTree(Bool_t fill=kTRUE){
    fFillTree=fill;
  }  
  void SetReadMC(Bool_t optMC=kTRUE){
    fReadMC=optMC;
  }
  void SetUseMCtruthForPID(Bool_t opt=kTRUE){
    fUseMCId=opt;
  }
  void SetUseGenPtInPlots(Bool_t opt=kTRUE){
    fUseGenPt=opt;
  }
  void SetUsePhysicsSelection(Bool_t opt=kTRUE){
    fUsePhysSel=opt;
  }
  void SetTriggerMask(Int_t mask){
    fTriggerMask=mask;
  }
  void SetUsePileupCut(Bool_t opt=kTRUE){
    fUsePileupCut=kTRUE;
  }
  void SetTPCTrackCuts(AliESDtrackCuts* cuts){
    if(fTrCutsTPC) delete fTrCutsTPC;
    fTrCutsTPC=new AliESDtrackCuts(*cuts);
  }
  void SeMinNumOfTPCPIDclu(Int_t minc){
    fMinNumOfTPCPIDclu=minc;
  }
  void SetUseTOFbcSelection(Bool_t opt){
    fUseTOFbcSelection=opt;
  }
  void SetPtBinning(Int_t nbins, Double_t minpt, Double_t maxpt){
    fNPtBins=nbins; fMinPt=minpt; fMaxPt=maxpt;
  }
  void SetPhiBinning(Int_t nbins){
    fNPhiBins=nbins;
  }
  void SetEtaBinning(Int_t nbins){
    fNEtaBins=nbins;
  }

  AliESDtrackCuts* GetTrackCutObject() const {return fTrCutsTPC;}

 private:

  enum EVarsTree {kNumOfIntVar=14, kNumOfFloatVar=34};

  AliAnalysisTaskCheckESDTracks(const AliAnalysisTaskCheckESDTracks &source);
  AliAnalysisTaskCheckESDTracks& operator=(const AliAnalysisTaskCheckESDTracks &source);
  
  TList*  fOutput;                   //!<!  list of output histos

  TH1F* fHistNEvents;                //!<!  histo with N of events
  TH1F* fHistNTracks;                //!<!  histo with N of tracks
  TH1F* fHistNTracksBackg;           //!<!  histo with N of background tracks
  TH1F* fHistNTracksEmbed;           //!<!  histo with N of embedded tracks
  TH1F* fHistNV0Daughters;                //!<!  histo with N of V0-tracks
  TH1F* fHistNV0DaughtersBackg;           //!<!  histo with N of background V0-tracks
  TH1F* fHistNV0DaughtersEmbed;           //!<!  histo with N of embedded V0-tracks
  TH1F* fHistCheckK0SelBackg;     //!<!  histo to check K0s selection
  TH1F* fHistCheckK0SelEmbed;     //!<!  histo to check K0s selection
  TH1F* fHistCheckK0SelMixed;     //!<!  histo to check K0s selection
  
  TH1F* fHistNITSClu;             //!<!  histo with N of ITS clusters
  TH1F* fHistCluInITSLay;        //!<!  histo with cluters in ITS layers
  
  TH2F* fHistNtracksTPCselVsV0befEvSel;    //!<!  histo of tracks vs. centr.
  TH2F* fHistNtracksSPDanyVsV0befEvSel;    //!<!  histo of tracks vs. centr.
  TH2F* fHistNtracksTPCselVsV0aftEvSel;    //!<!  histo of tracks vs. centr.
  TH2F* fHistNtracksSPDanyVsV0aftEvSel;    //!<!  histo of tracks vs. centr.

  TH2F* fHistdEdxVsP[9];              //!<!  histo of dE/dx for hypos (all tracks)
  TH2F* fHistdEdxVsPTPCsel[9];        //!<!  histo of dE/dx for hypos (TPC cuts)
  TH2F* fHistdEdxVsPTPCselITSref[9];  //!<!  histo of dE/dx for hypos (ITSrefit)
  TH2F* fHistdEdxVsP0[9];              //!<!  histo of dE/dx for hypos (all tracks)
  TH2F* fHistdEdxVsPTPCsel0[9];        //!<!  histo of dE/dx for hypos (TPC cuts)
  TH2F* fHistdEdxVsPTPCselITSref0[9];  //!<!  histo of dE/dx for hypos (ITSrefit)
  TH2F* fHistCorrelHypo0HypoTPCsel;        //!<!  correl. f PID hypos in tracking steps
  TH2F* fHistCorrelHypo0HypoTPCselITSref;  //!<!  correl. f PID hypos in tracking steps

  TH2F* fHistnSigmaVsPdEdxTPCsel[9];  //!<!  histo of nSigma for particle species

  TH3F* fHistEtaPhiPtTPCsel;         //!<!  histo of eta,phi,pt (TPC cuts)
  TH3F* fHistEtaPhiPtTPCselITSref;   //!<!  histo of eta,phi,pt (ITSrefit)
  TH3F* fHistEtaPhiPtTPCselSPDany;   //!<!  histo of eta,phi,pt (ITSrefit+SPDany)
    
  TH3F* fHistEtaPhiPtPosChargeTPCsel;         //!<!  histo of eta,phi,pt (TPC cuts)
  TH3F* fHistEtaPhiPtPosChargeTPCselITSref;   //!<!  histo of eta,phi,pt (ITSrefit)
  TH3F* fHistEtaPhiPtPosChargeTPCselSPDany;   //!<!  histo of eta,phi,pt (ITSrefit+SPDany)
  TH3F* fHistEtaPhiPtNegChargeTPCsel;         //!<!  histo of eta,phi,pt (TPC cuts)
  TH3F* fHistEtaPhiPtNegChargeTPCselITSref;   //!<!  histo of eta,phi,pt (ITSrefit)
  TH3F* fHistEtaPhiPtNegChargeTPCselSPDany;   //!<!  histo of eta,phi,pt (ITSrefit+SPDany)  
  TH3F* fHistEtaPhiPositionPtPosChargeTPCsel;   //!<!  histo of eta,phi (position) at TPC inner
  TH3F* fHistEtaPhiPositionPtNegChargeTPCsel;   //!<!  histo of eta,phi (position) at TPC inner

  TH3F* fHistEtaPhiPtTPCselTOFbc;         //!<!  histo of eta,phi,pt (TPC cuts)
  TH3F* fHistEtaPhiPtTPCselITSrefTOFbc;   //!<!  histo of eta,phi,pt (ITSrefit)
  TH3F* fHistEtaPhiPtTPCselSPDanyTOFbc;   //!<!  histo of eta,phi,pt (ITSrefit+SPDany)

  TH3F* fHistEtaPhiPtInnerTPCsel;         //!<!  histo of eta,phi,pt (TPC cuts)
  TH3F* fHistEtaPhiPtInnerTPCselITSref;   //!<!  histo of eta,phi,pt (ITSrefit)
  TH3F* fHistEtaPhiPtInnerTPCselSPDany;   //!<!  histo of eta,phi,pt (ITSrefit+SPDany)
  TH3F* fHistEtaPhiPtInnerTPCselTOFbc;         //!<!  histo of eta,phi,pt (TPC cuts)
  TH3F* fHistEtaPhiPtInnerTPCselITSrefTOFbc;   //!<!  histo of eta,phi,pt (ITSrefit)
  TH3F* fHistEtaPhiPtInnerTPCselSPDanyTOFbc;   //!<!  histo of eta,phi,pt (ITSrefit+SPDany)

  TH2F* fHistNtrackeltsPtTPCsel;         //!<!  histo of eta,phi,pt (TPC cuts)
  TH2F* fHistNtrackeltsPtTPCselITSref;   //!<!  histo of eta,phi,pt (ITSrefit)
  TH2F* fHistNtrackeltsPtTPCselSPDany;   //!<!  histo of eta,phi,pt (ITSrefit+SPDany)
  TH2F* fHistNtrackeltsPtTPCselTOFbc;         //!<!  histo of eta,phi,pt (TPC cuts)
  TH2F* fHistNtrackeltsPtTPCselITSrefTOFbc;   //!<!  histo of eta,phi,pt (ITSrefit)
  TH2F* fHistNtrackeltsPtTPCselSPDanyTOFbc;   //!<!  histo of eta,phi,pt (ITSrefit+SPDany)

  TH3F* fHistTPCchi2PerClusPhiPtTPCsel;        //!<!  histo of chi2 vs. pt and phi;
  TH3F* fHistTPCchi2PerClusPhiPtTPCselITSref;  //!<!  histo of chi2 vs. pt and phi;
  TH3F* fHistTPCchi2PerClusPhiPtTPCselSPDany;  //!<!  histo of chi2 vs. pt and phi;

  TH3F* fHistSig1ptCovMatPhiPtTPCsel;        //!<!  histo of sigma 1/pt vs. pt and phi;
  TH3F* fHistSig1ptCovMatPhiPtTPCselITSref;  //!<!  histo of sigma 1/pt vs. pt and phi;
  TH3F* fHistSig1ptCovMatPhiPtTPCselSPDany;  //!<!  histo of sigma 1/pt vs. pt and phi;

  // Pi,K,p with good hypothesis in tracking
  TH3F* fHistEtaPhiPtGoodHypTPCsel[3];             //!<!  histo of eta,phi,pt (TPC cuts)
  TH3F* fHistEtaPhiPtGoodHypTPCselITSref[3];       //!<!  histo of eta,phi,pt (ITSrefit)
  TH3F* fHistEtaPhiPtGoodHypTPCselSPDany[3];       //!<!  histo of eta,phi,pt (ITSrefit+SPDany)
  TH3F* fHistEtaPhiPtInnerGoodHypTPCsel[3];        //!<!  histo of eta,phi,pt (TPC cuts)
  TH3F* fHistEtaPhiPtInnerGoodHypTPCselITSref[3];  //!<!  histo of eta,phi,pt (ITSrefit)
  TH3F* fHistEtaPhiPtInnerGoodHypTPCselSPDany[3];  //!<!  histo of eta,phi,pt (ITSrefit+SPDany)
  TH3F* fHistTPCchi2PerClusPhiPtGoodHypTPCsel[3];       //!<!  histo of chi2 vs. pt and phi;
  TH3F* fHistTPCchi2PerClusPhiPtGoodHypTPCselITSref[3]; //!<!  histo of chi2 vs. pt and phi;
  TH3F* fHistTPCchi2PerClusPhiPtGoodHypTPCselSPDany[3]; //!<!  histo of chi2 vs. pt and phi;
  TH2F* fHistdEdxVsPGoodHyp[3];                         //!<!  histo of dE/dx for protons
  TH3F* fHistImpParXYPtMulGoodHypTPCselSPDany[3];       //!<!  histo of impact parameter

  // Pi,K,p with bad hypothesis in tracking
  TH3F* fHistEtaPhiPtBadHypTPCsel[3];         //!<!  histo of eta,phi,pt (TPC cuts)
  TH3F* fHistEtaPhiPtBadHypTPCselITSref[3];   //!<!  histo of eta,phi,pt (ITSrefit)
  TH3F* fHistEtaPhiPtBadHypTPCselSPDany[3];   //!<!  histo of eta,phi,pt (ITSrefit+SPDany)
  TH3F* fHistEtaPhiPtInnerBadHypTPCsel[3];         //!<!  histo of eta,phi,pt (TPC cuts)
  TH3F* fHistEtaPhiPtInnerBadHypTPCselITSref[3];   //!<!  histo of eta,phi,pt (ITSrefit)
  TH3F* fHistEtaPhiPtInnerBadHypTPCselSPDany[3];   //!<!  histo of eta,phi,pt (ITSrefit+SPDany)
  TH3F* fHistTPCchi2PerClusPhiPtBadHypTPCsel[3];        //!<!  histo of chi2 vs. pt and phi;
  TH3F* fHistTPCchi2PerClusPhiPtBadHypTPCselITSref[3];  //!<!  histo of chi2 vs. pt and phi;
  TH3F* fHistTPCchi2PerClusPhiPtBadHypTPCselSPDany[3];  //!<!  histo of chi2 vs. pt and phi;
  TH2F* fHistdEdxVsPBadHyp[3];                //!<!  histo of dE/dx for protons
  TH3F* fHistImpParXYPtMulBadHypTPCselSPDany[3];   //!<!  histo of impact parameter

  TH2F* fHistPtResidVsPtTPConlyAll;                             //!<!  Pt residuals for TPC only tracks
  TH2F* fHistPtResidVsPtTPCselAll;                              //!<!  Pt residuals for TPC+ITS track with TPC only cuts
  TH2F* fHistPtResidVsPtTPCselITSrefAll;                        //!<!  Pt residuals for ITS+TPC tracks 
  TH2F* fHistPtResidVsPtTPCselGoodHyp[AliPID::kSPECIESC];       //!<!  Pt residuals for TPC only tracks tracked with good mass hypothesis (for each species)
  TH2F* fHistPtResidVsPtTPCselBadHyp[AliPID::kSPECIESC];        //!<!  Pt residuals for TPC only tracks tracked with bad mass hypothesis (for each species)
  TH2F* fHistPtResidVsPtTPCselITSrefGoodHyp[AliPID::kSPECIESC]; //!<!  Pt residuals for ITS+TPC tracks tracked with good mass hypothesis (for each species)
  TH2F* fHistPtResidVsPtTPCselITSrefBadHyp[AliPID::kSPECIESC];  //!<!  Pt residuals for ITS+TPC tracks tracked with bad mass hypothesis (for each species)
  TH2F* fHistOneOverPtResidVsPtTPConlyAll;                             //!<!  1/Pt residuals for TPC only tracks
  TH2F* fHistOneOverPtResidVsPtTPCselAll;                              //!<!  1/Pt residuals for TPC+ITS track with TPC only cuts
  TH2F* fHistOneOverPtResidVsPtTPCselITSrefAll;                        //!<!  1/Pt residuals for ITS+TPC tracks
  TH2F* fHistOneOverPtResidVsPtTPCselGoodHyp[AliPID::kSPECIESC];       //!<!  1/Pt residuals for TPC only tracks tracked with good mass hypothesis (for each species)
  TH2F* fHistOneOverPtResidVsPtTPCselBadHyp[AliPID::kSPECIESC];        //!<!  1/Pt residuals for TPC only tracks tracked with bad mass hypothesis (for each species)
  TH2F* fHistOneOverPtResidVsPtTPCselITSrefGoodHyp[AliPID::kSPECIESC]; //!<!  1/Pt residuals for ITS+TPC tracks tracked with good mass hypothesis (for each species)
  TH2F* fHistOneOverPtResidVsPtTPCselITSrefBadHyp[AliPID::kSPECIESC];  //!<!  1/Pt residuals for ITS+TPC tracks tracked with bad mass hypothesis (for each species)
  TH2F* fHistPzResidVsPtTPCselAll;                       //!<!  Pz residuals for TPC only tracks tracked with good mass hypothesis
  TH2F* fHistPzResidVsPtTPCselITSrefAll;                 //!<!  Pz residuals for ITS+TPC tracks tracked with good mass hypothesis
  TH2F* fHistPzResidVsEtaTPCselAll;                      //!<!  Pz residuals for TPC only tracks tracked with good mass hypothesis
  TH2F* fHistPzResidVsEtaTPCselITSrefAll;                //!<!  Pz residuals for ITS+TPC tracks tracked with good mass hypothesis

  TH3F* fHistEtaPhiPtTPCselITSrefGood;        //!<!  histo of eta,phi,pt - good MC tracks
  TH3F* fHistEtaPhiPtTPCselITSrefFake;        //!<!  histo of eta,phi,pt - fake MC tracks
  TH3F* fHistImpParXYPtMulTPCselSPDanyGood;   //!<!  histo of impact parameter (pion)
  TH3F* fHistImpParXYPtMulTPCselSPDanyFake;   //!<!  histo of impact parameter (pion)
  TH3F* fHistImpParXYPtMulTPCselSPDanyPrim;   //!<!  histo of impact parameter (pion)
  TH3F* fHistImpParXYPtMulTPCselSPDanySecDec;   //!<!  histo of impact parameter (pion)
  TH3F* fHistImpParXYPtMulTPCselSPDanySecMat;   //!<!  histo of impact parameter (pion)

  TH3F* fHistInvMassK0s;
  TH3F* fHistInvMassLambda;
  TH3F* fHistInvMassAntiLambda;
  TH3F* fHistInvMassLambdaGoodHyp;     //!<!  histo of lambda inv mass
  TH3F* fHistInvMassAntiLambdaGoodHyp; //!<!  histo of lambdabar inv mass
  TH3F* fHistInvMassLambdaBadHyp;      //!<!  histo of lambda inv mass
  TH3F* fHistInvMassAntiLambdaBadHyp;  //!<!  histo of lambdabar inv mass

  Bool_t   fFillTree;          // flag to control fill of tree
  TTree*   fTrackTree;         //!<! output tree
  Float_t* fTreeVarFloat;      //!<! variables to be written to the tree
  Int_t*   fTreeVarInt;        //!<! variables to be written to the tree


  AliESDtrackCuts* fTrCutsTPC;        // TPC track cuts
  Int_t   fMinNumOfTPCPIDclu;  // cut on min. of TPC clust for PID
  Bool_t  fUseTOFbcSelection;  // flag use/not use TOF for pileup rejection
  Bool_t  fUsePhysSel;         // flag use/not use phys sel
  Bool_t  fUsePileupCut;       // flag use/not use phys pileup cut
  Int_t   fTriggerMask;        // mask used in physics selection
  Int_t fNEtaBins;             // number of eta intervals in histos
  Int_t fNPhiBins;             // number of phi intervals in histos
  Int_t fNPtBins;              // number of pt intervals in histos
  Double_t fMinPt;             // minimum pt for histos
  Double_t fMaxPt;             // maximum pt for histos
  Bool_t  fReadMC;             // flag read/not-read MC truth info
  Bool_t  fUseMCId;            // flag use/not-use MC identity for PID
  Bool_t  fUseGenPt;           // flag for reco/gen pt in plots

  ClassDef(AliAnalysisTaskCheckESDTracks,19);
};


#endif
