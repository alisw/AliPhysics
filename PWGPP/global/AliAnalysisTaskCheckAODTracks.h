#ifndef ALIANALYSISTASKCHECKTRACKS
#define ALIANALYSISTASKCHECKTRACKS

/* Copyright(c) 1998-2012, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// Class AliAnalysisTaskCheckAODTracks
// AliAnalysisTaskSE to extract QA and performance histos for tracks
// 
//
// Author:
//          F. Prino, prino@to.infn.it
//          
//*************************************************************************

class TList;
class TTree;
class TH1F;
class TH2F;
class TH3F;
class TString;
class AliAODTrack;
class AliESDVertex;

#include "AliESDtrackCuts.h"
#include "AliAnalysisTaskSE.h"
#include "AliPID.h"

class AliAnalysisTaskCheckAODTracks : public AliAnalysisTaskSE {

 public:
  
  AliAnalysisTaskCheckAODTracks();
  virtual ~AliAnalysisTaskCheckAODTracks();

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
  void SetCentralityInterval(Double_t minc, Double_t maxc, TString estim="V0M"){
    fSelectOnCentrality=kTRUE;
    fMinCentrality=minc;
    fMaxCentrality=maxc;
    fCentrEstimator=estim.Data();
  }
  void SetUsePileupCut(Bool_t opt=kTRUE){
    fUsePileupCut=opt;
  }
  void SetTPCTrackCuts(AliESDtrackCuts* cuts){
    if(fTrCutsTPC) delete fTrCutsTPC;
    fTrCutsTPC=new AliESDtrackCuts(*cuts);
  }
  void SeMinNumOfTPCPIDclu(Int_t minc){
    fMinNumOfTPCPIDclu=minc;
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
  void SetUpperMultiplicity(Double_t maxMult){
    fMaxMult=maxMult;
  }
  void SetRequireITSrefitForV0Daughters(Bool_t opt){
    if(opt) fRequireITSforV0dau |= (1<<kBitRequireITSrefit);
    else fRequireITSforV0dau &= ~(1<<kBitRequireITSrefit);
  }
  void SetRequireSPDanyForV0Daughters(Bool_t opt){
    if(opt) fRequireITSforV0dau |= (1<<kBitRequireSPDany);
    else fRequireITSforV0dau &= ~(1<<kBitRequireSPDany);
  }
  AliESDtrackCuts* GetTPCTrackCuts(){return fTrCutsTPC;}

  Bool_t ConvertAndSelectAODTrack(AliAODTrack* aTrack, const AliESDVertex vESD, Double_t magField);



 private:

  enum EVarsTree {kNumOfIntVar=12, kNumOfFloatVar=27};
  enum EITSRequirements {kBitRequireITSrefit=0, kBitRequireSPDany=1};
  enum EFiltBits {kNumOfFilterBits=12};

  AliAnalysisTaskCheckAODTracks(const AliAnalysisTaskCheckAODTracks &source);
  AliAnalysisTaskCheckAODTracks& operator=(const AliAnalysisTaskCheckAODTracks &source);
  
  TList*  fOutput;                   //!<!  list of output histos

  TH1F* fHistNEvents;                  //!<!  histo with N of events  
  TH1F* fHistNTracks;                  //!<!  histo with N of tracks
  TH2F* fHistNTracksVsTPCclusters;     //!<! histos of track-cluster correlations
  TH2F* fHistNTracksVsITSclusters;     //!<! histos of track-cluster correlations
  TH2F* fHistITSclustersVsTPCclusters; //!<! histos of track-cluster correlations
  TH2F* fHistNTracksFB4VsTPCclusters;  //!<! histos of track-cluster correlations
  TH2F* fHistNTracksFB4VsITSclusters;  //!<! histos of track-cluster correlations

  TH2D* fHistFilterBits;             //!<!  histo of fieter bits

  TH1F* fHistITSnClusTPCsel;        //!<! histo of ITS clusters
  TH1F* fHistITSnClusITSsa;         //!<! histo of ITS clusters
  TH1F* fHistITSnClusITSPureSA;     //!<! histo of ITS clusters
  TH1F* fHistITSCluInLayTPCsel;     //!<! histo of track pts in ITS layers
  TH1F* fHistITSCluInLayITSsa;      //!<! histo of track pts in ITS layers
  TH1F* fHistITSCluInLayITSPureSA;  //!<! histo of track pts in ITS layers

  TH2F* fHistNtracksFb4VsV0befEvSel;    //!<!  histo of tracks vs. centr.
  TH2F* fHistNtracksFb5VsV0befEvSel;    //!<!  histo of tracks vs. centr.
  TH2F* fHistNtracksFb4VsV0aftEvSel;    //!<!  histo of tracks vs. centr.
  TH2F* fHistNtracksFb5VsV0aftEvSel;    //!<!  histo of tracks vs. centr.

  TH1F* fHistNtracksFb0; //!<! histo of track multipl.
  TH1F* fHistNtracksFb1; //!<! histo of track multipl.
  TH1F* fHistNtracksFb4; //!<! histo of track multipl.
  TH1F* fHistNtracksFb5; //!<! histo of track multipl.
  TH1F* fHistNtracksFb6; //!<! histo of track multipl.
  TH1F* fHistNtracksFb7; //!<! histo of track multipl.
  TH1F* fHistNtracksFb8; //!<! histo of track multipl.

  TH3F* fHistEtaPhiPtTPCsel;           //!<!  histo of eta,phi,pt (TPC cuts)
  TH3F* fHistEtaPhiPtTPCselITSref;     //!<!  histo of eta,phi,pt (ITSrefit)
  TH3F* fHistEtaPhiPtTPCselSPDany;     //!<!  histo of eta,phi,pt (ITSrefit+SPDany)

  TH3F* fHistEtaPhiPtPosChargeTPCsel;         //!<!  histo of eta,phi,pt (TPC cuts)
  TH3F* fHistEtaPhiPtPosChargeTPCselITSref;   //!<!  histo of eta,phi,pt (ITSrefit)
  TH3F* fHistEtaPhiPtPosChargeTPCselSPDany;   //!<!  histo of eta,phi,pt (ITSrefit+SPDany)
  TH3F* fHistEtaPhiPtNegChargeTPCsel;         //!<!  histo of eta,phi,pt (TPC cuts)
  TH3F* fHistEtaPhiPtNegChargeTPCselITSref;   //!<!  histo of eta,phi,pt (ITSrefit)
  TH3F* fHistEtaPhiPtNegChargeTPCselSPDany;   //!<!  histo of eta,phi,pt (ITSrefit+SPDany)

  TH3F* fHistEtaPhiPtTPCselTOFbc;         //!<!  histo of eta,phi,pt (TPC cuts)
  TH3F* fHistEtaPhiPtTPCselITSrefTOFbc;   //!<!  histo of eta,phi,pt (ITSrefit)
  TH3F* fHistEtaPhiPtTPCselSPDanyTOFbc;   //!<!  histo of eta,phi,pt (ITSrefit+SPDany)
  TH3F* fHistEtaPhiPtTPCselSPDanyTOFpid;  //!<!  histo of eta,phi,pt (ITSrefit+SPDany) for tracks with TOF PID

  TH3F* fHistTPCchi2PerClusPhiPtTPCsel;        //!<!  histo of chi2 vs. pt and phi;
  TH3F* fHistTPCchi2PerClusPhiPtTPCselITSref;  //!<!  histo of chi2 vs. pt and phi;
  TH3F* fHistTPCchi2PerClusPhiPtTPCselSPDany;  //!<!  histo of chi2 vs. pt and phi;

  TH3F* fHistSig1ptCovMatPhiPtTPCsel;        //!<!  histo of sigma 1/pt vs. pt and phi;
  TH3F* fHistSig1ptCovMatPhiPtTPCselITSref;  //!<!  histo of sigma 1/pt vs. pt and phi;
  TH3F* fHistSig1ptCovMatPhiPtTPCselSPDany;  //!<!  histo of sigma 1/pt vs. pt and phi;

  TH3F* fHistImpParXYPtMulPionTPCselSPDany;    //!<!  histo of impact parameter (pion)
  TH3F* fHistImpParXYPtMulKaonTPCselSPDany;    //!<!  histo of impact parameter (kaon)
  TH3F* fHistImpParXYPtMulProtonTPCselSPDany;  //!<!  histo of impact parameter (proton)

  TH3F* fHistEtaPhiPtFiltBit[kNumOfFilterBits];               //!<!  histo of Eta,phi,pt per filter bit
  TH3F* fHistImpParXYPtMulFiltBit[kNumOfFilterBits];          //!<!  histo of impact parameter per filter bit
  TH2F* fHistITScluPtFiltBit[kNumOfFilterBits];               //!<!  histo of n. ITS clus per filter bit
  TH2F* fHistSPDcluPtFiltBit[kNumOfFilterBits];               //!<!  histo of n. SPD clus per filter bit
  TH2F* fHistTPCcluPtFiltBit[kNumOfFilterBits];               //!<!  histo of n. TPC clus per filter bit
  TH2F* fHistTPCcrrowsPtFiltBit[kNumOfFilterBits];               //!<!  histo of n. TPC crossed rows per filter bit
  TH2F* fHistTPCCrowOverFindPtFiltBit[kNumOfFilterBits];      //!<!  histo of crossedrows/findable per filter bit
  TH2F* fHistTPCChi2clusPtFiltBit[kNumOfFilterBits];           //!<!  histo of TPC chi2 per filter bit
  TH2F* fHistChi2TPCConstrVsGlobPtFiltBit[kNumOfFilterBits];  //!<!  histo of golden chi2 per filter bit
  TH2F* fHistSig1ptCovMatPtFiltBit[kNumOfFilterBits];           //!<!  histo of sig1pt per filter bit

  TH2F* fHistPtResidVsPtTPCselAll;                       //!<!  Pt residuals for TPC only tracks tracked 
  TH2F* fHistPtResidVsPtTPCselITSrefAll;                 //!<!  Pt residuals for ITS+TPC tracks tracked 
  TH2F* fHistOneOverPtResidVsPtTPCselAll;                //!<!  1/Pt residuals for TPC only tracks tracked 
  TH2F* fHistOneOverPtResidVsPtTPCselITSrefAll;          //!<!  1/Pt residuals for ITS+TPC tracks tracked 
  TH2F* fHistPtResidVsPtTPCselPrim;                       //!<!  Pt residuals for TPC only tracks tracked 
  TH2F* fHistPtResidVsPtTPCselITSrefPrim;                 //!<!  Pt residuals for ITS+TPC tracks tracked 
  TH2F* fHistOneOverPtResidVsPtTPCselPrim;                //!<!  1/Pt residuals for TPC only tracks tracked 
  TH2F* fHistOneOverPtResidVsPtTPCselITSrefPrim;          //!<!  1/Pt residuals for ITS+TPC tracks tracked 
  TH2F* fHistPtResidVsPtTPCselSecDec;                       //!<!  Pt residuals for TPC only tracks tracked 
  TH2F* fHistPtResidVsPtTPCselITSrefSecDec;                 //!<!  Pt residuals for ITS+TPC tracks tracked 
  TH2F* fHistOneOverPtResidVsPtTPCselSecDec;                //!<!  1/Pt residuals for TPC only tracks tracked 
  TH2F* fHistOneOverPtResidVsPtTPCselITSrefSecDec;          //!<!  1/Pt residuals for ITS+TPC tracks tracked 
  TH2F* fHistPtResidVsPtTPCselSecMat;                       //!<!  Pt residuals for TPC only tracks tracked 
  TH2F* fHistPtResidVsPtTPCselITSrefSecMat;                 //!<!  Pt residuals for ITS+TPC tracks tracked 
  TH2F* fHistOneOverPtResidVsPtTPCselSecMat;                //!<!  1/Pt residuals for TPC only tracks tracked 
  TH2F* fHistOneOverPtResidVsPtTPCselITSrefSecMat;          //!<!  1/Pt residuals for ITS+TPC tracks tracked 
  TH2F* fHistPtResidVsPtTPCsel[AliPID::kSPECIESC];       //!<!  Pt residuals for TPC only tracks tracked  (for each species)
  TH2F* fHistPtResidVsPtTPCselITSref[AliPID::kSPECIESC]; //!<!  Pt residuals for ITS+TPC tracks tracked  (for each species)
  TH2F* fHistOneOverPtResidVsPtTPCsel[AliPID::kSPECIESC];       //!<!  Pt residuals for TPC only tracks tracked  (for each species)
  TH2F* fHistOneOverPtResidVsPtTPCselITSref[AliPID::kSPECIESC]; //!<!  Pt residuals for ITS+TPC tracks tracked  (for each species)
  TH2F* fHistPzResidVsPtTPCselAll;                       //!<!  Pz residuals for TPC only tracks tracked 
  TH2F* fHistPzResidVsPtTPCselITSrefAll;                 //!<!  Pz residuals for ITS+TPC tracks tracked 
  TH2F* fHistPzResidVsEtaTPCselAll;                      //!<!  Pz residuals for TPC only tracks tracked 
  TH2F* fHistPzResidVsEtaTPCselITSrefAll;                //!<!  Pz residuals for ITS+TPC tracks tracked 

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

  Bool_t   fFillTree;          // flag to control fill of tree
  TTree*   fTrackTree;         //!<! output tree
  Float_t* fTreeVarFloat;      //!<! variables to be written to the tree
  Int_t*   fTreeVarInt;        //!<! variables to be written to the tree


  AliESDtrackCuts* fTrCutsTPC; // TPC track cuts
  Int_t   fMinNumOfTPCPIDclu;  // cut on min. of TPC clust for PID
  Bool_t  fUsePhysSel;         // flag use/not use phys sel
  Bool_t  fUsePileupCut;       // flag use/not use phys pileup cut
  Int_t   fTriggerMask;        // mask used in physics selection
  Bool_t fSelectOnCentrality;  // flag to activeta cut on centrality
  Double_t fMinCentrality;     // centrality: lower limit
  Double_t fMaxCentrality;     // centrality: upper limit
  TString fCentrEstimator;     // centrality: estimator
  Int_t fNEtaBins;             // number of eta intervals in histos
  Int_t fNPhiBins;             // number of phi intervals in histos
  Int_t fNPtBins;              // number of pt intervals in histos
  Double_t fMinPt;             // minimum pt for histos
  Double_t fMaxPt;             // maximum pt for histos
  Double_t fMaxMult;           // upper limit of multiplicity plots
  Int_t   fRequireITSforV0dau; // ITSrefit/SPDany requests for V0 daughters
  Bool_t  fReadMC;             // flag read/not-read MC truth info
  Bool_t  fUseMCId;            // flag use/not-use MC identity for PID
  Bool_t  fUseGenPt;           // flag for reco/gen pt in plots

  ClassDef(AliAnalysisTaskCheckAODTracks,19);
};


#endif
