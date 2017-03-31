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
class AliESDEvent;

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
  void SetUsePhysicsSelection(Bool_t opt=kTRUE){
    fUsePhysSel=opt;
  }
  void SetTriggerMask(Int_t mask){
    fTriggerMask=mask;
  }
  void SetTPCTrackCuts(AliESDtrackCuts* cuts){
    if(fTrCutsTPC) delete fTrCutsTPC;
    fTrCutsTPC=new AliESDtrackCuts(*cuts);
  }
  void SeMinNumOfTPCPIDclu(Int_t minc){
    fMinNumOfTPCPIDclu=minc;
  }

 private:

  enum EVarsTree {kNumOfIntVar=11, kNumOfFloatVar=27};
  enum EFiltBits {kNumOfFilterBits=12};

  AliAnalysisTaskCheckAODTracks(const AliAnalysisTaskCheckAODTracks &source);
  AliAnalysisTaskCheckAODTracks& operator=(const AliAnalysisTaskCheckAODTracks &source);
  
  TList*  fOutput;                   //!<!  list of output histos

  TH1F* fHistNEvents;                //!<!  histo with N of events  
  TH1F* fHistNTracks;                //!<!  histo with N of tracks
  TH2D* fHistFilterBits;             //!<!  histo of fieter bits

  TH2F* fHistNtracksFb4VsV0befEvSel;    //!<!  histo of tracks vs. centr.
  TH2F* fHistNtracksFb5VsV0befEvSel;    //!<!  histo of tracks vs. centr.
  TH2F* fHistNtracksFb4VsV0aftEvSel;    //!<!  histo of tracks vs. centr.
  TH2F* fHistNtracksFb5VsV0aftEvSel;    //!<!  histo of tracks vs. centr.

  TH3F* fHistEtaPhiPtTPCsel;           //!<!  histo of eta,phi,pt (TPC cuts)
  TH3F* fHistEtaPhiPtTPCselITSref;     //!<!  histo of eta,phi,pt (ITSrefit)
  TH3F* fHistEtaPhiPtTPCselSPDany;     //!<!  histo of eta,phi,pt (ITSrefit+SPDany)
  TH3F* fHistEtaPhiPtTPCselTOFbc;         //!<!  histo of eta,phi,pt (TPC cuts)
  TH3F* fHistEtaPhiPtTPCselITSrefTOFbc;   //!<!  histo of eta,phi,pt (ITSrefit)
  TH3F* fHistEtaPhiPtTPCselSPDanyTOFbc;   //!<!  histo of eta,phi,pt (ITSrefit+SPDany)

  TH3F* fHistTPCchi2PerClusPhiPtTPCsel;        //!<!  histo of chi2 vs. pt and phi;
  TH3F* fHistTPCchi2PerClusPhiPtTPCselITSref;  //!<!  histo of chi2 vs. pt and phi;
  TH3F* fHistTPCchi2PerClusPhiPtTPCselSPDany;  //!<!  histo of chi2 vs. pt and phi;

  TH3F* fHistImpParXYPtMulPionTPCselSPDany;    //!<!  histo of impact parameter (pion)
  TH3F* fHistImpParXYPtMulKaonTPCselSPDany;    //!<!  histo of impact parameter (kaon)
  TH3F* fHistImpParXYPtMulProtonTPCselSPDany;  //!<!  histo of impact parameter (proton)

  TH3F* fHistEtaPhiPtFiltBit[kNumOfFilterBits];               //!<!  histo of Eta,phi,pt per filter bit
  TH3F* fHistImpParXYPtMulFiltBit[kNumOfFilterBits];          //!<!  histo of impact parameter per filter bit
  TH2F* fHistITScluPtFiltBit[kNumOfFilterBits];               //!<!  histo of n. ITS clus per filter bit
  TH2F* fHistSPDcluPtFiltBit[kNumOfFilterBits];               //!<!  histo of n. SPD clus per filter bit
  TH2F* fHistTPCcluPtFiltBit[kNumOfFilterBits];               //!<!  histo of n. TPC clus per filter bit
  TH2F* fHistTPCcrrowsPtFiltBit[kNumOfFilterBits];               //!<!  histo of n. TPC clus per filter bit
  TH2F* fHistTPCCrowOverFindPtFiltBit[kNumOfFilterBits];      //!<!  histo of n. TPC clus per filter bit
  TH2F* fHistTPCChi2ndfPtFiltBit[kNumOfFilterBits];           //!<!  histo of n. TPC clus per filter bit
  TH2F* fHistChi2TPCConstrVsGlobPtFiltBit[kNumOfFilterBits];  //!<!  histo of n. TPC clus per filter bit

  TH2F* fHistPtResidVsPtTPCselAll;                       //!<!  Pt residuals for TPC only tracks tracked with good mass hypothesis
  TH2F* fHistPtResidVsPtTPCselITSrefAll;                 //!<!  Pt residuals for ITS+TPC tracks tracked with good mass hypothesis
  TH2F* fHistPtResidVsPtTPCsel[AliPID::kSPECIESC];       //!<!  Pt residuals for TPC only tracks tracked with good mass hypothesis (for each species)
  TH2F* fHistPtResidVsPtTPCselITSref[AliPID::kSPECIESC]; //!<!  Pt residuals for ITS+TPC tracks tracked with good mass hypothesis (for each species)

  TH3F* fHistEtaPhiPtTPCselITSrefGood;        //!<!  histo of eta,phi,pt - good MC tracks
  TH3F* fHistEtaPhiPtTPCselITSrefFake;        //!<!  histo of eta,phi,pt - fake MC tracks
  TH3F* fHistImpParXYPtMulTPCselSPDanyGood;   //!<!  histo of impact parameter (pion)
  TH3F* fHistImpParXYPtMulTPCselSPDanyFake;   //!<!  histo of impact parameter (pion)
  TH3F* fHistImpParXYPtMulTPCselSPDanyPrim;   //!<!  histo of impact parameter (pion)
  TH3F* fHistImpParXYPtMulTPCselSPDanySecDec;   //!<!  histo of impact parameter (pion)
  TH3F* fHistImpParXYPtMulTPCselSPDanySecMat;   //!<!  histo of impact parameter (pion)

  TH2F* fHistInvMassK0s;
  TH3F* fHistInvMassLambda;
  TH3F* fHistInvMassAntiLambda;

  Bool_t   fFillTree;          // flag to control fill of tree
  TTree*   fTrackTree;         //!<! output tree
  Float_t* fTreeVarFloat;      //!<! variables to be written to the tree
  Int_t*   fTreeVarInt;        //!<! variables to be written to the tree


  AliESDtrackCuts* fTrCutsTPC; // TPC track cuts
  Int_t   fMinNumOfTPCPIDclu;  // cut on min. of TPC clust for PID
  Bool_t  fUsePhysSel;         // flag use/not use phys sel
  Int_t   fTriggerMask;        // mask used in physics selection
  Int_t fNPtBins;              // number of pt intervals in histos
  Double_t fMinPt;             // minimum pt for histos
  Double_t fMaxPt;             // maximum pt for histos
  Bool_t  fReadMC;             // flag read/not-read MC truth info
  Bool_t  fUseMCId;            // flag use/not-use MC identity for PID

  ClassDef(AliAnalysisTaskCheckAODTracks,4);
};


#endif
