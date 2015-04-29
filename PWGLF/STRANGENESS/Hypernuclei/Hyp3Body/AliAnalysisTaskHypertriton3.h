#ifndef ALIANALYSISTASKHYPERTRITON3_H
#define ALIANALYSISTASKHYPERTRITON3_H


/**************************************************************************
 *                                                                        *                                                                         
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/



///////////////////////////////////////////////////////////////////////////
// AliAnalysisTaskHypertriton3 class
// analysis task for the study of the production of hypertriton
// which decays in 3 prongs: d+p+pi^-
// This task is optimized for ESDs.root
//
// Author: 
// S. Trogolo, trogolo@to.infn.it
///////////////////////////////////////////////////////////////////////////

class TChain;
class TF1;
class TH1F;
class TH2F;
class TList;
class TLegend;
class TParticle;
class TRandom3;
class TString;
class TTree;

class AliESDEvent;
class AliESDtrackCuts;
class AliESDv0Cuts;
class AliESDVertex;
class AliMCEvent;
class AliMCEventHandler;
class AliPIDResponse;
class AliVertexerTracks;

#include "AliAnalysisTaskSE.h" 

class AliAnalysisTaskHypertriton3 : public AliAnalysisTaskSE {

 public:
  AliAnalysisTaskHypertriton3();
  virtual ~AliAnalysisTaskHypertriton3();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t*);
  virtual void   Terminate(Option_t*);

  void SetReadMC(Bool_t flag = kTRUE) {fMC = flag;}
  void SetFillTree(Bool_t outTree = kFALSE) {fFillTree = outTree;}

  Double_t GetDCAcut(Int_t part, Double_t dca) const;

 private:

  AliESDEvent       *fESDevent; // ESD event
  AliESDtrackCuts   *fESDtrackCuts; // track cuts
  AliESDtrackCuts   *fESDtrackCutsV0; // cut variables for V0's
  AliESDVertex      *fPrimaryVertex;

  AliPIDResponse    *fPIDResponse;  

  //Variables
  Bool_t             fMC; // variables for MC selection
  Bool_t             fFillTree;
  Float_t            fCentrality; // Centrality class
  Float_t            fCentralityPercentile; // Centrality percentile

  //Output list
  TList             *fOutput;

  //Histograms
  TH1F              *fHistCount;
  TH1F              *fHistParticle;
  TH1F              *fHistCentralityClass;
  TH1F              *fHistCentralityPercentile;
  TH1F              *fHistTrigger;
  TH1F              *fHistMultiplicity;
  TH1F              *fHistZPrimaryVtx;
  TH1F              *fHistXPrimaryVtx;
  TH1F              *fHistYPrimaryVtx;
  TH1F              *fHistChi2perTPCcluster;


  //PID
  //--> TPC
  TH2F              *fHistTPCpid;
  TH2F              *fHistTPCdeusignal;
  TH2F              *fHistTPCprosignal;
  TH2F              *fHistTPCpionsignal;
  TH2F              *fHistTPCantideusignal;
  TH2F              *fHistTPCantiprosignal;
  TH2F              *fHistTPCpionplussignal;

  //--> TOF
  TH2F              *fHistTOFsignal;
  TH2F              *fHistTOFdeusignal;
  TH2F              *fHistTOFprosignal;
  TH2F              *fHistTOFantideusignal;
  TH2F              *fHistTOFantiprosignal;
  TH1F              *fHistTOFdeumass;
  TH1F              *fHistTOFpromass;

  //Candidate combination
  TH1F              *fHistpionTPCcls;
  TH1F              *fHistpTpion;
  TH2F              *fHistCorrDCAdprimary;
  TH2F              *fHistCorrDCApprimary;
  TH2F              *fHistCorrDCApiprimary;
  TH1F              *fHistDCApiprimary;
  TH1F              *fHistDCAdeupro;
  TH1F              *fHistDCApiondeu;	  
  TH1F              *fHistDCApionpro;
  TH1F              *fHistZDecayVtx;
  TH1F              *fHistXDecayVtx;
  TH1F              *fHistYDecayVtx;
  TH1F              *fHistDecayLengthH3L;
  TH1F              *fHistDCAXYdeuvtx;
  TH1F              *fHistDCAZdeuvtx;
  TH1F              *fHistDCAXYprovtx;
  TH1F              *fHistDCAZprovtx;
  TH1F              *fHistDCAXYpionvtx;
  TH1F              *fHistDCAZpionvtx;
  TH1F              *fHistCosPointingAngle;
  TH1F              *fHistMassHypertriton;
  TH1F              *fHistpionTPCclsMCt;
  TH1F              *fHistpTpionMCt;
  TH1F              *fHistpTproMCt;
  TH1F              *fHistpTdeuMCt;
  TH1F              *fHistDCAdeuproMCt;
  TH1F              *fHistDCApiondeuMCt;
  TH1F              *fHistDCApionproMCt;
  TH2F              *fHistCorrDCAdprimaryMCt;
  TH2F              *fHistCorrDCApprimaryMCt;
  TH2F              *fHistCorrDCApiprimaryMCt;
  TH1F              *fHistDCApiprimaryMCt;
  TH1F              *fHistDCApprimaryMCt;
  TH1F              *fHistDCAdprimaryMCt;
  TH1F              *fHistDCAXYdeuvtxMCt;
  TH1F              *fHistDCAZdeuvtxMCt;
  TH1F              *fHistDCAXYprovtxMCt;
  TH1F              *fHistDCAZprovtxMCt;
  TH1F              *fHistDCAXYpionvtxMCt;
  TH1F              *fHistDCAZpionvtxMCt;
  TH1F              *fHistZDecayVtxMCt;
  TH1F              *fHistXDecayVtxMCt;
  TH1F              *fHistYDecayVtxMCt;
  TH1F              *fHistMassHypertritonMCt;

  //TTree
  TTree             *fTTree;
  // Deuteron
  Float_t           fTchi2NDFdeu;
  UShort_t          fTPCclsdeu;
  UShort_t          fTPCclsPIDdeu;
  Float_t           fTpTPCdeu;
  Float_t           fTpTdeu;
  Float_t           fTpdeu;
  Float_t           fTTPCnsigmadeu;
  Float_t           fTTOFmassdeu;
  Float_t           fTDCAXYdeuprvtx;
  Float_t           fTDCAZdeuprvtx;
  // Proton
  Float_t           fTchi2NDFpro;
  UShort_t          fTPCclspro;
  UShort_t          fTPCclsPIDpro;
  Float_t           fTpTPCpro;
  Float_t           fTpTpro;
  Float_t           fTppro;
  Float_t           fTTPCnsigmapro;
  Float_t           fTTOFmasspro;
  Float_t           fTDCAXYproprvtx;
  Float_t           fTDCAZproprvtx;
  // Pion
  Float_t           fTchi2NDFpion;
  UShort_t          fTPCclspion;
  UShort_t          fTPCclsPIDpion;
  Float_t           fTpTPCpion;
  Float_t           fTpTpion;
  Float_t           fTppion;
  Float_t           fTTPCnsigmapion;
  Float_t           fTDCAXYpioprvtx;
  Float_t           fTDCAZpioprvtx;
  // Triplet
  Float_t           fTDCAdp;
  Float_t           fTDCAdpi;
  Float_t           fTDCAppi;
  
  Float_t           fTDCAXYdvtx;
  Float_t           fTDCAZdvtx;
  Float_t           fTDCAXYpvtx;
  Float_t           fTDCAZpvtx;
  Float_t           fTDCAXYpivtx;
  Float_t           fTDCAZpivtx;
  
  Float_t           fTDecayLength;
  Float_t           fTCosPA;
  Float_t           fTInvariantMass;

  
  AliAnalysisTaskHypertriton3(const AliAnalysisTaskHypertriton3&); // not implemented
  AliAnalysisTaskHypertriton3& operator=(const AliAnalysisTaskHypertriton3&); // not implemented
  
  ClassDef(AliAnalysisTaskHypertriton3, 1); // analysisclass
  
};

#endif
