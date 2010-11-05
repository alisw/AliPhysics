#ifndef ALIANALYSISTASKSPDDNDETA_H
#define ALIANALYSISTASKSPDDNDETA_H

///////////////////////////////////////////////////////////////////////////
// Class AliAnalysisTaskSPDdNdEta                                        //
// Analysis task for dN/dEta reconstruction with the SPD                 //
//                                                                       //
// Author:  M. Nicassio (INFN Bari)                                      //
// Contact: Maria.Nicassio@ba.infn.it, Domenico.Elia@ba.infn.it          //
///////////////////////////////////////////////////////////////////////////

class TH1F; 
class TH2F;
class TH3F;
class AliESDEvent;
class TList;

#include "AliTriggerAnalysis.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskSPDdNdEta : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskSPDdNdEta(const char *name = "AliAnalysisTaskSPDdNdEta");
  virtual ~AliAnalysisTaskSPDdNdEta(); 
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  enum MCCentralityBin{kcentral=1,kbin3to6=2,kbin6to9=3,kbin9to12=4,kbin12to15=5,kperipheral=6,kall=7};

  void SetReadMC(Bool_t readmc = kFALSE) { fUseMC = readmc; }
  void SetTrigger(AliTriggerAnalysis::Trigger trigger) { fTrigger = trigger; }
  void SetReadPbPb(Bool_t pbpb = kFALSE) { fPbPb = pbpb; } 
  void SetReadTrackRefs(Bool_t readtr = kFALSE) { fTR = readtr; }
  void SetRecoTracklets(Bool_t recotracklets = kFALSE) { fRecoTracklets = recotracklets; }
  void SetMCCentralityBin(MCCentralityBin mccentrbin) {fMCCentralityBin=mccentrbin;}
  void SetCentralityLowLim(Float_t centrlowlim) {fCentrLowLim=centrlowlim;}
  void SetCentralityUpLim(Float_t centruplim) {fCentrUpLim=centruplim;}
  void SetCentralityEst(TString centrest) {fCentrEst=centrest;}
 
  Bool_t IsDetectablePrimary(Int_t nref, AliMCParticle* mcpart);
  Bool_t IsDetectedPrimary(Int_t nref, AliMCParticle* mcpart, Int_t Layer);

 protected:
  AliESDEvent *fmyESD;             // ! ESD object 
  TList *fOutput;                  // ! output list send on output slot 1 
  
  MCCentralityBin fMCCentralityBin; // to select MC centrality bin in which corrections are calculated
  Float_t fCentrLowLim;             // to select centrality bin on data
  Float_t fCentrUpLim;              // to select centrality bin on data
  TString fCentrEst;                // to select centrality estimator


  Bool_t fUseMC;                   // flag to enable the calculation of correction histograms
  Bool_t fPbPb;                    // flag to analyze PbPb data 
  AliTriggerAnalysis::Trigger fTrigger;  
  Bool_t fTR;                      // to read track references and calculate factors of track to particle correction 
  Bool_t fRecoTracklets;           // flag to recostruct tracklets

  TH2F        *fHistSPDRAWMultvsZ;          // ! data to be corrected 
  TH2F        *fHistSPDRAWMultvsZTriggCentrEvts; // ! data to be corrected
  TH2F        *fHistSPDRAWMultvsZCentrEvts; // ! data to be corrected

  TH2F        *fHistSPDRAWEtavsZ;           // ! data to be corrected

  TH1F        *fHistSPDmultEtacut;          // ! cluster inner layer and tracklet check histos
  TH1F        *fHistSPDmult;                // ! cluster inner layer and tracklet check histos
  TH1F        *fHistSPDmultcl1;
  TH1F        *fHistSPDeta;                 // ! cluster inner layer and tracklet check histos
  TH1F        *fHistSPDphi;                 // ! cluster inner layer and tracklet check histos
  TH1F        *fHistSPDtheta;               // ! cluster inner layer and tracklet check histos
  TH1F        *fHistSPDdePhi;               // ! cluster inner layer and tracklet check histos
  TH2F        *fHistSPDphivsSPDeta;         // ! cluster inner layer and tracklet check histos
  TH1F        *fHistSPDdeTheta;             // ! cluster inner layer and tracklet check histos

  TH1F        *fHistSPDvtxAnalysis;         // ! SPD vertex distributions
  TH2F        *fHistSPDdePhideTheta;        // ! histogram for combinatorial background studies

  TH2F* fHistBkgCorrDen;             // ! track level correction histograms
  TH2F* fHistBkgCorrDenPrimGen;      // ! track level correction histograms
  TH2F* fHistBkgCorrNum;             // ! track level correction histograms
  TH2F* fHistAlgEffNum;              // ! track level correction histograms
  TH2F* fHistNonDetectableCorrDen;   // ! track level correction histograms

  TH2F* fHistNonDetectableCorrNum;   // ! track level correction histograms
  TH2F* fHistAllPrimaries;           // ! track level correction histograms
  TH2F* fHistTrackCentrEvts;         // ! track level correction histograms
  TH2F* fHistTrackTrigCentrEvts;     // ! track level correction histograms
 
  TH2F* fHistAllEvts;                // ! event level correction histograms
  TH2F* fHistCentrEvts;              // ! event level correction histograms
  TH2F* fHistTrigCentrEvts;          // ! event level correction histograms
  TH2F* fHistSelEvts;                // ! event level correction histograms

  TH1F* fHistMCmultEtacut;           // ! MC distributions
  TH2F* fHistMCmultEtacutvsSPDmultEtacut; // ! MC distributions

  TH1F* fHistMCvtxx;                 // ! MC vertex
  TH1F* fHistMCvtxy;                 // ! MC vertex
  TH1F* fHistMCvtxz;                 // ! MC vertex

  TH2F* fHistRecvsGenImpactPar;      // ! impact parameter correlation (ZDC est vs gen)
  TH1F* fHistMCNpart;                // ! distribution of number of participants from MC 

  TH1F* fHistdPhiPP;                 // ! tracklet check histo
  TH1F* fHistdPhiSS;                 // ! tracklet check histo
  TH1F* fHistdPhiComb;               // ! tracklet check histo

  TH2F* fHistDeVtx;                  // ! check histo

 private:    
  AliAnalysisTaskSPDdNdEta(const AliAnalysisTaskSPDdNdEta&); // not implemented
  AliAnalysisTaskSPDdNdEta& operator=(const AliAnalysisTaskSPDdNdEta&); // not implemented 
  
  ClassDef(AliAnalysisTaskSPDdNdEta, 3);  
};

#endif
