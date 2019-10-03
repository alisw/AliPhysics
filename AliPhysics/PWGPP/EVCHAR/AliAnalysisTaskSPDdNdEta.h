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
class AliTrackletAlg;
class AliMCParticle;

#include "AliAnalysisTaskSE.h"
#include "AliTriggerAnalysis.h"

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
  void SetReadTrackRefs(Bool_t readtr = kFALSE) { fTR = readtr; }
  void SetRecoTracklets(Bool_t recotracklets = kFALSE) { fRecoTracklets = recotracklets; }
  void SetMCCentralityBin(MCCentralityBin mccentrbin) {fMCCentralityBin=mccentrbin;}
  void SetCentralityLowLim(Float_t centrlowlim) {fCentrLowLim=centrlowlim;}
  void SetCentralityUpLim(Float_t centruplim) {fCentrUpLim=centruplim;}
//  void SetCentralityEst(TString centrest) {fCentrEst=centrest;}
  void SetCentralityEst(Bool_t centrest) {fCentrEst=centrest;}
  void SetMinClusterMultLay2(Int_t minClMultLay2=0) {fMinClMultLay2=minClMultLay2;}
  void SetMaxClusterMultLay2(Int_t maxClMultLay2=0) {fMaxClMultLay2=maxClMultLay2;}
  void SetMinV0Mult(Int_t minV0Mult=0) {fMinMultV0=minV0Mult;}
  void SetVertexRange(Float_t vl=7.) {fVtxLim=vl;}
  void SetPartSpecies(Bool_t partsp = kFALSE)  {fPartSpecies=partsp;}

  void SetPhiWindow(Float_t w=0.08) {fPhiWindow=w;}
  void SetThetaWindow(Float_t w=0.025) {fThetaWindow=w;}
  void SetPhiShift(Float_t w=0.0045) {fPhiShift=w;}
  void SetRemoveClustersFromOverlaps(Bool_t b = kFALSE) {fRemoveClustersFromOverlaps = b;}
  void SetPhiOverlapCut(Float_t w=0.005) {fPhiOverlapCut=w;}
  void SetZetaOverlapCut(Float_t w=0.05) {fZetaOverlapCut=w;}
  void SetPhiRotationAngle(Float_t w=0.0) {fPhiRotationAngle=w;}
  void SetPhiWindowAna(Float_t w=0.08) {fPhiWindowAna=w;}
 
  Bool_t IsDetectablePrimary(Int_t nref, AliMCParticle* mcpart);
  Bool_t IsDetectedPrimary(Int_t nref, AliMCParticle* mcpart, Int_t Layer);

 protected:
  AliESDEvent *fmyESD;             // ! ESD object 
  TList *fOutput;                  // ! output list send on output slot 1 

  Bool_t fUseMC;                   // flag to enable the calculation of correction histograms
  AliTriggerAnalysis::Trigger fTrigger;  
  Bool_t fTR;                      // to read track references and calculate factors of track to particle correction 
  Bool_t fRecoTracklets;           // flag to recostruct tracklets
  
  MCCentralityBin fMCCentralityBin; // to select MC centrality bin in which corrections are calculated
  Float_t fCentrLowLim;             // to select centrality bin on data
  Float_t fCentrUpLim;              // to select centrality bin on data
//  TString fCentrEst;                 // to select centrality estimator
  Bool_t fCentrEst;                 // to select centrality estimator

  Int_t fMinClMultLay2;             // to select multiplicity class
  Int_t fMaxClMultLay2;             // to select multiplicity class
  Int_t fMinMultV0;                 // to select centrality class 
  Float_t fVtxLim;                  // to select vertex range
  Bool_t fPartSpecies;              // to fill correction matrices for each part species (for syst studies)
 
  Float_t       fPhiWindow;                    // Search window in phi
  Float_t       fThetaWindow;                  // Search window in theta
  Float_t       fPhiShift;                     // Phi shift reference value (at 0.5 T)
  Bool_t        fRemoveClustersFromOverlaps;   // Option to skip clusters in the overlaps
  Float_t       fPhiOverlapCut;                // Fiducial window in phi for overlap cut
  Float_t       fZetaOverlapCut;               // Fiducial window in eta for overlap cut
  Float_t       fPhiRotationAngle;             // Angle to rotate the inner layer cluster for combinatorial reco only
  Float_t       fPhiWindowAna;                 // Final analysis tracklet definition window in phi

  AliTrackletAlg *fMultReco;       // tracklet reconstruction class


  TH1F        *fV0Ampl;                     // ! V0 amplitudes to cut on centrality
  TH2F        *fHistSPDRAWMultvsZ;          // ! data to be corrected 
  TH2F        *fHistSPDRAWMultvsZTriggCentrEvts; // ! data to be corrected
  TH2F        *fHistSPDRAWMultvsZCentrEvts; // ! data to be corrected

  TH2F        *fHistSPDRAWEtavsZ;           // ! data to be corrected

  TH1F        *fHistSPDmultEtacut;          // ! cluster inner layer and tracklet check histos
  TH1F        *fHistSPDmult;                // ! cluster inner layer and tracklet check histos
  TH1F        *fHistSPDmultcl1;             // ! cluster inner layer and tracklet check histos
  TH1F        *fHistSPDmultcl2;             // ! cluster inner layer and tracklet check histos
  TH2F        *fHistSPDmultcl1vscl2;        // ! cluster inner layer and tracklet check histos
  TH2F        *fHistSPDmultvscl1;           // ! cluster inner layer and tracklet check histos
  TH2F        *fHistSPDmultvscl2;           // ! cluster inner layer and tracklet check histos

  TH1F        *fHistSPDeta;                 // ! cluster inner layer and tracklet check histos
  TH1F        *fHistSPDphi;                 // ! cluster inner layer and tracklet check histos
  TH1F        *fHistSPDtheta;               // ! cluster inner layer and tracklet check histos
  TH1F        *fHistSPDdePhi;               // ! cluster inner layer and tracklet check histos
  TH2F        *fHistSPDphivsSPDeta;         // ! cluster inner layer and tracklet check histos
  TH1F        *fHistSPDdeTheta;             // ! cluster inner layer and tracklet check histos
  
  TH1F        *fHistSPDvtx;                 // ! SPD vertex distributions
  TH1F        *fHistSPDvtxAnalysis;         // ! SPD vertex distributions
  TH2F        *fHistSPDdePhideTheta;        // ! histogram for combinatorial background studies

  TH1F        *fHistSPDphicl1;              // ! cluster inner layer and tracklet check histos
  TH1F        *fHistSPDphicl2;              // ! cluster inner layer and tracklet check histos

  TH2F* fHistBkgCorrDen;                    // ! track level correction histograms
  TH2F* fHistReconstructedProtons;            // ! track level correction histograms
  TH2F* fHistReconstructedKaons;            // ! track level correction histograms
  TH2F* fHistReconstructedPions;            // ! track level correction histograms
  TH2F* fHistReconstructedOthers;           // ! track level correction histograms
  TH2F* fHistReconstructedSec;        // ! track level correction histograms

  TH2F* fHistBkgCorrDenPrimGen;      // ! track level correction histograms
  TH2F* fHistBkgCombLabels;          // ! track level correction histograms
  TH2F* fHistBkgCorrNum;             // ! track level correction histograms
  TH2F* fHistAlgEffNum;              // ! track level correction histograms
  TH2F* fHistNonDetectableCorrDen;   // ! track level correction histograms

  TH2F* fHistNonDetectableCorrNum;   // ! track level correction histograms
  TH2F* fHistProtons;                // ! track level correction histograms
  TH2F* fHistKaons;                  // ! track level correction histograms
  TH2F* fHistPions;                  // ! track level correction histograms
  TH2F* fHistOthers;                 // ! track level correction histograms
  TH2F* fHistAllPrimaries;           // ! track level correction histograms
  TH2F* fHistTrackCentrEvts;         // ! track level correction histograms
  TH2F* fHistTrackTrigCentrEvts;     // ! track level correction histograms
 
  TH2F* fHistAllEvts;                // ! event level correction histograms
  TH2F* fHistCentrEvts;              // ! event level correction histograms
  TH2F* fHistTrigCentrEvts;          // ! event level correction histograms
  TH2F* fHistSelEvts;                // ! event level correction histograms

  TH1F* fHistMCmultEtacut;                // ! MC distributions
  TH2F* fHistMCmultEtacutvsSPDmultEtacut; // ! MC distributions
  TH2F* fHistMCmultEtacutvsSPDmultcl1;    // ! MC distributions
  TH2F* fHistMCmultEtacutvsSPDmultcl2;    // ! MC distributions

  TH1F* fHistMCvtxx;                 // ! MC vertex
  TH1F* fHistMCvtxy;                 // ! MC vertex
  TH1F* fHistMCvtxz;                 // ! MC vertex

  TH2F* fHistRecvsGenImpactPar;      // ! impact parameter correlation (ZDC est vs gen)
  TH1F* fHistMCNpart;                // ! distribution of number of participants from MC 

  TH2F* fHistdPhidThetaPP;           // ! tracklet check histo
  TH2F* fHistdPhidThetaSS;           // ! tracklet check histo
  TH2F* fHistdPhidThetaComb;         // ! tracklet check histo

  TH2F* fHistDeVtx;                  // ! check histo

 private:    
  AliAnalysisTaskSPDdNdEta(const AliAnalysisTaskSPDdNdEta&); // not implemented
  AliAnalysisTaskSPDdNdEta& operator=(const AliAnalysisTaskSPDdNdEta&); // not implemented 
  
  ClassDef(AliAnalysisTaskSPDdNdEta, 5);  
};
#endif
