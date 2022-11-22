#ifndef AliAnalysisTaskGammaPHOSPP_cxx
#define AliAnalysisTaskGammaPHOSPP_cxx

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id: AliAnalysisTaskPi0.h 55336 2012-03-25 16:41:06Z kharlov $ */

// Analysis task for pi0 and eta meson analysis in pp collisions
// Authors: Yuri Kharlov, Dmitri Peressounko

class TObjArray;
class TH1F;
class TH2I;
class TH2F;
class TH3F;
class AliESDtrackCuts;
class AliTriggerAnalysis;
class AliAODCaloCluster ;
class AliAODTrack;
class AliPHOSAodCluster ;
class AliPHOSCalibData ;
class AliPHOSGeometry ;
class AliAODEvent;

#include "TH2I.h"
#include "AliAnalysisTaskSE.h"
#include "AliCaloPhoton.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"
#include "AliAODMCParticle.h"
#include "AliLog.h"

class AliAnalysisTaskGammaPHOSPP : public AliAnalysisTaskSE 
{
public:
  AliAnalysisTaskGammaPHOSPP(const char *name = "AliAnalysisTaskGammaPHOSPP");
  virtual ~AliAnalysisTaskGammaPHOSPP() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  //virtual void   Terminate(Option_t *);
  virtual void FinishTaskOutput();
  void SetBCgap(Double_t bcgap) {fBCgap = bcgap;}   
  void SetRecalib (Int_t mod, Double_t recalib) {  
      if (mod < 1 || mod>5) AliFatal(Form("Wrong module number: %d",mod)) ;     
      else fRecalib[mod-1] = recalib ;   
  }

  /*
  void SetPHOSBadMap(Int_t mod,TH2I * h) {
    if (fPHOSBadMap[mod]) 
      delete fPHOSBadMap[mod] ;
    fPHOSBadMap[mod] = new TH2I(*h) ;
    printf("Set %s \n", fPHOSBadMap[mod]->GetName());  
  }
  */  

  static  Bool_t PythiaInfoFromFile(TString currFile, Float_t & xsec, Float_t & trials) ;
    
private:
  AliAnalysisTaskGammaPHOSPP(const AliAnalysisTaskGammaPHOSPP&); 
  AliAnalysisTaskGammaPHOSPP& operator=(const AliAnalysisTaskGammaPHOSPP&); 
  
  Bool_t AcceptEvent(AliAODEvent *event);
  Int_t  GetEventCentrality(AliAODEvent *event);
  void SelectCluster(AliAODCaloCluster *clu1);
  void FillOnePhotonHistograms(AliCaloPhoton *ph);
  void FillTwoPhotonHistograms();
  void MixPhotons();

 // Bool_t IsGoodChannel(const char * det, Int_t mod,Int_t ix, Int_t iz);
  void FillHistogram(const char * key,Double_t x) const ; //Fill 1D histogram witn name key
  void FillHistogram(const char * key,Double_t x, Double_t y) const ; //Fill 2D histogram witn name key
  void FillHistogram(const char * key,Double_t x, Double_t y, Double_t z) const ; //Fill 3D histogram witn name key
  void PHOSvsEMCALClusters();
  void AnalyzeCells();
  //Bool_t TestLambda(Double_t l1,Double_t l2, Double_t R) ;
  void ProcessMC();
  
  Double_t TestGammaPt(AliCaloPhoton *ph);
  Int_t GetPrimaryLabel(AliVCluster *clu);
  Int_t GetPrimaryLabelAtVertex(AliVCluster *clu);
  Int_t TestTrack(AliAODTrack *track);
  Int_t TestBC(Double_t tof) ;

  Double_t NonlinearCorrection(Double_t en);
  void TestMatchingTrackPID(AliVCluster *clu, Double_t pt);
  
  TString              fCurrFileName;      // current file path name
  Bool_t               fCheckMCCrossSection; // retrieve from the pyxsec.root file only if requested
  Float_t              fAvgTrials;         // avg trials

  TF1   *fWeightFunction;
  Double_t Weight(AliAODMCParticle *particle);
  Bool_t PhotonWithinPeak(Double_t Minv, Double_t pt);

  void GammaEfficiencies();
  void CutEfficiencies();

  void AddQAHistograms();
  void AddClusterHistograms();
  void AddMassHistograms();
  void AddMCHistograms();
  void AddTrackHistograms();

private:
 // AliAODtrackCuts *fAODtrackCuts; // Track cut
  AliESDtrackCuts *fESDtrackCuts; // Track cut
  THashList * fOutputContainer;       //final histogram container
  THashList * fOutputContainer2;       //final histogram container
  TList * fPHOSEvents[10][2] ;    //Container for PHOS photons
  TClonesArray * fPHOSEvent ;     //PHOS photons in current event
  AliPIDResponse *fPIDResponse;    //!

  Int_t fnCINT1B;           // Number of CINT1B triggers
  Int_t fnCINT1A;           // Number of CINT1A triggers
  Int_t fnCINT1C;           // Number of CINT1C triggers
  Int_t fnCINT1E;           // Number of CINT1E triggers

  Double_t fRecalib[5];     // Correction for abs.calibration per module
  Double_t fBCgap;          // time gap between BC in seconds
  Double_t fVtx0[3];
  Double_t fVtx5[3];

  std::vector<TString> fCuts = {"all", "cpv", "disp", "both"};
  
 // TH2I *fPHOSBadMap[6] ;    //Container for PHOS bad channels map

  TH1F *  fh1Xsec ;                       //! Xsec pythia
  TH1F *  fh1Trials ;                     //! trials pythia

  Bool_t Notify();

  AliPHOSGeometry  *fPHOSGeo;  // PHOS geometry
  Int_t fEventCounter;         // number of analyzed events
  Int_t fAllEventCounter;      // number of analyzed events
  Int_t fInPHOS;               // number of PHOS clusters
  Int_t fInPHOS2;              // number of PHOS clusters
  Int_t fEventCentrality;      // number of analyzed events
  
  Bool_t fEventVtxExist;  
  Bool_t fEventVtxZ10cm;
  Bool_t fEventPileup;
  Bool_t fEventV0AND; 

  AliAODEvent *fEvent;
  TClonesArray *fMCArray;
  AliTriggerAnalysis *fTriggerAnalysis; //! Trigger Analysis for Normalisation

  ClassDef(AliAnalysisTaskGammaPHOSPP, 2); // PHOS analysis task
};

#endif
