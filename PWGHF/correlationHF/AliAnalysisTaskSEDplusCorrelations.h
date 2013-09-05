#ifndef ALIANALYSISTASKSEDPLUSCORRELATIONS_H
#define ALIANALYSISTASKSEDPLUSCORRELATIONS_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

//*************************************************************************
// Class AliAnalysisTaskSEDplusCorrelations

// Authors:
// Sadhana Dash (correlation)
//*************************************************************************

#include <TROOT.h>
#include <TSystem.h>
#include <TNtuple.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3D.h>
#include <THnSparse.h>
#include <TArrayD.h>

#include "AliRDHFCutsDplustoKpipi.h"
#include "AliHFAssociatedTrackCuts.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisVertexingHF.h"
#include "AliEventPoolManager.h"
#include "AliNormalizationCounter.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliHFCorrelator.h"

class TParticle ;
class TClonesArray ;
class AliAODMCParticle;
class AliAODEvent;
class AliVParticle;
class TObjArray;
class AliEventPoolManager;
class AliESDEvent;



class AliAnalysisTaskSEDplusCorrelations : public AliAnalysisTaskSE
{
 public:

  AliAnalysisTaskSEDplusCorrelations();
  AliAnalysisTaskSEDplusCorrelations(const char *name, AliRDHFCutsDplustoKpipi* analysiscuts, AliHFAssociatedTrackCuts* assotrcuts);
  virtual ~AliAnalysisTaskSEDplusCorrelations();

  void SetReadMC(Bool_t readMC=kTRUE){fReadMC=readMC;}
  void SetEventMix(Bool_t mixing){fMixing=mixing;}
  
  void CreateCorrelationObjs();

  //  standard way used in Dplus task  
  void SetMassLimits(Float_t range);
  void SetMassLimits(Float_t lowlimit, Float_t uplimit);
  void SetBinWidth(Float_t w);

  void SetUseBit(Bool_t dols=kTRUE){fUseBit=dols;}

  void SetCorrelator(Int_t l) {fSelect = l;} // select 1 for hadrons, 2 for Kaons, 3 for Kzeros
  void SetUseDisplacement(Int_t m) {fDisplacement=m;} // select 0 for no displ, 1 for abs displ, 2 for d0/sigma_d0
  void SetSystem(Bool_t system){fSystem=system;} // select between pp (kFALSE) or PbPb (kTRUE)

  void SetUseReconstruction(Bool_t reco){fReco = reco;}

  void SetTrigEfficiency(Bool_t trigeff = kTRUE) {fTrig = trigeff;}

  void FillCorrelations(Double_t ptTrack,Double_t mass, Double_t deltaPhi, Double_t deltaEta, Int_t ind, Int_t sel, Double_t eweight) const;

  void FillMCTruthCorrelations(Double_t ptTrack, Double_t deltaPhi, Double_t deltaEta, Int_t ind, Int_t mcSource, Int_t origDplus, Int_t sel) const;
 
  void FillMCRCorrelations(Double_t ptTrack,Double_t mass, Double_t deltaPhi, Double_t deltaEta, Int_t ind, Int_t mcSource, Int_t origDplus,Double_t eweight) const;
  Bool_t IsDDaughter(AliAODMCParticle* d, AliAODMCParticle* track, TClonesArray* mcArray) const;
  
  Int_t CheckTrackOrigin(TClonesArray* arrayMC, AliAODMCParticle* mcPartCandidate) const ;		
  
  Int_t CheckOrigin(TClonesArray* arrayMC, AliAODMCParticle* mcPartCandidate) const ;		
  



  
  
  Float_t GetUpperMassLimit(){return fUpmasslimit;}
  Float_t GetLowerMassLimit(){return fLowmasslimit;}
  Int_t GetNBinsPt(){return fNPtBins;}
  Float_t GetBinWidth(){return fBinWidth;}
  Int_t GetNBinsHistos();
  
    // Implementation of interface methods
  
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

      
 private:
  AliAnalysisTaskSEDplusCorrelations(const AliAnalysisTaskSEDplusCorrelations &source);
  AliAnalysisTaskSEDplusCorrelations& operator=(const AliAnalysisTaskSEDplusCorrelations& source); 

  

  Int_t GetHistoIndex(Int_t iPtBin) const { return iPtBin*3;}
  
  Int_t GetSignalHistoIndex(Int_t iPtBin) const { return iPtBin*3+1;}
   
 
  enum {kMaxPtBins=20};
  
  //For dplus efficiency
 

  TList  *fOutput; //! list send on output slot 0
 
  AliHFCorrelator* fCorrelator; // object for correlations
 
  Int_t fSelect; // select what to correlate with a Dplus 1-chargedtracks,2-chargedkaons,3-k0s
  Int_t fDisplacement; // set 0 for no displacement cut, 1 for absolute d0, 2 for d0/sigma_d0
  TH1F *fHistNEvents; //!hist. for No. of events
  
  Bool_t fTrig;  // flag for using trig eff         

  TH2D *fEventMix; //!hist. for event mixing
   
  TH1F *fMassHistK0S[3*kMaxPtBins]; //!hist. for inv mass (LC)
  
  TH1F *fLeadPt[3*kMaxPtBins]; //!hist. for D- inv mass (TC)

  TH1F *fPtSig[3*kMaxPtBins]; //!hist. for D- inv mass (TC)
  
  TH1F *fMassHist[3*kMaxPtBins]; //!hist. for inv mass (LC)
  TH1F *fMassHistTC[3*kMaxPtBins]; //!hist. for inv mass (TC)
  TH1F *fMassHistOrigC[3*kMaxPtBins]; //!hist. for inv mass (TC)
  TH1F *fMassHistOrigB[3*kMaxPtBins]; //!hist. for inv mass (TC)
  TH1F *fMassHistMC[3*kMaxPtBins]; //!hist. for inv mass (TC)
 

 
  Float_t fUpmasslimit;  //upper inv mass limit for histos
  Float_t fLowmasslimit; //lower inv mass limit for histos
  
  Int_t fNPtBins; //Number of Pt Bins
  
  Float_t fBinWidth;//width of one bin in output histos

  TList *fListCuts; //list of cuts
  
  //TList *fListCutsAsso; //list of cuts

  AliRDHFCutsDplustoKpipi *fRDCutsAnalysis; //Cuts for Analysis
  
  AliHFAssociatedTrackCuts *fCuts; 
  
  AliNormalizationCounter *fCounter;//!Counter for normalization
  
  Double_t fArrayBinLimits[kMaxPtBins+1]; //limits for the Pt bins
  
  Bool_t fReadMC;    //flag for access to MC
  //  Bool_t fUseStrangeness;//flag to enhance strangeness in MC to fit to data
  Bool_t fUseBit;      // flag to use bitmask
  Bool_t fMixing;      // flag to use bitmask
  
  Bool_t fSystem; //
  Bool_t fReco; // use reconstruction or MC truth
  ClassDef(AliAnalysisTaskSEDplusCorrelations,3); // AliAnalysisTaskSE for the MC association of heavy-flavour decay candidates
};

#endif
