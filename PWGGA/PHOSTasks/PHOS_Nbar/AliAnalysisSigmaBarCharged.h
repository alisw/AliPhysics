#ifndef AliAnalysisSigmaBarCharged_cxx
#define AliAnalysisSigmaBarCharged_cxx

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

// Analysis task for extracting spectra of resonances in PHOS plus tracks
// Authors: Pavel Gordeev, Dmitri Peresunko
// 01-Feb-2022

class TObjArray;
class TH1F;
class TH2I;
class TH2F;
class TH3F;
class AliESDtrackCuts;
class AliTriggerAnalysis;
class AliPIDResponse;
class AliESDtrack;
class AliAODTrack;
class AliESDVertex;
class AliAODVertex;
class AliExternalTrackParam;
class AliVertexerTracks;
class AliCascadeVertexer;

#include "AliAnalysisTaskSE.h"

class AliAnalysisSigmaBarCharged : public AliAnalysisTaskSE
{
 public:
  AliAnalysisSigmaBarCharged(const char* name = "AliAnalysisSigmaBarCharged");
  AliAnalysisSigmaBarCharged(const AliAnalysisSigmaBarCharged&);
  AliAnalysisSigmaBarCharged& operator=(const AliAnalysisSigmaBarCharged&); 
  virtual ~AliAnalysisSigmaBarCharged();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t* option);
  virtual void Terminate(Option_t*);
  // Analyze MC or real data
  void SetMC(Bool_t isMC) { fIsMC = isMC; }

  void SetClusterTOF(Int_t option) { fPHOSClusterTOFOption = option; }
  void SetMinNbarEnergy(Float_t e) { fCluNbarMinE = e; }
  void SetCPACut(Float_t c = 0.9) { fCPACut = c; }
  void SetCPVCut(Float_t sigma = 4) { fCPVCut = sigma; }
  void SetDispCut(Float_t antiph = 2, Float_t a = -1, Float_t b = 3.5)
  {
    fDispCut = antiph;
    fDispA = a;
    fDispB = b;
  }
  void SetOptimalDepth(Float_t d = 10) { fOptDepth = d; }
  void SetTrackBits(int bitmap = 21) { fTracksBits = bitmap; }

 protected:
  void SelectSigma();
  Double_t RealRes(Double_t x);
  Bool_t TestDCAXY(Double_t x, Double_t y);
  Bool_t TestDCAZ(Double_t x, Double_t y);
  Bool_t DCADaugPlus(Double_t x, Double_t y);
  Bool_t DCADaugMinus(Double_t x, Double_t y);

  void FillHistogram(const char* key, Double_t x) const;                         // Fill 1D histogram witn name key
  void FillHistogram(const char* key, Double_t x, Double_t y) const;             // Fill 2D histogram witn name key
  void FillHistogram(const char* key, Double_t x, Double_t y, Double_t z) const; // Fill 3D histogram witn name key

 private:
  THashList *   fOutputContainer ;  //!List of output histograms
  AliPIDResponse* fPIDResponse; //! PID response
  TClonesArray* fGamma;         //! List of selected photons
  TClonesArray* fPi;            //! List of selected photons
  Int_t fCentBin;               //! current centrality bin
  Int_t fNCenBin;               // Number of centrality bins
  TArrayI fCenBinEdges;         // Centrality binning
  TList* fPHOSEvents[20][20];   //! Previous events for mixing
  TList* fCurrentMixedList;     //! list of previous evetns for given centrality
  Double_t fvtx5[3];
  Float_t fCluMinECut = 0.3;    // Minimal cluster energy cut
  Float_t fCluTimeCut = 25.e-9; // Time cut
  Float_t fCluNbarMinE = 1.;    // Minimal energy for antineutron selection in GeV
  Float_t fOptDepth = 10.;      // Optimal shower depth for nbar momentum reconstruction
  Float_t fCPACut = 0.9;        // Cut on cos of pointing angle
  Float_t fCPVCut = 4;          // CPV cut on clusters in sigma
  Float_t fDispCut = 2;         // Dispersion antiphoton cut for clusters in sigma
  Float_t fDispA = -1;          // Dispersion cut for clusters: M20>=[A]*M02+B
  Float_t fDispB = 3.5;         // Dispersion cut for clusters: M20>=A*M02+[B]
  UInt_t fTracksBits = 21;      // Cut for track->GetFilterMap()
  Bool_t fIsMC;
  Int_t fPHOSClusterTOFOption;
  Int_t fPHOSClusterEnergyOption;

  TClonesArray* fStack;
  AliAODEvent * fEvent;

  ClassDef(AliAnalysisSigmaBarCharged, 1);
};
#endif
