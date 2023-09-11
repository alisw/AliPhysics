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
class AliCaloPhoton;
class AliAnalysisUtils;

#include "AliAnalysisTaskSE.h"

class AliAnalysisSigmaBarCharged : public AliAnalysisTaskSE {
public:
  AliAnalysisSigmaBarCharged(const char *name = "AliAnalysisSigmaBarCharged");
  AliAnalysisSigmaBarCharged(const AliAnalysisSigmaBarCharged &);
  AliAnalysisSigmaBarCharged &operator=(const AliAnalysisSigmaBarCharged &);
  virtual ~AliAnalysisSigmaBarCharged();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);
  // Analyze MC or real data
  void SetMC(Bool_t isMC) { fIsMC = isMC; }
  void SetAdditionHist(Bool_t AdditionHist) { fAdditionHist = AdditionHist; }
  void SetInvMassHist(Bool_t InvMassHist) { fInvMassHist = InvMassHist; }

  void SetClusterTOF(Int_t option) { fPHOSClusterTOFOption = option; }
  void SetMinNbarEnergy(Float_t e) { fCluNbarMinE = e; }
  void SetTOFCut(Float_t tof = 150.e-9) { fCluTimeCut = tof; }
  void SetCPACut(Float_t cpaplus = 0., Float_t cpaminus = 0.) {
    fCPAplusCut = cpaplus;
    fCPAminusCut = cpaminus;
  }
  void SetDCAdaugCut(Float_t dcadaugplus = 0.06, Float_t dcadaugminus = 0.06) {
    fDCAdaugplusCut = dcadaugplus;
    fDCAdaugminusCut = dcadaugminus;
  }
  void SetRADCut(Float_t radplus = 0.25, Float_t radminus = 0.15) {
    fRADplusCut = radplus;
    fRADminusCut = radminus;
  }
  void SetCPVCut(Float_t sigma = 10.) { fCPVCut = sigma; }
  void SetDispCut(Float_t antiph = 4.) { fDispCut = antiph; }
  void SetNcellCut(Int_t Ncell = 7) { fNcellCut = Ncell; }
  void SetNTPCclusters(Int_t TPCclust = 60) { fNTPCclusters = TPCclust; }
  void SetTPCsigmas(Float_t TPCsigmas = 3.) { fTPCsigmas = TPCsigmas; }
  void SetTrackEta(Float_t eta = 0.8) { fTrackEta = eta; }
  void SetTrackBits(Int_t bitmap = 4) { fTracksBits = bitmap; }

protected:
  void SelectSigma();
  Double_t RealRes(Double_t x);
  Double_t RealResv2(Double_t x);
  Bool_t TestDCADaug(Double_t x, Double_t y, Double_t z);
  Bool_t TestRAD(Double_t x, Double_t y, Double_t z);
  void PropagateToDCACurvedBachelor(Double_t *decayparams, AliCaloPhoton *v,
                                    Double_t v0vtx[3], AliExternalTrackParam *t,
                                    Double_t b);
  Float_t GetDCASigmabar(Double_t xyz[3], Double_t pxpypz[3],
                         Double_t bachelor[3]) const;
  void Evaluate(const Double_t *h, Double_t t, Double_t r[3], Double_t g[3],
                Double_t gg[3]);
  Double_t GetCosineOfPointingAngle(Double_t mom[3], Double_t pos[3],
                                    Double_t refPointX, Double_t refPointY,
                                    Double_t refPointZ) const;

  void FillHistogram(const char *key,
                     Double_t x) const; // Fill 1D histogram witn name key
  void FillHistogram(const char *key, Double_t x,
                     Double_t y) const; // Fill 2D histogram witn name key
  void FillHistogram(const char *key, Double_t x, Double_t y,
                     Double_t z) const; // Fill 3D histogram witn name key

private:
  THashList *fOutputContainer;  //! List of output histograms
  AliPIDResponse *fPIDResponse; //! PID response
  AliAnalysisUtils *fUtils;
  TClonesArray *fGamma;       //! List of selected photons
  TList *fPHOSEvents[10][10]; //! Previous events for mixing
  TList *fCurrentMixedList;   //! Previous events for mixing
  Double_t fvtx5[3];
  Double_t fCentrality;          // Centrality
  Int_t fCentBin;                // Centrality
  Float_t fCluTimeCut = 150.e-9; // Time cut
  Float_t fCluNbarMinE = 0.6; // Minimal energy for antineutron selection in GeV
  Float_t fCPAplusCut = 0.;   // Cut on cos of pointing angle AntiSigmaPlus
  Float_t fCPAminusCut = 0.;  // Cut on cos of pointing angle AntiSigmaMinus
  Float_t fDCAdaugplusCut = 0.06;  // Cut on DCA daughters AntiSigmaPlus
  Float_t fDCAdaugminusCut = 0.06; // Cut on DCA daughters AntiSigmaMinus
  Float_t fRADplusCut =
      0.25; // Cut on distance between prim and sec vertexes for AntiSigmaPlus
  Float_t fRADminusCut =
      0.15; // Cut on distance between prim and sec vertexes for AntiSigmaMinus
  Float_t fCPVCut = 10.;    // CPV cut on clusters in sigma
  Float_t fDispCut = 4.;    // Dispersion antiphoton cut for clusters in sigma
  Int_t fNcellCut = 7;      // Dispersion antiphoton cut for clusters in sigma
  Int_t fTracksBits = 4;    // FilterBit
  Float_t fTrackEta = 0.8;  // Pseudurapidity
  Int_t fNTPCclusters = 60; // Number of TPC clusters
  Float_t fTPCsigmas = 3.;  // TPC sigma for pion PID
  Bool_t fIsMC;
  Bool_t fAdditionHist;
  Bool_t fInvMassHist;
  Int_t fPHOSClusterTOFOption;

  TClonesArray *fStack;
  AliAODEvent *fEvent;

  ClassDef(AliAnalysisSigmaBarCharged, 2);
};
#endif
