#ifndef ALIEMCALTRACKINGQATASK_H
#define ALIEMCALTRACKINGQATASK_H

#include "AliAnalysisTaskEmcal.h"

class AliParticleContainer;
class THnSparse;
class TH3;

class AliEmcalTrackingQATask : public AliAnalysisTaskEmcal {

 public:
  AliEmcalTrackingQATask();
  AliEmcalTrackingQATask(const char *name); 
  virtual ~AliEmcalTrackingQATask();

  void                   UserCreateOutputObjects();
  void                   SetGeneratorLevelName(const char* name);
  void                   SetDetectorLevelName(const char* name);
  void                   SetDoSigma1OverPt(Bool_t s)      {fDoSigma1OverPt     = s; }
  void                   SetDoSigmaPtOverPtGen(Bool_t s)  {fDoSigmaPtOverPtGen = s; }

 protected:
  Bool_t                 FillHistograms()                               ;
  void                   ExecOnce()                                     ;
  void                   GenerateHistoBins()                            ;
  void                   AllocateDetectorLevelTHnSparse()               ;
  void                   AllocateGeneratorLevelTHnSparse()              ;
  void                   AllocateMatchedParticlesTHnSparse()            ;
  void                   FillDetectorLevelTHnSparse(Double_t cent, Double_t trackEta, Double_t trackPhi, Double_t trackPt, 
                                                    Double_t sigma1OverPt, Int_t mcGen, Byte_t trackType);
  void                   FillGeneratorLevelTHnSparse(Double_t cent, Double_t partEta, Double_t partPhi, Double_t partPt, Int_t mcGen, Byte_t findable);
  void                   FillMatchedParticlesTHnSparse(Double_t cent, Double_t partEta, Double_t partPhi, Double_t partPt,
						       Double_t trackEta, Double_t trackPhi, Double_t trackPt, Byte_t trackType);

  // Task configuration
  Bool_t                fDoSigma1OverPt        ; //  add sigma(1/pt), if false add sigma(pt)/pt instead
  Bool_t                fDoSigmaPtOverPtGen    ; //  MC: if true do sigma((ptgen - ptdet) / ptgen), otherwise do sigma((ptgen - ptdet) / ptdet)

  // Service fields (non-streamed)
  AliMCParticleContainer* fGeneratorLevel      ; //! generator level container
  AliTrackContainer*    fDetectorLevel         ; //! detector level container
  Int_t                 fNPtHistBins           ; //! number of pt bins
  Double_t*             fPtHistBins            ; //! pt bins
  Int_t                 fNEtaHistBins          ; //! number of eta bins
  Double_t*             fEtaHistBins           ; //! eta bins
  Int_t                 fNPhiHistBins          ; //! number of phi bins
  Double_t*             fPhiHistBins           ; //! phi bins
  Int_t                 fNCentHistBins         ; //! number of cent bins
  Double_t*             fCentHistBins          ; //! cent bins
  Int_t                 fNPtRelDiffHistBins    ; //! number of pt relative difference bins
  Double_t*             fPtRelDiffHistBins     ; //! pt relative difference bins
  Int_t                 fNPtResHistBins        ; //! number of pt res bins
  Double_t*             fPtResHistBins         ; //! pt res bins
  Double_t*             f1OverPtResHistBins    ; //! 1/pt res bins
  Int_t                 fN1OverPtResHistBins   ; //! number of 1/pt res bins
  Int_t                 fNIntegerHistBins      ; //! number of integer bins
  Double_t*             fIntegerHistBins       ; //! integer bins


  // Histograms
  THnSparse*            fTracks                ; //! all tracks
  THnSparse*            fParticlesPhysPrim     ; //! all physical primary particles
  THnSparse*            fParticlesMatched      ; //! primary particles matched to detector level tracks
  
 private:
  AliEmcalTrackingQATask(const AliEmcalTrackingQATask&);            // not implemented
  AliEmcalTrackingQATask &operator=(const AliEmcalTrackingQATask&); // not implemented

  ClassDef(AliEmcalTrackingQATask, 3) // Track QA task (efficiency and pt resolution)
};
#endif
