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
  void                   SetSelectHIJING(Bool_t s)  {fSelectHIJING=s;}

 protected:
  Bool_t                 FillHistograms()                               ;
  void                   ExecOnce()                                     ;
  void                   GenerateHistoBins()                            ;
  void                   AllocateDetectorLevelTHnSparse()                    ;
  void                   AllocateGeneratorLevelTHnSparse()                              ;
  void                   AllocateFindableParticlesTHnSparse()           ;
  void                   AllocateMatchedParticlesTHnSparse()            ;
  void                   FillDetectorLevelTHnSparse(Double_t cent, Double_t trackEta, Double_t trackPhi, Double_t trackPt, 
                                                    Double_t sigma1OverPt, Int_t mcGen, Byte_t trackType);
  void                   FillGeneratorLevelTHnSparse(Double_t cent, Double_t partEta, Double_t partPhi, Double_t partPt, Int_t mcGen, Byte_t findable);
  void                   FillMatchedParticlesTHnSparse(Double_t cent, Double_t partEta, Double_t partPhi, Double_t partPt,
						       Double_t trackEta, Double_t trackPhi, Double_t trackPt, Byte_t trackType);

  // Task configuration
  Bool_t                fSelectHIJING          ; //  select HIJING particles

  // Service fields (non-streamed)
  AliParticleContainer* fGeneratorLevel        ; //! generator level container
  AliParticleContainer* fDetectorLevel         ; //! detector level container
  Int_t                 fNPtHistBins           ; //! number of pt bins
  Double_t*             fPtHistBins            ; //! pt bins
  Int_t                 fNEtaHistBins          ; //! number of eta bins
  Double_t*             fEtaHistBins           ; //! eta bins
  Int_t                 fNPhiHistBins          ; //! number of phi bins
  Double_t*             fPhiHistBins           ; //! phi bins
  Int_t                 fNCentHistBins         ; //! number of cent bins
  Double_t*             fCentHistBins          ; //! cent bins
  Int_t                 fNPtResHistBins        ; //! number of pt res bins
  Double_t*             fPtResHistBins         ; //! 1/pt res bins
  Double_t*             f1OverPtResHistBins    ; //! pt res bins
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

  ClassDef(AliEmcalTrackingQATask, 1) // Track QA task (efficiency and pt resolution)
};
#endif
