#ifndef ALIEMCALTRACKINGQATASK_H
#define ALIEMCALTRACKINGQATASK_H

#if !(defined(__CINT__) || defined(__MAKECINT__))
#include <tuple>
#endif

#include <vector>

#include "AliAnalysisTaskEmcalLight.h"

class THnSparse;

/**
 * @class AliEmcalTrackingQATask
 * @brief Tracking QA task
 *
 * Performs tracking QA: efficiency and momentum resolution
 *
 * Based on code in AliAnalysisTaskEMCALClusterize.
 *
 * @author Salvatore Aiola <salvatore.aiola@yale.edu>, Yale University
 * @date Mar 8 2018
 */
class AliEmcalTrackingQATask : public AliAnalysisTaskEmcalLight {

 public:
  AliEmcalTrackingQATask();
  AliEmcalTrackingQATask(const char *name); 
  virtual ~AliEmcalTrackingQATask();

  void                   UserCreateOutputObjects();
  void                   SetDoSigma1OverPt(Bool_t s)      {fDoSigma1OverPt     = s; }
  void                   SetDoSigmaPtOverPtGen(Bool_t s)  {fDoSigmaPtOverPtGen = s; }
  void                   SetDoSeparateTRDrefit(Bool_t s)  {fDoSeparateTRDrefit = s; }
  void                   SetUseTRDUpdateFlag(Bool_t s)    {fUseTRDUpdateFlag   = s; }
  void                   SetQOverPtShift(Double_t s)      {fUseQOverPtShift = true; fQOverPtShift = s; }

  static AliEmcalTrackingQATask* AddTaskTrackingQA(Bool_t isMC, const char *suffix = "");

 protected:
  Bool_t                 FillHistograms()                               ;
  void                   ExecOnce();
  void                   GenerateHistoBins()                            ;
  void                   AllocateDetectorLevelTHnSparse()               ;
  void                   AllocateGeneratorLevelTHnSparse()              ;
  void                   AllocateMatchedParticlesTHnSparse()            ;
  void                   FillDetectorLevelTHnSparse(Double_t cent, Double_t trackEta, Double_t trackPhi, Double_t trackPt, 
                                                    Double_t sigma1OverPt, Int_t mcGen, Byte_t trackType, Double_t trackCharge);
  void                   FillGeneratorLevelTHnSparse(Double_t cent, Double_t partEta, Double_t partPhi, Double_t partPt, Int_t mcGen, Byte_t findable, Double_t partCharge);
  void                   FillMatchedParticlesTHnSparse(Double_t cent, Double_t partEta, Double_t partPhi, Double_t partPt,
                                                       Double_t trackEta, Double_t trackPhi, Double_t trackPt, Byte_t trackType, Double_t trackCharge);
#if !(defined(__CINT__) || defined(__MAKECINT__))
  THnSparse* GenerateTHnSparse(const char* name, const std::vector<std::tuple<std::string, std::vector<Double_t>::iterator, std::vector<Double_t>::iterator>>& axis);
#endif

  // Task configuration
  Bool_t                  fDoSigma1OverPt        ; ///<  add sigma(1/pt), if false add sigma(pt)/pt instead
  Bool_t                  fDoSigmaPtOverPtGen    ; ///<  MC: if true do sigma((ptgen - ptdet) / ptgen), otherwise do sigma((ptgen - ptdet) / ptdet)
  Bool_t                  fDoSeparateTRDrefit    ; ///<  Separate tracks into tracks with TRD refit and 4 tracklets (gold) or not (sub-gold)
  Bool_t                  fUseTRDUpdateFlag      ; ///<  Use TRD update flag to separate tracks into gold and semi-gold
  Bool_t                  fUseQOverPtShift;      ; ///<  Apply Q/pt shift
  Double_t                fQOverPtShift;         ; ///<  Q/pt shift

  // Service fields (non-streamed)
  Bool_t                  fIsEsd                 ; //!<! whether it is ESD data
  AliMCParticleContainer* fGeneratorLevel        ; //!<! generator level container
  AliTrackContainer*      fDetectorLevel         ; //!<! detector level container
  std::vector<Double_t>   fPtHistBins            ; //!<! pt bins
  std::vector<Double_t>   fEtaHistBins           ; //!<! eta bins
  std::vector<Double_t>   fPhiHistBins           ; //!<! phi bins
  std::vector<Double_t>   fCentHistBins          ; //!<! cent bins
  std::vector<Double_t>   fPtRelDiffHistBins     ; //!<! pt relative difference bins
  std::vector<Double_t>   fPtResHistBins         ; //!<! pt res bins
  std::vector<Double_t>   f1OverPtResHistBins    ; //!<! 1/pt res bins
  std::vector<Double_t>   fIntegerHistBins       ; //!<! integer bins
  std::vector<Double_t>   fChargeHistBins;       ; //!<! charge bins

  // Histograms
  THnSparse*              fTracks                ; //! all tracks
  THnSparse*              fParticlesPhysPrim     ; //! all physical primary particles
  THnSparse*              fParticlesMatched      ; //! primary particles matched to detector level tracks
  
 private:
  AliEmcalTrackingQATask(const AliEmcalTrackingQATask&);            // not implemented
  AliEmcalTrackingQATask &operator=(const AliEmcalTrackingQATask&); // not implemented

  ClassDef(AliEmcalTrackingQATask, 5) // Track QA task (efficiency and pt resolution)
};
#endif
