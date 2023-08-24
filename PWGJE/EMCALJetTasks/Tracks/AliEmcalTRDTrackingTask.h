#ifndef AliEmcalTRDTrackingTask_H
#define AliEmcalTRDTrackingTask_H

#if !(defined(__CINT__) || defined(__MAKECINT__))
#include <tuple>
#endif

#include <vector>

#include "AliAnalysisTaskEmcalLight.h"

class THnSparse;

/**
 * @class AliEmcalTRDTrackingTask
 * @brief Tracking QA task
 *
 * Performs tracking QA: efficiency and momentum resolution
 *
 * Based on code in AliAnalysisTaskEMCALClusterize.
 *
 * @author Salvatore Aiola <salvatore.aiola@yale.edu>, Yale University
 * @date Mar 8 2018
 */
class AliEmcalTRDTrackingTask : public AliAnalysisTaskEmcalLight {

 public:
  AliEmcalTRDTrackingTask();
  AliEmcalTRDTrackingTask(const char *name); 
  virtual ~AliEmcalTRDTrackingTask();

  void                   UserCreateOutputObjects();

  static AliEmcalTRDTrackingTask* AddTaskTRDTracking(const char *suffix = "");

 protected:
  Bool_t                 FillHistograms()                               ;
  void                   ExecOnce();
  void                   GenerateHistoBins()                            ;
  void                   AllocateDetectorLevelTHnSparse()               ;
  void                   FillDetectorLevelTHnSparse(Double_t trackPt, Double_t sigma1OverPt, Byte_t trackType, Int_t ntracklets, Double_t trackCharge);
#if !(defined(__CINT__) || defined(__MAKECINT__))
  THnSparse* GenerateTHnSparse(const char* name, const std::vector<std::tuple<std::string, std::vector<Double_t>::iterator, std::vector<Double_t>::iterator>>& axis);
#endif


  // Service fields (non-streamed)
  Bool_t                  fIsEsd                 ; //!<! whether it is ESD data
  AliMCParticleContainer* fGeneratorLevel        ; //!<! generator level container
  AliTrackContainer*      fDetectorLevel         ; //!<! detector level container
  std::vector<Double_t>   fPtHistBins            ; //!<! pt bins
  std::vector<Double_t>   fPtResHistBins         ; //!<! pt res bins
  std::vector<Double_t>   fIntegerHistBins       ; //!<! integer bins
  std::vector<Double_t>   fChargeHistBins;       ; //!<! charge bins

  // Histograms
  THnSparse*              fTracks                ; //! all tracks
  
 private:
  AliEmcalTRDTrackingTask(const AliEmcalTRDTrackingTask&);            // not implemented
  AliEmcalTRDTrackingTask &operator=(const AliEmcalTRDTrackingTask&); // not implemented

  ClassDef(AliEmcalTRDTrackingTask, 5) // Track QA task (efficiency and pt resolution)
};
#endif
