#ifndef ALIANALYSISTASKMTRRESPONSE_H
#define ALIANALYSISTASKMTRRESPONSE_H

/// \class AliAnalysisTaskMTRResponse
/// \brief Get the Lpt/Apt response from data and MC
///
/// The analyses data and MC and provides the trigger response
///
/// \author Diego Stocco <dstocco@cern.ch>, Subatech
/// \date June 3, 2016


#include "AliAnalysisTaskSE.h"
#include "AliMuonEventCuts.h"
#include "AliMuonTrackCuts.h"

class AliMergeableCollection;

class AliAnalysisTaskMTRResponse : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskMTRResponse();
  AliAnalysisTaskMTRResponse ( const char *name );
  virtual ~AliAnalysisTaskMTRResponse();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void NotifyRun();
  virtual void FinishTaskOutput();
  virtual void Terminate(Option_t *option);

  /// Get muon event cuts for muon triggered events
  AliMuonEventCuts* GetMuonEventCutsMuTrig() { return &fMuonEventCutsMuTrig; }
  /// Get muon event cuts for MB events
  AliMuonEventCuts* GetMuonEventCutsMB() { return &fMuonEventCutsMB; }
  /// Get muon track cuts
  AliMuonTrackCuts* GetMuonTrackCuts() { return &fMuonTrackCuts; }

  void SetMuonEventCuts ( AliMuonEventCuts* muonEventCuts, const char* mbClassPattern, const char* muLptClassPattern = 0x0 );
  /// Set muon track cuts
  void SetMuonTrackCuts ( AliMuonTrackCuts* muonTrackCuts ) { fMuonTrackCuts = *muonTrackCuts; }

 private:

  AliAnalysisTaskMTRResponse(const AliAnalysisTaskMTRResponse&);
  AliAnalysisTaskMTRResponse& operator=(const AliAnalysisTaskMTRResponse&);

  TH2* GetHisto ( TString identifier, Int_t imatch, Bool_t perBoard, Bool_t makeIt = kFALSE );

  TObjArray* fMatchTrigKeys;   ///< Match trigger names
  AliMuonEventCuts fMuonEventCutsMuTrig;  ///< Muon event cuts
  AliMuonEventCuts fMuonEventCutsMB;  ///< Muon event cuts
  AliMuonTrackCuts fMuonTrackCuts;  ///< Muon track cuts
  AliMergeableCollection* fMergeableCollection; //!<! collection of mergeable objects
 
  ClassDef(AliAnalysisTaskMTRResponse, 1); // MTR performance analysis
};

#endif
