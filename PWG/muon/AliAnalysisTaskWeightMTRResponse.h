#ifndef ALIANALYSISTASKWEIGHTMTRRESPONSE_H
#define ALIANALYSISTASKWEIGHTMTRRESPONSE_H

/// \class AliAnalysisTaskWeightMTRResponse
/// \brief Get the Lpt/Apt response from data and MC
///
/// The analyses data and MC and provides the trigger response
///
/// \author Diego Stocco <dstocco@cern.ch>, Subatech
/// \date June 3, 2016


#include "AliAnalysisTaskSE.h"
#include "AliMTRParameterizedResponse.h"
#include "AliMuonTrackCuts.h"
#include "AliUtilityMuonAncestor.h"
#include "AliUtilityDimuonSource.h"

class AliMergeableCollection;
class AliVParticle;

class AliAnalysisTaskWeightMTRResponse : public AliAnalysisTaskSE {
 public:
  class AliTrackMore;
  AliAnalysisTaskWeightMTRResponse();
  AliAnalysisTaskWeightMTRResponse ( const char *name, AliMTRParameterizedResponse* response, Int_t responseType );
  virtual ~AliAnalysisTaskWeightMTRResponse();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void NotifyRun();
  virtual void Terminate(Option_t *option);

  /// Get muon parameterized response
  AliMTRParameterizedResponse* GetMTRParameterizedResponse() { return &fResponse; }

  /// Set muon track cuts
  void SetMTRParameterizedResponse ( const AliMTRParameterizedResponse* response ) { fResponse = *response; }

  /// Get muon track cuts
  AliMuonTrackCuts* GetMuonTrackCuts() { return &fMuonTrackCuts; }

  /// Set muon track cuts
  void SetMuonTrackCuts ( const AliMuonTrackCuts* muonTrackCuts ) { fMuonTrackCuts = *muonTrackCuts; }

 private:

  AliAnalysisTaskWeightMTRResponse(const AliAnalysisTaskWeightMTRResponse&);
  AliAnalysisTaskWeightMTRResponse& operator=(const AliAnalysisTaskWeightMTRResponse&);

  TH2* GetHisto ( const char* source, Bool_t mcResponse, Bool_t perBoard, Bool_t useFit, Int_t makeIt = 0 );

  TString GetSingleMuSource ( AliVParticle* track );

  AliMTRParameterizedResponse fResponse; ///< MTR parameterized response
  AliMuonTrackCuts fMuonTrackCuts;  ///< Muon track cuts
  AliUtilityMuonAncestor fUtilityMuonAncestor; //!<! Muon ancestor finder utility
  AliUtilityDimuonSource fUtilityDimuonSource; //!<! Dimuon source finder utility
  AliMergeableCollection* fMergeableCollection; //!<! collection of mergeable objects
  Int_t fResponseType; ///< response type

  ClassDef(AliAnalysisTaskWeightMTRResponse, 1); // MTR performance analysis
};

class AliAnalysisTaskWeightMTRResponse::AliTrackMore : public TObject
{
public:
  AliTrackMore(AliVParticle* track);
  virtual ~AliTrackMore();

  // void SetPassTrigClassCut ( Int_t itrig ) { fTrigClassCut |= (1<<itrig); }
  // Bool_t MatchTrigPtForTrigClass ( Int_t itrig ) const { return (fTrigClassCut >> itrig) & 0x1; }

  void SetParticleType ( Int_t itype ) { fParticleType = itype; }
  Int_t GetParticleType () const { return fParticleType; }

  void SetAncestor ( Int_t ancestor ) { fAncestor = ancestor; }
  Int_t GetAncestor () const { return fAncestor; }

  void SetHistory ( TString history ) { fHistory = history; }
  TString GetHistory () const { return fHistory; }

  void SetWgt ( Double_t wgt, Bool_t isPerBoard, Bool_t isMC, Bool_t isFit );
  Double_t GetWgt ( Bool_t isPerBoard, Bool_t isMC, Bool_t isFit );

  AliVParticle* GetTrack() { return fTrack; }

private:
  AliVParticle* fTrack; // AliVParticle NOT OWNER
  // UInt_t fTrigClassCut; // Trigger class cut
  Int_t fParticleType; // Particle type
  Int_t fAncestor; // Ancestor index
  Double_t fWgts[8]; // Weight
  TString fHistory; // Track history
};

#endif
