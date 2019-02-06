#ifndef ALITASKMUONTRACKSMEARINGQA_H
#define ALITASKMUONTRACKSMEARINGQA_H

/* $Id$ */ 

///
/// \class AliTaskMuonTrackSmearingQA
/// \brief Task to control the smearing of the muon track parameter according to resolution
///
/// \author Philippe Pillot <pillot@subatech.in2p3.fr>, Subatech
/// \date Nov 8, 2017

#include "AliAnalysisTaskSE.h"
#include "AliMuonTrackSmearing.h"
#include "AliMuonTrackCuts.h"
class TRootIOCtor;

class AliTaskMuonTrackSmearingQA : public AliAnalysisTaskSE {
 public:
  AliTaskMuonTrackSmearingQA(TRootIOCtor* ioCtor);
  AliTaskMuonTrackSmearingQA(const char *name, Int_t chosenFunc);
  virtual ~AliTaskMuonTrackSmearingQA();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *);

  // set cuts to select tracks to be considered
  void SetMuonTrackCuts(AliMuonTrackCuts &trackCuts);
  
  /// get canvas containing histograms for generated tracks
  TCanvas* GetGen() {return fcGen;}
  /// get canvas containing histograms for smeared tracks
  TCanvas* GetRec() {return fcRec;}
  /// get canvas containing histograms for smeared over generated ratios
  TCanvas* GetRat() {return fcRat;}
  /// get canvas containing histograms for track resolution
  TCanvas* GetResVsP() {return fcResVsP;}

 private:
  AliTaskMuonTrackSmearingQA(const AliTaskMuonTrackSmearingQA&);
  AliTaskMuonTrackSmearingQA& operator=(const AliTaskMuonTrackSmearingQA&);

  enum eGenRecList {
    kPt  = 0,
    kP   = 1,
    kEta = 2,
    kY   = 3,
    kPhi = 4
  };
  
  enum eResList {
    kResPAtVtxVsPIn02degMC = 0,
    kResPAtVtxVsPIn23deg   = 1,
    kResPAtVtxVsPIn310deg  = 2,
    kResPtAtVtxVsPt        = 3,
    kResSlopeXAtVtxVsP     = 4,
    kResSlopeYAtVtxVsP     = 5,
    kResEtaAtVtxVsP        = 6,
    kResPhiAtVtxVsP        = 7
  };
  
  TObjArray* fGenList; //!< List of histograms for generated tracks
  TObjArray* fRecList; //!< List of histograms for smeared tracks
  TObjArray* fResList; //!< List of histograms for track resolution

  AliMuonTrackCuts     *fMuonTrackCuts;     ///< cuts to select tracks to be considered
  
  TCanvas *fcGen;    //!< generated tracks
  TCanvas *fcRec;    //!< smeared tracks
  TCanvas *fcRat;    //!< smeared over generated ratios
  TCanvas *fcResVsP; //!< track resolution

  ClassDef(AliTaskMuonTrackSmearingQA, 1);
};

//________________________________________________________________________
inline void AliTaskMuonTrackSmearingQA::SetMuonTrackCuts(AliMuonTrackCuts &trackCuts)
{
  /// set cuts to select tracks to be considered
  delete fMuonTrackCuts;
  fMuonTrackCuts = new AliMuonTrackCuts(trackCuts);
}


#endif
