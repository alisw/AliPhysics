#ifndef ALIANALYSISTASKPIDFLUCTUATION_H
#define ALIANALYSISTASKPIDFLUCTUATION_H

/* 
 * updated by Zubayer Ahammed to instruct CINT about streaming
 * of data members.
 * Event by event PID fluctuation analysis
 * author: Roberto Preghenella (R+)
 * email:  preghenella@bo.infn.it
 *
 */
 
class TList;
class TH1F;
class TH2F;
class TH3F;
class AliVEvent;
class AliVTrack;
class AliESDtrackCuts;
class AliPIDResponse;

#include "AliAnalysisTaskSE.h"
#include "AliPID.h"
#include "THnSparse.h"

class AliAnalysisTaskPIDFluctuation : 
public AliAnalysisTaskSE 
{
  
 public:
  
  AliAnalysisTaskPIDFluctuation(const Char_t *name = "PIDFluctuation"); // default constructor
  virtual ~AliAnalysisTaskPIDFluctuation(); // default destructor
  
  void UserCreateOutputObjects();
  void UserExec(Option_t *option);

  /* setters */
  void SetPIDMethod(Int_t value) {fPIDMethod = value;}; // set PID method
  void SetESDtrackCuts(AliESDtrackCuts *value) {fESDtrackCuts = value;}; // set ESD track cuts
  void SetAODfilterBit(Int_t value) {fAODfilterBit = value;}; // set AOD filter bit
  void SetEtaRange(Float_t etaMin, Float_t etaMax) {fEtaMin = etaMin; fEtaMax = etaMax;}; // set eta range
  void SetPtRange(Float_t ptMin, Float_t ptMax) {fPtMin = ptMin; fPtMax = ptMax;}; // set pt range

  static const Int_t kNCentralityBins = 10; // N centrality bins

  enum EEventCounter_t {
    kAllEvents,
    kPhysicsSelection,
    kPrimaryVertex,
    kPrimaryVertexSPD,
    kVertexAccepted,
    kGoodCentrality,
    kAcceptedEvents,
    kNEventCounters
  };

  enum ESparseData_t {
    kCent_V0M,
    kCent_TRK,
    kNch,
    kNch_plus,
    kNch_minus,
    kNpi,
    kNpi_plus,
    kNpi_minus,
    kNka,
    kNka_plus,
    kNka_minus,
    kNpr,
    kNpr_plus,
    kNpr_minus,
    kNSparseData
  };

  enum EPIDMethod_t {
    kTPCTOF,
    kTPConly,
    kTOFonly,
    kNPIDMethods
  };

  static void MeasureNuDyn(const Char_t *filename, Int_t i1, Int_t i2, Int_t centralityEstimator = kCent_V0M);

 private:

  AliAnalysisTaskPIDFluctuation(const AliAnalysisTaskPIDFluctuation &); // not implemented
  AliAnalysisTaskPIDFluctuation &operator=(const AliAnalysisTaskPIDFluctuation &); // not implemented
  
  /*** event and track selection ***/
  Bool_t AcceptEvent(AliVEvent *event) const; // accept event
  Bool_t AcceptTrack(AliVTrack *track) const; // accept track

  /*** PID functions ***/
  Bool_t HasTPCPID(AliVTrack *track) const; // has TPC PID
  Bool_t HasTOFPID(AliVTrack *track) const; // has TOF PID
  Double_t MakeTPCPID(AliVTrack *track, Double_t *nSigma) const; // make TPC PID
  Double_t MakeTOFPID(AliVTrack *track, Double_t *nSigma) const; // make TOF PID
  void MakePID(AliVTrack *track, Bool_t *pidFlag, Float_t centrality) const; // make PID
  Bool_t InitPID(AliVEvent *event); // init PID

  /*** PID objects and flags ***/
  Int_t fPIDMethod; // PID method
  AliESDtrackCuts *fESDtrackCuts; // ESD track cuts
  Int_t fAODfilterBit; // AOD filter bit
  Float_t fEtaMin; // eta min
  Float_t fEtaMax; // eta max
  Float_t fPtMin; // pt min
  Float_t fPtMax; // pt max
  AliPIDResponse *fPID; //! PID

  /*** PID histos ***/
  TList *fHistoList; //! histo list
  TH1F *fHistoEventCounter; //! event counter
  TH2F *fHistoAcceptedTracks; //! accepted tracks
  TH2F *fHistoTOFMatchedTracks; //! TOF-matched tracks
  TH3F *fHistoTPCdEdx; //! TPC dEdx
  TH3F *fHistoTPCdEdx_inclusive; //! TPC dEdx
  TH3F *fHistoTOFbeta; //! TOF beta
  TH3F *fHistoTPCdEdx_selected[AliPID::kSPECIES]; //! TPC dEdx
  TH3F *fHistoTOFbeta_selected[AliPID::kSPECIES]; //! TOF beta
  TH3F *fHistoNSigmaTPC[AliPID::kSPECIES]; //! nsigma TPC
  TH3F *fHistoNSigmaTOF[AliPID::kSPECIES]; //! nsigma TOF

  /*** correlation histos */
  THnSparseI *fHistoCorrelation; // correlation THnSparse

  /*** labels, names and titles ***/
  static const Char_t *fgkEventCounterName[kNEventCounters]; // event couter name
  static const Char_t *fgkEventCounterTitle[kNEventCounters]; // event couter title
  static const Char_t *fgkSparseDataName[kNSparseData]; // sparse data name
  static const Char_t *fgkSparseDataTitle[kNSparseData]; // sparse data title

    
  ClassDef(AliAnalysisTaskPIDFluctuation, 1);
};

#endif /* ALIANALYSISTASKPIDFLUCTUATION_H */
