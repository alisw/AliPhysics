#ifndef AliAnalysisTaskCVEUtil_cxx
#define AliAnalysisTaskCVEUtil_cxx
#include <TH3.h>
#include <memory>
#include "AliAnalysisTaskSE.h"
#include "AliAODTrack.h"
#include "AliEventCuts.h"
#include "AliPIDResponse.h"
#include "AliMultSelection.h"
#include "TList.h"
#include "TProfile.h"

class AliAnalysisTaskCVEUtil : public AliAnalysisTaskSE
{
 public:
  AliAnalysisTaskCVEUtil();
  AliAnalysisTaskCVEUtil(const char* name);
  virtual ~AliAnalysisTaskCVEUtil();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t* option);
  virtual void Terminate(Option_t*);

 private:
  ////////////////////////
  // Procedural function
  ////////////////////////
  bool LoopTracks();

  ////////////////////////
  // Functional function
  ////////////////////////
  // Pile-up
  bool RejectEvtTFFit(float centSPD0);
  // Track
  bool AcceptAODTrack(AliAODTrack* track);
  bool CheckPIDofParticle(AliAODTrack* ftrack, int pidToCheck);
  // Plane
  // Get DCA
  bool GetDCA(float &dcaxy, float &dcaz, AliAODTrack* track);

  //////////////////////
  // Switch           //
  //////////////////////
  bool isNarrowDcaCuts768;
  bool isProtonCustomizedDCACut;
  bool isTightPileUp;

  //////////////////////
  // Cuts and options //
  //////////////////////
  // Global
  TString fTrigger; //
  TString fPeriod;  // period
  // Event
  float fVzCut;      // vz cut
  // Track
  int fFilterBit;          // AOD filter bit selection
  int fNclsCut;            // ncls cut for all tracks
  float fChi2Max;          // upper limmit for chi2
  float fChi2Min;          // lower limmit for chi2
  // PID
  float fNSigmaTPC;
  float fNSigmaRMS;

  ///////////////////The following files are from the data//////////////////////////////////
  /////////////
  // Handles //
  /////////////
  AliAODEvent* fAOD;            // aod Event
  AliPIDResponse* fPIDResponse; // PID Handler
  AliMultSelection* fMultSel;

  ////////////////////////////////
  // Global Variables from data //
  ////////////////////////////////
  std::array<double, 3> fVertex;
  int fRunNum;       // runnumber
  int fOldRunNum;
  int fRunNumBin;    // runnumer bin; 10:139510...; 11:170387...; 15HIR:246994...
  float fCent;

  ///////////////////The following files are read from external sources////////////////////

  ////////////////////////
  // Pile up Function
  ////////////////////////
  std::unique_ptr<TF1> fSPDCutPU;
  std::unique_ptr<TF1> fV0CutPU;
  std::unique_ptr<TF1> fCenCutLowPU;
  std::unique_ptr<TF1> fCenCutHighPU;
  std::unique_ptr<TF1> fMultCutPU;

  ///////////////////The following files will be saved//////////////////////////////////
  //////////////
  // QA Plots //
  //////////////
  TList* fQAList;
  //Proton QA
  TH2D* h2ProtonPtDcaXY;
  TH2D* h2ProtonPtDcaZ;
  TH2D* h2AntiProtonPtDcaXY;
  TH2D* h2AntiProtonPtDcaZ;

  AliAnalysisTaskCVEUtil(const AliAnalysisTaskCVEUtil&);
  AliAnalysisTaskCVEUtil& operator=(const AliAnalysisTaskCVEUtil&);

  ClassDef(AliAnalysisTaskCVEUtil, 1);
};

#endif
