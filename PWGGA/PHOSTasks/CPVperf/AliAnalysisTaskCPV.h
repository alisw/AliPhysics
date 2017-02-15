#ifndef AliAnalysisTaskCPV_cxx
#define AliAnalysisTaskCPV_cxx

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id$ */

// Analysis task for CPV performance via correlation studies
// between CPV clusters and projection of global tracks to CPV
// Author: Sergey Evdokimov

class TObjArray;
class TH1F;
class TH2I;
class TH2F;
class TH3F;
class AliStack ;
class AliESDtrackCuts;
class AliESDtrack ;
class AliESDCaloCluster ;
class AliPHOSGeometry;
class AliTriggerAnalysis;
class AliESDEvent ;
class AliPIDResponse;
class AliPIDCombined;
class AliPHOSCPVGeometry;

#include "TH2I.h"
#include "AliAnalysisTaskSE.h"
#include "AliLog.h"
#include "AliPHOSCPVGeometry.h"

class AliAnalysisTaskCPV : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskCPV(const char *name = "AliAnalysisTaskCPV");
  virtual ~AliAnalysisTaskCPV() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  void SetRCPV(Double_t rCPV) {frCPV = rCPV;}
  Double_t GetRCPV() {return frCPV;}
  
private:
  AliAnalysisTaskCPV(const AliAnalysisTaskCPV&); // not implemented
  AliAnalysisTaskCPV& operator=(const AliAnalysisTaskCPV&); // not implemented
  void FillHistogram(const char * key,Double_t x) const ; //Fill 1D histogram witn name key
  void FillHistogram(const char * key,Double_t x, Double_t y) const ; //Fill 2D histogram witn name key
  void FillHistogram(const char * key,Double_t x, Double_t y, Double_t z) const ; //Fill 3D histogram witn name key
  AliESDtrackCuts* ITSTPCTrackCuts();
  Bool_t IsHotZone(Double_t x, Double_t z);


private:
  TList * fOutputContainer1;    // final histogram container
  TList * fOutputContainer2;    // final histogram container

  AliPHOSGeometry  *fPHOSGeo;  // PHOS geometry
  Int_t fEventCounter;         // number of analyzed events
  AliESDtrackCuts *fITSTPCTrackCuts;
  AliPIDResponse* fPIDResponse;  // PID response object
  AliPIDCombined* fPIDCombined;  // PID combined object
  AliPHOSCPVGeometry* fGeometryCPV;  // CPV geometry object
  Double_t frCPV; // Radial distance from IT to CPV

  ClassDef(AliAnalysisTaskCPV, 1); // PHOS analysis task
};

#endif
