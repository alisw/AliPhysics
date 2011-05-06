// This class contains common elements for the analysis in the First Physics Group.
// based on the code by Arvinder Palaha
// written by Andras Agocs and Anton Alkin

#ifndef ALIANALYSISTASKFIRSTPHYSICS_H
#define ALIANALYSISTASKFIRSTPHYSICS_H

class TH1D;
class TH2D;
class TList;
class AliESDtrackCuts;
class AliESDEvent;

#include "AliAnalysisTaskSE.h"
#include "AliTriggerAnalysis.h"

#include <string>

class AliAnalysisTaskFirstPhysics : public AliAnalysisTaskSE {
 public:
  // indexing various track cuts
  enum {
    kTrackCutQGlo,
    kTrackCutQITS,
    kTrackCutDCAwSPD,
    kTrackCutDCAwoSPD,
    kTrackCutNoSPD,
    kTrackCutTPConly,
    knTrackCuts}; // this must always the last

  // recognised MC process types
  enum ProcessType {
    kProcSD1, // single diffractive AB->XB
    kProcSD2, // single diffractive AB->AX
    kProcDD, // double diffractive
    kProcEL, // elastic
    kProcCD, // central diffractive
    kProcND, // non-diffractive
    kProcIndef};

  AliAnalysisTaskFirstPhysics(const char *name = "You should have given a name to this analysis");
  virtual ~AliAnalysisTaskFirstPhysics();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);

  void SetCutTrackPt(Double_t min, Double_t max) {fCutTrackPtMin = min; fCutTrackPtMax = max;}
  Double_t GetCutTrackPtMin() const {return fCutTrackPtMin;}
  Double_t GetCutTrackPtMax() const {return fCutTrackPtMax;}
  void SetCutEta(Double_t x) {fCutEta = x;}
  Double_t GetCutEta() const {return fCutEta;}
  void SetCutVertexZ(Double_t x) {fCutVertexZ = x;}
  Double_t GetCutVertexZ() const {return fCutVertexZ;}

 protected:
  void PrepareOutputList(); // create fOutput
  void PrepareDefaultTrackCuts(); // create cut objects
  bool PrepareMCInfo(); // check whether MC info is available and read process type
  bool GetESDEvent(); // sets fESD to the current ESD
  bool CheckVertex(); // checks for an appropriate vertex
  bool CheckVertexMC(); // checks for an appropriate vertex in the MC truth

  TH1D* UserHisto1d(const char *name, const char *title, const char *xlabel, Int_t nbinsx, Double_t xlow, Double_t xup);
  TH2D* UserHisto2d(const char *name, const char *title, const char *xlabel, Int_t nbinsx, Double_t xlow, Double_t xup, const char *ylabel, Int_t nbinsy, Double_t ylow, Double_t yup);
  bool GetHisto1FromOutput(const char *name, TH1D *&h) const; // read a histogram from fOutput; use this in Terminate()
  bool GetHisto2FromOutput(const char *name, TH2D *&h) const; // read a 2d histogram from fOutput

  AliESDEvent *fESD; //! the ESD information of the event
  AliMCEvent *fMCEvent; //! the MC information is available
  TList *fOutput; // Output list
  bool fbReadMC; //! indicates if MC information could be read; see PrepareMCInfo()
  ProcessType fMCProcessType; //! indicates the process type used in MC
  AliESDtrackCuts *fTrackCuts[knTrackCuts]; // Track cuts
  AliTriggerAnalysis* fTrigger; //!

 private:
  // simplest cut parameters
  Double_t fCutTrackPtMin;
  Double_t fCutTrackPtMax;
  Double_t fCutEta;
  Double_t fCutVertexZ;

  AliAnalysisTaskFirstPhysics(const AliAnalysisTaskFirstPhysics&); // not implemented
  AliAnalysisTaskFirstPhysics& operator=(const AliAnalysisTaskFirstPhysics&); // not implemented

  ClassDef(AliAnalysisTaskFirstPhysics, 1);
};

#endif

