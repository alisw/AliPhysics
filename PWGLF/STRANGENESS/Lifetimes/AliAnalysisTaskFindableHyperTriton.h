#ifndef AliAnalysisTaskFindableHyperTriton_H
#define AliAnalysisTaskFindableHyperTriton_H

#include <string>
#include <vector>

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliESDtrack.h"

class TTree;

struct Event {
  float fRecVertex[3] = {-999.f, -999.f, -999.f};
  float fTrueVertex[3] = {-999.f, -999.f, -999.f};
  float fMagField = 0.f;
  float fMultiplicity = -1.f;
};

struct FindableHyperTriton {
  FindableHyperTriton() = default;
  std::vector<AliESDtrack> fTracks = {}; //
  std::vector<int> fPDG = {};
  float fDecayVertex[3] = {0.f, 0.f, 0.f};
  float fMomentum[3] = {0.f, 0.f, 0.f};
  float fDeltaT = 0.f;
  unsigned char fFoundTracks = 0u;
  bool fIsPositive = false;
};


class AliAnalysisTaskFindableHyperTriton : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskFindableHyperTriton(std::string name = "FindableHyperTritonTask");
  virtual ~AliAnalysisTaskFindableHyperTriton();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t* option);
  virtual void Terminate(Option_t*) {}

  static AliAnalysisTaskFindableHyperTriton* AddTask(std::string tskname = "FindableHyperTriton", std::string suffix = "");

  AliEventCuts fEventCuts;  /// Event cuts class
  bool fBreakOnMultiBody = true;

 private:
  TTree* fTree  = nullptr;    //! Output Tree

  Event fEventSummary;
  std::vector<FindableHyperTriton> f2Body;
  std::vector<FindableHyperTriton> f3Body;

  AliAnalysisTaskFindableHyperTriton(
      const AliAnalysisTaskFindableHyperTriton&);  // not implemented
  AliAnalysisTaskFindableHyperTriton& operator=(
      const AliAnalysisTaskFindableHyperTriton&);  // not implemented

  ClassDef(AliAnalysisTaskFindableHyperTriton, 1);
};

#endif
