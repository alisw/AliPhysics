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

template<int N>
struct FindableHyperTriton {
  AliESDtrack* fTracks[N] = {nullptr}; //->
  int fPDG[N] = {0};           ///TODO: check if the daughters are sorted according to their mass
  float fDecayVertex[3] = {0.f, 0.f, 0.f};
  float fMomentum[3] = {0.f, 0.f, 0.f};
  float fDeltaT = 0.f;
  unsigned char fFoundTracks = 0u;

  void SetSign(bool pos) { fFoundTracks |= pos ? BIT(7) : 0; }
  bool IsPositive() const { return fFoundTracks & BIT(7); }
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
  std::vector<FindableHyperTriton<2>> f2Body;
  std::vector<FindableHyperTriton<3>> f3Body;

  AliAnalysisTaskFindableHyperTriton(
      const AliAnalysisTaskFindableHyperTriton&);  // not implemented
  AliAnalysisTaskFindableHyperTriton& operator=(
      const AliAnalysisTaskFindableHyperTriton&);  // not implemented

  ClassDef(AliAnalysisTaskFindableHyperTriton, 1);
};

#endif
