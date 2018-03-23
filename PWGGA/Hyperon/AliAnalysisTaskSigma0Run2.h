#ifndef AliAnalysisTaskSigma0Run2_H
#define AliAnalysisTaskSigma0Run2_H

// ROOT includes
#include "TChain.h"
#include "TClonesArray.h"
#include "TList.h"
#include "TString.h"

// AliROOT includes
#include "AliAnalysisTaskSE.h"
#include "AliConvEventCuts.h"
#include "AliMCEvent.h"
#include "AliStack.h"

#include "AliSigma0EventCuts.h"
#include "AliSigma0ParticleBase.h"
#include "AliSigma0ParticlePhotonMother.h"
#include "AliSigma0ParticleV0.h"
#include "AliSigma0PhotonCuts.h"
#include "AliSigma0PhotonMotherCuts.h"
#include "AliSigma0SingleParticleCuts.h"
#include "AliSigma0V0Cuts.h"
#include "AliSigma0EventContainer.h"

// forward delcarations
class AliVParticle;
class AliVTrack;

using std::vector;
#include <deque>

class AliAnalysisTaskSigma0Run2 : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskSigma0Run2();
  AliAnalysisTaskSigma0Run2(const char *name);
  virtual ~AliAnalysisTaskSigma0Run2() {}

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);

  AliSigma0SingleParticleCuts *GetSingleParticleCuts() {
    return fSingleParticleCuts;
  }

  void SetV0ReaderName(TString name) { fV0ReaderName = name; }
  void SetIsHeavyIon(bool isHeavyIon) { fIsHeavyIon = isHeavyIon; }
  void SetIsMC(bool isMC) { fIsMC = isMC; }
  void SetIsQA(bool isQA) { fIsQA = isQA; }
  void SetEventCuts(AliSigma0EventCuts *cuts) { fEventCuts = cuts; }
  void SetProtonCuts(AliSigma0SingleParticleCuts *cuts) {
    fSingleParticleCutsProton = cuts;
  }

  void SetSingleParticleCuts(AliSigma0SingleParticleCuts *cuts) {
    fSingleParticleCuts = cuts;
  }
  void SetV0Cuts(AliSigma0V0Cuts *cuts) { fV0Cuts = cuts; }
  void SetV0LambdaCuts(AliSigma0V0Cuts *cuts) { fV0LambdaCuts = cuts; }
  void SetPhotonCuts(AliSigma0PhotonCuts *cuts) { fPhotonCuts = cuts; }
  void SetPhotonMotherCuts(AliSigma0PhotonMotherCuts *cuts) {
    fPhotonMotherCuts = cuts;
  }
  void SetEventContainer(AliSigma0EventContainer *cont) {
    fEventContainer = cont;
  }

 private:
  AliAnalysisTaskSigma0Run2(const AliAnalysisTaskSigma0Run2 &task);
  AliAnalysisTaskSigma0Run2 &operator=(const AliAnalysisTaskSigma0Run2 &task);

  AliVEvent *fInputEvent;    // current event
  AliMCEvent *fMCEvent;      // corresponding MC event
  AliStack *fMCStack;        // stack belonging to MC event
  AliESDpid *fESDpid;        // class for Track PID calculation
  AliV0ReaderV1 *fV0Reader;  //! basic photon Selection Task
  TString fV0ReaderName;     //

  TList *fQA;

  AliSigma0EventCuts *fEventCuts;
  AliSigma0SingleParticleCuts *fSingleParticleCuts;
  AliSigma0SingleParticleCuts *fSingleParticleCutsProton;
  AliSigma0V0Cuts *fV0Cuts;
  AliSigma0V0Cuts *fV0LambdaCuts;
  AliSigma0PhotonCuts *fPhotonCuts;
  AliSigma0PhotonMotherCuts *fPhotonMotherCuts;
  AliSigma0EventContainer *fEventContainer;

  std::vector<AliSigma0ParticleBase> fElectron;  //!
  std::vector<AliSigma0ParticleBase> fPositron;  //!

  std::vector<AliSigma0ParticleBase> fProton;  //!
  std::vector<AliSigma0ParticleBase> fAntiProton;  //!

  std::vector<AliSigma0ParticleV0> fLambda;      //!
  std::vector<AliSigma0ParticleV0> fAntiLambda;  //!

  std::vector<AliSigma0ParticleV0> fLambdaFemto;      //!
  std::vector<AliSigma0ParticleV0> fAntiLambdaFemto;  //!


  std::vector<AliSigma0ParticlePhotonMother> fSigma;      //!
  std::vector<AliSigma0ParticlePhotonMother> fAntiSigma;  //!

  std::vector<AliAODConversionPhoton> fPhoton;  //!

  bool fIsMC;        //
  bool fIsHeavyIon;  //
  bool fIsQA;        //

  TList *fOutputContainer;  //!

  // Histograms
  // =====================================================================

  ClassDef(AliAnalysisTaskSigma0Run2, 1)
};
#endif
