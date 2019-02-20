///
/// \file PWGCF/FEMTOSCOPY/AliFemtoUser/AliAnalysisV0Efficiency.h
///

#ifndef ALIANALYSISV0EFFICIENCY_H
#define ALIANALYSISV0EFFICIENCY_H


#include <iostream>
#include <math.h>
#include <vector>

#include "TChain.h"
#include "TFile.h"
#include "TCollection.h"
#include "TKey.h"
#include "TObject.h"
#include "TObjString.h"
#include "TList.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TROOT.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAODMCParticle.h"
#include "AliAODv0.h"
#include "AliAODRecoDecay.h"
#include "AliAODMCHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliAODpidUtil.h"
#include "AliCentrality.h"
#include "AliESDtrack.h"

using std::vector;

class AliAnalysisV0Efficiency : public AliAnalysisTaskSE {
public:

  enum McOriginType {kUnassigned=0, kPrimary=1, kSig0=2, kXi0=3, kXiC=4, kSigSt0=5, kSigStP=6, kSigStM=7, kKSt0=8, kKStC=9, kOther=10, kFake=11, kMcOriginTypeMax=12};

  AliAnalysisV0Efficiency();
  AliAnalysisV0Efficiency(const char *name, bool aIgnoreInjectedV0s);
  virtual ~AliAnalysisV0Efficiency();

  AliAnalysisV0Efficiency(const AliAnalysisV0Efficiency &aV0Eff);
  AliAnalysisV0Efficiency &operator=(const AliAnalysisV0Efficiency &aV0Eff);


  void MyInit();
  virtual void UserCreateOutputObjects();



  bool IsCorrectEventTrigger();
  static void SetBinLabels(TH1* aHist);
  void FillMCTruthHistogram(int aPID, const AliAODMCParticle* aPart, const AliAODMCParticle* aMother, TH1* aHist, bool aFillFake=false);
  void FillParticleOriginHistogram(int aPID, const AliAODMCParticle* aPart, const AliAODMCParticle* aMother, TH1* aHist);
  void FillHistograms(int aPID, const AliAODMCParticle* aPart, const AliAODMCParticle* aMother, TH1* aMCTruthHist, TH1* aParticleOriginHist, bool aFillFake=false);
  static int GetNumberOfLastHijingLabel(const AliAODEvent *aEvent);
  bool IsInjected(const AliAODMCParticle* aMCv0, TClonesArray *mcArray, int aNumberOfLastHijingLabel);

  bool ContainsSharedDaughters(vector<vector<int> > &aDaughtersCollection, int aPosLabel, int aNegLabel);
  void ExtractOriginalParticles(const AliAODEvent *aEvent);

  double GetNSigmaTOF(const AliAODTrack *aPart, AliPID::EParticleType aParticleType);
  bool IsPionNSigma(const AliAODTrack *aPart);
  bool IsProtonNSigma(const AliAODTrack *aPart);

  bool IsMisIDK0s(int aParticleType, const AliAODv0 *aV0, const AliAODTrack *aPartPos, const AliAODTrack *aPartNeg);
  bool IsMisIDLambda(const AliAODv0 *aV0, const AliAODTrack *aPartPos, const AliAODTrack *aPartNeg);
  bool IsMisIDAntiLambda(const AliAODv0 *aV0, const AliAODTrack *aPartPos, const AliAODTrack *aPartNeg);


  bool V0PassBasicCuts(const AliAODv0* aV0, const AliAODVertex *aAodvertex);
  bool V0PassAllCuts(int aParticleType, const AliAODv0 *aV0, const AliAODTrack *aPartPos, const AliAODTrack *aPartNeg, const AliAODVertex *aAodvertex);

  void ExtractV0FinderParticles(const AliAODEvent *aEvent);
  void ExtractAll(const AliAODEvent *aEvent);

  virtual void Exec(Option_t *option);
  virtual void Terminate(Option_t *);


  TH1F* GetMCTruthOfOriginalParticles_Lam();
  TH1F* GetMCTruthOfOriginalParticles_ALam();
  TH1F* GetMCTruthOfOriginalParticles_K0s();

private:
  AliAODEvent    *fAOD; //! AOD object
  TList          *fOutputList; //! Compact Output list
  AliAODpidUtil  *fpidAOD; //!

protected:

  int fEventCount; 
  bool fIgnoreInjectedV0s;

  TH1F *fMCTruthOfOriginalParticles_Lam; //!
  TH1F *fMCTruthOfV0FinderParticles_Lam; //!
  TH1F *fMCTruthOfReconstructedParticles_Lam; //!

  TH1F *fParticleOriginOfOriginalParticles_Lam; //!
  TH1F *fParticleOriginOfV0FinderParticles_Lam; //!
  TH1F *fParticleOriginOfReconstructedParticles_Lam; //!

  //----------

  TH1F *fMCTruthOfOriginalParticles_ALam; //!
  TH1F *fMCTruthOfV0FinderParticles_ALam; //!
  TH1F *fMCTruthOfReconstructedParticles_ALam; //!

  TH1F *fParticleOriginOfOriginalParticles_ALam; //!
  TH1F *fParticleOriginOfV0FinderParticles_ALam; //!
  TH1F *fParticleOriginOfReconstructedParticles_ALam; //!

  //----------

  TH1F *fMCTruthOfOriginalParticles_K0s; //!
  TH1F *fMCTruthOfV0FinderParticles_K0s; //!
  TH1F *fMCTruthOfReconstructedParticles_K0s; //!

  TH1F *fParticleOriginOfOriginalParticles_K0s; //!
  TH1F *fParticleOriginOfV0FinderParticles_K0s; //!
  TH1F *fParticleOriginOfReconstructedParticles_K0s; //!

  //----------

  TH1F *fReconstructedPurityAid_Lam;
  TH1F *fReconstructedPurityAid_ALam;
  TH1F *fReconstructedPurityAid_K0s;

  //----------

  bool fRemoveMisidentified;

  ClassDef(AliAnalysisV0Efficiency, 1);

};

inline TH1F* AliAnalysisV0Efficiency::GetMCTruthOfOriginalParticles_Lam() {return fMCTruthOfOriginalParticles_Lam;}
inline TH1F* AliAnalysisV0Efficiency::GetMCTruthOfOriginalParticles_ALam() {return fMCTruthOfOriginalParticles_ALam;}
inline TH1F* AliAnalysisV0Efficiency::GetMCTruthOfOriginalParticles_K0s() {return fMCTruthOfOriginalParticles_K0s;}

#endif
