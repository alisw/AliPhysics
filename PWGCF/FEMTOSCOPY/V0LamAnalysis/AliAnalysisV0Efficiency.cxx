///
/// \file AliAnalysisV0Efficiency.cxx
///
#include <cassert>
#include "AliAnalysisV0Efficiency.h"


  /// \cond CLASSIMP
  ClassImp(AliAnalysisV0Efficiency);
  /// \endcond

using namespace std;


//____________________________
AliAnalysisV0Efficiency::AliAnalysisV0Efficiency():
  AliAnalysisTaskSE(),

  fAOD(nullptr),
  fOutputList(nullptr),
  fpidAOD(nullptr),

  fEventCount(0),
  fIgnoreInjectedV0s(false), 

  fMCTruthOfOriginalParticles_Lam(nullptr),
  fMCTruthOfV0FinderParticles_Lam(nullptr),
  fMCTruthOfReconstructedParticles_Lam(nullptr),

  fParticleOriginOfOriginalParticles_Lam(nullptr),
  fParticleOriginOfV0FinderParticles_Lam(nullptr),
  fParticleOriginOfReconstructedParticles_Lam(nullptr),

  fMCTruthOfOriginalParticles_ALam(nullptr),
  fMCTruthOfV0FinderParticles_ALam(nullptr),
  fMCTruthOfReconstructedParticles_ALam(nullptr),

  fParticleOriginOfOriginalParticles_ALam(nullptr),
  fParticleOriginOfV0FinderParticles_ALam(nullptr),
  fParticleOriginOfReconstructedParticles_ALam(nullptr),

  fMCTruthOfOriginalParticles_K0s(nullptr),
  fMCTruthOfV0FinderParticles_K0s(nullptr),
  fMCTruthOfReconstructedParticles_K0s(nullptr),

  fParticleOriginOfOriginalParticles_K0s(nullptr),
  fParticleOriginOfV0FinderParticles_K0s(nullptr),
  fParticleOriginOfReconstructedParticles_K0s(nullptr),

  fReconstructedPurityAid_Lam(nullptr),
  fReconstructedPurityAid_ALam(nullptr),
  fReconstructedPurityAid_K0s(nullptr),

  fRemoveMisidentified(true)
{
}



//____________________________
AliAnalysisV0Efficiency::AliAnalysisV0Efficiency(const char *name, bool aIgnoreInjectedV0s):
  AliAnalysisTaskSE(name),

  fAOD(nullptr),
  fOutputList(nullptr),
  fpidAOD(nullptr),

  fEventCount(0),
  fIgnoreInjectedV0s(aIgnoreInjectedV0s), 

  fMCTruthOfOriginalParticles_Lam(nullptr),
  fMCTruthOfV0FinderParticles_Lam(nullptr),
  fMCTruthOfReconstructedParticles_Lam(nullptr),

  fParticleOriginOfOriginalParticles_Lam(nullptr),
  fParticleOriginOfV0FinderParticles_Lam(nullptr),
  fParticleOriginOfReconstructedParticles_Lam(nullptr),

  fMCTruthOfOriginalParticles_ALam(nullptr),
  fMCTruthOfV0FinderParticles_ALam(nullptr),
  fMCTruthOfReconstructedParticles_ALam(nullptr),

  fParticleOriginOfOriginalParticles_ALam(nullptr),
  fParticleOriginOfV0FinderParticles_ALam(nullptr),
  fParticleOriginOfReconstructedParticles_ALam(nullptr),

  fMCTruthOfOriginalParticles_K0s(nullptr),
  fMCTruthOfV0FinderParticles_K0s(nullptr),
  fMCTruthOfReconstructedParticles_K0s(nullptr),

  fParticleOriginOfOriginalParticles_K0s(nullptr),
  fParticleOriginOfV0FinderParticles_K0s(nullptr),
  fParticleOriginOfReconstructedParticles_K0s(nullptr),

  fReconstructedPurityAid_Lam(nullptr),
  fReconstructedPurityAid_ALam(nullptr),
  fReconstructedPurityAid_K0s(nullptr),

  fRemoveMisidentified(true)
{
  // default constructor
  DefineOutput(1, TList::Class());
}



//________________________________________________________________________
AliAnalysisV0Efficiency::~AliAnalysisV0Efficiency()
{
  if(fOutputList) delete fOutputList;
}


//________________________________________________________________________
AliAnalysisV0Efficiency::AliAnalysisV0Efficiency(const AliAnalysisV0Efficiency &aV0Eff) :
  AliAnalysisTaskSE(aV0Eff),
  fEventCount(aV0Eff.fEventCount),
  fIgnoreInjectedV0s(aV0Eff.fIgnoreInjectedV0s),
  fRemoveMisidentified(aV0Eff.fRemoveMisidentified)
{
  if(aV0Eff.fAOD) fAOD = new AliAODEvent(*(aV0Eff.fAOD));
  else fAOD = nullptr;

  if(aV0Eff.fOutputList) 
  {
    fOutputList = new TList();
    TListIter next(aV0Eff.fOutputList);
    while(TObject* obj = next()) fOutputList->Add(obj);
  }
  else fOutputList = nullptr;

  if(aV0Eff.fpidAOD) fpidAOD = new AliAODpidUtil(*(aV0Eff.fpidAOD));
  else fpidAOD = nullptr;

  //-----

  if(aV0Eff.fMCTruthOfOriginalParticles_Lam) fMCTruthOfOriginalParticles_Lam = new TH1F(*(aV0Eff.fMCTruthOfOriginalParticles_Lam));
  else fMCTruthOfOriginalParticles_Lam = nullptr;

  if(aV0Eff.fMCTruthOfV0FinderParticles_Lam) fMCTruthOfV0FinderParticles_Lam = new TH1F(*(aV0Eff.fMCTruthOfV0FinderParticles_Lam));
  else fMCTruthOfV0FinderParticles_Lam = nullptr;

  if(aV0Eff.fMCTruthOfReconstructedParticles_Lam) fMCTruthOfReconstructedParticles_Lam = new TH1F(*(aV0Eff.fMCTruthOfReconstructedParticles_Lam));
  else fMCTruthOfReconstructedParticles_Lam = nullptr;


  if(aV0Eff.fParticleOriginOfOriginalParticles_Lam) fParticleOriginOfOriginalParticles_Lam = new TH1F(*(aV0Eff.fParticleOriginOfOriginalParticles_Lam));
  else fParticleOriginOfOriginalParticles_Lam = nullptr;

  if(aV0Eff.fParticleOriginOfV0FinderParticles_Lam) fParticleOriginOfV0FinderParticles_Lam = new TH1F(*(aV0Eff.fParticleOriginOfV0FinderParticles_Lam));
  else fParticleOriginOfV0FinderParticles_Lam = nullptr;

  if(aV0Eff.fParticleOriginOfReconstructedParticles_Lam) fParticleOriginOfReconstructedParticles_Lam = new TH1F(*(aV0Eff.fParticleOriginOfReconstructedParticles_Lam));
  else fParticleOriginOfReconstructedParticles_Lam = nullptr;

  //-----

  if(aV0Eff.fMCTruthOfOriginalParticles_ALam) fMCTruthOfOriginalParticles_ALam = new TH1F(*(aV0Eff.fMCTruthOfOriginalParticles_ALam));
  else fMCTruthOfOriginalParticles_ALam = nullptr;

  if(aV0Eff.fMCTruthOfV0FinderParticles_ALam) fMCTruthOfV0FinderParticles_ALam = new TH1F(*(aV0Eff.fMCTruthOfV0FinderParticles_ALam));
  else fMCTruthOfV0FinderParticles_ALam = nullptr;

  if(aV0Eff.fMCTruthOfReconstructedParticles_ALam) fMCTruthOfReconstructedParticles_ALam = new TH1F(*(aV0Eff.fMCTruthOfReconstructedParticles_ALam));
  else fMCTruthOfReconstructedParticles_ALam = nullptr;


  if(aV0Eff.fParticleOriginOfOriginalParticles_ALam) fParticleOriginOfOriginalParticles_ALam = new TH1F(*(aV0Eff.fParticleOriginOfOriginalParticles_ALam));
  else fParticleOriginOfOriginalParticles_ALam = nullptr;

  if(aV0Eff.fParticleOriginOfV0FinderParticles_ALam) fParticleOriginOfV0FinderParticles_ALam = new TH1F(*(aV0Eff.fParticleOriginOfV0FinderParticles_ALam));
  else fParticleOriginOfV0FinderParticles_ALam = nullptr;

  if(aV0Eff.fParticleOriginOfReconstructedParticles_ALam) fParticleOriginOfReconstructedParticles_ALam = new TH1F(*(aV0Eff.fParticleOriginOfReconstructedParticles_ALam));
  else fParticleOriginOfReconstructedParticles_ALam = nullptr;

  //-----

  if(aV0Eff.fMCTruthOfOriginalParticles_K0s) fMCTruthOfOriginalParticles_K0s = new TH1F(*(aV0Eff.fMCTruthOfOriginalParticles_K0s));
  else fMCTruthOfOriginalParticles_K0s = nullptr;

  if(aV0Eff.fMCTruthOfV0FinderParticles_K0s) fMCTruthOfV0FinderParticles_K0s = new TH1F(*(aV0Eff.fMCTruthOfV0FinderParticles_K0s));
  else fMCTruthOfV0FinderParticles_K0s = nullptr;

  if(aV0Eff.fMCTruthOfReconstructedParticles_K0s) fMCTruthOfReconstructedParticles_K0s = new TH1F(*(aV0Eff.fMCTruthOfReconstructedParticles_K0s));
  else fMCTruthOfReconstructedParticles_K0s = nullptr;


  if(aV0Eff.fParticleOriginOfOriginalParticles_K0s) fParticleOriginOfOriginalParticles_K0s = new TH1F(*(aV0Eff.fParticleOriginOfOriginalParticles_K0s));
  else fParticleOriginOfOriginalParticles_K0s = nullptr;

  if(aV0Eff.fParticleOriginOfV0FinderParticles_K0s) fParticleOriginOfV0FinderParticles_K0s = new TH1F(*(aV0Eff.fParticleOriginOfV0FinderParticles_K0s));
  else fParticleOriginOfV0FinderParticles_K0s = nullptr;

  if(aV0Eff.fParticleOriginOfReconstructedParticles_K0s) fParticleOriginOfReconstructedParticles_K0s = new TH1F(*(aV0Eff.fParticleOriginOfReconstructedParticles_K0s));
  else fParticleOriginOfReconstructedParticles_K0s = nullptr;

  //-----

  if(aV0Eff.fReconstructedPurityAid_Lam) fReconstructedPurityAid_Lam = new TH1F(*(aV0Eff.fReconstructedPurityAid_Lam));
  else fReconstructedPurityAid_Lam = nullptr;

  if(aV0Eff.fReconstructedPurityAid_ALam) fReconstructedPurityAid_ALam = new TH1F(*(aV0Eff.fReconstructedPurityAid_ALam));
  else fReconstructedPurityAid_ALam = nullptr;

  if(aV0Eff.fReconstructedPurityAid_K0s) fReconstructedPurityAid_K0s = new TH1F(*(aV0Eff.fReconstructedPurityAid_K0s));
  else fReconstructedPurityAid_K0s = nullptr;

}

//________________________________________________________________________
AliAnalysisV0Efficiency& AliAnalysisV0Efficiency::operator=(const AliAnalysisV0Efficiency &aV0Eff)
{
  //assignment operator
  if (this == &aV0Eff) return *this;

  AliAnalysisTaskSE::operator=(aV0Eff);

  fEventCount = aV0Eff.fEventCount;
  fIgnoreInjectedV0s = aV0Eff.fIgnoreInjectedV0s;
  fRemoveMisidentified = aV0Eff.fRemoveMisidentified;

  //-----

  if(aV0Eff.fAOD) fAOD = new AliAODEvent(*(aV0Eff.fAOD));
  else fAOD = nullptr;

  if(aV0Eff.fOutputList) 
  {
    fOutputList = new TList();
    TListIter next(aV0Eff.fOutputList);
    while(TObject* obj = next()) fOutputList->Add(obj);
  }
  else fOutputList = nullptr;

  if(aV0Eff.fpidAOD) fpidAOD = new AliAODpidUtil(*(aV0Eff.fpidAOD));
  else fpidAOD = nullptr;

  //-----

  if(aV0Eff.fMCTruthOfOriginalParticles_Lam) fMCTruthOfOriginalParticles_Lam = new TH1F(*(aV0Eff.fMCTruthOfOriginalParticles_Lam));
  else fMCTruthOfOriginalParticles_Lam = nullptr;

  if(aV0Eff.fMCTruthOfV0FinderParticles_Lam) fMCTruthOfV0FinderParticles_Lam = new TH1F(*(aV0Eff.fMCTruthOfV0FinderParticles_Lam));
  else fMCTruthOfV0FinderParticles_Lam = nullptr;

  if(aV0Eff.fMCTruthOfReconstructedParticles_Lam) fMCTruthOfReconstructedParticles_Lam = new TH1F(*(aV0Eff.fMCTruthOfReconstructedParticles_Lam));
  else fMCTruthOfReconstructedParticles_Lam = nullptr;


  if(aV0Eff.fParticleOriginOfOriginalParticles_Lam) fParticleOriginOfOriginalParticles_Lam = new TH1F(*(aV0Eff.fParticleOriginOfOriginalParticles_Lam));
  else fParticleOriginOfOriginalParticles_Lam = nullptr;

  if(aV0Eff.fParticleOriginOfV0FinderParticles_Lam) fParticleOriginOfV0FinderParticles_Lam = new TH1F(*(aV0Eff.fParticleOriginOfV0FinderParticles_Lam));
  else fParticleOriginOfV0FinderParticles_Lam = nullptr;

  if(aV0Eff.fParticleOriginOfReconstructedParticles_Lam) fParticleOriginOfReconstructedParticles_Lam = new TH1F(*(aV0Eff.fParticleOriginOfReconstructedParticles_Lam));
  else fParticleOriginOfReconstructedParticles_Lam = nullptr;

  //-----

  if(aV0Eff.fMCTruthOfOriginalParticles_ALam) fMCTruthOfOriginalParticles_ALam = new TH1F(*(aV0Eff.fMCTruthOfOriginalParticles_ALam));
  else fMCTruthOfOriginalParticles_ALam = nullptr;

  if(aV0Eff.fMCTruthOfV0FinderParticles_ALam) fMCTruthOfV0FinderParticles_ALam = new TH1F(*(aV0Eff.fMCTruthOfV0FinderParticles_ALam));
  else fMCTruthOfV0FinderParticles_ALam = nullptr;

  if(aV0Eff.fMCTruthOfReconstructedParticles_ALam) fMCTruthOfReconstructedParticles_ALam = new TH1F(*(aV0Eff.fMCTruthOfReconstructedParticles_ALam));
  else fMCTruthOfReconstructedParticles_ALam = nullptr;


  if(aV0Eff.fParticleOriginOfOriginalParticles_ALam) fParticleOriginOfOriginalParticles_ALam = new TH1F(*(aV0Eff.fParticleOriginOfOriginalParticles_ALam));
  else fParticleOriginOfOriginalParticles_ALam = nullptr;

  if(aV0Eff.fParticleOriginOfV0FinderParticles_ALam) fParticleOriginOfV0FinderParticles_ALam = new TH1F(*(aV0Eff.fParticleOriginOfV0FinderParticles_ALam));
  else fParticleOriginOfV0FinderParticles_ALam = nullptr;

  if(aV0Eff.fParticleOriginOfReconstructedParticles_ALam) fParticleOriginOfReconstructedParticles_ALam = new TH1F(*(aV0Eff.fParticleOriginOfReconstructedParticles_ALam));
  else fParticleOriginOfReconstructedParticles_ALam = nullptr;

  //-----

  if(aV0Eff.fMCTruthOfOriginalParticles_K0s) fMCTruthOfOriginalParticles_K0s = new TH1F(*(aV0Eff.fMCTruthOfOriginalParticles_K0s));
  else fMCTruthOfOriginalParticles_K0s = nullptr;

  if(aV0Eff.fMCTruthOfV0FinderParticles_K0s) fMCTruthOfV0FinderParticles_K0s = new TH1F(*(aV0Eff.fMCTruthOfV0FinderParticles_K0s));
  else fMCTruthOfV0FinderParticles_K0s = nullptr;

  if(aV0Eff.fMCTruthOfReconstructedParticles_K0s) fMCTruthOfReconstructedParticles_K0s = new TH1F(*(aV0Eff.fMCTruthOfReconstructedParticles_K0s));
  else fMCTruthOfReconstructedParticles_K0s = nullptr;


  if(aV0Eff.fParticleOriginOfOriginalParticles_K0s) fParticleOriginOfOriginalParticles_K0s = new TH1F(*(aV0Eff.fParticleOriginOfOriginalParticles_K0s));
  else fParticleOriginOfOriginalParticles_K0s = nullptr;

  if(aV0Eff.fParticleOriginOfV0FinderParticles_K0s) fParticleOriginOfV0FinderParticles_K0s = new TH1F(*(aV0Eff.fParticleOriginOfV0FinderParticles_K0s));
  else fParticleOriginOfV0FinderParticles_K0s = nullptr;

  if(aV0Eff.fParticleOriginOfReconstructedParticles_K0s) fParticleOriginOfReconstructedParticles_K0s = new TH1F(*(aV0Eff.fParticleOriginOfReconstructedParticles_K0s));
  else fParticleOriginOfReconstructedParticles_K0s = nullptr;

  //-----

  if(aV0Eff.fReconstructedPurityAid_Lam) fReconstructedPurityAid_Lam = new TH1F(*(aV0Eff.fReconstructedPurityAid_Lam));
  else fReconstructedPurityAid_Lam = nullptr;

  if(aV0Eff.fReconstructedPurityAid_ALam) fReconstructedPurityAid_ALam = new TH1F(*(aV0Eff.fReconstructedPurityAid_ALam));
  else fReconstructedPurityAid_ALam = nullptr;

  if(aV0Eff.fReconstructedPurityAid_K0s) fReconstructedPurityAid_K0s = new TH1F(*(aV0Eff.fReconstructedPurityAid_K0s));
  else fReconstructedPurityAid_K0s = nullptr;

  return *this;
}


//________________________________________________________________________
void AliAnalysisV0Efficiency::MyInit()
{
  AliAODInputHandler *aodH = dynamic_cast<AliAODInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  fpidAOD = aodH->GetAODpidUtil();
}

//________________________________________________________________________
void AliAnalysisV0Efficiency::UserCreateOutputObjects()
{
  // Create output histograms.
  // Histograms are added to fOutputList
  // When fOutputList is deleted, it automatically cleans up all
  // associated histograms
  // cout<<"Create histograms"<<endl;
  fOutputList = new TList();
  fOutputList->SetOwner();
  MyInit();// Initialize my settings


  fMCTruthOfOriginalParticles_Lam = new TH1F("fMCTruthOfOriginalParticles_Lam", "fMCTruthOfOriginalParticles_Lam",
                                         kMcOriginTypeMax+1, 0, kMcOriginTypeMax+1);
  fMCTruthOfV0FinderParticles_Lam = new TH1F("fMCTruthOfV0FinderParticles_Lam", "fMCTruthOfV0FinderParticles_Lam",
                                         kMcOriginTypeMax+1, 0, kMcOriginTypeMax+1);
  fMCTruthOfReconstructedParticles_Lam = new TH1F("fMCTruthOfReconstructedParticles_Lam", "fMCTruthOfReconstructedParticles_Lam",
                                         kMcOriginTypeMax+1, 0, kMcOriginTypeMax+1);


  fParticleOriginOfOriginalParticles_Lam = new TH1F("fParticleOriginOfOriginalParticles_Lam", "fParticleOriginOfOriginalParticles_Lam",
                                         6000, 0.0, 6000.0);
  fParticleOriginOfV0FinderParticles_Lam = new TH1F("fParticleOriginOfV0FinderParticles_Lam", "fParticleOriginOfV0FinderParticles_Lam",
                                         6000, 0.0, 6000.0);
  fParticleOriginOfReconstructedParticles_Lam = new TH1F("fParticleOriginOfReconstructedParticles_Lam", "fParticleOriginOfReconstructedParticles_Lam",
                                         6000, 0.0, 6000.0);

  //----------

  fMCTruthOfOriginalParticles_ALam = new TH1F("fMCTruthOfOriginalParticles_ALam", "fMCTruthOfOriginalParticles_ALam",
                                         kMcOriginTypeMax+1, 0, kMcOriginTypeMax+1);
  fMCTruthOfV0FinderParticles_ALam = new TH1F("fMCTruthOfV0FinderParticles_ALam", "fMCTruthOfV0FinderParticles_ALam",
                                         kMcOriginTypeMax+1, 0, kMcOriginTypeMax+1);
  fMCTruthOfReconstructedParticles_ALam = new TH1F("fMCTruthOfReconstructedParticles_ALam", "fMCTruthOfReconstructedParticles_ALam",
                                         kMcOriginTypeMax+1, 0, kMcOriginTypeMax+1);


  fParticleOriginOfOriginalParticles_ALam = new TH1F("fParticleOriginOfOriginalParticles_ALam", "fParticleOriginOfOriginalParticles_ALam",
                                         6000, 0.0, 6000.0);
  fParticleOriginOfV0FinderParticles_ALam = new TH1F("fParticleOriginOfV0FinderParticles_ALam", "fParticleOriginOfV0FinderParticles_ALam",
                                         6000, 0.0, 6000.0);
  fParticleOriginOfReconstructedParticles_ALam = new TH1F("fParticleOriginOfReconstructedParticles_ALam", "fParticleOriginOfReconstructedParticles_ALam",
                                         6000, 0.0, 6000.0);

  //----------

  fMCTruthOfOriginalParticles_K0s = new TH1F("fMCTruthOfOriginalParticles_K0s", "fMCTruthOfOriginalParticles_K0s",
                                         kMcOriginTypeMax+1, 0, kMcOriginTypeMax+1);
  fMCTruthOfV0FinderParticles_K0s = new TH1F("fMCTruthOfV0FinderParticles_K0s", "fMCTruthOfV0FinderParticles_K0s",
                                         kMcOriginTypeMax+1, 0, kMcOriginTypeMax+1);
  fMCTruthOfReconstructedParticles_K0s = new TH1F("fMCTruthOfReconstructedParticles_K0s", "fMCTruthOfReconstructedParticles_K0s",
                                         kMcOriginTypeMax+1, 0, kMcOriginTypeMax+1);


  fParticleOriginOfOriginalParticles_K0s = new TH1F("fParticleOriginOfOriginalParticles_K0s", "fParticleOriginOfOriginalParticles_K0s",
                                         6000, 0.0, 6000.0);
  fParticleOriginOfV0FinderParticles_K0s = new TH1F("fParticleOriginOfV0FinderParticles_K0s", "fParticleOriginOfV0FinderParticles_K0s",
                                         6000, 0.0, 6000.0);
  fParticleOriginOfReconstructedParticles_K0s = new TH1F("fParticleOriginOfReconstructedParticles_K0s", "fParticleOriginOfReconstructedParticles_K0s",
                                         6000, 0.0, 6000.0);


  //----------

  SetBinLabels(fMCTruthOfOriginalParticles_Lam);
  SetBinLabels(fMCTruthOfV0FinderParticles_Lam);
  SetBinLabels(fMCTruthOfReconstructedParticles_Lam);

  SetBinLabels(fMCTruthOfOriginalParticles_ALam);
  SetBinLabels(fMCTruthOfV0FinderParticles_ALam);
  SetBinLabels(fMCTruthOfReconstructedParticles_ALam);

  SetBinLabels(fMCTruthOfOriginalParticles_K0s);
  SetBinLabels(fMCTruthOfV0FinderParticles_K0s);
  SetBinLabels(fMCTruthOfReconstructedParticles_K0s);

  //-----
/*
  fOutputList->Add(fMCTruthOfOriginalParticles_Lam);
  fOutputList->Add(fMCTruthOfV0FinderParticles_Lam);
  fOutputList->Add(fMCTruthOfReconstructedParticles_Lam);

  fOutputList->Add(fMCTruthOfOriginalParticles_ALam);
  fOutputList->Add(fMCTruthOfV0FinderParticles_ALam);
  fOutputList->Add(fMCTruthOfReconstructedParticles_ALam);

  fOutputList->Add(fMCTruthOfOriginalParticles_K0s);
  fOutputList->Add(fMCTruthOfV0FinderParticles_K0s);
  fOutputList->Add(fMCTruthOfReconstructedParticles_K0s);
*/

  fOutputList->Add(fMCTruthOfOriginalParticles_Lam);
  fOutputList->Add(fMCTruthOfOriginalParticles_ALam);
  fOutputList->Add(fMCTruthOfOriginalParticles_K0s);

  fOutputList->Add(fMCTruthOfV0FinderParticles_Lam);
  fOutputList->Add(fMCTruthOfV0FinderParticles_ALam);
  fOutputList->Add(fMCTruthOfV0FinderParticles_K0s);

  fOutputList->Add(fMCTruthOfReconstructedParticles_Lam);
  fOutputList->Add(fMCTruthOfReconstructedParticles_ALam);
  fOutputList->Add(fMCTruthOfReconstructedParticles_K0s);

  //----------

  fOutputList->Add(fParticleOriginOfOriginalParticles_Lam);
  fOutputList->Add(fParticleOriginOfOriginalParticles_ALam);
  fOutputList->Add(fParticleOriginOfOriginalParticles_K0s);

  fOutputList->Add(fParticleOriginOfV0FinderParticles_Lam);
  fOutputList->Add(fParticleOriginOfV0FinderParticles_ALam);
  fOutputList->Add(fParticleOriginOfV0FinderParticles_K0s);

  fOutputList->Add(fParticleOriginOfReconstructedParticles_Lam);
  fOutputList->Add(fParticleOriginOfReconstructedParticles_ALam);
  fOutputList->Add(fParticleOriginOfReconstructedParticles_K0s);

  //----------------------------------------------------

  fReconstructedPurityAid_Lam = new TH1F("fReconstructedPurityAid_Lam", "fReconstructedPurityAid_Lam",
                                         100, 1.080683, 1.150683);

  fReconstructedPurityAid_ALam = new TH1F("fReconstructedPurityAid_ALam", "fReconstructedPurityAid_ALam",
                                         100, 1.080683, 1.150683);

  fReconstructedPurityAid_K0s = new TH1F("fReconstructedPurityAid_K0s", "fReconstructedPurityAid_K0s",
                                         100, 0.427614, 0.567614);

  fOutputList->Add(fReconstructedPurityAid_Lam);
  fOutputList->Add(fReconstructedPurityAid_ALam);
  fOutputList->Add(fReconstructedPurityAid_K0s);

  //----------------------------------------------------
  PostData(1, fOutputList);
}



//________________________________________________________________________s
bool AliAnalysisV0Efficiency::IsCorrectEventTrigger()
{
  //Pick out Central, SemiCentral, and MinBias events.  False if not using one of those event triggers.
  Bool_t isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & (AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral));
  return isSelected;
}


//____________________________
void AliAnalysisV0Efficiency::SetBinLabels(TH1* aHist)
{
  aHist->GetXaxis()->SetBinLabel(1, "Unassigned");
  aHist->GetXaxis()->SetBinLabel(2, "Primary");
  aHist->GetXaxis()->SetBinLabel(3, "#Sigma^{0}");
  aHist->GetXaxis()->SetBinLabel(4, "#Xi^{0}");
  aHist->GetXaxis()->SetBinLabel(5, "#Xi^{-}");
  aHist->GetXaxis()->SetBinLabel(6, "#Sigma^{*0}");
  aHist->GetXaxis()->SetBinLabel(7, "#Sigma^{*+}");
  aHist->GetXaxis()->SetBinLabel(8, "#Sigma^{*-}");
  aHist->GetXaxis()->SetBinLabel(9, "K^{*0}");
  aHist->GetXaxis()->SetBinLabel(10, "K^{*ch}");
  aHist->GetXaxis()->SetBinLabel(11, "Other");
  aHist->GetXaxis()->SetBinLabel(12, "Fake");
  aHist->GetXaxis()->SetBinLabel(13, "W.E.");

}

//____________________________
void AliAnalysisV0Efficiency::FillMCTruthHistogram(int aPID, const AliAODMCParticle* aPart, const AliAODMCParticle* aMother, TH1* aHist, bool aFillFake)
{
  //TODO it appears a particle can have a mother, but still be considered primary
  //     i.e. calling aPart->IsPhysicalPrimary() can still be true
  //     This appears to be the case when the mother decays via EM or strong

  if(aPart->GetPdgCode()==aPID)
  {
    if(aPart->GetMother() <= 0) aHist->Fill(kPrimary);
    else
    {
      int tMotherPDG = fabs(aMother->GetPdgCode());  //Note the fabs!

      //It appears K0Short(=310) will have a parent K0(=311)
      if(aPID==310 && tMotherPDG==311) aHist->Fill(kPrimary);

      else if(aMother->GetCalcMass() < aPart->GetCalcMass()) aHist->Fill(kUnassigned);  //Mother should be heavier than daughter!
      else if(tMotherPDG==3212) aHist->Fill(kSig0);
      else if(tMotherPDG==3322) aHist->Fill(kXi0);
      else if(tMotherPDG==3312) aHist->Fill(kXiC);
      else if(tMotherPDG==3214) aHist->Fill(kSigSt0);
      else if(tMotherPDG==3224) aHist->Fill(kSigStP);
      else if(tMotherPDG==3114) aHist->Fill(kSigStM);
      else if(tMotherPDG==313) aHist->Fill(kKSt0);
      else if(tMotherPDG==323) aHist->Fill(kKStC);
      else aHist->Fill(kOther);
    }
  }
  else if(aFillFake) aHist->Fill(kFake);
}


//____________________________
void AliAnalysisV0Efficiency::FillParticleOriginHistogram(int aPID, const AliAODMCParticle* aPart, const AliAODMCParticle* aMother, TH1* aHist)
{
  if(aPart->GetPdgCode()==aPID)
  {
    if(aPart->GetMother() > -1) aHist->Fill(TMath::Abs(aMother->GetPdgCode()));
    else aHist->Fill(0.);
  }
}

//____________________________
void AliAnalysisV0Efficiency::FillHistograms(int aPID, const AliAODMCParticle* aPart, const AliAODMCParticle* aMother, TH1* aMCTruthHist, TH1* aParticleOriginHist, bool aFillFake)
{
  FillMCTruthHistogram(aPID, aPart, aMother, aMCTruthHist, aFillFake);
  FillParticleOriginHistogram(aPID, aPart, aMother, aParticleOriginHist);
}


//____________________________
int AliAnalysisV0Efficiency::GetNumberOfLastHijingLabel(const AliAODEvent *aEvent)
{
  int tNumberOfLastHijingLabel = -1;
  AliAODMCHeader *mcHeader = (AliAODMCHeader*)aEvent->FindListObject(AliAODMCHeader::StdBranchName());
  if(mcHeader){
    // Get the iterator on the list of cocktail headers
    TIter next(mcHeader->GetCocktailHeaders());
    // Loop over the cocktail headers
    while (const TObject *obj=next()){
      // Check whether it's a Hijing header
      const AliGenHijingEventHeader* hijingHeader = dynamic_cast<const AliGenHijingEventHeader*>(obj);
      if(hijingHeader) {
	tNumberOfLastHijingLabel=hijingHeader->NProduced()-1;
      } // End of found the hijing header
    } // End of loop over cocktail headers
  } // End of MC header exists

  return tNumberOfLastHijingLabel;
}

//____________________________
bool AliAnalysisV0Efficiency::IsInjected(const AliAODMCParticle* aMCv0, TClonesArray *mcArray, int aNumberOfLastHijingLabel)
{
  //Reject injected particles, if desired.  
  //Injected particles have a label greater than tNumberOfLastHijingLabel.  Many secondary particles also have a label greater than tNumberOfLastHijingLabel.  
  //So first check if a particle has a parent.  If it does, check if that parent is "original".  If not original, reject that particle.  
  //If a particle has no parent and it has a label greater than numberOfLastHijingLabel, reject it

  bool aIsInjected = false;
  if(aMCv0->GetMother() > -1) //V0 has a mother
  {
    AliAODMCParticle *tMotherOfV0 = (AliAODMCParticle*)(mcArray->At(aMCv0->GetMother()));
    if(tMotherOfV0->GetLabel() > aNumberOfLastHijingLabel) aIsInjected=true;
  }
  else //V0 has no mother
  {
    if(aMCv0->GetLabel() > aNumberOfLastHijingLabel) aIsInjected=true;
  }

  return aIsInjected;
}

//____________________________
bool AliAnalysisV0Efficiency::ContainsSharedDaughters(vector<vector<int> > &aDaughtersCollection, int aPosLabel, int aNegLabel)
{
  bool aPosDuplicate=false;
  bool aNegDuplicate=false;

  for(int i=0; i<(int)aDaughtersCollection.size(); i++)
  {
    if(aDaughtersCollection[i][0]==aPosLabel) aPosDuplicate=true;
    if(aDaughtersCollection[i][1]==aNegLabel) aNegDuplicate=true;
  }

  if(!aPosDuplicate && !aNegDuplicate) return false;
  else return true;
}

//____________________________
void AliAnalysisV0Efficiency::ExtractOriginalParticles(const AliAODEvent *aEvent)
{
  TClonesArray *mcP = NULL;
  mcP = (TClonesArray *) aEvent->FindListObject(AliAODMCParticle::StdBranchName());
  if (!mcP) 
  {
    cout << "AOD MC information requested, but no particle array found!" << endl;
    return;
  }

  int tNumberOfLastHijingLabel = GetNumberOfLastHijingLabel(aEvent);
  if(tNumberOfLastHijingLabel <= 0) return;

  vector<vector<int> > tDaughterLabelsLam(0), tDaughterLabelsALam(0), tDaughterLabelsK0s(0);

  for(int i=0; i < mcP->GetEntriesFast(); i++)
  {

    const AliAODMCParticle *tPart = (AliAODMCParticle*)mcP->At(i);
    if(!tPart) continue;
    if(tPart->GetNDaughters() != 2) continue;
//    int tPartPID = tPart->GetPdgCode();

    bool tIsInjected = IsInjected(tPart, mcP, tNumberOfLastHijingLabel);
    AliAODMCParticle *tMother = NULL;
    if(tPart->GetMother() > -1)  //MC particle has a mother 
    {
      tMother = (AliAODMCParticle *)mcP->At(tPart->GetMother());
    }



    const AliAODMCParticle *mcParticlePos = (AliAODMCParticle*)(mcP->At(tPart->GetDaughterLabel(0))),
                           *mcParticleNeg = (AliAODMCParticle*)(mcP->At(tPart->GetDaughterLabel(1)));

    // The daughter info must exist for both
    if(mcParticlePos==NULL || mcParticleNeg==NULL) continue;

    if (mcParticlePos->Charge() == mcParticleNeg->Charge()) continue;

    // ensure that pos and neg are pointing to the correct children
    if (mcParticlePos->Charge() < 0 && mcParticleNeg->Charge() > 0) 
    {
      const AliAODMCParticle *tmp = mcParticlePos;
      mcParticlePos = mcParticleNeg;
      mcParticleNeg = tmp;
    }

    const int motherOfPosID = mcParticlePos->GetMother(),
              motherOfNegID = mcParticleNeg->GetMother();

    // If both daughter tracks refer to the same mother, we can continue
    if(motherOfPosID < 0 || motherOfPosID != motherOfNegID) continue;


    if(fabs(mcParticlePos->Eta()) > 0.8) continue;
    if(fabs(mcParticleNeg->Eta()) > 0.8) continue;

    if(mcParticlePos->Pt() < 0.16) continue;
    if(mcParticleNeg->Pt() < 0.16) continue;

    if(mcParticlePos->GetPdgCode() == 2212 && mcParticlePos->Pt() < 0.5) continue;
    if(mcParticlePos->GetPdgCode() == -2212 && mcParticlePos->Pt() < 0.3) continue;

    if(mcParticleNeg->GetPdgCode() == 2212 && mcParticleNeg->Pt() < 0.5) continue;
    if(mcParticleNeg->GetPdgCode() == -2212 && mcParticleNeg->Pt() < 0.3) continue;

    if(fIgnoreInjectedV0s && tIsInjected) continue;

    if(!ContainsSharedDaughters(tDaughterLabelsLam, mcParticlePos->GetLabel(), mcParticleNeg->GetLabel()))
    {
      FillHistograms(3122, tPart, tMother, fMCTruthOfOriginalParticles_Lam, fParticleOriginOfOriginalParticles_Lam); 
    }
    if(!ContainsSharedDaughters(tDaughterLabelsALam, mcParticlePos->GetLabel(), mcParticleNeg->GetLabel())) 
    {
      FillHistograms(-3122, tPart, tMother, fMCTruthOfOriginalParticles_ALam, fParticleOriginOfOriginalParticles_ALam); 
    }
    if(!ContainsSharedDaughters(tDaughterLabelsK0s, mcParticlePos->GetLabel(), mcParticleNeg->GetLabel()))  
    {
      FillHistograms(310, tPart, tMother, fMCTruthOfOriginalParticles_K0s, fParticleOriginOfOriginalParticles_K0s); 
    }

    if(tPart->GetPdgCode()==3122)  tDaughterLabelsLam.push_back(vector<int>{mcParticlePos->GetLabel(), mcParticleNeg->GetLabel()});
    if(tPart->GetPdgCode()==-3122) tDaughterLabelsALam.push_back(vector<int>{mcParticlePos->GetLabel(), mcParticleNeg->GetLabel()});
    if(tPart->GetPdgCode()==310)   tDaughterLabelsK0s.push_back(vector<int>{mcParticlePos->GetLabel(), mcParticleNeg->GetLabel()});
  }

}

//____________________________
double AliAnalysisV0Efficiency::GetNSigmaTOF(const AliAODTrack *aPart, AliPID::EParticleType aParticleType)
{
  double tNSigmaTOF = -1234;
  double tProbMis = 1.0;
  if(((aPart->GetStatus() & AliVTrack::kTOFout) == AliVTrack::kTOFout) && ((aPart->GetStatus() & AliVTrack::kTIME) == AliVTrack::kTIME)) 
  {
    tProbMis = fpidAOD->GetTOFMismatchProbability(aPart);
  }

  if(!(((aPart->GetStatus() & AliVTrack::kTOFout) == AliVTrack::kTOFout) && ((aPart->GetStatus() & AliVTrack::kTIME) == AliVTrack::kTIME)) || tProbMis > 0.01) 
  {
    tNSigmaTOF = -1000;
  }
  else
  {
    if(((aPart->GetStatus() & AliVTrack::kTOFout) == AliVTrack::kTOFout) && ((aPart->GetStatus() & AliVTrack::kTIME) == AliVTrack::kTIME) && tProbMis < 0.01) 
    {
      tNSigmaTOF = fpidAOD->NumberOfSigmasTOF(aPart, aParticleType);
    }
  }
  return tNSigmaTOF;
}


//____________________________
bool AliAnalysisV0Efficiency::IsPionNSigma(const AliAODTrack *aPart)
{
  double tNSigmaTPCPi = fpidAOD->NumberOfSigmasTPC(aPart, AliPID::kPion);
//  double tNSigmaTOFPi = fpidAOD->NumberOfSigmasTOF(aPart, AliPID::kPion);
  double tNSigmaTOFPi = GetNSigmaTOF(aPart, AliPID::kPion);
  double tPt = aPart->Pt();


  if (tPt < 0.5) 
  {
    if (TMath::Abs(tNSigmaTPCPi) < 3.) return true;
  } 
  else 
  {
    if (tNSigmaTOFPi < -999.) 
    {
      if (TMath::Abs(tNSigmaTPCPi) < 3.) return true;
    } 
    else 
    {
      if (TMath::Abs(tNSigmaTPCPi) < 3. && TMath::Abs(tNSigmaTOFPi) < 3.) return true;
    }
  }

  return false;
}


//____________________________
bool AliAnalysisV0Efficiency::IsProtonNSigma(const AliAODTrack *aPart)
{
  double tNSigmaTPCP = fpidAOD->NumberOfSigmasTPC(aPart, AliPID::kProton);
//  double tNSigmaTOFP = fpidAOD->NumberOfSigmasTOF(aPart, AliPID::kProton);
  double tNSigmaTOFP = GetNSigmaTOF(aPart, AliPID::kProton);
  double tPt = aPart->Pt();

  if (tPt < 0.8) 
  {
    if (TMath::Abs(tNSigmaTPCP) < 3.) return true;
  } 
  else 
  {
    if (tNSigmaTOFP < -999.) 
    {
      if (TMath::Abs(tNSigmaTPCP) < 3.) return true;
    } 
    else 
    {
      if (TMath::Abs(tNSigmaTPCP) < 3. && TMath::Abs(tNSigmaTOFP) < 3.) return true;
    }
  }

  return false;
}

//____________________________
bool AliAnalysisV0Efficiency::IsMisIDK0s(int aParticleType, const AliAODv0 *aV0, const AliAODTrack *aPartPos, const AliAODTrack *aPartNeg)
{
  double fInvMassRejectK0sMin=0.488614, fInvMassRejectK0sMax=0.506614;

  if(aV0->MassK0Short()<fInvMassRejectK0sMin || aV0->MassK0Short()>fInvMassRejectK0sMax) return false;
  if(!IsPionNSigma(aPartPos)) return false;
  if(!IsPionNSigma(aPartNeg)) return false;

  // At this point, the V0 passes as both a (Anti)Lambda (particle of interest) and a K0Short (MisID to be removed)
  // For now, use mass to decide winner
  double tK0ShortMass = 0.497614, tLambdaMass = 1.115683;
  if((aParticleType==0) && (TMath::Abs(aV0->MassLambda()-tLambdaMass) < TMath::Abs(aV0->MassK0Short()-tK0ShortMass))) return false;
  else if((aParticleType==1) && (TMath::Abs(aV0->MassAntiLambda()-tLambdaMass) < TMath::Abs(aV0->MassK0Short()-tK0ShortMass))) return false;
  else {}

  return true;


}

//____________________________
bool AliAnalysisV0Efficiency::IsMisIDLambda(const AliAODv0 *aV0, const AliAODTrack *aPartPos, const AliAODTrack *aPartNeg)
{
  double fInvMassRejectLambdaMin=1.106683, fInvMassRejectLambdaMax=1.124683;

  if(aV0->MassLambda()<fInvMassRejectLambdaMin || aV0->MassLambda()>fInvMassRejectLambdaMax) return false;
  if(!IsProtonNSigma(aPartPos)) return false;
  if(!IsPionNSigma(aPartNeg)) return false;

  // At this point, the V0 passes as both a K0Short (particle of interest) and a Lambda (MisID to be removed)
  // For now, use mass to decide winner
  double tK0ShortMass = 0.497614, tLambdaMass = 1.115683;
  if(TMath::Abs(aV0->MassK0Short()-tK0ShortMass) < TMath::Abs(aV0->MassLambda()-tLambdaMass)) return false;

  return true;
}

//____________________________
bool AliAnalysisV0Efficiency::IsMisIDAntiLambda(const AliAODv0 *aV0, const AliAODTrack *aPartPos, const AliAODTrack *aPartNeg)
{
  double fInvMassRejectAntiLambdaMin=1.106683, fInvMassRejectAntiLambdaMax=1.124683;

  if(aV0->MassAntiLambda()<fInvMassRejectAntiLambdaMin || aV0->MassAntiLambda()>fInvMassRejectAntiLambdaMax) return false;
  if(!IsPionNSigma(aPartPos)) return false;
  if(!IsProtonNSigma(aPartNeg)) return false;

  // At this point, the V0 passes as both a K0Short (particle of interest) and a AntiLambda (MisID to be removed)
  // For now, use mass to decide winner
  double tK0ShortMass = 0.497614, tLambdaMass = 1.115683;
  if(TMath::Abs(aV0->MassK0Short()-tK0ShortMass) < TMath::Abs(aV0->MassAntiLambda()-tLambdaMass)) return false;

  return true;
}



//____________________________
bool AliAnalysisV0Efficiency::V0PassBasicCuts(const AliAODv0* aV0, const AliAODVertex *aAodvertex)
{
  //Make sure the V0 satisifies a few basic criteria
  if(aV0->GetNDaughters() != 2)                           return false;
  if(aV0->GetNProngs() != 2)                              return false;
  if(aV0->GetCharge() != 0)                              return false;
  if(aV0->ChargeProng(0) == aV0->ChargeProng(1))         return false;

  double tPrimVertex[3];
  aAodvertex->GetPosition(tPrimVertex);
  if(aV0->CosPointingAngle(tPrimVertex) < 0.98)          return false;


  //Now look at daughter tracks
  if(!(AliAODTrack*)aV0->GetDaughterLabel(0))                 return false;
  if(!(AliAODTrack*)aV0->GetDaughterLabel(1))                 return false;

  // The V0 has passed the most basic track cuts.
  return true;
}

//____________________________
bool AliAnalysisV0Efficiency::V0PassAllCuts(int aParticleType, const AliAODv0 *aV0, const AliAODTrack *aPartPos, const AliAODTrack *aPartNeg, const AliAODVertex *aAodvertex)
{
  //aParticleType==0 -> Lambda; aParticleType==1 -> AntiLambda; aParticleType==2 -> K0Short

  if(!V0PassBasicCuts(aV0, aAodvertex)) return false;

  //Cut values set for Lambda analysis, i.e. aParticleType=0
  double fInvMassLambdaMin = 1.111883;
  double fInvMassLambdaMax = 1.119483;
  double fInvMassK0sMin = 0.483937;
  double fInvMassK0sMax = 0.517937;
  
  double fMinDcaDaughterPosToVert = 0.1;
  double fMinDcaDaughterNegToVert = 0.3;
  double fMaxDcaV0Daughters = 0.4;
  double fMaxDcaV0 = 0.5;
  double fMinDcaV0 = 0.0;
  double fMaxDecayLength = 60.;
  

  double fMinCosPointingAngle = 0.9993;
  double fEta = 0.8;
  double fPtMin = 0.4;
  double fPtMax = 100.;
  bool fOnFlyStatus = false;

  float fMaxEtaDaughters = 0.8;
  int fTPCNclsDaughters = 80.;
  int fNdofDaughters = 10.;
  unsigned long fStatusDaughters = AliESDtrack::kTPCrefit;
  float fPtMinPosDaughter = 0.5;
  float fPtMaxPosDaughter = 99.;
  float fPtMinNegDaughter = 0.16;
  float fPtMaxNegDaughter = 99.;

  if(aParticleType==1)
  {
    fMinDcaDaughterPosToVert = 0.3;
    fMinDcaDaughterNegToVert = 0.1;

    fPtMinPosDaughter = 0.16;
    fPtMinNegDaughter = 0.3;
  }
  if(aParticleType==2)
  {
    fMinDcaDaughterPosToVert = 0.3;
    fMinDcaDaughterNegToVert = 0.3;

    fMaxDcaV0Daughters = 0.3;
    fMaxDcaV0 = 0.3;
    fMaxDecayLength = 30.;

    fPtMin = 0.2;
    fPtMinPosDaughter = 0.15;
    fPtMinNegDaughter = 0.15;
  }


  if(TMath::Abs(aV0->PseudoRapV0()) > fEta) return false;

  if (aV0->Pt() < fPtMin || fPtMax < aV0->Pt()) return false;
  if (TMath::Abs(aPartPos->Eta()) > fMaxEtaDaughters) return false;
  if (TMath::Abs(aPartNeg->Eta()) > fMaxEtaDaughters) return false;

  if (aPartPos->Pt() < fPtMinPosDaughter) return false;
  if (aPartNeg->Pt() < fPtMinNegDaughter) return false;
  if (aPartPos->Pt() > fPtMaxPosDaughter) return false;
  if (aPartNeg->Pt() > fPtMaxNegDaughter) return false;
  

  //quality cuts
  if (aV0->GetOnFlyStatus() != fOnFlyStatus) return false;
  if (aPartNeg->GetStatus() == 999 || aPartPos->GetStatus() == 999) return false;
  if (aPartPos->GetTPCNcls() < fTPCNclsDaughters) return false;
  if (aPartNeg->GetTPCNcls() < fTPCNclsDaughters) return false;


  if (aPartPos->Chi2perNDF() > fNdofDaughters) return false;
  if (aPartNeg->Chi2perNDF() > fNdofDaughters) return false;
  if (!(aPartNeg->GetStatus() & fStatusDaughters)) return false;
  if (!(aPartPos->GetStatus() & fStatusDaughters)) return false;
  
  //fiducial volume radius
  if(aV0->RadiusV0()<0.0 || aV0->RadiusV0()>99999.0) 
    return false;

  //DCA between daughter particles
  if (TMath::Abs(aV0->DcaV0Daughters()) > fMaxDcaV0Daughters)
    return false;
  
  //DCA of daughters to primary vertex
  if (TMath::Abs(aV0->DcaPosToPrimVertex()) < fMinDcaDaughterPosToVert || TMath::Abs(aV0->DcaNegToPrimVertex()) < fMinDcaDaughterNegToVert)
    return false;

  //DCA V0 to prim vertex
  if (TMath::Abs(aV0->DcaV0ToPrimVertex()) > fMaxDcaV0 || TMath::Abs(aV0->DcaV0ToPrimVertex()) < fMinDcaV0)
    return false;


  double tPrimVertex[3];
  aAodvertex->GetPosition(tPrimVertex);

  //cos pointing angle
  if (aV0->CosPointingAngle(tPrimVertex) < fMinCosPointingAngle)
    return false;

  //decay length
  if (aV0->DecayLengthV0(tPrimVertex) > fMaxDecayLength)
    return false;




  bool pid_check = false;
  // Looking for lambdas = proton + pim
  if (aParticleType == 0) {
    if (IsProtonNSigma(aPartPos)) //proton
    {
      if (IsPionNSigma(aPartNeg)) //pion
      {
        pid_check = true;
        if(fRemoveMisidentified)
        {
          if(IsMisIDK0s(aParticleType, aV0, aPartPos, aPartNeg))
          {
            return false;
          }
        }
        fReconstructedPurityAid_Lam->Fill(aV0->MassLambda());
        //invariant mass lambda
        if (aV0->MassLambda() < fInvMassLambdaMin || aV0->MassLambda() > fInvMassLambdaMax) return false;    
      }
    }
  }

  //Looking for antilambdas =  antiproton + pip
  else if (aParticleType == 1) {
    if (IsProtonNSigma(aPartNeg)) //proton
    {
      if (IsPionNSigma(aPartPos)) //pion
      {
        pid_check = true;
        if(fRemoveMisidentified)
        {
          if(IsMisIDK0s(aParticleType, aV0, aPartPos, aPartNeg))
          {
            return false;
          }
        }
        fReconstructedPurityAid_ALam->Fill(aV0->MassAntiLambda());
        //invariant mass antilambda
        if (aV0->MassAntiLambda() < fInvMassLambdaMin || aV0->MassAntiLambda() > fInvMassLambdaMax) return false;
      }
    }
  }

  //Looking for K0s = pip + pim
  else if (aParticleType == 2) {
    if (IsPionNSigma(aPartNeg)) //pion-
    {
      if (IsPionNSigma(aPartPos)) //pion+
      {
        pid_check = true;
        if(fRemoveMisidentified)
        {
          if(IsMisIDLambda(aV0, aPartPos, aPartNeg))
          {
            return false;
          }
          else if(IsMisIDAntiLambda(aV0, aPartPos, aPartNeg))
          {
            return false;
          }
        }
        fReconstructedPurityAid_K0s->Fill(aV0->MassK0Short());
        //invariant mass K0s
        if (aV0->MassK0Short() < fInvMassK0sMin || aV0->MassK0Short() > fInvMassK0sMax) return false;
      }
    }
  }

  if(!pid_check) return false;
  return true;
}

//____________________________
void AliAnalysisV0Efficiency::ExtractV0FinderParticles(const AliAODEvent *aEvent)
{
  TClonesArray *mcP = NULL;
  mcP = (TClonesArray *) aEvent->FindListObject(AliAODMCParticle::StdBranchName());
  if (!mcP) 
  {
    cout << "AOD MC information requested, but no particle array found!" << endl;
    return;
  }
  //-----
  int tNumberOfLastHijingLabel = GetNumberOfLastHijingLabel(aEvent);
  if(tNumberOfLastHijingLabel <= 0) return;
  //-----
  const AliAODVertex *aodvertex = (AliAODVertex*)aEvent->GetPrimaryVertex();

  vector<vector<int> > tDaughterLabelsRecoLam(0), tDaughterLabelsRecoALam(0), tDaughterLabelsRecoK0s(0);
  vector<vector<int> > tDaughterLabelsV0Lam(0), tDaughterLabelsV0ALam(0), tDaughterLabelsV0K0s(0);

  for(int i = 0; i < aEvent->GetNumberOfV0s(); i++)
  {
    AliAODv0* tAODv0 = aEvent->GetV0(i);
    if(!tAODv0) continue;

    if(!V0PassBasicCuts(tAODv0, aodvertex)) continue;
    //-----------------------------------------------------------------

    bool tIsInjected = false;

    AliAODTrack* daughterTrackPos = (AliAODTrack*)tAODv0->GetDaughterLabel(0);
    AliAODTrack* daughterTrackNeg = (AliAODTrack*)tAODv0->GetDaughterLabel(1);

    daughterTrackPos->SetAODEvent(aEvent);
    daughterTrackNeg->SetAODEvent(aEvent);

    if (daughterTrackPos == NULL || daughterTrackNeg == NULL) continue;     // daughter tracks must exist
    if (daughterTrackNeg->Charge() == daughterTrackPos->Charge()) continue; // and have different charge

    // ensure that pos and neg are pointing to the correct children
    if (daughterTrackPos->Charge() < 0 && daughterTrackNeg->Charge() > 0) 
    {
      AliAODTrack *tmp = daughterTrackPos;
      daughterTrackPos = daughterTrackNeg;
      daughterTrackNeg = tmp;
    }

    const AliAODMCParticle *mcParticlePos=nullptr, *mcParticleNeg=nullptr;
    const AliAODMCParticle *tMCv0=nullptr;
    const AliAODMCParticle *tMother = nullptr;

    if (daughterTrackPos->GetLabel() <= 0 || daughterTrackNeg->GetLabel() <= 0) continue;
    mcParticlePos = (AliAODMCParticle*)(mcP->At(daughterTrackPos->GetLabel())),
    mcParticleNeg = (AliAODMCParticle*)(mcP->At(daughterTrackNeg->GetLabel()));
    if ((mcParticlePos == NULL) || (mcParticleNeg == NULL)) continue;

    // Get the mother ID of the two daughters
    const int motherOfPosID = mcParticlePos->GetMother(),
              motherOfNegID = mcParticleNeg->GetMother();

    // If both daughter tracks refer to the same mother, we can continue
    if(motherOfPosID < 0 || motherOfPosID != motherOfNegID) continue;

    // Our V0 particle
    tMCv0 = (AliAODMCParticle*)(mcP->At(motherOfPosID));
    if(!tMCv0) continue;

    tIsInjected = IsInjected(tMCv0, mcP, tNumberOfLastHijingLabel);
    if(tMCv0->GetMother() > -1) tMother = (AliAODMCParticle*)(mcP->At(tMCv0->GetMother())); //V0 has mother


    //-----------------------------------------------------------------
    if(fIgnoreInjectedV0s && tIsInjected) continue;

    if(fabs(mcParticlePos->Eta()) > 0.8) continue;
    if(fabs(mcParticleNeg->Eta()) > 0.8) continue;

    if(mcParticlePos->Pt() < 0.16) continue;
    if(mcParticleNeg->Pt() < 0.16) continue;

    if(mcParticlePos->GetPdgCode() == 2212 && mcParticlePos->Pt() < 0.5) continue;
    if(mcParticlePos->GetPdgCode() == -2212 && mcParticlePos->Pt() < 0.3) continue;

    if(mcParticleNeg->GetPdgCode() == 2212 && mcParticleNeg->Pt() < 0.5) continue;
    if(mcParticleNeg->GetPdgCode() == -2212 && mcParticleNeg->Pt() < 0.3) continue;




    //Fill fMCTruthOfV0FinderParticles histograms
    if(!ContainsSharedDaughters(tDaughterLabelsV0Lam, mcParticlePos->GetLabel(), mcParticleNeg->GetLabel()))  
    {
      FillHistograms( 3122, tMCv0, tMother, fMCTruthOfV0FinderParticles_Lam, fParticleOriginOfV0FinderParticles_Lam); 
    }
    if(!ContainsSharedDaughters(tDaughterLabelsV0ALam, mcParticlePos->GetLabel(), mcParticleNeg->GetLabel())) 
    {
      FillHistograms(-3122, tMCv0, tMother, fMCTruthOfV0FinderParticles_ALam, fParticleOriginOfV0FinderParticles_ALam); 
    }
    if(!ContainsSharedDaughters(tDaughterLabelsV0K0s, mcParticlePos->GetLabel(), mcParticleNeg->GetLabel()))  
    {
      FillHistograms(  310, tMCv0, tMother, fMCTruthOfV0FinderParticles_K0s, fParticleOriginOfV0FinderParticles_K0s); 
    }

    if(tMCv0->GetPdgCode()==3122)  tDaughterLabelsV0Lam.push_back(vector<int>{mcParticlePos->GetLabel(), mcParticleNeg->GetLabel()});
    if(tMCv0->GetPdgCode()==-3122) tDaughterLabelsV0ALam.push_back(vector<int>{mcParticlePos->GetLabel(), mcParticleNeg->GetLabel()});
    if(tMCv0->GetPdgCode()==310)   tDaughterLabelsV0K0s.push_back(vector<int>{mcParticlePos->GetLabel(), mcParticleNeg->GetLabel()});
    //-----------------------------------------------------------------

    if(V0PassAllCuts(0, tAODv0, daughterTrackPos, daughterTrackNeg, aodvertex)) 
    {
      if(!ContainsSharedDaughters(tDaughterLabelsRecoLam, mcParticlePos->GetLabel(), mcParticleNeg->GetLabel()))
      {
        FillHistograms( 3122, tMCv0, tMother, fMCTruthOfReconstructedParticles_Lam, fParticleOriginOfReconstructedParticles_Lam, true);
        tDaughterLabelsRecoLam.push_back(vector<int>{mcParticlePos->GetLabel(), mcParticleNeg->GetLabel()});
      }
    }
    if(V0PassAllCuts(1, tAODv0, daughterTrackPos, daughterTrackNeg, aodvertex)) 
    {
      if(!ContainsSharedDaughters(tDaughterLabelsRecoALam, mcParticlePos->GetLabel(), mcParticleNeg->GetLabel()))
      {
        FillHistograms(-3122, tMCv0, tMother, fMCTruthOfReconstructedParticles_ALam, fParticleOriginOfReconstructedParticles_ALam, true);
        tDaughterLabelsRecoALam.push_back(vector<int>{mcParticlePos->GetLabel(), mcParticleNeg->GetLabel()});
      }
    }
    if(V0PassAllCuts(2, tAODv0, daughterTrackPos, daughterTrackNeg, aodvertex)) 
    {
      if(!ContainsSharedDaughters(tDaughterLabelsRecoK0s, mcParticlePos->GetLabel(), mcParticleNeg->GetLabel()))
      {
        FillHistograms(  310, tMCv0, tMother, fMCTruthOfReconstructedParticles_K0s, fParticleOriginOfReconstructedParticles_K0s, true);
        tDaughterLabelsRecoK0s.push_back(vector<int>{mcParticlePos->GetLabel(), mcParticleNeg->GetLabel()});
      }
    }
  }


}


//____________________________
void AliAnalysisV0Efficiency::ExtractAll(const AliAODEvent *aEvent)
{
  ExtractOriginalParticles(aEvent);
  ExtractV0FinderParticles(aEvent);
}


//________________________________________________________________________
void AliAnalysisV0Efficiency::Exec(Option_t *)
{
  // Main loop
  // Called for each event
  if(!IsCorrectEventTrigger()) return;

  fEventCount++;
  fAOD = dynamic_cast<AliAODEvent*> (InputEvent());
  if (!fAOD) {Printf("ERROR: fAOD not available"); return;}

  //Centrality selection
/*
  AliCentrality *centrality = fAOD->GetCentrality();
  float centralityPercentile = centrality->GetCentralityPercentile("V0M");
  if(centralityPercentile > 10.) return;
*/

  ExtractAll(fAOD);

  // Post output data.
  PostData(1, fOutputList);
}


//________________________________________________________________________
void AliAnalysisV0Efficiency::Terminate(Option_t *)
{
  // Called once at the end of the query
  cout<<"Done"<<endl;

}

