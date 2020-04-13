#include "AliAnalysisTaskEventCutsValidation.h"

// ROOT includes
#include <TChain.h>
#include <TH2F.h>
#include <TList.h>

// ALIROOT includes
#include "AliAODEvent.h"
#include "AliAODHeader.h"
#include "AliAODTrack.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliVEvent.h"

///\cond CLASSIMP
ClassImp(AliAnalysisTaskEventCutsValidation);
///\endcond

/// Standard and default constructor of the class.
///
/// \param taskname Name of the task
/// \param partname Name of the analysed particle
///
AliAnalysisTaskEventCutsValidation::AliAnalysisTaskEventCutsValidation(bool storeCuts, TString taskname) :
  AliAnalysisTaskSE(taskname.Data()),
  fFillTree{false},
  fEventCut(false),
  fList(nullptr),
  fStoreCuts(storeCuts)
{
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());

  if (storeCuts) DefineOutput(2, AliEventCuts::Class());
}

/// Standard destructor
///
AliAnalysisTaskEventCutsValidation::~AliAnalysisTaskEventCutsValidation() {
  if (fList) delete fList;
  if (fTree) delete fList;
}

/// This function creates all the histograms and all the objects in general used during the analysis
/// \return void
///
void AliAnalysisTaskEventCutsValidation::UserCreateOutputObjects() {

  fList = new TList();
  fList->SetOwner(true);
  fEventCut.AddQAplotsToList(fList,true);

  fTree = new TTree("EventSummary","Event Summary");
  fTree->Branch("Event",&fCurrentEvent);

  PostData(1,fList);
  if (fStoreCuts) PostData(2,&fEventCut);
  if (fFillTree) PostData(3, fTree);
}

/// This is the function that is evaluated for each event. The analysis code stays here.
///
/// \param options Deprecated parameter
/// \return void
///
void AliAnalysisTaskEventCutsValidation::UserExec(Option_t *) {
  AliVEvent* ev = InputEvent();
  bool acc = fEventCut.AcceptEvent(ev);


  PostData(1,fList);
  if (fStoreCuts) PostData(2,&fEventCut);
  if (acc) {
    fCurrentEvent.trigger = ev->GetTriggerMask();
    fCurrentEvent.x = fEventCut.GetPrimaryVertex()->GetX();
    fCurrentEvent.y = fEventCut.GetPrimaryVertex()->GetY();
    fCurrentEvent.z = fEventCut.GetPrimaryVertex()->GetZ();
    fCurrentEvent.v0m = fEventCut.GetCentrality(0);
    fCurrentEvent.cl0 = fEventCut.GetCentrality(1);
    fCurrentEvent.esd = fEventCut.fContainer.fMultESD;
    fCurrentEvent.fb32 = fEventCut.fContainer.fMultTrkFB32;
    fCurrentEvent.fb32acc = fEventCut.fContainer.fMultTrkFB32Acc;
    fCurrentEvent.fb32tof = fEventCut.fContainer.fMultTrkFB32TOF;
    fCurrentEvent.tpc = fEventCut.fContainer.fMultTrkTPC;
    fCurrentEvent.tpcOut = fEventCut.fContainer.fMultTrkTPCout;
    fCurrentEvent.multvzero = fEventCut.fContainer.fMultVZERO;
    if (fFillTree) {
      fTree->Fill();
      PostData(3,fTree);
    }
  }
  return;
}

/// Merge the output. Called once at the end of the query.
///
/// \return void
///
void AliAnalysisTaskEventCutsValidation::Terminate(Option_t *) {
  return;
}

