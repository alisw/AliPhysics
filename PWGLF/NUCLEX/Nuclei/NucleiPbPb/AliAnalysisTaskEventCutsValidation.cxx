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
}

/// This function creates all the histograms and all the objects in general used during the analysis
/// \return void
///
void AliAnalysisTaskEventCutsValidation::UserCreateOutputObjects() {

  fList = new TList();
  fList->SetOwner(true);
  fEventCut.AddQAplotsToList(fList,true);
  fEventCut.fUseVariablesCorrelationCuts = true;

  PostData(1,fList);
  if (fStoreCuts) PostData(2,&fEventCut);
}

/// This is the function that is evaluated for each event. The analysis code stays here.
///
/// \param options Deprecated parameter
/// \return void
///
void AliAnalysisTaskEventCutsValidation::UserExec(Option_t *) {
  AliVEvent* ev = InputEvent();
  fEventCut.AcceptEvent(ev);
  PostData(1,fList);
  if (fStoreCuts) PostData(2,&fEventCut);
  return;
}

/// Merge the output. Called once at the end of the query.
///
/// \return void
///
void AliAnalysisTaskEventCutsValidation::Terminate(Option_t *) {
  return;
}

