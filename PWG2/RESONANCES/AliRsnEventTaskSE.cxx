//
// Class AliRsnEventTaskSE
//
// TODO
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#include <Riostream.h>
#include <TList.h>

#include "AliLog.h"

#include "AliAnalysisManager.h"

#include "AliRsnEvent.h"
#include "AliRsnEventFunction.h"
#include "AliRsnEventTaskSE.h"

ClassImp(AliRsnEventTaskSE)

//________________________________________________________________________
AliRsnEventTaskSE::AliRsnEventTaskSE(const char * name) :
  AliRsnAnalysisTaskSEBase(name),
  fEventFunctions("AliRsnEventFunction", 0),
  fOutList(0x0)
{
//
// Default constructor
//

  InitIOVars();
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliRsnEventTaskSE::~AliRsnEventTaskSE()
{
//
// Destructor
//
}

//________________________________________________________________________
void AliRsnEventTaskSE::InitIOVars()
{
//
// Init input output values
//

  AliRsnAnalysisTaskSEBase::InitIOVars();
  fOutList = 0;
}

//________________________________________________________________________
void AliRsnEventTaskSE::UserCreateOutputObjects()
{
//
// UserCreateOutputObjects() of AliAnalysisTaskSE
//

  OpenFile(0);
  fOutList = new TList();
  fOutList->SetOwner();
  fOutList->SetName("EventFunctions");

  AliRsnEventFunction *fcn = 0;
  TObjArrayIter next(&fEventFunctions);

  while ( (fcn = (AliRsnEventFunction*)next()) )
  {
    fcn->Init(fOutList);
  }

  AliInfo("*************************** List");
  fOutList->Print();
  AliInfo("*********************** End List");
}

//________________________________________________________________________
void AliRsnEventTaskSE::UserExec(Option_t *)
{
//
// UserExec() of AliAnalysisTaskSE
//

  if (fEntry++ % 100 == 0) cout << "[" << GetName() << "]: event " << fEntry-1 << endl;

  AliRsnEvent *curEvent = GetRsnEventFromInputType();
  if (!curEvent) return;

  AliRsnEventFunction *fcn = 0;
  TObjArrayIter next(&fEventFunctions);

  while ( (fcn = (AliRsnEventFunction*)next()) )
  {
    fcn->Fill(curEvent);
  }

  PostData(1, fOutList);
}

//________________________________________________________________________
void AliRsnEventTaskSE::Terminate(Option_t *)
{
//
// Terminate() of AliAnalysisTask
//

  fOutList = dynamic_cast<TList*>(GetOutputData(1));
  if (!fOutList) { AliError("*** Output list not available ***"); return; }
  fOutList->Print();
}

//________________________________________________________________________
void AliRsnEventTaskSE::AddEventFunction(AliRsnEventFunction *fcn)
{
  Int_t size = fEventFunctions.GetEntries();
  new (fEventFunctions[size]) AliRsnEventFunction(*fcn);
  AliInfo(Form("Function name: %s", fcn->GetFcnName().Data()));
}
