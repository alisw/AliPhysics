//
// Class AliRsnAnalysisTaskSEBase
//
// TODO
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>

#include "AliLog.h"

#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"

#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliRsnEvent.h"
#include "AliMCEvent.h"

#include "AliRsnMCInfo.h"
#include "AliRsnDaughter.h"
#include "AliRsnAnalysisTaskSEBase.h"


ClassImp(AliRsnAnalysisTaskSEBase)

//_____________________________________________________________________________
AliRsnAnalysisTaskSEBase::AliRsnAnalysisTaskSEBase(const char *name) :
    AliAnalysisTaskSE(name),
    fUseAutoHandler(kTRUE),
    fReader(),
    //fPID(),
    fAnalysisMgr(0x0) // pointer to current AnalysisMgr
{
//
// Default constructor
//
  InitIOVars();
  DefineInput(0, TChain::Class());
}

//_____________________________________________________________________________
void AliRsnAnalysisTaskSEBase::InitIOVars()
{
//
// Initial values for constructor
//

  fUseAutoHandler = kFALSE;

  Int_t i;
  for (i = 0; i < 2; i++)
  {
    fChain[i] = 0;
    fRSN[i] = 0;
    fRsnESD[i] = 0;
    fRsnMC[i] = 0;
    fRsnAOD[i] = 0;
    fRsnESDEH[i] = 0;
    fRsnMCEH[i] = 0;
    fRsnAODEH[i] = 0;
    fInputType[i] = kRSN;
  }

  fAnalysisMgr = 0;
}

//_____________________________________________________________________________
void AliRsnAnalysisTaskSEBase::LocalInit()
{
//
// LocalInit()
//
  AliAnalysisTaskSE::LocalInit();
}

//_____________________________________________________________________________
Bool_t AliRsnAnalysisTaskSEBase::Notify()
{
//
// Notify()
//

  return AliAnalysisTaskSE::Notify();
}

//_____________________________________________________________________________
void AliRsnAnalysisTaskSEBase::SetInputType
(EInputType type, AliAnalysisManager* am, Bool_t autohandler, Short_t inputIndex)
{
//
// Sets input type.
// When autohandler is kTRUE handlers are created and connected
// to the passed AliAnalysisManager if it exists.
// The internal AliAnalysisManager object is redirected to the passed one.
//

  fInputType[inputIndex] = type;
  fAnalysisMgr = am;

  UseAutoHandler(autohandler);
  if (!fUseAutoHandler) return;

  if (!fAnalysisMgr)
  {
    AliWarning(Form("fAnalysisMgr is %p and fUseAutoHandler is %d", fAnalysisMgr, fUseAutoHandler));
    return;
  }

  switch (fInputType[inputIndex])
  {
    case kAOD:
      fRsnAODEH[0] = new AliAODInputHandler();
      if (fRsnAODEH[0]) fAnalysisMgr->SetInputEventHandler(fRsnAODEH[0]);
      break;
    case kESD:
    case kESDTPC:
      fRsnESDEH[0] = new AliESDInputHandler();
      if (fRsnESDEH[0]) fAnalysisMgr->SetInputEventHandler(fRsnESDEH[0]);
      break;
    case kESDMC:
    case kESDMCTPC:
    case kMC:
      fRsnESDEH[0] = new AliESDInputHandler();
      fRsnMCEH[0] = new AliMCEventHandler();
      if ((fRsnESDEH[0]) && (fRsnMCEH[0]))
      {
        fAnalysisMgr->SetInputEventHandler(fRsnESDEH[0]);
        fAnalysisMgr->SetMCtruthEventHandler(fRsnMCEH[0]);
      }
      break;
    case kRSN:
      break;
    default:
      AliError("Type not supported ...");
      break;
  }

  // check if the TPC only is used
  if (fInputType[inputIndex] == kESDTPC || fInputType[inputIndex] == kESDMCTPC) fReader.SetTPCOnly();
}

//_____________________________________________________________________________
void AliRsnAnalysisTaskSEBase::ConnectInputData(Option_t *)
{
//
// ConnectInputData() for AliAnalysisTaskSE
// just define myTask->SetInputType ( AliRsnAnalysisTaskSEBase::kRSN ); for Rsn input
// just define myTask->SetInputType ( AliRsnAnalysisTaskSEBase::kESDMC ); for ESD and MC input
// just define myTask->SetInputType ( AliRsnAnalysisTaskSEBase::kAOD ); for Rsn input
//

  if (fInputType[0] != kRSN) AliAnalysisTaskSE::ConnectInputData();

  // connects input handlers according to the type of analysis being done
  ConnectInputDataByInputType(fInputType[0], 0);
}

//_____________________________________________________________________________
void AliRsnAnalysisTaskSEBase::ConnectInputDataByInputType
(EInputType type, Short_t inputIndex)
{
//
// Connect input data dependint on the input type used.
//
  AliDebug(AliLog::kDebug, "<-");

  switch (type)
  {
    case kAOD:
      ConnectAOD(inputIndex);
      break;
    case kESD:
    case kESDTPC:
      ConnectESD(inputIndex);
      break;
    case kESDMC:
    case kESDMCTPC:
    case kMC:
      ConnectESDMC(inputIndex);
      break;
    case kRSN:
      ConnectRSN(inputIndex);
      break;
    default:
      AliError("Type not supported ...");
      break;
  }

  AliDebug(AliLog::kDebug, "->");
}
//_____________________________________________________________________________
void AliRsnAnalysisTaskSEBase::ConnectRSN(Short_t inputIndex)
{
//
// Connect input data by RSN input type
//

  AliDebug(AliLog::kDebug, "<-");

  char ** address = (char **) GetBranchAddress(inputIndex, "rsnEvents");
  if (address)
  {
    fRSN[inputIndex] = (AliRsnEvent*)(*address);
  }
  else
  {
    fRSN[inputIndex] = 0;
    SetBranchAddress(inputIndex, "rsnEvents", &fRSN[inputIndex]);
  }

  AliDebug(AliLog::kDebug, "->");
}

//_____________________________________________________________________________
void AliRsnAnalysisTaskSEBase::ConnectESD(Short_t inputIndex)
{
//
// Connect input data by ESD input type
//

  AliDebug(AliLog::kDebug, "<-");
  fRsnESD[inputIndex] = (AliESDEvent*)fInputEvent;
  AliDebug(AliLog::kDebug, "->");

}

//_____________________________________________________________________________
void AliRsnAnalysisTaskSEBase::ConnectESDMC(Short_t inputIndex)
{
//
// Connect input data by ESDMC input type
//

  AliDebug(AliLog::kDebug, "<-");
  fRsnESD[inputIndex] = (AliESDEvent*)fInputEvent;
  //fRSNMC[inputIndex] = fMCEvent;
  AliDebug(AliLog::kDebug, "->");
}

//_____________________________________________________________________________
void AliRsnAnalysisTaskSEBase::ConnectAOD(Short_t inputIndex)
{
//
// Connect input data by AOD input type
//

  AliDebug(AliLog::kDebug, "<-");

  TTree* tree = dynamic_cast<TTree*>(GetInputData(inputIndex));
  if (!tree) { AliError("Could not read chain from input slot 0");}
  else
  {
    AliAODInputHandler *aodH = dynamic_cast<AliAODInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    if (!aodH) { AliError("Could not get AODInputHandler"); }
    else fRsnAOD[inputIndex] = aodH->GetEvent();
  }

  AliDebug(AliLog::kDebug, "->");
}

//_____________________________________________________________________________
AliRsnEvent * AliRsnAnalysisTaskSEBase::GetRsnEventFromInputType(const Short_t & index)
{
//
// Gets Event from input type
//

  switch (fInputType[index])
  {
    case kAOD:      return GetRsnFromAOD(index);
    case kESD:      return GetRsnFromESD(index);
    case kESDMC:    return GetRsnFromESDMC(index);
    case kESDTPC:   return GetRsnFromESD(index);
    case kESDMCTPC: return GetRsnFromESDMC(index);
    case kMC:       return GetRsnFromMC(index);
    case kRSN:      return GetRsnFromRSN(index);
    default:
      AliError("Type not supported ...");
      return (AliRsnEvent*) 0x0;
  }

  return (AliRsnEvent*) 0x0;
}

//_____________________________________________________________________________
AliRsnEvent * AliRsnAnalysisTaskSEBase::GetRsnFromAOD(const Short_t & index)
{
//
// Gets RSN event from AOD
//

  if (!fRsnAOD[index])
  {
    AliError("fRsnAOD not available.");
    return (AliRsnEvent *) 0x0;
  }

  if (!fRSN[0])
  {
    fRSN[0] = new AliRsnEvent();
    fRSN[0]->SetName("rsnEvents");
    fRSN[0]->Init();
  }

  // clear pevious event
  fRSN[0]->Clear();
  if (!fReader.FillFromAOD(fRSN[0], fRsnAOD[index])) return (AliRsnEvent*) 0x0;
  //if (!fPID.Process(fRSN[0])) AliWarning("Failed PID");

  return (AliRsnEvent*) fRSN[0];
}

//_____________________________________________________________________________
AliRsnEvent * AliRsnAnalysisTaskSEBase::GetRsnFromESD(const Short_t & index)
{
//
// Gets RSN event from ESD
//

  if (!fRsnESD[index])
  {
    AliError("fRsnESD not available.");
    return (AliRsnEvent *) 0x0;
  }

  if (!fRSN[index])
  {
    fRSN[index] = new AliRsnEvent();
    fRSN[index]->SetName("rsnEvents");
    fRSN[index]->Init();
  }

  // clear pevious event
  fRSN[index]->Clear();

  if (!fReader.FillFromESD(fRSN[index], fRsnESD[index], 0x0)) return (AliRsnEvent*) 0x0;

  //if (!fPID.Process(fRSN[index]))
  //{
  //  AliWarning("Failed PID");
  //  return (AliRsnEvent*) 0x0;
  //}

  return fRSN[index];
}

//_____________________________________________________________________________
AliRsnEvent * AliRsnAnalysisTaskSEBase::GetRsnFromMC(const Short_t & index)
{
//
// Gets RSN event from ESD
//

  if (!fRsnESD[index])
  {
    AliError("fRsnESD not available.");
    return (AliRsnEvent *) 0x0;
  }

  if (!fRSN[index])
  {
    fRSN[index] = new AliRsnEvent();
    fRSN[index]->SetName("rsnEvents");
    fRSN[index]->Init();
  }

  // clear pevious event
  fRSN[index]->Clear();
  fRsnMC[index] = MCEvent();

  if (!fRsnMC[index]) return (AliRsnEvent *) 0x0;
  if (!fReader.FillFromMC(fRSN[index], fRsnMC[index])) return (AliRsnEvent*) 0x0;
  //fPID.Process(fRSN[index]);
  return fRSN[index];
}

//_____________________________________________________________________________
AliRsnEvent * AliRsnAnalysisTaskSEBase::GetRsnFromESDMC(const Short_t & index)
{
//
// Gets RSN event from ESD and MC
//

  if (!fRsnESD[index])
  {
    AliError("fRsnESD not available.");
    return (AliRsnEvent *) 0x0;
  }

  if (!fRSN[index])
  {
    fRSN[index] = new AliRsnEvent();
    fRSN[index]->SetName("rsnEvents");
    fRSN[index]->Init();
  }

  // clear pevious event
  fRSN[index]->Clear();
  fRsnMC[index] = MCEvent();

  if (!fRsnMC[index]) return (AliRsnEvent *) 0x0;
  if (!fReader.FillFromESD(fRSN[index], fRsnESD[index], fRsnMC[index])) return (AliRsnEvent*) 0x0;
  //if (!fPID.Process(fRSN[index]))
  //{
  //  AliWarning("Failed PID");
  //  return (AliRsnEvent*) 0x0;
  //}

  return fRSN[index];
}

//_____________________________________________________________________________
AliRsnEvent * AliRsnAnalysisTaskSEBase::GetRsnFromRSN(const Short_t & index)
{
//
// Gets RSN event from RSN
// not fully implemented yet
//
  AliRsnEvent *event = fRSN[index];
  return event;
}
