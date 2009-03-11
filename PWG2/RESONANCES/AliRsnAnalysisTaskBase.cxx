//
// Class AliRsnAnalysisTaskBase
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
#include "AliAnalysisTask.h"

#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliRsnEvent.h"
#include "AliMCEvent.h"

#include "AliRsnAnalysisTaskBase.h"


ClassImp(AliRsnAnalysisTaskBase)

//________________________________________________________________________
AliRsnAnalysisTaskBase::AliRsnAnalysisTaskBase(const char *name)
    : AliAnalysisTask(name, ""),
    fNumOfEvents(100),
    fUseAutoHandler(100),
    fReader(),
    fPID(),
    fAnalysisMgr(0x0)
{
//=========================================================
// Default constructor
//=========================================================
  InitIOVars();
  DefineInput(0, TChain::Class());
//   DefineInput ( 1, AliRsnReader::Class() );
//   DefineInput ( 2, AliRsnPID::Class() );
}

//________________________________________________________________________
void AliRsnAnalysisTaskBase::InitIOVars()
{
//=========================================================
// Initial values for constructor
//=========================================================
  AliDebug(AliLog::kDebug, "<-");

  fNumOfEvents=0;
  fUseAutoHandler = kFALSE;
  for (Int_t i=0;i<2;i++)
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

  fAnalysisMgr=0;

  AliDebug(AliLog::kDebug, "->");
}

//________________________________________________________________________
void AliRsnAnalysisTaskBase::LocalInit()
{
//=========================================================
// LocalInit()
//=========================================================
}

//________________________________________________________________________
Bool_t AliRsnAnalysisTaskBase::Notify()
{
//=========================================================
// Notify()
//=========================================================

  return AliAnalysisTask::Notify();
}

//_____________________________________________________________________________
void AliRsnAnalysisTaskBase::SetInputType(EInputType type, AliAnalysisManager* am, Bool_t autohandler, Short_t inputIndex)
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
//________________________________________________________________________
void AliRsnAnalysisTaskBase::ConnectInputData(Option_t *)
{
//=========================================================
// ConectInputData() for AliAnalysisTask
// just define myTask->SetInputType ( AliRsnAnalysisTaskBase::kRSN ); for Rsn input
// just define myTask->SetInputType ( AliRsnAnalysisTaskBase::kESDMC ); for ESD and MC input
// just define myTask->SetInputType ( AliRsnAnalysisTaskBase::kAOD ); for Rsn input
//=========================================================

  ConnectInputDataByInputType(fInputType[0],0);
}

void AliRsnAnalysisTaskBase::ConnectInputDataByInputType(EInputType type ,Short_t inputIndex)
{
//=========================================================
// Connect input data by input type
//=========================================================

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
void AliRsnAnalysisTaskBase::ConnectRSN(Short_t inputIndex)
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
void AliRsnAnalysisTaskBase::ConnectESD(Short_t inputIndex)
{
//
// Connect input data by ESD input type
//
  AliDebug(AliLog::kDebug, "<-");
  TTree* tree = dynamic_cast<TTree*>(GetInputData(inputIndex));
  if (!tree) { AliError("Could not read chain from input slot 0"); }
  else
  {
    // Disable all branches, we want to process only MC
//     tree->SetBranchStatus("*", kFALSE);
//     tree->SetBranchStatus("fTracks.*", kTRUE);

    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

    if (!esdH) { AliError("Could not get ESDInputHandler"); }
    else
      fRsnESD[inputIndex] = esdH->GetEvent();
  }
  AliDebug(AliLog::kDebug, "->");
}

//_____________________________________________________________________________
void AliRsnAnalysisTaskBase::ConnectESDMC(Short_t inputIndex)
{
//
// Connect input data by ESDMC input type
//

  AliDebug(AliLog::kDebug, "<-");
  TTree* tree = dynamic_cast<TTree*>(GetInputData(inputIndex));
  if (!tree) { AliError("Could not read chain from input slot 0"); }
  else
  {
    // Disable all branches, we want to process only MC
//     tree->SetBranchStatus("*", kFALSE);
//     tree->SetBranchStatus("fTracks.*", kTRUE);

    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

    if (!esdH) { AliError("Could not get ESDInputHandler"); }
    else
      fRsnESD[inputIndex] = esdH->GetEvent();
  }
  AliDebug(AliLog::kDebug, "->");
}

//_____________________________________________________________________________
void AliRsnAnalysisTaskBase::ConnectAOD(Short_t inputIndex)
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
AliRsnEvent * AliRsnAnalysisTaskBase::GetRsnEventFromInputType(const Short_t & index)
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
AliRsnEvent * AliRsnAnalysisTaskBase::GetRsnFromAOD(const Short_t & index)
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
  if (!fPID.Process(fRSN[0])) AliWarning("Failed PID");

  return (AliRsnEvent*) fRSN[0];
}

//_____________________________________________________________________________
AliRsnEvent * AliRsnAnalysisTaskBase::GetRsnFromESD(const Short_t & index)
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

  if (!fPID.Process(fRSN[index]))
  {
    AliWarning("Failed PID");
    return (AliRsnEvent*) 0x0;
  }

  return fRSN[index];
}

//_____________________________________________________________________________
AliRsnEvent * AliRsnAnalysisTaskBase::GetRsnFromMC(const Short_t & index)
{
//
// Gets RSN event from ESD
//

  AliMCEventHandler* mcHandler = dynamic_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  if (!mcHandler) { AliError("Could not retrieve MC event handler"); return (AliRsnEvent *) 0x0; }

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

  fRsnMC[index] =  mcHandler->MCEvent();

  if (!fRsnMC[index]) return (AliRsnEvent *) 0x0;
  if (!fReader.FillFromMC(fRSN[index], fRsnMC[index])) return (AliRsnEvent*) 0x0;
  fPID.Process(fRSN[index]);
  return fRSN[index];
}

//_____________________________________________________________________________
AliRsnEvent * AliRsnAnalysisTaskBase::GetRsnFromESDMC(const Short_t & index)
{
//
// Gets RSN event from ESD and MC
//

  AliMCEventHandler* mcHandler = dynamic_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  if (!mcHandler) { AliError("Could not retrieve MC event handler"); return (AliRsnEvent *) 0x0; }

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
  fRsnMC[index] =  mcHandler->MCEvent();

  if (!fRsnMC[index]) return (AliRsnEvent *) 0x0;
  if (!fReader.FillFromESD(fRSN[index], fRsnESD[index], fRsnMC[index])) return (AliRsnEvent*) 0x0;
  if (!fPID.Process(fRSN[index]))
  {
    AliWarning("Failed PID");
    return (AliRsnEvent*) 0x0;
  }

  return fRSN[index];
}

//_____________________________________________________________________________
AliRsnEvent * AliRsnAnalysisTaskBase::GetRsnFromRSN(const Short_t & index)
{
//
// Gets RSN event from RSN
// not fully implemented yet
//
  AliRsnEvent *event = fRSN[index];
  return event;
}
