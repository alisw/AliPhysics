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
    fRsnInput(),
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
    fESD[i] = 0;
    fMC[i] = 0;
    fAOD[i] = 0;
    fESDEH[i] = 0;
    fMCEH[i] = 0;
    fAODEH[i] = 0;
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

//________________________________________________________________________
void AliRsnAnalysisTaskBase::SetInputType(EInputType type,AliAnalysisManager* am,Bool_t autohandler, Short_t inputIndex)
{
//=========================================================
// Sets input type.
// When autohandler is kTRUE handlers are created and sets
// to AliAnalysisManager (fAnalysisMgr) if exists
//=========================================================

  AliDebug(AliLog::kDebug, "<-");

  fInputType[inputIndex] = type;
  fAnalysisMgr = am;

  UseAutoHandler(autohandler);

  if (!fUseAutoHandler) return;

  if (!fAnalysisMgr)
  {
    AliWarning(Form("fAnalysisMgr is %p and fUseAutoHandler is %d",fAnalysisMgr,fUseAutoHandler));
    return;
  }

  switch (fInputType[inputIndex])
  {
    case kAOD:
    {
      if (fAnalysisMgr)
      {
        fAODEH[0] = new AliAODInputHandler();
        if (fAODEH[0])
          fAnalysisMgr->SetInputEventHandler(fAODEH[0]);
      }
      break;
    }
    case kESD:
    {
      if (fAnalysisMgr)
      {
        fESDEH[0] = new AliESDInputHandler();
        if (fESDEH[0])
          fAnalysisMgr->SetInputEventHandler(fESDEH[0]);
      }

      break;
    }
    case kESDMC:
    {
      if (fAnalysisMgr)
      {
        fESDEH[0] = new AliESDInputHandler();
        fMCEH[0] = new AliMCEventHandler();
        if ((fESDEH[0]) && (fMCEH[0]))
        {
          fAnalysisMgr->SetInputEventHandler(fESDEH[0]);
          fAnalysisMgr->SetMCtruthEventHandler(fMCEH[0]);
        }
      }
      break;
    }
    case kMC:
      AliError("Not Implemented Yet ...");
      break;
    case kRSN:
    {
      AliError("Not Implemented Yet ...");
      break;
    }
    default:
      AliError("Type not supported ...");
      break;
  }
  AliDebug(AliLog::kDebug, "->");
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
    {
      ConnectAOD(inputIndex);
      break;
    }
    case kESD:
    {
      ConnectESD(inputIndex);
      break;
    }
    case kESDMC:
    {
      ConnectESDMC(inputIndex);
      break;
    }
    case kMC:
      AliError("Not Implemented Yet ...");
      break;
    case kRSN:
    {
      ConnectRSN(inputIndex);
      break;
    }
    default:
      AliError("Type not supported ...");
      break;
  }
  AliDebug(AliLog::kDebug, "->");
}
//________________________________________________________________________
void AliRsnAnalysisTaskBase::ConnectRSN(Short_t inputIndex)
{
//=========================================================
// Connect input data by RSN input type
//=========================================================

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

void AliRsnAnalysisTaskBase::ConnectESD(Short_t inputIndex)
{
//=========================================================
// Connect input data by ESD input type
//=========================================================

  AliDebug(AliLog::kDebug, "<-");

  TTree* tree = dynamic_cast<TTree*>(GetInputData(inputIndex));
  if (!tree) { AliError("Could not read chain from input slot 0"); }
  else
  {
    // Disable all branches, we want to process only MC
    tree->SetBranchStatus("*", kFALSE);
    tree->SetBranchStatus("fTracks.*", kTRUE);

    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

    if (!esdH) { AliError("Could not get ESDInputHandler"); }
    else
      fESD[inputIndex] = esdH->GetEvent();
  }
  AliDebug(AliLog::kDebug, "->");

}
//________________________________________________________________________
void AliRsnAnalysisTaskBase::ConnectESDMC(Short_t inputIndex)
{
//=========================================================
// Connect input data by ESDMC input type
//=========================================================

  AliDebug(AliLog::kDebug, "<-");


  TTree* tree = dynamic_cast<TTree*>(GetInputData(inputIndex));
  if (!tree) { AliError("Could not read chain from input slot 0"); }
  else
  {
    // Disable all branches, we want to process only MC
    tree->SetBranchStatus("*", kFALSE);
    tree->SetBranchStatus("fTracks.*", kTRUE);

    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

    if (!esdH) { AliError("Could not get ESDInputHandler"); }
    else
      fESD[inputIndex] = esdH->GetEvent();
  }
  AliDebug(AliLog::kDebug, "->");

}
//________________________________________________________________________
void AliRsnAnalysisTaskBase::ConnectAOD(Short_t inputIndex)
{
//=========================================================
// Connect input data by AOD input type
//=========================================================

  AliDebug(AliLog::kDebug, "<-");

  TTree* tree = dynamic_cast<TTree*>(GetInputData(inputIndex));
  if (!tree) { AliError("Could not read chain from input slot 0");}
  else
  {
    AliAODInputHandler *aodH = dynamic_cast<AliAODInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

    if (!aodH) { AliError("Could not get AODInputHandler"); }
    else
    {
      fAOD[inputIndex] = aodH->GetEvent();
    }
  }
  AliDebug(AliLog::kDebug, "->");
}

//________________________________________________________________________
AliRsnEvent * AliRsnAnalysisTaskBase::GetRsnEventFromInputType(const Short_t & index)
{
//=========================================================
// Gets Evetn from input type
//=========================================================

  switch (fInputType[index])
  {
    case kAOD:
    {
      return GetRsnFromAOD(index);
      break;
    }
    case kESD:
    {
      AliWarning("Not Implemented Yet ...");
      return GetRsnFromESD(index);
      break;
    }
    case kESDMC:
    {
      return GetRsnFromESDMC(index);
      break;
    }
    case kMC:
      AliWarning("Not Implemented Yet ...");
      return (AliRsnEvent*) 0x0;
      break;
    case kRSN:
    {
      return GetRsnFromRSN(index);
      break;
    }
    default:
      AliError("Type not supported ...");
      return (AliRsnEvent*) 0x0;
      break;
  }
  return (AliRsnEvent*) 0x0;
}



//________________________________________________________________________
AliRsnEvent * AliRsnAnalysisTaskBase::GetRsnFromAOD(const Short_t & index)
{
//=========================================================
// Gets RSN event from AOD
//=========================================================

  if (!fAOD[index]) { AliError("fAOD not available."); return (AliRsnEvent *) 0x0; }

  if (!fRSN[0])
  {
    fRSN[0] = new AliRsnEvent();
    fRSN[0]->SetName("rsnEvents");
    fRSN[0]->Init();
  }
  // clear pevious event
  fRSN[0]->Clear();

  if (!fReader.Fill(fRSN[0], (AliVEvent*)  fAOD[index]))
  {
    return (AliRsnEvent*) 0x0;
  };

  if (!fPID.Process(fRSN[0])) AliWarning("Failed PID");

  return (AliRsnEvent*) fRSN[0];

}
//________________________________________________________________________
AliRsnEvent * AliRsnAnalysisTaskBase::GetRsnFromESD(const Short_t & index)
{
//=========================================================
// Gets RSN event from ESD
//=========================================================

  if (!fESD[index]) { AliError("fESD not available."); return (AliRsnEvent *) 0x0; }

  if (!fRSN[index])
  {
    fRSN[index] = new AliRsnEvent();
    fRSN[index]->SetName("rsnEvents");
    fRSN[index]->Init();
  }
  // clear pevious event
  fRSN[index]->Clear();

  if (!fReader.Fill(fRSN[index], (AliVEvent*) fESD[index]))
  {
    return (AliRsnEvent*) 0x0;
  };

  if (!fPID.Process(fRSN[index])) AliWarning("Failed PID");

  return fRSN[index];
}
//________________________________________________________________________
AliRsnEvent * AliRsnAnalysisTaskBase::GetRsnFromESDMC(const Short_t & index)
{
//=========================================================
// Gets RSN event from ESD and MC
//=========================================================

  if (!fESD[index]) { AliError("fESD not available."); return (AliRsnEvent *) 0x0; }
  AliMCEventHandler* mcHandler = dynamic_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  if (!mcHandler) { AliError("Could not retrieve MC event handler"); return (AliRsnEvent *) 0x0; }

  if (!fRSN[index])
  {
    fRSN[index] = new AliRsnEvent();
    fRSN[index]->SetName("rsnEvents");
    fRSN[index]->Init();
  }
  // clear pevious event
  fRSN[index]->Clear();

  fMC[index] = mcHandler->MCEvent();

  if (!fMC[index]) return (AliRsnEvent *) 0x0;

  Bool_t success = fReader.FillFromESD(fRSN[index], fESD[index], fMC[index]);
  if (!success)
  {
    AliInfo("Failed filling");
    return (AliRsnEvent*) 0x0;
  };

  if (!fPID.Process(fRSN[index]))
  {
    AliWarning("Failed PID");
    return (AliRsnEvent*) 0x0;
  }

  return fRSN[index];
}
//________________________________________________________________________
AliRsnEvent * AliRsnAnalysisTaskBase::GetRsnFromRSN(const Short_t & index)
{
//=========================================================
// Gets RSN event from RSN
// not fully implemented yet
//=========================================================

  AliRsnEvent *event = fRSN[index];
//   if ( fRsnEventBuffer->GetDeleteBufferWhenReset() == kTRUE )
//   {
//     event = ( AliRsnEvent * ) fRSN[index]->Clone();
//   }
//   AliInfo ( Form ( "%p %p",event,fRSN[index] ) );
  return event;
}

//________________________________________________________________________
void AliRsnAnalysisTaskBase::UseAutoHandler(const Bool_t & theValue)
{
//=========================================================
// Sets should create handlers
//=========================================================
  fUseAutoHandler = theValue;
}
