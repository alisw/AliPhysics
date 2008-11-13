//
// Class AliRsnAnalysisSE
//
// TODO
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#include <TSystem.h>
#include <TFile.h>
#include <TFolder.h>
#include <TROOT.h>

#include "AliLog.h"

#include "AliAnalysisManager.h"
#include "AliRsnPairMgr.h"
#include "AliRsnEventBuffer.h"

#include "AliMCEventHandler.h"
#include "AliESDEvent.h"

#include "AliRsnAnalysisSE.h"

ClassImp(AliRsnAnalysisSE)

//________________________________________________________________________
AliRsnAnalysisSE::AliRsnAnalysisSE(const char * name)
    : AliRsnAnalysisTaskSEBase(name),
    fPairMgrs(0),fOutList(0x0),fRsnEventBuffer(0x0),
    fNumOfEventsInBuffer(1000)
{
//=========================================================
// Default constructor
//=========================================================

  InitIOVars();
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliRsnAnalysisSE::~AliRsnAnalysisSE()
{
//=========================================================
// Destructor
//=========================================================
}

//________________________________________________________________________
void AliRsnAnalysisSE::InitIOVars()
{
//=========================================================
// Init input output values
//=========================================================

  AliDebug(AliLog::kDebug, "<-");
  AliRsnAnalysisTaskSEBase::InitIOVars();

  fRsnEventBuffer = 0;
  fOutList = 0;

  AliDebug(AliLog::kDebug, "->");
}

//________________________________________________________________________
void AliRsnAnalysisSE::UserCreateOutputObjects()
{
//=========================================================
// UserCreateOutputObjects() of AliAnalysisTaskSE
//=========================================================

  AliDebug(AliLog::kDebug, "<-");
//   fPID.DumpPriors();
  OpenFile(0);
  fOutList = new TList();
  fOutList->SetOwner();
  AliRsnPair *def=0;
  AliRsnPairMgr *mgr=0;
  TList *listTmp;
  for (Int_t iMgr=0 ;iMgr< fPairMgrs.GetEntries();iMgr++)
  {
    mgr = (AliRsnPairMgr *) fPairMgrs.At(iMgr);
    if (!mgr) continue;
    listTmp = new TList();
    listTmp->SetName(mgr->GetName());
    for (Int_t i=0;i< mgr->GetPairs()->GetEntriesFast();i++)
    {
      def = (AliRsnPair *) mgr->GetPairs()->At(i);
      if (def)
      {
        def->Init();
        listTmp->Add(def->GenerateHistograms(mgr->GetName()));
        //def->GenerateHistograms(mgr->GetName(), listTmp);
        //def->Print();
      }
    }
    fOutList->Add(listTmp);
  }

  fRsnEventBuffer = new AliRsnEventBuffer(fNumOfEventsInBuffer);
//   fRsnEventBuffer = new AliRsnEventBuffer ( 10000 ,kFALSE );
  AliDebug(AliLog::kDebug, "->");

}

//________________________________________________________________________
void AliRsnAnalysisSE::UserExec(Option_t *)
{
//=========================================================
// UserExec() of AliAnalysisTaskSE
//=========================================================

  if (fEntry++%1000==0)
    AliInfo(Form("Event %d",fEntry-1));

  AliRsnEvent *curEvent = GetRsnEventFromInputType();
  if (!curEvent) return;

  ProcessEventAnalysis(curEvent);
  PostEventProcess();

  PostData(1, fOutList);
}

//________________________________________________________________________
void AliRsnAnalysisSE::Terminate(Option_t *)
{
//=========================================================
// Terminate() of AliAnalysisTask
//=========================================================

  fOutList = dynamic_cast<TList*>(GetOutputData(1));
  if (!fOutList) { AliError(" fOutList not available"); return; }
  fOutList->Print();
}

//________________________________________________________________________
void AliRsnAnalysisSE::ProcessEventAnalysis(AliRsnEvent *curEvent)
{
//=========================================================
// Process of one event
//=========================================================


  // Adds event to Event Buffer
  fRsnEventBuffer->AddEvent(curEvent);

  // loop over all Pair managers
  AliRsnPairMgr *mgr=0;
  for (Int_t iMgr=0 ;iMgr< fPairMgrs.GetEntries();iMgr++)
  {
    mgr = (AliRsnPairMgr *) fPairMgrs.At(iMgr);
    AliRsnPair *pair=0;
    for (Int_t i=0;i< mgr->GetPairs()->GetEntriesFast();i++)
    {
      pair = (AliRsnPair *) mgr->GetPairs()->At(i);
      pair->ProcessPair(fRsnEventBuffer);
    }
  }
}

//________________________________________________________________________
void AliRsnAnalysisSE::PostEventProcess(const Short_t & index)
{
//=========================================================
// Post process of one event
//=========================================================

  switch (fInputType[index])
  {
    case kAOD:
      break;
    case kESD:
      break;
    case kESDMC:
      break;
    case kMC:
      break;
    case kRSN:
    {
      if (fRsnEventBuffer->GetDeleteBufferWhenReset() == kFALSE)
      {
        fRSN[index] = (AliRsnEvent*) fRsnEventBuffer->GetNextEvent();
        SetBranchAddress(0 , "rsnEvents", &fRSN[index]);
      }
      break;
    }
    default:
      break;
  }

}

void  AliRsnAnalysisSE::AddPairMgr(AliRsnPairMgr * pairmgr)
{
  fPairMgrs.Add(pairmgr);
}


void AliRsnAnalysisSE::AddPairMgrFromConfig(TString configfile)
{
  gROOT->LoadMacro(configfile.Data());

  configfile.ReplaceAll(".C","");

  AliRsnPairMgr *mgrRsn = (AliRsnPairMgr *) gROOT->ProcessLine(Form("%s();",configfile.Data()));
  if (!mgrRsn) return;

  fPairMgrs.Add(mgrRsn);
}
