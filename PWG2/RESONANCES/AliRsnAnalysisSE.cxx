//
// Class AliRsnAnalysisSE
//
// TODO
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#include <Riostream.h>

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
AliRsnAnalysisSE::AliRsnAnalysisSE(const char * name, Int_t bufferSize) :
  AliRsnAnalysisTaskSEBase(name),
  fDoesMixing(kFALSE),
  fMixingNum(0),
  fMixingCut(0x0),
  fPairMgrs(0),
  fOutList(0x0),
  fBuffer(0x0),
  fBufferSize(bufferSize)
{
//
// Default constructor
//

  InitIOVars();
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliRsnAnalysisSE::~AliRsnAnalysisSE()
{
//
// Destructor
//
}

//________________________________________________________________________
void AliRsnAnalysisSE::InitIOVars()
{
//
// Init input output values
//

  AliRsnAnalysisTaskSEBase::InitIOVars();

  fBuffer = 0;
  fOutList = 0;
}

//________________________________________________________________________
void AliRsnAnalysisSE::UserCreateOutputObjects()
{
//
// UserCreateOutputObjects() of AliAnalysisTaskSE
//

  OpenFile(0);
  fOutList = new TList();
  fOutList->SetOwner();
  AliRsnPair *def=0;
  AliRsnPairMgr *mgr=0;
  TList *listTmp;
  fDoesMixing = kFALSE;
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
        listTmp->Add(def->GenerateHistograms(mgr->GetName()));
        if (def->IsMixed()) fDoesMixing = kTRUE;
      }
    }
    fOutList->Add(listTmp);
  }

  TH1I *hUsed = new TH1I("hRsnUsed", "skipped and used events in this analysis", 2, 0, 2);
  fOutList->Add(hUsed);

  fBuffer = new AliRsnEventBuffer(fBufferSize);
}

//________________________________________________________________________
void AliRsnAnalysisSE::UserExec(Option_t *)
{
//
// UserExec() of AliAnalysisTaskSE
//

  static UInt_t eventID = 0;

  if (fEntry++ % 100 == 0) cout << "[" << GetName() << "] : processing entry " << fEntry-1 << endl;

  TH1I *h = (TH1I*)fOutList->FindObject("hRsnUsed");

  AliRsnEvent *curEvent = GetRsnEventFromInputType();
  if (curEvent) {
    curEvent->SetUniqueID(eventID++);
    ProcessEventAnalysis(curEvent);
    PostEventProcess();
    h->Fill(1);
  }
  else {
    h->Fill(0);
  }

  PostData(1, fOutList);
}

//________________________________________________________________________
void AliRsnAnalysisSE::Terminate(Option_t *)
{
//
// Terminate() of AliAnalysisTask
//

  fOutList = dynamic_cast<TList*>(GetOutputData(1));
  if (!fOutList) { AliError("At end of analysis, output list is NULL"); return; }
  //fOutList->Print();
}

//________________________________________________________________________
void AliRsnAnalysisSE::ProcessEventAnalysis(AliRsnEvent *curEvent)
{
//
// Process of one event
//

  // Adds event to Event Buffer
  fBuffer->AddEvent(curEvent);
  Int_t index = fBuffer->GetEventsBufferIndex();

  Int_t nmatches;
  TArrayI matched(0);
  if (fDoesMixing) matched = FindGoodMatches(index, nmatches);

  // loop over all Pair managers
  AliRsnPairMgr *mgr=0;
  for (Int_t iMgr=0 ;iMgr< fPairMgrs.GetEntries();iMgr++)
  {
    mgr = (AliRsnPairMgr *) fPairMgrs.At(iMgr);
    AliRsnPair *pair=0;
    for (Int_t i=0;i< mgr->GetPairs()->GetEntriesFast();i++)
    {
      pair = (AliRsnPair *) mgr->GetPairs()->At(i);
      if (!pair->IsMixed()) {
        pair->ProcessPair(curEvent, 0);
      }
      else {
        Int_t i, iev;
        for (i = 0; i < matched.GetSize(); i++) {
          iev = matched[i];
          if (iev < 0) continue;
          AliRsnEvent *evmatch = fBuffer->GetEvent(iev);
          pair->ProcessPair(curEvent, evmatch);
          if (!pair->IsPairEqual()) pair->ProcessPair(evmatch, curEvent);
        }
      }
    }
  }
}

//________________________________________________________________________
void AliRsnAnalysisSE::PostEventProcess(const Short_t & index)
{
//
// Post process of one event
//

  if (fInputType[index] != kRSN) return;

  if (fBuffer->GetDeleteBufferWhenReset() == kFALSE) {
    fRSN[index] = (AliRsnEvent*) fBuffer->GetNextEvent();
    SetBranchAddress(0 , "rsnEvents", &fRSN[index]);
  }
}

//________________________________________________________________________
void  AliRsnAnalysisSE::AddPairMgr(AliRsnPairMgr * pairmgr)
{
  fPairMgrs.Add(pairmgr);
}

//________________________________________________________________________
void AliRsnAnalysisSE::AddPairMgrFromConfig(TString configfile,TString analysisName)
{
  gROOT->LoadMacro(configfile.Data());

  configfile.ReplaceAll(".C","");

  analysisName.ReplaceAll("_","-");
  
  AliRsnPairMgr *mgrRsn = (AliRsnPairMgr *) gROOT->ProcessLine(Form("%s(\"%s\");", configfile.Data(),analysisName.Data()));
  if (!mgrRsn) return;

  fPairMgrs.Add(mgrRsn);
}

//________________________________________________________________________
TArrayI AliRsnAnalysisSE::FindGoodMatches(Int_t iRef, Int_t &foundMatches)
{
  // initialize the output array to the size of required mixed events
  // and initialize all members to -1
  Int_t i;
  TArrayI matched(fMixingNum);
  for (i = 0; i < fMixingNum; i++) matched[i] = -1;
  foundMatches = 0;

  // starts from the position behind the reference index
  // and goes backward; if it reaches the value 0, stops
  AliRsnEvent *refEvent = fBuffer->GetEvent(iRef);
  if (!refEvent) return matched;
  AliRsnEvent *matchEvent = 0x0;
  Int_t checkIndex;
  for (checkIndex = iRef - 1; ; checkIndex--) {
    if (checkIndex < 0) checkIndex = fBuffer->GetEventsBufferSize() - 1;
    if (checkIndex == iRef) break;
    matchEvent = fBuffer->GetEvent(checkIndex);
    if (!matchEvent) continue;
    if (fMixingCut) {
      if (!fMixingCut->IsSelected(AliRsnCut::kMixEvent, refEvent, matchEvent)) continue;
    }
    // assign to current array slot the matched event
    // and increment current slot and stops if it exceeds array size
    matched[foundMatches++] = checkIndex;
    if (foundMatches >= fMixingNum) break;
  }

  // returns the current index value,
  // which is also the number of matched events found
  return matched;
}
