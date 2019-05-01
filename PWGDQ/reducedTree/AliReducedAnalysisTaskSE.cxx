/*
 * **********************************************************
 * Virtual class for processing trees of AliReducedEventInfo
 * Authors: Ionut-Cristian Arsene, iarsene@cern.ch, i.c.arsene@fys.uio.no
 *                Jacobus Onderwaater, j.onderwaater@gsi.de, jacobus.onderwaater@cern.ch
 * Creation date: 2015/10/01
 *********************************************************
 */

#include "AliReducedAnalysisTaskSE.h"
#include "AliReducedEventInfo.h"

ClassImp(AliReducedAnalysisTaskSE);


//___________________________________________________________________________
AliReducedAnalysisTaskSE::AliReducedAnalysisTaskSE() :
  TObject(),
  fName(""),
  fTitle(""),
  fEvent(0x0),
  fFilteredTree(0x0),
  fActiveBranches(""),
  fInactiveBranches(""),
  fFilteredEvent(0x0),
  fFilteredTreeWritingOption(kBaseEventsWithBaseTracks),
  fProcessMCInfo(kFALSE),
  fEventCounter(0)
{
  //
  // default constructor
  //
  for(Int_t i=0; i<AliReducedVarManager::kNVars; ++i)
    fValues[i] = 0.0;
}


//___________________________________________________________________________
AliReducedAnalysisTaskSE::AliReducedAnalysisTaskSE(const Char_t* name, const Char_t* title) :
  TObject(),
  //fHistosManager(new AliHistogramManager("Histogram Manager", AliReducedVarManager::kNVars)),
  fName(name),
  fTitle(title),
  fEvent(0x0),
  fFilteredTree(0x0),
  fActiveBranches(""),
  fInactiveBranches(""),
  fFilteredEvent(0x0),
  fFilteredTreeWritingOption(kBaseEventsWithBaseTracks),
  fProcessMCInfo(kFALSE),
  fEventCounter(0)
{
  //
  // named constructor
  //
  for(Int_t i=0; i<AliReducedVarManager::kNVars; ++i)
    fValues[i] = 0.0;
}


//___________________________________________________________________________
AliReducedAnalysisTaskSE::~AliReducedAnalysisTaskSE() 
{
  //
  // destructor
  //
}

//___________________________________________________________________________
void AliReducedAnalysisTaskSE::InitFilteredTree() {
   //
   //
   //
   if(fFilteredTree) return; //already initialised
   fFilteredTree = new TTree("DstTree","Reduced ESD/AOD information");
   
   switch(fFilteredTreeWritingOption) {
      case kBaseEventsWithBaseTracks:
         fFilteredEvent = new AliReducedBaseEvent("DstEvent", AliReducedBaseEvent::kUseBaseTracks);
         break;
      case kBaseEventsWithFullTracks:
         fFilteredEvent = new AliReducedBaseEvent("DstEvent", AliReducedBaseEvent::kUseReducedTracks);
         break;
      case kFullEventsWithBaseTracks:
         fFilteredEvent = new AliReducedEventInfo("DstEvent", AliReducedBaseEvent::kUseBaseTracks);   
         break;
      case kFullEventsWithFullTracks:
         fFilteredEvent = new AliReducedEventInfo("DstEvent", AliReducedBaseEvent::kUseReducedTracks);   
         break;
      default:
         break;
   };
   
   fFilteredTree->Branch("Event",&fFilteredEvent,16000,99);
   
   // if user set active branches
   TObjArray* aractive=fActiveBranches.Tokenize(";");
   if(aractive->GetEntries()>0) {fFilteredTree->SetBranchStatus("*", 0);}
   for(Int_t i=0; i<aractive->GetEntries(); i++){
      fFilteredTree->SetBranchStatus(aractive->At(i)->GetName(), 1);
   }
   
   // if user set inactive branches
   TObjArray* arinactive=fInactiveBranches.Tokenize(";");
   for(Int_t i=0; i<arinactive->GetEntries(); i++){
      fFilteredTree->SetBranchStatus(arinactive->At(i)->GetName(), 0);
   }
}

//___________________________________________________________________________
void AliReducedAnalysisTaskSE::Init() {
   //
   // initialization (typically called in AliAnalysisTask::UserCreateOutputObjects())
   //
}

//___________________________________________________________________________
void AliReducedAnalysisTaskSE::Process() {
   //
   // process a given event (typically called in AliAnalysisTask::UserExec())
   //
}

//___________________________________________________________________________
void AliReducedAnalysisTaskSE::Finish() {
   //
   // finish, to be executed after all events were processed
   //
}
