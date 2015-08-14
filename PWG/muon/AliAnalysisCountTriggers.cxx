#include "AliAnalysisCountTriggers.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TString.h"
#include "AliCounterCollection.h"
#include "AliAnalysisManager.h"
#include "AliLog.h"
#include "AliESDEvent.h"

ClassImp(AliAnalysisCountTriggers)

//_____________________________________________________________________________
AliAnalysisCountTriggers::AliAnalysisCountTriggers()
: AliAnalysisTaskSE("AliAnalysisCountTriggers"),
fEventCounters(0),
fHTriggerMask(0)
{
  /// default ctor. This task has two output objects : one counter collection
  /// and one histogram
  DefineOutput(1,AliCounterCollection::Class());
  DefineOutput(2,TH1I::Class());
}

//_____________________________________________________________________________
AliAnalysisCountTriggers::~AliAnalysisCountTriggers()
{
  /// dtor

  if (fEventCounters && ! AliAnalysisManager::GetAnalysisManager()->IsProofMode())
  {
    delete fEventCounters;
  }
  if (fHTriggerMask && ! AliAnalysisManager::GetAnalysisManager()->IsProofMode())
  {
    delete fHTriggerMask;
  }
}

//_____________________________________________________________________________
AliVEvent*
AliAnalysisCountTriggers::Event() const
{
  // some const-dirty-dancing
  return const_cast<AliAnalysisCountTriggers*>(this)->InputEvent();
}

//_____________________________________________________________________________
void
AliAnalysisCountTriggers::Terminate(Option_t *)
{
  /// Called once at the end of the query
  /// Just a simple printout of the stat we analysed
  
  fEventCounters = dynamic_cast<AliCounterCollection*>(GetOutputData(1));
  
  if (!fEventCounters)
  {
    AliError("Could not find back counters in output...");
  }
  else
  {
    fEventCounters->Print("trigger/run");
  }
  
}

//_____________________________________________________________________________
void AliAnalysisCountTriggers::UserExec(Option_t* /*opt*/)
{
  /// Executed at each event
  
  TString firedTriggerClasses(Event()->GetFiredTriggerClasses());

  TObjArray* a = firedTriggerClasses.Tokenize(" ");
  TIter next(a);
  TObjString* s;

  // ANY will give the total number of analyzed events
  fEventCounters->Count(Form("trigger:ANY/run:%d",Event()->GetRunNumber()));

  while ( ( s = static_cast<TObjString*>(next())) )
  {
    fEventCounters->Count(Form("trigger:%s/run:%d",s->String().Data(),Event()->GetRunNumber()));
  }
  
  delete a;
  
  // fill the trigger mask histogram. One bin per trigger bit.
  // will obviously not work if running on several runs which do
  // not share a common trigger configuration...
  
  if ( Event()->IsA() == AliESDEvent::Class() )
  {
    AliESDEvent* e = static_cast<AliESDEvent*>(Event());
    
    for ( ULong64_t i = 0; i < 50; ++i )
    {
      ULong64_t low = ( 1ULL << i );
      ULong64_t high = ( 1ULL << (i+50) );
      
      if ( e->GetTriggerMask() & low ) { fHTriggerMask->Fill(i*1.0); }
      if ( e->GetTriggerMaskNext50() & high ) { fHTriggerMask->Fill(i*1.0+50.0); }
      
      
    }
  }

  // Post output data.
  PostData(1, fEventCounters);
  PostData(2, fHTriggerMask);
}

//_____________________________________________________________________________
void AliAnalysisCountTriggers::UserCreateOutputObjects()
{
  /// Create histograms
  /// Called once
  
  OpenFile(1);
  
  fEventCounters = new AliCounterCollection("CC");

  // initialize event counters

  fEventCounters->AddRubric("trigger", 100);
  
  fEventCounters->AddRubric("run", 1000000);
  
  // finally end the counters initialization
  fEventCounters->Init();
  
  fHTriggerMask = new TH1I("hTriggerMask","hTriggerMask",100,-0.5,99.5);

  // Post output data.
  PostData(1,fEventCounters);
  PostData(2,fHTriggerMask);
}
