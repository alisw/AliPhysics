#include "AliAnalysisTaskOnTheFlyAliAODDimuon.h"

#include "TClonesArray.h"
#include "AliAODDimuon.h"
#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "AliLog.h"
#include "AliCodeTimer.h"

#include <cassert>

//_________________________________________________________________________________________________
AliAnalysisTaskOnTheFlyAliAODDimuon::AliAnalysisTaskOnTheFlyAliAODDimuon() : AliAnalysisTaskSE("OnTheFlyAliAODDimuon"),
fDimuons(0x0)
{
  /// ctor
}

//_________________________________________________________________________________________________
AliAnalysisTaskOnTheFlyAliAODDimuon::~AliAnalysisTaskOnTheFlyAliAODDimuon()
{
  /// dtor
  /// Here we assume the ownership of fDimuons is now with (some) AODEvent, so we
  /// do _not_ delete it here
}

//_________________________________________________________________________________________________
void AliAnalysisTaskOnTheFlyAliAODDimuon::UserExec(Option_t* opt)
{
  /// Loop over muon tracks in the AOD event and build dimuon objects from them
  
  AliCodeTimerAuto("",0);
  
  AliAODEvent* event = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!event)
  {
    AliError("Could not get input event !!!");
    return;
  }
  
  if ( !event->GetDimuons() )
  {
	fDimuons = new TClonesArray("AliAODDimuon");
	AliInfo(Form("(Re)Creating on-the-fly dimuons array (%p)",fDimuons));
	fDimuons->SetName("dimuons"); // this name can not be random. See AliAODEvent::CreateStdContent
	fDimuons->SetOwner(kTRUE);
    AliInfo(Form("Adding on-the-fly dimuons array (%p) to AODEvent %p",fDimuons,event));
    event->AddObject(fDimuons); // by adding the object to the event, we transfer its ownership, so we won't take care of its deletion ourselves
    event->GetStdContent();
  }
  else
  {
    assert(fDimuons==event->GetDimuons());
  }
  
  if ( ! event->GetTracks() ) return;
  
  TClonesArray& tracks = *(event->GetTracks());

  Int_t nTracks = tracks.GetLast()+1;
  
  Int_t n(0);

  assert(fDimuons->GetLast()==-1);
  
  for ( Int_t i = 0; i < nTracks; ++i )
  {
    AliAODTrack* ti = static_cast<AliAODTrack*>(tracks[i]);
    if ( !ti->IsMuonTrack() ) continue;
    
    for ( Int_t j = i+1; j < nTracks; ++j )
    {
      AliAODTrack* tj = static_cast<AliAODTrack*>(tracks[j]);
      if ( !tj->IsMuonTrack() ) continue;

      new((*fDimuons)[n++]) AliAODDimuon(ti,tj);
    }
  }
}
