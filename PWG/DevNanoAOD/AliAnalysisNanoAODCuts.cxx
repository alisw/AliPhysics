#include "AliAnalysisNanoAODCuts.h"
#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "AliNanoAODHeader.h"
#include "AliNanoAODTrack.h"

ClassImp(AliAnalysisNanoAODTrackCuts)



AliAnalysisNanoAODTrackCuts::AliAnalysisNanoAODTrackCuts():
  AliAnalysisCuts(), fBitMask(1) 
{
  // default ctor 
}

Bool_t AliAnalysisNanoAODTrackCuts::IsSelected(TObject* obj)
{
  // Returns true if the track is good!
  AliAODTrack* track = dynamic_cast<AliAODTrack*>(obj);

  if (track->TestFilterBit(fBitMask)) return kTRUE;
  return kFALSE;  

}

AliAnalysisNanoAODEventCuts::AliAnalysisNanoAODEventCuts():
  AliAnalysisCuts(),  fVertexRange(10)
{
  // default ctor   
}

Bool_t AliAnalysisNanoAODEventCuts::IsSelected(TObject* obj)
{
  // Returns true if the object is a primary vertex
  AliAODEvent * evt = dynamic_cast<AliAODEvent*>(obj);
  
  AliAODVertex * vertex = evt->GetPrimaryVertex();
  if (!vertex) return kFALSE;
  if (TMath::Abs(vertex->GetZ()) > fVertexRange) return kFALSE;
  return kTRUE;
  
}
