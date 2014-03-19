#include "AliAnalysisNanoAODCuts.h"
#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "AliNanoAODHeader.h"
#include "AliNanoAODTrack.h"
#include <iomanip>

ClassImp(AliAnalysisNanoAODTrackCuts)
ClassImp(AliAnalysisNanoAODEventCuts)
ClassImp(AliNanoAODSimpleSetter)


AliAnalysisNanoAODTrackCuts::AliAnalysisNanoAODTrackCuts():
AliAnalysisCuts(), fBitMask(1), fMinPt(0), fMaxEta(10)
{
  // default ctor 
}

Bool_t AliAnalysisNanoAODTrackCuts::IsSelected(TObject* obj)
{
  // Returns true if the track is good!
  AliAODTrack* track = dynamic_cast<AliAODTrack*>(obj);
  
  
  if (!track->TestFilterBit(fBitMask))    return kFALSE;
  if (track->Pt() < fMinPt)               return kFALSE;
  if (TMath::Abs(track->Eta()) > fMaxEta) return kFALSE; 
  
  return kTRUE;  

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


void AliNanoAODSimpleSetter::SetNanoAODHeader(const AliAODEvent * event   , AliNanoAODHeader * head  ) {
  // Set custom nano aod vars
  Double_t centr    = event->GetHeader()->GetCentralityP()->GetCentralityPercentile("V0M");
  Double_t magfield = event->GetHeader()->GetMagneticField();
  head->SetVar(0, centr);
  head->SetVar(1, magfield);

}
