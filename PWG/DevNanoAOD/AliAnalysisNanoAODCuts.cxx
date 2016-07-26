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
  AliAnalysisCuts(), 
  fVertexRange(-1),
  fTrackCut(0),
  fMinMultiplicity(-1),
  fMaxMultiplicity(-1)
{
  // default ctor   
}

Bool_t AliAnalysisNanoAODEventCuts::IsSelected(TObject* obj)
{
  // Returns true if object accepted on the event level
  
  AliAODEvent * evt = dynamic_cast<AliAODEvent*>(obj);
  
  if (fVertexRange > 0)
  {
    AliAODVertex * vertex = evt->GetPrimaryVertex();
    if (!vertex) 
      return kFALSE;
    
    if (vertex->GetNContributors() < 1) 
    {
      // SPD vertex cut
      vertex = evt->GetPrimaryVertexSPD();    
      if (!vertex || vertex->GetNContributors() < 1) 
        return kFALSE;
    }    
    
    if (TMath::Abs(vertex->GetZ()) > fVertexRange) 
      return kFALSE;
  }
  
  if (fTrackCut != 0)
  {
    Int_t trackCount = 0;
    for (Int_t j=0; j<evt->GetNumberOfTracks(); j++)
      if (fTrackCut->IsSelected(evt->GetTrack(j)))
        trackCount++;
      
    if (fMinMultiplicity > 0 && trackCount < fMinMultiplicity)
      return kFALSE;
    if (fMaxMultiplicity > 0 && trackCount > fMaxMultiplicity)
      return kFALSE;
  }
      
  return kTRUE;
}


void AliNanoAODSimpleSetter::SetNanoAODHeader(const AliAODEvent * event   , AliNanoAODHeader * head  ) {

  AliAODHeader * header = dynamic_cast<AliAODHeader*>(event->GetHeader());
  if (!header) AliFatal("Not a standard AOD");
  // Set custom nano aod vars
  Double_t centr    = header->GetCentralityP()->GetCentralityPercentile("V0M");
  Double_t magfield = header->GetMagneticField();
  head->SetVar(0, centr);
  head->SetVar(1, magfield);

}
