
#include <TClonesArray.h>

#include "AliAnalysisNanoAODCutsJet.h"
#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "AliNanoAODHeader.h"
#include "AliNanoAODTrack.h"
#include "AliMultSelection.h"
#include "AliEventplane.h"
#include "AliAnalysisNanoAODCuts.h"
#include <iomanip>

ClassImp(AliNanoAODSimpleSetterJet)

AliNanoAODSimpleSetterJet::AliNanoAODSimpleSetterJet():
  AliNanoAODCustomSetter(),
  fArrayPythia(0x0),
  fArrayPythiaName("")
{
// default const
}
void AliNanoAODSimpleSetterJet::SetNanoAODTrack (const AliAODTrack * aodTrack, AliNanoAODTrack * spTrack){

    static  Int_t hybGlob  = AliNanoAODTrackMapping::GetInstance()->GetVarIndex("cstIsGlobalHybrid");
    static  Int_t inIsPyt  = AliNanoAODTrackMapping::GetInstance()->GetVarIndex("cstIsPythiaTrack");
    
    if (aodTrack-> IsHybridGlobalConstrainedGlobal()) spTrack->SetVar(hybGlob,1.);
    else spTrack->SetVar(hybGlob,0.);

    Int_t nTracksPythia = fArrayPythia->GetEntries();

    spTrack->SetVar(inIsPyt,0.);
    
    for(Int_t iTrackPyth = 0; iTrackPyth<nTracksPythia; iTrackPyth++){
        AliAODTrack *pythTrack = (AliAODTrack*)fArrayPythia->At(iTrackPyth);
	if((pythTrack->Pt()==aodTrack->Pt())&&(pythTrack->Eta()==aodTrack->Eta())&&(pythTrack->Phi()==aodTrack->Phi())) {
	  spTrack->SetVar(inIsPyt,1.);
	  return;
	}
    }


}

void AliNanoAODSimpleSetterJet::SetNanoAODHeader(const AliAODEvent * event   , AliNanoAODHeader * head , TString varListHeader  ) {

  AliNanoAODSimpleSetter simpleSetter;
  simpleSetter.SetNanoAODHeader(event,head,varListHeader);

  fArrayPythia = static_cast<TClonesArray*> (event->FindListObject(fArrayPythiaName.Data()));

  AliAODHeader *header = (AliAODHeader*)event->GetHeader();

  static Int_t indexEPV0 = head->GetVarIndex("cstEvPlaneV0");
  static Int_t indexEPV0A = head->GetVarIndex("cstEvPlaneV0A");
  static Int_t indexEPV0C = head->GetVarIndex("cstEvPlaneV0C");

  AliEventplane *evPlane = header->GetEventplaneP();
  Double_t EpV0 = evPlane->GetEventplane("V0" ,event);
  Double_t EpV0A = evPlane->GetEventplane("V0A",event);
  Double_t EpV0C = evPlane->GetEventplane("V0C",event);

  head->SetVar(indexEPV0,EpV0);
  head->SetVar(indexEPV0A,EpV0A);
  head->SetVar(indexEPV0C,EpV0C);

}
