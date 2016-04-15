#include <iostream>
#include <vector>

#include "THnSparse.h"
#include "THn.h"
#include "TList.h"
#include "TMath.h"
#include "TVectorF.h"
#include "TMatrix.h"
#include "TRandom3.h"

#include "AliAODEvent.h"
#include "AliAODMCParticle.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"
#include "AliLog.h"
#include "AliMCEvent.h"
#include "AliMultSelection.h"
#include "AliStack.h"
#include "AliVEvent.h"

#include "AliAnalysisTaskC2Base.h"
#include "AliAnalysisC2Utils.h"
#include "AliAnalysisC2Settings.h"

using std::cout;
using std::endl;
  
//________________________________________________________________________
AliAnalysisTaskC2Base::AliAnalysisTaskC2Base()
  : AliAnalysisTaskSE(),
    fOutputList(0),
    fDiscardedEvents(0),
    fDiscardedTracks(0),
    fSettings()
{
}

//________________________________________________________________________
AliAnalysisTaskC2Base::AliAnalysisTaskC2Base(const char *name)
  : AliAnalysisTaskSE(name),
    fOutputList(0),
    fDiscardedEvents(0),
    fDiscardedTracks(0),
    fSettings()
{
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
void AliAnalysisTaskC2Base::UserCreateOutputObjects()
{
  // Setup output list
  this->fOutputList = new TList();
  this->fOutputList->SetOwner();
  this->fDiscardedEvents = new TH1F("discardedEvents", "discardedEvents",
				    cDiscardEventReasons::nDiscardEventReasons,
				    0,
				    cDiscardEventReasons::nDiscardEventReasons);
  this->fDiscardedEvents->GetXaxis()->SetBinLabel(cDiscardEventReasons::zvtxPosition + 1,
						  "zvtx position");
  this->fDiscardedEvents->GetXaxis()->SetBinLabel(cDiscardEventReasons::invalidxVertex + 1,
						  "invalid vertex");
  this->fDiscardedEvents->GetXaxis()->SetBinLabel(cDiscardEventReasons::noTracks + 1,
						  "no Track array");
  this->fDiscardedEvents->GetXaxis()->SetBinLabel(cDiscardEventReasons::noTracksInPtRegion + 1,
						  "no track in pt interval");
  this->fOutputList->Add(this->fDiscardedEvents);

  this->fDiscardedTracks = new TH1F("discardedTracks", "discardedTracks",
				    cDiscardTrackReasons::nDiscardTrackReasons,
				    0,
				    cDiscardTrackReasons::nDiscardTrackReasons);
  this->fDiscardedTracks->GetXaxis()->SetBinLabel(cDiscardTrackReasons::dca + 1,
						  "DCA");
  this->fDiscardedTracks->GetXaxis()->SetBinLabel(cDiscardTrackReasons::etaAcceptance + 1,
						  "#eta acceptance");
  this->fDiscardedTracks->GetXaxis()->SetBinLabel(cDiscardTrackReasons::neutralCharge + 1,
						  "neutral");
  this->fDiscardedTracks->GetXaxis()->SetBinLabel(cDiscardTrackReasons::notAODPrimary + 1,
						  "!kPrimary");
  this->fDiscardedTracks->GetXaxis()->SetBinLabel(cDiscardTrackReasons::notHybridGCG + 1,
						  "!HybridGCG");
  this->fDiscardedTracks->GetXaxis()->SetBinLabel(cDiscardTrackReasons::notMCPrimary + 1,
						  "!PhysicalPrimary");
  this->fOutputList->Add(this->fDiscardedTracks);
  cout << this->fOutputList << endl;
  AliLog::SetGlobalLogLevel(AliLog::kError);
  PostData(1, fOutputList);
}

Bool_t AliAnalysisTaskC2Base::IsValidEvent()
{
  Double_t zvtxMax = 10;
  if (!this->InputEvent()->GetPrimaryVertex() ||
      !this->InputEvent()->GetPrimaryVertex()->GetZ())
    {
      this->fDiscardedEvents->Fill(cDiscardEventReasons::invalidxVertex);
      return false;
    }
  else if (TMath::Abs(this->InputEvent()->GetPrimaryVertex()->GetZ()) > zvtxMax)
    {
      this->fDiscardedEvents->Fill(cDiscardEventReasons::zvtxPosition);
      return false;
    }
  return true;
}

Bool_t AliAnalysisTaskC2Base::IsValidParticle(AliVParticle *particle)
{
  AliAODTrack* aodTrack = dynamic_cast< AliAODTrack* >(particle);
  AliAODMCParticle* mcParticle = dynamic_cast< AliAODMCParticle* >(particle);
  if (particle->Charge() == 0)
    return false;

  // restric eta region; assuming symmetric eta range and equal range for all histograms
  // TODO: This should be done more transparently, instead of just using one of the hists!
  else if (particle->Eta() < this->fSettings.fEtaAcceptanceLowEdge ||
	   particle->Eta() > this->fSettings.fEtaAcceptanceUpEdge) {
    this->fDiscardedTracks->Fill(cDiscardTrackReasons::etaAcceptance);
    return false;
  }
  if ((this->fSettings.kMCTRUTH == this->fSettings.fDataType) && !(mcParticle->IsPhysicalPrimary())) {
    this->fDiscardedTracks->Fill(cDiscardTrackReasons::notMCPrimary);
    return false;
  }
  if ((this->fSettings.kRECON == this->fSettings.fDataType) && !aodTrack->IsHybridGlobalConstrainedGlobal()){
    this->fDiscardedTracks->Fill(cDiscardTrackReasons::notHybridGCG);
    return false;
  }
  if ((this->fSettings.kRECON == this->fSettings.fDataType) && (aodTrack->GetType() != aodTrack->kPrimary)){
    this->fDiscardedTracks->Fill(cDiscardTrackReasons::notAODPrimary);
    return false;
  }
  if (this->fSettings.kRECON == this->fSettings.fDataType) {
    Double_t dcaTang;
    Double_t dcaLong;
    AliAnalysisC2Utils::GetDCA(dcaTang, dcaLong, aodTrack);
    if (TMath::Abs(dcaTang) > this->fSettings.fMaxDcaTang ||
	TMath::Abs(dcaLong) > this->fSettings.fMaxDcaLong){
      this->fDiscardedTracks->Fill(cDiscardTrackReasons::dca); \
      return false;
    }
  }
  return true;
}
