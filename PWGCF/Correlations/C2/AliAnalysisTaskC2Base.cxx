#include <iostream>
#include <vector>
#include <algorithm>

#include "TH1F.h"
#include "TAxis.h"
#include "THn.h"
#include "TList.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TVectorF.h"

#include "AliAnalysisUtils.h"
#include "AliAODEvent.h"
#include "AliAODITSsaTrackCuts.h"
#include "AliAODMCParticle.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"
#include "AliLog.h"
#include "AliMCEvent.h"
#include "AliMultSelection.h"
#include "AliStack.h"
#include "AliVEvent.h"
#include "AliVVZERO.h"

#include "AliAnalysisTaskC2Base.h"
#include "AliAnalysisC2Utils.h"
#include "AliAnalysisC2Settings.h"

using std::cout;
using std::endl;
  
//________________________________________________________________________
AliAnalysisTaskC2Base::AliAnalysisTaskC2Base()
  : AliAnalysisTaskSE(),
    fOutputList(0),
    fitssatrackcuts(0),
    fDiscardedEvents(0),
    fDiscardedTracks(0),
    fSettings()
{
}

//________________________________________________________________________
AliAnalysisTaskC2Base::AliAnalysisTaskC2Base(const char *name)
  : AliAnalysisTaskSE(name),
    fOutputList(0),
    fitssatrackcuts(0),
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
  TAxis *discardedEvtsAx = this->fDiscardedEvents->GetXaxis();
  discardedEvtsAx->SetBinLabel(cDiscardEventReasons::_eventIsValid + 1, "event is valid");
  discardedEvtsAx->SetBinLabel(cDiscardEventReasons::badRun + 1, "bad run");
  discardedEvtsAx->SetBinLabel(cDiscardEventReasons::invalidxVertex + 1, "invalid vertex");
  discardedEvtsAx->SetBinLabel(cDiscardEventReasons::isIncomplete + 1, "incomplete");
  discardedEvtsAx->SetBinLabel(cDiscardEventReasons::isOutOfBunchPileup + 1, "out of bunch PU");
  discardedEvtsAx->SetBinLabel(cDiscardEventReasons::noMultSelectionObject + 1, "!multSelection");
  discardedEvtsAx->SetBinLabel(cDiscardEventReasons::noTracks + 1, "!trackArray");
  discardedEvtsAx->SetBinLabel(cDiscardEventReasons::noTracksInPtRegion + 1, "no track in pt interval");
  discardedEvtsAx->SetBinLabel(cDiscardEventReasons::noTrigger + 1, "no trigger fired");
  discardedEvtsAx->SetBinLabel(cDiscardEventReasons::notV0AND + 1, "not V0AND");
  discardedEvtsAx->SetBinLabel(cDiscardEventReasons::spdFastOr + 1, "SPD fast OR");
  discardedEvtsAx->SetBinLabel(cDiscardEventReasons::spdPipeup + 1, "SPD PU");
  discardedEvtsAx->SetBinLabel(cDiscardEventReasons::spdVertexContributors + 1, "SPD vertex contr");
  discardedEvtsAx->SetBinLabel(cDiscardEventReasons::tklClusterCut + 1, "tkl cluster cut");
  discardedEvtsAx->SetBinLabel(cDiscardEventReasons::v0asymmetryCut + 1, "V0 asymmetry cut");
  discardedEvtsAx->SetBinLabel(cDiscardEventReasons::zvtxPosition + 1, "zvtx position");
  this->fOutputList->Add(this->fDiscardedEvents);

  this->fDiscardedTracks = new TH1F("discardedTracks", "discardedTracks",
				    cDiscardTrackReasons::nDiscardTrackReasons,
				    0,
				    cDiscardTrackReasons::nDiscardTrackReasons);
  TAxis *discardedTracksAx = this->fDiscardedTracks->GetXaxis();
  discardedTracksAx->SetBinLabel(cDiscardTrackReasons::_trackIsValid + 1, "Track is valid");
  discardedTracksAx->SetBinLabel(cDiscardTrackReasons::dca + 1, "DCA");
  discardedTracksAx->SetBinLabel(cDiscardTrackReasons::etaAcceptance + 1, "#eta acceptance");
  discardedTracksAx->SetBinLabel(cDiscardTrackReasons::failedFilterBits + 1, "failed filter bits");
    discardedTracksAx->SetBinLabel(cDiscardTrackReasons::failedITSCut + 1, "failed ITS cut");
  discardedTracksAx->SetBinLabel(cDiscardTrackReasons::neutralCharge + 1, "neutral");
  discardedTracksAx->SetBinLabel(cDiscardTrackReasons::notAODPrimary + 1, "not aod kPrimary");
  discardedTracksAx->SetBinLabel(cDiscardTrackReasons::notHybridGCG + 1, "!HybridGCG");
  discardedTracksAx->SetBinLabel(cDiscardTrackReasons::notMCPrimary + 1, "not MC PhysicalPrimary");
  this->fOutputList->Add(this->fDiscardedTracks);

  this->fitssatrackcuts = AliAODITSsaTrackCuts::GetStandardAODITSsaTrackCuts2015();

  AliLog::SetGlobalLogLevel(AliLog::kError);
  PostData(1, fOutputList);
}

Bool_t AliAnalysisTaskC2Base::IsValidEvent()
{
  // skip runs with long tails in V0C
  std::vector< Int_t > badRuns = {225611, 225609, 225589};
  
  if (!this->InputEvent()->GetPrimaryVertex() || !this->InputEvent()->GetPrimaryVertex()->GetZ()) {
    this->fDiscardedEvents->Fill(cDiscardEventReasons::invalidxVertex);
    return false;
  }
  if (this->InputEvent()->GetPrimaryVertex()->GetZ() < this->fSettings.fZVtxAcceptanceLowEdge
      || this->InputEvent()->GetPrimaryVertex()->GetZ() > this->fSettings.fZVtxAcceptanceUpEdge) {
    this->fDiscardedEvents->Fill(cDiscardEventReasons::zvtxPosition);
    return false;
  }
  if(std::find(badRuns.begin(), badRuns.end(), InputEvent()->GetRunNumber()) != badRuns.end()) {
    this->fDiscardedEvents->Fill(cDiscardEventReasons::badRun);
    return false;
  }
  if ((InputEvent()->GetVZEROData()->GetV0CDecision() != 1) ||
	   (InputEvent()->GetVZEROData()->GetV0ADecision() != 1)) {
    this->fDiscardedEvents->Fill(cDiscardEventReasons::notV0AND);
    return false;
  }
  //incomplete events
  if (InputEvent()->GetHeader()->GetL0TriggerInputs() == 0 && !fMCEvent) {
    this->fDiscardedEvents->Fill(cDiscardEventReasons::isIncomplete);
    return false;
  }
  //V0 asymmetry cut
  if (this->IsAsymmetricV0()) {
    this->fDiscardedEvents->Fill(cDiscardEventReasons::v0asymmetryCut);
    return false;
  }
  //trigger check; TODO: ask Katarina about "fMBtrigger" in her code!
  if (TString firedTrigger = InputEvent()->GetFiredTriggerClasses()) {
    // Make choice base on type of MB trigger?
    if (!firedTrigger.Contains(this->fSettings.fOfflineTrigger)) {
      this->fDiscardedEvents->Fill(cDiscardEventReasons::noTrigger);
      return false;
    }
  }
  //SPD vtx contributors (only aod events)
  if (AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(fInputEvent)) {
    if (aodEvent->GetPrimaryVertexSPD()->GetNContributors() < 1) {
      this->fDiscardedEvents->Fill(cDiscardEventReasons::spdVertexContributors);
      return false;
    }
  }
  //out-of-bunch 11 BC
  if (IsOutOfBunchPileup()) {
    this->fDiscardedEvents->Fill(cDiscardEventReasons::isOutOfBunchPileup);
    return false;
  }
  // Online-offline SPD fastor
  // Taken from Leonardo/Katarina; TODO: Cross-check and understand these values!
  if (AliVMultiplicity* mult = fInputEvent->GetMultiplicity()) {
    if(mult->GetFastOrFiredChipMap().CountBits(400)
       <= -20.589 + 0.73664 * mult->GetFiredChipMap().CountBits(400)) {
      this->fDiscardedEvents->Fill(cDiscardEventReasons::spdFastOr);
      return false;
    }
  }
  // SPD pileup
  if (AliAnalysisUtils* utils = new AliAnalysisUtils()) {
    utils->SetMinPlpContribSPD(3);
    utils->SetMinPlpContribMV(3);
    if (utils->IsPileUpSPD(InputEvent())) {
      this->fDiscardedEvents->Fill(cDiscardEventReasons::spdPipeup);
      return false;
    }
  }
  AliMultSelection *multSel =
    dynamic_cast< AliMultSelection* >(InputEvent()->FindListObject("MultSelection"));
  // no multSelection object found
  if (!multSel){
    this->fDiscardedEvents->Fill(cDiscardEventReasons::noMultSelectionObject);
    return false;
  }
  // tkl-cluster cut
  Float_t multTKL = multSel->GetEstimator("SPDTracklets")->GetValue();
  Int_t nITScluster = InputEvent()->GetNumberOfITSClusters(0) + InputEvent()->GetNumberOfITSClusters(1);
  if (nITScluster > 64+4*multTKL) {
    this->fDiscardedEvents->Fill(cDiscardEventReasons::tklClusterCut);
    return false;
  }
  this->fDiscardedEvents->Fill(cDiscardEventReasons::_eventIsValid);
  return true;
}

Bool_t AliAnalysisTaskC2Base::IsValidParticle(AliVParticle *particle) {
  AliAODTrack* aodTrack = dynamic_cast< AliAODTrack* >(particle);
  AliAODMCParticle* mcParticle = dynamic_cast< AliAODMCParticle* >(particle);

  if (particle->Charge() == 0){
    this->fDiscardedTracks->Fill(cDiscardTrackReasons::neutralCharge);
    return false;
  }
  if (particle->Eta() < this->fSettings.fEtaAcceptanceLowEdge ||
      particle->Eta() > this->fSettings.fEtaAcceptanceUpEdge) {
    this->fDiscardedTracks->Fill(cDiscardTrackReasons::etaAcceptance);
    return false;
  }
  // MC specific cuts
  if (this->fSettings.kMCTRUTH == this->fSettings.fDataType) {
    if (!mcParticle->IsPhysicalPrimary()) {
      this->fDiscardedTracks->Fill(cDiscardTrackReasons::notMCPrimary);
      return false;
    }
  }
  // Reconstruction specific cuts
  if (this->fSettings.kRECON == this->fSettings.fDataType) {
    if(!this->fitssatrackcuts->AcceptTrack(aodTrack)) {
      this->fDiscardedTracks->Fill(cDiscardTrackReasons::failedITSCut);
      return false;
    }
    // if (!aodTrack->IsHybridGlobalConstrainedGlobal()) {
    //   this->fDiscardedTracks->Fill(cDiscardTrackReasons::notHybridGCG);
    //   return false;
    // }
    if (aodTrack->GetType() != aodTrack->kPrimary) {
      this->fDiscardedTracks->Fill(cDiscardTrackReasons::notAODPrimary);
      return false;
    }
    {
      Double_t dcaTang;
      Double_t dcaLong;
      AliAnalysisC2Utils::GetDCA(dcaTang, dcaLong, aodTrack);
      if (TMath::Abs(dcaTang) > this->fSettings.fMaxDcaTang ||
	  TMath::Abs(dcaLong) > this->fSettings.fMaxDcaLong){
	this->fDiscardedTracks->Fill(cDiscardTrackReasons::dca);
	return false;
      }
    }
  }
  this->fDiscardedTracks->Fill(cDiscardTrackReasons::_trackIsValid);
  return true;
}

void AliAnalysisTaskC2Base::SetupEventForBase() {
    // Pass the current event's vertext to the ITS cut class
    fitssatrackcuts->ExtractAndSetPrimaryVertex(this->InputEvent());
}

Bool_t AliAnalysisTaskC2Base::IsAsymmetricV0() {
  // Lambda: Get the multiplicyt summed over all sectors of a V0C ring
  AliVVZERO* vzero = fInputEvent->GetVZEROData();
  auto getV0CRingMult = [vzero](Int_t iring) -> Float_t {
    Float_t mult = 0;
    for (int isector=0; isector<8; isector++) {
	mult += vzero->GetMultiplicityV0C(iring*8 + isector);
      }
    return mult;
  };
  Float_t multTotalV0A = vzero->GetMTotV0A();
  Float_t multTotalV0C = vzero->GetMTotV0C();
  Float_t multV0C_ring012 = getV0CRingMult(0) + getV0CRingMult(1) + getV0CRingMult(2);
  Float_t multV0C_ring3 = getV0CRingMult(3);

  if (multTotalV0C < (330 + 100 * TMath::Power(multTotalV0A, 0.2))) // V0A vs V0C asymmetry
    {
      if (multV0C_ring012 < 160 || //reject high activity in V0C012
	  multV0C_ring3 > 12*TMath::Power(0.01*(multV0C_ring012 - 160), 1.7)) // asymmetry V0C3 vs V0C012
	{
	  return false;
	}
    }
  return true;
}

Bool_t AliAnalysisTaskC2Base::IsOutOfBunchPileup(){
   // do not apply this cut in LHC15f - use spdOnline-spdOffline cut instead
  if(InputEvent()->GetRunNumber() < 227730)
    return false;

  std::vector< Int_t > irInt1BitNumbers = {87, 88, 89, 93, 94, 95, 96, 97, 98, 99, 100, 101};
  for(auto const& bitNumber: irInt1BitNumbers) {
    if (InputEvent()->GetHeader()->GetIRInt1InteractionMap().TestBitNumber(bitNumber)) {
      return true;
    }
  }
  return false;
}
