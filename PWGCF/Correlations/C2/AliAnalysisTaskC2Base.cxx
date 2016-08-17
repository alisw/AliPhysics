#include <iostream>
#include <vector>
#include <algorithm>

#include "TH1F.h"
#include "TH2F.h"
#include "TAxis.h"
#include "THn.h"
#include "TList.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TVectorF.h"

#include "AliAODEvent.h"
#include "AliAODForwardMult.h"
#include "AliAODITSsaTrackCuts.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisUtils.h"
#include "AliGenEventHeader.h"
#include "AliHeader.h"
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
    fCachedValues(),
    fDiscardedEvents(0),
    fDiscardedTracks(0),
    fOutputList(0),
    fSettings(),
    fValidTracks(0),
    fitssatrackcuts(0),
    fmultDistribution(0),
    fetaVsZvtx(0)
{
}

//________________________________________________________________________
AliAnalysisTaskC2Base::AliAnalysisTaskC2Base(const char *name)
  : AliAnalysisTaskSE(name),
    fCachedValues(),
    fDiscardedEvents(0),
    fDiscardedTracks(0),
    fOutputList(0),
    fSettings(),
    fValidTracks(0),
    fitssatrackcuts(0),
    fmultDistribution(0),
    fetaVsZvtx(0)
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
  discardedEvtsAx->SetBinLabel(cDiscardEventReasons::invalidxVertex + 1, "invalid vertex");
  discardedEvtsAx->SetBinLabel(cDiscardEventReasons::isIncomplete + 1, "incomplete");
  discardedEvtsAx->SetBinLabel(cDiscardEventReasons::isOutOfBunchPileup + 1, "out of bunch PU");
  discardedEvtsAx->SetBinLabel(cDiscardEventReasons::multEstimatorNotAvailable + 1, "!multEstimator");
  discardedEvtsAx->SetBinLabel(cDiscardEventReasons::noMultSelectionObject + 1, "!multSelection");
  discardedEvtsAx->SetBinLabel(cDiscardEventReasons::noForwardMultObj + 1, "!ForwardMult");
  discardedEvtsAx->SetBinLabel(cDiscardEventReasons::noEntriesInFMD + 1, "!Entries in FMD");
  discardedEvtsAx->SetBinLabel(cDiscardEventReasons::MeanMult0 + 1, "MeanMult0");
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

  this->fmultDistribution = new TH1F("multDistribution", "multDistribution", 210, 0, 210);
  this->fOutputList->Add(this->fmultDistribution);

  this->fetaVsZvtx = new TH2F("etaVsZvtx", "etaVsZvtx",
  			      200,
  			      -4,
  			      6,
			      100,
  			      this->fSettings.fZVtxAcceptanceLowEdge,
  			      this->fSettings.fZVtxAcceptanceUpEdge
			      );
  this->fOutputList->Add(this->fetaVsZvtx);
  
  this->fitssatrackcuts = AliAODITSsaTrackCuts::GetStandardAODITSsaTrackCuts2015();

  AliLog::SetGlobalLogLevel(AliLog::kError);
  PostData(1, fOutputList);
}

Bool_t AliAnalysisTaskC2Base::IsValidEvent()
{
  AliAODForwardMult* aodForward =
    dynamic_cast<AliAODForwardMult*>(fInputEvent->FindListObject("Forward"));
  if (!aodForward) {
    this->fDiscardedEvents->Fill(cDiscardEventReasons::noForwardMultObj);
    return false;
  }
  if (aodForward->GetHistogram().GetEntries() == 0) {
    this->fDiscardedEvents->Fill(cDiscardEventReasons::noEntriesInFMD);
    return false;
  }
  
  if (!this->InputEvent()->GetPrimaryVertex() || !this->InputEvent()->GetPrimaryVertex()->GetZ()) {
    this->fDiscardedEvents->Fill(cDiscardEventReasons::invalidxVertex);
    return false;
  }
  if (this->InputEvent()->GetPrimaryVertex()->GetZ() < this->fSettings.fZVtxAcceptanceLowEdge
      || this->InputEvent()->GetPrimaryVertex()->GetZ() > this->fSettings.fZVtxAcceptanceUpEdge) {
    this->fDiscardedEvents->Fill(cDiscardEventReasons::zvtxPosition);
    return false;
  }
  if (!this->GetAllTracks() || this->GetAllTracks()->GetSize() == 0) {
    this->fDiscardedEvents->Fill(cDiscardEventReasons::noTracks);
    return false;
  }

  // are we using a multSelector from the multSelection framework? Is it present?
  if (this->fSettings.fMultEstimator != this->fSettings.fMultEstimatorValidTracks) {
    AliMultSelection *multSel =
      dynamic_cast< AliMultSelection* >(InputEvent()->FindListObject("MultSelection"));
    // no multSelection object found
    if (!multSel){
      this->fDiscardedEvents->Fill(cDiscardEventReasons::noMultSelectionObject);
      return false;
    }
    if (!multSel->GetEstimator(this->fSettings.fMultEstimator)){
      this->fDiscardedEvents->Fill(cDiscardEventReasons::multEstimatorNotAvailable);
      return false;
    }
    if (multSel->GetEstimator(this->fSettings.fMultEstimator)->GetMean() == 0){
      this->fDiscardedEvents->Fill(cDiscardEventReasons::MeanMult0);
      return false;
    }
  }

  // ITSsa specific checks
  if (this->fSettings.fIsITSsa) {
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
      if (!firedTrigger.Contains(this->fSettings.fTrigger)) {
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
    // tkl-cluster cut
    // Float_t multTKL = multSel->GetEstimator("SPDTracklets")->GetValue();
    // Int_t nITScluster = (InputEvent()->GetNumberOfITSClusters(0)
    // 			 + InputEvent()->GetNumberOfITSClusters(1));
    // if (nITScluster > 64+4*multTKL) {
    //   this->fDiscardedEvents->Fill(cDiscardEventReasons::tklClusterCut);
    //   return false;
    // }
  }
  this->fDiscardedEvents->Fill(cDiscardEventReasons::_eventIsValid);
  return true;
}

Bool_t AliAnalysisTaskC2Base::IsValidParticle(AliVParticle *particle) {
  if (particle->Charge() == 0){
    this->fDiscardedTracks->Fill(cDiscardTrackReasons::neutralCharge);
    return false;
  }
  // if (particle->Eta() < this->fSettings.fEtaAcceptanceLowEdge ||
  //     particle->Eta() > this->fSettings.fEtaAcceptanceUpEdge) {
  //   this->fDiscardedTracks->Fill(cDiscardTrackReasons::etaAcceptance);
  //   return false;
  // }
  // MC specific cuts
  if (this->fSettings.kMCTRUTH == this->fSettings.fDataType) {
    AliAODMCParticle* mcParticle = dynamic_cast< AliAODMCParticle* >(particle);
    if (!mcParticle->IsPhysicalPrimary()) {
      this->fDiscardedTracks->Fill(cDiscardTrackReasons::notMCPrimary);
      return false;
    }
  }
  // Reconstruction specific cuts
  if (this->fSettings.kRECON == this->fSettings.fDataType) {
    AliAODTrack* aodTrack = dynamic_cast< AliAODTrack* >(particle);
    // ITSsa specific cuts
    if(this->fSettings.fIsITSsa) {
      if (!this->fitssatrackcuts->AcceptTrack(aodTrack)) {
	this->fDiscardedTracks->Fill(cDiscardTrackReasons::failedITSCut);
	return false;
      }
    }
    // Normal (ie., not ITSsa cuts)
    if (!this->fSettings.fIsITSsa) {
      // IsHybridGlobalConstrainedGlobal() seems to be equivalent to TestFilterBit(768) == false
      // at least on LHC10h AOD160
      // if (!aodTrack->IsHybridGlobalConstrainedGlobal()) {
      //   this->fDiscardedTracks->Fill(cDiscardTrackReasons::notHybridGCG);
      //   return false;
      // }
      if (aodTrack->TestFilterBit(768) == false) {
	this->fDiscardedTracks->Fill(cDiscardTrackReasons::failedFilterBits);
	return false;
      }
    }

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
  // reset event status
  for (Int_t i = 0; i < cCachedValues::nCachedValues; i++) {
    this->fCachedValues[i] = false;
  }

  // Pass the current event's vertext to the ITS cut class
  fitssatrackcuts->ExtractAndSetPrimaryVertex(this->InputEvent());
}

TClonesArray* AliAnalysisTaskC2Base::GetAllTracks() {
  // Yes, the following is realy "aodEvent" not mcEvent :P
  AliAODEvent* aodEvent = dynamic_cast< AliAODEvent* >(this->InputEvent());
  TClonesArray* tracksArray = (this->fSettings.kMCTRUTH == this->fSettings.fDataType)
    ? dynamic_cast<TClonesArray*>(aodEvent->GetList()->FindObject(AliAODMCParticle::StdBranchName()))
    : aodEvent->GetTracks();
  return tracksArray;
}

std::vector< AliAnalysisC2NanoTrack > AliAnalysisTaskC2Base::GetFMDhits() {
  // Relies on the event being vaild (no extra checks if object exists done here)
  AliAODForwardMult* aodForward =
    static_cast<AliAODForwardMult*>(fInputEvent->FindListObject("Forward"));
  // Shape of d2Ndetadphi: 200, -4, 6, 20, 0, 2pi
  const TH2D& d2Ndetadphi = aodForward->GetHistogram();
  Double_t fmdSignalThreshold = 0.01;
  Int_t nEta = d2Ndetadphi.GetXaxis()->GetNbins();
  Int_t nPhi = d2Ndetadphi.GetYaxis()->GetNbins();
  std::vector< AliAnalysisC2NanoTrack > ret_vector;
  for (Int_t iEta = 1; iEta <= nEta; iEta++) {
    Int_t valid = Int_t(d2Ndetadphi.GetBinContent(iEta, 0));
    if (!valid) {
       // No data expected for this eta
      continue;
    }
    Float_t eta = d2Ndetadphi.GetXaxis()->GetBinCenter(iEta);
    for (Int_t iPhi = 1; iPhi <= nPhi; iPhi++) {
      Float_t signal = d2Ndetadphi.GetBinContent(iEta, iPhi);
      if(signal > fmdSignalThreshold) {
	Float_t phi = d2Ndetadphi.GetYaxis()->GetBinCenter(iPhi);
	ret_vector.push_back(AliAnalysisC2NanoTrack(eta, phi, 0/*pt*/));
      }
    }
  }
  return ret_vector;
}

std::vector< AliAnalysisC2NanoTrack > &AliAnalysisTaskC2Base::GetValidTracks() {
  if (this->fCachedValues[cCachedValues::validTracks]) {
    return this->fValidTracks;
  }
  this->fValidTracks.clear();
  TIter nextTrack(this->GetAllTracks());
  while (TObject* obj = nextTrack()){
    // The naming around VTrack, VParticle, AODTrack mcParticle is a mess!
    // Take-away message: They all derive from AliVParticle one way or another.
    AliVParticle* particle = static_cast< AliVParticle* >(obj);
    if (!this->IsValidParticle(particle)){
      continue;
    }
    AliAnalysisC2NanoTrack tmp_track = {
      particle->Eta(),
      particle->Phi(),
      particle->Pt()
    };

    // The following is a bit ugly. We do not have the pt information
    // available in the reconstructed data. In order to compare apples
    // with apples in a MC colsure test, it is necessary to blind the
    // pt value in for MC truth in the FMD region. TODO: create a
    // switch for this in the settings; rewrite this nicer...
    if (this->fSettings.kMCTRUTH == this->fSettings.fDataType &&
	std::abs(tmp_track.eta) > 1.5) { // 1.5 is between its and fmd1/2
      tmp_track.pt = 0;
    }
    this->fValidTracks.push_back(tmp_track);
  }
  // Append the fmd hits to this vector if we are looking at reconstructed data,
  // All hits on the FMD (above the internally used threshold) are "valid"
  if (this->fSettings.kRECON == this->fSettings.fDataType) {
    auto fmdhits = this->GetFMDhits();
    this->fValidTracks.insert(this->fValidTracks.end(), fmdhits.begin(), fmdhits.end());
  }
  // fill the valid tracks into the appropriate QA histograms
  Float_t zvtx = this->InputEvent()->GetPrimaryVertex()->GetZ();
  for (auto &t: this->fValidTracks) {
    this->fetaVsZvtx->Fill(t.eta, zvtx);
  }
  // Mark the valid tracks as cached:
  this->fCachedValues[cCachedValues::validTracks] = true;
  return this->fValidTracks;
}

Float_t AliAnalysisTaskC2Base::GetEventClassifierValue() {
  if (this->fSettings.fMultEstimator == this->fSettings.fMultEstimatorValidTracks){
    // 7.049 seems to be the mean for 16j; this is an ugly hack!
    return this->GetValidTracks().size() / 7.049;
  }
  else {
    // Event is invalid if no multselection is present; ie. tested in IsValidEvent() already
    AliMultEstimator *multEstimator =
      (dynamic_cast< AliMultSelection* >(this->InputEvent()->FindListObject("MultSelection")))
      ->GetEstimator(this->fSettings.fMultEstimator);
    const Float_t multiplicity = multEstimator->GetValue() / multEstimator->GetMean();
    return multiplicity;
  }
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
