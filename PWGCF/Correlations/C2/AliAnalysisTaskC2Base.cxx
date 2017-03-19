#include <algorithm>
#include <iostream>
#include <random>
#include <vector>

#include "TH1F.h"
#include "THn.h"
#include "TMath.h"

#include "AliAODEvent.h"
#include "AliAODForwardMult.h"
#include "AliAODHandler.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisUtils.h"
#include "AliGenEventHeader.h"
#include "AliHeader.h"
#include "AliInputEventHandler.h"
#include "AliMultSelection.h"
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
    fmultDistribution(0),
    fEtaPhiZvtx_max_res(0)
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
    fmultDistribution(0),
    fEtaPhiZvtx_max_res(0)
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
  discardedEvtsAx->SetBinLabel(cDiscardEventReasons::SPDClusterVsTrackletBG + 1, "SPDClstrs vs trkl BG");
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
  discardedEvtsAx->SetBinLabel(cDiscardEventReasons::isMB + 1, "IsMB");
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
  discardedTracksAx->SetBinLabel(cDiscardTrackReasons::neutralCharge + 1, "neutral");
  discardedTracksAx->SetBinLabel(cDiscardTrackReasons::notHybridGCG + 1, "!HybridGCG");
  discardedTracksAx->SetBinLabel(cDiscardTrackReasons::notMCPrimary + 1, "not MC PhysicalPrimary");
  this->fOutputList->Add(this->fDiscardedTracks);

  this->fmultDistribution = new TH1F("multDistribution", "multDistribution", 210, 0, 210);
  this->fOutputList->Add(this->fmultDistribution);

  {
    const Int_t ndims = 3;
    const Int_t nbins[ndims] = {200, 20, 100};
    const Double_t lowerBounds[ndims] = {-4.0, 0, -10.0};
    const Double_t upperBounds[ndims] = {6.0, 2*TMath::Pi(), 10.0};
    this->fEtaPhiZvtx_max_res = new THnF("etaPhiZvtx_max_res",
					 "etaPhiZvtx_max_res;#eta;#phi;#zvtx;",
					 ndims,
					 nbins,
					 lowerBounds,
					 upperBounds
					 );
    this->fOutputList->Add(this->fEtaPhiZvtx_max_res);
  }
  AliLog::SetGlobalLogLevel(AliLog::kError);
  PostData(1, fOutputList);
}

Bool_t AliAnalysisTaskC2Base::IsValidEvent()
{
  // if (AliAnalysisC2Utils::EventFitsTrigger(AliVEvent::kINT7)) {
  //   this->fDiscardedEvents->Fill(cDiscardEventReasons::isMB);
  //   return false;
  // }
  //trigger check; TODO: ask Katarina about "fMBtrigger" in her code!
  if (TString firedTrigger = InputEvent()->GetFiredTriggerClasses()) {
    // Make choice base on type of MB trigger?
    if (!firedTrigger.Contains(this->fSettings.fTrigger_str)) {
      this->fDiscardedEvents->Fill(cDiscardEventReasons::noTrigger);
      return false;
    }
    // Fish the trigger mask of this event out of the input handler and compare it to the
    // desired one
    UInt_t fSelectMask= fInputHandler->IsEventSelected();
    if (fSelectMask & this->fSettings.fTriggerMask) {
      this->fDiscardedEvents->Fill(cDiscardEventReasons::noTrigger);
      return false;
    }
  }

  if (this->fSettings.fUseFMD) {
    AliAODForwardMult* aodForward =
      dynamic_cast<AliAODForwardMult*>(fInputEvent->FindListObject("Forward"));
    if (!aodForward) {
      this->fDiscardedEvents->Fill(cDiscardEventReasons::noForwardMultObj);
      return false;
    }
    if (aodForward->GetHistogram().GetIntegral() <= 0) {
      this->fDiscardedEvents->Fill(cDiscardEventReasons::noEntriesInFMD);
      return false;
    }
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
  // if (this->IsAsymmetricV0()) {
  //   this->fDiscardedEvents->Fill(cDiscardEventReasons::v0asymmetryCut);
  //   return false;
  // }
  //SPD vtx contributors (only aod events)
  if (AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(fInputEvent)) {
    if (aodEvent->GetPrimaryVertexSPD()->GetNContributors() < 1) {
      this->fDiscardedEvents->Fill(cDiscardEventReasons::spdVertexContributors);
      return false;
    }
  }
  // Aliphysics Analysis utils as it it used by several checks below
  AliAnalysisUtils utils;
  //out-of-bunch 11 BC
  if (utils.IsOutOfBunchPileUp(this->InputEvent())) {
    this->fDiscardedEvents->Fill(cDiscardEventReasons::isOutOfBunchPileup);
    return false;
  }
  if (utils.IsSPDClusterVsTrackletBG(this->InputEvent())) {
    this->fDiscardedEvents->Fill(cDiscardEventReasons::SPDClusterVsTrackletBG);
    return false;
  }
  // SPD pileup
  if (utils.IsPileUpSPD(InputEvent())) {
    this->fDiscardedEvents->Fill(cDiscardEventReasons::spdPipeup);
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
}

TClonesArray* AliAnalysisTaskC2Base::GetAllTracks() {
  // If we are dealing with an ESD event, we have to have an AOD handler as well!
  // We get all the particles/tracks from this AOD handler. This function makes no asumptions
  // If a track is valid (ie. caused a hit on the FMD) or not. We have to check that later!
  AliAODEvent* aodEvent = dynamic_cast< AliAODEvent* >(this->InputEvent());
  // Yes, the following is realy "aodEvent" not mcEvent :P
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
  Int_t nEta = d2Ndetadphi.GetXaxis()->GetNbins();
  Int_t nPhi = d2Ndetadphi.GetYaxis()->GetNbins();
  std::vector< AliAnalysisC2NanoTrack > ret_vector;
  // minimum pt value for this analysis. This pt value is used for the FMD tracks
  // The reason being that otherwise a pt bin with pt<0 has to exist which is a wast of memory
  // if we want to also use ITS tracks with real pt values
  Double_t analysis_pt_min = this->fSettings.fPtBinEdges.at(0);
  for (Int_t iEta = 1; iEta <= nEta; iEta++) {
    Int_t valid = Int_t(d2Ndetadphi.GetBinContent(iEta, 0));
    if (!valid) {
       // No data expected for this eta
      continue;
    }
    Float_t eta = d2Ndetadphi.GetXaxis()->GetBinCenter(iEta);
    for (Int_t iPhi = 1; iPhi <= nPhi; iPhi++) {
      // Bin content is most likely number of particles!
      Float_t mostProbableN = d2Ndetadphi.GetBinContent(iEta, iPhi);
      if (mostProbableN > 0) {
	Float_t phi = d2Ndetadphi.GetYaxis()->GetBinCenter(iPhi);
	ret_vector.push_back(AliAnalysisC2NanoTrack(eta, phi, analysis_pt_min, mostProbableN));
      }
    }
  }
  // Above, we looped through the histogram in an orderly
  // fashion. This means that the eta and phi values are in increasing
  // orders. For eta this is not an issue. We will just always have
  // the trigger (or associated) with a larger eta. In a eta1, eta2
  // plots, this will show by all values being in one triangle. That
  // is ok. However, for phi this is not ok if we require
  // eta1==eta2. In that case, phi2 > phi1. This leads to a very
  // shifted correlation function. Especially if the eta bin width is
  // larger than the resolution. Bottom line, I suffle this vector
  // here to avoid this mess.
  std::random_device rd;
  std::default_random_engine engine{rd()};
  std::shuffle(std::begin(ret_vector), std::end(ret_vector), engine);
  return ret_vector;
}

std::vector< AliAnalysisC2NanoTrack > AliAnalysisTaskC2Base::GetSPDhits() {
  AliAODEvent* aodEvent = dynamic_cast< AliAODEvent* >(this->InputEvent());
  AliAODTracklets* aodTracklets = aodEvent->GetTracklets();
  Double_t analysis_pt_min = this->fSettings.fPtBinEdges.at(0);
  std::vector< AliAnalysisC2NanoTrack > ret_vector;
  for (Int_t i = 0; i < aodTracklets->GetNumberOfTracklets(); i++) {
    auto eta = aodTracklets->GetEta(i);
    auto phi = aodTracklets->GetPhi(i);
    ret_vector.push_back(AliAnalysisC2NanoTrack(eta, phi, analysis_pt_min, 1));
  }
  return ret_vector;
}

std::vector< AliAnalysisC2NanoTrack > AliAnalysisTaskC2Base::GetValidCentralTracks() {
  std::vector< AliAnalysisC2NanoTrack > ret_vector;
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
      particle->Pt(),
      1 /* track weight; tracks in FMD have a probability assigned */ 
    };
    ret_vector.push_back(tmp_track);
  }
  return ret_vector;
}


std::vector< AliAnalysisC2NanoTrack > &AliAnalysisTaskC2Base::GetValidTracks() {
  if (this->fCachedValues[cCachedValues::validTracks]) {
    return this->fValidTracks;
  }
  this->fValidTracks.clear();

  // Append central tracks
  // auto centralTracks = this->GetValidCentralTracks();
  // this->fValidTracks.insert(this->fValidTracks.end(), centralTracks.begin(), centralTracks.end());

  // Append central tracklets
  if (this->fSettings.kRECON == this->fSettings.fDataType) {
    auto spdhits = this->GetSPDhits();
    this->fValidTracks.insert(this->fValidTracks.end(), spdhits.begin(), spdhits.end());
  
    // Append the fmd hits to this vector if we are looking at reconstructed data,
    // All hits on the FMD (above the internally used threshold) are "valid"
    if (this->fSettings.fUseFMD) {
      auto fmdhits = this->GetFMDhits();
      this->fValidTracks.insert(this->fValidTracks.end(), fmdhits.begin(), fmdhits.end());
    }
  } else {
    // Sorry future me: GetValidCentralTracks returns all tracks in case of MC truth :P
    auto allMCtracks = this->GetValidCentralTracks();
    this->fValidTracks.insert(this->fValidTracks.end(), allMCtracks.begin(), allMCtracks.end());
  }
  // fill the valid tracks into the appropriate QA histograms
  Float_t zvtx = this->InputEvent()->GetPrimaryVertex()->GetZ();
  for (auto &t: this->fValidTracks) {
    const Double_t stuffing[] = {t.eta, t.phi, zvtx};
    this->fEtaPhiZvtx_max_res->Fill(stuffing, t.weight);
  }
  // Mark the valid tracks as cached:
  this->fCachedValues[cCachedValues::validTracks] = true;
  return this->fValidTracks;
}

Float_t AliAnalysisTaskC2Base::GetEventClassifierValue() {
  if (this->fSettings.fMultEstimator == this->fSettings.fMultEstimatorValidTracks){
    return this->GetValidTracks().size();
  }
  else {
    // Event is invalid if no multselection is present; ie. tested in IsValidEvent() already
    AliMultEstimator *multEstimator =
      (dynamic_cast< AliMultSelection* >(this->InputEvent()->FindListObject("MultSelection")))
      ->GetEstimator(this->fSettings.fMultEstimator);
    const Float_t multiplicity = ((Float_t)multEstimator->GetValue()) / multEstimator->GetMean();
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

Bool_t AliAnalysisTaskC2Base::IsAODdataset() {
  AliAODHandler *aodHandler =
    dynamic_cast< AliAODHandler* >(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (aodHandler) {
    return true;
  }
  else {
    return false;
  }
}
