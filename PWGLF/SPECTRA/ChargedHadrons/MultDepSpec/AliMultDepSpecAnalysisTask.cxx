#include "AliMultDepSpecAnalysisTask.h"

/// \cond CLASSIMP
ClassImp(AliMultDepSpecAnalysisTask);
/// \endcond

/***************************************************************************//**
 * ROOT I/O Constructor.
 ******************************************************************************/
 AliMultDepSpecAnalysisTask::AliMultDepSpecAnalysisTask() : AliAnalysisTaskSE(),
   //General member variables
   fOutputList(nullptr),
   fEventCuts(),
   fESDtrackCuts(nullptr),
   fRand(nullptr),
   fMCSpectraWeights(nullptr),
   fCutMode(100),
   //Toggles
   fIsESD(kTRUE),
   fIsMC(kFALSE),
   fUseCent(kFALSE),
   fUseZDCCut(kFALSE),
   fOverridePbPbEventCuts(kFALSE),
   fMCUseDataDrivenCorrections(kFALSE),
   fMCSecScalingSysFlag(0),
   // Cut Parameters
   fTriggerMask(AliVEvent::kMB | AliVEvent::kINT7),
   fMinEta(-10),
   fMaxEta(10),
   fMinPt(0.0),
   fMaxPt(50.0),
   fMaxZv(10.0),
   fMinCent(-1000.0),
   fMaxCent(1000.0),
   //Arrays for Binning
   fBinsEventCuts(nullptr),
   fBinsMult(nullptr),
   fBinsCent(nullptr),
   fBinsPt(nullptr),
   fBinsEta(nullptr),
   fBinsZv(nullptr),
   fBinsPtReso(nullptr),
   //Event-Histograms
   fHistEventSelection(nullptr),
   fHistEvents(nullptr),
   fHistTracks(nullptr),
   fHistRelPtReso(nullptr),
   fHistMCEventEfficiency(nullptr),
   fHistMCEventEfficiencyScaled(nullptr),
   fHistMCRelPtReso(nullptr),
   fHistMCMultCorrelMatrix(nullptr),
   fHistMCPtCorrelMatrix(nullptr),
   fHistMCEtaCorrelMatrix(nullptr),
   fHistMCPrimTrue(nullptr),
   fHistMCPrimMeas(nullptr),
   fHistMCSecMeas(nullptr),
   fHistMCEdgeContam(nullptr),
   fHistMCDoubleCountig(nullptr),
   fHistMCEventsScaled(nullptr),
   fHistMCTracksScaled(nullptr),
   fHistMCMultCorrelMatrixScaled(nullptr),
   fHistMCPrimTrueScaled(nullptr),
   fHistMCPrimMeasScaled(nullptr),
   fHistMCSecMeasScaled(nullptr),
   fHistMCEdgeContamScaled(nullptr),
   fHistMCMultMeasScaleEffect(nullptr),
   fHistMCMultTrueScaleEffect(nullptr),
   // transient event and track properties
   fEvent(nullptr),
   fMCEvent(nullptr),
   fCent(0),
   fMultMeas(0),
   fMultTrue(0),
   fMultMeasScaled(0),
   fMultTrueScaled(0),
   fRunNumber(0),
   fEventNumberInFile(0),
   fTimeStamp(0),
   fPt(0),
   fEta(0),
   fSigmaPt(0),
   fMCPt(0),
   fMCEta(0),
   fMCLabel(0),
   fIsParticleInAcceptance(kFALSE),
   fMCIsPhysicalPrimary(kFALSE),
   fMCIsCharged(kFALSE),
   fMCIsChargedPrimary(kFALSE),
   fMCIsChargedSecondary(kFALSE),
   fMCParticleWeight(1.0),
   fMCSecScaleWeight(1.0),
   fNRepetitions(1),
   fUseRandomSeed(kFALSE)
{
  // ROOT IO constructor, don't allocate memory here!
}


/***************************************************************************//**
 * Constructor.
 ******************************************************************************/
AliMultDepSpecAnalysisTask::AliMultDepSpecAnalysisTask(const char* name) : AliAnalysisTaskSE(name),
  //General member variables
  fOutputList(nullptr),
  fEventCuts(),
  fESDtrackCuts(nullptr),
  fRand(nullptr),
  fMCSpectraWeights(nullptr),
  fCutMode(100),
  //Toggles
  fIsESD(kTRUE),
  fIsMC(kFALSE),
  fUseCent(kFALSE),
  fUseZDCCut(kFALSE),
  fOverridePbPbEventCuts(kFALSE),
  fMCUseDataDrivenCorrections(kFALSE),
  fMCSecScalingSysFlag(0),
  // Cut Parameters
  fTriggerMask(AliVEvent::kMB | AliVEvent::kINT7),
  fMinEta(-10),
  fMaxEta(10),
  fMinPt(0.0),
  fMaxPt(50.0),
  fMaxZv(10.0),
  fMinCent(-1000.0),
  fMaxCent(1000.0),
  //Arrays for Binning
  fBinsEventCuts(nullptr),
  fBinsMult(nullptr),
  fBinsCent(nullptr),
  fBinsPt(nullptr),
  fBinsEta(nullptr),
  fBinsZv(nullptr),
  fBinsPtReso(nullptr),
  //Event-Histograms
  fHistEventSelection(nullptr),
  fHistEvents(nullptr),
  fHistTracks(nullptr),
  fHistRelPtReso(nullptr),
  fHistMCEventEfficiency(nullptr),
  fHistMCEventEfficiencyScaled(nullptr),
  fHistMCRelPtReso(nullptr),
  fHistMCMultCorrelMatrix(nullptr),
  fHistMCPtCorrelMatrix(nullptr),
  fHistMCEtaCorrelMatrix(nullptr),
  fHistMCPrimTrue(nullptr),
  fHistMCPrimMeas(nullptr),
  fHistMCSecMeas(nullptr),
  fHistMCEdgeContam(nullptr),
  fHistMCDoubleCountig(nullptr),
  fHistMCEventsScaled(nullptr),
  fHistMCTracksScaled(nullptr),
  fHistMCMultCorrelMatrixScaled(nullptr),
  fHistMCPrimTrueScaled(nullptr),
  fHistMCPrimMeasScaled(nullptr),
  fHistMCSecMeasScaled(nullptr),
  fHistMCEdgeContamScaled(nullptr),
  fHistMCMultMeasScaleEffect(nullptr),
  fHistMCMultTrueScaleEffect(nullptr),
  // transient event and track properties
  fEvent(nullptr),
  fMCEvent(nullptr),
  fCent(0),
  fMultMeas(0),
  fMultTrue(0),
  fMultMeasScaled(0),
  fMultTrueScaled(0),
  fRunNumber(0),
  fEventNumberInFile(0),
  fTimeStamp(0),
  fPt(0),
  fEta(0),
  fSigmaPt(0),
  fMCPt(0),
  fMCEta(0),
  fMCLabel(0),
  fIsParticleInAcceptance(kFALSE),
  fMCIsPhysicalPrimary(kFALSE),
  fMCIsCharged(kFALSE),
  fMCIsChargedPrimary(kFALSE),
  fMCIsChargedSecondary(kFALSE),
  fMCParticleWeight(1.0),
  fMCSecScaleWeight(1.0),
  fNRepetitions(1),
  fUseRandomSeed(kFALSE)
{
  // Set default binning
  Double_t binsEventCutsDefault[8] = {-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.6};
  fBinsEventCuts = new TArrayD(8, binsEventCutsDefault);

  Double_t binsMultDefault[2] = {0., 10000.};
  Double_t binsCentDefault[2] = {0., 100.};
  Double_t binsPtDefault[53] = {0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.5,5.0,5.5,6.0,6.5,7.0,8.0,9.0,10.0,20.0,30.0,40.0,50.0,60.0};
  Double_t binsEtaDefault[19] = {-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9};
  Double_t binsZvDefault[13] = {-30.,-25.,-20.,-15.,-10.,-5.,0.,5.,10.,15.,20.,25.,30.};

  // binning for relative pT resolution
  const Int_t nBinsPtReso = 300;
  Double_t binsPtReso[nBinsPtReso+1];
  SetFixedBinEdges(binsPtReso, 0., 0.3, nBinsPtReso);
  SetBinsPtReso(nBinsPtReso, binsPtReso);

  SetBinsMult(1, binsMultDefault);
  SetBinsCent(1, binsCentDefault);
  SetBinsPt(52, binsPtDefault);
  SetBinsEta(18, binsEtaDefault);
  SetBinsZv(12, binsZvDefault);

  DefineOutput(1, TList::Class());
}

/***************************************************************************//**
 * Function executed once before the event loop. Create histograms here.
 ******************************************************************************/
void AliMultDepSpecAnalysisTask::UserCreateOutputObjects(){

  OpenFile(1, "recreate");
  fOutputList = new TList();
  fOutputList->SetOwner();

  //TList* eventSelection = new TList();
  //eventSelection->SetName("eventSelection");
  //fEventCuts.AddQAplotsToList(eventSelection);
  //fOutputList->Add(eventSelection);

  fHistEventSelection = CreateHistogram("fHistEventSelection", {"eventcuts"});
  fOutputList->Add(fHistEventSelection);

  fHistEvents = CreateHistogram("fHistEvents", {"mult_meas", "cent"});
  fOutputList->Add(fHistEvents);
  fHistTracks = CreateHistogram("fHistTracks", {"pt_meas", "eta_meas", "mult_meas", "cent"});
  fOutputList->Add(fHistTracks);

  fHistRelPtReso = CreateHistogram("fHistRelPtReso", {"sigmapt", "pt_meas", "cent"});
  fOutputList->Add(fHistRelPtReso);

  if(fIsMC)
  {
    fHistMCEventEfficiency = CreateHistogram("fHistMCEventEfficiency", {"eventcuts", "mult_true"});
    fOutputList->Add(fHistMCEventEfficiency);

    fHistMCEventEfficiencyScaled = CreateHistogram("fHistMCEventEfficiencyScaled", {"eventcuts", "mult_true"});
    fOutputList->Add(fHistMCEventEfficiencyScaled);

    fHistMCRelPtReso = CreateHistogram("fHistMCRelPtReso", {"deltapt", "pt_meas", "cent"});
    fOutputList->Add(fHistMCRelPtReso);
    fHistMCMultCorrelMatrix = CreateHistogram("fHistMCMultCorrelMatrix", {"mult_meas", "mult_true"});
    fOutputList->Add(fHistMCMultCorrelMatrix);
    fHistMCPtCorrelMatrix = CreateHistogram("fHistMCPtCorrelMatrix", {"pt_meas", "pt_true"});
    fOutputList->Add(fHistMCPtCorrelMatrix);
    fHistMCEtaCorrelMatrix = CreateHistogram("fHistMCEtaCorrelMatrix", {"eta_meas", "eta_true"});
    fOutputList->Add(fHistMCEtaCorrelMatrix);
    fHistMCPrimTrue = CreateHistogram("fHistMCPrimTrue", {"pt_true", "eta_true", "mult_true", "cent"});
    fOutputList->Add(fHistMCPrimTrue);
    fHistMCPrimMeas = CreateHistogram("fHistMCPrimMeas", {"pt_true", "eta_true", "mult_true", "cent"});
    fOutputList->Add(fHistMCPrimMeas);
    fHistMCSecMeas = CreateHistogram("fHistMCSecMeas", {"pt_true", "eta_true", "mult_true", "cent"});
    fOutputList->Add(fHistMCSecMeas);
    fHistMCEdgeContam = CreateHistogram("fHistMCEdgeContam", {"pt_meas", "eta_meas", "mult_true", "cent"});
    fOutputList->Add(fHistMCEdgeContam);
    fHistMCDoubleCountig = CreateLogHistogram("fHistMCDoubleCountig");
    fOutputList->Add(fHistMCDoubleCountig);

    // scaled histos
    fHistMCEventsScaled = CreateHistogram("fHistMCEventsScaled", {"mult_meas", "cent"});
    fOutputList->Add(fHistMCEventsScaled);
    fHistMCTracksScaled = CreateHistogram("fHistMCTracksScaled", {"pt_meas", "eta_meas", "mult_meas", "cent"});
    fOutputList->Add(fHistMCTracksScaled);
    fHistMCMultCorrelMatrixScaled = CreateHistogram("fHistMCMultCorrelMatrixScaled", {"mult_meas", "mult_true"});
    fOutputList->Add(fHistMCMultCorrelMatrixScaled);
    fHistMCPrimTrueScaled = CreateHistogram("fHistMCPrimTrueScaled", {"pt_true", "eta_true", "mult_true", "cent"});
    fOutputList->Add(fHistMCPrimTrueScaled);
    fHistMCPrimMeasScaled = CreateHistogram("fHistMCPrimMeasScaled", {"pt_true", "eta_true", "mult_true", "cent"});
    fOutputList->Add(fHistMCPrimMeasScaled);
    fHistMCSecMeasScaled = CreateHistogram("fHistMCSecMeasScaled", {"pt_true", "eta_true", "mult_true", "cent"});
    fOutputList->Add(fHistMCSecMeasScaled);
    fHistMCEdgeContamScaled = CreateHistogram("fHistMCEdgeContamScaled", {"pt_meas", "eta_meas", "mult_true", "cent"});
    fOutputList->Add(fHistMCEdgeContamScaled);

    fHistMCMultMeasScaleEffect = CreateHistogram("fHistMCMultMeasScaleEffect", {"mult_true", "mult_true"});
    fOutputList->Add(fHistMCMultMeasScaleEffect);
    fHistMCMultTrueScaleEffect = CreateHistogram("fHistMCMultTrueScaleEffect", {"mult_meas", "mult_meas"});
    fOutputList->Add(fHistMCMultTrueScaleEffect);
  }

  // override event automatic event selection settings
  if(fOverridePbPbEventCuts)
  {
    fEventCuts.SetManualMode();
    fEventCuts.SetupLHC15o(); // first set default values and then override
    fEventCuts.SetCentralityRange(fMinCent, fMaxCent);
    fEventCuts.fUseEstimatorsCorrelationCut = false;
  }
  //fEventCuts.SetMaxVertexZposition(fMaxZv); // has no effect in automatic mode...
  fEventCuts.OverrideAutomaticTriggerSelection(fTriggerMask);

  if(fIsESD) InitESDTrackCuts();

  fRand = new TRandom3();

  PostData(1, fOutputList);
}

/// Destructor
AliMultDepSpecAnalysisTask::~AliMultDepSpecAnalysisTask(){
  if(fESDtrackCuts){delete fESDtrackCuts; fESDtrackCuts = nullptr;}
  if(fRand){delete fRand; fRand = nullptr;}
}

/***************************************************************************//**
 * Function executed for each event.
 ******************************************************************************/
void AliMultDepSpecAnalysisTask::UserExec(Option_t *){
  if(!InitEvent()) return;
  FillEventHistos();
  LoopMeas();
  if(fIsMC) LoopTrue();
  PostData(1, fOutputList);
}

/***************************************************************************//**
 * Function executed after all events were processed.
 ******************************************************************************/
void AliMultDepSpecAnalysisTask::Terminate(Option_t*)
{

}

/***************************************************************************//**
 * Function to get data-driven secondary scaling weights.
 ******************************************************************************/
Double_t AliMultDepSpecAnalysisTask::GetSecScalingFactor(AliMCParticle* particle)
{
  return AlidNdPtTools::MCScalingFactor(particle, fMCEvent, fMCSecScalingSysFlag); // syst variations -1, +1
}

/***************************************************************************//**
 * Function to get data-driven particle composition weights.
 ******************************************************************************/
Double_t AliMultDepSpecAnalysisTask::GetParticleWeight(AliMCParticle* particle)
{
  if(!fMCSpectraWeights) {return 1.0;}
  else return fMCSpectraWeights->GetMCSpectraWeight(particle->Particle(), fMCEvent);
}

/***************************************************************************//**
 * Initialize event quantities. Sets fEvent, fMCEvent, fMeasMult, fTrueMult
 ******************************************************************************/
Bool_t AliMultDepSpecAnalysisTask::InitEvent()
{
    fEvent = InputEvent();
    if (!fEvent) {AliError("fEvent not available\n"); return kFALSE;}

    if(fIsMC){
      fMCEvent = MCEvent();
      if (!fMCEvent) {AliError("fMCEvent not available\n"); return kFALSE;}
    }

    fCent = ((fBinsCent->GetSize()-1) > 1) ? GetCentrality(fEvent) : 50;
    // event info for random seeds
    fRunNumber = fEvent->GetRunNumber();
    AliESDEvent* esdEvent = dynamic_cast<AliESDEvent*>(fEvent);
    fEventNumberInFile = esdEvent->GetEventNumberInFile();
    fTimeStamp = esdEvent->GetTimeStamp();
    //esdEvent->GetHeader()->GetEventIdAsLong();

    if(fUseZDCCut)
    {
      AliESDZDC* zdc = esdEvent->GetESDZDC();
      Double_t znaEnergy = zdc->GetZNAEnergy();
      Double_t zncEnergy = zdc->GetZNCEnergy();
      if ( (znaEnergy <  3503.0) || (zncEnergy <  3609.0) ) {return kFALSE;}
    }
    Bool_t acceptEvent = fEventCuts.AcceptEvent(fEvent);

    LoopMeas(kTRUE); // set measured multiplicity fMeasMult, fMeasMultScaled
    if(fIsMC) LoopTrue(kTRUE); // set true multiplicity fTrueMult, fTrueMultScaled

    // fill histograms for event selection and Nch dependent efficiency
    std::array <AliEventCuts::NormMask, 5> norm_masks {
      AliEventCuts::kAnyEvent,
      AliEventCuts::kTriggeredEvent,
      AliEventCuts::kPassesNonVertexRelatedSelections,
      AliEventCuts::kHasReconstructedVertex,
      AliEventCuts::kPassesAllCuts // => acceptEvent = kTRUE
    };
    for (int iC = 0; iC < 5; ++iC) {
      if (fEventCuts.CheckNormalisationMask(norm_masks[iC])) {
        FillHisto(fHistEventSelection, {Double_t(iC)});
        if(fIsMC){
          FillHisto(fHistMCEventEfficiency, {Double_t(iC), fMultTrue});
          FillHisto(fHistMCEventEfficiencyScaled, {Double_t(iC), fMultTrueScaled});
        }
      }
    }
    // now count only the events which also contribute to the measurement
    if(acceptEvent){
      if (fMultMeas > 0) FillHisto(fHistEventSelection, {5.0});
      if(fIsMC){
        if (fMultMeas > 0) FillHisto(fHistMCEventEfficiency, {5.0, fMultTrue});
        if (fMultMeasScaled > 0) FillHisto(fHistMCEventEfficiencyScaled, {5.0, fMultTrueScaled});
      }
    }
    // additional info: triggered and vertex in acceptacne
    if(fEventCuts.CheckNormalisationMask(AliEventCuts::kTriggeredEvent) && fEventCuts.PassedCut(AliEventCuts::kVertexPosition))
    {
      FillHisto(fHistEventSelection, {6.0});
      if(fIsMC){
        FillHisto(fHistMCEventEfficiency, {6.0, fMultTrue});
        FillHisto(fHistMCEventEfficiencyScaled, {6.0, fMultTrueScaled});
      }
    }

    return acceptEvent;
}

/***************************************************************************//**
 * Fill event histograms. Event related members are set.
 ******************************************************************************/
void AliMultDepSpecAnalysisTask::FillEventHistos()
{
  FillHisto(fHistEvents, {fMultMeas, fCent});
  if(fIsMC){
    FillHisto(fHistMCMultCorrelMatrix, {fMultMeas, fMultTrue});
    // scaled histos
    FillHisto(fHistMCEventsScaled, {fMultMeasScaled, fCent});
    FillHisto(fHistMCMultCorrelMatrixScaled, {fMultMeasScaled, fMultTrueScaled});

    FillHisto(fHistMCMultMeasScaleEffect, {fMultMeas, fMultMeasScaled});
    FillHisto(fHistMCMultTrueScaleEffect, {fMultTrue, fMultTrueScaled});
  }
}

/***************************************************************************//**
 * Fill track histograms. Track related members are set.
 ******************************************************************************/
void AliMultDepSpecAnalysisTask::FillMeasTrackHistos()
{
  FillHisto(fHistTracks, {fPt, fEta, fMultMeas, fCent});
  FillHisto(fHistRelPtReso, {fSigmaPt, fPt, fCent});
}

/***************************************************************************//**
 * Fill scaled track histograms. Track related members are set.
 ******************************************************************************/
void AliMultDepSpecAnalysisTask::FillMeasScaledTrackHistos()
{
  FillHisto(fHistMCTracksScaled, {fPt, fEta, fMultMeasScaled, fCent});
}

/***************************************************************************//**
 * Fill measured particle histograms. Track and mc particle related members are set.
 ******************************************************************************/
void AliMultDepSpecAnalysisTask::FillMeasParticleHistos()
{
  if(!fIsParticleInAcceptance)
  {
    FillHisto(fHistMCEdgeContam, {fPt, fEta, fMultTrue, fCent});
  }
  else
  {
    FillHisto(fHistMCRelPtReso, {TMath::Abs(fPt - fMCPt)/fPt, fPt, fCent});

    if(fMCIsChargedPrimary)
    {
      FillHisto(fHistMCPtCorrelMatrix, {fPt, fMCPt});
      FillHisto(fHistMCEtaCorrelMatrix, {fEta, fMCEta});
      FillHisto(fHistMCPrimMeas, {fMCPt, fMCEta, fMultTrue, fCent});
    }else{
      FillHisto(fHistMCSecMeas, {fMCPt, fMCEta, fMultTrue, fCent});
    }
  }
}

/***************************************************************************//**
 * Fill measured scaled particle histograms. Track and mc particle related members are set.
 ******************************************************************************/
void AliMultDepSpecAnalysisTask::FillMeasScaledParticleHistos()
{
  if(!fIsParticleInAcceptance)
  {
    FillHisto(fHistMCEdgeContamScaled, {fPt, fEta, fMultTrueScaled, fCent});
  }
  else
  {
    if(fMCIsChargedPrimary)
    {
      FillHisto(fHistMCPrimMeasScaled, {fMCPt, fMCEta, fMultTrueScaled, fCent});
    }else{
      FillHisto(fHistMCSecMeasScaled, {fMCPt, fMCEta, fMultTrueScaled, fCent});
    }
  }
}

/***************************************************************************//**
 * Fill generated particle histograms. MC particle related members are set.
 ******************************************************************************/
void AliMultDepSpecAnalysisTask::FillTrueParticleHistos()
{
  if(fMCIsChargedPrimary && fIsParticleInAcceptance) FillHisto(fHistMCPrimTrue, {fMCPt, fMCEta, fMultTrue, fCent});
}

/***************************************************************************//**
 * Fill scaled generated particle histograms. MC particle related members are set.
 ******************************************************************************/
void AliMultDepSpecAnalysisTask::FillTrueScaledParticleHistos()
{
  if(fMCIsChargedPrimary && fIsParticleInAcceptance) FillHisto(fHistMCPrimTrueScaled, {fMCPt, fMCEta, fMultTrueScaled, fCent});
}


/***************************************************************************//**
 * Loop over measured tracks. Can either count multiplicity or fill histograms.
 ******************************************************************************/
void AliMultDepSpecAnalysisTask::LoopMeas(Bool_t count)
{
  if(count) {fMultMeas = 0; fMultMeasScaled = 0;}
  AliVTrack* track = nullptr;
  AliMCParticle* particle = nullptr;
  std::vector<Int_t> mcLableLedger;
  for (Int_t i = 0; i < fEvent->GetNumberOfTracks(); i++){
    track = fEvent->GetVTrack(i);
    // Set fPt, fEta, fSigmapt; Check if track in kin range and has good quality
    if(!InitTrack(track)) continue;

    // initialize particle properties
    if(fIsMC)
    {
      particle  = (AliMCParticle*)fMCEvent->GetTrack(TMath::Abs(track->GetLabel()));
      // Set fMCPt, fMCEta, fMCIsChargedPrimary; Check if particle in kin range
      if(!InitParticle(particle)) continue;
      // Control hist to check if one particle results in multiple tracks
      if(!count)
      {
        if(std::find(mcLableLedger.begin(), mcLableLedger.end(), fMCLabel) != mcLableLedger.end()){
            FillLogHisto(fHistMCDoubleCountig, "doubleCountedTracks");
        } else {
            mcLableLedger.push_back(fMCLabel);
        }
      }
    }

    if(count)
    {
      fMultMeas++;
      if(fIsMC) fMultMeasScaled += fNRepetitions;
    }
    else
    {
      FillMeasTrackHistos();
      if(fIsMC)
      {
        FillMeasParticleHistos();
        // mc scaled
        for(Int_t i = 0; i < fNRepetitions; i++)
        {
          FillMeasScaledTrackHistos();
          FillMeasScaledParticleHistos();
        }
      }
    }
  }
}

/***************************************************************************//**
 * Loop over generated mc particles. Can either count multiplicity or fill histograms.
 ******************************************************************************/
void AliMultDepSpecAnalysisTask::LoopTrue(Bool_t count)
{
  if(count) {fMultTrue = 0; fMultTrueScaled = 0;}
  AliMCParticle* particle = nullptr;
  for(Int_t i = 0; i < fMCEvent->GetNumberOfTracks(); i++) {
    particle  = (AliMCParticle*)fMCEvent->GetTrack(i);
    // Sets fMCPt, fMCEta, ... and checks if particle in kin range
    if(!InitParticle(particle)) continue;

    // mc truth
    if(count)
    {
      if(fMCIsChargedPrimary && fIsParticleInAcceptance){
         fMultTrue++;
         fMultTrueScaled += fNRepetitions;
      }
    }else{
      FillTrueParticleHistos();
      // mc scaled
      for(Int_t i = 0; i < fNRepetitions; i++)
      {
        FillTrueScaledParticleHistos();
      }
    }
  }
}



/***************************************************************************//**
 * Initializes track properties and returns false if track bad or out of range.
 ******************************************************************************/
Bool_t AliMultDepSpecAnalysisTask::InitTrack(AliVTrack* track)
{
  if(!track) {AliError("Track not available\n"); return kFALSE;}
  fPt = track->Pt();
  fEta = track->Eta();
  fSigmaPt = 1./TMath::Abs(dynamic_cast<AliESDtrack*>(track)->GetSigned1Pt())*TMath::Sqrt(dynamic_cast<AliESDtrack*>(track)->GetSigma1Pt2());


  if(fPt  <= fMinPt  + PRECISION)   return kFALSE;
  if(fPt  >= fMaxPt  - PRECISION)   return kFALSE;
  if(fEta <= fMinEta + PRECISION)   return kFALSE;
  if(fEta >= fMaxEta - PRECISION)   return kFALSE;
  if(!AcceptTrackQuality(track))    return kFALSE;

  return kTRUE;
}

/***************************************************************************//**
 * Initializes particle properties and returns false if out of range.
 ******************************************************************************/
Bool_t AliMultDepSpecAnalysisTask::InitParticle(AliMCParticle* particle)
{
  if(!particle) {AliError("Particle not available\n"); return kFALSE;}
  fMCPt = particle->Pt();
  fMCEta = particle->Eta();

  fIsParticleInAcceptance = kTRUE;

  if((fMCPt  <= fMinPt  + PRECISION)  || (fMCPt  >= fMaxPt  - PRECISION) ||
  (fMCEta <= fMinEta + PRECISION) || (fMCEta >= fMaxEta - PRECISION))
    fIsParticleInAcceptance = kFALSE;

  fMCLabel = particle->GetLabel();
  fMCIsPhysicalPrimary = (fMCEvent->IsPhysicalPrimary(fMCLabel)) ? kTRUE : kFALSE;
  fMCIsCharged = ((TMath::Abs(particle->Charge()) > 0.01)) ? kTRUE : kFALSE;
  fMCIsChargedPrimary = fMCIsCharged && fMCIsPhysicalPrimary;
  fMCIsChargedSecondary = fMCIsCharged && !fMCIsPhysicalPrimary;

  fMCParticleWeight = 1.0;
  fMCSecScaleWeight = 1.0;
  fNRepetitions = 1;

  if(fMCUseDataDrivenCorrections)
  {
    if(fMCIsChargedPrimary)
    {
      fMCParticleWeight = GetParticleWeight(particle);
      fNRepetitions = GetNRepetitons(fMCParticleWeight);
    }
    else if(fMCIsChargedSecondary)
    {
      fMCSecScaleWeight = GetSecScalingFactor(particle);
      fNRepetitions = GetNRepetitons(fMCSecScaleWeight);
    }
  }
  return kTRUE;
}

/***************************************************************************//**
 * Decide how often to repeat particle in MC to match data.
 ******************************************************************************/
Int_t AliMultDepSpecAnalysisTask::GetNRepetitons(Double_t scalingFactor)
{
  Int_t nRepetitions = (Int_t)scalingFactor;
  Double_t rest = scalingFactor - nRepetitions;

  fRand->SetSeed(GetSeed());
  nRepetitions += (fRand->Rndm() <= rest) ? 1 : 0;

  return nRepetitions;
}

/***************************************************************************//**
 * Define random (but reproducable) seed.
 ******************************************************************************/
UInt_t AliMultDepSpecAnalysisTask::GetSeed()
{
    if (fUseRandomSeed) { return 0; }

    UInt_t seed = fEventNumberInFile;
    seed <<= 7;
    seed += fRunNumber;
    seed <<= 7;
    seed += fMCLabel;
    seed <<= 7;
    seed += fTimeStamp;
    return seed;
}

/***************************************************************************//**
 * Function to select tracks with required quality.
 ******************************************************************************/
Bool_t AliMultDepSpecAnalysisTask::AcceptTrackQuality(AliVTrack* track){
  if(fIsESD){
    AliESDtrack* esdTrack = dynamic_cast<AliESDtrack*> (track);
    if(!fESDtrackCuts->AcceptTrack(esdTrack)) return kFALSE;
  }
  return kTRUE;
}

/***************************************************************************//**
 * Function to obtain V0M centrality.
 ******************************************************************************/
Double_t AliMultDepSpecAnalysisTask::GetCentrality(AliVEvent* event)
{
  AliMultSelection* multSelection = (AliMultSelection*) fEvent->FindListObject("MultSelection");
  if(!multSelection){AliError("No MultSelection found!"); return 999;}
  return multSelection->GetMultiplicityPercentile("V0M");
}

/***************************************************************************//**
 * Function to initialize the ESD track cuts object.
 ******************************************************************************/
void AliMultDepSpecAnalysisTask::InitESDTrackCuts(){

  if(fESDtrackCuts) delete fESDtrackCuts;
  fESDtrackCuts = new AliESDtrackCuts("AliESDtrackCuts");
  if(!fESDtrackCuts) {AliError("fESDtrackCuts not available"); return;}

  fESDtrackCuts->SetRequireTPCRefit(kTRUE);
  fESDtrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  fESDtrackCuts->SetMaxChi2PerClusterTPC(4.0);
  fESDtrackCuts->SetMaxFractionSharedTPCClusters(0.4);
  fESDtrackCuts->SetRequireITSRefit(kTRUE);
  fESDtrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  fESDtrackCuts->SetMaxChi2PerClusterITS(36.);
  fESDtrackCuts->SetDCAToVertex2D(kFALSE);
  fESDtrackCuts->SetRequireSigmaToVertex(kFALSE);
  fESDtrackCuts->SetMaxDCAToVertexZ(2.0);
  fESDtrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01"); // 7*(0.0026+0.0050/pt^1.01)
  fESDtrackCuts->SetAcceptKinkDaughters(kFALSE);
  fESDtrackCuts->SetMaxChi2TPCConstrainedGlobal(36.); // tpcc cut
  fESDtrackCuts->SetCutGeoNcrNcl(3,130,1.5,0.85,0.7); // Geometrical-Length Cut

  // cut-variations:
  if(fCutMode == 101) {fESDtrackCuts->SetMaxChi2PerClusterITS(25.);}
  if(fCutMode == 102) {fESDtrackCuts->SetMaxChi2PerClusterITS(49.);}

  if(fCutMode == 103) {fESDtrackCuts->SetMaxChi2PerClusterTPC(3.0); }
  if(fCutMode == 104) {fESDtrackCuts->SetMaxChi2PerClusterTPC(5.0); }

  if(fCutMode == 105) {fESDtrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.7);}
  if(fCutMode == 106) {fESDtrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.9);}

  if(fCutMode == 107) {fESDtrackCuts->SetMaxFractionSharedTPCClusters(0.2);}
  if(fCutMode == 108) {fESDtrackCuts->SetMaxFractionSharedTPCClusters(1.0);}

  if(fCutMode == 109) {fESDtrackCuts->SetMaxChi2TPCConstrainedGlobal(25.);}
  if(fCutMode == 110) {fESDtrackCuts->SetMaxChi2TPCConstrainedGlobal(49.);}

  if(fCutMode == 111) {fESDtrackCuts->SetMaxDCAToVertexXYPtDep("0.0104+0.0200/pt^1.01");}
  if(fCutMode == 112) {fESDtrackCuts->SetMaxDCAToVertexXYPtDep("0.0260+0.0500/pt^1.01");}

  if(fCutMode == 113) {fESDtrackCuts->SetMaxDCAToVertexZ(1.0);}
  if(fCutMode == 114) {fESDtrackCuts->SetMaxDCAToVertexZ(5.0);}

  if(fCutMode == 115) {fESDtrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kOff);}

  if(fCutMode == 116) {fESDtrackCuts->SetCutGeoNcrNcl(3,120,1.5,0.85,0.7);}
  if(fCutMode == 117) {fESDtrackCuts->SetCutGeoNcrNcl(3,140,1.5,0.85,0.7);}

  if(fCutMode == 118) {fESDtrackCuts->SetCutGeoNcrNcl(4,130,1.5,0.85,0.7);}
  if(fCutMode == 119) {fESDtrackCuts->SetCutGeoNcrNcl(2,130,1.5,0.85,0.7);}
}


/***************************************************************************//**
 * Function to get array of equidistant bin edges between lower and upper edge.
 ******************************************************************************/
void AliMultDepSpecAnalysisTask::SetFixedBinEdges(Double_t* array, Double_t lowerEdge, Double_t upperEdge, Int_t nBins){
  for(Int_t i = 0; i <= nBins; i++){
    array[i] = lowerEdge + i*(upperEdge - lowerEdge)/nBins;
  }
}

/***************************************************************************//**
 * Function to create THnSparseF histogram with the specified axes.
 ******************************************************************************/
THnSparseF* AliMultDepSpecAnalysisTask::CreateHistogram(const string& name, const vector<string>& axes){
  Int_t nAxes = axes.size();
  if(nAxes > MAX_HISTO_DIM) return nullptr;

  Int_t nBins[MAX_HISTO_DIM] = {0};
  Double_t lowerBounds[MAX_HISTO_DIM] = {0.0};
  Double_t upperBounds[MAX_HISTO_DIM] = {0.0};

  string title = name + " [";
  // first figure out number of bins and dimensions
  for(Int_t i = 0; i < nAxes; i++){
    TArrayD* binEdges = GetBinEdges(axes[i]);
    nBins[i] = binEdges->GetSize()-1;
    lowerBounds[i] = binEdges->GetAt(0);
    upperBounds[i] = binEdges->GetAt(binEdges->GetSize()-1);
    title += axes[i];
    if(i < nAxes-1) title += " : "; else title += "]";
  }
  // create histogram
  THnSparseF* histogram = new THnSparseF(name.c_str(), title.c_str(), nAxes, nBins, lowerBounds, upperBounds);

  // set histogram axes
  for(Int_t i = 0; i < nAxes; i++){
    TArrayD* binEdges = GetBinEdges(axes[i]);
    histogram->SetBinEdges(i, binEdges->GetArray());
    histogram->GetAxis(i)->SetTitle(GetAxisTitle(axes[i]).c_str());
    histogram->GetAxis(i)->SetName((std::to_string(i) + "-" + axes[i]).c_str());
  }
  histogram->Sumw2();
  return histogram;
}

/***************************************************************************//**
 * Function to obtain the correct binning for the respective axis.
 ******************************************************************************/
TArrayD* AliMultDepSpecAnalysisTask::GetBinEdges(const string& axisName){
       if(axisName.find("sigmapt") != string::npos) return fBinsPtReso;
  else if(axisName.find("deltapt") != string::npos) return fBinsPtReso;
  else if(axisName.find("pt") != string::npos)      return fBinsPt;
  else if(axisName.find("eta") != string::npos)     return fBinsEta;
  else if(axisName.find("mult") != string::npos)    return fBinsMult;
  else if(axisName.find("cent") != string::npos)    return fBinsCent;
  else if(axisName.find("zv") != string::npos)      return fBinsZv;
  else if(axisName.find("eventcuts") != string::npos)      return fBinsEventCuts;
  else return nullptr;
}

/***************************************************************************//**
 * Function to get the correct title for each histogram axis.
 ******************************************************************************/
string AliMultDepSpecAnalysisTask::GetAxisTitle(const string& axisName){
       if(axisName == "pt")         return "#it{p}_{T} (GeV/#it{c})";
  else if(axisName == "deltapt")    return "#Delta(#it{p}_{T}) / #it{p}^{ meas}_{T}";
  else if(axisName == "mult")       return "Multiplicity";
  else if(axisName == "cent")       return "Centrality (%)";
  else if(axisName == "eta_meas")   return "#eta^{ meas}";
  else if(axisName == "eta_true")   return "#eta^{ true}";
  else if(axisName == "pt_meas")    return "#it{p}^{ meas}_{T} (GeV/#it{c})";
  else if(axisName == "pt_true")    return "#it{p}^{ true}_{T} (GeV/#it{c})";
  else if(axisName == "mult_meas")  return "#it{N}^{ meas}_{ch}";
  else if(axisName == "mult_true")  return "#it{N}^{ true}_{ch}";
  else if(axisName == "sigmapt")    return "#sigma(#it{p}^{ meas}_{T}) / #it{p}^{ meas}_{T}";
  else if(axisName == "eventcuts")    return "[No cuts; Trigger selection; Event selection; Vertex reconstruction and quality; Vertex position; Track selection; triggered and vertex position]";
  else                              return "dummyTitle";
}

/***************************************************************************//**
 * Function to fill a histogram.
 ******************************************************************************/
void AliMultDepSpecAnalysisTask::FillHisto(THnSparseF* histo, const array<Double_t, MAX_HISTO_DIM>& values){
  histo->Fill(values.data());
}

/***************************************************************************//**
 * Function to create a log histogram.
 ******************************************************************************/
TH1D* AliMultDepSpecAnalysisTask::CreateLogHistogram(const string& name){
  return new TH1D(name.c_str(), name.c_str(), 1, 0, 1);
}

/***************************************************************************//**
 * Function to fill a log histogram.
 ******************************************************************************/
void AliMultDepSpecAnalysisTask::FillLogHisto(TH1D* logHist, const string& entry){
  logHist->Fill(entry.c_str(), 1);
}

/***************************************************************************//**
 * Function to set variable binning for multiplicity.
 ******************************************************************************/
void AliMultDepSpecAnalysisTask::SetBinsMult(vector<Int_t> multSteps, vector<Int_t> multBinWidth)
{
    if(multSteps.size() != multBinWidth.size())
    {
      AliError("SetBinsMult:: Vectors need to have same size!");
      return;
    }

    Int_t nMultSteps = multSteps.size();
    Int_t nBinsMult = 1; // for mult=0 bin
    for(Int_t multBins : multSteps) nBinsMult += multBins;
    Double_t* multBinEdges = new Double_t [nBinsMult+1]; // edges need one more

    multBinEdges[0] = -0.5;
    multBinEdges[1] = 0.5;
    Int_t startBin = 1;
    Int_t endBin = 1;
    for(Int_t multStep = 0; multStep < nMultSteps; multStep++){
      endBin += multSteps[multStep];
      for(Int_t multBin = startBin; multBin < endBin; multBin++)  {
        multBinEdges[multBin+1] = multBinEdges[multBin] + multBinWidth[multStep];
      }
      startBin = endBin;
    }
    SetBinsMult(nBinsMult, multBinEdges);
}

/***************************************************************************//**
 * Function to set multiplicity binning for expected maximum multiplicity maxMult.
 * Single multiplicity steps are used for maxMult < 500.
 * For larger multiplicities binning is adjusted such that the maximum of bins
 * is 500 and the first 100 bins are in single multiplicity steps.
 ******************************************************************************/
void AliMultDepSpecAnalysisTask::SetBinsMult(Int_t maxMult)
{
  // for more than 500 bins output becomes large and unfolding takes too long
  if(maxMult > MAX_ALLOWED_MULT_BINS)
  {
    // use single multiplicity for first 100 bins
    // calculate appropriate bin with for the rest
    Int_t nSingleBins = 100;
    Int_t remainingBins = MAX_ALLOWED_MULT_BINS - nSingleBins;

    // increase max mult in case uneven bin width is not possible
    Int_t stepWidth2 = 1;
    while(remainingBins * stepWidth2 < (maxMult - nSingleBins)) stepWidth2 += 2;
    SetBinsMult({100, remainingBins}, {1, stepWidth2});
  }
  else
  {
    SetBinsMult({maxMult}, {1});
  }
}

/***************************************************************************//**
 * Function to add this task to a train.
 ******************************************************************************/
AliMultDepSpecAnalysisTask* AliMultDepSpecAnalysisTask::AddTaskMultDepSpec(TString controlstring, Int_t cutModeLow, Int_t cutModeHigh, Bool_t useDataDrivenCorrections, string pccTrainOutputPath, Int_t pccSysFlag, Int_t secSysFlag)
{
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskMultDepSpec", "No analysis manager found.");
    return nullptr;
  }
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskMultDepSpec", "No input event handler found.");
    return nullptr;
  }
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  Bool_t isMC = (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler() != nullptr);


  // Default cut settings:
  UInt_t triggerMask = AliVEvent::kMB | AliVEvent::kINT7;

  Double_t cutVertexZ = 10.0;

  Double_t cutCentLow = -1000.0;
  Double_t cutCentHigh = 1000.0;

  Double_t cutPtLow = 0.15;
  Double_t cutPtHigh = 50.0;

  Double_t cutEtaLow = -0.8;
  Double_t cutEtaHigh = 0.8;

  Int_t maxMult = 100;
  string colsys = "pp";

  Bool_t useCent = kFALSE;
  if(controlstring.Contains("useCent")) useCent = kTRUE;
  Double_t centBinEdges[9] = {0., 5., 10., 20., 40., 60., 80., 90., 100.};

  Bool_t overridePbPbEventCuts = kFALSE;
  Bool_t useZDC = kFALSE;


  // colison system specific settings
  if(controlstring.Contains("pp")){
    colsys = "pp";
    maxMult = 100;
  }
  else if(controlstring.Contains("pPb"))  {
    colsys = "pPb";
    maxMult = 200;
  }
  else if(controlstring.Contains("XeXe")){
    colsys = "XeXe";
    maxMult = 3700;
    overridePbPbEventCuts = kTRUE;
    useZDC = kTRUE;
  }
  else if(controlstring.Contains("PbPb")) {
    colsys = "PbPb";
    maxMult = 4500;
    overridePbPbEventCuts = kTRUE;
    useZDC = kTRUE;
  }

  if(controlstring.Contains("autoEventCutsPbPb")) overridePbPbEventCuts = kFALSE;
  if(controlstring.Contains("noZDC")) useZDC = kFALSE;


  AliMCSpectraWeights* pccWeights = nullptr;
  if(isMC && pccTrainOutputPath != ""){
    pccWeights = new AliMCSpectraWeights(colsys.c_str(), "fMCSpectraWeights", (AliMCSpectraWeights::SysFlag) pccSysFlag);
    pccWeights->SetMCSpectraFile(pccTrainOutputPath.c_str());
    pccWeights->Init();
  }

  AliMultDepSpecAnalysisTask* returnTask = nullptr;

  char taskName[100] = "";

  for(Int_t cutMode = cutModeLow; cutMode <= cutModeHigh; cutMode++){
    sprintf(taskName, "multDepSpec_%s_cutMode_%d", colsys.c_str(), cutMode);

    AliMultDepSpecAnalysisTask* task = new AliMultDepSpecAnalysisTask(taskName);
    if(cutMode == cutModeLow) returnTask = task; // return one of the tasks

    task->SetSecScalingSysFlag(secSysFlag); //-1, 0, 1

    task->SetUseDataDrivenCorrections(useDataDrivenCorrections);
    task->SetCutMode(cutMode);
    task->SetTriggerMask(triggerMask);

    task->SetIsMC(isMC);
    if(type.Contains("ESD")) task->SetUseESD();
    else task->SetUseAOD();

    task->SetBinsMult(maxMult);
    task->SetUseZDCCut(useZDC);
    task->SetOverridePbPbEventCuts(overridePbPbEventCuts);

    if(useCent){
      task->SetMinCent(cutCentLow);
      task->SetMaxCent(cutCentHigh);
      task->SetBinsCent(8, centBinEdges);
    }

    // kinematic cuts:
    task->SetMinEta(cutEtaLow);
    task->SetMaxEta(cutEtaHigh);
    task->SetMinPt(cutPtLow);
    task->SetMaxPt(cutPtHigh);
    task->SetMaxZv(cutVertexZ);

    task->SetMCSpectraWeights(pccWeights);

    // hang task in train
    mgr->AddTask(task);
    mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task, 1, mgr->CreateContainer(taskName, TList::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root"));
  }
  return returnTask;
}
