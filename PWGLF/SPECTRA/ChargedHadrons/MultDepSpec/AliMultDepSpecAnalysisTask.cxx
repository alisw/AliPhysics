#include "AliMultDepSpecAnalysisTask.h"
/// \cond CLASSIMP
ClassImp(AliMultDepSpecAnalysisTask);
/// \endcond

 //****************************************************************************************
 /**
  * ROOT I/O Constructor.
  */
 //****************************************************************************************
AliMultDepSpecAnalysisTask::AliMultDepSpecAnalysisTask() : AliAnalysisTaskSE(),
  //General member variables
  fOutputList(nullptr),
  fEventCuts(),
  fTrackCuts(nullptr),
  fRand(nullptr),
  fTrainMetadata(""),
  //Toggles
  fIsESD(true),
  fIsMC(false),
  fUseZDCCut(false),
  fOverridePbPbEventCuts(false),
  fMCUseDDC(false),
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
  fHistMCRelPtReso(nullptr),
  fHistMCMultCorrelMatrix(nullptr),
  fHistMCPtCorrelMatrix(nullptr),
  fHistMCEtaCorrelMatrix(nullptr),
  fHistMCPrimTrue(nullptr),
  fHistMCPrimMeas(nullptr),
  fHistMCSecMeas(nullptr),
  // transient event and track properties
  fEvent(nullptr),
  fMCEvent(nullptr),
  fMultMeas(0),
  fMultTrue(0),
  fRunNumber(0),
  fEventNumber(0),
  fTimeStamp(0),
  fPt(0),
  fEta(0),
  fSigmaPt(0),
  fMCPt(0),
  fMCEta(0),
  fMCLabel(0),
  fIsParticleInAcceptance(false),
  fMCIsChargedPrimary(false),
  fMCIsChargedSecDecay(false),
  fMCIsChargedSecMat(false),
  fMCIsChargedSecondary(false),
  fMCParticleWeight(1.0),
  fMCSecScaleWeight(1.0),
  fNRepetitions(1),
  fUseRandomSeed(false)
{
  // ROOT IO constructor, don't allocate memory here!
}

//****************************************************************************************
/**
 * Constructor.
 */
//****************************************************************************************
AliMultDepSpecAnalysisTask::AliMultDepSpecAnalysisTask(const char* name) : AliAnalysisTaskSE(name),
  //General member variables
  fOutputList(nullptr),
  fEventCuts(),
  fTrackCuts(nullptr),
  fRand(nullptr),
  fTrainMetadata(""),
  //Toggles
  fIsESD(true),
  fIsMC(false),
  fUseZDCCut(false),
  fOverridePbPbEventCuts(false),
  fMCUseDDC(false),
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
  fHistMCRelPtReso(nullptr),
  fHistMCMultCorrelMatrix(nullptr),
  fHistMCPtCorrelMatrix(nullptr),
  fHistMCEtaCorrelMatrix(nullptr),
  fHistMCPrimTrue(nullptr),
  fHistMCPrimMeas(nullptr),
  fHistMCSecMeas(nullptr),
  // transient event and track properties
  fEvent(nullptr),
  fMCEvent(nullptr),
  fMultMeas(0),
  fMultTrue(0),
  fRunNumber(0),
  fEventNumber(0),
  fTimeStamp(0),
  fPt(0),
  fEta(0),
  fSigmaPt(0),
  fMCPt(0),
  fMCEta(0),
  fMCLabel(0),
  fIsParticleInAcceptance(false),
  fMCIsChargedPrimary(false),
  fMCIsChargedSecDecay(false),
  fMCIsChargedSecMat(false),
  fMCIsChargedSecondary(false),
  fMCParticleWeight(1.0),
  fMCSecScaleWeight(1.0),
  fNRepetitions(1),
  fUseRandomSeed(false)
{
  // Set default binning
  double binsEventCutsDefault[8] = {-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.6};
  fBinsEventCuts = new TArrayD(8, binsEventCutsDefault);

  double binsMultDefault[2] = {0., 10000.};
  double binsPtDefault[53] = {0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.5,5.0,5.5,6.0,6.5,7.0,8.0,9.0,10.0,20.0,30.0,40.0,50.0,60.0};
  double binsEtaDefault[19] = {-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9};
  double binsZvDefault[13] = {-30.,-25.,-20.,-15.,-10.,-5.,0.,5.,10.,15.,20.,25.,30.};

  // binning for relative pT resolution
  const int nBinsPtReso = 300;
  double binsPtReso[nBinsPtReso+1];
  SetFixedBinEdges(binsPtReso, 0., 0.3, nBinsPtReso);
  SetBinsPtReso(nBinsPtReso, binsPtReso);

  SetBinsMult(1, binsMultDefault);
  SetBinsPt(52, binsPtDefault);
  SetBinsEta(18, binsEtaDefault);
  SetBinsZv(12, binsZvDefault);

  DefineOutput(1, TList::Class());
}

//****************************************************************************************
/**
 * Function executed once before the event loop. Create histograms here.
 */
//****************************************************************************************
void AliMultDepSpecAnalysisTask::UserCreateOutputObjects(){
  
  OpenFile(1, "recreate");
  fOutputList = new TList();
  fOutputList->SetOwner();

  // some meta info histograms filled only once for each train
  fHistEventSelection = CreateHistogram("fHistEventSelection", {"eventcuts"});
  fOutputList->Add(fHistEventSelection);

  //fHistRuns = CreateHistogram("fHistRuns", {"runs"});
  //fOutputList->Add(fHistRuns);
  
  fHistEvents = CreateHistogram("fHistEvents", {"mult_meas"});
  fOutputList->Add(fHistEvents);
  fHistTracks = CreateHistogram("fHistTracks", {"pt_meas", "eta_meas", "mult_meas"});
  fOutputList->Add(fHistTracks);

  fHistRelPtReso = CreateHistogram("fHistRelPtReso", {"sigmapt", "pt_meas"});
  fOutputList->Add(fHistRelPtReso);

  if(fIsMC)
  {
    fHistMCEventEfficiency = CreateHistogram("fHistMCEventEfficiency", {"eventcuts", "mult_true"});
    fOutputList->Add(fHistMCEventEfficiency);

    fHistMCRelPtReso = CreateHistogram("fHistMCRelPtReso", {"deltapt", "pt_meas"});
    fOutputList->Add(fHistMCRelPtReso);
    fHistMCMultCorrelMatrix = CreateHistogram("fHistMCMultCorrelMatrix", {"mult_meas", "mult_true"});
    fOutputList->Add(fHistMCMultCorrelMatrix);
    fHistMCPtCorrelMatrix = CreateHistogram("fHistMCPtCorrelMatrix", {"pt_meas", "pt_true"});
    fOutputList->Add(fHistMCPtCorrelMatrix);
    fHistMCEtaCorrelMatrix = CreateHistogram("fHistMCEtaCorrelMatrix", {"eta_meas", "eta_true"});
    fOutputList->Add(fHistMCEtaCorrelMatrix);
    fHistMCPrimTrue = CreateHistogram("fHistMCPrimTrue", {"pt_true", "eta_true", "mult_true"});
    fOutputList->Add(fHistMCPrimTrue);
    fHistMCPrimMeas = CreateHistogram("fHistMCPrimMeas", {"pt_true", "eta_true", "mult_true"});
    fOutputList->Add(fHistMCPrimMeas);
    fHistMCSecMeas = CreateHistogram("fHistMCSecMeas", {"pt_true", "eta_true", "mult_true"});
    fOutputList->Add(fHistMCSecMeas);
  }

  // override event automatic event selection settings
  // this is needed because AliEventCuts by default throws away all events beyond 90% cent
  if(fOverridePbPbEventCuts)
  {
    fEventCuts.SetManualMode();
    fEventCuts.SetupLHC15o(); // first set default values and then override
    fEventCuts.SetCentralityRange(fMinCent, fMaxCent);
    fEventCuts.fUseEstimatorsCorrelationCut = false;
  }
  //fEventCuts.SetMaxVertexZposition(fMaxZv); // has no effect in automatic mode...
  fEventCuts.OverrideAutomaticTriggerSelection(fTriggerMask);

  fRand = new TRandom3();

  PostData(1, fOutputList);
}

//****************************************************************************************
/**
 * Destructor
 */
//****************************************************************************************
AliMultDepSpecAnalysisTask::~AliMultDepSpecAnalysisTask(){
  if(fRand){delete fRand; fRand = nullptr;}
}

//****************************************************************************************
/**
 * Function executed for each event.
 */
//****************************************************************************************
void AliMultDepSpecAnalysisTask::UserExec(Option_t *){
  if(!InitEvent()) return;
  FillEventHistos();
  LoopMeas();
  if(fIsMC) LoopTrue();
  PostData(1, fOutputList);
}

//****************************************************************************************
/**
 * Function executed after all events were processed.
 */
//****************************************************************************************
void AliMultDepSpecAnalysisTask::Terminate(Option_t*)
{
}

//****************************************************************************************
/**
 * Function to get data-driven secondary scaling weights.
 */
//****************************************************************************************
double AliMultDepSpecAnalysisTask::GetSecScalingFactor(AliVParticle* particle)
{
  // FIXME: this is still missing
  return 1.0;
}

//****************************************************************************************
/**
 * Function to get data-driven particle composition weights.
 */
//****************************************************************************************
double AliMultDepSpecAnalysisTask::GetParticleWeight(AliVParticle* particle)
{
  // FIXME: this will be similar to AliMultSelection
  return 1.0;
}

//****************************************************************************************
/**
 * Initialize event quantities. Sets fEvent, fMCEvent, fMeasMult, fTrueMult
 */
//****************************************************************************************
bool AliMultDepSpecAnalysisTask::InitEvent()
{
    fEvent = InputEvent();
    if (!fEvent) {AliError("fEvent not available\n"); return false;}

    if(fIsMC){
      fMCEvent = MCEvent();
      if (!fMCEvent) {AliError("fMCEvent not available\n"); return false;}
    }

    //v0-mult: fEvent->GetVZEROData()->GetMTotV0A() + fEvent->GetVZEROData()->GetMTotV0C();
    
    // event info for random seeds
    fRunNumber = fEvent->GetRunNumber();
    fTimeStamp = fEvent->GetTimeStamp();
    fEventNumber = fEvent->GetHeader()->GetEventIdAsLong();
  
    if(fUseZDCCut)
    {
      double znaEnergy = fEvent->GetZDCN2Energy();
      double zncEnergy = fEvent->GetZDCN1Energy();
      if ( (znaEnergy <  3503.0) || (zncEnergy <  3609.0) ) {return false;}
    }
    bool acceptEvent = fEventCuts.AcceptEvent(fEvent);

    LoopMeas(true); // set measured multiplicity fMeasMult
    if(fIsMC) LoopTrue(true); // set true multiplicity fTrueMult

    // fill histograms for event selection and Nch dependent efficiency
    std::array <AliEventCuts::NormMask, 5> norm_masks {
      AliEventCuts::kAnyEvent,
      AliEventCuts::kTriggeredEvent,
      AliEventCuts::kPassesNonVertexRelatedSelections,
      AliEventCuts::kHasReconstructedVertex,
      AliEventCuts::kPassesAllCuts // => acceptEvent = true
    };
    for (int iC = 0; iC < 5; ++iC) {
      if (fEventCuts.CheckNormalisationMask(norm_masks[iC])) {
        FillHisto(fHistEventSelection, {double(iC)});
        if(fIsMC){
          FillHisto(fHistMCEventEfficiency, {double(iC), fMultTrue});
        }
      }
    }
    // now count only the events which also contribute to the measurement
    if(acceptEvent){
      if (fMultMeas > 0) FillHisto(fHistEventSelection, {5.0});
      if(fIsMC){
        if (fMultMeas > 0) FillHisto(fHistMCEventEfficiency, {5.0, fMultTrue});
      }
    }
    // additional info: triggered and vertex in acceptacne
    if(fEventCuts.CheckNormalisationMask(AliEventCuts::kTriggeredEvent) && fEventCuts.PassedCut(AliEventCuts::kVertexPosition))
    {
      FillHisto(fHistEventSelection, {6.0});
      if(fIsMC){
        FillHisto(fHistMCEventEfficiency, {6.0, fMultTrue});
      }
    }

    return acceptEvent;
}

//****************************************************************************************
/**
 * Fill event histograms. Event related members are set.
 */
//****************************************************************************************
void AliMultDepSpecAnalysisTask::FillEventHistos()
{
  FillHisto(fHistEvents, {fMultMeas});
  if(fIsMC){
    FillHisto(fHistMCMultCorrelMatrix, {fMultMeas, fMultTrue});
  }
}

//****************************************************************************************
/**
 * Fill track histograms. Track related members are set.
 */
//****************************************************************************************
void AliMultDepSpecAnalysisTask::FillMeasTrackHistos()
{
  FillHisto(fHistTracks, {fPt, fEta, fMultMeas});
  FillHisto(fHistRelPtReso, {fSigmaPt, fPt});
}

//****************************************************************************************
/**
 * Fill measured particle histograms. Track and mc particle related members are set.
 */
//****************************************************************************************
void AliMultDepSpecAnalysisTask::FillMeasParticleHistos()
{
  if(fIsParticleInAcceptance)
  {
    FillHisto(fHistMCRelPtReso, {TMath::Abs(fPt - fMCPt)/fPt, fPt});

    if(fMCIsChargedPrimary)
    {
      FillHisto(fHistMCPtCorrelMatrix, {fPt, fMCPt});
      FillHisto(fHistMCEtaCorrelMatrix, {fEta, fMCEta});
      FillHisto(fHistMCPrimMeas, {fMCPt, fMCEta, fMultTrue});
    }else{
      FillHisto(fHistMCSecMeas, {fMCPt, fMCEta, fMultTrue});
    }
  }
}

//****************************************************************************************
/**
 * Fill generated particle histograms. MC particle related members are set.
 */
//****************************************************************************************
void AliMultDepSpecAnalysisTask::FillTrueParticleHistos()
{
  if(fMCIsChargedPrimary && fIsParticleInAcceptance) FillHisto(fHistMCPrimTrue, {fMCPt, fMCEta, fMultTrue});
}

//****************************************************************************************
/**
 * Loop over measured tracks. Can either count multiplicity or fill histograms.
 */
//****************************************************************************************
void AliMultDepSpecAnalysisTask::LoopMeas(bool count)
{
  if(count) {fMultMeas = 0;}
  AliVTrack* track = nullptr;
  for (int i = 0; i < fEvent->GetNumberOfTracks(); i++){
    track = dynamic_cast<AliVTrack*>(fEvent->GetTrack(i));
    // Set fPt, fEta, fSigmapt; Check if track in kin range and has good quality
    if(!InitTrack(track)) continue;

    // initialize particle properties
    if(fIsMC)
    {
      // mc lable corresponding to measured track
      // negative lable indicates bad quality track - use it anyway
      int mcLable = TMath::Abs(track->GetLabel());

      // Set fMCPt, fMCEta, fMCIsChargedPrimary; Stop computation if it is not charged primary or secondary
      if(fIsESD){
        if(!InitParticle((AliMCParticle*)fMCEvent->GetTrack(mcLable))) continue;
      }else{
        if(!InitParticle((AliAODMCParticle*)fMCEvent->GetTrack(mcLable))) continue;
      }
    }

    if(count)
    {
      fMultMeas += fNRepetitions;
    }
    else
    {
      for(int i = 0; i < fNRepetitions; i++)
      {
        FillMeasTrackHistos();
        if(fIsMC) FillMeasParticleHistos();
      }
    }
  }
}

//****************************************************************************************
/**
 * Loop over generated mc particles. Can either count multiplicity or fill histograms.
 */
//****************************************************************************************
void AliMultDepSpecAnalysisTask::LoopTrue(bool count)
{
  if(count) {fMultTrue = 0;}

  for(int i = 0; i < fMCEvent->GetNumberOfTracks(); i++) {
    // Sets fMCPt, fMCEta, ... and checks if particle in kin range
    if(fIsESD){
      if(!InitParticle((AliMCParticle*)fMCEvent->GetTrack(i))) continue;
    }else{
      if(!InitParticle((AliAODMCParticle*)fMCEvent->GetTrack(i))) continue;
    }

    // mc truth
    if(count)
    {
      if(fMCIsChargedPrimary && fIsParticleInAcceptance){
         fMultTrue += fNRepetitions;
      }
    }else{
      for(int i = 0; i < fNRepetitions; i++)
      {
        FillTrueParticleHistos();
      }
    }
  }
}

//****************************************************************************************
/**
 * Initializes track properties and returns false if track bad or out of range.
 */
//****************************************************************************************
bool AliMultDepSpecAnalysisTask::InitTrack(AliVTrack* track)
{
  if(!track) {AliError("Track not available\n"); return false;}
  fPt = track->Pt();
  fEta = track->Eta();
  if(fPt  <= fMinPt  + PRECISION)   return false;
  if(fPt  >= fMaxPt  - PRECISION)   return false;
  if(fEta <= fMinEta + PRECISION)   return false;
  if(fEta >= fMaxEta - PRECISION)   return false;
  if(!AcceptTrackQuality(track))    return false;
  fNRepetitions = 1;

  if(fIsESD){
    fSigmaPt = fPt * TMath::Sqrt(dynamic_cast<AliESDtrack*>(track)->GetSigma1Pt2());
  }else{
    // for AODs this is only possible with massive overhead
    // (cov matrix entries defined in pxpypz space need to be converted back to sigma 1/pt)
    AliAODTrack* aodTrack = dynamic_cast<AliAODTrack*>(track);
    double cov[21] = {0,}, pxpypz[3] = {0,}, xyz[3] = {0,};
    AliExternalTrackParam exParam;
    aodTrack->GetCovMatrix(cov);
    aodTrack->PxPyPz(pxpypz);
    aodTrack->GetXYZ(xyz);
    exParam.Set(xyz,pxpypz,cov,aodTrack->Charge());
    fSigmaPt = fPt * TMath::Sqrt(exParam.GetSigma1Pt2());
  }
  return true;
}

//****************************************************************************************
/**
 * Initializes particle properties and returns false if out of range.
 * Works for AliMCParticles (ESD) and AliAODMCParticles (AOD).
 */
//****************************************************************************************
template<typename Particle_t>
bool AliMultDepSpecAnalysisTask::InitParticle(Particle_t* particle)
{
  if(!particle) {AliError("Particle not available\n"); return false;}

  bool isCharged = ((TMath::Abs(particle->Charge()) > 0.01)) ? true : false;
  fMCIsChargedPrimary = isCharged && particle->IsPhysicalPrimary();
  fMCIsChargedSecDecay = isCharged && particle->IsSecondaryFromWeakDecay();
  fMCIsChargedSecMat = isCharged && particle->IsSecondaryFromMaterial();
  fMCIsChargedSecondary = fMCIsChargedSecDecay || fMCIsChargedSecMat;
  // not interested in anything non-final or non-charged
  if(!(fMCIsChargedPrimary || fMCIsChargedSecondary)) return false;

  fMCPt = particle->Pt();
  fMCEta = particle->Eta();

  fIsParticleInAcceptance = true;

  if((fMCPt  <= fMinPt  + PRECISION)  || (fMCPt  >= fMaxPt  - PRECISION) ||
  (fMCEta <= fMinEta + PRECISION) || (fMCEta >= fMaxEta - PRECISION))
    fIsParticleInAcceptance = false;

  fMCLabel = particle->GetLabel();
    
  fMCParticleWeight = 1.0;
  fMCSecScaleWeight = 1.0;
  fNRepetitions = 1;

  if(fMCUseDDC)
  {
    if(fMCIsChargedPrimary)
    {
      fMCParticleWeight = GetParticleWeight((AliVParticle*)particle);
      fNRepetitions = GetNRepetitons(fMCParticleWeight);
    }
    else if(fMCIsChargedSecondary)
    {
      fMCSecScaleWeight = 1.;
      fNRepetitions = 1;
    }
  }
  return true;
}

//****************************************************************************************
/**
 * Decide how often to repeat particle in MC to match data.
 */
//****************************************************************************************
int AliMultDepSpecAnalysisTask::GetNRepetitons(double scalingFactor)
{
  int nRepetitions = (int)scalingFactor;
  double rest = scalingFactor - nRepetitions;

  fRand->SetSeed(GetSeed());
  nRepetitions += (fRand->Rndm() <= rest) ? 1 : 0;
  return nRepetitions;
}

//****************************************************************************************
/**
 * Define random (but reproducable) seed.
 */
//****************************************************************************************
unsigned long AliMultDepSpecAnalysisTask::GetSeed()
{
    if (fUseRandomSeed) { return 0; }

    unsigned long seed = fEventNumber;
    seed <<= 5;
    seed += fRunNumber;
    seed <<= 2;
    seed += fMCLabel;
    seed <<= 3;
    seed += fTimeStamp;
    return seed;
}

//****************************************************************************************
/**
 * Function to select tracks with required quality.
 */
//****************************************************************************************
bool AliMultDepSpecAnalysisTask::AcceptTrackQuality(AliVTrack* track){
  if(fTrackCuts)
  {
    return fTrackCuts->IsSelected(track);
  }else{
    AliError("No track cuts are defined!"); return true;
  }
}

//****************************************************************************************
/**
 * Function to obtain centrality.
 */
//****************************************************************************************
double AliMultDepSpecAnalysisTask::GetCentrality(AliVEvent* event)
{
  AliMultSelection* multSelection = (AliMultSelection*) fEvent->FindListObject("MultSelection");
  if(!multSelection){AliError("No MultSelection found!"); return 999;}
  return multSelection->GetMultiplicityPercentile("V0M");
}

//****************************************************************************************
/**
 * Function to get array of equidistant bin edges between lower and upper edge.
 */
//****************************************************************************************
void AliMultDepSpecAnalysisTask::SetFixedBinEdges(double* array, double lowerEdge, double upperEdge, int nBins){
  for(int i = 0; i <= nBins; i++){
    array[i] = lowerEdge + i*(upperEdge - lowerEdge)/nBins;
  }
}

//****************************************************************************************
/**
 * Function to create THnSparseF histogram with the specified axes.
 */
//****************************************************************************************
THnSparseF* AliMultDepSpecAnalysisTask::CreateHistogram(const string& name, const vector<string>& axes){
  int nAxes = axes.size();
  if(nAxes > MAX_HISTO_DIM) return nullptr;

  int nBins[MAX_HISTO_DIM] = {0};
  double lowerBounds[MAX_HISTO_DIM] = {0.0};
  double upperBounds[MAX_HISTO_DIM] = {0.0};

  string title = name + " [";
  // first figure out number of bins and dimensions
  for(int i = 0; i < nAxes; i++){
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
  for(int i = 0; i < nAxes; i++){
    TArrayD* binEdges = GetBinEdges(axes[i]);
    histogram->SetBinEdges(i, binEdges->GetArray());
    histogram->GetAxis(i)->SetTitle(GetAxisTitle(axes[i]).c_str());
    histogram->GetAxis(i)->SetName((std::to_string(i) + "-" + axes[i]).c_str());
  }
  //histogram->Sumw2(); // NEVER do this for simple weightless counting!!
  return histogram;
}

//****************************************************************************************
/**
 * Function to obtain the correct binning for the respective axis.
 */
//****************************************************************************************
TArrayD* AliMultDepSpecAnalysisTask::GetBinEdges(const string& axisName){
       if(axisName.find("sigmapt") != string::npos)     return fBinsPtReso;
  else if(axisName.find("deltapt") != string::npos)     return fBinsPtReso;
  else if(axisName.find("pt") != string::npos)          return fBinsPt;
  else if(axisName.find("eta") != string::npos)         return fBinsEta;
  else if(axisName.find("mult") != string::npos)        return fBinsMult;
  else if(axisName.find("zv") != string::npos)          return fBinsZv;
  else if(axisName.find("eventcuts") != string::npos)   return fBinsEventCuts;
  else                                                  return nullptr;
}

//****************************************************************************************
/**
 * Function to get the correct title for each histogram axis.
 */
//****************************************************************************************
string AliMultDepSpecAnalysisTask::GetAxisTitle(const string& axisName){
       if(axisName == "pt")           return "#it{p}_{T} (GeV/#it{c})";
  else if(axisName == "deltapt")      return "#Delta(#it{p}_{T}) / #it{p}^{ meas}_{T}";
  else if(axisName == "mult")         return "Multiplicity";
  else if(axisName == "eta_meas")     return "#eta^{ meas}";
  else if(axisName == "eta_true")     return "#eta";
  else if(axisName == "pt_meas")      return "#it{p}^{ meas}_{T} (GeV/#it{c})";
  else if(axisName == "pt_true")      return "#it{p}_{T} (GeV/#it{c})";
  else if(axisName == "mult_meas")    return "#it{N}^{ meas}_{ch}";
  else if(axisName == "mult_true")    return "#it{N}_{ch}";
  else if(axisName == "sigmapt")      return "#sigma(#it{p}^{ meas}_{T}) / #it{p}^{ meas}_{T}";
  else if(axisName == "eventcuts")    return "[No cuts; Trigger selection; Event selection; Vertex reconstruction and quality; Vertex position; Track selection; triggered and vertex position]";
  else                                return "dummyTitle";
}

//****************************************************************************************
/**
 * Function to fill a histogram.
 */
//****************************************************************************************
void AliMultDepSpecAnalysisTask::FillHisto(THnSparseF* histo, const array<double, MAX_HISTO_DIM>& values){
  histo->Fill(values.data());
}

//****************************************************************************************
/**
 * Function to create a log histogram.
 */
//****************************************************************************************
TH1D* AliMultDepSpecAnalysisTask::CreateLogHistogram(const string& name){
  return new TH1D(name.c_str(), name.c_str(), 1, 0, 1);
}

//****************************************************************************************
/**
 * Function to fill a log histogram.
 */
//****************************************************************************************
void AliMultDepSpecAnalysisTask::FillLogHisto(TH1D* logHist, const string& entry){
  logHist->Fill(entry.c_str(), 1);
}

//****************************************************************************************
/**
 * Function to set variable binning for multiplicity.
 */
//****************************************************************************************
void AliMultDepSpecAnalysisTask::SetBinsMult(vector<int> multSteps, vector<int> multBinWidth)
{
    if(multSteps.size() != multBinWidth.size())
    {
      AliError("SetBinsMult:: Vectors need to have same size!");
      return;
    }

    int nMultSteps = multSteps.size();
    int nBinsMult = 1; // for mult=0 bin
    for(int multBins : multSteps) nBinsMult += multBins;
    double* multBinEdges = new double [nBinsMult+1]; // edges need one more

    multBinEdges[0] = -0.5;
    multBinEdges[1] = 0.5;
    int startBin = 1;
    int endBin = 1;
    for(int multStep = 0; multStep < nMultSteps; multStep++){
      endBin += multSteps[multStep];
      for(int multBin = startBin; multBin < endBin; multBin++)  {
        multBinEdges[multBin+1] = multBinEdges[multBin] + multBinWidth[multStep];
      }
      startBin = endBin;
    }
    SetBinsMult(nBinsMult, multBinEdges);
}

//****************************************************************************************
/**
 * Function to set multiplicity binning for expected maximum multiplicity maxMult.
 * Single multiplicity steps are used for maxMult < 500.
 * For larger multiplicities binning is adjusted such that the maximum of bins
 * is 500 and the first 100 bins are in single multiplicity steps.
 */
//****************************************************************************************
void AliMultDepSpecAnalysisTask::SetBinsMult(int maxMult)
{
  // for more than 500 bins output becomes large and unfolding takes too long
  if(maxMult > MAX_ALLOWED_MULT_BINS)
  {
    // use single multiplicity for first 100 bins
    // calculate appropriate bin with for the rest
    int nSingleBins = 100;
    int remainingBins = MAX_ALLOWED_MULT_BINS - nSingleBins;

    // increase max mult in case uneven bin width is not possible
    int stepWidth2 = 1;
    while(remainingBins * stepWidth2 < (maxMult - nSingleBins)) stepWidth2 += 2;
    SetBinsMult({100, remainingBins}, {1, stepWidth2});
  }
  else
  {
    SetBinsMult({maxMult}, {1});
  }
}

//****************************************************************************************
/**
 * Function to add this task to a train.
 */
//****************************************************************************************
AliMultDepSpecAnalysisTask* AliMultDepSpecAnalysisTask::AddTaskMultDepSpec
 (string dataSet, TString options, int cutModeLow, int cutModeHigh, bool isMC)
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
  bool isAOD = false;
  if(type.Contains("AOD")){
    isAOD = true;
  }
  else{
    // for ESDs isMC can be determined automatically
    isMC = (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler() != nullptr);
  }
  string mode = (isMC) ? "MC" : "Data";
  
  AliMultDepSpecAnalysisTask* returnTask = nullptr;
  char taskName[100] = "";
  for(int cutMode = cutModeLow; cutMode <= cutModeHigh; cutMode++){
    sprintf(taskName, "%s_%s_cutMode_%d", dataSet.c_str(), mode.c_str(), cutMode);
    AliMultDepSpecAnalysisTask* task = new AliMultDepSpecAnalysisTask(taskName);
    if(!task->InitTask(isMC, isAOD, dataSet, options, cutMode))
    {
      delete task;
      task = nullptr;
      continue;
    }
    if(!returnTask) returnTask = task; // return one of the tasks
    task->SaveTrainMetadata();
    
    // hang task in train
    mgr->AddTask(task);
    mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task, 1, mgr->CreateContainer(taskName, TList::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root"));
  }
  return returnTask;
}


//****************************************************************************************
/**
 * Function to extract metadata string for current train run.
 */
//****************************************************************************************
void AliMultDepSpecAnalysisTask::SaveTrainMetadata()
{
   // Save train metadata
   string trainID     = (gSystem->Getenv("TRAIN_ID")) ? gSystem->Getenv("TRAIN_ID") : "";
   string trainRun    = (gSystem->Getenv("TRAIN_RUN_ID")) ? gSystem->Getenv("TRAIN_RUN_ID") : "";
   string aliPhysTag  = (gSystem->Getenv("ALIROOT_VERSION")) ? gSystem->Getenv("ALIROOT_VERSION") : "";
   if(aliPhysTag.find("::") != string::npos) aliPhysTag = aliPhysTag.substr(aliPhysTag.find("::")+2);
   std::map<string, string> trainIdNames =
   {
     // ESD trains
     {"36", "LF_pp"},
     {"37", "LF_pp_MC"},
     {"51", "LF_pPb"},
     {"53", "LF_pPb_MC"},
     {"26", "LF_PbPb"},
     {"27", "LF_PbPb_MC"},
     // AOD trains
     {"39", "LF_pp_AOD"},
     {"38", "LF_pp_MC_AOD"},
     {"72", "LF_pPb_AOD"},
     {"73", "LF_pPb_MC_AOD"},
     {"20", "LF_PbPb_AOD"},
     {"21", "LF_PbPb_MC_AOD"},
   };
   string trainName =  (trainIdNames.find(trainID) == trainIdNames.end()) ? "-" : trainIdNames[trainID];
   fTrainMetadata = trainName + "#" + trainRun + " @ " + aliPhysTag;
}


//****************************************************************************************
/**
 * Initialize task for specific cut mode. Creates ESD track cuts and particle compositon objects wit correct settings.
 * Call this funktion in the AddTask function. The resulting settings will be streamed.
 * In AODs the track selections are based on the best possible conversion of the specified esd cuts (see AliESDTrackCuts::AcceptVTrack()).
 * Due to the nature of the AOD format it is not possible to have an exact 1-to-1conversion for all of the cuts (e.g. geom length cut).
 * In addition, some cut variables are not available in early AOD productions: e.g. the golden chi2 cut can only be applied since August 2016.
 */
//****************************************************************************************
bool AliMultDepSpecAnalysisTask::InitTask(bool isMC, bool isAOD, string dataSet, TString options, int cutMode)
{
  if((!isMC && (cutMode == 99 || cutMode > 119)) || cutMode < 99 || cutMode > 123) return false;
  SetIsMC(isMC);
  SetIsAOD(isAOD);
  if(!SetupTask(dataSet, options)) return false;
  string colSys = dataSet.substr(0, dataSet.find("_"));
  
  if(fTrackCuts) delete fTrackCuts;
  fTrackCuts = new AliESDtrackCuts("AliESDtrackCuts");

  // these are the default track selections (cutMode == 100)
  fTrackCuts->SetRequireTPCRefit(true);
  fTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  fTrackCuts->SetMaxChi2PerClusterTPC(4.0); // might be lower for PbPb2018
  fTrackCuts->SetMaxFractionSharedTPCClusters(0.4);
  fTrackCuts->SetRequireITSRefit(true);
  fTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  fTrackCuts->SetMaxChi2PerClusterITS(36.);
  fTrackCuts->SetDCAToVertex2D(false);
  fTrackCuts->SetRequireSigmaToVertex(false);
  fTrackCuts->SetMaxDCAToVertexZ(2.0);
  fTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01"); // 7 sigma cut, dataset dependent
  fTrackCuts->SetAcceptKinkDaughters(false);
  fTrackCuts->SetMaxChi2TPCConstrainedGlobal(36.); // golden chi2 cut
  fTrackCuts->SetCutGeoNcrNcl(3,130,1.5,0.85,0.7); // geometrical length cut

  // cut-variations:
  if(cutMode == 101) {fTrackCuts->SetMaxChi2PerClusterITS(25.);}
  if(cutMode == 102) {fTrackCuts->SetMaxChi2PerClusterITS(49.);}

  if(cutMode == 103) {fTrackCuts->SetMaxChi2PerClusterTPC(3.0); }
  if(cutMode == 104) {fTrackCuts->SetMaxChi2PerClusterTPC(5.0); }

  if(cutMode == 105) {fTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.7);}
  if(cutMode == 106) {fTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.9);}

  if(cutMode == 107) {fTrackCuts->SetMaxFractionSharedTPCClusters(0.2);}
  if(cutMode == 108) {fTrackCuts->SetMaxFractionSharedTPCClusters(1.0);}

  if(cutMode == 109) {fTrackCuts->SetMaxChi2TPCConstrainedGlobal(25.);}
  if(cutMode == 110) {fTrackCuts->SetMaxChi2TPCConstrainedGlobal(49.);}

  if(cutMode == 111) {fTrackCuts->SetMaxDCAToVertexXYPtDep("0.0104+0.0200/pt^1.01");}
  if(cutMode == 112) {fTrackCuts->SetMaxDCAToVertexXYPtDep("0.0260+0.0500/pt^1.01");}

  if(cutMode == 113) {fTrackCuts->SetMaxDCAToVertexZ(1.0);}
  if(cutMode == 114) {fTrackCuts->SetMaxDCAToVertexZ(5.0);}

  if(cutMode == 115) {fTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kOff);}

  if(cutMode == 116) {fTrackCuts->SetCutGeoNcrNcl(3,120,1.5,0.85,0.7);}
  if(cutMode == 117) {fTrackCuts->SetCutGeoNcrNcl(3,140,1.5,0.85,0.7);}

  if(cutMode == 118) {fTrackCuts->SetCutGeoNcrNcl(4,130,1.5,0.85,0.7);}
  if(cutMode == 119) {fTrackCuts->SetCutGeoNcrNcl(2,130,1.5,0.85,0.7);}
  
  // in MC we always apply data driven corrections to account for wrong particle composition in the generator
  // cutMode 99 is for crosschecks without any data driven corrections (not part of systematics)
  if(isMC && cutMode != 99)
  {
    // TODO: use enums for this once they are available
    int pccMode = 0;        // 0 = default
    int secScalingMode = 0; // 0 = default

    // systematic variations related to the data driven corrections
    // the particle composition framework picks a new random systematic setup per event
    // do this multiple times to have a better feeling for the systematics
    if(cutMode >= 120 && cutMode <= 123) {pccMode = 1;}
    //if(cutMode >= 124 && cutMode <= 127) {secScalingMode = 1;} // maybe here up or down specifically?

    //SetParticleCompositionMode(pccMode);      // TODO: add this once PCC is ready
    //SetSecondaryScalingMode(secScalingMode);  // TODO: add this once sec scaling is ready

    SetUseDataDrivenCorrections();
  }
  return true;
}

//****************************************************************************************
/**
 * Apply default and data set specific settings like triggers and mulitplicity binning. Override defaults with user options.
 */
//****************************************************************************************
bool AliMultDepSpecAnalysisTask::SetupTask(string dataSet, TString options)
{
  vector<string> dataSets =
  {
    "pp_2TeV",
    "pp_5TeV",
    "pp_7TeV",
    "pp_13TeV",
    "pp_13TeV_trig",
    "pPb_5TeV",
    "pPb_8TeV",
    "XeXe_5TeV",
    "PbPb_2TeV",
    "PbPb_5TeV",
  };
  if(std::find(dataSets.begin(), dataSets.end(), dataSet) == dataSets.end())
  {
    AliError("Settings for specified dataset are not defined!\n");
    return false;
  }
  
  // Default cut settings:
  unsigned int triggerMask = AliVEvent::kMB | AliVEvent::kINT7; // central, semicentral??
  double cutVertexZ = 10.0;
  double cutCentLow = -1000.0;
  double cutCentHigh = 1000.0;
  double cutPtLow = 0.15;
  double cutPtHigh = 50.0;
  double cutEtaLow = -0.8;
  double cutEtaHigh = 0.8;
  int maxMult = 100;

  bool overridePbPbEventCuts = false;
  bool useZDC = false;

  // colison system specific settings
  if(dataSet.find("pp") != string::npos)
  {
    maxMult = 100;
    if(options.Contains("triggered")){
      maxMult = 300;
      triggerMask = AliVEvent::kHighMultV0;
    }
  }
  else if(dataSet.find("pPb") != string::npos)
  {
    maxMult = 200;
  }
  else if(dataSet.find("XeXe") != string::npos){
    maxMult = 3700;
    overridePbPbEventCuts = true;
    useZDC = true;
  }
  else if(dataSet.find("PbPb") != string::npos) {
    maxMult = 4500;
    overridePbPbEventCuts = true;
    useZDC = true;
  }

  if(options.Contains("autoEventCutsPbPb")) overridePbPbEventCuts = false;
  if(options.Contains("noZDC")) useZDC = false;


  // now apply settings
  SetTriggerMask(triggerMask);

  SetBinsMult(maxMult);
  SetUseZDCCut(useZDC);
  SetOverridePbPbEventCuts(overridePbPbEventCuts);

  SetMinCent(cutCentLow);
  SetMaxCent(cutCentHigh);

  // kinematic cuts:
  SetMinEta(cutEtaLow);
  SetMaxEta(cutEtaHigh);
  SetMinPt(cutPtLow);
  SetMaxPt(cutPtHigh);
  SetMaxZv(cutVertexZ);
  
  return true;
}
