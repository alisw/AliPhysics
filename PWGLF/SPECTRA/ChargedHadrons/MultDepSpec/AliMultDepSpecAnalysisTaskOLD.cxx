#include "AliMultDepSpecAnalysisTaskOLD.h"

using std::array;
using std::string;
using std::vector;

//**************************************************************************************************
/**
 * ROOT I/O Constructor.
 */
//**************************************************************************************************
AliMultDepSpecAnalysisTaskOLD::AliMultDepSpecAnalysisTaskOLD() : AliAnalysisTaskSE()
{
  // ROOT IO constructor, don't allocate memory here!
}

//**************************************************************************************************
/**
 * Constructor.
 */
//**************************************************************************************************
AliMultDepSpecAnalysisTaskOLD::AliMultDepSpecAnalysisTaskOLD(const char* name) : AliAnalysisTaskSE(name)
{
  DefineOutput(1, TList::Class());
}

//**************************************************************************************************
/**
 * Define axis that can be used in the histograms.
 */
//**************************************************************************************************
void AliMultDepSpecAnalysisTaskOLD::SetAxis(unsigned int dim, const std::string name,
                                         const std::string title,
                                         const std::vector<double>& binEdges, int nBins)
{
  if(binEdges.size() != 2 && nBins != 0)
  {
    AliError("Specifying the number of bins is only required for fixed bin widths.");
    nBins = 0;
  }
  fAxes[dim] = { name, title, binEdges, nBins };
}

//**************************************************************************************************
/**
 * Define default axis properties .
 */
//**************************************************************************************************
void AliMultDepSpecAnalysisTaskOLD::DefineDefaultAxes(int maxMult)
{
  std::vector<double> ptBins
    = { 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4,  0.45, 0.5,  0.55, 0.6, 0.65, 0.7, 0.75,
        0.8, 0.85, 0.9, 0.95, 1.0, 1.1,  1.2,  1.3,  1.4,  1.5,  1.6, 1.7,  1.8, 1.9,
        2.0, 2.2,  2.4, 2.6,  2.8, 3.0,  3.2,  3.4,  3.6,  3.8,  4.0, 4.5,  5.0, 5.5,
        6.0, 6.5,  7.0, 8.0,  9.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0 };

  if(maxMult > MAX_ALLOWED_MULT_BINS)
  {
    SetAxis(mult_meas, "mult_meas", "#it{N}^{ meas}_{ch}", GetMultBinEdges(maxMult));
    SetAxis(mult_true, "mult_true", "#it{N}_{ch}", GetMultBinEdges(maxMult));
  }
  else
  {
    int nBinsMult = maxMult + 1;
    SetAxis(mult_meas, "mult_meas", "#it{N}^{ meas}_{ch}", { -0.5, nBinsMult - 0.5 }, nBinsMult);
    SetAxis(mult_true, "mult_true", "#it{N}_{ch}", { -0.5, nBinsMult - 0.5 }, nBinsMult);
  }

  int nBinsRelPtReso = 300; // FIXME: too much!
  int nBinsEta = 18;        // FIXME: too much! (maybe remove dimension)

  SetAxis(pt_meas, "pt_meas", "#it{p}_{T} (GeV/#it{c})", ptBins);
  SetAxis(eta_meas, "eta_meas", "#eta^{ meas}", { -0.9, 0.9 }, nBinsEta);
  SetAxis(sigma_pt, "sigma_pt", "#sigma(#it{p}^{ meas}_{T}) / #it{p}^{ meas}_{T}",
          {
            0.,
            0.3,
          },
          nBinsRelPtReso);
  SetAxis(event_cuts, "event_cuts",
          "[No cuts; Trigger selection; Event selection; Vertex reconstruction and quality; Vertex "
          "position; Track selection; triggered and vertex position]",
          { -0.5, 6.5 }, 7);
  SetAxis(zv, "zv", "z vertex position", { -30., 30. }, 12);

  // TODO: add only for mc
  SetAxis(pt_true, "pt_true", "#it{p}_{T} (GeV/#it{c})", ptBins);
  SetAxis(eta_true, "eta_true", "#eta", { -0.9, 0.9 }, nBinsEta);
  SetAxis(delta_pt, "delta_pt", "#Delta(#it{p}_{T}) / #it{p}^{ meas}_{T}",
          {
            0.,
            0.3,
          },
          nBinsRelPtReso);
}

//**************************************************************************************************
/**
 * Function executed once before the event loop. Create histograms here.
 */
//**************************************************************************************************
void AliMultDepSpecAnalysisTaskOLD::BookHistograms()
{
  BookHistogram(fHistEventSelection, "fHistEventSelection", { event_cuts });
  BookHistogram(fHistEvents, "fHistEvents", { mult_meas });
  BookHistogram(fHistTracks, "fHistTracks", { pt_meas, eta_meas, mult_meas });
  BookHistogram(fHistRelPtReso, "fHistRelPtReso", { sigma_pt, pt_meas });

  if(fIsMC)
  {
    BookHistogram(fHistMCEventEfficiency, "fHistMCEventEfficiency", { event_cuts, mult_true });
    BookHistogram(fHistMCRelPtReso, "fHistMCRelPtReso", { delta_pt, pt_meas });
    BookHistogram(fHistMCMultCorrelMatrix, "fHistMCMultCorrelMatrix", { mult_meas, mult_true });
    BookHistogram(fHistMCPtCorrelMatrix, "fHistMCPtCorrelMatrix", { pt_meas, pt_true });
    BookHistogram(fHistMCEtaCorrelMatrix, "fHistMCEtaCorrelMatrix", { eta_meas, eta_true });
    BookHistogram(fHistMCPrimTrue, "fHistMCPrimTrue", { pt_true, eta_true, mult_true });
    BookHistogram(fHistMCPrimMeas, "fHistMCPrimMeas", { pt_true, eta_true, mult_true });
    BookHistogram(fHistMCSecMeas, "fHistMCSecMeas", { pt_true, eta_true, mult_true });
  }

  // check required memory
  if(true)
  {
    double requiredMemory = fHistEventSelection.GetSize() + fHistEvents.GetSize()
                            + fHistTracks.GetSize() + fHistRelPtReso.GetSize()
                            + fHistMCEventEfficiency.GetSize() + fHistMCRelPtReso.GetSize()
                            + fHistMCMultCorrelMatrix.GetSize() + fHistMCPtCorrelMatrix.GetSize()
                            + fHistMCEtaCorrelMatrix.GetSize() + fHistMCPrimTrue.GetSize()
                            + fHistMCPrimMeas.GetSize() + fHistMCSecMeas.GetSize();

    AliError(Form("Estimated memory usage of histograms: %.2f MiB.", requiredMemory / 1048576));
    // Max allowed Memory per train job: 8 GiB
  }
}

//**************************************************************************************************
/**
 * Function executed once before the event loop.
 */
//**************************************************************************************************
void AliMultDepSpecAnalysisTaskOLD::UserCreateOutputObjects()
{

  OpenFile(1, "recreate");
  fOutputList = new TList();
  fOutputList->SetOwner();

  // save train metadata
  fOutputList->Add(fHistTrainInfo.GenerateHist("trainInfo"));
  fOutputList->Add(fHistTrainInfoUE.GenerateHist("trainInfoUE"));

  // book user histograms
  BookHistograms();

  // by default AliEventCuts selects MB and centrality triggered events
  // to avoid confusing results, override this manually
  fEventCuts.OverrideAutomaticTriggerSelection(fTriggerMask);

  // fEventCuts.AddQAplotsToList(fOutputList); // better put in sub-list!
  fRand = new TRandom3();

  PostData(1, fOutputList);
}

//**************************************************************************************************
/**
 * Destructor
 */
//**************************************************************************************************
AliMultDepSpecAnalysisTaskOLD::~AliMultDepSpecAnalysisTaskOLD()
{
  if(fRand)
  {
    delete fRand;
    fRand = nullptr;
  }
}

//**************************************************************************************************
/**
 * Function executed for each event.
 */
//**************************************************************************************************
void AliMultDepSpecAnalysisTaskOLD::UserExec(Option_t*)
{
  if(!InitEvent()) return;
  FillEventHistos();
  LoopMeas();
  if(fIsMC) LoopTrue();
  PostData(1, fOutputList);
}

//**************************************************************************************************
/**
 * Function executed after all events were processed.
 */
//**************************************************************************************************
void AliMultDepSpecAnalysisTaskOLD::Terminate(Option_t*)
{
}

//**************************************************************************************************
/**
 * Function to get data-driven secondary scaling weights.
 */
//**************************************************************************************************
double AliMultDepSpecAnalysisTaskOLD::GetSecScalingFactor(AliVParticle* particle)
{
  // FIXME: this is still missing
  return 1.0;
}

//**************************************************************************************************
/**
 * Function to get data-driven particle composition weights.
 */
//**************************************************************************************************
double AliMultDepSpecAnalysisTaskOLD::GetParticleWeight(AliVParticle* particle)
{
  // FIXME: this will be similar to AliMultSelection
  return 1.0;
}

//**************************************************************************************************
/**
 * Initialize event quantities. Sets fEvent, fMCEvent, fMeasMult, fTrueMult
 */
//**************************************************************************************************
bool AliMultDepSpecAnalysisTaskOLD::InitEvent()
{
  fEvent = InputEvent();
  if(!fEvent)
  {
    AliError("fEvent not available\n");
    return false;
  }

  if(fIsMC)
  {
    fMCEvent = MCEvent();
    if(!fMCEvent)
    {
      AliError("fMCEvent not available\n");
      return false;
    }
  }

  // v0-mult: fEvent->GetVZEROData()->GetMTotV0A() + fEvent->GetVZEROData()->GetMTotV0C();

  // event info for random seeds
  fRunNumber = fEvent->GetRunNumber();
  fTimeStamp = fEvent->GetTimeStamp();
  fEventNumber = fEvent->GetHeader()->GetEventIdAsLong();

  // default check if event is accepted
  bool acceptEvent = fEventCuts.AcceptEvent(fEvent);

  // in PbPb a centrality cut at 90 is always applied by the AliEventCuts (because centrality is
  // only provided up to 90%) to include also more peripheral events, re-do the selection without
  // centrality related cuts for those events this way the event will still be removed if it does
  // not pass the other event quality requirements but the good events can populate the low
  // multiplicities
  if(!acceptEvent && fIncludePeripheralEvents)
  {
    fCent = GetCentrality(fEvent); // TODO: get rid of useless member copying -> InitCent()

    // only consider those events which possibliy were removed by max centrality cut
    // std::cout << "Peripheral event, cent: " << fCent << std::endl;
    if(fCent < 0. || fCent > 90.)
    {
      // since we already ran AcceptEvent once, it is already configured correctly for this
      // event/run/period we only need to override the centrality related cuts before re-running
      fEventCuts.SetManualMode(true);
      fEventCuts.SetCentralityRange(-1000.0, 1000.0);
      fEventCuts.fUseEstimatorsCorrelationCut = false;

      if(fUseZDCCut) // TODO: check if we can do this also in PbPb_2TeV and 2018 PbPb
      {
        double znaEnergy = fEvent->GetZDCN2Energy();
        double zncEnergy = fEvent->GetZDCN1Energy();
        if((znaEnergy < 3503.0) || (zncEnergy < 3609.0))
        {
          return false;
        }
      }

      fIsAcceptedPeripheralEvent = fEventCuts.AcceptEvent(fEvent);
      acceptEvent = fIsAcceptedPeripheralEvent;

      // reset to automatic mode for the next event!
      fEventCuts.SetManualMode(false);
    }
  }

  /*
   PILEUP considerations:

   // for PbPb2018 if one wants to cut out pileup events: (one of the two or both)
   fEventCuts.fUseStrongVarCorrelationCut  = true;  // cut on the V0 multiplicity vs number of
   TPCout traks and reject ~30% of events fEventCuts.fUseVariablesCorrelationCuts = true;  // cut on
   the correlation between the number of SDD+SSD clusters and the number of TPC clusters (may
   introduce some small non-uniformity in the V0M centrality distribution)

   bool fIsPuleupParticle = false;
   if(fRunNumber >= 295581)
   {
    if(fIsESD)
    {
      fIsPuleupParticle = AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(Int_t index,
   AliMCEvent* mcEv) // esd
    }
    else
    {
      fIsPuleupParticle = AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(Int_t index,
   AliAODMCHeader* aodMCHeader, TClonesArray *arrayMC) // aod
    }

   }
   2018 MC per particle mc track removal:

   */

  LoopMeas(true);           // set measured multiplicity fMeasMult
  if(fIsMC) LoopTrue(true); // set true multiplicity fTrueMult

  // fill histograms for event selection and Nch dependent efficiency
  std::array<AliEventCuts::NormMask, 5> norm_masks{
    AliEventCuts::kAnyEvent, AliEventCuts::kTriggeredEvent,
    AliEventCuts::kPassesNonVertexRelatedSelections, AliEventCuts::kHasReconstructedVertex,
    AliEventCuts::kPassesAllCuts // => acceptEvent = true
  };
  for(int iC = 0; iC < 5; ++iC)
  {
    if(fEventCuts.CheckNormalisationMask(norm_masks[iC]))
    {
      fHistEventSelection.Fill(iC);
      if(fIsMC)
      {
        fHistMCEventEfficiency.Fill(iC, fMultTrue);
      }
    }
  }
  // now count only the events which also contribute to the measurement
  if(acceptEvent)
  {
    if(fMultMeas > 0) fHistEventSelection.Fill(5.0);
    if(fIsMC)
    {
      if(fMultMeas > 0) fHistMCEventEfficiency.Fill(5.0, fMultTrue);
    }
  }
  // additional info: triggered and vertex in acceptacne
  if(fEventCuts.CheckNormalisationMask(AliEventCuts::kTriggeredEvent)
     && fEventCuts.PassedCut(AliEventCuts::kVertexPosition))
  {
    fHistEventSelection.Fill(6.0);
    if(fIsMC)
    {
      fHistMCEventEfficiency.Fill(6.0, fMultTrue);
    }
  }
  return acceptEvent;
}

//**************************************************************************************************
/**
 * Fill event histograms. Event related members are set.
 */
//**************************************************************************************************
void AliMultDepSpecAnalysisTaskOLD::FillEventHistos()
{
  // TODO: add this to only in first wagon of train!
  // TODO: find a clever way to log also events per run
  // fill meta train info once per computing job
  if(fIsFirstEventInJob)
  {
    string trainInfo = AliMultSelectionTask::GetSystemTypeByRunNumber(fRunNumber).Data();
    trainInfo += ", ";
    trainInfo += (fIsMC) ? "MC" : "data";
    trainInfo += ", ";
    trainInfo += AliMultSelectionTask::GetPeriodNameByRunNumber(fRunNumber).Data();
    trainInfo += ", ";
    trainInfo += fTrainMetadata;
    fHistTrainInfo.Fill(trainInfo.data());
    // fIsFirstEventInJob = false;
  }

  fHistEvents.Fill(fMultMeas);
  if(fIsMC)
  {
    fHistMCMultCorrelMatrix.Fill(fMultMeas, fMultTrue);
  }
}

//**************************************************************************************************
/**
 * Fill track histograms. Track related members are set.
 */
//**************************************************************************************************
void AliMultDepSpecAnalysisTaskOLD::FillMeasTrackHistos()
{
  fHistTracks.Fill(fPt, fEta, fMultMeas);
  fHistRelPtReso.Fill(fSigmaPt, fPt);
}

//**************************************************************************************************
/**
 * Fill measured particle histograms. Track and mc particle related members are set.
 */
//**************************************************************************************************
void AliMultDepSpecAnalysisTaskOLD::FillMeasParticleHistos()
{
  fHistMCRelPtReso.Fill(TMath::Abs(fPt - fMCPt) / fPt, fPt);

  if(fMCIsChargedPrimary)
  {
    fHistMCPtCorrelMatrix.Fill(fPt, fMCPt);
    fHistMCEtaCorrelMatrix.Fill(fEta, fMCEta);
    fHistMCPrimMeas.Fill(fMCPt, fMCEta, fMultTrue);
  }
  else
  {
    fHistMCSecMeas.Fill(fMCPt, fMCEta, fMultTrue);
  }
}

//**************************************************************************************************
/**
 * Fill generated particle histograms. MC particle related members are set.
 */
//**************************************************************************************************
void AliMultDepSpecAnalysisTaskOLD::FillTrueParticleHistos()
{
  if(fMCIsChargedPrimary) fHistMCPrimTrue.Fill(fMCPt, fMCEta, fMultTrue);
}

//**************************************************************************************************
/**
 * Loop over measured tracks. Can either count multiplicity or fill histograms.
 */
//**************************************************************************************************
void AliMultDepSpecAnalysisTaskOLD::LoopMeas(bool count)
{
  if(count)
  {
    fMultMeas = 0;
    fSumPtMeas = 0;
  }
  AliVTrack* track = nullptr;
  for(int i = 0; i < fEvent->GetNumberOfTracks(); i++)
  {
    track = dynamic_cast<AliVTrack*>(fEvent->GetTrack(i));
    // Set track properties and check if track is in kin range and has good quality
    if(!InitTrack(track)) continue;
    if(!SelectTrack(count)) continue; // apply additional selection criteria on track sample

    // initialize particle properties
    if(fIsMC)
    {
      // mc lable corresponding to measured track
      // negative lable indicates bad quality track - use it anyway
      int mcLable = TMath::Abs(track->GetLabel());

      // Set mc particle properties and check if it is charged prim/sec and in kin range
      if(fIsESD)
      {
        if(!InitParticle((AliMCParticle*)fMCEvent->GetTrack(mcLable))) continue;
      }
      else
      {
        if(!InitParticle((AliAODMCParticle*)fMCEvent->GetTrack(mcLable))) continue;
      }
      if(!SelectParticle(count)) continue; // apply additional selection criteria on particle sample
    }

    if(count)
    {
      fMultMeas += fNRepetitions;
      fSumPtMeas += fPt;
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

//**************************************************************************************************
/**
 * Loop over generated mc particles. Can either count multiplicity or fill histograms.
 */
//**************************************************************************************************
void AliMultDepSpecAnalysisTaskOLD::LoopTrue(bool count)
{
  if(count)
  {
    fMultTrue = 0;
    fSumPtTrue = 0;
  }

  for(int i = 0; i < fMCEvent->GetNumberOfTracks(); i++)
  {
    // Sets fMCPt, fMCEta, ... and checks if particle in kin range
    if(fIsESD)
    {
      if(!InitParticle((AliMCParticle*)fMCEvent->GetTrack(i))) continue;
    }
    else
    {
      if(!InitParticle((AliAODMCParticle*)fMCEvent->GetTrack(i))) continue;
    }
    if(!SelectParticle(count)) continue; // apply additional selection criteria on particle sample

    // mc truth
    if(count)
    {
      if(fMCIsChargedPrimary)
      {
        fMultTrue += fNRepetitions;
        fSumPtTrue += fMCPt;
      }
    }
    else
    {
      for(int i = 0; i < fNRepetitions; i++)
      {
        FillTrueParticleHistos();
      }
    }
  }
}

//**************************************************************************************************
/**
 * Initializes track properties and returns false if track is not available, has bad quality or is
 * not in kinematic range.
 */
//**************************************************************************************************
bool AliMultDepSpecAnalysisTaskOLD::InitTrack(AliVTrack* track)
{
  if(!track)
  {
    AliError("Track not available\n");
    return false;
  }
  fPt = track->Pt();
  fEta = track->Eta();
  fPhi = track->Phi();
  if((fPt <= fMinPt + PRECISION) || (fPt >= fMaxPt - PRECISION) || (fEta <= fMinEta + PRECISION)
     || (fEta >= fMaxEta - PRECISION))
  {
    return false;
  }
  if(!AcceptTrackQuality(track)) return false;
  fNRepetitions = 1;

  if(fIsESD)
  {
    fSigmaPt = fPt * TMath::Sqrt(dynamic_cast<AliESDtrack*>(track)->GetSigma1Pt2());
  }
  else
  {
    // for AODs this is only possible with massive overhead
    // (cov matrix entries defined in pxpypz space need to be converted back to sigma 1/pt)
    AliAODTrack* aodTrack = dynamic_cast<AliAODTrack*>(track);
    double cov[21] = {0,}, pxpypz[3] = {0,}, xyz[3] = {0,};
    AliExternalTrackParam exParam;
    aodTrack->GetCovMatrix(cov);
    aodTrack->PxPyPz(pxpypz);
    aodTrack->GetXYZ(xyz);
    exParam.Set(xyz, pxpypz, cov, aodTrack->Charge());
    fSigmaPt = fPt * TMath::Sqrt(exParam.GetSigma1Pt2());
  }
  return true;
}

//**************************************************************************************************
/**
 * Initializes particle properties and returns false if not charged primary or secondary out if
 * particle is of kinematic range range. Works for AliMCParticles (ESD) and AliAODMCParticles (AOD).
 */
//**************************************************************************************************
template <typename Particle_t>
bool AliMultDepSpecAnalysisTaskOLD::InitParticleBase(Particle_t* particle)
{
  if(!particle)
  {
    AliError("Particle not available\n");
    return false;
  }

  fMCPt = particle->Pt();
  fMCEta = particle->Eta();
  fMCPhi = particle->Phi();

  if((fMCPt <= fMinPt + PRECISION) || (fMCPt >= fMaxPt - PRECISION)
     || (fMCEta <= fMinEta + PRECISION) || (fMCEta >= fMaxEta - PRECISION))
  {
    return false;
  }

  bool isCharged = ((TMath::Abs(particle->Charge()) > 0.01)) ? true : false;
  fMCIsChargedPrimary = isCharged && particle->IsPhysicalPrimary();
  fMCIsChargedSecDecay = isCharged && particle->IsSecondaryFromWeakDecay();
  fMCIsChargedSecMat = isCharged && particle->IsSecondaryFromMaterial();
  fMCIsChargedSecondary = fMCIsChargedSecDecay || fMCIsChargedSecMat;

  // not interested in anything non-final or non-charged
  if(!(fMCIsChargedPrimary || fMCIsChargedSecondary)) return false;

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

//**************************************************************************************************
/**
 * Decide how often to repeat particle in MC to match data.
 */
//**************************************************************************************************
int AliMultDepSpecAnalysisTaskOLD::GetNRepetitons(double scalingFactor)
{
  int nRepetitions = (int)scalingFactor;
  double rest = scalingFactor - nRepetitions;

  fRand->SetSeed(GetSeed());
  nRepetitions += (fRand->Rndm() <= rest) ? 1 : 0;
  return nRepetitions;
}

//**************************************************************************************************
/**
 * Define random (but reproducable) seed.
 */
//**************************************************************************************************
unsigned long AliMultDepSpecAnalysisTaskOLD::GetSeed()
{
  if(fUseRandomSeed)
  {
    return 0;
  }

  unsigned long seed = fEventNumber;
  seed <<= 5;
  seed += fRunNumber;
  seed <<= 2;
  seed += fMCLabel;
  seed <<= 3;
  seed += fTimeStamp;
  return seed;
}

//**************************************************************************************************
/**
 * Function to select tracks with required quality.
 */
//**************************************************************************************************
bool AliMultDepSpecAnalysisTaskOLD::AcceptTrackQuality(AliVTrack* track)
{
  if(fTrackCuts)
  {
    return fTrackCuts->IsSelected(track);
  }
  else
  {
    AliError("No track cuts are defined!");
    return true;
  }
}

//**************************************************************************************************
/**
 * Function to obtain centrality.
 */
//**************************************************************************************************
double AliMultDepSpecAnalysisTaskOLD::GetCentrality(AliVEvent* event)
{
  AliMultSelection* multSelection = (AliMultSelection*)fEvent->FindListObject("MultSelection");
  if(!multSelection)
  {
    // check for legacy cent selection used in old PbPb data
    AliCentrality* centSelection = event->GetCentrality();
    if(!centSelection)
    {
      AliError("Centrality wagon not found!");
      return 999.;
    }
    return centSelection->GetCentralityPercentile("V0M");
  }
  return multSelection->GetMultiplicityPercentile("V0M");

  // TODO: this would be shorter:
  // AliMultSelectionBase::GetMultiplicityPercentileWithFallback(fEvent, "V0M");
}

//**************************************************************************************************
/**
 * Function to set variable binning for multiplicity.
 */
//**************************************************************************************************
std::vector<double> AliMultDepSpecAnalysisTaskOLD::GetMultBinEdges(vector<int> multSteps,
                                                                vector<int> multBinWidth)
{
  if(multSteps.size() != multBinWidth.size())
  {
    AliFatal("Vectors need to have same size!");
    return {};
  }

  int nMultSteps = multSteps.size();
  int nBinsMult = 1; // for mult=0 bin
  for(int multBins : multSteps)
    nBinsMult += multBins;

  std::vector<double> multBinEdges;
  multBinEdges.resize(nBinsMult + 1); // edges need one more

  multBinEdges[0] = -0.5;
  multBinEdges[1] = 0.5;
  int startBin = 1;
  int endBin = 1;
  for(int multStep = 0; multStep < nMultSteps; multStep++)
  {
    endBin += multSteps[multStep];
    for(int multBin = startBin; multBin < endBin; ++multBin)
    {
      multBinEdges[multBin + 1] = multBinEdges[multBin] + multBinWidth[multStep];
    }
    startBin = endBin;
  }
  return multBinEdges;
}

//**************************************************************************************************
/**
 * Function to set multiplicity binning for expected maximum multiplicity maxMult.
 * Single multiplicity steps are used for maxMult < 500.
 * For larger multiplicities binning is adjusted such that the maximum of bins
 * is 500 and the first 100 bins are in single multiplicity steps.
 */
//**************************************************************************************************
std::vector<double> AliMultDepSpecAnalysisTaskOLD::GetMultBinEdges(int maxMult)
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
    while(remainingBins * stepWidth2 < (maxMult - nSingleBins))
      stepWidth2 += 2;
    return GetMultBinEdges({ 100, remainingBins }, { 1, stepWidth2 });
  }
  else
  {
    return GetMultBinEdges({ maxMult }, { 1 });
  }
}

//**************************************************************************************************
/**
 * Function to add this task to a train.
 */
//**************************************************************************************************
AliMultDepSpecAnalysisTaskOLD* AliMultDepSpecAnalysisTaskOLD::AddTaskMultDepSpec(
  const string& dataSet, int cutModeLow, int cutModeHigh, TString options, bool isMC)
{
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr)
  {
    ::Error("AddTaskMultDepSpec", "No analysis manager found.");
    return nullptr;
  }
  if(!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskMultDepSpec", "No input event handler found.");
    return nullptr;
  }
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  bool isAOD = false;
  if(type.Contains("AOD"))
  {
    isAOD = true;
  }
  else
  {
    // for ESDs isMC can be determined automatically
    isMC = (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler() != nullptr);
  }
  string mode = (isMC) ? "MC" : "Data";

  AliMultDepSpecAnalysisTaskOLD* returnTask = nullptr;
  char taskName[100] = "";
  for(int cutMode = cutModeLow; cutMode <= cutModeHigh; cutMode++)
  {
    sprintf(taskName, "%s_%s_cutMode_%d", dataSet.data(), mode.data(), cutMode);
    AliMultDepSpecAnalysisTaskOLD* task = new AliMultDepSpecAnalysisTaskOLD(taskName);
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
    mgr->ConnectOutput(
      task, 1,
      mgr->CreateContainer(taskName, TList::Class(), AliAnalysisManager::kOutputContainer,
                           "AnalysisResults.root"));
  }
  return returnTask;
}

//**************************************************************************************************
/**
 * Function to extract metadata string for current train run.
 */
//**************************************************************************************************
void AliMultDepSpecAnalysisTaskOLD::SaveTrainMetadata()
{
  // Save train metadata
  string trainID = (gSystem->Getenv("TRAIN_ID")) ? gSystem->Getenv("TRAIN_ID") : "";
  string trainRun = (gSystem->Getenv("TRAIN_RUN_ID")) ? gSystem->Getenv("TRAIN_RUN_ID") : "";
  string aliPhysTag
    = (gSystem->Getenv("ALIROOT_VERSION")) ? gSystem->Getenv("ALIROOT_VERSION") : "";
  if(aliPhysTag.find("::") != string::npos)
    aliPhysTag = aliPhysTag.substr(aliPhysTag.find("::") + 2);
  std::map<string, string> trainIdNames = {
    // ESD trains
    { "36", "LF_pp" },
    { "37", "LF_pp_MC" },
    { "51", "LF_pPb" },
    { "53", "LF_pPb_MC" },
    { "26", "LF_PbPb" },
    { "27", "LF_PbPb_MC" },
    // AOD trains
    { "39", "LF_pp_AOD" },
    { "38", "LF_pp_MC_AOD" },
    { "72", "LF_pPb_AOD" },
    { "73", "LF_pPb_MC_AOD" },
    { "20", "LF_PbPb_AOD" },
    { "21", "LF_PbPb_MC_AOD" },
  };
  string trainName
    = (trainIdNames.find(trainID) == trainIdNames.end()) ? "-" : trainIdNames[trainID];
  fTrainMetadata = trainName + "#" + trainRun + " @ " + aliPhysTag;
}

//**************************************************************************************************
/**
 * Initialize task for specific cut mode. Creates ESD track cuts and particle compositon objects wit
 * correct settings. Call this funktion in the AddTask function. The resulting settings will be
 * streamed. In AODs the track selections are based on the best possible conversion of the specified
 * esd cuts (see AliESDTrackCuts::AcceptVTrack()). Due to the nature of the AOD format it is not
 * possible to have an exact 1-to-1conversion for all of the cuts (e.g. geom length cut). In
 * addition, some cut variables are not available in early AOD productions: e.g. the golden chi2 cut
 * can only be applied since August 2016.
 */
//**************************************************************************************************
bool AliMultDepSpecAnalysisTaskOLD::InitTask(bool isMC, bool isAOD, string dataSet, TString options,
                                          int cutMode)
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
  fTrackCuts->SetMaxChi2PerClusterTPC(4.0); // lower for PbPb2018 -> 2.5
  fTrackCuts->SetMaxFractionSharedTPCClusters(0.4);
  fTrackCuts->SetRequireITSRefit(true);
  fTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
  fTrackCuts->SetMaxChi2PerClusterITS(36.);
  fTrackCuts->SetDCAToVertex2D(false);
  fTrackCuts->SetRequireSigmaToVertex(false);
  fTrackCuts->SetMaxDCAToVertexZ(2.0);
  fTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01"); // 7 sigma cut, dataset dependent
  fTrackCuts->SetAcceptKinkDaughters(false);
  fTrackCuts->SetMaxChi2TPCConstrainedGlobal(36.);     // golden chi2 cut
  fTrackCuts->SetCutGeoNcrNcl(3, 130, 1.5, 0.85, 0.7); // geometrical length cut

  // cut-variations:
  if(cutMode == 101)
  {
    fTrackCuts->SetMaxChi2PerClusterITS(25.);
  }
  if(cutMode == 102)
  {
    fTrackCuts->SetMaxChi2PerClusterITS(49.);
  }

  if(cutMode == 103)
  {
    fTrackCuts->SetMaxChi2PerClusterTPC(3.0);
  }
  if(cutMode == 104)
  {
    fTrackCuts->SetMaxChi2PerClusterTPC(5.0);
  }

  if(cutMode == 105)
  {
    fTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.7);
  }
  if(cutMode == 106)
  {
    fTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.9);
  }

  if(cutMode == 107)
  {
    fTrackCuts->SetMaxFractionSharedTPCClusters(0.2);
  }
  if(cutMode == 108)
  {
    fTrackCuts->SetMaxFractionSharedTPCClusters(1.0);
  }

  if(cutMode == 109)
  {
    fTrackCuts->SetMaxChi2TPCConstrainedGlobal(25.);
  }
  if(cutMode == 110)
  {
    fTrackCuts->SetMaxChi2TPCConstrainedGlobal(49.);
  }

  if(cutMode == 111)
  {
    fTrackCuts->SetMaxDCAToVertexXYPtDep("0.0104+0.0200/pt^1.01");
  }
  if(cutMode == 112)
  {
    fTrackCuts->SetMaxDCAToVertexXYPtDep("0.0260+0.0500/pt^1.01");
  }

  if(cutMode == 113)
  {
    fTrackCuts->SetMaxDCAToVertexZ(1.0);
  }
  if(cutMode == 114)
  {
    fTrackCuts->SetMaxDCAToVertexZ(5.0);
  }

  if(cutMode == 115)
  {
    fTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kOff);
  }

  if(cutMode == 116)
  {
    fTrackCuts->SetCutGeoNcrNcl(3, 120, 1.5, 0.85, 0.7);
  }
  if(cutMode == 117)
  {
    fTrackCuts->SetCutGeoNcrNcl(3, 140, 1.5, 0.85, 0.7);
  }

  if(cutMode == 118)
  {
    fTrackCuts->SetCutGeoNcrNcl(4, 130, 1.5, 0.85, 0.7);
  }
  if(cutMode == 119)
  {
    fTrackCuts->SetCutGeoNcrNcl(2, 130, 1.5, 0.85, 0.7);
  }

  // in MC we always apply data driven corrections to account for wrong particle composition in the
  // generator cutMode 99 is for crosschecks without any data driven corrections (not part of
  // systematics)
  if(isMC && cutMode != 99)
  {
    // TODO: use enums for this once they are available
    int pccMode = 0;        // 0 = default
    int secScalingMode = 0; // 0 = default

    // systematic variations related to the data driven corrections
    // the particle composition framework picks a new random systematic setup per event
    // do this multiple times to have a better feeling for the systematics
    if(cutMode >= 120 && cutMode <= 123)
    {
      pccMode = 1;
    }
    // if(cutMode >= 124 && cutMode <= 127) {secScalingMode = 1;} // maybe here up or down
    // specifically?

    // SetParticleCompositionMode(pccMode);      // TODO: add this once PCC is ready
    // SetSecondaryScalingMode(secScalingMode);  // TODO: add this once sec scaling is ready

    SetUseDataDrivenCorrections();
  }
  return true;
}

//**************************************************************************************************
/**
 * Apply default and data set specific settings like triggers and mulitplicity binning. Override
 * defaults with user options.
 */
//**************************************************************************************************
bool AliMultDepSpecAnalysisTaskOLD::SetupTask(string dataSet, TString options)
{
  vector<string> dataSets = {
    "pp_2TeV",        "pp_5TeV",        "pp_7TeV",        "pp_13TeV",
    "pp_13TeV_trig",  "pPb_5TeV",       "pPb_8TeV",       "XeXe_5TeV",
    "PbPb_2TeV",      "PbPb_5TeV",      "XeXe_5TeV_semi", "PbPb_2TeV_semi",
    "PbPb_5TeV_semi", "XeXe_5TeV_cent", "PbPb_2TeV_cent", "PbPb_5TeV_cent",
  };

  if(std::find(dataSets.begin(), dataSets.end(), dataSet) == dataSets.end())
  {
    AliError("Settings for specified dataset are not defined!\n");
    return false;
  }

  // Default cut settings:
  unsigned int triggerMask = AliVEvent::kAnyINT;
  // unsigned int triggerMask = AliVEvent::kMB | AliVEvent::kINT7;
  double cutPtLow = 0.15;
  double cutPtHigh = 50.0;
  double cutEtaLow = -0.8;
  double cutEtaHigh = 0.8;
  int maxMult = 100;

  bool includePeripheralEvents = false;
  bool useZDC = false;

  // colison system specific settings
  if(dataSet.find("pp") != string::npos)
  {
    maxMult = 100;
    if(dataSet.find("trig") != string::npos)
    {
      triggerMask = AliVEvent::kHighMultV0;
    }
  }
  else if(dataSet.find("pPb") != string::npos)
  {
    maxMult = 180;
  }
  else if(dataSet.find("PbPb") != string::npos || dataSet.find("XeXe") != string::npos)
  {
    maxMult = 4500;
    includePeripheralEvents = true;
    useZDC = true;

    if(dataSet.find("semi") != string::npos)
    {
      triggerMask = AliVEvent::kSemiCentral;
    }
    if(dataSet.find("cent") != string::npos)
    {
      triggerMask = AliVEvent::kCentral;
    }
  }
  if(dataSet.find("XeXe") != string::npos)
  {
    maxMult = 3700;
  }

  if(options.Contains("noPeripheral")) includePeripheralEvents = false;
  if(options.Contains("noZDC")) useZDC = false;

  // now apply settings
  SetTriggerMask(triggerMask);

  // SetBinsMult(maxMult);
  SetUseZDCCut(useZDC);
  SetIncludePeripheralEvents(includePeripheralEvents);

  // kinematic cuts:
  SetMinEta(cutEtaLow);
  SetMaxEta(cutEtaHigh);
  SetMinPt(cutPtLow);
  SetMaxPt(cutPtHigh);

  DefineDefaultAxes(maxMult);
  return true;
}
