#include "AliMultDepSpecAnalysisTaskUE.h"

using std::array;
using std::string;
using std::vector;

//**************************************************************************************************
/**
 * ROOT I/O Constructor.
 */
//**************************************************************************************************
AliMultDepSpecAnalysisTaskUE::AliMultDepSpecAnalysisTaskUE()
  : AliMultDepSpecAnalysisTaskOLD()
{
  // ROOT IO constructor, don't allocate memory here!
}

//**************************************************************************************************
/**
 * Constructor.
 */
//**************************************************************************************************
AliMultDepSpecAnalysisTaskUE::AliMultDepSpecAnalysisTaskUE(const char* name)
  : AliMultDepSpecAnalysisTaskOLD(name)
{
}

//**************************************************************************************************
/**
 * Destructor
 */
//**************************************************************************************************
AliMultDepSpecAnalysisTaskUE::~AliMultDepSpecAnalysisTaskUE()
{
  // AliMultDepSpecAnalysisTaskOLD::~AliMultDepSpecAnalysisTaskOLD();
}

//**************************************************************************************************
/**
 * Define default axis properties .
 */
//**************************************************************************************************
void AliMultDepSpecAnalysisTaskUE::DefineDefaultAxes(int maxMult)
{
  // add dimensions used in base class
  AliMultDepSpecAnalysisTaskOLD::DefineDefaultAxes(maxMult);

  // additional dimensions for UE studies
  std::vector<double> ptBins
    = { 0.1, 0.15, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0,  3.2,  3.4,  3.6,
        3.8,  4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0 };

  std::vector<double> sumptBins
    = { 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 3.2, 3.4, 3.6, 3.8, 4.0,
        4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0 };

  SetAxis(pt_lead_meas, "pt_lead_meas", "#it{p}_{T,Lead} (GeV/#it{c})", ptBins);
  SetAxis(sum_pt, "sum_pt", "sum #it{p}_{T} (GeV/#it{c})", sumptBins);
  SetAxis(phi_lead_meas, "phi_lead_meas", "#it{#it{#phi}_{Lead}}", { 0, 2 * TMath::Pi() }, 180);
  SetAxis(delta_pt_lead, "delta_pt_lead", "#Delta(#it{p}_{T,Lead})", { 0., 2. }, 40);
  SetAxis(delta_phi_lead, "delta_phi_lead", "#Delta(#it{#phi}_{Lead})", { 0, TMath::Pi() / 4 },
          45);
  SetAxis(delta_phi, "delta_phi", "#it{#phi}_{Lead} - #it{#phi}_{Track}", { 0, 2 * TMath::Pi() },
          180);
}

//**************************************************************************************************
/**
 * Function executed once before the event loop. Create histograms here.
 */
//**************************************************************************************************
void AliMultDepSpecAnalysisTaskUE::BookHistograms()
{
  // book all default histograms
  AliMultDepSpecAnalysisTaskOLD::BookHistograms();

  // book UE specific histograms
  BookHistogram(fHistLeadPt, "fHistLeadPt", { pt_lead_meas });
  BookHistogram(fHistLeadPhi, "fHistLeadPhi", { phi_lead_meas });
  BookHistogram(fHistLeadPtPhi, "fHistLeadPtPhi", { pt_lead_meas, phi_lead_meas });
  BookHistogram(fHistDeltaPhi, "fHistDeltaPhi", { delta_phi });
  BookHistogram(fHistPtLeadCutLoss, "fHistPtLeadCutLoss", { pt_meas });
  BookHistogram(fHistPlateau, "fHistPlateau", { pt_lead_meas, mult_meas });
  BookHistogram(fHistSumPt, "fHistSumPt", { sum_pt, sum_pt });
  BookHistogram(fHistMultCorr, "fHistMultCorr", { mult_meas, mult_meas });
  BookHistogram(fHistMeanSumPt, "fHistMeanSumPt", { sum_pt });
  BookHistogram(fHistMeanMult, "fHistMeanMult", { mult_meas });
  BookHistogram(fHistMeanSumPtMultTot, "fHistMeanSumPtMultTot", { mult_meas});
  BookHistogram(fHistSumPtMult, "fHistSumPtMult", { mult_meas, sum_pt });


  if(fIsMC)
  {
    BookHistogram(fHistMCResoPtLead, "fHistMCResoPtLead", { delta_pt_lead });
    BookHistogram(fHistMCResoPhiLead, "fHistMCResoPhiLead", { delta_phi_lead });
    BookHistogram(fHistMCPlateau, "fHistMCPlateau", { pt_lead_meas, mult_meas });
    BookHistogram(fHistMCSumPt, "fHistMCSumPt", { sum_pt, sum_pt });
    BookHistogram(fHistMCMultCorr, "fHistMCMultCorr", { mult_meas, mult_meas });
    BookHistogram(fHistMCMeanSumPt, "fHistMCMeanSumPt", { sum_pt });
    BookHistogram(fHistMCMeanMult, "fHistMCMeanMult", { mult_meas });
    BookHistogram(fHistMCMeanSumPtMultTot, "fHistMCMeanSumPtMultTot", { mult_meas});
    BookHistogram(fHistMCSumPtMult, "fHistMCSumPtMult", { mult_meas, sum_pt });
  }

  // check additional required memory
  double requiredMemory = fHistLeadPt.GetSize() + fHistLeadPhi.GetSize() + fHistDeltaPhi.GetSize()
                          + fHistPtLeadCutLoss.GetSize() + fHistPlateau.GetSize()
                          + fHistMCResoPtLead.GetSize() + fHistMCResoPhiLead.GetSize()
                          + fHistSumPt.GetSize() + fHistMultCorr.GetSize()
                          + fHistMeanSumPt.GetSize() + fHistMeanMult.GetSize()
                          + fHistMCMeanSumPt.GetSize() + fHistMCMeanMult.GetSize()
                          + fHistMCPlateau.GetSize() + fHistMCSumPt.GetSize() + fHistMCMultCorr.GetSize()
                          + fHistLeadPtPhi.GetSize() + fHistMeanSumPtMultTot.GetSize()
                          + fHistMCMeanSumPtMultTot.GetSize() + fHistSumPtMult.GetSize()
                          + fHistMCSumPtMult.GetSize();

  AliError(Form("Estimated additional memory usage of histograms from UE analysis: %.2f MiB.",
                requiredMemory / 1048576));
}

//**************************************************************************************************
/**
 * Initialize event quantities. Sets fEvent, fMCEvent, fMeasMult, fTrueMult
 */
//**************************************************************************************************
bool AliMultDepSpecAnalysisTaskUE::InitEvent()
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
  FindLeadingTrack(); // sets fPhiLead, fMCPhilead, fPtLead, fMCPtLead
  if(fPtLead < fPtLeadCut)
  {
    if(fMultMeas != 0) fHistPtLeadCutLoss.Fill(fPtLead);
    return false;
  }
  // if(fIsMC && (fPtLead < fPtLeadCut))
  // {
  //   if(fMultMeas != 0) fHistPtLeadCutLoss.Fill(fPtLead);
  //   return false;
  // }
  return AliMultDepSpecAnalysisTaskOLD::InitEvent();
}

//**************************************************************************************************
/**
 * Finds the leading track.
 */
//**************************************************************************************************
void AliMultDepSpecAnalysisTaskUE::FindLeadingTrack()
{
  fPtLead = 0;
  fPhiLead = 0;
  fSumPtMeasTot = 0;
  fSumPtTrueTot = 0;
  fMultMeasTot = 0;
  fMultTrueTot = 0;
  if(fIsMC)
  {
    fMCPtOfLead = 0;
    fMCPhiOfLead = 0;
  }

  AliVTrack* track = nullptr;
  for(int i = 0; i < fEvent->GetNumberOfTracks(); i++)
  {
    track = dynamic_cast<AliVTrack*>(fEvent->GetTrack(i));
    // Set fPt, fEta, fSigmapt; Check if track in kin range and has good quality

    if(!InitTrack(track)) continue;
    // initialize particle properties
    if(fIsMC)
    {
      // mc lable corresponding to measured track
      // negative lable indicates bad quality track - use it anyway
      int mcLable = TMath::Abs(track->GetLabel());

      // Set fMCPt, fMCEta, fMCIsChargedPrimary; Stop computation if it is not charged primary or
      // secondary
      if(fIsESD)
      {
        if(!InitParticle((AliMCParticle*)fMCEvent->GetTrack(mcLable))) continue;
      }
      else
      {
        if(!InitParticle((AliAODMCParticle*)fMCEvent->GetTrack(mcLable))) continue;
      }
    }
    // Find track with highest Pt
    if(fPt > fPtLead)
    {
      fPtLead  = fPt;
      fPhiLead = fPhi;
      if(fIsMC)
      {
        fMCPtOfLead  = fMCPt;
        fMCPhiOfLead = fMCPhi;
      }
    }

    fSumPtMeasTot += fPt;
    fMultMeasTot  += 1;

    if (fIsMC)
    {
      if (fMCIsChargedPrimary) fSumPtTrueTot += fMCPt;
      if (fMCIsChargedPrimary) fMultTrueTot  += 1;
    }

  }
}

//**************************************************************************************************
/**
 * Fill event histograms. Event related members are set.
 */
//**************************************************************************************************
void AliMultDepSpecAnalysisTaskUE::FillEventHistos()
{
  AliMultDepSpecAnalysisTaskOLD::FillEventHistos();
  if(fIsFirstEventInJob){
     fHistTrainInfoUE.Fill(Form("cos min: %.1f - cos max %.1f - ptleadcut %.2f", fCosPhiMin, fCosPhiMax, fPtLeadCut));
     fIsFirstEventInJob = false;
  }

  fHistSumPt.Fill(fSumPtMeasTot, fSumPtMeas);
  fHistMultCorr.Fill(fMultMeasTot, fMultMeas);
  fHistSumPtMult.Fill(fMultMeas, fSumPtMeas);

  fHistMeanSumPt.FillWeight(fSumPtMeas, fSumPtMeasTot);
  fHistMeanMult.FillWeight(fMultMeas, fMultMeasTot);
  fHistMeanSumPtMultTot.FillWeight(fSumPtMeas, fMultMeasTot);

  if(fMultMeas == 0) return; // ensure that number of tracks is non-zero
  fHistLeadPt.Fill(fPtLead);
  fHistLeadPhi.Fill(fPhiLead);
  fHistPlateau.Fill(fPtLead, fMultMeas);
  fHistLeadPtPhi.Fill(fPtLead, fPhiLead);

  if(fIsMC)
  {
    fHistMCSumPt.Fill(fSumPtTrueTot, fSumPtTrue);
    fHistMCMultCorr.Fill(fMultTrueTot, fMultTrue);
    fHistMCSumPtMult.Fill(fMultTrue, fSumPtTrue);

    fHistMCMeanSumPt.FillWeight(fSumPtTrue, fSumPtTrueTot);
    fHistMCMeanMult.FillWeight(fMultTrue, fMultTrueTot); //first weight then position
    fHistMCMeanSumPtMultTot.FillWeight(fSumPtTrue, fMultTrueTot);

    if (fMultTrue == 0) return;
    fHistMCPlateau.Fill(fMCPtOfLead, fMultTrue);
    // resolution of measurement:
    fHistMCResoPtLead.Fill(fabs(fPtLead - fMCPtOfLead));
    fHistMCResoPhiLead.Fill(fabs(fPhiLead - fMCPhiOfLead));
  }
}

//**************************************************************************************************
/**
 * Fill track histograms. Track related members are set.
 */
//**************************************************************************************************
void AliMultDepSpecAnalysisTaskUE::FillMeasTrackHistos()
{
  AliMultDepSpecAnalysisTaskOLD::FillMeasTrackHistos();

  if (fMultMeas == 0) return;
  if (fPhiLead-fPhi >= 0) fHistDeltaPhi.Fill(fPhiLead-fPhi);
  else fHistDeltaPhi.Fill(2*TMath::Pi() + fPhiLead-fPhi);
}

//**************************************************************************************************
/**
 * Fill measured particle histograms. Track and mc particle related members are set.
 */
//**************************************************************************************************
void AliMultDepSpecAnalysisTaskUE::FillMeasParticleHistos()
{
  AliMultDepSpecAnalysisTaskOLD::FillMeasParticleHistos();
}

//**************************************************************************************************
/**
 * Fill generated particle histograms. MC particle related members are set.
 */
//**************************************************************************************************
void AliMultDepSpecAnalysisTaskUE::FillTrueParticleHistos()
{
  AliMultDepSpecAnalysisTaskOLD::FillTrueParticleHistos();
}

//**************************************************************************************************
/**
 * Initializes track properties and returns false if track bad or out of range.
 */
//**************************************************************************************************
bool AliMultDepSpecAnalysisTaskUE::InitTrack(AliVTrack* track)
{
  bool isValidTrack = AliMultDepSpecAnalysisTaskOLD::InitTrack(track);
  // extract additional UE track infos
  // fDeltaLeadingPt = fLeadingPt - fPt;
  // fDeltaLeadingPhi = fLeadingPhi - fPhi;
  // fIsInTransverseRegion = ....

  return isValidTrack;
}

//**************************************************************************************************
/**
 * Initializes particle properties and returns false if out of range.
 * Different functions are called for AliMCParticles (ESD) and AliAODMCParticles (AOD).
 */
//**************************************************************************************************
bool AliMultDepSpecAnalysisTaskUE::InitParticle(AliMCParticle* particle)
{
  bool isValidParticle = AliMultDepSpecAnalysisTaskOLD::InitParticle(particle);

  return isValidParticle;
}

bool AliMultDepSpecAnalysisTaskUE::InitParticle(AliAODMCParticle* particle)
{
  bool isValidParticle = AliMultDepSpecAnalysisTaskOLD::InitParticle(particle);

  return isValidParticle;
}

//**************************************************************************************************
/**
 * User defined track selection criteria beyond good quality and fulfilling acceptance.that are
 * applied before counting multiplicities and filling histograms (after InitTrack() was called).
 */
//**************************************************************************************************
bool AliMultDepSpecAnalysisTaskUE::SelectTrack(bool count)
{

  return (fCosPhiMin <= TMath::Cos(fPhiLead - fPhi)) && (TMath::Cos(fPhiLead - fPhi) <= fCosPhiMax);

}

//**************************************************************************************************
/**
 * User defined particle selection criteria that are applied before counting multiplicities and
 * filling histograms (after InitParticle() was called).
 */
//**************************************************************************************************
bool AliMultDepSpecAnalysisTaskUE::SelectParticle(bool count)
{

   return ((fCosPhiMin <= TMath::Cos(fMCPhiOfLead - fMCPhi)) && (TMath::Cos(fMCPhiOfLead - fMCPhi) <= fCosPhiMax));

}

//**************************************************************************************************
/**
 * Function to hang an instance of this task in a LEGO train.
 */
//**************************************************************************************************
AliMultDepSpecAnalysisTaskUE* AliMultDepSpecAnalysisTaskUE::AddTaskMultDepSpecUE(
  const string& dataSet, int cutModeLow, int cutModeHigh, TString options, bool isMC,
  AliMultDepSpecAnalysisTaskUE::PhiRange phiRange, double ptLeadMIN,  const char* suffix)
{
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr)
  {
    ::Error("AddTaskMultDepSpecUE", "No analysis manager found.");
    return nullptr;
  }
  if(!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskMultDepSpecUE", "No input event handler found.");
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

  AliMultDepSpecAnalysisTaskUE* returnTask = nullptr;
  char taskName[100] = "";
  for(int cutMode = cutModeLow; cutMode <= cutModeHigh; cutMode++)
  {
    sprintf(taskName, "%s_%s_cutMode_%d_phiRange%d_%s", dataSet.data(), mode.data(), cutMode, phiRange, suffix);
    AliMultDepSpecAnalysisTaskUE* task = new AliMultDepSpecAnalysisTaskUE(taskName);
    if(!task->InitTask(phiRange, isMC, isAOD, dataSet, options, cutMode, ptLeadMIN))
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
 * Initialize task.
 */
//**************************************************************************************************
bool AliMultDepSpecAnalysisTaskUE::InitTask(AliMultDepSpecAnalysisTaskUE::PhiRange phiRange, bool isMC, bool isAOD, string dataSet,
                                            TString options, int cutMode, float ptLeadMIN)
{
  switch(phiRange){

    case full:
      SetPhiRange(-1., 1.);
      break;
    case toward:
      SetPhiRange(.5, 1.);
      break;
    case away:
      SetPhiRange(-1., -.5);
      break;
    case transverse:
      SetPhiRange(-.5, .5);
      break;
    std::cout << "----------------" << std::endl;
    std::cout << "cosphi:" << fCosPhiMin << " | " << fCosPhiMax << std::endl;
    std::cout << "----------------" << std::endl;

  }

  SetPtLeadCut(ptLeadMIN);
  return AliMultDepSpecAnalysisTaskOLD::InitTask(isMC, isAOD, dataSet, options, cutMode);
}
