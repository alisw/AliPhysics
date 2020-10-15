#include <iostream>
#include "AlidNdPtTools.h"
#include "AliAnalysisManager.h"
#include "AliLog.h"
#include "AliESDtrack.h"
#include "AliAnalysisTaskCutStudies.h"

//**************************************************************************************************
/**
 * Default constructor.
 */
//**************************************************************************************************
AliAnalysisTaskCutStudies::AliAnalysisTaskCutStudies() : AliAnalysisTaskMKBase()
{
}

//**************************************************************************************************
/**
 * Named constructor.
 */
//**************************************************************************************************
AliAnalysisTaskCutStudies::AliAnalysisTaskCutStudies(const char* name) : AliAnalysisTaskMKBase(name)
{
}

//**************************************************************************************************
/**
 * Destructor.
 */
//**************************************************************************************************
AliAnalysisTaskCutStudies::~AliAnalysisTaskCutStudies()
{
}

//**************************************************************************************************
/**
 * Add output to this task.
 */
//**************************************************************************************************
void AliAnalysisTaskCutStudies::AddOutput()
{
  std::vector<double> centBins = { 0., 30., 60., 90. };

  std::vector<double> ptBins = { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7,  0.8,  0.9, 1.0,
                                 1.1, 1.2, 1.3, 1.4, 1.5, 2.0, 5.0, 10.0, 20.0, 50.0 };
  const int nCuts = 5;
  Hist::Axis cutAxis = { "cut", "cut setting", { -0.5, nCuts - 0.5 }, nCuts };
  Hist::Axis centAxis = { "cent", "centrality", centBins };
  Hist::Axis ptAxis = { "pt", "#it{p}_{T} (GeV/c)", ptBins };
  Hist::Axis etaAxis = { "eta", "#eta", { -0.8, 0.8 }, 2 };
  Hist::Axis phiAxis = { "phi", "#phi", { 0., 2. * M_PI }, 4 }; // 36 to see tpc sectors

  std::vector<Hist::Axis> defaultAxes = { cutAxis, centAxis, ptAxis, etaAxis, phiAxis };

  double requiredMemory = 0.;

  fHist_x.AddAxis("x", "x [cm]", 200, -0.4, 0.4);
  fOutputList->Add(fHist_x.GenerateHist("trackpar-x"));
  requiredMemory += fHist_x.GetSize();

  fHist_y.AddAxis("y", "y [cm]", 100, -4., 4.);
  fOutputList->Add(fHist_y.GenerateHist("trackpar-y"));
  requiredMemory += fHist_y.GetSize();

  fHist_z.AddAxis("z", "z [cm]", 100, -20., 20.);
  fOutputList->Add(fHist_z.GenerateHist("trackpar-z"));
  requiredMemory += fHist_z.GetSize();

  fHist_alpha.AddAxis("alpha", "#alpha (rad)", 100, -(M_PI + 0.01), (M_PI + 0.01));
  fOutputList->Add(fHist_alpha.GenerateHist("trackpar-alpha"));
  requiredMemory += fHist_alpha.GetSize();

  fHist_signed1Pt.AddAxis("signed1Pt", "q/p_{T}", 200, -8, 8);
  fOutputList->Add(fHist_signed1Pt.GenerateHist("trackpar-signed1Pt"));
  requiredMemory += fHist_signed1Pt.GetSize();

  fHist_snp.AddAxis("snp", "snp", 100, -1., 1.);
  fOutputList->Add(fHist_snp.GenerateHist("trackpar-snp"));
  requiredMemory += fHist_snp.GetSize();

  fHist_tgl.AddAxis("tgl", "tgl", 1000, -2, 2);
  fOutputList->Add(fHist_tgl.GenerateHist("trackpar-tgl"));
  requiredMemory += fHist_tgl.GetSize();

  fHist_dcaxy.AddAxis("dcaxy", "dca xy", 200, -3., 3.);
  fOutputList->Add(fHist_dcaxy.GenerateHist("trackpar-dcaXY"));
  requiredMemory += fHist_dcaxy.GetSize();

  fHist_dcaz.AddAxis("dcaz", "dca z", 200, -3., 3.);
  fOutputList->Add(fHist_dcaz.GenerateHist("trackpar-dcaZ"));
  requiredMemory += fHist_dcaz.GetSize();

  fHist_zInner.AddAxis("zInner", "z at inner param", 400, -85., 85.);
  fOutputList->Add(fHist_zInner.GenerateHist("trackpar-zInner"));
  requiredMemory += fHist_zInner.GetSize();

  fHist_flag.AddAxis("flag", "track flag", 64, -0.5, 63.5);
  fOutputList->Add(fHist_flag.GenerateHist("track-flags"));
  requiredMemory += fHist_flag.GetSize();

  fHist_pt.AddAxis("pt", "p_{T} [GeV/c]", 100, 0., 50.);
  fOutputList->Add(fHist_pt.GenerateHist("track-pt"));
  requiredMemory += fHist_pt.GetSize();

  fHist_eta.AddAxis("eta", "#eta", 101, -1.0, 1.0);
  fOutputList->Add(fHist_eta.GenerateHist("track-eta"));
  requiredMemory += fHist_eta.GetSize();

  fHist_phi.AddAxis("phi", "#phi [rad]", 100, 0., 2 * M_PI);
  fOutputList->Add(fHist_phi.GenerateHist("track-phi"));
  requiredMemory += fHist_phi.GetSize();

  fHist_itsFoundClusters.AddAxis("itsFoundClusters", "# clusters ITS", 8, -0.5, 7.5);
  fOutputList->Add(fHist_itsFoundClusters.GenerateHist("its-foundClusters"));
  requiredMemory += fHist_itsFoundClusters.GetSize();

  fHist_itsHits.AddAxis("itsHits", "layer ITS", 7, -0.5, 6.5);
  fOutputList->Add(fHist_itsHits.GenerateHist("its-hits"));
  requiredMemory += fHist_itsHits.GetSize();

  // TODO: add CHI2 / hit ITS (500, 0, 100) and number of hits SPD

  fHist_tpcSharedClusters.AddAxis("sharedClustersTPC", "# shared clusters TPC", 161, -0.5, 160.5);
  fOutputList->Add(fHist_tpcSharedClusters.GenerateHist("tpc-sharedClusters"));
  requiredMemory += fHist_tpcSharedClusters.GetSize();

  fHist_tpcFractionSharedClusters.AddAxis("tpcFractionSharedClusters",
                                          "fraction shared clusters TPC", 100, 0., 1.);
  fOutputList->Add(fHist_tpcFractionSharedClusters.GenerateHist("tpc-fractionSharedClusters"));
  requiredMemory += fHist_tpcFractionSharedClusters.GetSize();

  fHist_tpcFindableClusters.AddAxes(defaultAxes);
  fHist_tpcFindableClusters.AddAxis("findableClustersTPC", "# findable clusters TPC", 81, 79.5,
                                    160.5);
  fOutputList->Add(fHist_tpcFindableClusters.GenerateHist("tpc-findableClusters"));
  requiredMemory += fHist_tpcFindableClusters.GetSize();

  fHist_tpcFoundClusters.AddAxes(defaultAxes);
  fHist_tpcFoundClusters.AddAxis("foundClustersTPC", "# found clusters TPC", 81, 79.5, 160.5);
  fOutputList->Add(fHist_tpcFoundClusters.GenerateHist("tpc-foundClusters"));
  requiredMemory += fHist_tpcFoundClusters.GetSize();

  fHist_tpcCrossedRows.AddAxes(defaultAxes);
  fHist_tpcCrossedRows.AddAxis("crossedRowsTPC", "# crossed rows TPC", 81, 79.5, 160.5);
  fOutputList->Add(fHist_tpcCrossedRows.GenerateHist("tpc-crossedRows"));
  requiredMemory += fHist_tpcCrossedRows.GetSize();

  fHist_tpcCrossedRowsOverFindableClusters.AddAxes(defaultAxes);
  fHist_tpcCrossedRowsOverFindableClusters.AddAxis(
    "crossedRowsOverFindableTPC", "crossed rows / findable clusters TPC", 50, 0.7, 1.2);
  fOutputList->Add(
    fHist_tpcCrossedRowsOverFindableClusters.GenerateHist("tpc-crossedRowsOverFindable"));
  requiredMemory += fHist_tpcCrossedRowsOverFindableClusters.GetSize();

  fHist_tpcChi2PerCluster.AddAxes(defaultAxes);
  fHist_tpcChi2PerCluster.AddAxis("chi2PerClusterTPC", "chi2 / cluster TPC", 140, 0., 7.);
  fOutputList->Add(fHist_tpcChi2PerCluster.GenerateHist("tpc-chi2PerCluster"));
  requiredMemory += fHist_tpcChi2PerCluster.GetSize();

  fHist_tpcNClustersPID.AddAxes(defaultAxes);
  fHist_tpcNClustersPID.AddAxis("nClustersPIDTPC", "# clusters PID TPC", 91, 69.5, 160.5);
  fOutputList->Add(fHist_tpcNClustersPID.GenerateHist("tpc-nClustersPID"));
  requiredMemory += fHist_tpcNClustersPID.GetSize();

  fHist_tpcGoldenChi2.AddAxes(defaultAxes);
  fHist_tpcGoldenChi2.AddAxis("goldenChi2TPC", "chi2 global vs. TPC constrained", 41, -0.5, 40.5);
  fOutputList->Add(fHist_tpcGoldenChi2.GenerateHist("tpc-goldenChi2"));
  requiredMemory += fHist_tpcGoldenChi2.GetSize();

  fHist_tpcGeomLength.AddAxes(defaultAxes);
  fHist_tpcGeomLength.AddAxis("geomLengthTPC", "geometric length in TPC", 51, 111.5, 162.5);
  fOutputList->Add(fHist_tpcGeomLength.GenerateHist("tpc-geomLength"));
  requiredMemory += fHist_tpcGeomLength.GetSize();

  // tmp crosscheck histogram
  fHist_correlChi2GeomLength.AddAxis("geomLengthTPC", "geometric length in TPC", 51, 111.5, 162.5);
  fHist_correlChi2GeomLength.AddAxis("chi2PerClusterTPC", "chi2 / cluster TPC", 140, 0., 7.);
  fOutputList->Add(fHist_correlChi2GeomLength.GenerateHist("tpc-correlChi2GeomLength"));
  requiredMemory += fHist_correlChi2GeomLength.GetSize();

  AliError(Form("Estimated memory usage of histograms: %.2f MiB.", requiredMemory / 1048576));
}

//**************************************************************************************************
/**
 * Event selection.
 */
//**************************************************************************************************
Bool_t AliAnalysisTaskCutStudies::IsEventSelected()
{
  return fIsAcceptedAliEventCuts;
}

//**************************************************************************************************
/**
 * Analyse the event.
 */
//**************************************************************************************************
void AliAnalysisTaskCutStudies::AnaEvent()
{
  LoopOverAllTracks();
  if(fIsMC) LoopOverAllParticles();
}

//**************************************************************************************************
/**
 * Analyse the track.
 */
//**************************************************************************************************
void AliAnalysisTaskCutStudies::AnaTrack(Int_t flag)
{
  const int nCuts = 5;
  // if none of the track cuts lets this one pass, just skip the following calculations
  bool interestingTrack = false;
  for(int i = 0; i < nCuts; i++)
    interestingTrack = interestingTrack || fAcceptTrack[i];
  if(!interestingTrack) return;

  InitTrackQA();
  fZInner = (fESDTrack->GetInnerParam()) ? fESDTrack->GetInnerParam()->GetZ() : 999.;

  // fill some basic quantities only if they fulfil the strictest cut
  if(fAcceptTrack[0])
  {
    // track related properties
    fHist_x.Fill(fX);
    fHist_y.Fill(fY);
    fHist_z.Fill(fZ);
    fHist_alpha.Fill(fAlpha);
    fHist_signed1Pt.Fill(fSigned1Pt);
    fHist_snp.Fill(fSnp);
    fHist_tgl.Fill(fTgl);
    fHist_dcaxy.Fill(fDCAr);
    fHist_dcaz.Fill(fDCAz);
    for(unsigned int i = 0; i < 64; i++)
    {
      if(fFlags & (1 << i)) fHist_flag.Fill(i);
    }
    fHist_pt.Fill(fPt);
    fHist_eta.Fill(fEta);
    fHist_phi.Fill(fPhi);
    fHist_zInner.Fill(fZInner);

    // its related properties
    fHist_itsFoundClusters.Fill(fITSFoundClusters);
    for(unsigned int i = 0; i < 6; i++)
    {
      if(fITSClusterMap & (1 << i)) fHist_itsHits.Fill(i);
    }

    // tpc related properties
    fHist_tpcSharedClusters.Fill(fTPCSharedClusters);
    fHist_tpcFractionSharedClusters.Fill(fTPCFractionSharedClusters);

    fHist_correlChi2GeomLength.Fill(fTPCGeomLength, fTPCChi2PerCluster);
  }

  for(int cutID = 0; cutID < nCuts; ++cutID)
  {
    if(fAcceptTrack[cutID])
    {
      if(cutID == 2 && fTPCSignalN < 70)
        continue; // in cut mode 2 apply additional requirement on number of clusters for pid
      fHist_tpcFindableClusters.Fill(cutID, fMultPercentileV0M, fPt, fEta, fPhi,
                                     fTPCFindableClusters);
      fHist_tpcCrossedRows.Fill(cutID, fMultPercentileV0M, fPt, fEta, fPhi, fTPCCrossedRows);
      fHist_tpcCrossedRowsOverFindableClusters.Fill(cutID, fMultPercentileV0M, fPt, fEta, fPhi,
                                                    fTPCCrossedRowsOverFindableClusters);
      fHist_tpcGoldenChi2.Fill(cutID, fMultPercentileV0M, fPt, fEta, fPhi, fTPCGoldenChi2);

      fHist_tpcFoundClusters.Fill(cutID, fMultPercentileV0M, fPt, fEta, fPhi, fTPCFoundClusters);
      fHist_tpcChi2PerCluster.Fill(cutID, fMultPercentileV0M, fPt, fEta, fPhi, fTPCChi2PerCluster);
      fHist_tpcNClustersPID.Fill(cutID, fMultPercentileV0M, fPt, fEta, fPhi, fTPCSignalN);
      fHist_tpcGeomLength.Fill(cutID, fMultPercentileV0M, fPt, fEta, fPhi, fTPCGeomLength);
    }
  }
}

//**************************************************************************************************
/**
 * Analyse the MC track.
 */
//**************************************************************************************************
void AliAnalysisTaskCutStudies::AnaTrackMC(Int_t flag)
{
  // if (!fAcceptTrackM) return;
}

//**************************************************************************************************
/**
 * Analyse the MC particle.
 */
//**************************************************************************************************
void AliAnalysisTaskCutStudies::AnaParticleMC(Int_t flag)
{
  if(!fMCisPrim) return;
  if(!fMCIsCharged) return;
  if(TMath::Abs(fMCEta) > 0.8) return;
}

//**************************************************************************************************
/**
 * Add task of this kind to a train.
 */
//**************************************************************************************************
AliAnalysisTaskCutStudies* AliAnalysisTaskCutStudies::AddTaskCutStudies(const char* name)
{
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr)
  {
    ::Error("AddTaskCutStudies", "No analysis manager to connect to.");
    return nullptr;
  }

  if(!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskCutStudies", "This task requires an input event handler.");
    return nullptr;
  }

  AliAnalysisTaskCutStudies* task = new AliAnalysisTaskCutStudies(name);
  if(!task)
  {
    return nullptr;
  }

  task->SelectCollisionCandidates(AliVEvent::kAnyINT);

  task->SetSkipMCtruth();
  task->SetNeedEventMult();

  task->SetESDtrackCuts(0, GetCutSetting("chargedParticles"));
  task->SetESDtrackCuts(1, GetCutSetting("global"));
  task->SetESDtrackCuts(2, GetCutSetting("global-nodca-pidClCut"));
  task->SetESDtrackCuts(3, GetCutSetting("global-nodca"));
  task->SetESDtrackCuts(4, GetCutSetting("tpc-only"));

  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(
    task, 1,
    mgr->CreateContainer(name, TList::Class(), AliAnalysisManager::kOutputContainer,
                         "AnalysisResults.root"));

  return task;
}

// temporary function that handles cut setting creation
AliESDtrackCuts* AliAnalysisTaskCutStudies::GetCutSetting(const std::string& identifier)
{
  AliESDtrackCuts* trackCuts;

  if(identifier == "chargedParticles")
  {
    trackCuts = AlidNdPtTools::CreateESDtrackCuts("defaultEta08");
  }
  else if(identifier == "tpc-only") // Filter Bit 0
  {
    trackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
  }
  else if(identifier == "global-nodca") // Filter Bit 4 (standard cuts with very loose DCA)
  {
    trackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE);
    trackCuts->SetMaxDCAToVertexXY(2.4);
    trackCuts->SetMaxDCAToVertexZ(3.2);
    trackCuts->SetDCAToVertex2D(kTRUE);
  }
  else if(identifier == "global-nodca-pidClCut") // Filter Bit 4 with additional pid cluster cut <70
                                                 // (applied in AnaTrack)
  {
    trackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE);
    trackCuts->SetMaxDCAToVertexXY(2.4);
    trackCuts->SetMaxDCAToVertexZ(3.2);
    trackCuts->SetDCAToVertex2D(kTRUE);
  }
  else if(identifier == "global") // Filter Bit 5 (standard cuts with tight DCA cut)
  {
    trackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();
  }
  // always restrict eta range
  trackCuts->SetEtaRange(-0.8, 0.8);
  // open chi2 for all cut settings
  trackCuts->SetMaxChi2PerClusterTPC(1e10);

  return trackCuts;
}
