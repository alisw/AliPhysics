#include <iostream>
#include "AlidNdPtTools.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskCutStudies.h"
#include "AliLog.h"
#include "AliESDtrack.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskCutStudies)
/// \endcond

//****************************************************************************************
/**
 * Default constructor.
 */
//****************************************************************************************
AliAnalysisTaskCutStudies::AliAnalysisTaskCutStudies()
: AliAnalysisTaskMKBase(), fHist_x{}, fHist_y{}, fHist_z{}, fHist_alpha{}, fHist_signed1Pt{}, fHist_snp{}, fHist_tgl{},
fHist_dcaxy{}, fHist_dcaz{}, fHist_flag{}, fHist_pt{}, fHist_eta{}, fHist_phi{}, fHist_zInner{}, fHist_itsFoundClusters{},
fHist_itsChi2PerCluster{}, fHist_itsHits{}, fHist_tpcFindableClusters{}, fHist_tpcFoundClusters{},
fHist_tpcSharedClusters{}, fHist_tpcFractionSharedClusters{}, fHist_tpcCrossedRows{},
fHist_tpcCrossedRowsOverFindableClusters{}, fHist_tpcChi2PerCluster{}
{
}



//****************************************************************************************
/**
 * Named constructor.
 */
//****************************************************************************************
AliAnalysisTaskCutStudies::AliAnalysisTaskCutStudies(const char* name)
: AliAnalysisTaskMKBase(name), fHist_x{}, fHist_y{}, fHist_z{}, fHist_alpha{}, fHist_signed1Pt{}, fHist_snp{}, fHist_tgl{},
fHist_dcaxy{}, fHist_dcaz{}, fHist_flag{}, fHist_pt{}, fHist_eta{}, fHist_phi{}, fHist_zInner{}, fHist_itsFoundClusters{},
fHist_itsChi2PerCluster{}, fHist_itsHits{}, fHist_tpcFindableClusters{}, fHist_tpcFoundClusters{},
fHist_tpcSharedClusters{}, fHist_tpcFractionSharedClusters{}, fHist_tpcCrossedRows{},
fHist_tpcCrossedRowsOverFindableClusters{}, fHist_tpcChi2PerCluster{}
{
}

//****************************************************************************************
/**
 * Destructor.
 */
//****************************************************************************************
AliAnalysisTaskCutStudies::~AliAnalysisTaskCutStudies()
{
}

//****************************************************************************************
/**
 * Add output to this task.
 */
//****************************************************************************************
void AliAnalysisTaskCutStudies::AddOutput()
{
  std::vector<double> centBins =
  {
    0., 30., 60., 90.
  };

  std::vector<double> ptBins =
  {
    0.0, 0.1, 0.2,  0.3,  0.4,  0.5,  0.6, 0.7,
    0.8, 0.9, 1.0,  1.1,  1.2,  1.3,  1.4, 1.5,
    2.0, 5.0, 10.0, 20.0, 50.0
  };
  const int nCuts = 2;
  Axis cutAxis =  {"cut", "cut setting", {-0.5, nCuts - 0.5}, nCuts};
  Axis centAxis = {"cent", "centrality", centBins};
  Axis ptAxis =   {"pt", "p_{T} [GeV/c]", ptBins};
  Axis etaAxis =  {"eta", "#eta", {-0.8, 0.8}, 2};
  Axis phiAxis =  {"phi", "#phi", {0., 2.*M_PI}, 4}; // 36 to see tpc sectors

  std::vector<Axis> defaultAxes = {cutAxis, centAxis, ptAxis, etaAxis, phiAxis};
  
  double requiredMemory;
  
  fHist_x.AddAxis("x", "x [cm]", 200, -0.4, 0.4);
  fOutputList->Add(fHist_x.GenerateHist("trackpar-x"));
  requiredMemory += fHist_x.GetSize();

  fHist_y.AddAxis("y", "y [cm]", 100, -4., 4.);
  fOutputList->Add(fHist_y.GenerateHist("trackpar-y"));
  requiredMemory += fHist_y.GetSize();

  fHist_z.AddAxis("z", "z [cm]", 100, -20., 20.);
  fOutputList->Add(fHist_z.GenerateHist("trackpar-z"));
  requiredMemory += fHist_z.GetSize();

  fHist_alpha.AddAxis("alpha", "#alpha [rad]", 100, -(M_PI + 0.01), (M_PI + 0.01));
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

  fHist_tpcFractionSharedClusters.AddAxis("tpcFractionSharedClusters", "fraction shared clusters TPC",  100, 0., 1.);
  fOutputList->Add(fHist_tpcFractionSharedClusters.GenerateHist("tpc-fractionSharedClusters"));
  requiredMemory += fHist_tpcFractionSharedClusters.GetSize();

  fHist_tpcFindableClusters.AddAxes(defaultAxes);
  fHist_tpcFindableClusters.AddAxis("findableClustersTPC", "# findable clusters TPC", 81, 79.5, 160.5);
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
  fHist_tpcCrossedRowsOverFindableClusters.AddAxis("crossedRowsOverFindableTPC", "crossed rows / findable clusters TPC",  50, 0.7, 1.2);
  fOutputList->Add(fHist_tpcCrossedRowsOverFindableClusters.GenerateHist("tpc-crossedRowsOverFindable"));
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

  
  AliError(Form("Estimated memory usage of histograms: %.0f Bytes (%f MiB)", requiredMemory, requiredMemory/1048576));
}

//****************************************************************************************
/**
 * Event selection.
 */
//****************************************************************************************
Bool_t AliAnalysisTaskCutStudies::IsEventSelected()
{
  return fIsAcceptedAliEventCuts;
}

//****************************************************************************************
/**
 * Analyse the event.
 */
//****************************************************************************************
void AliAnalysisTaskCutStudies::AnaEvent()
{
   LoopOverAllTracks();
   if (fIsMC) LoopOverAllParticles();
}

//****************************************************************************************
/**
 * Analyse the track.
 */
//****************************************************************************************
void AliAnalysisTaskCutStudies::AnaTrack(Int_t flag)
{
  if (!fAcceptTrackM) return;
  InitTrackQA();
  fZInner = (fESDTrack->GetInnerParam()) ? fESDTrack->GetInnerParam()->GetZ() : 999.;
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
  for(unsigned int i = 0; i < 64; i++) {
    if(fFlags & (1 << i))
      fHist_flag.Fill(i);
  }
  fHist_pt.Fill(fPt);
  fHist_eta.Fill(fEta);
  fHist_phi.Fill(fPhi);
  fHist_zInner.Fill(fZInner);

  // its related properties
  fHist_itsFoundClusters.Fill(fITSFoundClusters);
  for(unsigned int i = 0; i < 6; i++) {
    if(fITSClusterMap & (1 << i))
      fHist_itsHits.Fill(i);
  }

  // tpc related properties
  fHist_tpcSharedClusters.Fill(fTPCSharedClusters);
  fHist_tpcFractionSharedClusters.Fill(fTPCFractionSharedClusters);
  
  //cutID 0 == master cut
  for(int cutID = 0; cutID < 2; ++cutID)
  {
    if(fAcceptTrack[cutID])
    {
      fHist_tpcFindableClusters.Fill(cutID, fMultPercentileV0M, fPt, fEta, fPhi, fTPCFindableClusters);
      fHist_tpcCrossedRows.Fill     (cutID, fMultPercentileV0M, fPt, fEta, fPhi, fTPCCrossedRows);
      fHist_tpcCrossedRowsOverFindableClusters.Fill(cutID, fMultPercentileV0M, fPt, fEta, fPhi, fTPCCrossedRowsOverFindableClusters);
      fHist_tpcGoldenChi2.Fill      (cutID, fMultPercentileV0M, fPt, fEta, fPhi, fTPCGoldenChi2);

      fHist_tpcFoundClusters.Fill (cutID, fMultPercentileV0M, fPt, fEta, fPhi, fTPCFoundClusters);
      fHist_tpcChi2PerCluster.Fill(cutID, fMultPercentileV0M, fPt, fEta, fPhi, fTPCChi2PerCluster);
      fHist_tpcNClustersPID.Fill  (cutID, fMultPercentileV0M, fPt, fEta, fPhi, fTPCSignalN);
      fHist_tpcGeomLength.Fill    (cutID, fMultPercentileV0M, fPt, fEta, fPhi, fTPCGeomLength);
    }
  }
}

//****************************************************************************************
/**
 * Analyse the MC track.
 */
//****************************************************************************************
void AliAnalysisTaskCutStudies::AnaTrackMC(Int_t flag)
{
    if (!fAcceptTrackM) return;

}

//****************************************************************************************
/**
 * Analyse the MC particle.
 */
//****************************************************************************************
void AliAnalysisTaskCutStudies::AnaParticleMC(Int_t flag)
{            
    if (!fMCisPrim) return;    
    if (!fMCIsCharged) return;    
    if (TMath::Abs(fMCEta) > 0.8) return;    
  
}

//****************************************************************************************
/**
 * Add task of this kind to a train.
 */
//****************************************************************************************
AliAnalysisTaskCutStudies* AliAnalysisTaskCutStudies::AddTaskCutStudies(const char* name)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
      ::Error("AddTaskCutStudies", "No analysis manager to connect to.");
      return nullptr;
  }

  if (!mgr->GetInputEventHandler()) {
      ::Error("AddTaskCutStudies", "This task requires an input event handler.");
      return nullptr;
  }

  AliAnalysisTaskCutStudies *task = new AliAnalysisTaskCutStudies(name);
  if (!task) { return nullptr; }

  task->SelectCollisionCandidates(AliVEvent::kAnyINT);
  
  task->SetSkipMCtruth();
  task->SetNeedEventMult();

  AliESDtrackCuts* baseLineCut = AlidNdPtTools::CreateESDtrackCuts("defaultEta08");
  baseLineCut->SetMaxChi2PerClusterTPC(1e10);

  AliESDtrackCuts* defaultCut = AlidNdPtTools::CreateESDtrackCuts("defaultEta08");

  //AliESDtrackCuts* widerChi2Cut = AlidNdPtTools::CreateESDtrackCuts("defaultEta08");
  //widerChi2Cut->SetMaxChi2PerClusterTPC(6);

  //AliESDtrackCuts* tighterChi2Cut = AlidNdPtTools::CreateESDtrackCuts("defaultEta08");
  //tighterChi2Cut->SetMaxChi2PerClusterTPC(3);

  
  task->SetESDtrackCutsM(baseLineCut);
  task->SetESDtrackCuts(0, baseLineCut);
  task->SetESDtrackCuts(1, defaultCut);
  //task->SetESDtrackCuts(2, widerChi2Cut);
  //task->SetESDtrackCuts(3, tighterChi2Cut);

  
  
  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer(name, TList::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root"));

  return task;
}
