/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* AliAnalysisTaskAO2Dconverter
 *
 * Convert Run 2 ESDs to Run 3 AODs (AO2D.root).
 */

#include <TFile.h>
#include <TDirectory.h>
#include <TChain.h>
#include <TTree.h>
#include <TMath.h>
#include <TTimeStamp.h>
#include <TSystem.h>
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliEMCALGeometry.h"
#include "AliAnalysisTaskAO2Dconverter.h"
#include "AliVHeader.h"
#include "COMMON/MULTIPLICITY/AliMultSelection.h"

#include "AliESDCaloCells.h"
#include "AliESDCaloTrigger.h"
#include "AliESDHeader.h"
#include "AliESDtrack.h"
#include "AliESDMuonTrack.h"
#include "AliESDMuonCluster.h"
#include "AliESDZDC.h"
#include "AliESDVZERO.h"
#include "AliESDv0.h"
#include "AliESDcascade.h"

#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliPIDResponse.h"

#include "AliGenCocktailEventHeader.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliGenEpos3EventHeader.h"
#include "AliGenEposEventHeader.h"
#include "AliGenEventHeader.h"
#include "AliGenEventHeaderTunedPbPb.h"
#include "AliGenGeVSimEventHeader.h"
#include "AliGenHepMCEventHeader.h"
#include "AliGenHerwigEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenToyEventHeader.h"

#include "AliMathBase.h"
#include "AliLog.h"

ClassImp(AliAnalysisTaskAO2Dconverter);

const TString AliAnalysisTaskAO2Dconverter::TreeName[kTrees] = { "O2collision", "DbgEventExtra", "O2track", "O2calo",  "O2calotrigger", "O2muon", "O2muoncluster", "O2zdc", "O2fv0a", "O2fv0c", "O2ft0", "O2fdd", "O2v0", "O2cascade", "O2tof", "O2mcparticle", "O2mccollision", "O2mctracklabel", "O2mccalolabel", "O2mccollisionlabel", "O2bc" };

const TString AliAnalysisTaskAO2Dconverter::TreeTitle[kTrees] = { "Collision tree", "Collision extra", "Barrel tracks", "Calorimeter cells", "Calorimeter triggers", "MUON tracks", "MUON clusters", "ZDC", "FV0A", "FV0C", "FT0", "FDD", "V0s", "Cascades", "TOF hits", "Kinematics", "MC collisions", "MC track labels", "MC calo labels", "MC collision labels", "BC info" };

const TClass* AliAnalysisTaskAO2Dconverter::Generator[kGenerators] = { AliGenEventHeader::Class(), AliGenCocktailEventHeader::Class(), AliGenDPMjetEventHeader::Class(), AliGenEpos3EventHeader::Class(), AliGenEposEventHeader::Class(), AliGenEventHeaderTunedPbPb::Class(), AliGenGeVSimEventHeader::Class(), AliGenHepMCEventHeader::Class(), AliGenHerwigEventHeader::Class(), AliGenHijingEventHeader::Class(), AliGenPythiaEventHeader::Class(), AliGenToyEventHeader::Class() };

namespace
{
  // Helper function
  ULong64_t GetGlobalBC(AliVHeader *header)
  {
    return ((ULong64_t)header->GetBunchCrossNumber() +
	    (ULong64_t)header->GetOrbitNumber() * 3564 +
	    (ULong64_t)header->GetPeriodNumber() * 16777216 * 3564);
  }

  // Initialize the precision masks used to truncate the corresponding float data members
  // By default no truncation
  
  UInt_t mCollisionPosition = 0xFFFFFFFF;    // Position in x,y,z
  UInt_t mCollisionPositionCov = 0xFFFFFFFF; // Covariance matrix and chi2

  UInt_t mTrackX =  0xFFFFFFFF;
  UInt_t mTrackAlpha = 0xFFFFFFFF;
  UInt_t mtrackSnp = 0xFFFFFFFF;
  UInt_t mTrackTgl = 0xFFFFFFFF;
  UInt_t mTrack1Pt = 0xFFFFFFFF; // Including the momentun at the inner wall of TPC
  UInt_t mTrackCovDiag = 0xFFFFFFFF; // Including the chi2
  UInt_t mTrackCovOffDiag = 0xFFFFFFFF;
  UInt_t mTrackSignal = 0xFFFFFFFF; // PID signals and track length
  UInt_t mTrackPosEMCAL = 0xFFFFFFFF;

  UInt_t mTracklets = 0xFFFFFFFF; // tracklet members

  UInt_t mMcParticleW   = 0xFFFFFFFF; // Precision for weight
  UInt_t mMcParticlePos = 0xFFFFFFFF; // Precision for (x,y,z,t)
  UInt_t mMcParticleMom = 0xFFFFFFFF; // Precision for (Px,Py,Pz,E)

  UInt_t mCaloAmp = 0xFFFFFFFF;
  UInt_t mCaloTime = 0xFFFFFFFF;

  UInt_t mMuonTr1P = 0xFFFFFFFF;
  UInt_t mMuonTrThetaX = 0xFFFFFFFF;
  UInt_t mMuonTrThetaY = 0xFFFFFFFF;
  UInt_t mMuonTrZmu = 0xFFFFFFFF;
  UInt_t mMuonTrBend = 0xFFFFFFFF;
  UInt_t mMuonTrNonBend = 0xFFFFFFFF;
  UInt_t mMuonTrCov = 0xFFFFFFFF; // Covariance matrix and chi2

  UInt_t mMuonCl = 0xFFFFFFFF; // Position and charge
  UInt_t mMuonClErr = 0xFFFFFFFF;

  UInt_t mV0Time = 0xFFFFFFFF;
  UInt_t mADTime = 0xFFFFFFFF;
  UInt_t mT0Time = 0xFFFFFFFF;
  UInt_t mV0Amplitude = 0xFFFFFFFF;
  UInt_t mADAmplitude = 0xFFFFFFFF;
  UInt_t mT0Amplitude = 0xFFFFFFFF;
  
  // No compression for ZDC for the moment

} // namespace

AliAnalysisTaskAO2Dconverter::AliAnalysisTaskAO2Dconverter(const char* name)
    : AliAnalysisTaskSE(name)
    , fTrackFilter(Form("AO2Dconverter%s", name), Form("fTrackFilter%s", name))
    , fEventCuts{}
    , collision()
    , eventextra()
    , bc()
    , tracks()
    , mccollision()
    , mctracklabel()
    , mccalolabel()
    , mccollisionlabel()
    , mcparticle()
#ifdef USE_TOF_CLUST
    , tofClusters()
#endif
    , calo()
    , calotrigger()
    , muons()
    , mucls()
    , zdc()
    , fv0a()
    , fv0c()
    , ft0()
    , fdd()
    , v0s()
    , cascs()
{
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  for (Int_t i = 0; i < kTrees; i++) {
    fTreeStatus[i] = kTRUE;
  }
} // AliAnalysisTaskAO2Dconverter::AliAnalysisTaskAO2Dconverter(const char* name)

AliAnalysisTaskAO2Dconverter::~AliAnalysisTaskAO2Dconverter()
{
  fOutputList->Delete();
  delete fOutputList;
} // AliAnalysisTaskAO2Dconverter::~AliAnalysisTaskAO2Dconverter()

void AliAnalysisTaskAO2Dconverter::UserCreateOutputObjects()
{
  // Setting active/inactive containers based on the TaskMode
  switch (fTaskMode) {
  case kStandard:
    DisableTree(kMcParticle);
    DisableTree(kMcCollision);
    DisableTree(kMcTrackLabel);
    DisableTree(kMcCaloLabel);
    DisableTree(kMcCollisionLabel);
    break;
  default:
    break;
  }

  // Set the truncation
  if (fTruncate) {
    mCollisionPosition = 0xFFFFFFF0; // 19 bits mantissa
    mCollisionPositionCov = 0xFFFFE000; // 10 bits mantissa

    mTrackX =  0xFFFFFFF0; // 19 bits
    mTrackAlpha = 0xFFFFFFF0; // 19 bits
    mtrackSnp = 0xFFFFFF00; // 15 bits
    mTrackTgl = 0xFFFFFF00; // 15 bits
    mTrack1Pt = 0xFFFFFC00; // 13 bits
    mTrackCovDiag = 0xFFFFFF00; // 15 bits
    mTrackCovOffDiag = 0xFFFF0000; // 7 bits
    mTrackSignal = 0xFFFFFF00; // 15 bits
    mTrackPosEMCAL = 0xFFFFFF00; // 15 bits;

    mTracklets = 0xFFFFFF00; // 15 bits

    mMcParticleW   = 0xFFFFFFF0; // 19 bits
    mMcParticlePos = 0xFFFFFFF0; // 19 bits
    mMcParticleMom = 0xFFFFFFF0; // 19 bits

    mCaloAmp = 0xFFFFFF00; // 15 bits
    mCaloTime = 0xFFFFFF00; // 15 bits

    mMuonTr1P = 0xFFFFFC00; // 13 bits 
    mMuonTrThetaX = 0xFFFFFF00; // 15 bits
    mMuonTrThetaY = 0xFFFFFF00; // 15 bits
    mMuonTrZmu = 0xFFFFFFF0; // 19 bits
    mMuonTrBend = 0xFFFFFFF0; // 19 bits
    mMuonTrNonBend = 0xFFFFFFF0; // 19 bits
    mMuonTrCov = 0xFFFF0000; // 7 bits

    mMuonCl = 0xFFFFFF00; // 15 bits
    mMuonClErr = 0xFFFF0000; // 7 bits
    
    mV0Time = 0xFFFFF000; // 11 bits
    mADTime = 0xFFFFF000; // 11 bits
    mT0Time = 0xFFFFFF00; // 15 bits
    mV0Amplitude = 0xFFFFF000; // 11 bits
    mADAmplitude = 0xFFFFF000; // 11 bits
    mT0Amplitude = 0xFFFFF000; // 11 bits
  }
  
  // create output objects
  OpenFile(1); // Here we have the histograms
  /// Option compress is used to specify the compression level and algorithm:
  ///
  ///     compress = 100 * algorithm + level
  ///
  /// Level | Explanation
  /// ------|-------------
  /// 0   | objects written to this file will not be compressed.
  /// 1   | minimal compression level but fast.
  /// ... | ....
  /// 9   | maximal compression level but slower and might use more memory.
  /// algorithm = 1 : ZLIB compression algorithm is used (default)
  /// algorithm = 2 : LZMA compression algorithm is used
  /// algorithm = 4 : LZ4  compression algorithm is used
  /// algorithm = 5 : ZSTD compression algorithm is used
  /// So fCompress = 409 is LZ4 algorithm level 9


  fOutputFile = TFile::Open("AO2D.root","RECREATE", "O2 AOD", fCompress); // File to store the trees of time frames
  fOutputFile->Print();

  // create the list of output histograms
  fOutputList = new TList();
  fOutputList->SetOwner();

  // Add centrality histogram
  fCentralityHist = new TH1F("centrality", TString::Format("Centrality %s", fCentralityMethod.Data()),
                             100, 0.0, 100.0);
  fCentralityINT7 = new TH1F("centralityINT7", TString::Format("Centrality %s INT7", fCentralityMethod.Data()),
                             100, 0.0, 100.0);
  fHistPileupEvents = new TH1I("puEvents", "Pileup events", 2, 0, 2);
  fHistPileupEvents->SetStats(0);

  fOutputList->Add(fCentralityHist);
  fOutputList->Add(fCentralityINT7);
  fOutputList->Add(fHistPileupEvents);
  if (fSkipTPCPileup || fSkipPileup || fUseEventCuts) fEventCuts.AddQAplotsToList(fOutputList);
  if (fSkipTPCPileup) fEventCuts.SetRejectTPCPileupWithITSTPCnCluCorr(true);

  PostData(1, fOutputList);
} // void AliAnalysisTaskAO2Dconverter::UserCreateOutputObjects()

void AliAnalysisTaskAO2Dconverter::UserExec(Option_t *)
{
  // Initialisation

  const char *kPileupRejType[2] = {"PU_rej", "PU_TPC_rej"};

  fESD = dynamic_cast<AliESDEvent *>(InputEvent());
  if (!fESD) {
    ::Fatal("AliAnalysisTaskAO2Dconverter::UserExec", "Something is wrong with the event handler");
  }

  // We can use event cuts to avoid cases where we have zero reconstructed tracks
  bool skip_event = false;
  if (fUseEventCuts || fSkipPileup || fSkipTPCPileup) {
    skip_event = !fEventCuts.AcceptEvent(fESD) && fUseEventCuts;
  }

  // Skip pileup events if requested
  if (fSkipPileup && !fEventCuts.PassedCut(AliEventCuts::kPileUp)) {
    fHistPileupEvents->Fill(kPileupRejType[0], 1);
    skip_event = true;
  }

  // Skip TPC pileup events if requested, but don't skip event since it may affect physics
  if (fSkipTPCPileup && !fEventCuts.PassedCut(AliEventCuts::kTPCPileUp)) {
    fHistPileupEvents->Fill(kPileupRejType[1], 1);
    //skip_event = true;
  }
  
  if (fTaskMode == kStandard)
    if (fESD->GetHeader()->GetEventType() != 7) // check for PHYSICS events
      skip_event = true;

  if (skip_event) {
    return;
  }
  
  if (!fTfInitialized) {
    ULong64_t tfId = GetGlobalBC(fESD->GetHeader());
    if (tfId == 0) {
      // The time stamp of the event is not set, for example in MC
      // Try the AliEn job ID
      TString alienPID(gSystem->Getenv("ALIEN_PROC_ID"));
      if (alienPID.Length() > 0) {
        tfId = alienPID.Atoll() * 1000 + fTFCount;
      } else {
        // Fallback:
        // Use the time period from 2020/11/01 in seconds (similar to what we did in the DAQ LDCs)
        TTimeStamp ts0(2020,11,1,0,0,0);
        TTimeStamp ts1;
        tfId = ts1.GetSec() - ts0.GetSec();
        tfId *= 1000 + fTFCount;
      }
    }
    // Make globally unique TF id (across all data-taking)
    // The longest fill in Run 2 was 38 hours, which needs 43 bits. We reserve the values up to 1e13 which corresponds to 69 hours.
    // Run numbers in Run 2 were < 300k, which needs 19 bits
    // To make the number human-readable, we avoid a bit shift, but use a multiplication
    if (fESD->GetRunNumber() > 0)
      tfId += (ULong64_t) fESD->GetRunNumber() * 10000000000000L;
    
    InitTF(tfId);
  }
  // Get multiplicity selection
  AliMultSelection *multSelection = (AliMultSelection*) fESD->FindListObject("MultSelection");
  if (!multSelection)
    AliFatal("MultSelection not found in input event");

  float centrality = multSelection->GetMultiplicityPercentile(fCentralityMethod);

  // Selection of events with at least two contributors (GMI)
  // Can this be done using the physics selection? (PH)

  const AliESDVertex * pvtx = fESD->GetPrimaryVertex();
  if (!pvtx) {
    ::Fatal("AliAnalysisTaskAO2Dconverter::UserExec", "Vertex not defined");
  }
  TString title=pvtx->GetTitle();
  
  // bypass vertex selection for muon triggers
  Bool_t applyVertexSelection = kTRUE;
  if (fESD->GetFiredTriggerClasses().Contains("-MUON")) applyVertexSelection = kFALSE; // MUON cluster
  if (fESD->GetFiredTriggerClasses().Contains("-MUFAST")) applyVertexSelection = kFALSE; // MUFAST cluster
  if (fESD->GetFiredTriggerClasses().Contains("CMUP")) applyVertexSelection = kFALSE; // MUON UPC including semiforward

  if (applyVertexSelection) {
    if(pvtx->IsFromVertexer3D() || pvtx->IsFromVertexerZ()) return;
    if(pvtx->GetNContributors()<2) return;
  }

  // Fill centrality QA plots
  fCentralityHist->Fill(centrality);
  if ((fInputHandler->IsEventSelected() & AliVEvent::kINT7) != 0)
    fCentralityINT7->Fill(centrality);

  // Now fill the content of the TF
  FillEventInTF();

  // Finish the current TF and initialize a new one, if the size is above the limit
  if (fBytes > fMaxBytes) {
    AliInfo(Form("Total size of output trees: %lu bytes\n", fBytes));
    fBytes = 0; // Reset the byte counter
    fTfInitialized = false;
    FinishTF();
  }

  //---------------------------------------------------------------------------
  //Posting data
  PostData(1, fOutputList);
} // void AliAnalysisTaskAO2Dconverter::UserExec(Option_t *)

void AliAnalysisTaskAO2Dconverter::FinishTaskOutput()
{
  // called at the end of the event loop on the worker
  FinishTF();
  fOutputFile->Write(); // Do not close the file since this is then re-opened and overwritten by the framework
  AliInfo(Form("Total size of output trees: %lu bytes\n", fBytes));
}

void AliAnalysisTaskAO2Dconverter::Terminate(Option_t *)
{
  // called at the END of the analysis AFTER merging. In grid this is NOT called on the workers
}

AliAnalysisTaskAO2Dconverter *AliAnalysisTaskAO2Dconverter::AddTask(TString suffix)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    return nullptr;
  }
  // get the input event handler, again via a static method.
  // this handler is part of the managing system and feeds events
  // to your task
  if (!mgr->GetInputEventHandler())
  {
    return nullptr;
  }
  // by default, a file is open for writing. here, we get the filename
  TString fileName = "qa.root";
  if (!suffix.IsNull())
    fileName += ":" + suffix; // create a subfolder in the file
  // now we create an instance of your task
  AliAnalysisTaskAO2Dconverter *task = new AliAnalysisTaskAO2Dconverter((TString("AO2Dfriend") + suffix).Data());
  if (!task)
    return nullptr;
  // add your task to the manager
  mgr->AddTask(task);
  // your task needs input: here we connect the manager to your task
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  // same for the output
  mgr->ConnectOutput(task, 1, mgr->CreateContainer("QAlist", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  // for (Int_t i = 0; i < kTrees; i++)
  //   mgr->ConnectOutput(task, 2 + i, mgr->CreateContainer(TreeName[i], TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  // in the end, this macro returns a pointer to your task. this will be convenient later on
  // when you will run your analysis in an analysis train on grid
  return task;
} // AliAnalysisTaskAO2Dconverter *AliAnalysisTaskAO2Dconverter::AddTask(TString suffix)

////////////////////////////////////////////////////////////
TTree* AliAnalysisTaskAO2Dconverter::CreateTree(TreeIndex t)
{
  if (!fTreeStatus[t]) return 0x0;
  // Create the tree in the corresponding (TF) directory
  if (!fOutputDir) AliFatal("No Root subdir|");
  fOutputDir->cd();
  AliInfo(Form("Creating tree %s\n", TreeName[t].Data()));
  fTree[t] = new TTree(TreeName[t], TreeTitle[t]);
  fTree[t]->SetAutoFlush(0);
  return fTree[t];
} // TTree* AliAnalysisTaskAO2Dconverter::CreateTree(TreeIndex t)

void AliAnalysisTaskAO2Dconverter::Prune()
{
  if (fPruneList.IsNull() || fPruneList.IsWhitespace())
    return;
  TObjArray* arr = fPruneList.Tokenize(" ");
  for (Int_t i = 0; i < arr->GetEntries(); i++) {
    Bool_t found = kFALSE;
    for (Int_t j = 0; j < kTrees; j++) {
      TObjArray* branches = fTree[j]->GetListOfBranches();
      for (Int_t k = 0; k < branches->GetEntries(); k++) {
        TString bname = branches->At(k)->GetName();
        if (!bname.EqualTo(arr->At(i)->GetName()))
          continue;
        fTree[j]->SetBranchStatus(bname, 0);
        found = kTRUE;
      }
    }
    if (!found)
      AliFatal(Form("Did not find Branch %s", arr->At(i)->GetName()));
  }
  fPruneList = "";
} // void AliAnalysisTaskAO2Dconverter::Prune()

void AliAnalysisTaskAO2Dconverter::FillTree(TreeIndex t)
{
  if (!fTreeStatus[t]) return;
  Int_t nbytes = fTree[t]->Fill();
  if (nbytes > 0) fBytes += nbytes;
} // void AliAnalysisTaskAO2Dconverter::FillTree(TreeIndex t)

void AliAnalysisTaskAO2Dconverter::WriteTree(TreeIndex t)
{
  if (!fTreeStatus[t]) return;
  // Write the tree in the corrsponding (TF) directory
  if (!fOutputDir) AliFatal("No Root subdir|");
  fOutputDir->cd();
  AliInfo(Form("Writing tree %s\n", TreeName[t].Data()));
  fTree[t]->Write();
} // void AliAnalysisTaskAO2Dconverter::WriteTree(TreeIndex t)

void AliAnalysisTaskAO2Dconverter::InitTF(ULong64_t tfId)
{
  // Reset the event count
  fTFCount++;
  fEventCount = 0;
  fTfInitialized = true;
  
  // Reset the offsets
  fOffsetMuTrackID = 0;
  fOffsetTrackID = 0;
  fOffsetV0ID = 0;
  fOffsetLabel = 0;

  // Reset the content of eventextra
  for (auto i = 0; i < kTrees; ++i) {
     eventextra.fStart[i] = 0;
     eventextra.fNentries[i] = 0;
  }

  // Create the output directory for the current time frame
  fOutputDir = fOutputFile->mkdir(Form("TF_%llu", tfId));


  // Associate branches for fEventTree
  TTree* tEvents = CreateTree(kEvents);
  if (fTreeStatus[kEvents]) {
    tEvents->Branch("fBCsID", &collision.fBCsID, "fBCsID/I");
    tEvents->Branch("fPosX", &collision.fPosX, "fPosX/F");
    tEvents->Branch("fPosY", &collision.fPosY, "fPosY/F");
    tEvents->Branch("fPosZ", &collision.fPosZ, "fPosZ/F");
    tEvents->Branch("fCovXX", &collision.fCovXX, "fCovXX/F");
    tEvents->Branch("fCovXY", &collision.fCovXY, "fCovXY/F");
    tEvents->Branch("fCovXZ", &collision.fCovXZ, "fCovXZ/F");
    tEvents->Branch("fCovYY", &collision.fCovYY, "fCovYY/F");
    tEvents->Branch("fCovYZ", &collision.fCovYZ, "fCovYZ/F");
    tEvents->Branch("fCovZZ", &collision.fCovZZ, "fCovZZ/F");
    tEvents->Branch("fChi2", &collision.fChi2, "fChi2/F");
    tEvents->Branch("fNumContrib", &collision.fN, "fNumContrib/i");
    tEvents->Branch("fCollisionTime", &collision.fCollisionTime, "fCollisionTime/F");
    tEvents->Branch("fCollisionTimeRes", &collision.fCollisionTimeRes, "fCollisionTimeRes/F");
    tEvents->Branch("fCollisionTimeMask", &collision.fCollisionTimeMask, "fCollisionTimeMask/b");
    tEvents->SetBasketSize("*", fBasketSizeEvents);
  }
  
  // Extra information for debugging for event table
  TTree* tEventsExtra = CreateTree(kEventsExtra);
  if (fTreeStatus[kEventsExtra]) {
    TString sstart = TString::Format("fStart[%d]/I", kTrees);
    TString sentries = TString::Format("fNentries[%d]/I", kTrees);
    tEventsExtra->Branch("fStart", eventextra.fStart, sstart.Data());
    tEventsExtra->Branch("fNentries", eventextra.fNentries, sentries.Data());
    tEventsExtra->SetBasketSize("*", fBasketSizeEvents);
  }

  // Associate branches for fEventTree
  TTree* tBC = CreateTree(kBC);
  if (fTreeStatus[kBC]) {
    tBC->Branch("fRunNumber", &bc.fRunNumber, "fRunNumber/I");
    tBC->Branch("fGlobalBC", &bc.fGlobalBC, "fGlobalBC/l");
    tBC->Branch("fTriggerMask", &bc.fTriggerMask, "fTriggerMask/l");
    tBC->SetBasketSize("*", fBasketSizeEvents);
  }
  
  // Associate branches for fTrackTree
  TTree* tTracks = CreateTree(kTracks);
  if (fTreeStatus[kTracks]) {
    tTracks->Branch("fCollisionsID", &tracks.fCollisionsID, "fCollisionsID/I");
    tTracks->Branch("fTrackType", &tracks.fTrackType, "fTrackType/b");
    //    tTracks->Branch("fTOFclsIndex", &tracks.fTOFclsIndex, "fTOFclsIndex/I");
    //    tTracks->Branch("fNTOFcls", &tracks.fNTOFcls, "fNTOFcls/I");
    tTracks->Branch("fX", &tracks.fX, "fX/F");
    tTracks->Branch("fAlpha", &tracks.fAlpha, "fAlpha/F");
    tTracks->Branch("fY", &tracks.fY, "fY/F");
    tTracks->Branch("fZ", &tracks.fZ, "fZ/F");
    tTracks->Branch("fSnp", &tracks.fSnp, "fSnp/F");
    tTracks->Branch("fTgl", &tracks.fTgl, "fTgl/F");
    tTracks->Branch("fSigned1Pt", &tracks.fSigned1Pt, "fSigned1Pt/F");
    // Modified covariance matrix
    tTracks->Branch("fSigmaY", &tracks.fSigmaY, "fSigmaY/F");
    tTracks->Branch("fSigmaZ", &tracks.fSigmaZ, "fSigmaZ/F");
    tTracks->Branch("fSigmaSnp", &tracks.fSigmaSnp, "fSigmaSnp/F");
    tTracks->Branch("fSigmaTgl", &tracks.fSigmaTgl, "fSigmaTgl/F");
    tTracks->Branch("fSigma1Pt", &tracks.fSigma1Pt, "fSigma1Pt/F");
    tTracks->Branch("fRhoZY", &tracks.fRhoZY, "fRhoZY/B");
    tTracks->Branch("fRhoSnpY", &tracks.fRhoSnpY, "fRhoSnpY/B");
    tTracks->Branch("fRhoSnpZ", &tracks.fRhoSnpZ, "fRhoSnpZ/B");
    tTracks->Branch("fRhoTglY", &tracks.fRhoTglY, "fRhoTglY/B");
    tTracks->Branch("fRhoTglZ", &tracks.fRhoTglZ, "fRhoTglZ/B");
    tTracks->Branch("fRhoTglSnp", &tracks.fRhoTglSnp, "fRhoTglSnp/B");
    tTracks->Branch("fRho1PtY", &tracks.fRho1PtY, "fRho1PtY/B");
    tTracks->Branch("fRho1PtZ", &tracks.fRho1PtZ, "fRho1PtZ/B");
    tTracks->Branch("fRho1PtSnp", &tracks.fRho1PtSnp, "fRho1PtSnp/B");
    tTracks->Branch("fRho1PtTgl", &tracks.fRho1PtTgl, "fRho1PtTgl/B");
    //
    tTracks->Branch("fTPCInnerParam", &tracks.fTPCinnerP, "fTPCInnerParam/F");
    tTracks->Branch("fFlags", &tracks.fFlags, "fFlags/i");
    tTracks->Branch("fITSClusterMap", &tracks.fITSClusterMap, "fITSClusterMap/b");
    tTracks->Branch("fTPCNClsFindable", &tracks.fTPCNClsFindable, "fTPCNClsFindable/b");
    tTracks->Branch("fTPCNClsFindableMinusFound",&tracks.fTPCNClsFindableMinusFound, "fTPCNClsFindableMinusFound/B");
    tTracks->Branch("fTPCNClsFindableMinusCrossedRows", &tracks.fTPCNClsFindableMinusCrossedRows, "fTPCNClsFindableMinusCrossedRows/B");
    tTracks->Branch("fTPCNClsShared", &tracks.fTPCNClsShared, "fTPCNClsShared/b");
    tTracks->Branch("fTRDPattern", &tracks.fTRDPattern, "fTRDPattern/b");
    tTracks->Branch("fITSChi2NCl", &tracks.fITSChi2NCl, "fITSChi2NCl/F");
    tTracks->Branch("fTPCChi2NCl", &tracks.fTPCChi2NCl, "fTPCChi2NCl/F");
    tTracks->Branch("fTRDChi2", &tracks.fTRDChi2, "fTRDChi2/F");
    tTracks->Branch("fTOFChi2", &tracks.fTOFChi2, "fTOFChi2/F");
    tTracks->Branch("fTPCSignal", &tracks.fTPCSignal, "fTPCSignal/F");
    tTracks->Branch("fTRDSignal", &tracks.fTRDSignal, "fTRDSignal/F");
    tTracks->Branch("fTOFSignal", &tracks.fTOFSignal, "fTOFSignal/F");
    tTracks->Branch("fLength", &tracks.fLength, "fLength/F");
    tTracks->Branch("fTOFExpMom", &tracks.fTOFExpMom, "fTOFExpMom/F");
    tTracks->Branch("fTrackEtaEMCAL", &tracks.fTrackEtaEMCAL, "fTrackEtaEMCAL/F");
    tTracks->Branch("fTrackPhiEMCAL", &tracks.fTrackPhiEMCAL, "fTrackPhiEMCAL/F");
    tTracks->SetBasketSize("*", fBasketSizeTracks);
  }

  // Associate branches for Calo
  TTree* tCalo = CreateTree(kCalo);
  if (fTreeStatus[kCalo]) {
    tCalo->Branch("fBCsID", &calo.fBCsID, "fBCsID/I");
    tCalo->Branch("fCellNumber", &calo.fCellNumber, "fCellNumber/S");
    tCalo->Branch("fAmplitude", &calo.fAmplitude, "fAmplitude/F");
    tCalo->Branch("fTime", &calo.fTime, "fTime/F");
    tCalo->Branch("fCellType", &calo.fCellType, "fCellType/B");
    tCalo->Branch("fCaloType", &calo.fCaloType, "fCaloType/B");
    tCalo->SetBasketSize("*", fBasketSizeEvents);
  }

  TTree *tCaloTrigger = CreateTree(kCaloTrigger);
  if (fTreeStatus[kCaloTrigger]) {
    tCaloTrigger->Branch("fBCsID", &calotrigger.fBCsID, "fBCsID/I");
    tCaloTrigger->Branch("fFastOrAbsID", &calotrigger.fFastOrAbsID, "fFastOrAbsID/S");
    tCaloTrigger->Branch("fL0Amplitude", &calotrigger.fL0Amplitude, "fL0Amplitude/F");
    tCaloTrigger->Branch("fL1TimeSum", &calotrigger.fL1TimeSum, "fL1TimeSum/F");
    tCaloTrigger->Branch("fNL0Times", &calotrigger.fNL0Times, "fNL0Times/B");
    tCaloTrigger->Branch("fTriggerBits", &calotrigger.fTriggerBits, "fTriggerBits/I");
    tCaloTrigger->Branch("fCaloType", &calotrigger.fCaloType, "fCaloType/B");
    tCaloTrigger->SetBasketSize("*", fBasketSizeEvents);
  }

  // Associuate branches for MUON tracks
  TTree* tMuon = CreateTree(kMuon);
  if (fTreeStatus[kMuon]) {
    tMuon->Branch("fBCsID", &muons.fBCsID, "fBCsID/I");
//    tMuon->Branch("fClusterIndex", &muons.fClusterIndex, "fClusterIndex/I");
//    tMuon->Branch("fNclusters", &muons.fNclusters, "fNclusters/I");
    tMuon->Branch("fInverseBendingMomentum", &muons.fInverseBendingMomentum, "fInverseBendingMomentum/F");
    tMuon->Branch("fThetaX", &muons.fThetaX, "fThetaX/F");
    tMuon->Branch("fThetaY", &muons.fThetaY, "fThetaY/F");
    tMuon->Branch("fZMu", &muons.fZMu, "fZMu/F");
    tMuon->Branch("fBendingCoor", &muons.fBendingCoor, "fBendingCoor/F");
    tMuon->Branch("fNonBendingCoor", &muons.fNonBendingCoor, "fNonBendingCoor/F");
    tMuon->Branch("fCovariances", muons.fCovariances, "fCovariances[15]/F");
    tMuon->Branch("fChi2", &muons.fChi2, "fChi2/F");
    tMuon->Branch("fChi2MatchTrigger", &muons.fChi2MatchTrigger, "fChi2MatchTrigger/F");
    tMuon->SetBasketSize("*", fBasketSizeEvents);
  }

  // Associate branches for MUON tracks
  TTree* tMuonCls = CreateTree(kMuonCls);
  if (fTreeStatus[kMuonCls]) {
    tMuonCls->Branch("fMuonsID",&mucls.fMuonsID,"fMuonsID/I");
    tMuonCls->Branch("fX",&mucls.fX,"fX/F");
    tMuonCls->Branch("fY",&mucls.fY,"fY/F");
    tMuonCls->Branch("fZ",&mucls.fZ,"fZ/F");
    tMuonCls->Branch("fErrX",&mucls.fErrX,"fErrX/F");
    tMuonCls->Branch("fErrY",&mucls.fErrY,"fErrY/F");
    tMuonCls->Branch("fCharge",&mucls.fCharge,"fCharge/F");
    tMuonCls->Branch("fChi2",&mucls.fChi2,"fChi2/F");
    tMuonCls->SetBasketSize("*", fBasketSizeEvents);
  }

  // Associuate branches for ZDC
  TTree* tZdc = CreateTree(kZdc);
  if (fTreeStatus[kZdc]) {
    tZdc->Branch("fBCsID",           &zdc.fBCsID          , "fBCsID/I");
    tZdc->Branch("fEnergyZEM1",      &zdc.fEnergyZEM1     , "fEnergyZEM1/F");
    tZdc->Branch("fEnergyZEM2",      &zdc.fEnergyZEM2     , "fEnergyZEM2/F");
    tZdc->Branch("fEnergyCommonZNA", &zdc.fEnergyCommonZNA, "fEnergyCommonZNA/F");
    tZdc->Branch("fEnergyCommonZNC", &zdc.fEnergyCommonZNC, "fEnergyCommonZNC/F");
    tZdc->Branch("fEnergyCommonZPA", &zdc.fEnergyCommonZPA, "fEnergyCommonZPA/F");
    tZdc->Branch("fEnergyCommonZPC", &zdc.fEnergyCommonZPC, "fEnergyCommonZPC/F");
    tZdc->Branch("fEnergySectorZNA", &zdc.fEnergySectorZNA, "fEnergySectorZNA[4]/F");
    tZdc->Branch("fEnergySectorZNC", &zdc.fEnergySectorZNC, "fEnergySectorZNC[4]/F");
    tZdc->Branch("fEnergySectorZPA", &zdc.fEnergySectorZPA, "fEnergySectorZPA[4]/F");
    tZdc->Branch("fEnergySectorZPC", &zdc.fEnergySectorZPC, "fEnergySectorZPC[4]/F");
    tZdc->Branch("fTimeZEM1",        &zdc.fTimeZEM1       , "fTimeZEM1/F");
    tZdc->Branch("fTimeZEM2",        &zdc.fTimeZEM2       , "fTimeZEM2/F");
    tZdc->Branch("fTimeZNA",         &zdc.fTimeZNA        , "fTimeZNA/F");
    tZdc->Branch("fTimeZNC",         &zdc.fTimeZNC        , "fTimeZNC/F");
    tZdc->Branch("fTimeZPA",         &zdc.fTimeZPA        , "fTimeZPA/F");
    tZdc->Branch("fTimeZPC",         &zdc.fTimeZPC        , "fTimeZPC/F");
    tZdc->SetBasketSize("*", fBasketSizeEvents);
  }  

  // Associate branches for V0A
  TTree* tFV0A = CreateTree(kFV0A);
  if (fTreeStatus[kFV0A]) {
    tFV0A->Branch("fBCsID", &fv0a.fBCsID, "fBCsID/I");
    tFV0A->Branch("fAmplitude", fv0a.fAmplitude, "fAmplitude[48]/F");
    tFV0A->Branch("fTime", &fv0a.fTime, "fTime/F");
    tFV0A->Branch("fTriggerMask", &fv0a.fTriggerMask, "fTriggerMask/b");
    tFV0A->SetBasketSize("*", fBasketSizeEvents);
  }

  // Associate branches for V0C
  TTree* tFV0C = CreateTree(kFV0C);
  if (fTreeStatus[kFV0C]) {
    tFV0C->Branch("fBCsID", &fv0c.fBCsID, "fBCsID/I");
    tFV0C->Branch("fAmplitude", fv0c.fAmplitude, "fAmplitude[32]/F");
    tFV0C->Branch("fTime", &fv0c.fTime, "fTime/F");
    tFV0C->SetBasketSize("*", fBasketSizeEvents);
  }

  // Associate branches for FT0
  TTree* tFT0 = CreateTree(kFT0);
  if (fTreeStatus[kFT0]) {
    tFT0->Branch("fBCsID", &ft0.fBCsID, "fBCsID/I");
    tFT0->Branch("fAmplitudeA", ft0.fAmplitudeA, "fAmplitudeA[96]/F");
    tFT0->Branch("fAmplitudeC", ft0.fAmplitudeC, "fAmplitudeC[112]/F");
    tFT0->Branch("fTimeA", &ft0.fTimeA, "fTimeA/F");
    tFT0->Branch("fTimeC", &ft0.fTimeC, "fTimeC/F");
    tFT0->Branch("fTriggerMask", &ft0.fTriggerMask, "fTriggerMask/b");
    tFT0->SetBasketSize("*", fBasketSizeEvents);
  }
  
  // Associate branches for FDD (AD)
  TTree* tFDD = CreateTree(kFDD);
  if (fTreeStatus[kFDD]) {
    tFDD->Branch("fBCsID", &fdd.fBCsID, "fBCsID/I");
    tFDD->Branch("fAmplitudeA", fdd.fAmplitudeA, "fAmplitudeA[4]/F");
    tFDD->Branch("fAmplitudeC", fdd.fAmplitudeC, "fAmplitudeC[4]/F");
    tFDD->Branch("fTimeA", &fdd.fTimeA, "fTimeA/F");
    tFDD->Branch("fTimeC", &fdd.fTimeC, "fTimeC/F");
    tFDD->Branch("fTriggerMask", &fdd.fTriggerMask, "fTriggerMask/b");
    tFDD->SetBasketSize("*", fBasketSizeEvents);
  }
  
  // Associuate branches for V0s
  TTree* tV0s = CreateTree(kV0s);
  if (fTreeStatus[kV0s]) {
    tV0s->Branch("fPosTrackID", &v0s.fPosTrackID, "fPosTrackID/I");
    tV0s->Branch("fNegTrackID", &v0s.fNegTrackID, "fNegTrackID/I");
    tV0s->SetBasketSize("*", fBasketSizeTracks);
  }

  // Associuate branches for cascades
  TTree* tCascades = CreateTree(kCascades);
  if (fTreeStatus[kCascades]) {
    tCascades->Branch("fV0sID", &cascs.fV0sID, "fV0sID/I");
    tCascades->Branch("fTracksID", &cascs.fTracksID, "fTracksID/I");
    tCascades->SetBasketSize("*", fBasketSizeTracks);
  }

#ifdef USE_TOF_CLUST
  // Associate branches for TOF
  TTree* TOF = CreateTree(kTOF);
  if (fTreeStatus[kTOF]) {
    TOF->Branch("fTOFChannel", &tofClusters.fTOFChannel, "fTOFChannel/I");
    TOF->Branch("fTOFncls", &tofClusters.fTOFncls, "fTOFncls/S");
    TOF->Branch("fDx", &tofClusters.fDx, "fDx/F");
    TOF->Branch("fDz", &tofClusters.fDz, "fDz/F");
    TOF->Branch("fToT", &tofClusters.fToT, "fToT/F");
    TOF->SetBasketSize("*", fBasketSizeEvents);
  }
#else
  DisableTree(kTOF);
#endif

  if (fTaskMode == kMC) {
    TTree * tMCvtx = CreateTree(kMcCollision);
    if(fTreeStatus[kMcCollision]) {
      tMCvtx->Branch("fBCsID", &mccollision.fBCsID, "fBCsID/I");
      tMCvtx->Branch("fGeneratorsID", &mccollision.fGeneratorsID, "fGeneratorsID/S");
      tMCvtx->Branch("fPosX", &mccollision.fPosX, "fPosX/F");
      tMCvtx->Branch("fPosY", &mccollision.fPosY, "fPosY/F");
      tMCvtx->Branch("fPosZ", &mccollision.fPosZ, "fPosZ/F");
      tMCvtx->Branch("fT", &mccollision.fT, "fT/F");
      tMCvtx->Branch("fWeight", &mccollision.fWeight, "fWeight/F");
      tMCvtx->Branch("fImpactParameter", &mccollision.fImpactParameter, "fImpactParameter/F");
      tMCvtx->SetBasketSize("*", fBasketSizeEvents);
    }

    // Associate branches for Kinematics
    TTree* Kinematics = CreateTree(kMcParticle);
    if (fTreeStatus[kMcParticle]) {
      Kinematics->Branch("fMcCollisionsID", &mcparticle.fMcCollisionsID, "fMcCollisionsID/I");

      Kinematics->Branch("fPdgCode", &mcparticle.fPdgCode, "fPdgCode/I");
      Kinematics->Branch("fStatusCode", &mcparticle.fStatusCode, "fStatusCode/I");
      Kinematics->Branch("fFlags", &mcparticle.fFlags, "fFlags/b");
      
      Kinematics->Branch("fMother0", &mcparticle.fMother0, "fMother0/I");
      Kinematics->Branch("fMother1", &mcparticle.fMother1, "fMother1/I");
      Kinematics->Branch("fDaughter0", &mcparticle.fDaughter0, "fDaughter0/I");
      Kinematics->Branch("fDaughter1", &mcparticle.fDaughter1, "fDaughter1/I");
      Kinematics->Branch("fWeight", &mcparticle.fWeight, "fWeight/F");
      
      Kinematics->Branch("fPx", &mcparticle.fPx, "fPx/F");
      Kinematics->Branch("fPy", &mcparticle.fPy, "fPy/F");
      Kinematics->Branch("fPz", &mcparticle.fPz, "fPz/F");
      Kinematics->Branch("fE", &mcparticle.fE, "fE/F");
      
      Kinematics->Branch("fVx", &mcparticle.fVx, "fVx/F");
      Kinematics->Branch("fVy", &mcparticle.fVy, "fVy/F");
      Kinematics->Branch("fVz", &mcparticle.fVz, "fVz/F");
      Kinematics->Branch("fVt", &mcparticle.fVt, "fVt/F");
      Kinematics->SetBasketSize("*", fBasketSizeTracks);
    }

    // MC labels of each reconstructed track
    TTree* tLabels = CreateTree(kMcTrackLabel);
    if (fTreeStatus[kMcTrackLabel]) {
      tLabels->Branch("fLabel", &mctracklabel.fLabel, "fLabel/i");
      tLabels->Branch("fLabelMask", &mctracklabel.fLabelMask, "fLabelMask/s");
      tLabels->SetBasketSize("*", fBasketSizeTracks);
    }

    // MC labels of each reconstructed calo cluster
    TTree* tCaloLabels = CreateTree(kMcCaloLabel);
    if (fTreeStatus[kMcCaloLabel]) {
      tCaloLabels->Branch("fLabel", &mccalolabel.fLabel, "fLabel/i");
      tCaloLabels->Branch("fLabelMask", &mccalolabel.fLabelMask, "fLabelMask/s");
      tCaloLabels->SetBasketSize("*", fBasketSizeEvents);
    }

    // MC labels of each reconstructed calo cluster
    TTree* tCollisionLabels = CreateTree(kMcCollisionLabel);
    if (fTreeStatus[kMcCaloLabel]) {
      tCollisionLabels->Branch("fLabel", &mccollisionlabel.fLabel, "fLabel/i");
      tCollisionLabels->Branch("fLabelMask", &mccollisionlabel.fLabelMask, "fLabelMask/s");
      tCollisionLabels->SetBasketSize("*", fBasketSizeEvents);
    }
  }

  Prune(); //Removing all unwanted branches (if any)
} // void AliAnalysisTaskAO2Dconverter::InitTF(Int_t tfId)

void AliAnalysisTaskAO2Dconverter::FillEventInTF()
{
  // Event counter
  Int_t eventID = fEventCount++;

  // Primary vertex
  const AliESDVertex * pvtx = fESD->GetPrimaryVertex();

  // Configuration of the PID response
  AliPIDResponse* PIDResponse = (AliPIDResponse*)((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetPIDResponse();
  PIDResponse->SetTOFResponse(fESD, AliPIDResponse::kBest_T0);
  AliTOFPIDResponse & TOFResponse = PIDResponse->GetTOFResponse();

  // Configuration of the MC event (if needed)
  AliMCEvent* MCEvt = nullptr;
  if (fTaskMode == kMC) {
    AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()); //Get the MC handler

    if (!eventHandler) //Check on the MC handler
      AliFatal("Could not retrieve MC event handler");
    MCEvt = eventHandler->MCEvent(); //Get the MC Event

    if (!MCEvt) // Check on the MC Event
      AliFatal("Could not retrieve MC event");
    PIDResponse->SetCurrentMCEvent(MCEvt); //Set The PID response on the current MC event
  }

  //---------------------------------------------------------------------------
  // Collision data

  // Adjust start indices for this event in all trees by adding the number of entries of the previous event
  for (auto i = 0; i < kTrees; ++i)
     eventextra.fStart[i] += eventextra.fNentries[i];

  eventextra.fNentries[kEvents] = 1;  // one entry per vertex
  collision.fBCsID = eventID;
  collision.fPosX = AliMathBase::TruncateFloatFraction(pvtx->GetX(), mCollisionPosition);
  collision.fPosY = AliMathBase::TruncateFloatFraction(pvtx->GetY(), mCollisionPosition);
  collision.fPosZ = AliMathBase::TruncateFloatFraction(pvtx->GetZ(), mCollisionPosition);

  Double_t covmatrix[6];
  pvtx->GetCovMatrix(covmatrix);

  collision.fCovXX = AliMathBase::TruncateFloatFraction(covmatrix[0], mCollisionPositionCov);
  collision.fCovXY = AliMathBase::TruncateFloatFraction(covmatrix[1], mCollisionPositionCov);
  collision.fCovXZ = AliMathBase::TruncateFloatFraction(covmatrix[2], mCollisionPositionCov);
  collision.fCovYY = AliMathBase::TruncateFloatFraction(covmatrix[3], mCollisionPositionCov);
  collision.fCovYZ = AliMathBase::TruncateFloatFraction(covmatrix[4], mCollisionPositionCov);
  collision.fCovZZ = AliMathBase::TruncateFloatFraction(covmatrix[5], mCollisionPositionCov);

  collision.fChi2 = AliMathBase::TruncateFloatFraction(pvtx->GetChi2(), mCollisionPositionCov);
  collision.fN = (pvtx->GetNDF()+3)/2;

  Float_t eventTime[10];
  Float_t eventTimeRes[10];
  Double_t eventTimeWeight[10];
  
  for (Int_t i = 0; i < TOFResponse.GetNmomBins(); i++) {
    if (i >= 10)
      AliFatal("Index is too high!");
    Float_t mom = (TOFResponse.GetMinMom(i) + TOFResponse.GetMaxMom(i)) / 2.f;
    eventTime[i] = TOFResponse.GetStartTime(mom)*1.e-3; // ps to ns
    eventTimeRes[i] = TOFResponse.GetStartTimeRes(mom)*1.e-3; // ps to ns
    eventTimeWeight[i] = 1./(eventTimeRes[i]*eventTimeRes[i]);

    //PH The part below is just a place holder
    if (TOFResponse.GetStartTimeMask(mom) & 0x1)
      SETBIT(collision.fCollisionTimeMask, 0);
    else
      CLRBIT(collision.fCollisionTimeMask, 0);
    //
    if (TOFResponse.GetStartTimeMask(mom) & 0x2)
      SETBIT(collision.fCollisionTimeMask, 1);
    else
      CLRBIT(collision.fCollisionTimeMask, 1);
    //
    if (TOFResponse.GetStartTimeMask(mom) & 0x3)
      SETBIT(collision.fCollisionTimeMask, 2);
    else
      CLRBIT(collision.fCollisionTimeMask, 2);
  }

  // Recalculate unique event time and its resolution
  collision.fCollisionTime = AliMathBase::TruncateFloatFraction(TMath::Mean(10,eventTime,eventTimeWeight), mCollisionPosition); // Weighted mean of times per momentum interval
  collision.fCollisionTimeRes = AliMathBase::TruncateFloatFraction(TMath::Sqrt(9./10.)*TMath::Mean(10,eventTimeRes), mCollisionPositionCov); // PH bad approximation

  //---------------------------------------------------------------------------
  // BC data
  
  bc.fRunNumber = fESD->GetRunNumber();
  
  ULong64_t evtid = GetGlobalBC(fESD->GetHeader());
  if(!evtid){
    evtid = (ULong64_t(fESD->GetTimeStamp())<<32) + ULong64_t((fESD->GetNumberOfTPCClusters()<<5)|(fESD->GetNumberOfTPCTracks()));
  }
  bc.fGlobalBC = evtid;
  
  bc.fTriggerMask = fESD->GetTriggerMask();
  TString firedClasses = fESD->GetFiredTriggerClasses();
  if (firedClasses.Contains("CINT7-B-NOPF-CENTNOTRD")) bc.fTriggerMask |= 1ull << 50;
  if (firedClasses.Contains("CCUP8-B-NOPF-CENTNOTRD")) bc.fTriggerMask |= 1ull << 51;
  if (firedClasses.Contains("CCUP9-B-NOPF-CENTNOTRD")) bc.fTriggerMask |= 1ull << 52;
  if (firedClasses.Contains("CMUP10-B-NOPF-MUFAST"))   bc.fTriggerMask |= 1ull << 53;
  if (firedClasses.Contains("CMUP11-B-NOPF-MUFAST"))   bc.fTriggerMask |= 1ull << 54;
  if (firedClasses.Contains("CINT7-B-NOPF-MUFAST"))    bc.fTriggerMask |= 1ull << 55;
  if (firedClasses.Contains("CMSL7-B-NOPF-MUFAST"))    bc.fTriggerMask |= 1ull << 56;
  if (firedClasses.Contains("CMLL7-B-NOPF-MUFAST"))    bc.fTriggerMask |= 1ull << 57;
  if (firedClasses.Contains("CMUL7-B-NOPF-MUFAST"))    bc.fTriggerMask |= 1ull << 58;
  if (firedClasses.Contains("CMSH7-B-NOPF-MUFAST"))    bc.fTriggerMask |= 1ull << 59;
  
  FillTree(kBC);
  //---------------------------------------------------------------------------
  // MC kinematic tree
  // It has to be done before the reconstructed tracks,
  // since the cleaning requires reasignment of the labels
  Int_t nkine_filled = 0; // Number of kine tracks filled

  TArrayC toWrite;
  TArrayI kineIndex;
  if (MCEvt) {
    // Kinematics
    TParticle* particle = nullptr;
    Int_t nMCtracks = MCEvt->GetNumberOfTracks();
    Int_t nMCprim = MCEvt->GetNumberOfPrimaries();

    toWrite.Reset();
    toWrite.Set(nMCtracks);
    kineIndex.Reset();
    kineIndex.Set(nMCtracks);

    // For each reconstructed track keep the corresponding MC particle
    Int_t ntracks = fESD->GetNumberOfTracks();

    for (Int_t itrack=0; itrack < ntracks; ++itrack) {
      AliESDtrack *track = fESD->GetTrack(itrack);
      Int_t alabel = track->GetLabel();
      toWrite[TMath::Abs(alabel)] = 1;
    }

    // For each calo cluster keep the corresponding MC particle
    AliESDCaloCells *emcalCells = fESD->GetEMCALCells();
    Short_t nEmcalCells = emcalCells->GetNumberOfCells();
    for (Short_t ice = 0; ice < nEmcalCells; ++ice) {
      Short_t cellNumber;
      Double_t amplitude;
      Double_t time;
      Int_t mclabel;
      Double_t efrac;

      emcalCells->GetCell(ice, cellNumber, amplitude, time, mclabel, efrac);
      toWrite[TMath::Abs(mclabel)] = 1;
    }
    AliESDCaloCells *phosCells = fESD->GetPHOSCells();
    Short_t nPhosCells = phosCells->GetNumberOfCells();
    for (Short_t ice = 0; ice < nPhosCells; ++ice) {
      Short_t cellNumber;
      Double_t amplitude;
      Double_t time;
      Int_t mclabel;
      Double_t efrac;

      phosCells->GetCell(ice, cellNumber, amplitude, time, mclabel, efrac);
      toWrite[TMath::Abs(mclabel)] = 1;
    }

    // For each tracklet keep the corresponding MC particle
    AliMultiplicity *mult = fESD->GetMultiplicity();
    Int_t ntracklets = mult->GetNumberOfTracklets();

    for (Int_t itr = ntracklets; itr--;) {
      Int_t alabel = mult->GetLabel(itr, 0); // Take the label of the first layer
      toWrite[TMath::Abs(alabel)] = 1;
    }

    // For each MUON track keep the corresponding MC particle
    Int_t nmuons = fESD->GetNumberOfMuonTracks();

    for (Int_t imu=0; imu<nmuons; ++imu) {
      AliESDMuonTrack* mutrk = fESD->GetMuonTrack(imu);
      Int_t alabel = mutrk->GetLabel();
      toWrite[TMath::Abs(alabel)] = 1;
    }

    // Define which MC particles pass kinematic selection and update toWrite index
    for (Int_t i = 0; i < nMCtracks; ++i) {
      AliVParticle* vpt = MCEvt->GetTrack(i);
      particle = vpt->Particle();

      Float_t xv = particle->Vx();
      Float_t yv = particle->Vy();
      Float_t zv = particle->Vz();
      Float_t rv = TMath::Sqrt(xv * xv + yv * yv);

      // This part is taken from AliAnalysisTaskMCParticleFilter
      Bool_t write = kFALSE;

      if (i < nMCprim) {
	// Select all primary particles
	write = kTRUE;
      } else if (particle->GetUniqueID() == kPDecay) {
	// Particles from decay
	// Check that the decay chain ends at a primary particle
	AliVParticle* mother = vpt;
	Int_t imo = vpt->GetMother();
	while((imo >= nMCprim) && (mother->Particle()->GetUniqueID() == kPDecay)) {
	  mother =  (AliMCParticle*) MCEvt->GetTrack(imo);
	  imo =  mother->GetMother();
	}
	// Select according to pseudorapidity and production point of primary ancestor
	if (imo < nMCprim) write = kTRUE;
      } else if (particle->GetUniqueID() == kPPair) {
	// Now look for pair production
	Int_t imo = vpt->GetMother();
	if (imo < nMCprim) {
	  // Select, if the gamma is a primary
	  write = kTRUE;
	} else {
	  // Check if the gamma comes from the decay chain of a primary particle
	  AliMCParticle* mother =  (AliMCParticle*) MCEvt->GetTrack(imo);
	  imo = mother->GetMother();
	  while((imo >= nMCprim) && (mother->Particle()->GetUniqueID() == kPDecay)) {
	    mother =   (AliMCParticle*) MCEvt->GetTrack(imo);
	    imo =  mother->GetMother();
	  }
	  // Select according to pseudorapidity and production point
	  if (imo < nMCprim && Select(mother->Particle(), rv, zv))
	    write = kTRUE;
	}
      }
      if (toWrite[i] > 0 || write) {
	toWrite[i] = 1;
	kineIndex[i] = nkine_filled;
	nkine_filled++;
      }
      else {
	kineIndex[i] = -1;
      }
    }

    for (Int_t i = 0; i < nMCtracks; ++i) { //loop on primary MC tracks Before Event Selection
      AliVParticle* vpt = MCEvt->GetTrack(i);
      particle = vpt->Particle();

      mcparticle.fMcCollisionsID = eventID;

      //Get the kinematic values of the particles
      mcparticle.fPdgCode = particle->GetPdgCode();
      mcparticle.fStatusCode = particle->GetStatusCode();
      mcparticle.fFlags = 0;
      if (i >= MCEvt->Stack()->GetNprimary())
        mcparticle.fFlags |= MCParticleFlags::ProducedInTransport;
      mcparticle.fMother0 = vpt->GetMother();
      if (mcparticle.fMother0 > -1)
	mcparticle.fMother0 = kineIndex[mcparticle.fMother0] > -1 ? kineIndex[mcparticle.fMother0]+fOffsetLabel : -1;
      mcparticle.fMother1 = -1;
      mcparticle.fDaughter0 = particle->GetFirstDaughter();
      if (mcparticle.fDaughter0 > -1)
	mcparticle.fDaughter0 = kineIndex[mcparticle.fDaughter0] > -1 ? kineIndex[mcparticle.fDaughter0]+fOffsetLabel : -1;
      mcparticle.fDaughter1 = particle->GetLastDaughter();
      if (mcparticle.fDaughter1 > -1)
	mcparticle.fDaughter1 = kineIndex[mcparticle.fDaughter1] > -1 ? kineIndex[mcparticle.fDaughter1]+fOffsetLabel : -1;
      mcparticle.fWeight = AliMathBase::TruncateFloatFraction(particle->GetWeight(), mMcParticleW);

      mcparticle.fPx = AliMathBase::TruncateFloatFraction(particle->Px(), mMcParticleMom);
      mcparticle.fPy = AliMathBase::TruncateFloatFraction(particle->Py(), mMcParticleMom);
      mcparticle.fPz = AliMathBase::TruncateFloatFraction(particle->Pz(), mMcParticleMom);
      mcparticle.fE  = AliMathBase::TruncateFloatFraction(particle->Energy(), mMcParticleMom);

      mcparticle.fVx = AliMathBase::TruncateFloatFraction(particle->Vx(), mMcParticlePos);
      mcparticle.fVy = AliMathBase::TruncateFloatFraction(particle->Vy(), mMcParticlePos);
      mcparticle.fVz = AliMathBase::TruncateFloatFraction(particle->Vz(), mMcParticlePos);
      mcparticle.fVt = AliMathBase::TruncateFloatFraction(particle->T(), mMcParticlePos);

      if (toWrite[i]>0) {
	FillTree(kMcParticle);
      }
    }
  }
  eventextra.fNentries[kMcParticle] = nkine_filled;
  
  //---------------------------------------------------------------------------
  // Track data

  Int_t ntrk = fESD->GetNumberOfTracks();

  Int_t ntrk_filled = 0;     // total number of tracks filled per event
  Int_t ntofcls_filled = 0;  // total number of TOF clusters filled per event

  for (Int_t itrk = 0; itrk < ntrk; itrk++)
  {
    AliESDtrack *track = fESD->GetTrack(itrk);
//    if (!fTrackFilter.IsSelected(track))
//      continue;

    tracks.fCollisionsID = eventID;
    tracks.fTrackType = TrackTypeEnum::Run2GlobalTrack;

    tracks.fX = AliMathBase::TruncateFloatFraction(track->GetX(), mTrackX);
    tracks.fAlpha = AliMathBase::TruncateFloatFraction(track->GetAlpha(), mTrackAlpha);

    tracks.fY = track->GetY(); // no lossy compression
    tracks.fZ = track->GetZ();
    tracks.fSnp = AliMathBase::TruncateFloatFraction(track->GetSnp(), mtrackSnp);
    tracks.fTgl = AliMathBase::TruncateFloatFraction(track->GetTgl(), mTrackTgl);
    tracks.fSigned1Pt = AliMathBase::TruncateFloatFraction(track->GetSigned1Pt(), mTrack1Pt);

    // Modified covariance matrix
    // First sigmas on the diagonal
    tracks.fSigmaY = AliMathBase::TruncateFloatFraction(TMath::Sqrt(track->GetSigmaY2()), mTrackCovDiag);
    tracks.fSigmaZ = AliMathBase::TruncateFloatFraction(TMath::Sqrt(track->GetSigmaZ2()), mTrackCovDiag);
    tracks.fSigmaSnp = AliMathBase::TruncateFloatFraction(TMath::Sqrt(track->GetSigmaSnp2()), mTrackCovDiag);
    tracks.fSigmaTgl = AliMathBase::TruncateFloatFraction(TMath::Sqrt(track->GetSigmaTgl2()), mTrackCovDiag);
    tracks.fSigma1Pt = AliMathBase::TruncateFloatFraction(TMath::Sqrt(track->GetSigma1Pt2()), mTrackCovDiag);
    //
    tracks.fRhoZY = (Char_t)(128.*track->GetSigmaZY()/tracks.fSigmaZ/tracks.fSigmaY);
    tracks.fRhoSnpY = (Char_t)(128.*track->GetSigmaSnpY()/tracks.fSigmaSnp/tracks.fSigmaY);
    tracks.fRhoSnpZ = (Char_t)(128.*track->GetSigmaSnpZ()/tracks.fSigmaSnp/tracks.fSigmaZ);
    tracks.fRhoTglY = (Char_t)(128.*track->GetSigmaTglY()/tracks.fSigmaTgl/tracks.fSigmaY);
    tracks.fRhoTglZ = (Char_t)(128.*track->GetSigmaTglZ()/tracks.fSigmaTgl/tracks.fSigmaZ);
    tracks.fRhoTglSnp = (Char_t)(128.*track->GetSigmaTglSnp()/tracks.fSigmaTgl/tracks.fSigmaSnp);
    tracks.fRho1PtY = (Char_t)(128.*track->GetSigma1PtY()/tracks.fSigma1Pt/tracks.fSigmaY);
    tracks.fRho1PtZ = (Char_t)(128.*track->GetSigma1PtZ()/tracks.fSigma1Pt/tracks.fSigmaZ);
    tracks.fRho1PtSnp = (Char_t)(128.*track->GetSigma1PtSnp()/tracks.fSigma1Pt/tracks.fSigmaSnp);
    tracks.fRho1PtTgl = (Char_t)(128.*track->GetSigma1PtTgl()/tracks.fSigma1Pt/tracks.fSigmaTgl);

    const AliExternalTrackParam *intp = track->GetInnerParam();
    tracks.fTPCinnerP = AliMathBase::TruncateFloatFraction((intp ? intp->GetP() : 0), mTrack1Pt); // Set the momentum to 0 if the track did not reach TPC

    // Compressing and reassigned flags. Keeping only the ones we need.
    tracks.fFlags = 0x0;
    if (track->GetStatus() & AliVTrack::kITSrefit)
      tracks.fFlags |= TrackFlagsRun2Enum::ITSrefit;
    if (track->GetStatus() & AliVTrack::kTPCrefit)
      tracks.fFlags |= TrackFlagsRun2Enum::TPCrefit;

    // add status bit if golden chi2 cut was passed
    const AliESDVertex* vertex = (fESD->GetPrimaryVertex()) ? fESD->GetPrimaryVertex() : fESD->GetPrimaryVertexSPD();
    bool goldenChi2Status = (vertex) ? (track->GetChi2TPCConstrainedVsGlobal(vertex) > 0. && track->GetChi2TPCConstrainedVsGlobal(vertex) < 36.) : false;
    if (goldenChi2Status) 
      tracks.fFlags |= TrackFlagsRun2Enum::GoldenChi2;

    tracks.fITSClusterMap = track->GetITSClusterMap();
    tracks.fTPCNClsFindable = track->GetTPCNclsF();
    
    if ((int) tracks.fTPCNClsFindable - track->GetTPCNcls() >= -128)
      tracks.fTPCNClsFindableMinusFound = tracks.fTPCNClsFindable - track->GetTPCNcls();
    else
      tracks.fTPCNClsFindableMinusFound = -128;
    
    if ((int) tracks.fTPCNClsFindable - track->GetTPCCrossedRows() >= -128)
      tracks.fTPCNClsFindableMinusCrossedRows = tracks.fTPCNClsFindable - track->GetTPCCrossedRows();
    else
      tracks.fTPCNClsFindableMinusCrossedRows = -128;
    
    tracks.fTPCNClsShared = (track->GetTPCSharedMap()).CountBits();

    tracks.fTRDPattern = 0;
    for (int i=0;i<6;i++)
      if (track->GetTRDslice(i)>0)
        tracks.fTRDPattern |= 0x1<<i; // flag tracklet on this layer

    tracks.fITSChi2NCl = AliMathBase::TruncateFloatFraction((track->GetITSNcls() ? track->GetITSchi2() / track->GetITSNcls() : 0), mTrackCovOffDiag);
    tracks.fTPCChi2NCl = AliMathBase::TruncateFloatFraction((track->GetTPCNcls() ? track->GetTPCchi2() / track->GetTPCNcls() : 0), mTrackCovOffDiag);
    tracks.fTRDChi2 = AliMathBase::TruncateFloatFraction(track->GetTRDchi2(), mTrackCovOffDiag);
    tracks.fTOFChi2 = AliMathBase::TruncateFloatFraction(track->GetTOFchi2(), mTrackCovOffDiag);

    tracks.fTPCSignal = AliMathBase::TruncateFloatFraction(track->GetTPCsignal(), mTrackSignal);
    tracks.fTRDSignal = AliMathBase::TruncateFloatFraction(track->GetTRDsignal(), mTrackSignal);
    tracks.fTOFSignal = AliMathBase::TruncateFloatFraction(track->GetTOFsignal(), mTrackSignal);
    tracks.fLength = AliMathBase::TruncateFloatFraction(track->GetIntegratedLength(), mTrackSignal);

    // Speed of ligth in TOF units
    const Float_t cspeed = 0.029979246f;
    // PID hypothesis for the momentum extraction
    const AliPID::EParticleType tof_pid = AliPID::kPion;
    // Expected beta for such hypothesis
    const Float_t exp_beta =
        (track->GetIntegratedLength() /
         TOFResponse.GetExpectedSignal(track, tof_pid) / cspeed);

    tracks.fTOFExpMom = AliMathBase::TruncateFloatFraction(
        AliPID::ParticleMass(tof_pid) * exp_beta * cspeed /
            TMath::Sqrt(1. - (exp_beta * exp_beta)),
        mTrack1Pt);

    tracks.fTrackEtaEMCAL = AliMathBase::TruncateFloatFraction(track->GetTrackEtaOnEMCal(), mTrackPosEMCAL);
    tracks.fTrackPhiEMCAL = AliMathBase::TruncateFloatFraction(track->GetTrackPhiOnEMCal(), mTrackPosEMCAL);

    if (fTaskMode == kMC) {
      // Separate tables (trees) for the MC labels
      Int_t alabel = track->GetLabel();
      // Find the modified label
      Int_t klabel = kineIndex[TMath::Abs(alabel)];
      mctracklabel.fLabel = TMath::Abs(klabel) + fOffsetLabel;
      mctracklabel.fLabelMask = 0;
      // Use the ITS shared clusters to set the corresponding bits 0-6
      UChar_t itsMask = track->GetITSSharedMap() & 0x1F; // Normally only bits 0-5 are set in Run1/2
      mctracklabel.fLabelMask |= itsMask;
      // Use the number of TPC shared clusters as number of TPC mismatches
      // encode in bits 7-9 the values in the ranges 0, 1, 2-3, 4-7, 8-15, 16-31, 32-63, >64
      const TBits * tpcShared = track->GetTPCSharedMapPtr();
      UInt_t tpcCount = tpcShared->CountBits();
      UShort_t tpcMask = 0;
      while (tpcCount>0) {
	tpcCount = tpcCount >> 1;
	tpcMask++;
      }
      if (tpcMask>7) tpcMask = 7;
      mctracklabel.fLabelMask |= (tpcMask<<7);
      // TRD (bit 10)
      // We can also use labels per tracklet in the future
      Int_t trdLabel = track->GetTRDLabel();
      if (TMath::Abs(alabel)!=TMath::Abs(trdLabel)) mctracklabel.fLabelMask |= (0x1 << 10);
      // TOF (bit 11)
      Int_t tofLabel[3]={-1};
      track->GetTOFLabel(tofLabel);
      // Check if at least one of the TOF hits matches the track label
      if (!( TMath::Abs(alabel)==TMath::Abs(tofLabel[0])
	     || TMath::Abs(alabel)==TMath::Abs(tofLabel[1])
	     || TMath::Abs(alabel)==TMath::Abs(tofLabel[2])))
	mctracklabel.fLabelMask |= (0x1 << 11);

      if (alabel<0 || klabel<0) mctracklabel.fLabelMask |= (0x1 << 15);

      FillTree(kMcTrackLabel);
    }
  
#ifdef USE_TOF_CLUST
    tofClusters.fTOFncls = track->GetNTOFclusters();

    if (fTreeStatus[kTOF] && tofClusters.fTOFncls > 0) {
      Int_t* TOFclsIndex = track->GetTOFclusterArray(); //Index of the matchable cluster (there are fNTOFClusters of them)
      for (Int_t icls = 0; icls < tofClusters.fTOFncls; icls++) {
        AliESDTOFCluster* TOFcls = (AliESDTOFCluster*)fESD->GetESDTOFClusters()->At(TOFclsIndex[icls]);
        tofClusters.fToT = TOFcls->GetTOFsignalToT(0);
        tofClusters.fTOFChannel = TOFcls->GetTOFchannel();
        for (Int_t mtchbl = 0; mtchbl < TOFcls->GetNMatchableTracks(); mtchbl++) {
          if (TOFcls->GetTrackIndex(mtchbl) != track->GetID())
            continue;
          tofClusters.fDx = TOFcls->GetDx(mtchbl);
          tofClusters.fDz = TOFcls->GetDz(mtchbl);
          tofClusters.fLengthRatio = tracks.fLength > 0 ? TOFcls->GetLength(mtchbl) / tracks.fLength : -1;
          break;
        }
        FillTree(kTOF);
	if (fTreeStatus[kTOF]) ntofcls_filled++;
      }
    }
#endif

    // In case we need connection to clusters, activate next lines
    // tracks.fTOFclsIndex += tracks.fNTOFcls;
    // tracks.fNTOFcls = ntofcls_filled;
    FillTree(kTracks);
    if (fTreeStatus[kTracks]) ntrk_filled++;
  } // end loop on tracks

  eventextra.fNentries[kTOF]    = ntofcls_filled;

  AliMultiplicity *mlt = fESD->GetMultiplicity();
  Int_t Ntracklets = mlt->GetNumberOfTracklets();

  Int_t ntracklet_filled = 0;
  Float_t theta, phi, dphi, dphiS, dist, x, tgl, alpha;

  for (Int_t itr = Ntracklets; itr--;) {
    dphi   = mlt->GetDeltaPhi(itr);
    dist   = mlt->CalcDist(itr);
    
    // on-the-fly filtering based on parameters tuned in Run2
    dphiS  = TMath::Abs(dphi) - 0.0045; 
    if (dphi<0) dphiS = -dphiS;
    if (dist<1. && dphiS<0.06) {
      theta  = mlt->GetTheta(itr);
      phi = mlt->GetPhi(itr);
      tracks.fCollisionsID = eventID;
      tracks.fTrackType = TrackTypeEnum::Run2Tracklet;
      
      // inversion formulas for snp and alpha
      tracks.fSnp = 0.;
      alpha = phi;
      tracks.fAlpha = AliMathBase::TruncateFloatFraction(alpha, mTracklets);

      // inversion formulas for tgl
      x = (TMath::Tan(theta/2.)-1.) / (TMath::Tan(theta/2.)+1.);
      if (TMath::Log(TMath::Tan(theta/2)) >= 0)
        tgl = TMath::Sqrt((TMath::Power((1.+TMath::Power(x,2))/(1.-TMath::Power(x,2)),2))-1.);
      else 
        tgl = - TMath::Sqrt((TMath::Power((1.+TMath::Power(x,2))/(1.-TMath::Power(x,2)),2))-1.);
      tracks.fTgl = AliMathBase::TruncateFloatFraction(tgl, mTracklets);
    
      // set global track parameters to NAN
      tracks.fX = NAN;
      tracks.fY = NAN;
      tracks.fZ = NAN; 
      tracks.fSigned1Pt = NAN;
      tracks.fSigmaY = NAN;
      tracks.fSigmaZ = NAN;
      tracks.fSigmaSnp = NAN;
      tracks.fSigmaTgl = NAN;
      tracks.fSigma1Pt = NAN;
      tracks.fRhoZY = 0;
      tracks.fRhoSnpY = 0;
      tracks.fRhoSnpZ = 0;
      tracks.fRhoTglY = 0;
      tracks.fRhoTglZ = 0;
      tracks.fRhoTglSnp = 0;
      tracks.fRho1PtY = 0;
      tracks.fRho1PtZ = 0;
      tracks.fRho1PtSnp = 0;
      tracks.fRho1PtTgl = 0;
      tracks.fTPCinnerP = NAN; 
      tracks.fFlags = 0;
      tracks.fITSClusterMap = 0;
      tracks.fTPCNClsFindable = 0;
      tracks.fTPCNClsFindableMinusFound = 0;
      tracks.fTPCNClsFindableMinusCrossedRows = 0;
      tracks.fTPCNClsShared = 0;
      tracks.fTRDPattern = 0;
      tracks.fITSChi2NCl = NAN;
      tracks.fTPCChi2NCl = NAN;
      tracks.fTRDChi2 = NAN; 
      tracks.fTOFChi2 = NAN;
      tracks.fTPCSignal = NAN; 
      tracks.fTRDSignal = NAN;
      tracks.fTOFSignal = NAN;
      tracks.fLength = NAN;
      tracks.fTOFExpMom = NAN;

      if (fTaskMode == kMC) {
	// Separate tables (trees) for the MC labels: tracklets
	Int_t alabel = mlt->GetLabel(itr, 0); // Take the label of the first layer
	// Find the modified label
	Int_t klabel = kineIndex[TMath::Abs(alabel)];
	mctracklabel.fLabel = TMath::Abs(klabel) + fOffsetLabel;
	mctracklabel.fLabelMask = 0;
	// Mask fake tracklets
	if (alabel<0 || klabel<0) mctracklabel.fLabelMask |= (0x1 << 15);
	if (mlt->GetLabel(itr, 0) != mlt->GetLabel(itr, 1)) mctracklabel.fLabelMask |= (0x1 << 15);

	FillTree(kMcTrackLabel);
      }

      FillTree(kTracks);
      if (fTreeStatus[kTracks]) ntracklet_filled++;
    }
  } // end loop on tracklets
  eventextra.fNentries[kTracks] = ntrk_filled + ntracklet_filled; 
  
  //---------------------------------------------------------------------------
  // Calorimeter data

  const double kSecToNanoSec = 1e9;
  AliESDCaloCells *cells = fESD->GetEMCALCells();
  Short_t nCells = cells->GetNumberOfCells();
  Int_t ncalocells_filled = 0; // total number of calo cells filled per event
  for (Short_t ice = 0; ice < nCells; ++ice)
  {
    Short_t cellNumber;
    Double_t amplitude;
    Double_t time;
    Int_t mclabel;
    Double_t efrac;

    calo.fBCsID = eventID;
    
    cells->GetCell(ice, cellNumber, amplitude, time, mclabel, efrac);
    calo.fCellNumber = cellNumber;
    // Mimic run3 compression: Store only cells with energy larger than the threshold
    if(amplitude < fEMCALAmplitudeThreshold) continue;
    calo.fAmplitude = AliMathBase::TruncateFloatFraction(amplitude, mCaloAmp);
    calo.fTime = AliMathBase::TruncateFloatFraction(time * kSecToNanoSec, mCaloAmp);
    calo.fCaloType = cells->GetType(); // common for all cells
    calo.fCellType = cells->GetHighGain(ice) ? 1. : 0.; 
    FillTree(kCalo);
    if (fTreeStatus[kCalo]) ncalocells_filled++;
    if (fTaskMode == kMC) {
      // Find the modified label
      Int_t klabel = kineIndex[TMath::Abs(mclabel)];
      mccalolabel.fLabel = TMath::Abs(klabel) + fOffsetLabel;
      mccalolabel.fLabelMask = 0;
      if (mclabel<0 || klabel<0) mccalolabel.fLabelMask |= (0x1 << 15);

      FillTree(kMcCaloLabel);
    }
  } // end loop on calo cells
  eventextra.fNentries[kCalo] = ncalocells_filled;

  AliEMCALGeometry *geo = AliEMCALGeometry::GetInstanceFromRunNumber(fESD->GetRunNumber()); // Needed for EMCAL trigger mapping
  AliESDCaloTrigger *calotriggers = fESD->GetCaloTrigger("EMCAL");
  calotriggers->Reset();
  Int_t ncalotriggers_filled = 0; // total number of EMCAL triggers filled per event
  while(calotriggers->Next()){
    calotrigger.fBCsID = eventID;
    int col, row, fastorID;
    calotriggers->GetPosition(col, row);
    // filter null entries: they usually have negative entries and no trigger bits
    // in case of trigger bits the energy can be 0 or negative but the trigger position is marked
    int l1timesum, triggerbits;
    calotriggers->GetTriggerBits(triggerbits);
    calotriggers->GetL1TimeSum(l1timesum);
    if(!triggerbits && l1timesum <= 0) continue;
    // store trigger
    geo->GetTriggerMapping()->GetAbsFastORIndexFromPositionInEMCAL(col, row, fastorID);
    calotrigger.fFastOrAbsID = fastorID;
    calotriggers->GetAmplitude(calotrigger.fL0Amplitude);
    calotrigger.fL0Amplitude = AliMathBase::TruncateFloatFraction(calotrigger.fL0Amplitude, mCaloAmp);
    calotrigger.fL1TimeSum = AliMathBase::TruncateFloatFraction(l1timesum, mCaloAmp);
    calotriggers->GetTime(calotrigger.fL0Time);
    calotrigger.fL0Time = AliMathBase::TruncateFloatFraction(calotrigger.fL0Time, mCaloTime);
    calotriggers->GetTriggerBits(calotrigger.fTriggerBits);
    Int_t nL0times;
    calotriggers->GetNL0Times(nL0times);
    calotrigger.fNL0Times = nL0times;
    calotrigger.fTriggerBits = triggerbits;
    calotrigger.fCaloType = 1;
    FillTree(kCaloTrigger);
    if (fTreeStatus[kCaloTrigger]) ncalotriggers_filled++;
  }
  eventextra.fNentries[kCaloTrigger] = ncalotriggers_filled;

  cells = fESD->GetPHOSCells();
  nCells = cells->GetNumberOfCells();
  Int_t nphoscells_filled = 0;
  for (Short_t icp = 0; icp < nCells; ++icp)
  {
    Short_t cellNumber;
    Double_t amplitude;
    Double_t time;
    Int_t mclabel;
    Double_t efrac;

    calo.fBCsID = eventID;
    
    cells->GetCell(icp, cellNumber, amplitude, time, mclabel, efrac);
    calo.fCellNumber = cellNumber;
    calo.fAmplitude = AliMathBase::TruncateFloatFraction(amplitude, mCaloAmp);
    calo.fTime = AliMathBase::TruncateFloatFraction(time, mCaloTime);
    calo.fCellType = cells->GetHighGain(icp) ? 0. : 1.;     /// @TODO cell type value to be confirmed by PHOS experts
    calo.fCaloType = cells->GetType(); // common for all cells

    FillTree(kCalo);
    if (fTreeStatus[kCalo]) nphoscells_filled++;
    if (fTaskMode == kMC) {
      // Find the modified label
      Int_t klabel = kineIndex[TMath::Abs(mclabel)];
      mccalolabel.fLabel = TMath::Abs(klabel) + fOffsetLabel;
      mccalolabel.fLabelMask = 0;
      if (mclabel<0 || klabel<0) mccalolabel.fLabelMask |= (0x1 << 15);

      FillTree(kMcCaloLabel);
    }
  } // end loop on PHOS cells
  eventextra.fNentries[kCalo] = nphoscells_filled;

  //---------------------------------------------------------------------------
  // Muon tracks
  muons.fBCsID  = eventID;
  
  Int_t nmu = fESD->GetNumberOfMuonTracks();
  Int_t nmu_filled = 0;    // total number of muons filled per event
  Int_t nmucl_filled = 0;  // total number of clusters filled per event
  for (Int_t imu=0; imu<nmu; ++imu) {
    AliESDMuonTrack* mutrk = fESD->GetMuonTrack(imu);

    muons.fInverseBendingMomentum = AliMathBase::TruncateFloatFraction(mutrk->GetInverseBendingMomentum(), mMuonTr1P);
    muons.fThetaX = AliMathBase::TruncateFloatFraction(mutrk->GetThetaX(), mMuonTrThetaX);
    muons.fThetaY = AliMathBase::TruncateFloatFraction(mutrk->GetThetaY(), mMuonTrThetaY);
    muons.fZMu = AliMathBase::TruncateFloatFraction(mutrk->GetZ(), mMuonTrZmu);
    muons.fBendingCoor = AliMathBase::TruncateFloatFraction(mutrk->GetBendingCoor(), mMuonTrBend);
    muons.fNonBendingCoor = AliMathBase::TruncateFloatFraction(mutrk->GetNonBendingCoor(), mMuonTrNonBend);

    TMatrixD cov;
    mutrk->GetCovariances(cov);
    for (Int_t i = 0; i < 5; i++)
      for (Int_t j = 0; j <= i; j++)
	muons.fCovariances[i*(i+1)/2 + j] = AliMathBase::TruncateFloatFraction(cov(i,j), mMuonTrCov);

    muons.fChi2 = AliMathBase::TruncateFloatFraction(mutrk->GetChi2(), mMuonTrCov);
    muons.fChi2MatchTrigger = AliMathBase::TruncateFloatFraction(mutrk->GetChi2MatchTrigger(), mMuonTrCov);

    // Now MUON clusters for the current track
    Int_t muTrackID = fOffsetMuTrackID + imu;
    Int_t nmucl = mutrk->GetNClusters();
    for (Int_t imucl=0; imucl<nmucl; ++imucl){
      AliESDMuonCluster *muCluster = fESD->FindMuonCluster(mutrk->GetClusterId(imucl));
      mucls.fMuonsID = muTrackID;
      mucls.fX = AliMathBase::TruncateFloatFraction(muCluster->GetX(), mMuonCl);
      mucls.fY = AliMathBase::TruncateFloatFraction(muCluster->GetY(), mMuonCl);
      mucls.fZ = AliMathBase::TruncateFloatFraction(muCluster->GetZ(), mMuonCl);
      mucls.fErrX = AliMathBase::TruncateFloatFraction(muCluster->GetErrX(), mMuonClErr);
      mucls.fErrY = AliMathBase::TruncateFloatFraction(muCluster->GetErrY(), mMuonClErr);
      mucls.fCharge = AliMathBase::TruncateFloatFraction(muCluster->GetCharge(), mMuonCl);
      mucls.fChi2   = AliMathBase::TruncateFloatFraction(muCluster->GetChi2(), mMuonClErr);
      FillTree(kMuonCls);
      if (fTreeStatus[kMuonCls]) nmucl_filled++;
    } // End loop on muon clusters for the current muon track

    // In case we need connection to clusters, activate next lines
    // muons.fClusterIndex += muons.fNclusters;
    // muons.fNclusters = nmucl_filled;

    FillTree(kMuon);
    if (fTreeStatus[kMuon]) nmu_filled++;
  } // End loop on muon tracks
  eventextra.fNentries[kMuon] = nmu_filled;
  eventextra.fNentries[kMuonCls] = nmucl_filled;

  //---------------------------------------------------------------------------
  // ZDC
  AliESDZDC* esdzdc  =    fESD->GetESDZDC();
  zdc.fBCsID = eventID;
  // ZEM
  zdc.fEnergyZEM1      = esdzdc->GetZEM1Energy();
  zdc.fEnergyZEM2      = esdzdc->GetZEM2Energy();
  zdc.fEnergyCommonZNA = esdzdc->GetZNATowerEnergy()[0];
  zdc.fEnergyCommonZNC = esdzdc->GetZNCTowerEnergy()[0];
  zdc.fEnergyCommonZPA = esdzdc->GetZPATowerEnergy()[0];
  zdc.fEnergyCommonZPC = esdzdc->GetZPCTowerEnergy()[0];
  
  // ZDC (P,N) sectors
  for (Int_t ich=0; ich<4; ++ich) {
    zdc.fEnergySectorZNA[ich] = esdzdc->GetZNATowerEnergy()[ich+1];
    zdc.fEnergySectorZNC[ich] = esdzdc->GetZNCTowerEnergy()[ich+1];
    zdc.fEnergySectorZPA[ich] = esdzdc->GetZPATowerEnergy()[ich+1];
    zdc.fEnergySectorZPC[ich] = esdzdc->GetZPCTowerEnergy()[ich+1];
  }
  // ZDC TDC
  Bool_t isHitFlagFilled = fESD->GetRunNumber()>=208502;
  Bool_t isZNAhit  = isHitFlagFilled ? esdzdc->IsZNAhit() : 1;
  Bool_t isZNChit  = isHitFlagFilled ? esdzdc->IsZNChit() : 1;
  Bool_t isZPAhit  = isHitFlagFilled ? esdzdc->IsZPAhit() : 1;
  Bool_t isZPChit  = isHitFlagFilled ? esdzdc->IsZPChit() : 1;
  Bool_t isZEM1hit = isHitFlagFilled ? esdzdc->IsZEM1hit() : 1;
  Bool_t isZEM2hit = isHitFlagFilled ? esdzdc->IsZEM2hit() : 1;
  
  zdc.fTimeZNA  = 999.f;
  zdc.fTimeZNC  = 999.f;
  zdc.fTimeZPA  = 999.f;
  zdc.fTimeZPC  = 999.f;
  zdc.fTimeZEM1 = 999.f;
  zdc.fTimeZEM2 = 999.f;

  // Storing first ZDC hit in +/-12.5 ns around 0
  for (Int_t i=0;i<4;i++) {
    Float_t tZNA  = isZNAhit  ? esdzdc->GetZDCTDCCorrected(esdzdc->GetZNATDCChannel(),i)  : 999.f;
    Float_t tZNC  = isZNChit  ? esdzdc->GetZDCTDCCorrected(esdzdc->GetZNCTDCChannel(),i)  : 999.f;
    Float_t tZPA  = isZPAhit  ? esdzdc->GetZDCTDCCorrected(esdzdc->GetZPATDCChannel(),i)  : 999.f;
    Float_t tZPC  = isZPChit  ? esdzdc->GetZDCTDCCorrected(esdzdc->GetZPCTDCChannel(),i)  : 999.f;
    Float_t tZEM1 = isZEM1hit ? esdzdc->GetZDCTDCCorrected(esdzdc->GetZEM1TDCChannel(),i) : 999.f;
    Float_t tZEM2 = isZEM2hit ? esdzdc->GetZDCTDCCorrected(esdzdc->GetZEM2TDCChannel(),i) : 999.f;
    if (tZNA >-12.5 && tZNA <12.5 && zdc.fTimeZNA >998) zdc.fTimeZNA  = tZNA;
    if (tZNC >-12.5 && tZNC <12.5 && zdc.fTimeZNC >998) zdc.fTimeZNC  = tZNC;
    if (tZPA >-12.5 && tZPA <12.5 && zdc.fTimeZPA >998) zdc.fTimeZPA  = tZPA;
    if (tZPC >-12.5 && tZPC <12.5 && zdc.fTimeZPC >998) zdc.fTimeZPC  = tZPC;
    if (tZEM1>-12.5 && tZEM1<12.5 && zdc.fTimeZEM1>998) zdc.fTimeZEM1 = tZEM1;
    if (tZEM2>-12.5 && tZEM2<12.5 && zdc.fTimeZEM2>998) zdc.fTimeZEM2 = tZEM2;
  }
  
  FillTree(kZdc);
  if (fTreeStatus[kZdc]) eventextra.fNentries[kZdc] = 1;

  //---------------------------------------------------------------------------
  // V0A and V0C
  AliESDVZERO* vz = fESD->GetVZEROData();
  fv0a.fBCsID = eventID;
  fv0c.fBCsID = eventID;
  for (Int_t ich=0; ich<32; ++ich) fv0a.fAmplitude[ich] = AliMathBase::TruncateFloatFraction(vz->GetMultiplicityV0A(ich),mV0Amplitude);
  for (Int_t ich=0; ich<32; ++ich) fv0c.fAmplitude[ich] = AliMathBase::TruncateFloatFraction(vz->GetMultiplicityV0C(ich),mV0Amplitude);
  fv0a.fTime = AliMathBase::TruncateFloatFraction(vz->GetV0ATime(),mV0Time);
  fv0c.fTime = AliMathBase::TruncateFloatFraction(vz->GetV0CTime(),mV0Time);
  fv0a.fTriggerMask = 0; // not filled for the moment
  FillTree(kFV0A);
  FillTree(kFV0C);
  if (fTreeStatus[kFV0A]) eventextra.fNentries[kFV0A] = 1;
  if (fTreeStatus[kFV0C]) eventextra.fNentries[kFV0C] = 1;

  //---------------------------------------------------------------------------
  // FT0
  ft0.fBCsID = eventID;
  for (Int_t ich=0; ich<12; ++ich) ft0.fAmplitudeA[ich] = AliMathBase::TruncateFloatFraction(fESD->GetT0amplitude()[ich+12],mT0Amplitude);
  for (Int_t ich=0; ich<12; ++ich) ft0.fAmplitudeC[ich] = AliMathBase::TruncateFloatFraction(fESD->GetT0amplitude()[ich   ],mT0Amplitude);
  ft0.fTimeA = AliMathBase::TruncateFloatFraction(fESD->GetT0TOF(1)*1e-3,mT0Time); // ps to ns
  ft0.fTimeC = AliMathBase::TruncateFloatFraction(fESD->GetT0TOF(2)*1e-3,mT0Time); // ps to ns
  ft0.fTriggerMask = fESD->GetT0Trig();
  FillTree(kFT0);
  if (fTreeStatus[kFT0]) eventextra.fNentries[kFT0] = 1;
  
  //---------------------------------------------------------------------------
  // AD (FDD)
  AliESDAD* esdad = fESD->GetADData();
  fdd.fBCsID = eventID;
  for (Int_t ich=0; ich<4; ++ich) fdd.fAmplitudeA[ich] = 0; // not filled for the moment
  for (Int_t ich=0; ich<4; ++ich) fdd.fAmplitudeC[ich] = 0; // not filled for the moment
  fdd.fTimeA = AliMathBase::TruncateFloatFraction(esdad->GetADATime(),mADTime);
  fdd.fTimeC = AliMathBase::TruncateFloatFraction(esdad->GetADCTime(),mADTime);
  fdd.fTriggerMask = 0; // not filled for the moment
  FillTree(kFDD);
  if (fTreeStatus[kFDD]) eventextra.fNentries[kFDD] = 1;
  
  //---------------------------------------------------------------------------
  // V0s (Lambda and KS)
  Int_t nv0 = fESD->GetNumberOfV0s();
  Int_t nv0_filled = 0; // total number of v0's filled per event
  for (Int_t iv0=0; iv0<nv0; ++iv0) {
    AliESDv0 * v0 = fESD->GetV0(iv0);
    // select only "offline" V0s, skip the "on-the-fly" ones
    if (v0 && !v0->GetOnFlyStatus()) {
      Int_t pidx = v0->GetPindex();
      Int_t nidx = v0->GetNindex();
      v0s.fPosTrackID = TMath::Sign(TMath::Abs(pidx) + fOffsetTrackID, pidx); // Positive track ID
      v0s.fNegTrackID = TMath::Sign(TMath::Abs(nidx) + fOffsetTrackID, nidx); // Negative track ID
      FillTree(kV0s);
      if (fTreeStatus[kV0s]) nv0_filled++;
    }
  } // End loop on V0s
  eventextra.fNentries[kV0s] = nv0_filled;

  //---------------------------------------------------------------------------
  // Cascades
  // If we do not have V0s, we do not have cascades
  Int_t ncascades_filled = 0; // total number of cascades filled per event
  if (nv0>0) {
    // Combine the track indexes of V0 daughters in unique identifier
    ULong64_t * packedPosNeg = new ULong64_t[nv0];
    ULong64_t * sortedPosNeg = new ULong64_t[nv0];
    Int_t * sortIdx = new Int_t[nv0];

    //Don't forget that OTF V0s might exist
    Int_t nv0offline = 0;
    for (Int_t iv0=0; iv0<nv0; ++iv0) {
      AliESDv0 * v0 = fESD->GetV0(iv0);
      if (v0 && !v0->GetOnFlyStatus()){
          packedPosNeg[nv0offline++] = (((ULong64_t)(v0->GetPindex())) << 31) | ((ULong64_t)(v0->GetNindex()));
      }
    }
    TMath::Sort(nv0offline,packedPosNeg,sortIdx,kFALSE);
    for (Int_t iv0=0; iv0<nv0offline; ++iv0) {
      sortedPosNeg[iv0] = packedPosNeg[sortIdx[iv0]];
    }
  
    Int_t ncas = fESD->GetNumberOfCascades();
    for (Int_t icas=0; icas<ncas; ++icas) {
      AliESDcascade *cas = fESD->GetCascade(icas);
      // Select only cascades containing "offline" V0s
      if (cas && !cas->GetOnFlyStatus()) {
	// Find the identifier of the V0 using the indexes of its daughters
	ULong64_t currV0 = (((ULong64_t)(cas->GetPindex())) << 31) | ((ULong64_t)(cas->GetNindex()));
	// Use binary search in the sorted array
	Int_t v0idx = TMath::BinarySearch(nv0offline, sortedPosNeg, currV0);
	// Check if the match is exact
	if (sortedPosNeg[v0idx] == currV0) {
	  cascs.fV0sID = sortIdx[v0idx] + fOffsetV0ID;
	  cascs.fTracksID = cas->GetBindex() + fOffsetTrackID;
	  FillTree(kCascades);
	  if (fTreeStatus[kCascades]) ncascades_filled++;
	}
      }
    } // End loop on cascades

    delete [] packedPosNeg;
    delete [] sortedPosNeg;
    delete [] sortIdx;
  } // End if V0s
  eventextra.fNentries[kCascades] = ncascades_filled;
  
  //---------------------------------------------------------------------------
  // MC data (to be modified)


  if (MCEvt) {
    // MC vertex
    const AliVVertex* MCvtx = MCEvt->GetPrimaryVertex();
    if (!MCvtx) //Check on the MC vertex
      AliFatal("Could not retrieve MC vertex");

    mccollision.fBCsID = eventID;

    mccollision.fPosX = AliMathBase::TruncateFloatFraction(MCvtx->GetX(), mCollisionPosition);
    mccollision.fPosY = AliMathBase::TruncateFloatFraction(MCvtx->GetY(), mCollisionPosition);
    mccollision.fPosZ = AliMathBase::TruncateFloatFraction(MCvtx->GetZ(), mCollisionPosition);

    AliGenEventHeader* mcGenH = MCEvt->GenEventHeader();
    mccollision.fT = AliMathBase::TruncateFloatFraction(mcGenH->InteractionTime(), mCollisionPosition);
    mccollision.fWeight = AliMathBase::TruncateFloatFraction(mcGenH->EventWeight(), mCollisionPosition);

    // Impact parameter
    AliCollisionGeometry * cGeo = dynamic_cast<AliCollisionGeometry*>(mcGenH);
    mccollision.fImpactParameter = (cGeo ? cGeo->ImpactParameter() : -999.f);

    mccollision.fGeneratorsID = 0;
    for (Int_t gen = 0; gen < kGenerators; gen++) {
      if (mcGenH->InheritsFrom(Generator[gen]))
        SETBIT(mccollision.fGeneratorsID, gen);
      else
        CLRBIT(mccollision.fGeneratorsID, gen);
    }
    if (mcGenH->InheritsFrom(Generator[kAliGenCocktailEventHeader])) {
      TList* headers = ((AliGenCocktailEventHeader*)mcGenH)->GetHeaders();
      TListIter cocktail(headers);
      TObject *to = 0x0;
      while ((to = cocktail())) {
	if (mccollision.fImpactParameter < 0) {
	  // Change the impact parameter if not set
	  AliCollisionGeometry * toCGeo = dynamic_cast<AliCollisionGeometry*>(to);
	  mccollision.fImpactParameter = (toCGeo ? toCGeo->ImpactParameter() : -999.f);
	}
        for (Int_t gen = 0; gen < kGenerators; gen++) {
          if (to->InheritsFrom(Generator[gen]))
            SETBIT(mccollision.fGeneratorsID, gen);
        }
      }
    }
    mccollision.fImpactParameter = AliMathBase::TruncateFloatFraction(mccollision.fImpactParameter, mCollisionPosition);
    eventextra.fNentries[kMcCollision] = 1;
  } else {
    eventextra.fNentries[kMcCollision] = 0;
  }
  // Filling the tree of vertices has to be done last because it contains the
  // index data for the other trees
  FillTree(kMcCollision);

  // MC collision label
  mccollisionlabel.fLabel = eventID;
  mccollisionlabel.fLabelMask = 0;
  FillTree(kMcCollisionLabel);

  // We can fill now the vertex + indexing data
  FillTree(kEvents);

  //---------------------------------------------------------------------------
  // Update the offsets at the end of each collision    
  fOffsetLabel += nkine_filled; // Offset for the labels of the next event
  fOffsetTrackID += ntrk_filled + ntracklet_filled;
  fOffsetMuTrackID += nmu_filled;
  fOffsetV0ID += nv0_filled;
} // void AliAnalysisTaskAO2Dconverter::FillEventInTF()

void AliAnalysisTaskAO2Dconverter::FinishTF()
{
  // Write all trees
  for (Int_t i = 0; i < kTrees; i++)
    WriteTree((TreeIndex)i);
  // Remove trees
  for (Int_t i = 0; i < kTrees; i++)
    if (fTree[i]) {
      delete fTree[i];
      fTree[i] = 0x0;
    }
} // AliAnalysisTaskAO2Dconverter::FinishTF()

Bool_t AliAnalysisTaskAO2Dconverter::Select(TParticle* part, Float_t rv, Float_t zv)
{
  /// Selection accoring to eta of the mother and production point
  /// This has to be refined !!!!!!

  // Esp if we don't have collisison in the central barrel but e.g. beam gas
  //  return kTRUE;

  Float_t eta = part->Eta();

  // central barrel consider everything in the ITS...
  // large eta window for smeared vertex and SPD acceptance (2 at z = 0)
  // larger for V0s in the TPC
  //  if (TMath::Abs(eta) < 2.5 && rv < 250. && TMath::Abs(zv)<255)return kTRUE;

  if (TMath::Abs(eta) < 2.5 && rv < 170)return kTRUE;

  // Andreas' Cuts
  //  if (TMath::Abs(eta) < 1. && rv < 170)return kTRUE;



  // Muon arm
  if(eta > -4.2 && eta < -2.3 && zv > -500.)return kTRUE; // Muon arms

  // PMD acceptance 2.3 <= eta < = 3.5
  //  if(eta>2.0&&eta<3.7)return kTRUE;

  return kFALSE;

} // AliAnalysisTaskAO2Dconverter::Select(TParticle* part, Float_t rv, Float_t zv)

////////////////////////////////////////////////////////////
