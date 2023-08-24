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

#include <algorithm>
#include <TFile.h>
#include <TDirectory.h>
#include <TChain.h>
#include <TTree.h>
#include <TMath.h>
#include <Math/SMatrix.h>
#include <TTimeStamp.h>
#include <TClonesArray.h>
#include <TRandom.h>
#include <TH2.h>
#include <TSystem.h>
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliVVertex.h"
#include "AliESDVertex.h"
#include "AliAODVertex.h"
#include "AliCentrality.h"
#include "AliDAQ.h"
#include "AliGRPObject.h"
#include "AliVMultiplicity.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALTriggerDataGrid.h"
#include "AliEMCALTriggerPatchInfo.h"
#include "AliAnalysisTaskAO2Dconverter.h"
#include "AliVHeader.h"
#include "COMMON/MULTIPLICITY/AliMultSelection.h"
#include "limits.h"

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

#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"

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
#include "AliTriggerAnalysis.h"
#include "AliOADBContainer.h"
#include "AliMathBase.h"
#include "AliLog.h"
#include "AliProdInfo.h"

#include "RVersion.h"
#include "ARVersion.h"
// #include "APVersion.h"

ClassImp(AliAnalysisTaskAO2Dconverter);

const TString AliAnalysisTaskAO2Dconverter::TreeName[kTrees] = {
  "O2collision_001",
  "DbgEventExtra",
  "O2track",
  "O2trackcov",
  "O2trackextra",
  "O2fwdtrack",
  "O2fwdtrackcov",
  "O2calo",
  "O2calotrigger",
  "O2zdc",
  "O2fv0a",
  "O2fv0c",
  "O2ft0",
  "O2fdd_001",
  "O2v0_001",
  "O2cascade_001",
  "O2tof",
  "O2mcparticle_001",
  "O2mccollision",
  "O2mctracklabel",
  "O2mccalolabel",
  "O2mccollisionlabel",
  "O2bc",
  "O2run2bcinfo",
  "O2origin",
  "O2hmpid",
  "O2hf2prong",
  "O2hf3prong",
  "O2hfcascade",
  "O2hfdstar",
  "O2hepmcxsection",
  "O2hepmcpdfinfo",
  "O2hepmcheavyion"
};

const TString AliAnalysisTaskAO2Dconverter::TreeTitle[kTrees] = {
  "Collision tree",
  "Collision extra",
  "Barrel tracks Parameters",
  "Barrel tracks Covariance",
  "Barrel tracks Extra",
  "Forward tracks Parameters",
  "Forward tracks Covariances",
  "Calorimeter cells",
  "Calorimeter triggers",
  "ZDC",
  "FV0A",
  "FV0C",
  "FT0",
  "FDD",
  "V0s",
  "Cascades",
  "TOF hits",
  "Kinematics",
  "MC collisions",
  "MC track labels",
  "MC calo labels",
  "MC collision labels",
  "BC info",
  "Run 2 BC Info",
  "DF ids",
  "HMPID info",
  "HF 2 prong candidates",
  "HF 3 prong candidates",
  "HF cascade candidates",
  "HF D* candidates",
  "O2 HepMc Cross Sections",
  "O2 HepMc Pdf Info",
  "O2 HepMc Heavy Ion"
};

const TClass *AliAnalysisTaskAO2Dconverter::Generator[kGenerators] = {AliGenEventHeader::Class(), AliGenCocktailEventHeader::Class(), AliGenDPMjetEventHeader::Class(), AliGenEpos3EventHeader::Class(), AliGenEposEventHeader::Class(), AliGenEventHeaderTunedPbPb::Class(), AliGenGeVSimEventHeader::Class(), AliGenHepMCEventHeader::Class(), AliGenHerwigEventHeader::Class(), AliGenHijingEventHeader::Class(), AliGenPythiaEventHeader::Class(), AliGenToyEventHeader::Class()};

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

  UInt_t mTrackX = 0xFFFFFFFF;
  UInt_t mTrackAlpha = 0xFFFFFFFF;
  UInt_t mtrackSnp = 0xFFFFFFFF;
  UInt_t mTrackTgl = 0xFFFFFFFF;
  UInt_t mTrack1Pt = 0xFFFFFFFF;     // Including the momentun at the inner wall of TPC
  UInt_t mTrackCovDiag = 0xFFFFFFFF; // Including the chi2
  UInt_t mTrackCovOffDiag = 0xFFFFFFFF;
  UInt_t mTrackSignal = 0xFFFFFFFF; // PID signals and track length
  UInt_t mTrackPosEMCAL = 0xFFFFFFFF;

  UInt_t mTracklets = 0xFFFFFFFF; // tracklet members

  UInt_t mMcParticleW = 0xFFFFFFFF;   // Precision for weight
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
  : AliAnalysisTaskSE(name),
    fTrackFilter(Form("AO2Dconverter%s", name), Form("fTrackFilter%s", name)),
    fEventCuts{},
    fTriggerAnalysis(),
    collision(),
    eventextra(),
    bc(),
    run2bcinfo(),
    origin(),
    tracks(),
    hmpids(),
    mccollision(),
    mctracklabel(),
    mccalolabel(),
    mccollisionlabel(),
    mcparticle(),
    hepMcCrossSections(),
    hepMcPdfInfo(),
    hepMcHeavyIon(),
#ifdef USE_TOF_CLUST
    tofClusters(),
#endif
    calo(),
    calotrigger(),
    fwdtracks(),
    zdc(),
    fv0a(),
    fv0c(),
    ft0(),
    fdd(),
    v0s(),
    cascs(),
    hf2Prong(),
    hf3Prong(),
    hfCascades(),
    hfDStar(),
    fMetaData()
{
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  for (Int_t i = 0; i < kTrees; i++)
  {
    fTreeStatus[i] = kTRUE;
  }
} // AliAnalysisTaskAO2Dconverter::AliAnalysisTaskAO2Dconverter(const char* name)



AliAnalysisTaskAO2Dconverter::~AliAnalysisTaskAO2Dconverter()
{
  if (fOutputList) {
    fOutputList->Delete();
    delete fOutputList;
  }
} // AliAnalysisTaskAO2Dconverter::~AliAnalysisTaskAO2Dconverter()

void AliAnalysisTaskAO2Dconverter::NotifyRun(){
  const char* oadbfilename = Form("%s/COMMON/PHYSICSSELECTION/data/physicsSelection.root", AliAnalysisManager::GetOADBPath());
  TFile* oadb = TFile::Open(oadbfilename);
  if(!oadb->IsOpen()) AliFatal(Form("Cannot open OADB file %s", oadbfilename));
  AliOADBContainer* triggerContainer = (AliOADBContainer*) oadb->Get("trigAnalysis");
  if (!triggerContainer) AliFatal("Cannot fetch OADB container for trigger analysis");
  AliOADBTriggerAnalysis* oadbTriggerAnalysis = (AliOADBTriggerAnalysis*) triggerContainer->GetObject(fCurrentRunNumber, "Default");
  fTriggerAnalysis.SetParameters(oadbTriggerAnalysis);
  fTriggerAnalysis.ApplyPileupCuts(kTRUE);
  //read PHOS trigger bad map
  if (fUsePHOSBadMap){
    AliOADBContainer phosBadmapContainer(Form("phosTriggerBadMap"));
    phosBadmapContainer.InitFromFile(Form("%s/PHOS/PHOSTrigBadMaps.root", AliAnalysisManager::GetOADBPath()),
                                      "phosTriggerBadMap");
    TObjArray *maps = (TObjArray*)phosBadmapContainer.GetObject(fCurrentRunNumber,"phosTriggerBadMap");
    if(!maps){
      AliFatal(Form("Can not read PHOS Trigger Bad map for run %d. \n",fCurrentRunNumber)) ;
    }
    else{
      for(Int_t mod=0; mod<5;mod++){
        if(fPHOSBadMap[mod])
          delete fPHOSBadMap[mod] ;
        TH2I * h = (TH2I*)maps->At(mod) ;
        if(h)
          fPHOSBadMap[mod]=new TH2I(*h) ;
      }
    }
  }
  AliCDBManager *cdb = AliCDBManager::Instance();
  if(cdb->IsDefaultStorageSet() && cdb->GetRun() > 0) {
    // CDB manage initialized
    // Get GRP
    AliCDBEntry *en = cdb->Get("GRP/GRP/Data");
    if(en) {
      fGRP = static_cast<AliGRPObject *>(en->GetObject());
    }
  }
}


void AliAnalysisTaskAO2Dconverter::UserCreateOutputObjects()
{
  // Setting active/inactive containers based on the TaskMode
  switch (fTaskMode)
  {
  case kStandard:
    DisableTree(kMcParticle);
    DisableTree(kMcCollision);
    DisableTree(kMcTrackLabel);
    DisableTree(kMcCaloLabel);
    DisableTree(kMcCollisionLabel);
    DisableTree(kHepMcCrossSections);
    DisableTree(kHepMcPdfInfo);
    DisableTree(kHepMcHeavyIon);
    break;
  default:
    break;
  }

  // Set the truncation
  if (fTruncate)
  {
    mCollisionPosition = 0xFFFFFFF0;    // 19 bits mantissa
    mCollisionPositionCov = 0xFFFFE000; // 10 bits mantissa

    mTrackX = 0xFFFFFFF0;          // 19 bits
    mTrackAlpha = 0xFFFFFFF0;      // 19 bits
    mtrackSnp = 0xFFFFFF00;        // 15 bits
    mTrackTgl = 0xFFFFFF00;        // 15 bits
    mTrack1Pt = 0xFFFFFC00;        // 13 bits
    mTrackCovDiag = 0xFFFFFF00;    // 15 bits
    mTrackCovOffDiag = 0xFFFF0000; // 7 bits
    mTrackSignal = 0xFFFFFF00;     // 15 bits
    mTrackPosEMCAL = 0xFFFFFF00;   // 15 bits;

    mTracklets = 0xFFFFFF00; // 15 bits

    mMcParticleW = 0xFFFFFFF0;   // 19 bits
    mMcParticlePos = 0xFFFFFFF0; // 19 bits
    mMcParticleMom = 0xFFFFFFF0; // 19 bits

    mCaloAmp = 0xFFFFFF00;  // 15 bits
    mCaloTime = 0xFFFFFF00; // 15 bits

    mMuonTr1P = 0xFFFFFC00;      // 13 bits
    mMuonTrThetaX = 0xFFFFFF00;  // 15 bits
    mMuonTrThetaY = 0xFFFFFF00;  // 15 bits
    mMuonTrZmu = 0xFFFFFFF0;     // 19 bits
    mMuonTrBend = 0xFFFFFFF0;    // 19 bits
    mMuonTrNonBend = 0xFFFFFFF0; // 19 bits
    mMuonTrCov = 0xFFFF0000;     // 7 bits

    mMuonCl = 0xFFFFFF00;    // 15 bits
    mMuonClErr = 0xFFFF0000; // 7 bits

    mV0Time = 0xFFFFF000;      // 11 bits
    mADTime = 0xFFFFF000;      // 11 bits
    mT0Time = 0xFFFFFF00;      // 15 bits
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

  fOutputFile = TFile::Open("AO2D.root", "RECREATE", "O2 AOD", fCompress); // File to store the trees of time frames
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
  if (fSkipTPCPileup || fSkipPileup || fUseEventCuts)
    fEventCuts.AddQAplotsToList(fOutputList);
  if (fSkipTPCPileup)
    fEventCuts.SetRejectTPCPileupWithITSTPCnCluCorr(true);

  PostData(1, fOutputList);

  if (!fStoreHF) {
    DisableTree(kHF2Prong);
    DisableTree(kHF3Prong);
    DisableTree(kHFCascade);
    DisableTree(kHFDStar);
  }
} // void AliAnalysisTaskAO2Dconverter::UserCreateOutputObjects()

void AliAnalysisTaskAO2Dconverter::UserExec(Option_t *)
{
  // Initialisation

  const char *kPileupRejType[2] = {"PU_rej", "PU_TPC_rej"};

  fVEvent = InputEvent();
  fESD = dynamic_cast<AliESDEvent *>(fVEvent);
  fAOD = dynamic_cast<AliAODEvent *>(fVEvent);
  if (!fESD && !fAOD)
  {
    ::Fatal("AliAnalysisTaskAO2Dconverter::UserExec", "Something is wrong with the event handler");
  }

  // populate meta data only once
  if (fMetaData.GetEntries() == 0) {
    fMetaData.Add(new TObjString("DataType"), new TObjString((fTaskMode == kStandard) ? "RAW" : "MC"));
    fMetaData.Add(new TObjString("Run"), new TObjString("2"));
    TString converterVersion("aliroot ");
    converterVersion += ALIROOT_VERSION;
    converterVersion += ":";
    converterVersion += ALIROOT_REVISION;
    // TODO include above does not work?
    // converterVersion += "; aliphysics "
    // converterVersion += ALIPHYSICS_VERSION;
    // converterVersion += ":";
    // converterVersion += ALIPHYSICS_REVISION;
    converterVersion += "; root ";
    converterVersion += ROOT_RELEASE;
    fMetaData.Add(new TObjString("Run2ConverterVersion"), new TObjString(converterVersion));

    auto userInfo = fInputHandler->GetUserInfo();
    AliProdInfo prodInfo(userInfo);
    fMetaData.Add(new TObjString("ProducerAliRootVersion"), new TObjString(prodInfo.GetAlirootVersion()));
    fMetaData.Add(new TObjString("ProducerROOTVersion"), new TObjString(prodInfo.GetRootVersion()));
    fMetaData.Add(new TObjString("RecoPassName"), new TObjString(prodInfo.GetRecoPassName()));
    fMetaData.Add(new TObjString("AnchorProduction"), new TObjString(prodInfo.GetAnchorProduction()));
    fMetaData.Add(new TObjString("AnchorPassName"), new TObjString(prodInfo.GetAnchorPassName()));
    fMetaData.Add(new TObjString("LPMProductionTag"), new TObjString(prodInfo.GetTag(AliProdInfo::kProdTag)));
  }

  // This call is necessary to initialize event cuts according to the current run number
  bool alieventcut = fEventCuts.AcceptEvent(fVEvent);

  // In case of ESD we skip events like in the AOD filtering, for AOD this is not needed
  // We can use event cuts to avoid cases where we have zero reconstructed tracks
  bool skip_event = false;
  if (fESD && (fUseEventCuts || fSkipPileup || fSkipTPCPileup))
  {
    skip_event = !alieventcut && fUseEventCuts;
  }

  // Skip pileup events if requested
  if (fESD && fSkipPileup && !fEventCuts.PassedCut(AliEventCuts::kPileUp))
  {
    fHistPileupEvents->Fill(kPileupRejType[0], 1);
    skip_event = true;
  }

  // Check for TPC pileup events if requested, but don't skip event since it may affect physics
  if (fESD && fSkipTPCPileup && !fEventCuts.PassedCut(AliEventCuts::kTPCPileUp))
  {
    fHistPileupEvents->Fill(kPileupRejType[1], 1);
    //skip_event = true;
  }

  if (fTaskMode == kStandard)
    if (fESD && (fESD->GetHeader()->GetEventType() != 7)) // check for PHYSICS events
      skip_event = true;

  if (skip_event)
  {
    return;
  }

  if (!fTfInitialized)
  {
    ULong64_t tfId = GetGlobalBC(fVEvent->GetHeader());
    if (tfId == 0)
    {
      // The time stamp of the event is not set, for example in MC
      // Try the AliEn job ID
      TString alienPID(gSystem->Getenv("ALIEN_PROC_ID"));
      if (alienPID.Length() > 0)
      {
        tfId = alienPID.Atoll() * 1000 + fTFCount;
      }
      else
      {
        // Fallback:
        // Use the time period from 2020/11/01 in seconds (similar to what we did in the DAQ LDCs)
        TTimeStamp ts0(2020, 11, 1, 0, 0, 0);
        TTimeStamp ts1;
        tfId = ts1.GetSec() - ts0.GetSec();
        tfId *= 1000 + fTFCount;
      }
    }
    // Make globally unique TF id (across all data-taking)
    // The longest fill in Run 2 was 38 hours, which needs 43 bits. We reserve the values up to 1e13 which corresponds to 69 hours.
    // Run numbers in Run 2 were < 300k, which needs 19 bits
    // To make the number human-readable, we avoid a bit shift, but use a multiplication
    Int_t runNumber = fVEvent->GetRunNumber();
    if (runNumber > 0)
      tfId += (ULong64_t)(runNumber) * 10000000000000L;

    InitTF(tfId);
  }
  // Centrality QA
  Float_t centrality = -999;
  if (fESD) {
    // Get multiplicity selection
    AliMultSelection *multSelection = (AliMultSelection *)fESD->FindListObject("MultSelection");
    if (!multSelection)
      AliFatal("MultSelection not found in input event");

    centrality = multSelection->GetMultiplicityPercentile(fCentralityMethod);
  }
  else if (fAOD) {
    AliCentrality*  centralityObj = fAOD->GetCentrality();
    if (centralityObj) centrality = centralityObj->GetCentralityPercentile(fCentralityMethod);
  }

  // Fill centrality QA plots
  fCentralityHist->Fill(centrality);
  if ((fInputHandler->IsEventSelected() & AliVEvent::kINT7) != 0)
    fCentralityINT7->Fill(centrality);

  // Selection of events with at least two contributors (GMI)
  // Can this be done using the physics selection? (PH)

  const AliVVertex *pvtx = fVEvent->GetPrimaryVertex();
  if (!pvtx)
  {
    ::Fatal("AliAnalysisTaskAO2Dconverter::UserExec", "Vertex not defined");
  }

  // Now fill the content of the TF
  FillEventInTF();

  // Finish the current TF and initialize a new one, if the size is above the limit
  if (fBytes > fMaxBytes)
  {
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

  // write meta data
  fOutputFile->WriteObject(&fMetaData, "metaData");
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
  AliAnalysisTaskAO2Dconverter *task = new AliAnalysisTaskAO2Dconverter((TString("AO2D converter") + suffix).Data());
  if (!task)
    return nullptr;
  // add your task to the manager
  mgr->AddTask(task);
  // your task needs input: here we connect the manager to your task
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  // same for the output
  mgr->ConnectOutput(task, 1, mgr->CreateContainer("QAlist", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  // we need to register an unconnected output container to register the output file in the lego system
  mgr->CreateContainer("AO2D", TTree::Class(), AliAnalysisManager::kOutputContainer, "AO2D.root");
  // for (Int_t i = 0; i < kTrees; i++)
  //   mgr->ConnectOutput(task, 2 + i, mgr->CreateContainer(TreeName[i], TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  // in the end, this macro returns a pointer to your task. this will be convenient later on
  // when you will run your analysis in an analysis train on grid

  return task;
} // AliAnalysisTaskAO2Dconverter *AliAnalysisTaskAO2Dconverter::AddTask(TString suffix)

////////////////////////////////////////////////////////////
TTree *AliAnalysisTaskAO2Dconverter::CreateTree(TreeIndex t)
{
  if (!fTreeStatus[t])
    return 0x0;
  // Create the tree in the corresponding (TF) directory
  if (!fOutputDir)
    AliFatal("No Root subdir|");
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
  TObjArray *arr = fPruneList.Tokenize(" ");
  for (Int_t i = 0; i < arr->GetEntries(); i++)
  {
    Bool_t found = kFALSE;
    for (Int_t j = 0; j < kTrees; j++)
    {
      if (!fTreeStatus[j])
        continue;
      TObjArray *branches = fTree[j]->GetListOfBranches();
      for (Int_t k = 0; k < branches->GetEntries(); k++)
      {
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
  if (!fTreeStatus[t])
    return;
  Int_t nbytes = fTree[t]->Fill();
  if (nbytes > 0)
    fBytes += nbytes;
} // void AliAnalysisTaskAO2Dconverter::FillTree(TreeIndex t)

void AliAnalysisTaskAO2Dconverter::WriteTree(TreeIndex t)
{
  if (!fTreeStatus[t])
    return;
  // Write the tree in the corrsponding (TF) directory
  if (!fOutputDir)
    AliFatal("No Root subdir|");
  fOutputDir->cd();
  AliInfo(Form("Writing tree %s\n", TreeName[t].Data()));
  fTree[t]->Write();
} // void AliAnalysisTaskAO2Dconverter::WriteTree(TreeIndex t)

void AliAnalysisTaskAO2Dconverter::InitTF(ULong64_t tfId)
{
  Printf("Initializing TF %lld", tfId);

  // Reset the event count
  fTFCount++;
  fCollisionCount = 0;
  fBCCount = 0;
  fTfInitialized = true;

  // Reset the offsets
  fOffsetMuTrackID = 0;
  fOffsetTrack = 0;
  fOffsetV0 = 0;
  fOffsetLabel = 0;
  fOffsetHF2Prong = 0;

  // Reset the content of eventextra
  for (auto i = 0; i < kTrees; ++i)
  {
    eventextra.fStart[i] = 0;
    eventextra.fNentries[i] = 0;
  }

  // Create the output directory for the current time frame
  fOutputDir = fOutputFile->mkdir(Form("DF_%llu", tfId));

  // Associate branches for origin table
  TTree* tOrigin = CreateTree(kOrigin);
  if (fTreeStatus[kOrigin]) {
    tOrigin->Branch("fDataframeID", &origin.fDataframeID, "fDataframeID/l");
    origin.fDataframeID = tfId;
    FillTree(kOrigin);
  }

  // Associate branches for fEventTree
  TTree *tEvents = CreateTree(kEvents);
  if (fTreeStatus[kEvents])
  {
    tEvents->Branch("fIndexBCs", &collision.fIndexBCs, "fIndexBCs/I");
    tEvents->Branch("fPosX", &collision.fPosX, "fPosX/F");
    tEvents->Branch("fPosY", &collision.fPosY, "fPosY/F");
    tEvents->Branch("fPosZ", &collision.fPosZ, "fPosZ/F");
    tEvents->Branch("fCovXX", &collision.fCovXX, "fCovXX/F");
    tEvents->Branch("fCovXY", &collision.fCovXY, "fCovXY/F");
    tEvents->Branch("fCovXZ", &collision.fCovXZ, "fCovXZ/F");
    tEvents->Branch("fCovYY", &collision.fCovYY, "fCovYY/F");
    tEvents->Branch("fCovYZ", &collision.fCovYZ, "fCovYZ/F");
    tEvents->Branch("fCovZZ", &collision.fCovZZ, "fCovZZ/F");
    tEvents->Branch("fFlags", &collision.fFlags, "fFlags/s");
    tEvents->Branch("fChi2", &collision.fChi2, "fChi2/F");
    tEvents->Branch("fNumContrib", &collision.fN, "fNumContrib/s");
    tEvents->Branch("fCollisionTime", &collision.fCollisionTime, "fCollisionTime/F");
    tEvents->Branch("fCollisionTimeRes", &collision.fCollisionTimeRes, "fCollisionTimeRes/F");
    tEvents->SetBasketSize("*", fBasketSizeEvents);
  }

  // Extra information for debugging for event table
  TTree *tEventsExtra = CreateTree(kEventsExtra);
  if (fTreeStatus[kEventsExtra])
  {
    TString sstart = TString::Format("fStart[%d]/I", kTrees);
    TString sentries = TString::Format("fNentries[%d]/I", kTrees);
    tEventsExtra->Branch("fStart", eventextra.fStart, sstart.Data());
    tEventsExtra->Branch("fNentries", eventextra.fNentries, sentries.Data());
    tEventsExtra->SetBasketSize("*", fBasketSizeEvents);
  }

  // Associate branches for fEventTree
  TTree *tBC = CreateTree(kBC);
  if (fTreeStatus[kBC])
  {
    tBC->Branch("fRunNumber", &bc.fRunNumber, "fRunNumber/I");
    tBC->Branch("fGlobalBC", &bc.fGlobalBC, "fGlobalBC/l");
    tBC->Branch("fTriggerMask", &bc.fTriggerMask, "fTriggerMask/l");
    tBC->SetBasketSize("*", fBasketSizeEvents);
  }

  // Associate branches for Run 2 BC info
  TTree* tRun2BCInfo = CreateTree(kRun2BCInfo);
  if (fTreeStatus[kRun2BCInfo]) {
    tRun2BCInfo->Branch("fEventCuts", &run2bcinfo.fEventCuts, "fEventCuts/i");
    tRun2BCInfo->Branch("fTriggerMaskNext50", &run2bcinfo.fTriggerMaskNext50, "fTriggerMaskNext50/l");
    tRun2BCInfo->Branch("fL0TriggerInputMask", &run2bcinfo.fL0TriggerInputMask, "fL0TriggerInputMask/i");
    tRun2BCInfo->Branch("fSPDClustersL0", &run2bcinfo.fSPDClustersL0, "fSPDClustersL0/s");
    tRun2BCInfo->Branch("fSPDClustersL1", &run2bcinfo.fSPDClustersL1, "fSPDClustersL1/s");
    tRun2BCInfo->Branch("fSPDFiredChipsL0", &run2bcinfo.fSPDFiredChipsL0, "fSPDFiredChipsL0/s");
    tRun2BCInfo->Branch("fSPDFiredChipsL1", &run2bcinfo.fSPDFiredChipsL1, "fSPDFiredChipsL1/s");
    tRun2BCInfo->Branch("fSPDFiredFastOrL0", &run2bcinfo.fSPDFiredFastOrL0, "fSPDFiredFastOrL0/s");
    tRun2BCInfo->Branch("fSPDFiredFastOrL1", &run2bcinfo.fSPDFiredFastOrL1, "fSPDFiredFastOrL1/s");
    tRun2BCInfo->Branch("fV0TriggerChargeA", &run2bcinfo.fV0TriggerChargeA, "fV0TriggerChargeA/s");
    tRun2BCInfo->Branch("fV0TriggerChargeC", &run2bcinfo.fV0TriggerChargeC, "fV0TriggerChargeC/s");
    tRun2BCInfo->SetBasketSize("*", fBasketSizeEvents);
  }

  // Associate branches for the three track trees
  TTree* tTracks = CreateTree(kTracks);
  if (fTreeStatus[kTracks]) {
    tTracks->Branch("fIndexCollisions", &tracks.fIndexCollisions, "fIndexCollisions/I");
    tTracks->Branch("fTrackType", &tracks.fTrackType, "fTrackType/b");
    tTracks->Branch("fX", &tracks.fX, "fX/F");
    tTracks->Branch("fAlpha", &tracks.fAlpha, "fAlpha/F");
    tTracks->Branch("fY", &tracks.fY, "fY/F");
    tTracks->Branch("fZ", &tracks.fZ, "fZ/F");
    tTracks->Branch("fSnp", &tracks.fSnp, "fSnp/F");
    tTracks->Branch("fTgl", &tracks.fTgl, "fTgl/F");
    tTracks->Branch("fSigned1Pt", &tracks.fSigned1Pt, "fSigned1Pt/F");
    tTracks->SetBasketSize("*", fBasketSizeTracks);
  }

  TTree* tTracksCov = CreateTree(kTracksCov);
  if (fTreeStatus[kTracksCov]) {
    // Modified covariance matrix
    tTracksCov->Branch("fSigmaY", &tracks.fSigmaY, "fSigmaY/F");
    tTracksCov->Branch("fSigmaZ", &tracks.fSigmaZ, "fSigmaZ/F");
    tTracksCov->Branch("fSigmaSnp", &tracks.fSigmaSnp, "fSigmaSnp/F");
    tTracksCov->Branch("fSigmaTgl", &tracks.fSigmaTgl, "fSigmaTgl/F");
    tTracksCov->Branch("fSigma1Pt", &tracks.fSigma1Pt, "fSigma1Pt/F");
    tTracksCov->Branch("fRhoZY", &tracks.fRhoZY, "fRhoZY/B");
    tTracksCov->Branch("fRhoSnpY", &tracks.fRhoSnpY, "fRhoSnpY/B");
    tTracksCov->Branch("fRhoSnpZ", &tracks.fRhoSnpZ, "fRhoSnpZ/B");
    tTracksCov->Branch("fRhoTglY", &tracks.fRhoTglY, "fRhoTglY/B");
    tTracksCov->Branch("fRhoTglZ", &tracks.fRhoTglZ, "fRhoTglZ/B");
    tTracksCov->Branch("fRhoTglSnp", &tracks.fRhoTglSnp, "fRhoTglSnp/B");
    tTracksCov->Branch("fRho1PtY", &tracks.fRho1PtY, "fRho1PtY/B");
    tTracksCov->Branch("fRho1PtZ", &tracks.fRho1PtZ, "fRho1PtZ/B");
    tTracksCov->Branch("fRho1PtSnp", &tracks.fRho1PtSnp, "fRho1PtSnp/B");
    tTracksCov->Branch("fRho1PtTgl", &tracks.fRho1PtTgl, "fRho1PtTgl/B");
    tTracksCov->SetBasketSize("*", fBasketSizeTracks);
  }

  TTree* tTracksExtra = CreateTree(kTracksExtra);
  if (fTreeStatus[kTracksExtra]) {
    //Extra
    tTracksExtra->Branch("fTPCInnerParam", &tracks.fTPCinnerP, "fTPCInnerParam/F");
    tTracksExtra->Branch("fFlags", &tracks.fFlags, "fFlags/i");
    tTracksExtra->Branch("fITSClusterMap", &tracks.fITSClusterMap, "fITSClusterMap/b");
    tTracksExtra->Branch("fTPCNClsFindable", &tracks.fTPCNClsFindable, "fTPCNClsFindable/b");
    tTracksExtra->Branch("fTPCNClsFindableMinusFound",&tracks.fTPCNClsFindableMinusFound, "fTPCNClsFindableMinusFound/B");
    tTracksExtra->Branch("fTPCNClsFindableMinusCrossedRows", &tracks.fTPCNClsFindableMinusCrossedRows, "fTPCNClsFindableMinusCrossedRows/B");
    tTracksExtra->Branch("fTPCNClsShared", &tracks.fTPCNClsShared, "fTPCNClsShared/b");
    tTracksExtra->Branch("fTRDPattern", &tracks.fTRDPattern, "fTRDPattern/b");
    tTracksExtra->Branch("fITSChi2NCl", &tracks.fITSChi2NCl, "fITSChi2NCl/F");
    tTracksExtra->Branch("fTPCChi2NCl", &tracks.fTPCChi2NCl, "fTPCChi2NCl/F");
    tTracksExtra->Branch("fTRDChi2", &tracks.fTRDChi2, "fTRDChi2/F");
    tTracksExtra->Branch("fTOFChi2", &tracks.fTOFChi2, "fTOFChi2/F");
    tTracksExtra->Branch("fTPCSignal", &tracks.fTPCSignal, "fTPCSignal/F");
    tTracksExtra->Branch("fTRDSignal", &tracks.fTRDSignal, "fTRDSignal/F");
    // tTracksExtra->Branch("fTOFSignal", &tracks.fTOFSignal, "fTOFSignal/F");
    tTracksExtra->Branch("fLength", &tracks.fLength, "fLength/F");
    tTracksExtra->Branch("fTOFExpMom", &tracks.fTOFExpMom, "fTOFExpMom/F");
    tTracksExtra->Branch("fTrackEtaEMCAL", &tracks.fTrackEtaEMCAL, "fTrackEtaEMCAL/F");
    tTracksExtra->Branch("fTrackPhiEMCAL", &tracks.fTrackPhiEMCAL, "fTrackPhiEMCAL/F");
    tTracksExtra->Branch("fTrackTime", &tracks.fTrackTime, "fTrackTime/F");
    tTracksExtra->Branch("fTrackTimeRes", &tracks.fTrackTimeRes, "fTrackTimeRes/F");
    tTracksExtra->SetBasketSize("*", fBasketSizeTracks);
  }

  // Associate branches for Calo
  TTree *tCalo = CreateTree(kCalo);
  if (fTreeStatus[kCalo])
  {
    tCalo->Branch("fIndexBCs", &calo.fIndexBCs, "fIndexBCs/I");
    tCalo->Branch("fCellNumber", &calo.fCellNumber, "fCellNumber/S");
    tCalo->Branch("fAmplitude", &calo.fAmplitude, "fAmplitude/F");
    tCalo->Branch("fTime", &calo.fTime, "fTime/F");
    tCalo->Branch("fCellType", &calo.fCellType, "fCellType/B");
    tCalo->Branch("fCaloType", &calo.fCaloType, "fCaloType/B");
    tCalo->SetBasketSize("*", fBasketSizeEvents);
  }

  TTree *tCaloTrigger = CreateTree(kCaloTrigger);
  if (fTreeStatus[kCaloTrigger])
  {
    tCaloTrigger->Branch("fIndexBCs", &calotrigger.fIndexBCs, "fIndexBCs/I");
    tCaloTrigger->Branch("fFastOrAbsID", &calotrigger.fFastOrAbsID, "fFastOrAbsID/S");
    tCaloTrigger->Branch("fLnAmplitude", &calotrigger.fLnAmplitude, "fLnAmplitude/S");
    tCaloTrigger->Branch("fTriggerBits", &calotrigger.fTriggerBits, "fTriggerBits/I");
    tCaloTrigger->Branch("fCaloType", &calotrigger.fCaloType, "fCaloType/B");
    tCaloTrigger->SetBasketSize("*", fBasketSizeEvents);
  }

  // Associate FwdTrack branches for MUON tracks
  TTree *tFwdTrack = CreateTree(kFwdTrack);
  if (fTreeStatus[kFwdTrack])
  {
    tFwdTrack->Branch("fIndexCollisions", &fwdtracks.fIndexCollisions, "fIndexCollisions/I");
    tFwdTrack->Branch("fTrackType", &fwdtracks.fTrackType, "fTrackType/I");
    tFwdTrack->Branch("fX", &fwdtracks.fX, "fX/F");
    tFwdTrack->Branch("fY", &fwdtracks.fY, "fY/F");
    tFwdTrack->Branch("fZ", &fwdtracks.fZ, "fZ/F");
    tFwdTrack->Branch("fPhi", &fwdtracks.fPhi, "fPhi/F");
    tFwdTrack->Branch("fTgl", &fwdtracks.fTgl, "fTgl/F");
    tFwdTrack->Branch("fSigned1Pt", &fwdtracks.fSigned1Pt, "fSigned1Pt/F");
    tFwdTrack->Branch("fNClusters", &fwdtracks.fNClusters, "fNClusters/I");
    tFwdTrack->Branch("fPDca", &fwdtracks.fPDca, "fPDca/F");
    tFwdTrack->Branch("fRAtAbsorberEnd", &fwdtracks.fRAtAbsorberEnd, "fRAtAbsorberEnd/F");
    tFwdTrack->Branch("fChi2", &fwdtracks.fChi2, "fChi2/F");
    tFwdTrack->Branch("fChi2MatchMCHMID", &fwdtracks.fChi2MatchMCHMID, "fChi2MatchMCHMID/F");
    tFwdTrack->Branch("fChi2MatchMCHMFT", &fwdtracks.fChi2MatchMCHMFT, "fChi2MatchMCHMFT/F");
    tFwdTrack->Branch("fMatchScoreMCHMFT", &fwdtracks.fMatchScoreMCHMFT, "fMatchScoreMCHMFT/F");
    tFwdTrack->Branch("fTrackTime", &fwdtracks.fTrackTime, "fTrackTime/F");
    tFwdTrack->Branch("fTrackTimeRes", &fwdtracks.fTrackTimeRes, "fTrackTimeRes/F");
    tFwdTrack->Branch("fIndexMFTTracks", &fwdtracks.fIndexMFTTracks, "fIndexMFTTracks/I");
    tFwdTrack->Branch("fIndexFwdTracks_MatchMCHTrack", &fwdtracks.fIndexFwdTracks_MatchMCHTrack, "fIndexFwdTracks_MatchMCHTrack/I");
    tFwdTrack->Branch("fMCHBitMap", &fwdtracks.fMCHBitMap, "fMCHBitMap/s");
    tFwdTrack->Branch("fMIDBitMap", &fwdtracks.fMIDBitMap, "fMIDBitMap/s");
    tFwdTrack->Branch("fMIDBoards", &fwdtracks.fMIDBoards, "fMIDBoards/i");
    tFwdTrack->SetBasketSize("*", fBasketSizeEvents);
  }

  // Associate FwdTrack covariances for MUON tracks
  TTree *tFwdTrackCov = CreateTree(kFwdTrackCov);
  if (fTreeStatus[kFwdTrackCov])
  {
    tFwdTrackCov->Branch("fSigmaX", &fwdtracks.fSigmaX, "fSigmaX/F");
    tFwdTrackCov->Branch("fSigmaY", &fwdtracks.fSigmaY, "fSigmaY/F");
    tFwdTrackCov->Branch("fSigmaPhi", &fwdtracks.fSigmaPhi, "fSigmaPhi/F");
    tFwdTrackCov->Branch("fSigmaTgl", &fwdtracks.fSigmaTgl, "fSigmaTgl/F");
    tFwdTrackCov->Branch("fSigma1Pt", &fwdtracks.fSigma1Pt, "fSigma1Pt/F");
    tFwdTrackCov->Branch("fRhoXY", &fwdtracks.fRhoXY, "fRhoXY/B");
    tFwdTrackCov->Branch("fRhoPhiX", &fwdtracks.fRhoPhiX, "fRhoPhiX/B");
    tFwdTrackCov->Branch("fRhoPhiY", &fwdtracks.fRhoPhiY, "fRhoPhiY/B");
    tFwdTrackCov->Branch("fRhoTglX", &fwdtracks.fRhoTglX, "fRhoTglX/B");
    tFwdTrackCov->Branch("fRhoTglY", &fwdtracks.fRhoTglY, "fRhoTglY/B");
    tFwdTrackCov->Branch("fRhoTglPhi", &fwdtracks.fRhoTglPhi, "fRhoTglPhi/B");
    tFwdTrackCov->Branch("fRho1PtX", &fwdtracks.fRho1PtX, "fRho1PtX/B");
    tFwdTrackCov->Branch("fRho1PtY", &fwdtracks.fRho1PtY, "fRho1PtY/B");
    tFwdTrackCov->Branch("fRho1PtPhi", &fwdtracks.fRho1PtPhi, "fRho1PtPhi/B");
    tFwdTrackCov->Branch("fRho1PtTgl", &fwdtracks.fRho1PtTgl, "fRho1PtTgl/B");
    tFwdTrack->SetBasketSize("*", fBasketSizeEvents);
  }


  // Associuate branches for ZDC
  TTree *tZdc = CreateTree(kZdc);
  if (fTreeStatus[kZdc])
  {
    tZdc->Branch("fIndexBCs", &zdc.fIndexBCs, "fIndexBCs/I");
    tZdc->Branch("fEnergyZEM1", &zdc.fEnergyZEM1, "fEnergyZEM1/F");
    tZdc->Branch("fEnergyZEM2", &zdc.fEnergyZEM2, "fEnergyZEM2/F");
    tZdc->Branch("fEnergyCommonZNA", &zdc.fEnergyCommonZNA, "fEnergyCommonZNA/F");
    tZdc->Branch("fEnergyCommonZNC", &zdc.fEnergyCommonZNC, "fEnergyCommonZNC/F");
    tZdc->Branch("fEnergyCommonZPA", &zdc.fEnergyCommonZPA, "fEnergyCommonZPA/F");
    tZdc->Branch("fEnergyCommonZPC", &zdc.fEnergyCommonZPC, "fEnergyCommonZPC/F");
    tZdc->Branch("fEnergySectorZNA", &zdc.fEnergySectorZNA, "fEnergySectorZNA[4]/F");
    tZdc->Branch("fEnergySectorZNC", &zdc.fEnergySectorZNC, "fEnergySectorZNC[4]/F");
    tZdc->Branch("fEnergySectorZPA", &zdc.fEnergySectorZPA, "fEnergySectorZPA[4]/F");
    tZdc->Branch("fEnergySectorZPC", &zdc.fEnergySectorZPC, "fEnergySectorZPC[4]/F");
    tZdc->Branch("fTimeZEM1", &zdc.fTimeZEM1, "fTimeZEM1/F");
    tZdc->Branch("fTimeZEM2", &zdc.fTimeZEM2, "fTimeZEM2/F");
    tZdc->Branch("fTimeZNA", &zdc.fTimeZNA, "fTimeZNA/F");
    tZdc->Branch("fTimeZNC", &zdc.fTimeZNC, "fTimeZNC/F");
    tZdc->Branch("fTimeZPA", &zdc.fTimeZPA, "fTimeZPA/F");
    tZdc->Branch("fTimeZPC", &zdc.fTimeZPC, "fTimeZPC/F");
    tZdc->SetBasketSize("*", fBasketSizeEvents);
  }

  // Associate branches for V0A
  TTree *tFV0A = CreateTree(kFV0A);
  if (fTreeStatus[kFV0A])
  {
    tFV0A->Branch("fIndexBCs", &fv0a.fIndexBCs, "fIndexBCs/I");
    tFV0A->Branch("fChannel_size", &fv0a.fChannel_size, "fChannel_size/I");
    tFV0A->Branch("fChannel", fv0a.fChannel, "fChannel[fChannel_size]/b");
    tFV0A->Branch("fAmplitude_size", &fv0a.fAmplitude_size, "fAmplitude_size/I");
    tFV0A->Branch("fAmplitude", fv0a.fAmplitude, "fAmplitude[fAmplitude_size]/F");
    tFV0A->Branch("fTime", &fv0a.fTime, "fTime/F");
    tFV0A->Branch("fTriggerMask", &fv0a.fTriggerMask, "fTriggerMask/b");
    tFV0A->SetBasketSize("*", fBasketSizeEvents);
  }

  // Associate branches for V0C
  TTree *tFV0C = CreateTree(kFV0C);
  if (fTreeStatus[kFV0C])
  {
    tFV0C->Branch("fIndexBCs", &fv0c.fIndexBCs, "fIndexBCs/I");
    tFV0C->Branch("fChannel_size", &fv0c.fChannel_size, "fChannel_size/I");
    tFV0C->Branch("fChannel", fv0c.fChannel, "fChannel[fChannel_size]/b");
    tFV0C->Branch("fAmplitude_size", &fv0c.fAmplitude_size, "fAmplitude_size/I");
    tFV0C->Branch("fAmplitude", fv0c.fAmplitude, "fAmplitude[fAmplitude_size]/F");
    tFV0C->Branch("fTime", &fv0c.fTime, "fTime/F");
    tFV0C->SetBasketSize("*", fBasketSizeEvents);
  }

  // Associate branches for FT0
  TTree *tFT0 = CreateTree(kFT0);
  if (fTreeStatus[kFT0])
  {
    tFT0->Branch("fIndexBCs", &ft0.fIndexBCs, "fIndexBCs/I");
    tFT0->Branch("fChannelA_size", &ft0.fChannelA_size, "fChannelA_size/I");
    tFT0->Branch("fChannelA", ft0.fChannelA, "fChannelA[fChannelA_size]/b");
    tFT0->Branch("fAmplitudeA_size", &ft0.fAmplitudeA_size, "fAmplitudeA_size/I"); // will be removed with O2 improvement (one size field used for two VLAs)
    tFT0->Branch("fAmplitudeA", ft0.fAmplitudeA, "fAmplitudeA[fChannelA_size]/F");
    tFT0->Branch("fChannelC_size", &ft0.fChannelC_size, "fChannelC_size/I");
    tFT0->Branch("fChannelC", ft0.fChannelC, "fChannelC[fChannelC_size]/b");
    tFT0->Branch("fAmplitudeC_size", &ft0.fAmplitudeC_size, "fAmplitudeC_size/I"); // will be removed with O2 improvement (one size field used for two VLAs)
    tFT0->Branch("fAmplitudeC", ft0.fAmplitudeC, "fAmplitudeC[fChannelC_size]/F");
    tFT0->Branch("fTimeA", &ft0.fTimeA, "fTimeA/F");
    tFT0->Branch("fTimeC", &ft0.fTimeC, "fTimeC/F");
    tFT0->Branch("fTriggerMask", &ft0.fTriggerMask, "fTriggerMask/b");
    tFT0->SetBasketSize("*", fBasketSizeEvents);
  }

  // Associate branches for FDD (AD)
  TTree *tFDD = CreateTree(kFDD);
  if (fTreeStatus[kFDD])
  {
    tFDD->Branch("fIndexBCs", &fdd.fIndexBCs, "fIndexBCs/I");
    tFDD->Branch("fChargeA", fdd.fChargeA, "fChargeA[8]/I");
    tFDD->Branch("fChargeC", fdd.fChargeC, "fChargeC[8]/I");
    tFDD->Branch("fTimeA", &fdd.fTimeA, "fTimeA/F");
    tFDD->Branch("fTimeC", &fdd.fTimeC, "fTimeC/F");
    tFDD->Branch("fTriggerMask", &fdd.fTriggerMask, "fTriggerMask/b");
    tFDD->SetBasketSize("*", fBasketSizeEvents);
  }

  // Associuate branches for V0s
  TTree *tV0s = CreateTree(kV0s);
  if (fTreeStatus[kV0s])
  {
    tV0s->Branch("fIndexCollisions", &v0s.fIndexCollisions, "fIndexCollisions/I");
    tV0s->Branch("fIndexTracks_Pos", &v0s.fIndexTracksPos, "fIndexTracks_Pos/I");
    tV0s->Branch("fIndexTracks_Neg", &v0s.fIndexTracksNeg, "fIndexTracks_Neg/I");
    tV0s->SetBasketSize("*", fBasketSizeTracks);
  }

  // Associuate branches for cascades
  TTree *tCascades = CreateTree(kCascades);
  if (fTreeStatus[kCascades])
  {
    tCascades->Branch("fIndexCollisions", &cascs.fIndexCollisions, "fIndexCollisions/I");
    tCascades->Branch("fIndexV0s", &cascs.fIndexV0s, "fIndexV0s/I");
    tCascades->Branch("fIndexTracks", &cascs.fIndexTracks, "fIndexTracks/I");
    tCascades->SetBasketSize("*", fBasketSizeTracks);
  }

#ifdef USE_TOF_CLUST
  // Associate branches for TOF
  TTree *TOF = CreateTree(kTOF);
  if (fTreeStatus[kTOF])
  {
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

  // HMPID information
  TTree* tHMPID = CreateTree(kHMPID);
  if (fTreeStatus[kHMPID]) {
    tHMPID->Branch("fIndexTracks", &hmpids.fIndexTracks, "fIndexTracks/I");
    tHMPID->Branch("fHMPIDSignal", &hmpids.fHMPIDSignal, "fHMPIDSignal/F");
    tHMPID->Branch("fHMPIDDistance", &hmpids.fHMPIDDistance, "fHMPIDDistance/F");
    tHMPID->Branch("fHMPIDNPhotons", &hmpids.fHMPIDNPhotons, "fHMPIDNPhotons/S");
    tHMPID->Branch("fHMPIDQMip", &hmpids.fHMPIDQMip, "fHMPIDQMip/S");
    tHMPID->SetBasketSize("*", fBasketSizeTracks);
  }

  if (fTaskMode == kMC)
  {
    TTree *tMCvtx = CreateTree(kMcCollision);
    if (fTreeStatus[kMcCollision])
    {
      tMCvtx->Branch("fIndexBCs", &mccollision.fIndexBCs, "fIndexBCs/I");
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
    TTree *Kinematics = CreateTree(kMcParticle);
    if (fTreeStatus[kMcParticle])
    {
      Kinematics->Branch("fIndexMcCollisions", &mcparticle.fIndexMcCollisions, "fIndexMcCollisions/I");

      Kinematics->Branch("fPdgCode", &mcparticle.fPdgCode, "fPdgCode/I");
      Kinematics->Branch("fStatusCode", &mcparticle.fStatusCode, "fStatusCode/I");
      Kinematics->Branch("fFlags", &mcparticle.fFlags, "fFlags/b");

      Kinematics->Branch("fIndexArray_Mothers_size", &mcparticle.fIndexArray_Mothers_size, "fIndexArray_Mothers_size/I");
      Kinematics->Branch("fIndexArray_Mothers", &mcparticle.fIndexArray_Mothers, "fIndexArray_Mothers[fIndexArray_Mothers_size]/I");
      Kinematics->Branch("fIndexSlice_Daughters", &mcparticle.fIndexSlice_Daughters, "fIndexSlice_Daughters[2]/I");
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
    TTree *tLabels = CreateTree(kMcTrackLabel);
    if (fTreeStatus[kMcTrackLabel])
    {
      tLabels->Branch("fIndexMcParticles", &mctracklabel.fIndexMcParticles, "fIndexMcParticles/I");
      tLabels->Branch("fMcMask", &mctracklabel.fMcMask, "fMcMask/s");
      tLabels->SetBasketSize("*", fBasketSizeTracks);
    }

    // MC labels of each reconstructed calo cluster
    TTree *tCaloLabels = CreateTree(kMcCaloLabel);
    if (fTreeStatus[kMcCaloLabel])
    {
      tCaloLabels->Branch("fIndexMcParticles", &mccalolabel.fIndexMcParticles, "fIndexMcParticles/I");
      tCaloLabels->Branch("fMcMask", &mccalolabel.fMcMask, "fMcMask/s");
      tCaloLabels->SetBasketSize("*", fBasketSizeEvents);
    }

    // MC labels of each reconstructed calo cluster
    TTree *tCollisionLabels = CreateTree(kMcCollisionLabel);
    if (fTreeStatus[kMcCaloLabel])
    {
      tCollisionLabels->Branch("fIndexMcCollisions", &mccollisionlabel.fIndexMcCollisions, "fIndexMcCollisions/I");
      tCollisionLabels->Branch("fMcMask", &mccollisionlabel.fMcMask, "fMcMask/s");
      tCollisionLabels->SetBasketSize("*", fBasketSizeEvents);
    }

    // HepMC trees
    // Cross-sections
    TTree *tHepMcCrossSections = CreateTree(kHepMcCrossSections);
    if (fTreeStatus[kHepMcCrossSections])
    {
      tHepMcCrossSections->Branch("fIndexMcCollisions", &hepMcCrossSections.fIndexMcCollisions, "fIndexMcCollisions/I");
      // The generator may differ from the one in the MC collision or from the other HepMC tables in case of cocktail
      tHepMcCrossSections->Branch("fGeneratorsID", &hepMcCrossSections.fGeneratorsID, "fGeneratorsID/S");

      tHepMcCrossSections->Branch("fAccepted", &hepMcCrossSections.fAccepted, "fAccepted/l");
      tHepMcCrossSections->Branch("fAttempted", &hepMcCrossSections.fAttempted, "fAttempted/l");
      tHepMcCrossSections->Branch("fXsectGen", &hepMcCrossSections.fXsectGen, "fXsectGen/F");
      tHepMcCrossSections->Branch("fXsectErr", &hepMcCrossSections.fXsectErr, "fXsectErr/F");
      tHepMcCrossSections->Branch("fPtHard", &hepMcCrossSections.fPtHard, "fPtHard/F");
    }

    // PDF Information
    TTree *tHepMcPdfInfo = CreateTree(kHepMcPdfInfo);
    if (fTreeStatus[kHepMcPdfInfo])
    {
      tHepMcPdfInfo->Branch("fIndexMcCollisions", &hepMcPdfInfo.fIndexMcCollisions, "fIndexMcCollisions/I");
      // The generator may differ from the one in the MC collision or from the other HepMC tables in case of cocktail
      tHepMcPdfInfo->Branch("fGeneratorsID", &hepMcPdfInfo.fGeneratorsID, "fGeneratorsID/S");

      tHepMcPdfInfo->Branch("fId1", &hepMcPdfInfo.fId1, "fId1/I");
      tHepMcPdfInfo->Branch("fId2", &hepMcPdfInfo.fId2, "fId2/I");
      tHepMcPdfInfo->Branch("fPdfId1", &hepMcPdfInfo.fPdfId1, "fPdfId1/I");
      tHepMcPdfInfo->Branch("fPdfId2", &hepMcPdfInfo.fPdfId2, "fPdfId2/I");
      tHepMcPdfInfo->Branch("fX1", &hepMcPdfInfo.fX1, "fX1/F");
      tHepMcPdfInfo->Branch("fX2", &hepMcPdfInfo.fX2, "fX2/F");
      tHepMcPdfInfo->Branch("fScalePdf", &hepMcPdfInfo.fScalePdf, "fScalePdf/F");
      tHepMcPdfInfo->Branch("fPdf1", &hepMcPdfInfo.fPdf1, "fPdf1/F");
      tHepMcPdfInfo->Branch("fPdf2", &hepMcPdfInfo.fPdf2, "fPdf2/F");
    }

    // Heavy Ion information
    TTree *tHepMcHeavyIon = CreateTree(kHepMcHeavyIon);
    if (fTreeStatus[kHepMcHeavyIon])
    {
      tHepMcHeavyIon->Branch("fIndexMcCollisions", &hepMcHeavyIon.fIndexMcCollisions, "fIndexMcCollisions/I");
      // The generator may differ from the one in the MC collision or from the other HepMC tables in case of cocktail
      tHepMcHeavyIon->Branch("fGeneratorsID", &hepMcHeavyIon.fGeneratorsID, "fGeneratorsID/S");

      tHepMcHeavyIon->Branch("fNcollHard", &hepMcHeavyIon.fNcollHard, "fNcollHard/I");
      tHepMcHeavyIon->Branch("fNpartProj", &hepMcHeavyIon.fNpartProj, "fNpartProj/I");
      tHepMcHeavyIon->Branch("fNpartTarg", &hepMcHeavyIon.fNpartTarg, "fNpartTarg/I");
      tHepMcHeavyIon->Branch("fNcoll", &hepMcHeavyIon.fNcoll, "fNcoll/I");
      tHepMcHeavyIon->Branch("fNNwoundedCollisions", &hepMcHeavyIon.fNNwoundedCollisions, "fNNwoundedCollisions/I");
      tHepMcHeavyIon->Branch("fNwoundedNCollisions", &hepMcHeavyIon.fNwoundedNCollisions, "fNwoundedNCollisions/I");
      tHepMcHeavyIon->Branch("fNwoundedNwoundedCollisions", &hepMcHeavyIon.fNwoundedNwoundedCollisions, "fNwoundedNwoundedCollisions/I");
      tHepMcHeavyIon->Branch("fSpectatorNeutrons", &hepMcHeavyIon.fSpectatorNeutrons, "fSpectatorNeutrons/I");
      tHepMcHeavyIon->Branch("fSpectatorProtons", &hepMcHeavyIon.fSpectatorProtons, "fSpectatorProtons/I");
      tHepMcHeavyIon->Branch("fImpactParameter", &hepMcHeavyIon.fImpactParameter, "fImpactParameter/F");
      tHepMcHeavyIon->Branch("fEventPlaneAngle", &hepMcHeavyIon.fEventPlaneAngle, "fEventPlaneAngle/F");
      tHepMcHeavyIon->Branch("fEccentricity", &hepMcHeavyIon.fEccentricity, "fEccentricity/F");
      tHepMcHeavyIon->Branch("fSigmaInelNN", &hepMcHeavyIon.fSigmaInelNN, "fSigmaInelNN/F");
      tHepMcHeavyIon->Branch("fCentrality", &hepMcHeavyIon.fCentrality, "fCentrality/F");
    }
  } // End if(fTaskMode == kMC)

  if (fStoreHF) {
    // HF 2 prong candidates
    TTree *tHF2Prong = CreateTree(kHF2Prong);
    if (fTreeStatus[kHF2Prong])
    {
      tHF2Prong->Branch("fIndexTracks_0", &hf2Prong.fIndexTracks_0, "fIndexTracks_0/I");
      tHF2Prong->Branch("fIndexTracks_1", &hf2Prong.fIndexTracks_1, "fIndexTracks_1/I");
      tHF2Prong->Branch("fHFflag", &hf2Prong.fHFflag, "fHFflag/b");
      tHF2Prong->SetBasketSize("*", fBasketSizeEvents);
    }

    TTree *tHF3Prong = CreateTree(kHF3Prong);
    if (fTreeStatus[kHF3Prong])
    {
      tHF3Prong->Branch("fIndexTracks_0", &hf3Prong.fIndexTracks_0, "fIndexTracks_0/I");
      tHF3Prong->Branch("fIndexTracks_1", &hf3Prong.fIndexTracks_1, "fIndexTracks_1/I");
      tHF3Prong->Branch("fIndexTracks_2", &hf3Prong.fIndexTracks_2, "fIndexTracks_2/I");
      tHF3Prong->Branch("fHFflag", &hf3Prong.fHFflag, "fHFflag/b");
      tHF3Prong->SetBasketSize("*", fBasketSizeEvents);
    }

    TTree *tHFDStar = CreateTree(kHFDStar);
    if (fTreeStatus[kHFDStar])
    {
      tHFDStar->Branch("fIndexTracks_0", &hfDStar.fIndexTracks_0, "fIndexTracks_0/I");
      tHFDStar->Branch("fIndexHf2Prongs", &hfDStar.fIndexHf2Prongs, "fIndexHf2Prongs/I");
      tHFDStar->SetBasketSize("*", fBasketSizeEvents);
    }

    TTree *tHFCascade = CreateTree(kHFCascade);
    if (fTreeStatus[kHFCascade])
    {
      tHFCascade->Branch("fIndexV0s", &hfCascades.fIndexV0s, "fIndexV0s/I");
      tHFCascade->Branch("fIndexTracks_0", &hfCascades.fIndexTracks_0, "fIndexTracks_0/I");
      tHFCascade->SetBasketSize("*", fBasketSizeEvents);
    }
  }

  Prune(); //Removing all unwanted branches (if any)
} // void AliAnalysisTaskAO2Dconverter::InitTF(Int_t tfId)

void AliAnalysisTaskAO2Dconverter::FillEventInTF()
{
  // Configuration of the PID response
  AliPIDResponse *PIDResponse = (AliPIDResponse *)((AliInputEventHandler *)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetPIDResponse();
  if (!PIDResponse) {
    ((AliInputEventHandler *)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->CreatePIDResponse(fTaskMode==kMC);
    PIDResponse = (AliPIDResponse *)((AliInputEventHandler *)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetPIDResponse();
  }
  PIDResponse->SetTOFResponse(fVEvent, AliPIDResponse::kBest_T0);
  AliTOFPIDResponse &TOFResponse = PIDResponse->GetTOFResponse();

  // Configuration of the MC event (if needed)
  AliMCEvent *MCEvt = nullptr;
  // In case of AOD the access is via TClonesArray containing AliAODMCParticles, and AliAODMCHeader
  TClonesArray *MCArray = nullptr;
  AliAODMCHeader *MCHeader = nullptr;

  if (fTaskMode == kMC)
  {
    if (fESD) {
      AliMCEventHandler *eventHandler = dynamic_cast<AliMCEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()); //Get the MC handler

      if (!eventHandler) //Check on the MC handler
	AliFatal("Could not retrieve MC event handler");
      MCEvt = eventHandler->MCEvent(); //Get the MC Event

      if (!MCEvt) // Check on the MC Event
	AliFatal("Could not retrieve MC event");
      PIDResponse->SetCurrentMCEvent(MCEvt); //Set The PID response on the current MC event
    }
    else if (fAOD) {
      MCArray = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
      MCHeader = dynamic_cast<AliAODMCHeader*>(fAOD->FindListObject(AliAODMCHeader::StdBranchName()));
    }
  }

  // Adjust start indices for this event in all trees by adding the number of entries of the previous event
  for (auto i = 0; i < kTrees; ++i) {
    eventextra.fStart[i] += eventextra.fNentries[i];
    eventextra.fNentries[i] = 0;
  }

  // Decide if vertex should be written
  // Primary vertex
  const AliVVertex *pvtx = fVEvent->GetPrimaryVertex();
  UChar_t vertexType = 0;

  TString title(pvtx->GetTitle());
  if (pvtx->IsFromVertexer3D()) {
    vertexType = Run2Vertexer3D;
  } else if (pvtx->IsFromVertexerZ()) {
    vertexType = Run2VertexerZ;
  } else if (title.Contains("VertexerTracks")) {
    if (pvtx->GetNContributors() >= 2) {
      vertexType = Run2VertexerTracks;
      if (title.Contains("WithConstraint"))
        vertexType |= Run2VertexerTracksWithConstraint;
      if (title.Contains("OnlyFitter"))
        vertexType |= Run2VertexerTracksOnlyFitter;
      if (title.Contains("MV"))
        vertexType |= Run2VertexerTracksMultiVertex;
    }
  } else {
    AliWarning(Form("Unknown vertexer type: %s", title.Data()));
  }

  Bool_t fillCollision = kFALSE;
  if (vertexType > 0) {
    fillCollision = kTRUE;
    AliInfo(Form("Event with vertex %s and %d contributors: accepted as vertex type %d", title.Data(), pvtx->GetNContributors(), vertexType));
  } else {
    AliInfo(Form("Event with vertex %s and %d contributors: rejected", title.Data(), pvtx->GetNContributors()));
  }

  //---------------------------------------------------------------------------
  // Collision data
  if (fillCollision)
  {
    eventextra.fNentries[kEvents] = 1; // one entry per vertex
    collision.fIndexBCs = fBCCount;
    collision.fPosX = AliMathBase::TruncateFloatFraction(pvtx->GetX(), mCollisionPosition);
    collision.fPosY = AliMathBase::TruncateFloatFraction(pvtx->GetY(), mCollisionPosition);
    collision.fPosZ = AliMathBase::TruncateFloatFraction(pvtx->GetZ(), mCollisionPosition);

    Double_t covmatrix[6];
    pvtx->GetCovarianceMatrix(covmatrix);

    collision.fCovXX = AliMathBase::TruncateFloatFraction(covmatrix[0], mCollisionPositionCov);
    collision.fCovXY = AliMathBase::TruncateFloatFraction(covmatrix[1], mCollisionPositionCov);
    collision.fCovXZ = AliMathBase::TruncateFloatFraction(covmatrix[3], mCollisionPositionCov);
    collision.fCovYY = AliMathBase::TruncateFloatFraction(covmatrix[2], mCollisionPositionCov);
    collision.fCovYZ = AliMathBase::TruncateFloatFraction(covmatrix[4], mCollisionPositionCov);
    collision.fCovZZ = AliMathBase::TruncateFloatFraction(covmatrix[5], mCollisionPositionCov);

    collision.fFlags = vertexType;
    collision.fChi2 = AliMathBase::TruncateFloatFraction(pvtx->GetChi2(), mCollisionPositionCov);
    collision.fN = (pvtx->GetNContributors() > USHRT_MAX) ? USHRT_MAX : pvtx->GetNContributors();

    Float_t eventTime[10];
    Float_t eventTimeRes[10];
    Double_t eventTimeWeight[10];

    for (Int_t i = 0; i < TOFResponse.GetNmomBins(); i++)
    {
      if (i >= 10)
        AliFatal("Index is too high!");
      Float_t mom = (TOFResponse.GetMinMom(i) + TOFResponse.GetMaxMom(i)) / 2.f;
      eventTime[i] = TOFResponse.GetStartTime(mom) * 1.e-3;       // ps to ns
      eventTimeRes[i] = TOFResponse.GetStartTimeRes(mom) * 1.e-3; // ps to ns
      eventTimeWeight[i] = 1. / (eventTimeRes[i] * eventTimeRes[i]);
    }

    // Recalculate unique event time and its resolution
    collision.fCollisionTime = AliMathBase::TruncateFloatFraction(TMath::Mean(10, eventTime, eventTimeWeight), mCollisionPosition);                 // Weighted mean of times per momentum interval
    collision.fCollisionTimeRes = AliMathBase::TruncateFloatFraction(TMath::Sqrt(9. / 10.) * TMath::Mean(10, eventTimeRes), mCollisionPositionCov); // PH bad approximation

    FillTree(kEvents);
  }

  //---------------------------------------------------------------------------
  // BC data

  bc.fRunNumber = fVEvent->GetRunNumber();

  ULong64_t evtid = GetGlobalBC(fVEvent->GetHeader());
  bc.fGlobalBC = evtid;
  bc.fTriggerMask = fVEvent->GetTriggerMask();
  // NOTE upper 64 bit of trigger classes stored few lines below in run2bcinfo.fTriggerMaskNext50
  FillTree(kBC);

  //---------------------------------------------------------------------------
  // Run 2 BC information
  run2bcinfo.fTriggerMaskNext50 = fVEvent->GetTriggerMaskNext50();
  run2bcinfo.fL0TriggerInputMask = fVEvent->GetHeader()->GetL0TriggerInputs();
  run2bcinfo.fSPDClustersL0 = (fVEvent->GetNumberOfITSClusters(0) > USHRT_MAX) ? USHRT_MAX : fVEvent->GetNumberOfITSClusters(0);
  run2bcinfo.fSPDClustersL1 = (fVEvent->GetNumberOfITSClusters(1) > USHRT_MAX) ? USHRT_MAX : fVEvent->GetNumberOfITSClusters(1);
  const TBits& onMap = fInputEvent->GetMultiplicity()->GetFastOrFiredChips();
  const TBits& ofMap = fInputEvent->GetMultiplicity()->GetFiredChipMap();
  run2bcinfo.fSPDFiredFastOrL1 = onMap.CountBits(400);
  run2bcinfo.fSPDFiredChipsL1 = ofMap.CountBits(400);
  run2bcinfo.fSPDFiredFastOrL0 = onMap.CountBits(0)-run2bcinfo.fSPDFiredFastOrL1;
  run2bcinfo.fSPDFiredChipsL0 = ofMap.CountBits(0)-run2bcinfo.fSPDFiredChipsL1;
  AliVVZERO* vzero = fInputEvent->GetVZEROData();
  run2bcinfo.fV0TriggerChargeA = vzero->GetTriggerChargeA();
  run2bcinfo.fV0TriggerChargeC = vzero->GetTriggerChargeC();

  run2bcinfo.fEventCuts = 0;

  if (fESD) {
    // Get multiplicity selection
    AliMultSelection *multSelection = (AliMultSelection*) fESD->FindListObject("MultSelection");
    if (!multSelection)
      AliFatal("MultSelection not found in input event");

    if( multSelection->GetThisEventINELgtZERO() )
      SETBIT (run2bcinfo.fEventCuts, kINELgtZERO);

    if( multSelection->GetThisEventIsNotPileupInMultBins() )
      SETBIT (run2bcinfo.fEventCuts, kPileupInMultBins);

    if( multSelection->GetThisEventHasNoInconsistentVertices() )
      SETBIT (run2bcinfo.fEventCuts, kConsistencySPDandTrackVertices);

    if( multSelection->GetThisEventPassesTrackletVsCluster() )
      SETBIT (run2bcinfo.fEventCuts, kTrackletsVsClusters);

    if( fESD->GetPrimaryVertex()->GetNContributors()>0 )
      SETBIT (run2bcinfo.fEventCuts, kNonZeroNContribs);

    if( multSelection->GetThisEventIsNotIncompleteDAQ() )
      SETBIT (run2bcinfo.fEventCuts, kIncompleteDAQ);

    if (fEventCuts.PassedCut(AliEventCuts::kPileUp))
      SETBIT(run2bcinfo.fEventCuts, kPileUpMV);

    if (fEventCuts.PassedCut(AliEventCuts::kTPCPileUp))
      SETBIT(run2bcinfo.fEventCuts, kTPCPileUp);

    if (fEventCuts.PassedCut(AliEventCuts::kTimeRangeCut))
      SETBIT(run2bcinfo.fEventCuts, kTimeRangeCut);

    if (fEventCuts.PassedCut(AliEventCuts::kEMCALEDCut))
      SETBIT(run2bcinfo.fEventCuts, kEMCALEDCut);

    if (fEventCuts.PassedCut(AliEventCuts::kAllCuts))
      SETBIT(run2bcinfo.fEventCuts, kAliEventCutsAccepted);

    if (fUseTriggerAnalysis) {
      if (fTriggerAnalysis.IsSPDVtxPileup(fInputEvent))
        SETBIT(run2bcinfo.fEventCuts, kIsPileupFromSPD);

      if (fTriggerAnalysis.IsV0PFPileup(fInputEvent))
        SETBIT(run2bcinfo.fEventCuts, kIsV0PFPileup);

      if (fTriggerAnalysis.IsHVdipTPCEvent(fInputEvent))
        SETBIT(run2bcinfo.fEventCuts, kIsTPCHVdip);

      if (fTriggerAnalysis.IsLaserWarmUpTPCEvent(fInputEvent))
        SETBIT(run2bcinfo.fEventCuts, kIsTPCLaserWarmUp);

      if (fTriggerAnalysis.TRDTrigger(fInputEvent,AliTriggerAnalysis::kTRDHCO))
        SETBIT(run2bcinfo.fEventCuts, kTRDHCO);

      if (fTriggerAnalysis.TRDTrigger(fInputEvent,AliTriggerAnalysis::kTRDHJT))
        SETBIT(run2bcinfo.fEventCuts, kTRDHJT);

      if (fTriggerAnalysis.TRDTrigger(fInputEvent,AliTriggerAnalysis::kTRDHSE))
        SETBIT(run2bcinfo.fEventCuts, kTRDHSE);

      if (fTriggerAnalysis.TRDTrigger(fInputEvent,AliTriggerAnalysis::kTRDHQU))
        SETBIT(run2bcinfo.fEventCuts, kTRDHQU);

      if (fTriggerAnalysis.TRDTrigger(fInputEvent,AliTriggerAnalysis::kTRDHEE))
        SETBIT(run2bcinfo.fEventCuts, kTRDHEE);
    }
  }
  else {
    //PH AOD case: What should we use here?
  }

  FillTree(kRun2BCInfo);

  //---------------------------------------------------------------------------
  // MC kinematic tree
  // It has to be done before the reconstructed tracks,
  // since the cleaning requires reasignment of the labels
  Int_t nkine_filled = 0; // Number of kine tracks filled

  TArrayC toWrite;
  TArrayI kineIndex;
  if (MCEvt || MCArray)
  {
    // Kinematics
    TParticle *particle = nullptr;
    Int_t nMCtracks = MCEvt? MCEvt->GetNumberOfTracks() : MCArray->GetEntriesFast();

    toWrite.Reset();
    toWrite.Set(nMCtracks);
    kineIndex.Reset();
    kineIndex.Set(nMCtracks);

    Int_t nMCprim = 0;
    if (MCEvt) {
      nMCprim = MCEvt->GetNumberOfPrimaries();
    }
    else {
      for (Int_t it = 0; it < nMCtracks; it++) {
        AliAODMCParticle* mcAODTrack = dynamic_cast<AliAODMCParticle*>(MCArray->At(it));
        if (!mcAODTrack) {
	  Error("UserExec", "Could not receive MC track %d", it);
	  continue;
        }
	if (mcAODTrack->IsPrimary()) nMCprim++;
      }
    }

    // The selection below is needed only for ESD. In case of AOD we should just copy the content of the array
    if (MCEvt) {
      // For each reconstructed track keep the corresponding MC particle
      Int_t ntracks = fVEvent->GetNumberOfTracks();

      for (Int_t itrack = 0; itrack < ntracks; ++itrack)
      {
	AliVTrack *track = dynamic_cast<AliVTrack*>(fVEvent->GetTrack(itrack));
	Int_t alabel = track->GetLabel();
	toWrite[TMath::Abs(alabel)] = 1;
      }

      // For each calo cluster keep the corresponding MC particle
      AliVCaloCells *emcalCells = fVEvent->GetEMCALCells();
      Short_t nEmcalCells = emcalCells->GetNumberOfCells();
      for (Short_t ice = 0; ice < nEmcalCells; ++ice)
      {
	Short_t cellNumber;
	Double_t amplitude;
	Double_t time;
	Int_t mclabel;
	Double_t efrac;

	emcalCells->GetCell(ice, cellNumber, amplitude, time, mclabel, efrac);
	toWrite[TMath::Abs(mclabel)] = 1;
      }
      AliVCaloCells *phosCells = fVEvent->GetPHOSCells();
      Short_t nPhosCells = phosCells->GetNumberOfCells();
      for (Short_t ice = 0; ice < nPhosCells; ++ice)
      {
	Short_t cellNumber;
	Double_t amplitude;
	Double_t time;
	Int_t mclabel;
	Double_t efrac;

	phosCells->GetCell(ice, cellNumber, amplitude, time, mclabel, efrac);
        if(mclabel>=0) //label -1 means no primary
          toWrite[mclabel] = 1;
      }

      // For each tracklet keep the corresponding MC particle
      AliVMultiplicity *mult = fVEvent->GetMultiplicity();
      Int_t ntracklets = mult->GetNumberOfTracklets();

      for (Int_t itr = ntracklets; itr--;)
      {
	Int_t alabel = mult->GetLabel(itr, 0); // Take the label of the first layer
	toWrite[TMath::Abs(alabel)] = 1;
      }

      // For each MUON track keep the corresponding MC particle
      Int_t nmuons = fVEvent->GetNumberOfMuonTracks();

      for (Int_t imu = 0; imu < nmuons; ++imu)
      {
	Int_t alabel = 0;
	if (fESD) {
	  AliESDMuonTrack *mutrk = fESD->GetMuonTrack(imu);
	  alabel = mutrk->GetLabel();
	  toWrite[TMath::Abs(alabel)] = 1;
	}
      }
    }

    // Define which MC particles pass kinematic selection and update toWrite index
    for (Int_t i = 0; i < nMCtracks; ++i)
    {
      AliVParticle *vpt = MCEvt ? MCEvt->GetTrack(i) : nullptr;
      AliAODMCParticle * aodmcpt = MCArray ? dynamic_cast<AliAODMCParticle*>(MCArray->At(i)) : nullptr;
      particle = vpt ? vpt->Particle() : nullptr;
      //PH In case of AliAODMCParticle the method Particle id not reimplemented and returns 0

      Float_t xv = particle ? particle->Vx() : aodmcpt->Xv();
      Float_t yv = particle ? particle->Vy() : aodmcpt->Yv();
      Float_t zv = particle ? particle->Vz() : aodmcpt->Zv();
      Float_t rv = TMath::Sqrt(xv * xv + yv * yv);

      // This part is taken from AliAnalysisTaskMCParticleFilter
      // If AOD, write all MC particles, esle make selection
      Bool_t write = aodmcpt ? kTRUE : kFALSE;

      if (MCEvt) {
	if (i < nMCprim)
        {
	  // Select all primary particles
	  write = kTRUE;
	}
	else if (particle->GetUniqueID() == kPDecay)
        {
	  // Particles from decay
	  // Check that the decay chain ends at a primary particle
	  AliVParticle *mother = vpt;
	  Int_t imo = vpt->GetMother();
	  while ((imo >= nMCprim) && (mother->Particle()->GetUniqueID() == kPDecay))
          {
	    mother = MCEvt ? MCEvt->GetTrack(imo) : dynamic_cast<AliAODMCParticle*>(MCArray->At(imo));
	    imo = mother->GetMother();
         }
	  // Select according to pseudorapidity and production point of primary ancestor
	  if (imo < nMCprim)
	    write = kTRUE;
	}
	else if (particle->GetUniqueID() == kPPair)
        {
	  // Now look for pair production
	  Int_t imo = vpt->GetMother();
	  if (imo < nMCprim)
	  {
	    // Select, if the gamma is a primary
	    write = kTRUE;
	  }
	  else
          {
	    // Check if the gamma comes from the decay chain of a primary particle
	    AliVParticle *mother = MCEvt ? MCEvt->GetTrack(imo) : dynamic_cast<AliAODMCParticle*>(MCArray->At(imo));
	    imo = mother->GetMother();
	    while ((imo >= nMCprim) && (mother->Particle()->GetUniqueID() == kPDecay))
	    {
	      mother = MCEvt ? MCEvt->GetTrack(imo) : dynamic_cast<AliAODMCParticle*>(MCArray->At(imo));
	      imo = mother->GetMother();
	    }
	    // Select according to pseudorapidity and production point
	    if (imo < nMCprim && Select(mother->Particle(), rv, zv))
	      write = kTRUE;
	  }
	}
      }
      if (toWrite[i] > 0 || write)
      {
        toWrite[i] = 1;
        kineIndex[i] = nkine_filled;
        nkine_filled++;
      }
      else
      {
        kineIndex[i] = -1;
      }
    }

    for (Int_t i = 0; i < nMCtracks; ++i)
    { //loop on primary MC tracks Before Event Selection
      AliVParticle *vpt = MCEvt ? MCEvt->GetTrack(i) : nullptr;
      AliAODMCParticle * aodmcpt = MCArray ? dynamic_cast<AliAODMCParticle*>(MCArray->At(i)) : nullptr;
      particle = vpt ? vpt->Particle() : nullptr;

      mcparticle.fIndexMcCollisions = fBCCount;

      //Get the kinematic values of the particles
      mcparticle.fPdgCode = particle ? particle->GetPdgCode() : aodmcpt->GetPdgCode();
      mcparticle.fStatusCode = particle ? particle->GetStatusCode() : aodmcpt->MCStatusCode(); //PH possible issue
      mcparticle.fFlags = 0;
      //PH The part below probably needs modifications
      Int_t nPrim = MCEvt ? MCEvt->Stack()->GetNprimary() : nMCprim;
      if (i >= nPrim)
        mcparticle.fFlags |= MCParticleFlags::ProducedInTransport;

      if (MCEvt && MCEvt->IsPhysicalPrimary(i)) // ESD
        mcparticle.fFlags |= MCParticleFlags::PhysicalPrimary;
      if (aodmcpt && aodmcpt->IsPhysicalPrimary()) // AOD
        mcparticle.fFlags |= MCParticleFlags::PhysicalPrimary;

      if (MCEvt && MCEvt->IsFromBGEvent(i)) // ESD
        mcparticle.fFlags |= MCParticleFlags::FromBackgroundEvent;
      if (aodmcpt && aodmcpt->IsFromSubsidiaryEvent()) // AOD: PH Not sure if corresponds to BKG event
        mcparticle.fFlags |= MCParticleFlags::FromBackgroundEvent;

      mcparticle.fIndexArray_Mothers_size = 1;
      mcparticle.fIndexArray_Mothers[0] = particle ? particle->GetMother(0) : aodmcpt->GetMother();
      if (mcparticle.fIndexArray_Mothers[0] > -1)
        mcparticle.fIndexArray_Mothers[0] = kineIndex[mcparticle.fIndexArray_Mothers[0]] > -1 ? kineIndex[mcparticle.fIndexArray_Mothers[0]] + fOffsetLabel : -1;
      if (mcparticle.fIndexArray_Mothers[0] == -1)
        mcparticle.fIndexArray_Mothers_size = 0;

      mcparticle.fIndexSlice_Daughters[0] = particle ? particle->GetFirstDaughter() : aodmcpt->GetDaughterFirst();
      if (mcparticle.fIndexSlice_Daughters[0] > -1)
        mcparticle.fIndexSlice_Daughters[0] = kineIndex[mcparticle.fIndexSlice_Daughters[0]] > -1 ? kineIndex[mcparticle.fIndexSlice_Daughters[0]] + fOffsetLabel : -1;
      mcparticle.fIndexSlice_Daughters[1] = particle ? particle->GetLastDaughter() : aodmcpt->GetDaughterLast();
      if (mcparticle.fIndexSlice_Daughters[1] > -1)
        mcparticle.fIndexSlice_Daughters[1] = kineIndex[mcparticle.fIndexSlice_Daughters[1]] > -1 ? kineIndex[mcparticle.fIndexSlice_Daughters[1]] + fOffsetLabel : -1;
      if (mcparticle.fIndexSlice_Daughters[0] > -1 && mcparticle.fIndexSlice_Daughters[1] == -1)
        mcparticle.fIndexSlice_Daughters[1] = mcparticle.fIndexSlice_Daughters[0];
      if (mcparticle.fIndexSlice_Daughters[1] > -1 && mcparticle.fIndexSlice_Daughters[0] == -1)
        mcparticle.fIndexSlice_Daughters[0] = mcparticle.fIndexSlice_Daughters[1];
      mcparticle.fWeight = AliMathBase::TruncateFloatFraction(particle ? particle->GetWeight() : 1., mMcParticleW);

      mcparticle.fPx = AliMathBase::TruncateFloatFraction(particle ? particle->Px() : aodmcpt->Px(), mMcParticleMom);
      mcparticle.fPy = AliMathBase::TruncateFloatFraction(particle ? particle->Py() : aodmcpt->Py(), mMcParticleMom);
      mcparticle.fPz = AliMathBase::TruncateFloatFraction(particle ? particle->Pz() : aodmcpt->Pz(), mMcParticleMom);
      mcparticle.fE = AliMathBase::TruncateFloatFraction(particle ? particle->Energy() : aodmcpt->E(), mMcParticleMom);

      mcparticle.fVx = AliMathBase::TruncateFloatFraction(particle ? particle->Vx() : aodmcpt->Xv(), mMcParticlePos);
      mcparticle.fVy = AliMathBase::TruncateFloatFraction(particle ? particle->Vy() : aodmcpt->Yv(), mMcParticlePos);
      mcparticle.fVz = AliMathBase::TruncateFloatFraction(particle ? particle->Vz() : aodmcpt->Zv(), mMcParticlePos);
      mcparticle.fVt = AliMathBase::TruncateFloatFraction(particle ? particle->T() : aodmcpt->T(), mMcParticlePos);

      if (toWrite[i] > 0)
      {
        FillTree(kMcParticle);
      }
    }
  }
  eventextra.fNentries[kMcParticle] = nkine_filled;

  //---------------------------------------------------------------------------
  // Track data
  Int_t ntrk_filled = 0;    // total number of tracks filled per event
  Int_t ntracklet_filled = 0; // total number of tracklets filled per event
  if (fillCollision)
  {
    Int_t ntrk = fVEvent->GetNumberOfTracks();
    Int_t ntofcls_filled = 0; // total number of TOF clusters filled per event

    for (Int_t itrk = 0; itrk < ntrk; itrk++)
    {
      AliESDtrack *track = nullptr;
      Bool_t deleteTrack = kFALSE;
      if (fESD) {
	track = fESD->GetTrack(itrk);
	deleteTrack = kFALSE; // No need to delete the returned track (?)
      }
      else {
	AliVTrack * aodVTrack = fAOD->GetTrack(itrk);
	AliAODTrack * aodTrack = dynamic_cast<AliAODTrack *>(aodVTrack);
	if (!aodTrack) continue; // Should not happen
	// Skip MUON tracks and constrained tracks
	if (aodTrack->IsMuonTrack()) continue;
	if (aodTrack->IsTPCConstrained()) continue;
	if (aodTrack->IsGlobalConstrained()) continue;
	if (aodTrack->IsHybridGlobalConstrainedGlobal()) continue;
	if (aodTrack->IsHybridTPCConstrainedGlobal()) continue;

	track = new AliESDtrack(aodVTrack);
	deleteTrack = kTRUE; // Since we use new, we have to delete at the end
      }
      //    if (!fTrackFilter.IsSelected(track))
      //      continue;

      tracks.fIndexCollisions = fCollisionCount;
      tracks.fTrackType = TrackTypeEnum::Run2Track;

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
      tracks.fRhoZY = (Char_t)(128. * track->GetSigmaZY() / tracks.fSigmaZ / tracks.fSigmaY);
      tracks.fRhoSnpY = (Char_t)(128. * track->GetSigmaSnpY() / tracks.fSigmaSnp / tracks.fSigmaY);
      tracks.fRhoSnpZ = (Char_t)(128. * track->GetSigmaSnpZ() / tracks.fSigmaSnp / tracks.fSigmaZ);
      tracks.fRhoTglY = (Char_t)(128. * track->GetSigmaTglY() / tracks.fSigmaTgl / tracks.fSigmaY);
      tracks.fRhoTglZ = (Char_t)(128. * track->GetSigmaTglZ() / tracks.fSigmaTgl / tracks.fSigmaZ);
      tracks.fRhoTglSnp = (Char_t)(128. * track->GetSigmaTglSnp() / tracks.fSigmaTgl / tracks.fSigmaSnp);
      tracks.fRho1PtY = (Char_t)(128. * track->GetSigma1PtY() / tracks.fSigma1Pt / tracks.fSigmaY);
      tracks.fRho1PtZ = (Char_t)(128. * track->GetSigma1PtZ() / tracks.fSigma1Pt / tracks.fSigmaZ);
      tracks.fRho1PtSnp = (Char_t)(128. * track->GetSigma1PtSnp() / tracks.fSigma1Pt / tracks.fSigmaSnp);
      tracks.fRho1PtTgl = (Char_t)(128. * track->GetSigma1PtTgl() / tracks.fSigma1Pt / tracks.fSigmaTgl);

      const AliExternalTrackParam *intp = track->GetInnerParam();
      tracks.fTPCinnerP = AliMathBase::TruncateFloatFraction((intp ? intp->GetP() : 0), mTrack1Pt); // Set the momentum to 0 if the track did not reach TPC

      // Compressing and reassigned flags. Keeping only the ones we need.
      tracks.fFlags = 0x0;
      if (track->GetStatus() & AliVTrack::kITSrefit)
        tracks.fFlags |= TrackFlagsRun2Enum::ITSrefit;
      if (track->GetStatus() & AliVTrack::kTPCrefit)
        tracks.fFlags |= TrackFlagsRun2Enum::TPCrefit;

      // add status bit if golden chi2 cut was passed
      if (fESD) {
	const AliESDVertex *vertex = (fESD->GetPrimaryVertex()) ? fESD->GetPrimaryVertex() : fESD->GetPrimaryVertexSPD();
	bool goldenChi2Status = (vertex) ? (track->GetChi2TPCConstrainedVsGlobal(vertex) > 0. && track->GetChi2TPCConstrainedVsGlobal(vertex) < 36.) : false;
	if (goldenChi2Status)
	  tracks.fFlags |= TrackFlagsRun2Enum::GoldenChi2;
      }
      else {
	// FIXME: In case of AOD we suppose the cut has been applied
	  tracks.fFlags |= TrackFlagsRun2Enum::GoldenChi2;
      }

      // Uppermost 4 bits contain PID hypothesis used during tracking
      if (track->GetPIDForTracking() >= 0 && track->GetPIDForTracking() <= 15)
        tracks.fFlags |= track->GetPIDForTracking() << 28;

      tracks.fITSClusterMap = track->GetITSClusterMap();
      tracks.fTPCNClsFindable = track->GetTPCNclsF();

      if ((int)tracks.fTPCNClsFindable - track->GetTPCNcls() >= -128)
        tracks.fTPCNClsFindableMinusFound = tracks.fTPCNClsFindable - track->GetTPCNcls();
      else
        tracks.fTPCNClsFindableMinusFound = -128;

      if ((int)tracks.fTPCNClsFindable - track->GetTPCCrossedRows() >= -128)
        tracks.fTPCNClsFindableMinusCrossedRows = tracks.fTPCNClsFindable - track->GetTPCCrossedRows();
      else
        tracks.fTPCNClsFindableMinusCrossedRows = -128;

      tracks.fTPCNClsShared = (track->GetTPCSharedMap()).CountBits();

      tracks.fTRDPattern = 0;
      for (int i = 0; i < 6; i++)
        if (track->GetTRDslice(i) > 0)
          tracks.fTRDPattern |= 0x1 << i; // flag tracklet on this layer

      // Checking that the track has a TOF measurement matched
      const bool hasTOF = (track->GetStatus() & AliESDtrack::kTOFout) && (track->GetStatus() & AliESDtrack::kTIME);

      tracks.fITSChi2NCl = AliMathBase::TruncateFloatFraction((track->GetITSNcls() ? track->GetITSchi2() / track->GetITSNcls() : 0), mTrackCovOffDiag);
      tracks.fTPCChi2NCl = AliMathBase::TruncateFloatFraction((track->GetTPCNcls() ? track->GetTPCchi2() / track->GetTPCNcls() : 0), mTrackCovOffDiag);
      tracks.fTRDChi2 = AliMathBase::TruncateFloatFraction(track->GetTRDchi2(), mTrackCovOffDiag);
      tracks.fTOFChi2 = AliMathBase::TruncateFloatFraction(hasTOF ? sqrt(track->GetTOFsignalDx() * track->GetTOFsignalDx() + track->GetTOFsignalDz() * track->GetTOFsignalDz()) : -1.f, mTrackCovOffDiag);

      tracks.fTPCSignal = AliMathBase::TruncateFloatFraction(track->GetTPCsignal(), mTrackSignal);
      tracks.fTRDSignal = AliMathBase::TruncateFloatFraction(track->GetTRDsignal(), mTrackSignal);
      // tracks.fTOFSignal = AliMathBase::TruncateFloatFraction(hasTOF ? track->GetTOFsignal() : 0.f, mTrackSignal);
      tracks.fLength = AliMathBase::TruncateFloatFraction(track->GetIntegratedLength(), mTrackSignal);

      // Speed of ligth in TOF units
      constexpr Float_t cspeed = 0.029979246f;
      // PID hypothesis for the momentum extraction
      const AliPID::EParticleType tof_pid = AliPID::kPion;
      const float exp_time = TOFResponse.GetExpectedSignal(track, tof_pid);
      if (exp_time > 0) { // Check that the expected time is positive
        // Expected beta for such hypothesis
        const Float_t exp_beta = (track->GetIntegratedLength() / exp_time / cspeed);

        tracks.fTOFExpMom = exp_beta < 1.f ? AliMathBase::TruncateFloatFraction(
                                               AliPID::ParticleMass(tof_pid) * exp_beta * cspeed /
                                                 TMath::Sqrt(1. - (exp_beta * exp_beta)),
                                               mTrack1Pt)
                                           : 0.f;
      } else {
        tracks.fTOFExpMom = 0.f;
      }

      tracks.fTrackTime = 0.f;
      tracks.fTrackTimeRes = -1.f;
      if (hasTOF) {
        if (track->GetPIDForTracking() >= 0 && track->GetPIDForTracking() <= 15) {
          tracks.fTrackTimeRes = 200e-3;
          const float tofExpMom = tracks.fTOFExpMom / cspeed;
          const float length = tracks.fLength;
          const float massZ = AliPID::ParticleMassZ(track->GetPIDForTracking());
          const float energy = sqrt((massZ * massZ) + (tofExpMom * tofExpMom));
          const float exp = length * energy / (cspeed * tofExpMom);
          tracks.fTrackTime = (track->GetTOFsignal() - exp) * 1e-3; // tof time in ns, taken from the definition of Ruben scaled by 1000 to convert from mus to ns
        }
      }

      tracks.fTrackTime = AliMathBase::TruncateFloatFraction(tracks.fTrackTime, mTrackSignal);
      tracks.fTrackTimeRes = AliMathBase::TruncateFloatFraction(tracks.fTrackTimeRes, mTrackSignal);

      tracks.fTrackEtaEMCAL = AliMathBase::TruncateFloatFraction(track->GetTrackEtaOnEMCal(), mTrackPosEMCAL);
      tracks.fTrackPhiEMCAL = AliMathBase::TruncateFloatFraction(track->GetTrackPhiOnEMCal(), mTrackPosEMCAL);

      if (fTaskMode == kMC)
      {
        // Separate tables (trees) for the MC labels
        Int_t alabel = track->GetLabel();
        // Find the modified label
        Int_t klabel = kineIndex[TMath::Abs(alabel)];
        mctracklabel.fIndexMcParticles = TMath::Abs(klabel) + fOffsetLabel;
        mctracklabel.fMcMask = 0;
        // Use the ITS shared clusters to set the corresponding bits 0-6
        UChar_t itsMask = track->GetITSSharedMap() & 0x1F; // Normally only bits 0-5 are set in Run1/2
        mctracklabel.fMcMask |= itsMask;
        // Use the number of TPC shared clusters as number of TPC mismatches
        // encode in bits 7-9 the values in the ranges 0, 1, 2-3, 4-7, 8-15, 16-31, 32-63, >64
        const TBits *tpcShared = track->GetTPCSharedMapPtr();
        UInt_t tpcCount = tpcShared->CountBits();
        UShort_t tpcMask = 0;
        while (tpcCount > 0)
        {
          tpcCount = tpcCount >> 1;
          tpcMask++;
        }
        if (tpcMask > 7)
          tpcMask = 7;
        mctracklabel.fMcMask |= (tpcMask << 7);
        // TRD (bit 10)
        // We can also use labels per tracklet in the future
        Int_t trdLabel = track->GetTRDLabel();
        if (TMath::Abs(alabel) != TMath::Abs(trdLabel))
          mctracklabel.fMcMask |= (0x1 << 10);
        // TOF (bit 11)
        Int_t tofLabel[3] = {-1};
        track->GetTOFLabel(tofLabel);
        // Check if at least one of the TOF hits matches the track label
        if (!(TMath::Abs(alabel) == TMath::Abs(tofLabel[0]) || TMath::Abs(alabel) == TMath::Abs(tofLabel[1]) || TMath::Abs(alabel) == TMath::Abs(tofLabel[2])))
          mctracklabel.fMcMask |= (0x1 << 11);

        if (alabel < 0 || klabel < 0)
          mctracklabel.fMcMask |= (0x1 << 15);

        FillTree(kMcTrackLabel);
      }

  #ifdef USE_TOF_CLUST
      tofClusters.fTOFncls = track->GetNTOFclusters();
      // Only in case of ESD
      if (fESD && fTreeStatus[kTOF] && tofClusters.fTOFncls > 0)
      {
        Int_t *TOFclsIndex = track->GetTOFclusterArray(); //Index of the matchable cluster (there are fNTOFClusters of them)
        for (Int_t icls = 0; icls < tofClusters.fTOFncls; icls++)
        {
          AliESDTOFCluster *TOFcls = (AliESDTOFCluster *)fESD->GetESDTOFClusters()->At(TOFclsIndex[icls]);
          tofClusters.fToT = TOFcls->GetTOFsignalToT(0);
          tofClusters.fTOFChannel = TOFcls->GetTOFchannel();
          for (Int_t mtchbl = 0; mtchbl < TOFcls->GetNMatchableTracks(); mtchbl++)
          {
            if (TOFcls->GetTrackIndex(mtchbl) != track->GetID())
              continue;
            tofClusters.fDx = TOFcls->GetDx(mtchbl);
            tofClusters.fDz = TOFcls->GetDz(mtchbl);
            tofClusters.fLengthRatio = tracks.fLength > 0 ? TOFcls->GetLength(mtchbl) / tracks.fLength : -1;
            break;
          }
          FillTree(kTOF);
          if (fTreeStatus[kTOF])
            ntofcls_filled++;
        }
      }
  #endif

      // Fill information of the HMPID
      if (track->GetHMPIDsignal() > 0.) {
        if (track->GetHMPIDcluIdx() >= 0) {
          hmpids.fIndexTracks = ntrk_filled + fOffsetTrack;

          Float_t xPc = 0., yPc = 0.;
          Float_t xMip = 0., yMip = 0.;
          Float_t thetaTrk = 0., phiTrk = 0.;
          Int_t nPhot = 0, qMip = 0;

          track->GetHMPIDtrk(xPc, yPc, thetaTrk, phiTrk);
          track->GetHMPIDmip(xMip, yMip, qMip, nPhot);

          hmpids.fHMPIDSignal = AliMathBase::TruncateFloatFraction(track->GetHMPIDsignal(), mTrackSignal);
          hmpids.fHMPIDDistance = AliMathBase::TruncateFloatFraction(TMath::Sqrt((xPc - xMip) * (xPc - xMip) + (yPc - yMip) * (yPc - yMip)), mTrackSignal);
          hmpids.fHMPIDNPhotons = static_cast<Short_t>(nPhot);
          hmpids.fHMPIDQMip = static_cast<Short_t>(qMip);
          FillTree(kHMPID);
          if (fTreeStatus[kHMPID])
            eventextra.fNentries[kHMPID]++;
        }
      }

      // In case we need connection to clusters, activate next lines
      // tracks.fTOFclsIndex += tracks.fNTOFcls;
      // tracks.fNTOFcls = ntofcls_filled;
      FillTree(kTracks);
      FillTree(kTracksCov);
      FillTree(kTracksExtra);
      if (fTreeStatus[kTracks])
        ntrk_filled++;

      if (deleteTrack) delete track;
    } // end loop on tracks

    eventextra.fNentries[kTOF] = ntofcls_filled;

    AliMultiplicity *mlt = dynamic_cast<AliMultiplicity*>(fVEvent->GetMultiplicity());
    Int_t Ntracklets = mlt ? mlt->GetNumberOfTracklets() : 0;

    Float_t theta, phi, dphi, dphiS, dist, x, tgl, alpha;

    for (Int_t itr = Ntracklets; itr--;)
    {
      dphi = mlt->GetDeltaPhi(itr);
      dist = mlt->CalcDist(itr);

      // on-the-fly filtering based on parameters tuned in Run2
      dphiS = TMath::Abs(dphi) - 0.0045;
      if (dphi < 0)
        dphiS = -dphiS;
      if (dist < 1. && dphiS < 0.06)
      {
        theta = mlt->GetTheta(itr);
        phi = mlt->GetPhi(itr);
        tracks.fIndexCollisions = fCollisionCount;
        tracks.fTrackType = TrackTypeEnum::Run2Tracklet;

        // inversion formulas for snp and alpha
        tracks.fSnp = 0.;
        alpha = phi;
        tracks.fAlpha = AliMathBase::TruncateFloatFraction(alpha, mTracklets);

        // inversion formulas for tgl
        x = (TMath::Tan(theta / 2.) - 1.) / (TMath::Tan(theta / 2.) + 1.);
        if (TMath::Log(TMath::Tan(theta / 2)) >= 0)
          tgl = TMath::Sqrt((TMath::Power((1. + TMath::Power(x, 2)) / (1. - TMath::Power(x, 2)), 2)) - 1.);
        else
          tgl = -TMath::Sqrt((TMath::Power((1. + TMath::Power(x, 2)) / (1. - TMath::Power(x, 2)), 2)) - 1.);
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
        // tracks.fTOFSignal = NAN;
        tracks.fLength = NAN;
        tracks.fTOFExpMom = NAN;
        tracks.fTrackTime = NAN;
        tracks.fTrackTimeRes = NAN;

        if (fTaskMode == kMC)
        {
          // Separate tables (trees) for the MC labels: tracklets
          Int_t alabel = mlt->GetLabel(itr, 0); // Take the label of the first layer
          // Find the modified label
          Int_t klabel = kineIndex[TMath::Abs(alabel)];
          mctracklabel.fIndexMcParticles = TMath::Abs(klabel) + fOffsetLabel;
          mctracklabel.fMcMask = 0;
          // Mask fake tracklets
          if (alabel < 0 || klabel < 0)
            mctracklabel.fMcMask |= (0x1 << 15);
          if (mlt->GetLabel(itr, 0) != mlt->GetLabel(itr, 1))
            mctracklabel.fMcMask |= (0x1 << 15);

          FillTree(kMcTrackLabel);
        }

        FillTree(kTracks);
        FillTree(kTracksCov);
        FillTree(kTracksExtra);
        if (fTreeStatus[kTracks]) ntracklet_filled++;
      }
    } // end loop on tracklets
    eventextra.fNentries[kTracks] = ntrk_filled + ntracklet_filled;
    eventextra.fNentries[kTracksCov] = eventextra.fNentries[kTracks];
    eventextra.fNentries[kTracksExtra] = eventextra.fNentries[kTracks];
  }

  //---------------------------------------------------------------------------
  // Calorimeter data

  const double kSecToNanoSec = 1e9;
  AliVCaloCells *cells = fVEvent->GetEMCALCells();
  Short_t nCells = cells->GetNumberOfCells();
  Int_t ncalocells_filled = 0; // total number of calo cells filled per event
  for (Short_t ice = 0; ice < nCells; ++ice)
  {
    Short_t cellNumber;
    Double_t amplitude;
    Double_t time;
    Int_t mclabel;
    Double_t efrac;

    calo.fIndexBCs = fBCCount;

    cells->GetCell(ice, cellNumber, amplitude, time, mclabel, efrac);
    calo.fCellNumber = cellNumber;
    // Mimic run3 compression: Store only cells with energy larger than the threshold
    if (amplitude < fEMCALAmplitudeThreshold)
      continue;
    calo.fAmplitude = AliMathBase::TruncateFloatFraction(amplitude, mCaloAmp);
    calo.fTime = AliMathBase::TruncateFloatFraction(time * kSecToNanoSec, mCaloAmp);
    calo.fCaloType = cells->GetType(); // common for all cells
    calo.fCellType = cells->GetHighGain(ice) ? 1. : 0.;
    FillTree(kCalo);
    if (fTreeStatus[kCalo])
      ncalocells_filled++;
    if (fTaskMode == kMC)
    {
      // Find the modified label
      Int_t klabel = kineIndex[TMath::Abs(mclabel)];
      mccalolabel.fIndexMcParticles = TMath::Abs(klabel) + fOffsetLabel;
      mccalolabel.fMcMask = 0;
      if (mclabel < 0 || klabel < 0)
        mccalolabel.fMcMask |= (0x1 << 15);

      FillTree(kMcCaloLabel);
    }
  } // end loop on calo cells
  eventextra.fNentries[kCalo] = ncalocells_filled;

  bool fillEMCtriggers = true;
  bool fillEMCheaderMedian = true;
  Int_t ncalotriggers_filled = 0;
  TString triggerclasses = fInputEvent->GetFiredTriggerClasses();
  if(!(triggerclasses.Contains("CENT") || triggerclasses.Contains("ALL") || triggerclasses.Contains("-FAST-") || triggerclasses.Contains("CALO"))) {
    // Don't write EMCAL trigger table for clusters where EMCAL was neither trigger nor readout detector
    fillEMCtriggers = false;
  }
  if(fGRP) {
    // In case global run parameters are available write EMCAL trigger 
    // entries only in case EMCAL was a readout nor a trigger detector
    // in order to reduce the data volume used for header words
    TString detectorstring = AliDAQ::ListOfTriggeredDetectors(fGRP->GetDetectorMask());
    if(!detectorstring.Contains("EMCAL")) {
      fillEMCtriggers = false;
    } else {
      // Don't fill median header words for small system (STU not in median mode)
      TString beamtype = fGRP->GetBeamType();
      if((beamtype == "p-p" || beamtype == "p-A" || beamtype == "A-p")) { 
        fillEMCheaderMedian = false;
      }
    }
  }
  if(fillEMCtriggers) {
    // Trigger data for EMCAL:
    // - For full payload (monitoring) events store all non-0 L1 ADCs
    // - For regular events store non-0 L1 ADCs of the 3 leading Gamma patches
    auto geo = AliEMCALGeometry::GetInstance();
    if(!geo) {
      // singleton not yet initialized - initialize it for the first time
      // based on the current run number expecting data from different years
      // will not be mixed in the same job
      geo = AliEMCALGeometry::GetInstanceFromRunNumber(fVEvent->GetRunNumber()); // Needed for EMCAL trigger mapping
    }
    AliVCaloTrigger *calotriggers = fVEvent->GetCaloTrigger("EMCAL");
    Bool_t fullPayload = fEMCALReducedTriggerPayload ? gRandom->Uniform() < fFractionL1MonitorEventsEMCAL : true;
    calotrigger.fIndexBCs = fBCCount;
    calotrigger.fFastOrAbsID = 10001;
    calotrigger.fLnAmplitude = fullPayload ? 1 : 0;
    calotrigger.fTriggerBits = 0;
    calotrigger.fCaloType = 1;
    FillTree(kCaloTrigger);
    int nheader = 1;
    if(fillEMCheaderMedian) {
      // Median for EMCAL
      calotrigger.fFastOrAbsID = 10002;
      calotrigger.fLnAmplitude = calotriggers->GetMedian(0);
      FillTree(kCaloTrigger);
      // Median for DCAL
      calotrigger.fFastOrAbsID = 10003;
      calotrigger.fLnAmplitude = calotriggers->GetMedian(1);
      FillTree(kCaloTrigger);
      nheader += 2;
    }

    calotriggers->Reset();
    ncalotriggers_filled = nheader; // total number of EMCAL triggers filled per event, offset 1 or 3 for the global event properties for EMCAL
    int col, row, fastorID, l1timesums;
    if(fullPayload) {
      // Full payload - store all ADCs
      while (calotriggers->Next())
      {
        calotriggers->GetPosition(col, row);
        calotriggers->GetL1TimeSum(l1timesums);
        if(l1timesums <=0) 
          continue;
        geo->GetTriggerMapping()->GetAbsFastORIndexFromPositionInEMCAL(col, row, fastorID);
        calotrigger.fFastOrAbsID = fastorID;
        calotrigger.fLnAmplitude = l1timesums;
        FillTree(kCaloTrigger);
        if (fTreeStatus[kCaloTrigger])
          ncalotriggers_filled++;
      }
    } else {
      TClonesArray *emcalpatches = dynamic_cast<TClonesArray *>(fInputEvent->FindListObject("EmcalTriggers"));
      if(emcalpatches) {
        AliEMCALTriggerDataGrid<int> l1adcs;
        l1adcs.Allocate(48, 104);
        // Pre-sort the ADCs in the data grid in order to find the ADCs belonging to the 3 leading patches later
        while (calotriggers->Next())
        {
          calotriggers->GetPosition(col, row);
          calotriggers->GetL1TimeSum(l1timesums);
          if(l1timesums <=0) 
            continue;
          l1adcs(col, row) = l1timesums;
        }
        std::vector<AliEMCALTriggerPatchInfo *> allpatches;
        for(int ipatch = 0; ipatch < emcalpatches->GetEntries(); ipatch++) {
          AliEMCALTriggerPatchInfo *nextpatch = static_cast<AliEMCALTriggerPatchInfo *>(emcalpatches->At(ipatch));
          if(nextpatch->IsGammaLowRecalc()) allpatches.emplace_back(nextpatch);
        }
        // sort patches in descending order according to the patch ADC
        std::sort(allpatches.begin(), allpatches.end(), [](const AliEMCALTriggerPatchInfo *lhs, const AliEMCALTriggerPatchInfo *rhs) { return lhs->GetADCAmp() > rhs->GetADCAmp(); } );
        std::vector<int> fastOrIDsInTree;  // in order to avoid double counting
        int npatches = 0;
        for(auto patch : allpatches) {
          for(int icol = patch->GetColStart(); icol < patch->GetColStart() + patch->GetPatchSize(); icol++) {
            for(int irow = patch->GetRowStart(); irow < patch->GetRowStart() + patch->GetPatchSize(); irow++) {
              int adc = l1adcs(icol, irow);
              if(adc > 0) {
                geo->GetTriggerMapping()->GetAbsFastORIndexFromPositionInEMCAL(icol, irow, fastorID);
                if(std::find(fastOrIDsInTree.begin(), fastOrIDsInTree.end(), fastorID) == fastOrIDsInTree.end()) {
                  // FastOr not present in other (selected) trigger patch - add to tree
                  calotrigger.fFastOrAbsID = fastorID;
                  calotrigger.fLnAmplitude = adc;
                  FillTree(kCaloTrigger);
                  if (fTreeStatus[kCaloTrigger])
                    ncalotriggers_filled++;
                  fastOrIDsInTree.emplace_back(fastorID);
                }
              }
            }
          }
          npatches++;
          if(npatches == 3) {
            // select at maximum the 3 leading trigger patches
            break;
          }
        }
      } else {
        AliErrorStream() << "Needs EMCAL patches (from EMCAL trigger maker) for the selection of the FastORs belonging to the 3 leading gamma patches" << std::endl;
      }
    }
  }
  eventextra.fNentries[kCaloTrigger] = ncalotriggers_filled;

  //------PHOS trigger -----------
  const double mPHOSCalib = 0.005; // Mean PHOS calibration 5 MeV/ADC
  // add PHOS trigger digits to list of phos cells
  AliVCaloTrigger *phostriggers = fVEvent->GetCaloTrigger("PHOS");
  phostriggers->Reset();
  calotrigger.fIndexBCs = fBCCount;
  int relid[3];
  Float_t amplitude;
  while (phostriggers->Next())
  {
    //Write trigger digits to same stream as readout cells
    calo.fIndexBCs = fBCCount;
    int triggerbits;
    phostriggers->GetTriggerBits(triggerbits);
    int mod, absId;
    phostriggers->GetPosition(mod, absId);
    //here absId is normal Run2 readout absId
    //Remove noisy triggers
    Int_t phosmodulenumber = TMath:: Ceil(float(absId)/3584) ;
    if (fUsePHOSBadMap)
    {
      int id = absId - ( phosmodulenumber - 1 ) * 3584 ;
      int ix = (Int_t)TMath::Ceil( float(id) / 64 )  ;
      int iz = (Int_t)( id - ( ix - 1 ) * 64 ) ;
      if(fPHOSBadMap[phosmodulenumber] != nullptr && fPHOSBadMap[phosmodulenumber]->GetBinContent(ix,iz)>0) { //bad channel
        continue ;
      }
    }
    //transform to Run3 truID
    absId--;
    relid[0] = 4 - absId / 3584  ;  //Aliroot<->O2 module numbering
    absId = absId % 3584  ;  //module
    relid[1] = absId / 64  ; //x
    relid[2] = absId % 64  ; //z

    relid[0] = relid[0]*4 -2 + relid[1]/16 ;
    relid[1] = (relid[1]%16)/2 ;
    relid[2] /=2 ;

    Int_t truId= relid[0] * 224 + // the offset of PHOS modules
                 relid[1] +       // the offset along phi
                 relid[2] * 8;    // the offset along z
    // filter null entries: they usually have negative entries and no trigger bits
    // in case of trigger bits the energy can be 0 or negative but the trigger position is marked
    // store trigger
    calotrigger.fCaloType = 0 ; //PHOS
    calotrigger.fFastOrAbsID = truId;

    phostriggers->GetAmplitude(amplitude);
    calotrigger.fLnAmplitude = AliMathBase::TruncateFloatFraction(amplitude/mPHOSCalib, 0xFFF); //12 bit
    if(triggerbits==0){ //L0 trigger
      calotrigger.fTriggerBits =0 ; //0:L0, 1:L1
    }
    else{
      int timesum;
      phostriggers->GetL1TimeSum(timesum);
      calotrigger.fTriggerBits =1+timesum ; // 1,2,3:L1
    }
    FillTree(kCaloTrigger);
    if (fTreeStatus[kCaloTrigger])
      ncalotriggers_filled++;
  }
  eventextra.fNentries[kCaloTrigger] = ncalotriggers_filled;




  AliVCaloCells * phoscells = fVEvent->GetPHOSCells();
  Int_t nphoscells_filled = 0;
  Int_t nPHCells = phoscells->GetNumberOfCells();
  for (Short_t icp = 0; icp < nPHCells; ++icp)
  {
    Short_t cellNumber;
    Double_t amplitude,time;
    Int_t mclabel;
    Double_t efrac;

    calo.fIndexBCs = fBCCount;
    phoscells->GetCell(icp, cellNumber, amplitude, time, mclabel, efrac);
    //Run2: absId=1..4*56*64 ; Run3: absId = 32*56...4*56*64, module numbering is opposite
    int mod = cellNumber/3584;
    calo.fCellNumber = (4-mod)*3584 + cellNumber%3584 ;
    //Run3: uncalibrated amplitude in ADC counts
    // here we assume fixed calibration
    calo.fAmplitude = AliMathBase::TruncateFloatFraction(amplitude/mPHOSCalib, 0xFFF); //12 bit
    calo.fTime = AliMathBase::TruncateFloatFraction(time, 0x1FFF);  //13 bit
    calo.fCellType = phoscells->GetHighGain(icp) ? 0. : 1.;
    calo.fCaloType = phoscells->GetType();

    FillTree(kCalo);
    if (fTreeStatus[kCalo])
      nphoscells_filled++;
    if (fTaskMode == kMC)
    {
      // Find the modified label
      if(mclabel>=0){  //label -1 == no primary
        Int_t klabel = kineIndex[mclabel];
        mccalolabel.fIndexMcParticles = TMath::Abs(klabel) + fOffsetLabel;
        mccalolabel.fMcMask = 0;
        if ( klabel < 0)
          mccalolabel.fMcMask |= (0x1 << 15);

        FillTree(kMcCaloLabel);
      }
    }
  } // end loop on PHOS cells
  eventextra.fNentries[kCalo] = nphoscells_filled;

  //---------------------------------------------------------------------------
  // Muon tracks -> fwdtracks

  Int_t nmu = fVEvent->GetNumberOfMuonTracks();

  TRefArray muonTracks; // Needed for AOD
  if (fAOD) fAOD->GetMuonTracks(&muonTracks);

  Int_t nmu_filled = 0;   // total number of muons filled per event

  if (fESD) {
    for (Int_t imu = 0; imu < nmu; ++imu)
    {
      AliESDMuonTrack *mutrk = fESD->GetMuonTrack(imu);
      if (mutrk->GetInverseBendingMomentum()==FLT_MAX && mutrk->GetThetaX()==0 && mutrk->GetThetaY()==0 && mutrk->GetZ()==0) continue; //skip tracks whose parameter values are still at their default (not real muon tracks)
      fwdtracks = MUONtoFwdTrack(*mutrk);
      fwdtracks.fIndexCollisions = fCollisionCount;

      // Now MUON clusters for the current track
      // Int_t muTrackID = fOffsetMuTrackID + imu;
      Int_t nmucl = mutrk->GetNClusters();

      fwdtracks.fNClusters = nmucl;

      FillTree(kFwdTrack);
      FillTree(kFwdTrackCov);
      if (fTreeStatus[kFwdTrack])
      nmu_filled++;
    } // End loop on muon tracks
  }
  else {
    for (Int_t imu = 0; imu < nmu; ++imu) {
      //PH It seems the MUON information in the "standard" AOD is not really useful.
      AliAODTrack *mutrk = dynamic_cast<AliAODTrack*>(muonTracks.At(imu));

      if (!mutrk) continue; //PH This should not be needed!

      fwdtracks = MUONtoFwdTrack(*mutrk);
      fwdtracks.fIndexCollisions = fCollisionCount;

      // No MUON clusters for the AOD track

      FillTree(kFwdTrack);
      FillTree(kFwdTrackCov);
      if (fTreeStatus[kFwdTrack])
      nmu_filled++;
     } // End loop on muon AOD tracks
  }

  eventextra.fNentries[kFwdTrack] = nmu_filled;
  eventextra.fNentries[kFwdTrackCov] = nmu_filled;

  //---------------------------------------------------------------------------
  // ZDC
  AliESDZDC *esdzdc = fESD ? fESD->GetESDZDC() : nullptr;
  AliAODZDC *aodzdc = fAOD ? fAOD->GetZDCData() : nullptr;
  zdc.fIndexBCs = fBCCount;

  zdc.fTimeZNA = 999.f;
  zdc.fTimeZNC = 999.f;
  zdc.fTimeZPA = 999.f;
  zdc.fTimeZPC = 999.f;
  zdc.fTimeZEM1 = 999.f;
  zdc.fTimeZEM2 = 999.f;

  if (esdzdc) {
    // ZEM
    zdc.fEnergyZEM1 = esdzdc->GetZEM1Energy();
    zdc.fEnergyZEM2 = esdzdc->GetZEM2Energy();
    zdc.fEnergyCommonZNA = esdzdc->GetZNATowerEnergy()[0];
    zdc.fEnergyCommonZNC = esdzdc->GetZNCTowerEnergy()[0];
    zdc.fEnergyCommonZPA = esdzdc->GetZPATowerEnergy()[0];
    zdc.fEnergyCommonZPC = esdzdc->GetZPCTowerEnergy()[0];

    // ZDC (P,N) sectors
    for (Int_t ich = 0; ich < 4; ++ich) {
      zdc.fEnergySectorZNA[ich] = esdzdc->GetZNATowerEnergy()[ich + 1];
      zdc.fEnergySectorZNC[ich] = esdzdc->GetZNCTowerEnergy()[ich + 1];
      zdc.fEnergySectorZPA[ich] = esdzdc->GetZPATowerEnergy()[ich + 1];
      zdc.fEnergySectorZPC[ich] = esdzdc->GetZPCTowerEnergy()[ich + 1];
    }
    // ZDC TDC
    Bool_t isHitFlagFilled = fESD->GetRunNumber() >= 208502;
    Bool_t isZNAhit = isHitFlagFilled ? esdzdc->IsZNAhit() : 1;
    Bool_t isZNChit = isHitFlagFilled ? esdzdc->IsZNChit() : 1;
    Bool_t isZPAhit = isHitFlagFilled ? esdzdc->IsZPAhit() : 1;
    Bool_t isZPChit = isHitFlagFilled ? esdzdc->IsZPChit() : 1;
    Bool_t isZEM1hit = isHitFlagFilled ? esdzdc->IsZEM1hit() : 1;
    Bool_t isZEM2hit = isHitFlagFilled ? esdzdc->IsZEM2hit() : 1;

    // Storing first ZDC hit in +/-12.5 ns around 0
    for (Int_t i = 0; i < 4; i++) {
      Float_t tZNA = isZNAhit ? esdzdc->GetZDCTDCCorrected(esdzdc->GetZNATDCChannel(), i) : 999.f;
      Float_t tZNC = isZNChit ? esdzdc->GetZDCTDCCorrected(esdzdc->GetZNCTDCChannel(), i) : 999.f;
      Float_t tZPA = isZPAhit ? esdzdc->GetZDCTDCCorrected(esdzdc->GetZPATDCChannel(), i) : 999.f;
      Float_t tZPC = isZPChit ? esdzdc->GetZDCTDCCorrected(esdzdc->GetZPCTDCChannel(), i) : 999.f;
      Float_t tZEM1 = isZEM1hit ? esdzdc->GetZDCTDCCorrected(esdzdc->GetZEM1TDCChannel(), i) : 999.f;
      Float_t tZEM2 = isZEM2hit ? esdzdc->GetZDCTDCCorrected(esdzdc->GetZEM2TDCChannel(), i) : 999.f;
      if (tZNA > -12.5 && tZNA < 12.5 && zdc.fTimeZNA > 998)
	zdc.fTimeZNA = tZNA;
      if (tZNC > -12.5 && tZNC < 12.5 && zdc.fTimeZNC > 998)
	zdc.fTimeZNC = tZNC;
      if (tZPA > -12.5 && tZPA < 12.5 && zdc.fTimeZPA > 998)
	zdc.fTimeZPA = tZPA;
      if (tZPC > -12.5 && tZPC < 12.5 && zdc.fTimeZPC > 998)
	zdc.fTimeZPC = tZPC;
      if (tZEM1 > -12.5 && tZEM1 < 12.5 && zdc.fTimeZEM1 > 998)
	zdc.fTimeZEM1 = tZEM1;
      if (tZEM2 > -12.5 && tZEM2 < 12.5 && zdc.fTimeZEM2 > 998)
	zdc.fTimeZEM2 = tZEM2;
    }
  }
  else {
    // ZEM
    zdc.fEnergyZEM1 = aodzdc->GetZEM1Energy();
    zdc.fEnergyZEM2 = aodzdc->GetZEM2Energy();
    zdc.fEnergyCommonZNA = aodzdc->GetZNATowerEnergy()[0];
    zdc.fEnergyCommonZNC = aodzdc->GetZNCTowerEnergy()[0];
    zdc.fEnergyCommonZPA = aodzdc->GetZPATowerEnergy()[0];
    zdc.fEnergyCommonZPC = aodzdc->GetZPCTowerEnergy()[0];

    // ZDC (P,N) sectors
    for (Int_t ich = 0; ich < 4; ++ich) {
      zdc.fEnergySectorZNA[ich] = aodzdc->GetZNATowerEnergy()[ich + 1];
      zdc.fEnergySectorZNC[ich] = aodzdc->GetZNCTowerEnergy()[ich + 1];
      zdc.fEnergySectorZPA[ich] = aodzdc->GetZPATowerEnergy()[ich + 1];
      zdc.fEnergySectorZPC[ich] = aodzdc->GetZPCTowerEnergy()[ich + 1];
    }
    // ZDC TDC
    // Storing first ZDC hit in +/-12.5 ns around 0
    zdc.fTimeZNA = aodzdc->GetZNATime();
    zdc.fTimeZNC = aodzdc->GetZNCTime();
    zdc.fTimeZPA = aodzdc->GetZPATime();
    zdc.fTimeZPC = aodzdc->GetZPCTime();
    //PH ZEM1,2 times are not in AOD
  }

  FillTree(kZdc);
  if (fTreeStatus[kZdc])
    eventextra.fNentries[kZdc] = 1;

  //---------------------------------------------------------------------------
  // V0A and V0C
  AliVVZERO *vz = fVEvent->GetVZEROData();
  fv0a.fIndexBCs = fBCCount;
  fv0c.fIndexBCs = fBCCount;
  fv0a.fAmplitude_size = 0;
  fv0c.fAmplitude_size = 0;
  fv0a.fChannel_size = 0;
  fv0c.fChannel_size = 0;
  for (Int_t ich = 0; ich < 32; ++ich) {
    if (vz->GetMultiplicityV0A(ich) > 0) {
      fv0a.fAmplitude[fv0a.fAmplitude_size++] = AliMathBase::TruncateFloatFraction(vz->GetMultiplicityV0A(ich), mV0Amplitude);
      fv0a.fChannel[fv0a.fChannel_size++] = ich;
    }
  }
  for (Int_t ich = 0; ich < 32; ++ich) {
    if (vz->GetMultiplicityV0C(ich) > 0) {
      fv0c.fAmplitude[fv0c.fAmplitude_size++] = AliMathBase::TruncateFloatFraction(vz->GetMultiplicityV0C(ich), mV0Amplitude);
      fv0c.fChannel[fv0c.fChannel_size++] = ich;
    }
  }
  fv0a.fTime = AliMathBase::TruncateFloatFraction(vz->GetV0ATime(), mV0Time);
  fv0c.fTime = AliMathBase::TruncateFloatFraction(vz->GetV0CTime(), mV0Time);
  fv0a.fTriggerMask = 0; // not filled for the moment
  FillTree(kFV0A);
  FillTree(kFV0C);
  if (fTreeStatus[kFV0A])
    eventextra.fNentries[kFV0A] = 1;
  if (fTreeStatus[kFV0C])
    eventextra.fNentries[kFV0C] = 1;

  //---------------------------------------------------------------------------
  // FT0
  ft0.fIndexBCs = fBCCount;
  ft0.fChannelA_size = 0;
  ft0.fChannelC_size = 0;
  if (fESD) {
    for (Int_t ich = 0; ich < 12; ++ich) {
      if (fESD->GetT0amplitude()[ich + 12] > 0) {
        ft0.fAmplitudeA[ft0.fChannelA_size] = AliMathBase::TruncateFloatFraction(fESD->GetT0amplitude()[ich + 12], mT0Amplitude);
        ft0.fChannelA[ft0.fChannelA_size++] = ich;
      }
    }
    for (Int_t ich = 0; ich < 12; ++ich) {
      if (fESD->GetT0amplitude()[ich] > 0) {
        ft0.fAmplitudeC[ft0.fChannelC_size] = AliMathBase::TruncateFloatFraction(fESD->GetT0amplitude()[ich], mT0Amplitude);
        ft0.fChannelC[ft0.fChannelC_size++] = ich;
      }
    }
    ft0.fTimeA = AliMathBase::TruncateFloatFraction(fESD->GetT0TOF(1) * 1e-3, mT0Time); // ps to ns
    ft0.fTimeC = AliMathBase::TruncateFloatFraction(fESD->GetT0TOF(2) * 1e-3, mT0Time); // ps to ns
    ft0.fTriggerMask = fESD->GetT0Trig();
  }
  else {
    AliAODTZERO * aodtzero = fAOD->GetTZEROData();
    for (Int_t ich = 0; ich < 12; ++ich) {
      if (aodtzero->GetAmp(ich + 12) > 0) {
        ft0.fAmplitudeA[ft0.fChannelA_size] = AliMathBase::TruncateFloatFraction(aodtzero->GetAmp(ich + 12), mT0Amplitude);
        ft0.fChannelA[ft0.fChannelA_size++] = ich;
      }
    }
    for (Int_t ich = 0; ich < 12; ++ich) {
      if (aodtzero->GetAmp(ich) > 0) {
        ft0.fAmplitudeC[ft0.fChannelC_size] = AliMathBase::TruncateFloatFraction(aodtzero->GetAmp(ich), mT0Amplitude);
        ft0.fChannelC[ft0.fChannelC_size++] = ich;
      }
    }
    ft0.fTimeA = AliMathBase::TruncateFloatFraction(aodtzero->GetT0TOF(1) * 1e-3, mT0Time); // ps to ns
    ft0.fTimeC = AliMathBase::TruncateFloatFraction(aodtzero->GetT0TOF(2) * 1e-3, mT0Time); // ps to ns
    ft0.fTriggerMask = 0; // Not available in AOD
  }
  // will be removed with O2 improvement (one size field used for two VLAs)
  ft0.fAmplitudeA_size = ft0.fChannelA_size;
  ft0.fAmplitudeC_size = ft0.fChannelC_size;

  FillTree(kFT0);
  if (fTreeStatus[kFT0])
    eventextra.fNentries[kFT0] = 1;

  //---------------------------------------------------------------------------
  // AD (FDD)
  if (fVEvent->GetADData()) {
    fdd.fIndexBCs = fBCCount;
    if (fESD) {
      AliESDAD *esdad = fESD->GetADData();
      for (Int_t ich = 0; ich < 8; ++ich)
        fdd.fChargeA[ich] = int16_t(esdad->GetMultiplicityADA(ich));
      for (Int_t ich = 0; ich < 8; ++ich)
        fdd.fChargeC[ich] = int16_t(esdad->GetMultiplicityADC(ich));
      fdd.fTimeA = AliMathBase::TruncateFloatFraction(esdad->GetADATime(), mADTime);
      fdd.fTimeC = AliMathBase::TruncateFloatFraction(esdad->GetADCTime(), mADTime);
      fdd.fTriggerMask = 0; // not filled for the moment
    }
    else {
      AliAODAD *aodad = fAOD->GetADData();
      for (Int_t ich = 0; ich < 8; ++ich)
        fdd.fChargeA[ich] = int16_t(aodad->GetMultiplicityADA(ich));
      for (Int_t ich = 0; ich < 8; ++ich)
        fdd.fChargeC[ich] = int16_t(aodad->GetMultiplicityADC(ich));
      fdd.fTimeA = AliMathBase::TruncateFloatFraction(aodad->GetADATime(), mADTime);
      fdd.fTimeC = AliMathBase::TruncateFloatFraction(aodad->GetADCTime(), mADTime);
      fdd.fTriggerMask = 0; // not filled for the moment
    }
    FillTree(kFDD);
    if (fTreeStatus[kFDD])
      eventextra.fNentries[kFDD] = 1;
  }

  Int_t nv0_filled = 0; // total number of v0's filled per event
  Int_t n2prong_filled = 0;
  if (fillCollision)
  {
    //---------------------------------------------------------------------------
    // V0s (Lambda and K0S)
    Int_t nv0 = fVEvent->GetNumberOfV0s();
    Int_t* v0Lookup = new Int_t[nv0]; // needed for HF indices
    v0s.fIndexCollisions = fCollisionCount;
    for (Int_t iv0 = 0; iv0 < nv0; ++iv0) {
      v0Lookup[iv0] = -1; // not stored
      if (fESD) {
	AliESDv0 *v0 = fESD->GetV0(iv0);
	// select only "offline" V0s, skip the "on-the-fly" ones
	if (!v0 || v0->GetOnFlyStatus())
          continue;
        Int_t pidx = v0->GetPindex();
        Int_t nidx = v0->GetNindex();
        v0s.fIndexTracksPos = TMath::Abs(pidx) + fOffsetTrack; // Positive track ID
        v0s.fIndexTracksNeg = TMath::Abs(nidx) + fOffsetTrack; // Negative track ID
      }
      else {
	AliAODv0 *v0 = fAOD->GetV0(iv0);
	// select only "offline" V0s, skip the "on-the-fly" ones
	if (!v0 || v0->GetOnFlyStatus())
          continue;
        Int_t pidx = v0->GetPosID();
        Int_t nidx = v0->GetNegID();
        v0s.fIndexTracksPos = TMath::Abs(pidx) + fOffsetTrack; // Positive track ID
        v0s.fIndexTracksNeg = TMath::Abs(nidx) + fOffsetTrack; // Negative track ID
      }

      v0Lookup[iv0] = nv0_filled; // stored
      FillTree(kV0s);
      if (fTreeStatus[kV0s])
	nv0_filled++;
    } // End loop on V0s
    eventextra.fNentries[kV0s] = nv0_filled;

    //---------------------------------------------------------------------------
    // Cascades
    // If we do not have V0s, we do not have cascades
    Int_t ncascades_filled = 0; // total number of cascades filled per event
    if (nv0 > 0)
    {
      // Combine the track indexes of V0 daughters in unique identifier
      ULong64_t *packedPosNeg = new ULong64_t[nv0];
      ULong64_t *sortedPosNeg = new ULong64_t[nv0];
      Int_t *sortIdx = new Int_t[nv0];

      //Fill cascades in order
      Int_t ncas = fVEvent->GetNumberOfCascades();
      Int_t ncastosort = 0;
      ULong64_t *packedV0indices       = new ULong64_t[ncas];
      ULong64_t *packedbachelorindices = new ULong64_t[ncas];
      Int_t *sortV0Idx = new Int_t[ncas];

      //Don't forget that OTF V0s might exist
      Int_t nv0offline = 0;
      for (Int_t iv0 = 0; iv0 < nv0; ++iv0)
      {
	if (fESD) {
	  AliESDv0 *v0 = fESD->GetV0(iv0);
	  if (v0 && !v0->GetOnFlyStatus()) {
	    packedPosNeg[nv0offline++] = (((ULong64_t)(v0->GetPindex())) << 31) | ((ULong64_t)(v0->GetNindex()));
	  }
	}
	else {
	  AliAODv0 *v0 = fAOD->GetV0(iv0);
	  if (v0 && !v0->GetOnFlyStatus()) {
	    packedPosNeg[nv0offline++] = (((ULong64_t)(v0->GetPosID())) << 31) | ((ULong64_t)(v0->GetNegID()));
	  }
	}
      }
      TMath::Sort(nv0offline, packedPosNeg, sortIdx, kFALSE);
      for (Int_t iv0 = 0; iv0 < nv0offline; ++iv0)
      {
        sortedPosNeg[iv0] = packedPosNeg[sortIdx[iv0]];
      }

      for (Int_t icas = 0; icas < ncas; ++icas)
      {
	if (fESD) {
	  AliESDcascade *cas = fESD->GetCascade(icas);
	  // Select only cascades containing "offline" V0s
	  if (cas && !cas->GetOnFlyStatus()) {
	    // Find the identifier of the V0 using the indexes of its daughters
	    ULong64_t currV0 = (((ULong64_t)(cas->GetPindex())) << 31) | ((ULong64_t)(cas->GetNindex()));
	    // Use binary search in the sorted array
	    Int_t v0idx = TMath::BinarySearch(nv0offline, sortedPosNeg, currV0);
	    // Check if the match is exact
	    if (sortedPosNeg[v0idx] == currV0) {
        packedV0indices[ncastosort] = sortIdx[v0idx] + fOffsetV0;
        packedbachelorindices[ncastosort] = cas->GetBindex() + fOffsetTrack;
        ncastosort++;
	    }
	  }
	}
	else {
	  AliAODcascade *cas = fAOD->GetCascade(icas);
	  // Select only cascades containing "offline" V0s
	  if (cas && !cas->GetOnFlyStatus()) {
	    // Find the identifier of the V0 using the indexes of its daughters
	    ULong64_t currV0 = (((ULong64_t)(cas->GetPosID())) << 31) | ((ULong64_t)(cas->GetNegID()));
	    // Use binary search in the sorted array
	    Int_t v0idx = TMath::BinarySearch(nv0offline, sortedPosNeg, currV0);
	    // Check if the match is exact
	    if (sortedPosNeg[v0idx] == currV0) {
        packedV0indices[ncastosort] = sortIdx[v0idx] + fOffsetV0;
        packedbachelorindices[ncastosort] = cas->GetBachID() + fOffsetTrack;
        ncastosort++;
	    }
	  }
	}
      } // End loop on cascades

      //Sort cascades
      TMath::Sort(ncastosort, packedV0indices, sortV0Idx, kFALSE);
      //Fill cascades only after V0 sorting
      cascs.fIndexCollisions = fCollisionCount;
      for (Int_t icas = 0; icas < ncastosort; ++icas){
        //Fill tree only with ordered information
        cascs.fIndexV0s    = packedV0indices[sortV0Idx[icas]];
        cascs.fIndexTracks = packedbachelorindices[sortV0Idx[icas]];
        FillTree(kCascades);
        ncascades_filled++;
      }

      delete[] packedPosNeg;
      delete[] sortedPosNeg;
      delete[] sortIdx;

      delete[] packedV0indices;
      delete[] packedbachelorindices;
      delete[] sortV0Idx;

      // Store HF here when we have the V0 infos in place
      if (fStoreHF) {
        // Read input from HF task
        TTree *hf2ProngCandidateTree = (TTree*) fInputHandler->GetUserInfo()->FindObject("hf2ProngCandidateTree");
        //Printf("HF hf2ProngCandidateTree has %lld entries", hf2ProngCandidateTree->GetEntries());
        hf2ProngCandidateTree->SetBranchAddress("trackId0", &hf2Prong.fIndexTracks_0);
        hf2ProngCandidateTree->SetBranchAddress("trackId1", &hf2Prong.fIndexTracks_1);
        hf2ProngCandidateTree->SetBranchAddress("hfflag", &hf2Prong.fHFflag);

        for (int i=0; i<hf2ProngCandidateTree->GetEntries(); i++) {
          hf2ProngCandidateTree->GetEntry(i);
          hf2Prong.fIndexTracks_0 += fOffsetTrack;
          hf2Prong.fIndexTracks_1 += fOffsetTrack;
          FillTree(kHF2Prong);
          n2prong_filled++;
        }
        eventextra.fNentries[kHF2Prong] = n2prong_filled;

        TTree *hf3ProngCandidateTree = (TTree*) fInputHandler->GetUserInfo()->FindObject("hf3ProngCandidateTree");
        //Printf("HF hf3ProngCandidateTree has %lld entries", hf3ProngCandidateTree->GetEntries());
        hf3ProngCandidateTree->SetBranchAddress("trackId0", &hf3Prong.fIndexTracks_0);
        hf3ProngCandidateTree->SetBranchAddress("trackId1", &hf3Prong.fIndexTracks_1);
        hf3ProngCandidateTree->SetBranchAddress("trackId2", &hf3Prong.fIndexTracks_2);
        hf3ProngCandidateTree->SetBranchAddress("hfflag", &hf3Prong.fHFflag);

        for (int i=0; i<hf3ProngCandidateTree->GetEntries(); i++) {
          hf3ProngCandidateTree->GetEntry(i);
          hf3Prong.fIndexTracks_0 += fOffsetTrack;
          hf3Prong.fIndexTracks_1 += fOffsetTrack;
          hf3Prong.fIndexTracks_2 += fOffsetTrack;
          FillTree(kHF3Prong);
        }
        eventextra.fNentries[kHF3Prong] = hf3ProngCandidateTree->GetEntries();

        TTree *hfDstarCandidateTree = (TTree*) fInputHandler->GetUserInfo()->FindObject("hfDstarCandidateTree");
        //Printf("HF hfDstarCandidateTree has %lld entries", hfDstarCandidateTree->GetEntries());
        hfDstarCandidateTree->SetBranchAddress("trackSoftPi", &hfDStar.fIndexTracks_0);
        hfDstarCandidateTree->SetBranchAddress("trackD0", &hfDStar.fIndexHf2Prongs);

        for (int i=0; i<hfDstarCandidateTree->GetEntries(); i++) {
          hfDstarCandidateTree->GetEntry(i);
          hfDStar.fIndexTracks_0 += fOffsetTrack;
          hfDStar.fIndexHf2Prongs += fOffsetHF2Prong;
          FillTree(kHFDStar);
        }
        eventextra.fNentries[kHFDStar] = hfDstarCandidateTree->GetEntries();

        TTree *hfCascadeCandidateTree = (TTree*) fInputHandler->GetUserInfo()->FindObject("hfCascadeCandidateTree");
        //Printf("HF hfCascadeCandidateTree has %lld entries", hfCascadeCandidateTree->GetEntries());
        Int_t v0Index = -1;
        hfCascadeCandidateTree->SetBranchAddress("v0index", &v0Index);
        hfCascadeCandidateTree->SetBranchAddress("trackBachel", &hfCascades.fIndexTracks_0);

        Int_t nhfcascades_filled = 0;
        for (int i=0; i<hfCascadeCandidateTree->GetEntries(); i++) {
          hfCascadeCandidateTree->GetEntry(i);
          hfCascades.fIndexTracks_0 += fOffsetTrack;
          if (v0Lookup[v0Index] == -1) {
            AliWarning(Form("V0 %d is requested for a HF cascade but not planned to be stored. Skipping this HF cascade.", v0Index));
          } else {
            hfCascades.fIndexV0s = v0Lookup[v0Index] + fOffsetV0;
            FillTree(kHFCascade);
            nhfcascades_filled++;
          }
        }
        eventextra.fNentries[kHFCascade] = nhfcascades_filled;
      }

      delete[] v0Lookup;
    } // End if V0s
    eventextra.fNentries[kCascades] = ncascades_filled;
  }

  //---------------------------------------------------------------------------
  // MC data (to be modified)

  if (MCEvt || MCHeader)
  {
    // MC vertex
    const AliVVertex *MCvtx = MCEvt ? MCEvt->GetPrimaryVertex() : nullptr;
    if (!MCvtx && !MCHeader) //Check on the MC vertex
      AliFatal("Could not retrieve MC vertex");

    // Fill MC collision
    mccollision.fIndexBCs = fBCCount;

    mccollision.fPosX = AliMathBase::TruncateFloatFraction(MCvtx ? MCvtx->GetX() : MCHeader->GetVtxX(), mCollisionPosition);
    mccollision.fPosY = AliMathBase::TruncateFloatFraction(MCvtx ? MCvtx->GetY() : MCHeader->GetVtxY(), mCollisionPosition);
    mccollision.fPosZ = AliMathBase::TruncateFloatFraction(MCvtx ? MCvtx->GetZ() : MCHeader->GetVtxZ(), mCollisionPosition);

    AliGenEventHeader *mcGenH = MCEvt ? MCEvt->GenEventHeader()  : MCHeader->GetCocktailHeader(0); //PH Probably not OK
    mccollision.fT = AliMathBase::TruncateFloatFraction(mcGenH ? mcGenH->InteractionTime() : -999., mCollisionPosition);
    mccollision.fWeight = AliMathBase::TruncateFloatFraction(mcGenH ? mcGenH->EventWeight() : 1., mCollisionPosition);

    // Impact parameter
    AliCollisionGeometry *cGeo = dynamic_cast<AliCollisionGeometry *>(mcGenH);
    mccollision.fImpactParameter = (cGeo ? cGeo->ImpactParameter() : -999.f);

    mccollision.fGeneratorsID = 0;
    for (Int_t gen = 0; gen < kGenerators; gen++)
    {
      if (mcGenH && mcGenH->InheritsFrom(Generator[gen]))
        SETBIT(mccollision.fGeneratorsID, gen);
      else
        CLRBIT(mccollision.fGeneratorsID, gen);
    }
    if (mcGenH && mcGenH->InheritsFrom(Generator[kAliGenCocktailEventHeader]))
    {
      //PH This part probably needs modifications for the MC AOD
      TList *headers = ((AliGenCocktailEventHeader *)mcGenH)->GetHeaders();
      TListIter cocktail(headers);
      TObject *to = 0x0;
      while ((to = cocktail()))
      {
        if (mccollision.fImpactParameter < 0)
        {
          // Change the impact parameter if not set
          AliCollisionGeometry *toCGeo = dynamic_cast<AliCollisionGeometry *>(to);
          mccollision.fImpactParameter = (toCGeo ? toCGeo->ImpactParameter() : -999.f);
        }
        for (Int_t gen = 0; gen < kGenerators; gen++)
        {
          if (to->InheritsFrom(Generator[gen]))
            SETBIT(mccollision.fGeneratorsID, gen);
        }
      }
    }
    mccollision.fImpactParameter = AliMathBase::TruncateFloatFraction(mccollision.fImpactParameter, mCollisionPosition);
    eventextra.fNentries[kMcCollision] = 1;

    // HepMC inforrmation
    hepMcCrossSections.Reset();
    hepMcCrossSections.fIndexMcCollisions = fBCCount;
    hepMcCrossSections.Fill(mcGenH, this);

    hepMcPdfInfo.Reset();
    hepMcPdfInfo.fIndexMcCollisions = fBCCount;
    hepMcPdfInfo.Fill(mcGenH, this);

    hepMcHeavyIon.Reset();
    hepMcHeavyIon.fIndexMcCollisions = fBCCount;
    hepMcHeavyIon.Fill(mcGenH, this);
  }
  else
  {
    eventextra.fNentries[kMcCollision] = 0;
  }
  // Filling the tree of vertices has to be done last because it contains the
  // index data for the other trees
  FillTree(kMcCollision);

  // MC collision label
  // will be joined with Collisions, therefore fill it only when we fill Collisions
  if (fillCollision) {
    mccollisionlabel.fIndexMcCollisions = fBCCount;
    mccollisionlabel.fMcMask = 0;
    FillTree(kMcCollisionLabel);
  }

  //---------------------------------------------------------------------------
  // Update the offsets at the end of each collision
  fOffsetLabel += nkine_filled; // Offset for the labels of the next event
  fOffsetTrack += ntrk_filled + ntracklet_filled;
  fOffsetMuTrackID += nmu_filled;
  fOffsetV0 += nv0_filled;
  fOffsetHF2Prong += n2prong_filled;

  // Update event counters
  if (fillCollision)
    fCollisionCount++;
  fBCCount++;
} // void AliAnalysisTaskAO2Dconverter::FillEventInTF()

void AliAnalysisTaskAO2Dconverter::FinishTF()
{
  // Write all trees
  for (Int_t i = 0; i < kTrees; i++)
    WriteTree((TreeIndex)i);
  // Remove trees
  for (Int_t i = 0; i < kTrees; i++)
    if (fTree[i])
    {
      delete fTree[i];
      fTree[i] = 0x0;
    }
} // AliAnalysisTaskAO2Dconverter::FinishTF()

Bool_t AliAnalysisTaskAO2Dconverter::Select(TParticle *part, Float_t rv, Float_t zv)
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

  if (TMath::Abs(eta) < 2.5 && rv < 170)
    return kTRUE;

  // Andreas' Cuts
  //  if (TMath::Abs(eta) < 1. && rv < 170)return kTRUE;

  // Muon arm
  if (eta > -4.2 && eta < -2.3 && zv > -500.)
    return kTRUE; // Muon arms

  // PMD acceptance 2.3 <= eta < = 3.5
  //  if(eta>2.0&&eta<3.7)return kTRUE;

  return kFALSE;

} // AliAnalysisTaskAO2Dconverter::Select(TParticle* part, Float_t rv, Float_t zv)

//_________________________________________________________________________________________________
AliAnalysisTaskAO2Dconverter::FwdTrackPars AliAnalysisTaskAO2Dconverter::MUONtoFwdTrack(AliESDMuonTrack &MUONTrack) {
  // Convert MCH Track parameters and covarianceS matrix to the RUN3 Forward track coordinate system

  FwdTrackPars convertedTrack;

  //pxdca
  double pdca;
  pdca= MUONTrack.P() * MUONTrack.GetDCA();

  // Parameter conversion
  double alpha1, alpha3, alpha4, x2, x3, x4;

  alpha1 = TMath::Tan(MUONTrack.GetThetaX());
  alpha3 = TMath::Tan(MUONTrack.GetThetaY());
  alpha4 = MUONTrack.GetInverseBendingMomentum();

  x2 = TMath::ATan2(-alpha3, -alpha1);
  if (alpha3 != 0 || alpha1 != 0)
    x3 = -1. / TMath::Sqrt(alpha3 * alpha3 + alpha1 * alpha1);
  else
    x3 = -FLT_MAX;
  x4 = alpha4 * -x3 * TMath::Sqrt(1 + alpha3 * alpha3);

  // Set output parameters
  convertedTrack.fX = AliMathBase::TruncateFloatFraction(MUONTrack.GetNonBendingCoor(), mMuonTrNonBend);
  convertedTrack.fY = AliMathBase::TruncateFloatFraction(MUONTrack.GetBendingCoor(), mMuonTrBend);
  convertedTrack.fZ = AliMathBase::TruncateFloatFraction(MUONTrack.GetZ(), mMuonTrZmu);
  convertedTrack.fPhi = AliMathBase::TruncateFloatFraction(x2, mMuonTrThetaX);
  convertedTrack.fTgl = AliMathBase::TruncateFloatFraction(x3, mMuonTrThetaX);
  convertedTrack.fSigned1Pt = AliMathBase::TruncateFloatFraction(x4, mMuonTr1P);
  convertedTrack.fChi2 = AliMathBase::TruncateFloatFraction(MUONTrack.GetChi2(), mMuonTrCov);
  convertedTrack.fChi2MatchMCHMID = AliMathBase::TruncateFloatFraction(MUONTrack.GetChi2MatchTrigger(), mMuonTrCov);
  convertedTrack.fRAtAbsorberEnd = AliMathBase::TruncateFloatFraction(MUONTrack.GetRAtAbsorberEnd(), mMuonTrCov);
  convertedTrack.fPDca = AliMathBase::TruncateFloatFraction(pdca, mMuonTrCov);

  //MCH and MID words
  convertedTrack.fMCHBitMap = MUONTrack.GetMuonClusterMap();
  UInt_t midpattern = MUONTrack.GetHitsPatternInTrigChTrk();
  UShort_t midbitmap = 0;
  for (int icath=0; icath < 2; icath++){
    for (int ich=0; ich < 4; ich++){
      if (AliESDMuonTrack::IsChamberHit(midpattern,icath,ich)) midbitmap |= 1<<(ich+icath*4);
    }
  }
  convertedTrack.fMIDBitMap = midbitmap;
  UInt_t midboard = static_cast<UInt_t>(AliESDMuonTrack::GetCrossedBoard(midpattern));
  convertedTrack.fMIDBoards = (midboard << 24) | (midboard << 16) | (midboard  << 8) | midboard;


  // Covariances matrix conversion
  using SMatrix55Std = ROOT::Math::SMatrix<double, 5>;
  using SMatrix55Sym = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
  SMatrix55Std jacobian;
  SMatrix55Sym convertedCovariances;
  TMatrixD inputCovariances;
  MUONTrack.GetCovariances(inputCovariances);

  auto K = alpha1 * alpha1 + alpha3 * alpha3;
  auto K32 = K * TMath::Sqrt(K);
  auto L = TMath::Sqrt(alpha3 * alpha3 + 1);

  convertedCovariances(0, 0) = inputCovariances(0, 0);
  convertedCovariances(0, 1) = inputCovariances(0, 1);
  convertedCovariances(0, 2) = inputCovariances(0, 2);
  convertedCovariances(0, 3) = inputCovariances(0, 3);
  convertedCovariances(0, 4) = inputCovariances(0, 4);

  convertedCovariances(1, 1) = inputCovariances(1, 1);
  convertedCovariances(1, 2) = inputCovariances(1, 2);
  convertedCovariances(1, 3) = inputCovariances(1, 3);
  convertedCovariances(1, 4) = inputCovariances(1, 4);

  convertedCovariances(2, 2) = inputCovariances(2, 2);
  convertedCovariances(2, 3) = inputCovariances(2, 3);
  convertedCovariances(2, 4) = inputCovariances(2, 4);

  convertedCovariances(3, 3) = inputCovariances(3, 3);
  convertedCovariances(3, 4) = inputCovariances(3, 4);

  convertedCovariances(4, 4) = inputCovariances(4, 4);

  jacobian(0, 0) = 1;

  jacobian(1, 2) = 1;

  if(K32 != 0) {
     jacobian(2, 1) = -alpha3 / K;
     jacobian(2, 3) = alpha1 / K;

     jacobian(3, 1) = alpha1 / K32;
     jacobian(3, 3) = alpha3 / K32;

     jacobian(4, 1) = -alpha1 * alpha4 * L / K32;
     jacobian(4, 3) = alpha3 * alpha4 * (1 / (TMath::Sqrt(K) * L) - L / K32);
     jacobian(4, 4) = L / TMath::Sqrt(K);
  }
  else
  {
     jacobian(2, 1) = 0;
     jacobian(2, 3) = 0;

     jacobian(3, 1) = 0;
     jacobian(3, 3) = 0;

     jacobian(4, 1) = 0;
     jacobian(4, 3) = 0;
     jacobian(4, 4) = 0;
  }

  // jacobian*covariances*jacobian^T
  convertedCovariances = ROOT::Math::Similarity(jacobian, convertedCovariances);

  // Set output covariances
  convertedTrack.fSigmaX   = AliMathBase::TruncateFloatFraction(TMath::Sqrt(convertedCovariances(0,0)), mTrackCovDiag);
  convertedTrack.fSigmaY   = AliMathBase::TruncateFloatFraction(TMath::Sqrt(convertedCovariances(1,1)), mTrackCovDiag);
  convertedTrack.fSigmaPhi = AliMathBase::TruncateFloatFraction(TMath::Sqrt(convertedCovariances(2,2)), mTrackCovDiag);
  convertedTrack.fSigmaTgl = AliMathBase::TruncateFloatFraction(TMath::Sqrt(convertedCovariances(3,3)), mTrackCovDiag);
  convertedTrack.fSigma1Pt = AliMathBase::TruncateFloatFraction(TMath::Sqrt(convertedCovariances(4,4)), mTrackCovDiag);

  if(fwdtracks.fSigmaX != 0 && fwdtracks.fSigmaY != 0)
      convertedTrack.fRhoXY = (Char_t)(128. * convertedCovariances(0,1) / fwdtracks.fSigmaX / fwdtracks.fSigmaY);
  else
     convertedTrack.fRhoXY = 0;

  if(fwdtracks.fSigmaPhi != 0 && fwdtracks.fSigmaX != 0)
     convertedTrack.fRhoPhiX = (Char_t)(128. * convertedCovariances(0,2) / fwdtracks.fSigmaPhi / fwdtracks.fSigmaX);
  else
     convertedTrack.fRhoPhiX = 0;

  if(fwdtracks.fSigmaPhi != 0 && fwdtracks.fSigmaY != 0)
     convertedTrack.fRhoPhiY = (Char_t)(128. * convertedCovariances(1,2) / fwdtracks.fSigmaPhi / fwdtracks.fSigmaY);
  else
     convertedTrack.fRhoPhiY = 0;

  if(fwdtracks.fSigmaTgl != 0 && fwdtracks.fSigmaX != 0)
     convertedTrack.fRhoTglX = (Char_t)(128. * convertedCovariances(3,0) / fwdtracks.fSigmaTgl / fwdtracks.fSigmaX);
  else
     convertedTrack.fRhoTglX = 0;

  if(fwdtracks.fSigmaTgl != 0 && fwdtracks.fSigmaY != 0)
     convertedTrack.fRhoTglY = (Char_t)(128. * convertedCovariances(3,1) / fwdtracks.fSigmaTgl / fwdtracks.fSigmaY);
  else
     convertedTrack.fRhoTglY = 0;

  if(fwdtracks.fSigmaTgl != 0 && fwdtracks.fSigmaPhi != 0)
     convertedTrack.fRhoTglPhi = (Char_t)(128. * convertedCovariances(3,2) / fwdtracks.fSigmaTgl / fwdtracks.fSigmaPhi);
  else
     convertedTrack.fRhoTglPhi = 0;

  if(fwdtracks.fSigma1Pt != 0 && fwdtracks.fSigmaX != 0)
     convertedTrack.fRho1PtX = (Char_t)(128. * convertedCovariances(4,0) / fwdtracks.fSigma1Pt / fwdtracks.fSigmaX);
  else
     convertedTrack.fRho1PtX = 0;

  if(fwdtracks.fSigma1Pt != 0 && fwdtracks.fSigmaY != 0)
     convertedTrack.fRho1PtY = (Char_t)(128. * convertedCovariances(4,1) / fwdtracks.fSigma1Pt / fwdtracks.fSigmaY);
  else
     convertedTrack.fRho1PtY = 0;

  if(fwdtracks.fSigma1Pt != 0 && fwdtracks.fSigmaPhi != 0)
     convertedTrack.fRho1PtPhi = (Char_t)(128. * convertedCovariances(4,2) / fwdtracks.fSigma1Pt / fwdtracks.fSigmaPhi);
  else
     convertedTrack.fRho1PtPhi = 0;

  if(fwdtracks.fSigma1Pt != 0 && fwdtracks.fSigmaTgl != 0)
    convertedTrack.fRho1PtTgl = (Char_t)(128. * convertedCovariances(4,3) / fwdtracks.fSigma1Pt / fwdtracks.fSigmaTgl);
  else
    convertedTrack.fRho1PtTgl = 0;

  return convertedTrack;
}

//_________________________________________________________________________________________________
AliAnalysisTaskAO2Dconverter::FwdTrackPars AliAnalysisTaskAO2Dconverter::MUONtoFwdTrack(AliAODTrack &MUONTrack) {
  // Convert MCH Track parameters and covarianceS matrix to the RUN3 Forward track coordinate system

  FwdTrackPars convertedTrack;

  //pxdca
  double pdca;
  pdca= MUONTrack.P() * MUONTrack.DCA();

  // Parameter conversion
  double alpha1, alpha3, alpha4, x2, x3, x4;

  alpha1 = TMath::Tan(TMath::ACos(MUONTrack.Px()/MUONTrack.P()));
  alpha3 = TMath::Tan(TMath::ACos(MUONTrack.Py()/MUONTrack.P()));
  alpha4 = MUONTrack.OneOverPt();

  x2 = TMath::ATan2(-alpha3, -alpha1);
  if (alpha3 != 0 || alpha1 != 0)
    x3 = -1. / TMath::Sqrt(alpha3 * alpha3 + alpha1 * alpha1);
  else
    x3 = 0;
  x4 = alpha4 * -x3 * TMath::Sqrt(1 + alpha3 * alpha3);

  // Set output parameters
  convertedTrack.fX = AliMathBase::TruncateFloatFraction(MUONTrack.Xv(), mMuonTrNonBend);
  convertedTrack.fY = AliMathBase::TruncateFloatFraction(MUONTrack.Yv(), mMuonTrBend);
  convertedTrack.fZ = AliMathBase::TruncateFloatFraction(MUONTrack.Zv(), mMuonTrZmu);
  convertedTrack.fPhi = AliMathBase::TruncateFloatFraction(x2, mMuonTrThetaX);
  convertedTrack.fTgl = AliMathBase::TruncateFloatFraction(x3, mMuonTrThetaX);
  convertedTrack.fSigned1Pt = AliMathBase::TruncateFloatFraction(x4, mMuonTr1P);
  convertedTrack.fChi2 = AliMathBase::TruncateFloatFraction(MUONTrack.Chi2perNDF(), mMuonTrCov);
  convertedTrack.fChi2MatchMCHMID = AliMathBase::TruncateFloatFraction(MUONTrack.GetChi2MatchTrigger(), mMuonTrCov);
  convertedTrack.fRAtAbsorberEnd = AliMathBase::TruncateFloatFraction(MUONTrack.GetRAtAbsorberEnd(), mMuonTrCov);
  convertedTrack.fPDca = AliMathBase::TruncateFloatFraction(pdca, mMuonTrCov);

  //MCH and MID words
  convertedTrack.fMCHBitMap = MUONTrack.GetMUONClusterMap();
  UInt_t midpattern = MUONTrack.GetMUONTrigHitsMapTrk();
  UShort_t midbitmap = 0;
  for (int icath=0; icath < 2; icath++){
    for (int ich=0; ich < 4; ich++){
      if (AliESDMuonTrack::IsChamberHit(midpattern,icath,ich)) midbitmap |= 1<<(ich+icath*4);
    }
  }
  convertedTrack.fMIDBitMap = midbitmap;
  UInt_t midboard = static_cast<UInt_t>(AliESDMuonTrack::GetCrossedBoard(midpattern));
  convertedTrack.fMIDBoards = (midboard << 24) | (midboard << 16) | (midboard  << 8) | midboard;

  // Covariances matrix conversion
  using SMatrix55Std = ROOT::Math::SMatrix<double, 5>;
  using SMatrix55Sym = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
  SMatrix55Std jacobian;
  SMatrix55Sym convertedCovariances;
  TMatrixD inputCovariances;
  MUONTrack.GetCovMatrix(&inputCovariances);

  auto K = alpha1 * alpha1 + alpha3 * alpha3;
  auto K32 = K * TMath::Sqrt(K);
  auto L = TMath::Sqrt(alpha3 * alpha3 + 1);

  convertedCovariances(0, 0) = inputCovariances(0, 0);
  convertedCovariances(0, 1) = inputCovariances(0, 1);
  convertedCovariances(0, 2) = inputCovariances(0, 2);
  convertedCovariances(0, 3) = inputCovariances(0, 3);
  convertedCovariances(0, 4) = inputCovariances(0, 4);

  convertedCovariances(1, 1) = inputCovariances(1, 1);
  convertedCovariances(1, 2) = inputCovariances(1, 2);
  convertedCovariances(1, 3) = inputCovariances(1, 3);
  convertedCovariances(1, 4) = inputCovariances(1, 4);

  convertedCovariances(2, 2) = inputCovariances(2, 2);
  convertedCovariances(2, 3) = inputCovariances(2, 3);
  convertedCovariances(2, 4) = inputCovariances(2, 4);

  convertedCovariances(3, 3) = inputCovariances(3, 3);
  convertedCovariances(3, 4) = inputCovariances(3, 4);

  convertedCovariances(4, 4) = inputCovariances(4, 4);

  jacobian(0, 0) = 1;

  jacobian(1, 2) = 1;

  if(K32 != 0) {
     jacobian(2, 1) = -alpha3 / K;
     jacobian(2, 3) = alpha1 / K;

     jacobian(3, 1) = alpha1 / K32;
     jacobian(3, 3) = alpha3 / K32;

     jacobian(4, 1) = -alpha1 * alpha4 * L / K32;
     jacobian(4, 3) = alpha3 * alpha4 * (1 / (TMath::Sqrt(K) * L) - L / K32);
     jacobian(4, 4) = L / TMath::Sqrt(K);
  }
  else
  {
     jacobian(2, 1) = 0;
     jacobian(2, 3) = 0;

     jacobian(3, 1) = 0;
     jacobian(3, 3) = 0;

     jacobian(4, 1) = 0;
     jacobian(4, 3) = 0;
     jacobian(4, 4) = 0;
  }

  // jacobian*covariances*jacobian^T
  convertedCovariances = ROOT::Math::Similarity(jacobian, convertedCovariances);

  // Set output covariances
  convertedTrack.fSigmaX   = AliMathBase::TruncateFloatFraction(TMath::Sqrt(convertedCovariances(0,0)), mTrackCovDiag);
  convertedTrack.fSigmaY   = AliMathBase::TruncateFloatFraction(TMath::Sqrt(convertedCovariances(1,1)), mTrackCovDiag);
  convertedTrack.fSigmaPhi = AliMathBase::TruncateFloatFraction(TMath::Sqrt(convertedCovariances(2,2)), mTrackCovDiag);
  convertedTrack.fSigmaTgl = AliMathBase::TruncateFloatFraction(TMath::Sqrt(convertedCovariances(3,3)), mTrackCovDiag);
  convertedTrack.fSigma1Pt = AliMathBase::TruncateFloatFraction(TMath::Sqrt(convertedCovariances(4,4)), mTrackCovDiag);

  if(fwdtracks.fSigmaX != 0 && fwdtracks.fSigmaY != 0)
      convertedTrack.fRhoXY = (Char_t)(128. * convertedCovariances(0,1) / fwdtracks.fSigmaX / fwdtracks.fSigmaY);
  else
     convertedTrack.fRhoXY = 0;

  if(fwdtracks.fSigmaPhi != 0 && fwdtracks.fSigmaX != 0)
     convertedTrack.fRhoPhiX = (Char_t)(128. * convertedCovariances(0,2) / fwdtracks.fSigmaPhi / fwdtracks.fSigmaX);
  else
     convertedTrack.fRhoPhiX = 0;

  if(fwdtracks.fSigmaPhi != 0 && fwdtracks.fSigmaY != 0)
     convertedTrack.fRhoPhiY = (Char_t)(128. * convertedCovariances(1,2) / fwdtracks.fSigmaPhi / fwdtracks.fSigmaY);
  else
     convertedTrack.fRhoPhiY = 0;

  if(fwdtracks.fSigmaTgl != 0 && fwdtracks.fSigmaX != 0)
     convertedTrack.fRhoTglX = (Char_t)(128. * convertedCovariances(3,0) / fwdtracks.fSigmaTgl / fwdtracks.fSigmaX);
  else
     convertedTrack.fRhoTglX = 0;

  if(fwdtracks.fSigmaTgl != 0 && fwdtracks.fSigmaY != 0)
     convertedTrack.fRhoTglY = (Char_t)(128. * convertedCovariances(3,1) / fwdtracks.fSigmaTgl / fwdtracks.fSigmaY);
  else
     convertedTrack.fRhoTglY = 0;

  if(fwdtracks.fSigmaTgl != 0 && fwdtracks.fSigmaPhi != 0)
     convertedTrack.fRhoTglPhi = (Char_t)(128. * convertedCovariances(3,2) / fwdtracks.fSigmaTgl / fwdtracks.fSigmaPhi);
  else
     convertedTrack.fRhoTglPhi = 0;

  if(fwdtracks.fSigma1Pt != 0 && fwdtracks.fSigmaX != 0)
     convertedTrack.fRho1PtX = (Char_t)(128. * convertedCovariances(4,0) / fwdtracks.fSigma1Pt / fwdtracks.fSigmaX);
  else
     convertedTrack.fRho1PtX = 0;

  if(fwdtracks.fSigma1Pt != 0 && fwdtracks.fSigmaY != 0)
     convertedTrack.fRho1PtY = (Char_t)(128. * convertedCovariances(4,1) / fwdtracks.fSigma1Pt / fwdtracks.fSigmaY);
  else
     convertedTrack.fRho1PtY = 0;

  if(fwdtracks.fSigma1Pt != 0 && fwdtracks.fSigmaPhi != 0)
     convertedTrack.fRho1PtPhi = (Char_t)(128. * convertedCovariances(4,2) / fwdtracks.fSigma1Pt / fwdtracks.fSigmaPhi);
  else
     convertedTrack.fRho1PtPhi = 0;

  if(fwdtracks.fSigma1Pt != 0 && fwdtracks.fSigmaTgl != 0)
    convertedTrack.fRho1PtTgl = (Char_t)(128. * convertedCovariances(4,3) / fwdtracks.fSigma1Pt / fwdtracks.fSigmaTgl);
  else
    convertedTrack.fRho1PtTgl = 0;

  return convertedTrack;
}

//_________________________________________________________________________________________________
AliAnalysisTaskAO2Dconverter::MCGeneratorID AliAnalysisTaskAO2Dconverter::GetMCGeneratorID(AliGenEventHeader * genHeader) {
  // Not very elegant way of finding the generator's ID
  if (genHeader->InheritsFrom("AliGenDPMjetEventHeader")) {
    return kAliGenDPMjetEventHeader;
  }
  if (genHeader->InheritsFrom("AliGenEpos3EventHeader")) {
    return kAliGenEpos3EventHeader;
  }
  if (genHeader->InheritsFrom("AliGenEposEventHeader")) {
    return kAliGenEposEventHeader;
  }
  if (genHeader->InheritsFrom("AliGenEventHeaderTunedPbPb")) {
    return kAliGenEventHeaderTunedPbPb;
  }
  if (genHeader->InheritsFrom("AliGenGeVSimEventHeader")) {
    return kAliGenGeVSimEventHeader;
  }
  if (genHeader->InheritsFrom("AliGenHepMCEventHeader")) {
    return kAliGenHepMCEventHeader;
  }
  if (genHeader->InheritsFrom("AliGenHerwigEventHeader")) {
    return kAliGenHerwigEventHeader;
  }
  if (genHeader->InheritsFrom("AliGenHijingEventHeader")) {
    return kAliGenHijingEventHeader;
  }
  if (genHeader->InheritsFrom("AliGenPythiaEventHeader")) {
    return kAliGenPythiaEventHeader;
  }
  if (genHeader->InheritsFrom("AliGenToyEventHeader")) {
    return kAliGenToyEventHeader;
  }
  return kAliGenCocktailEventHeader;
}
//_________________________________________________________________________________________________
void AliAnalysisTaskAO2Dconverter::CrossSections::Print() {
  std::cout << "Content of CrossSections --------------->" << std::endl;
  std::cout << "fGeneratorsID = " << fGeneratorsID << std::endl;
  std::cout << "fAccepted = " << fAccepted << std::endl;
  std::cout << "fAttempted = " << fAttempted << std::endl;
  std::cout << "fXsectGen = " << fXsectGen << std::endl;
  std::cout << "fXsectErr = " << fXsectErr << std::endl;
  std::cout << "fPtHard = " << fPtHard << std::endl;
  std::cout << "<---------------" << std::endl;
}

//_________________________________________________________________________________________________
void AliAnalysisTaskAO2Dconverter::CrossSections::Reset() {
  fIndexMcCollisions = -1;
  fGeneratorsID = 0u;
  fAccepted = 0;
  fAttempted = 0;
  fXsectGen = -999.f;
  fXsectErr = -999.f;
  fPtHard = -999.f;
}

//_________________________________________________________________________________________________
void AliAnalysisTaskAO2Dconverter::CrossSections::Fill(AliGenEventHeader * genHeader, AliAnalysisTaskAO2Dconverter * task) {
  if (!genHeader) {
    // Protection, normally genHeader should always be set
    std::cerr << "No AliGenEventHeader!" << std::endl;
    return;
  }
  // Define the number of accepted events from the map
  const std::map<std::string, Float_t> weights = genHeader->GetEventWeights();
  fAccepted = weights.size();

  fGeneratorsID = AliAnalysisTaskAO2Dconverter::GetMCGeneratorID(genHeader);
  switch (fGeneratorsID) {
  case kAliGenDPMjetEventHeader: {
    AliGenDPMjetEventHeader * hepMcHeader = (AliGenDPMjetEventHeader*)genHeader;
    fAttempted = hepMcHeader->Trials();
    fXsectGen = -999.f;
    fXsectErr = -999.f;
    fPtHard = -999.f;
    break;
  }
  case kAliGenHepMCEventHeader: {
    AliGenHepMCEventHeader * hepMcHeader = (AliGenHepMCEventHeader*)genHeader;
    fAttempted = hepMcHeader->ntrials();
    fXsectGen = hepMcHeader->sigma_gen();
    fXsectErr = hepMcHeader->sigma_err();
    fPtHard = hepMcHeader->pthard();
    break;
  }
  case kAliGenHerwigEventHeader: {
    AliGenHerwigEventHeader * hepMcHeader = (AliGenHerwigEventHeader*)genHeader;
    fAttempted = hepMcHeader->Trials();
    fXsectGen = -999.f;
    fXsectErr = -999.f;
    fPtHard = hepMcHeader->GetPtHard();
    break;
  }
  case kAliGenHijingEventHeader: {
    AliGenHijingEventHeader * hepMcHeader = (AliGenHijingEventHeader*)genHeader;
    fAttempted = hepMcHeader->Trials();
    fXsectGen = -999.f;
    fXsectErr = -999.f;
    fPtHard = -999.f;
    break;
  }
  case kAliGenPythiaEventHeader: {
    AliGenPythiaEventHeader * hepMcHeader = (AliGenPythiaEventHeader*)genHeader;
    fAttempted = hepMcHeader->Trials();
    fXsectGen = hepMcHeader->GetXsection();
    fXsectErr = -999;
    fPtHard = hepMcHeader->GetPtHard();
    break;
  }
  default: {
    fAttempted = 0;
    fXsectGen = -999.f;
    fXsectErr = -999.f;
    fPtHard = -999.f;
    break;
  }
  }
  task->FillTree(kHepMcCrossSections);

  // Special case - AliGenCocktailEventHeader
  if (fGeneratorsID == kAliGenCocktailEventHeader) {
    TList* genList = ((AliGenCocktailEventHeader*)genHeader)->GetHeaders();
    TListIter next(genList);
    TObject * to;
    while((to=next())) {
      Fill((AliGenEventHeader *)to, task);
    }
  }
}

//_________________________________________________________________________________________________
void AliAnalysisTaskAO2Dconverter::PdfInfo::Print() {
  std::cout << "Content of PdfInfo --------------->" << std::endl;
  std::cout << "fGeneratorsID = " << fGeneratorsID << std::endl;
  std::cout << "fId1 = " << fId1 << std::endl;
  std::cout << "fId2 = " << fId2 << std::endl;
  std::cout << "fPdfId1 = " << fPdfId1 << std::endl;
  std::cout << "fPdfId2 = " << fPdfId2 << std::endl;
  std::cout << "fX1 = " << fX1 << std::endl;
  std::cout << "fX2 = " << fX2 << std::endl;
  std::cout << "fScalePdf = " << fScalePdf << std::endl;
  std::cout << "fPdf1 = " << fPdf1 << std::endl;
  std::cout << "fPdf2 = " << fPdf2 << std::endl;
  std::cout << "<---------------" << std::endl;
}

//_________________________________________________________________________________________________
void AliAnalysisTaskAO2Dconverter::PdfInfo::Reset(){
  fIndexMcCollisions = -1;
  fGeneratorsID = 0u;
  fId1 = 0;
  fId2 = 0;
  fPdfId1 = 0;
  fPdfId2 = 0;
  fX1 = 0.f;
  fX2 = 0.f;
  fScalePdf = 0.f;
  fPdf1 = 0.f;
  fPdf2 = 0.f;
}

//_________________________________________________________________________________________________
void AliAnalysisTaskAO2Dconverter::PdfInfo::Fill(AliGenEventHeader * genHeader, AliAnalysisTaskAO2Dconverter * task) {
  if (!genHeader) {
    // Protection, normally genHeader should always be set
    std::cerr << "No AliGenEventHeader!" << std::endl;
    return;
  }
  fGeneratorsID = AliAnalysisTaskAO2Dconverter::GetMCGeneratorID(genHeader);
  if (fGeneratorsID == kAliGenHepMCEventHeader) {
    AliGenHepMCEventHeader * hepMcHeader = (AliGenHepMCEventHeader*)genHeader;

    fId1 = hepMcHeader->id1();
    fId2 = hepMcHeader->id2();
    fPdfId1 = hepMcHeader->pdf_id1();
    fPdfId2 = hepMcHeader->pdf_id2();
    fX1 = hepMcHeader->x1();
    fX2 = hepMcHeader->x2();
    fScalePdf = hepMcHeader->scalePDF();
    fPdf1 = hepMcHeader->pdf1();
    fPdf2 = hepMcHeader->pdf2();
    task->FillTree(kHepMcPdfInfo);
    return;
  }
  // Special case - AliGenCocktailEventHeader
  if (fGeneratorsID == kAliGenCocktailEventHeader) {
    TList* genList = ((AliGenCocktailEventHeader*)genHeader)->GetHeaders();
    TListIter next(genList);
    TObject * to;
    while((to=next())) {
      Fill((AliGenEventHeader *)to, task);
    }
  }
}
//_________________________________________________________________________________________________
void AliAnalysisTaskAO2Dconverter::HeavyIon::Print() {
  std::cout << "Content of HeavyIon --------------->" << std::endl;
  std::cout << "fGeneratorsID = " << fGeneratorsID << std::endl;
  std::cout << "fNcollHard = " << fNcollHard << std::endl;
  std::cout << "fNpartProj = " << fNpartProj << std::endl;
  std::cout << "fNpartTarg = " << fNpartTarg << std::endl;
  std::cout << "fNcoll = " << fNcoll << std::endl;
  std::cout << "fNNwoundedCollisions = " << fNNwoundedCollisions << std::endl;
  std::cout << "fNwoundedNCollisions = " << fNwoundedNCollisions << std::endl;
  std::cout << "fNwoundedNwoundedCollisions = " << fNwoundedNwoundedCollisions << std::endl;
  std::cout << "fSpectatorNeutrons = " << fSpectatorNeutrons << std::endl;
  std::cout << "fSpectatorProtons = " << fSpectatorProtons << std::endl;
  std::cout << "fImpactParameter = " << fImpactParameter << std::endl;
  std::cout << "fEventPlaneAngle = " << fEventPlaneAngle << std::endl;
  std::cout << "fEccentricity = " << fEccentricity << std::endl;
  std::cout << "fSigmaInelNN = " << fSigmaInelNN << std::endl;
  std::cout << "fCentrality = " << fCentrality << std::endl;
  std::cout << "<---------------" << std::endl;
}

//_________________________________________________________________________________________________
void AliAnalysisTaskAO2Dconverter::HeavyIon::Reset(){
  fIndexMcCollisions = -1;
  fGeneratorsID = 0u;
  fNcollHard = -999;
  fNpartProj = -999;
  fNpartTarg = -999;
  fNcoll = -999;
  fNNwoundedCollisions = -999;
  fNwoundedNCollisions = -999;
  fNwoundedNwoundedCollisions = -999;
  fSpectatorNeutrons = -999;
  fSpectatorProtons = -999;
  fImpactParameter = -999.f;
  fEventPlaneAngle = 0.f;
  fEccentricity = 0.f;
  fSigmaInelNN = -999.f;
  fCentrality = -999.f;
}

//_________________________________________________________________________________________________
void AliAnalysisTaskAO2Dconverter::HeavyIon::Fill(AliGenEventHeader * genHeader, AliAnalysisTaskAO2Dconverter * task) {
  if (!genHeader) {
    // Protection, normally genHeader should always be set
    std::cerr << "No AliGenEventHeader!" << std::endl;
    return;
  }
  fGeneratorsID = AliAnalysisTaskAO2Dconverter::GetMCGeneratorID(genHeader);
  if (fGeneratorsID == kAliGenHepMCEventHeader) {
    // The information is already available
    AliGenHepMCEventHeader * hepMcHeader = (AliGenHepMCEventHeader*)genHeader;

    fNcollHard = hepMcHeader->Ncoll_hard();
    fNpartProj = hepMcHeader->Npart_proj();
    fNpartTarg = hepMcHeader->Npart_targ();
    fNcoll = hepMcHeader->Ncoll();
    fNNwoundedCollisions = hepMcHeader->N_Nwounded_collisions();
    fNwoundedNCollisions = hepMcHeader->Nwounded_N_collisions();
    fNwoundedNwoundedCollisions = hepMcHeader->Nwounded_Nwounded_collisions();
    fSpectatorNeutrons = hepMcHeader->spectator_neutrons();
    fSpectatorProtons = hepMcHeader->spectator_protons();
    fImpactParameter = hepMcHeader->impact_parameter();
    fEventPlaneAngle = hepMcHeader->event_plane_angle();
    fEccentricity = hepMcHeader->eccentricity();
    fSigmaInelNN = hepMcHeader->sigma_inel_NN();
    //PH fCentrality is not in HepMC
    task->FillTree(kHepMcHeavyIon);
    return;
  }
  if (genHeader->InheritsFrom("AliCollisionGeometry")) {
    // Extracting form AliCollisionInfo
    // Covers: AliGenDPMjetEventHeader, AliGenEpos3EventHeader,
    //         AliGenEposEventHeader, AliGenHijingEventHeader,
    //         AliGenHydjetEventHeader
    AliCollisionGeometry * hepMcHeader = (AliCollisionGeometry*)genHeader;

    fNcollHard = hepMcHeader->HardScatters();
    fNpartProj = hepMcHeader->ProjectileParticipants();
    fNpartTarg = hepMcHeader->TargetParticipants();
    fNcoll = hepMcHeader->NN();
    fNNwoundedCollisions = hepMcHeader->NNw();
    fNwoundedNCollisions = hepMcHeader->NwN();
    fNwoundedNwoundedCollisions = hepMcHeader->NwNw();
    fSpectatorNeutrons = hepMcHeader->ProjSpectatorsn() + hepMcHeader->TargSpectatorsn();
    fSpectatorProtons = hepMcHeader->ProjSpectatorsp() + hepMcHeader->TargSpectatorsp();
    fImpactParameter = hepMcHeader->ImpactParameter();
    fEventPlaneAngle = hepMcHeader->ReactionPlaneAngle();

    //PH fEccentricity, fSigmaInelNN, fCentrallity are not in AliCollisionGeometry
    task->FillTree(kHepMcHeavyIon);
    return;
  }
  if (fGeneratorsID == kAliGenEventHeaderTunedPbPb) {
    AliGenEventHeaderTunedPbPb * hepMcHeader = (AliGenEventHeaderTunedPbPb*)genHeader;
    fCentrality = hepMcHeader->GetCentrality();
    task->FillTree(kHepMcHeavyIon);
    return;
  }
  // Special case - AliGenCocktailEventHeader
  if (fGeneratorsID == kAliGenCocktailEventHeader) {
    TList* genList = ((AliGenCocktailEventHeader*)genHeader)->GetHeaders();
    TListIter next(genList);
    TObject * to;
    while((to=next())) {
      Fill((AliGenEventHeader *)to, task);
    }
  }
}

////////////////////////////////////////////////////////////
