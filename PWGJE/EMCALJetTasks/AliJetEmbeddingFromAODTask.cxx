// $Id: AliJetEmbeddingFromAODTask.cxx $
//
// Jet embedding from AOD task.
//
// Author: S.Aiola, C.Loizides

#include "AliJetEmbeddingFromAODTask.h"

// C++ standard library
#include <vector>

// ROOT
#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TGrid.h>
#include <TH2C.h>
#include <TList.h>
#include <TStreamerInfo.h>
#include <TRandom.h>
#include <TSystem.h>
#include <TLorentzVector.h>

// AliRoot
#include "AliVEvent.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliPicoTrack.h"
#include "AliVTrack.h"
#include "AliVCluster.h"
#include "AliVCaloCells.h"
#include "AliAODMCParticle.h"
#include "AliCentrality.h"
#include "AliVHeader.h"
#include "AliVVertex.h"
#include "AliAODHeader.h"
#include "AliFJWrapper.h"
#include "AliLog.h"
#include "AliInputEventHandler.h"

ClassImp(AliJetEmbeddingFromAODTask)

//________________________________________________________________________
AliJetEmbeddingFromAODTask::AliJetEmbeddingFromAODTask() : 
  AliJetModelBaseTask("AliJetEmbeddingFromAODTask"),
  fFileList(0),
  fRandomAccess(kFALSE),
  fAODTreeName(),
  fAODHeaderName(),
  fAODVertexName(),
  fAODTrackName(),
  fAODClusName(),
  fAODCellsName(),
  fAODMCParticlesName(),
  fMinCentrality(0),
  fMaxCentrality(10),
  fTriggerMask(AliVEvent::kAny),
  fZVertexCut(10),
  fJetMinPt(0),
  fJetMinEta(-0.5),
  fJetMaxEta(0.5),
  fJetMinPhi(-999),
  fJetMaxPhi(999),
  fJetConstituentMinPt(0),
  fJetRadius(0.4),
  fJetType(0),
  fJetAlgo(1),
  fJetParticleLevel(kTRUE),
  fIncludeNoITS(kTRUE),
  fUseNegativeLabels(kTRUE),
  fTrackEfficiency(1),
  fIsAODMC(kFALSE),
  fTotalFiles(2050),
  fAttempts(5),
  fEsdTreeMode(kFALSE),
  fCurrentFileID(0),
  fCurrentAODFileID(0),
  fCurrentAODFile(0),
  fPicoTrackVersion(0),
  fCurrentAODTree(0),
  fAODHeader(0),
  fAODVertex(0),
  fAODTracks(0),
  fAODClusters(0),
  fAODCaloCells(0),
  fAODMCParticles(0),
  fCurrentAODEntry(-1),
  fHistFileMatching(0),
  fHistAODFileError(0),
  fHistNotEmbedded(0),
  fHistEmbeddingQA(0)
{
  // Default constructor.
  SetSuffix("AODEmbedding");
  SetMarkMC(0);
  fAODfilterBits[0] = -1;
  fAODfilterBits[1] = -1;
  fEtaMin = -1;
  fEtaMax = 1;
  fPhiMin = -10;
  fPhiMax = 10;
  fPtMin = 0;
  fPtMax = 1000;
}

//________________________________________________________________________
AliJetEmbeddingFromAODTask::AliJetEmbeddingFromAODTask(const char *name, Bool_t drawqa) : 
  AliJetModelBaseTask(name, drawqa),
  fFileList(0),
  fRandomAccess(kFALSE),
  fAODTreeName("aodTree"),
  fAODHeaderName("header"),
  fAODVertexName("vertices"),
  fAODTrackName("tracks"),
  fAODClusName(""),
  fAODCellsName("emcalCells"),
  fAODMCParticlesName(AliAODMCParticle::StdBranchName()),
  fMinCentrality(0),
  fMaxCentrality(10),
  fTriggerMask(AliVEvent::kAny),
  fZVertexCut(10),
  fJetMinPt(0),
  fJetMinEta(-0.5),
  fJetMaxEta(0.5),
  fJetMinPhi(-999),
  fJetMaxPhi(999),
  fJetConstituentMinPt(0),
  fJetRadius(0.4),
  fJetType(0),
  fJetAlgo(1),
  fJetParticleLevel(kTRUE),
  fIncludeNoITS(kTRUE),
  fUseNegativeLabels(kTRUE),
  fTrackEfficiency(1),
  fIsAODMC(kFALSE),
  fTotalFiles(2050),
  fAttempts(5),
  fEsdTreeMode(kFALSE),
  fCurrentFileID(0),
  fCurrentAODFileID(0),
  fCurrentAODFile(0),
  fPicoTrackVersion(0),
  fCurrentAODTree(0),
  fAODHeader(0),
  fAODVertex(0),
  fAODTracks(0),
  fAODClusters(0),
  fAODCaloCells(0),  
  fAODMCParticles(0),
  fCurrentAODEntry(-1),
  fHistFileMatching(0),
  fHistAODFileError(0),
  fHistNotEmbedded(0),
  fHistEmbeddingQA(0)
{
  // Standard constructor.
  SetSuffix("AODEmbedding");
  SetMarkMC(0);
  fAODfilterBits[0] = -1;
  fAODfilterBits[1] = -1;
  fEtaMin = -1;
  fEtaMax = 1;
  fPhiMin = -10;
  fPhiMax = 10;
  fPtMin = 0;
  fPtMax = 1000;
}

//________________________________________________________________________
AliJetEmbeddingFromAODTask::~AliJetEmbeddingFromAODTask()
{
  // Destructor

  if (fCurrentAODFile) {
    fCurrentAODFile->Close();
    delete fCurrentAODFile;
  }
}

//________________________________________________________________________
void AliJetEmbeddingFromAODTask::UserCreateOutputObjects()
{
  if (!fQAhistos)
    return;

  AliJetModelBaseTask::UserCreateOutputObjects();
  
  fHistFileMatching = new TH2C("fHistFileMatching", "fHistFileMatching", fTotalFiles, -0.5, fTotalFiles - 0.5, fFileList->GetEntriesFast(), -0.5, fFileList->GetEntriesFast() -0.5);
  fHistFileMatching->GetXaxis()->SetTitle("File no. (PYTHIA)");
  fHistFileMatching->GetYaxis()->SetTitle("File no. (Embedded AOD)");
  fHistFileMatching->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistFileMatching);

  fHistAODFileError = new TH1C("fHistAODFileError", "fHistAODFileError", fFileList->GetEntriesFast(), -0.5, fFileList->GetEntriesFast() -0.5);
  fHistAODFileError->GetXaxis()->SetTitle("File no. (Embedded AOD)");
  fHistAODFileError->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistAODFileError);

  fHistNotEmbedded = new TH1C("fHistNotEmbedded", "fHistNotEmbedded", fTotalFiles, -0.5, fTotalFiles - 0.5);
  fHistNotEmbedded->GetXaxis()->SetTitle("File no. (PYTHIA)");
  fHistNotEmbedded->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistNotEmbedded);

  fHistEmbeddingQA = new TH1F("fHistEmbeddingQA", "fHistEmbeddingQA", 2, 0, 2);
  fHistEmbeddingQA->GetXaxis()->SetTitle("Event state");
  fHistEmbeddingQA->GetYaxis()->SetTitle("counts");
  fHistEmbeddingQA->GetXaxis()->SetBinLabel(1, "OK");
  fHistEmbeddingQA->GetXaxis()->SetBinLabel(2, "Not embedded");
  fOutput->Add(fHistEmbeddingQA);

  PostData(1, fOutput);
}

//________________________________________________________________________
Bool_t AliJetEmbeddingFromAODTask::UserNotify()
{
  TString path(fInputHandler->GetTree()->GetCurrentFile()->GetPath());
  if (path.EndsWith("/"))
    path.Remove(path.Length()-1);
  path.Remove(path.Last('/'));
  path.Remove(0, path.Last('/')+1);
  if (!path.IsDec()) {
    AliWarning(Form("Could not extract file number from path %s", path.Data()));
    return kTRUE;
  }
  fCurrentFileID = path.Atoi();
  if (!fRandomAccess) {
    fCurrentAODFileID = fFileList->GetEntriesFast() * fCurrentFileID / fTotalFiles-1;
    AliInfo(Form("Start embedding from file ID %d", fCurrentAODFileID));
  }
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliJetEmbeddingFromAODTask::ExecOnce() 
{
  if (fAODTreeName.Contains("aod"))
    fEsdTreeMode = kFALSE;
  else
    fEsdTreeMode = kTRUE;

  return AliJetModelBaseTask::ExecOnce();
}

//________________________________________________________________________
Bool_t AliJetEmbeddingFromAODTask::OpenNextFile() 
{
  if (fCurrentAODFile) {
    fCurrentAODFile->Close();
    delete fCurrentAODFile;
    fCurrentAODFile = 0;
  }

  Int_t i = 0;

  while ((!fCurrentAODFile || fCurrentAODFile->IsZombie()) && i < fAttempts) {
    if (i > 0 && fHistAODFileError) {
	fHistAODFileError->Fill(fCurrentAODFileID);
    }

    fCurrentAODFile = GetNextFile();
    i++;
  } 

  if (!fCurrentAODFile || fCurrentAODFile->IsZombie())
    return kFALSE;

  const TList *clist = fCurrentAODFile->GetStreamerInfoCache();
  if(clist) {
    TStreamerInfo *cinfo = static_cast<TStreamerInfo*>(clist->FindObject("AliPicoTrack"));
    if(cinfo) 
      fPicoTrackVersion = cinfo->GetClassVersion();
    else
      fPicoTrackVersion = 0;
  }

  fCurrentAODTree = static_cast<TTree*>(fCurrentAODFile->Get(fAODTreeName));
  if (!fCurrentAODTree)
    return kFALSE;

  if (!fAODHeaderName.IsNull()) 
    fCurrentAODTree->SetBranchAddress(fAODHeaderName, &fAODHeader);
  
  if (!fAODVertexName.IsNull()) 
    fCurrentAODTree->SetBranchAddress(fAODVertexName, &fAODVertex);
      
  if (!fAODTrackName.IsNull()) 
    fCurrentAODTree->SetBranchAddress(fAODTrackName, &fAODTracks);
  
  if (!fAODClusName.IsNull()) 
    fCurrentAODTree->SetBranchAddress(fAODClusName, &fAODClusters);
  
  if (!fAODCellsName.IsNull()) 
    fCurrentAODTree->SetBranchAddress(fAODCellsName, &fAODCaloCells);
  
  if (!fAODMCParticlesName.IsNull()) 
    fCurrentAODTree->SetBranchAddress(fAODMCParticlesName, &fAODMCParticles);
  
  if (fRandomAccess)
    fCurrentAODEntry = TMath::Nint(gRandom->Rndm()*fCurrentAODTree->GetEntries());
  else
    fCurrentAODEntry = 0;

  AliDebug(3,Form("Will start embedding from entry %d", fCurrentAODEntry));
  
  if (fHistFileMatching)
    fHistFileMatching->Fill(fCurrentFileID, fCurrentAODFileID-1);
  
  return kTRUE;
}

//________________________________________________________________________
TFile* AliJetEmbeddingFromAODTask::GetNextFile()
{
  if (fRandomAccess) 
    fCurrentAODFileID = TMath::Nint(gRandom->Rndm()*fFileList->GetEntriesFast());
  else
    fCurrentAODFileID++;
  
  if (fCurrentAODFileID >= fFileList->GetEntriesFast()) {
    AliError("No more file in the list!");
    return 0;
  }
  
  TObjString *objFileName = static_cast<TObjString*>(fFileList->At(fCurrentAODFileID));
  TString fileName(objFileName->GetString());

  if (fileName.BeginsWith("alien://") && !gGrid) {
    AliInfo("Trying to connect to AliEn ...");
    TGrid::Connect("alien://");
  }

  TString baseFileName(fileName);
  if (baseFileName.Contains(".zip#")) {
    Ssiz_t pos = baseFileName.Last('#');
    baseFileName.Remove(pos);
  }
  
  if (gSystem->AccessPathName(baseFileName)) {
    AliDebug(3,Form("File %s does not exist!", baseFileName.Data()));
    return 0;
  }

  AliDebug(3,Form("Trying to open file %s...", fileName.Data()));
  TFile *file = TFile::Open(fileName);

  if (!file)
    AliDebug(3,Form("Unable to open file: %s!", fileName.Data()));

  return file;
}

//________________________________________________________________________
Bool_t AliJetEmbeddingFromAODTask::GetNextEntry() 
{
  Int_t attempts = 0;

  do {
    if (!fCurrentAODFile || !fCurrentAODTree || fCurrentAODEntry >= fCurrentAODTree->GetEntries()) {
      if (!OpenNextFile())
	return kFALSE;
    }
    
    fCurrentAODTree->GetEntry(fCurrentAODEntry);
    fCurrentAODEntry++;
    attempts++;

    if (attempts == 1000) 
      AliWarning("After 1000 attempts no event has been accepted by the event selection (trigger, centrality...)!");

  } while (!IsAODEventSelected() && fCurrentAODTree);

  return (fCurrentAODTree!=0);
}

//________________________________________________________________________
Bool_t AliJetEmbeddingFromAODTask::IsAODEventSelected()
{
  // AOD event selection.
  
  if (!fEsdTreeMode && fAODHeader) {
    AliAODHeader *aodHeader = static_cast<AliAODHeader*>(fAODHeader);

    // Trigger selection
    if (fTriggerMask != AliVEvent::kAny) {
      UInt_t offlineTrigger = aodHeader->GetOfflineTrigger();
      
      if ((offlineTrigger & fTriggerMask) == 0) {
	AliDebug(2,Form("Event rejected due to physics selection. Event trigger mask: %d, trigger mask selection: %d.", 
			offlineTrigger, fTriggerMask));
	return kFALSE;
      }
    }
    
    // Centrality selection
    if (fMinCentrality >= 0) {
      AliCentrality *cent = aodHeader->GetCentralityP();
      Float_t centVal = cent->GetCentralityPercentile("V0M");
      if (centVal < fMinCentrality || centVal >= fMaxCentrality) {
	AliDebug(2,Form("Event rejected due to centrality selection. Event centrality: %f, centrality range selection: %f to %f", 
			centVal, fMinCentrality, fMaxCentrality));
	return kFALSE;
      }
    }
  }

  // Vertex selection
  if (fAODVertex) {
    Double_t vert[3]={0};
    ((AliVVertex*)fAODVertex->At(0))->GetXYZ(vert);
    if (TMath::Abs(vert[2]) > fZVertexCut) {
      AliDebug(2,Form("Event rejected due to Z vertex selection. Event Z vertex: %f, Z vertex cut: %f", 
		      vert[2], fZVertexCut));
      return kFALSE;
    }
  }

  // Jet selection
  if (fJetMinPt > 0) {
    TLorentzVector jet;

    if (fJetParticleLevel) {
      if (fAODMCParticles)
	jet = GetLeadingJet(fAODMCParticles);
      else {
	AliWarning("Particle level jets selected, but not MC particles found. The jet event selection will be skipped.");
	return kTRUE;
      }
    }
    else {
      if (fAODTracks || fAODClusters) 
	jet = GetLeadingJet(fAODTracks, fAODClusters);
      else {
	AliWarning("Detector level jets selected, but not tracks or clusters found. The jet event selection will be skipped.");
	return kTRUE;
      }
    }
    
    AliDebug(1, Form("Leading jet pt = %f", jet.Pt()));
    if (jet.Pt() < fJetMinPt)
      return kFALSE;
  }

  return kTRUE;
}

//________________________________________________________________________
void AliJetEmbeddingFromAODTask::Run() 
{
  if (!GetNextEntry()) {
    if (fHistNotEmbedded)
      fHistNotEmbedded->Fill(fCurrentFileID);
    if (fHistEmbeddingQA)
      fHistEmbeddingQA->Fill("Not embedded", 1);
    AliError("Unable to get the AOD event to embed. Nothing will be embedded.");
    return;
  }

  if (fHistEmbeddingQA)
    fHistEmbeddingQA->Fill("OK", 1);

  if (fOutMCParticles) {

    if (fCopyArray && fMCParticles)
      CopyMCParticles();

    if (fAODMCParticles) {
      AliDebug(3, Form("%d MC particles will be processed for embedding.", fAODMCParticles->GetEntriesFast()));
      for (Int_t i = 0; i < fAODMCParticles->GetEntriesFast(); i++) {
	AliAODMCParticle *part = static_cast<AliAODMCParticle*>(fAODMCParticles->At(i));
	if (!part) {
	  AliError(Form("Could not find MC particle %d in branch %s of tree %s!", i, fAODMCParticlesName.Data(), fAODTreeName.Data()));
	  continue;
	}
	
	AliDebug(3, Form("Processing MC particle %d with pT = %f, eta = %f, phi = %f", i, part->Pt(), part->Eta(), part->Phi()));
	
	if (!part->IsPhysicalPrimary()) 
	  continue;

	if (part->Pt() < fPtMin || part->Pt() > fPtMax ||
	    part->Eta() < fEtaMin || part->Eta() > fEtaMax ||
	    part->Phi() < fPhiMin || part->Phi() > fPhiMax)
	  continue;
	
	AddMCParticle(part, i);
	AliDebug(3, "Embedded!");
      }
    }
  }

  if (fOutTracks) {

    if (fCopyArray && fTracks)
      CopyTracks();

    AliDebug(3, Form("Start embedding with %d tracks.", fOutTracks->GetEntriesFast()));

    if (fAODTracks) {
      AliDebug(3, Form("%d tracks will be processed for embedding.", fAODTracks->GetEntriesFast()));
      for (Int_t i = 0; i < fAODTracks->GetEntriesFast(); i++) {
	AliVTrack *track = static_cast<AliVTrack*>(fAODTracks->At(i));
	if (!track) {
	  AliError(Form("Could not find track %d in branch %s of tree %s!", i, fAODTrackName.Data(), fAODTreeName.Data()));
	  continue;
	}
	
	AliDebug(3, Form("Processing track %d with pT = %f, eta = %f, phi = %f, label = %d", i, track->Pt(), track->Eta(), track->Phi(), track->GetLabel()));
	
	Int_t type = 0;
	Bool_t isEmc = kFALSE;
	if (!fEsdTreeMode) {
	  AliAODTrack *aodtrack = static_cast<AliAODTrack*>(track);
	  if (aodtrack->TestFilterBit(fAODfilterBits[0])) {
	    type = 0;
	  }
	  else if (aodtrack->TestFilterBit(fAODfilterBits[1])) {
	    if ((aodtrack->GetStatus()&AliESDtrack::kITSrefit)==0) {
	      if (fIncludeNoITS)
		type = 2;
	      else {
		AliDebug(3, "Track not embedded because ITS refit failed.");
		continue;
	    }
	    }
	    else {
	      type = 1;
	    }
	  }
	  else { /*not a good track*/
	    AliDebug(3, "Track not embedded because not an hybrid track.");
	    continue;
	  }

	  if (TMath::Abs(aodtrack->GetTrackEtaOnEMCal()) < 0.75 && 
	      aodtrack->GetTrackPhiOnEMCal() > 70 * TMath::DegToRad() && 
	      aodtrack->GetTrackPhiOnEMCal() < 190 * TMath::DegToRad())
	    isEmc = kTRUE;
	}
	else if (fPicoTrackVersion > 0) { /*not AOD mode, let's see if it is PicoTrack*/
	  AliPicoTrack *ptrack = dynamic_cast<AliPicoTrack*>(track);
	  if (ptrack) {
	    if (fPicoTrackVersion >= 3)
	      type = ptrack->GetTrackType();
	    else
	      type = ptrack->GetLabel();
	    isEmc = ptrack->IsEMCAL();
	  }
	}
	
	if (track->Pt() < fPtMin || track->Pt() > fPtMax ||
	    track->Eta() < fEtaMin || track->Eta() > fEtaMax ||
	    track->Phi() < fPhiMin || track->Phi() > fPhiMax) {
	  AliDebug(3, "Track not embedded because out of limits.");
	  continue;
	}
	
	if (fTrackEfficiency < 1) {
	  Double_t r = gRandom->Rndm();
	  if (fTrackEfficiency < r) {
	    AliDebug(3, "Track not embedded because of artificial inefiiciency.");
	    continue;
	  }
	}
	
	Int_t label = 0;
	if (fIsAODMC) {
	  if (fUseNegativeLabels)
	    label = track->GetLabel();
	  else 
	    label = TMath::Abs(track->GetLabel());
	  
	  if (label == 0) 
	    AliWarning(Form("%s: Track %d with label==0", GetName(), i));
	}

	AddTrack(track->Pt(), track->Eta(), track->Phi(), type, track->GetTrackEtaOnEMCal(), track->GetTrackPhiOnEMCal(), isEmc, label);
	AliDebug(3, "Track embedded!");
      }
    }
  }

  if (fOutClusters) {

    if (fCopyArray && fClusters)
      CopyClusters();

    if (fAODClusters) {
      for (Int_t i = 0; i < fAODClusters->GetEntriesFast(); i++) {
	AliVCluster *clus = static_cast<AliVCluster*>(fAODClusters->At(i));
	if (!clus) {
	  AliError(Form("Could not find cluster %d in branch %s of tree %s!", i, fAODClusName.Data(), fAODTreeName.Data()));
	  continue;
	}
	TLorentzVector vect;
	Double_t vert[3] = {0,0,0};
	clus->GetMomentum(vect,vert);

	if (vect.Pt() < fPtMin || vect.Pt() > fPtMax ||
	    vect.Eta() < fEtaMin || vect.Eta() > fEtaMax ||
	    vect.Phi() < fPhiMin || vect.Phi() > fPhiMax)
	  continue;

	Int_t label = 0;
	if (fIsAODMC) {
	  label = clus->GetLabel();
	  if (label <= 0) 
	    AliWarning(Form("%s: Clus %d with label<=0", GetName(), i));
	}

	AddCluster(clus);
      }
    }
  }

  if (fOutCaloCells) {

    Double_t totalEnergy = 0;
    Int_t totalCells = 0;

    if (fCaloCells)
      totalCells += fCaloCells->GetNumberOfCells();
    if (fAODCaloCells)
      totalCells += fAODCaloCells->GetNumberOfCells();

    SetNumberOfOutCells(totalCells);

    if (fCaloCells)
      CopyCells();

    if (fAODCaloCells) {
      for (Short_t i = 0; i < fAODCaloCells->GetNumberOfCells(); i++) {
	Int_t mclabel = 0;
	Double_t efrac = 0.;
	Double_t time = -1;
	Short_t cellNum = -1;
	Double_t amp = -1;
	
	fAODCaloCells->GetCell(i, cellNum, amp, time, mclabel, efrac);

	if (fIsAODMC) {
	  if (mclabel <= 0) 
	    AliWarning(Form("%s: Cell %d with label<=0", GetName(), i));
	}
	else {
	  mclabel = 0;
	}

	AliDebug(2,Form("Adding cell with amplitude %f, absolute ID %d, time %f, mc label %d", amp, cellNum, time, mclabel));
	AddCell(amp, cellNum, time, mclabel);
	totalEnergy += amp;
      }
    }

    AliDebug(2,Form("Added cells = %d (energy = %f), total cells = %d", fAddedCells, totalEnergy, totalCells));
  }
}

//________________________________________________________________________
TLorentzVector AliJetEmbeddingFromAODTask::GetLeadingJet(TClonesArray *tracks, TClonesArray *clusters)
{
  TString name("kt");
  fastjet::JetAlgorithm jalgo(fastjet::kt_algorithm);
  if (fJetAlgo == 1) {
    name  = "antikt";
    jalgo = fastjet::antikt_algorithm;
  }

  // setup fj wrapper
  AliFJWrapper fjw(name, name);
  fjw.SetAreaType(fastjet::active_area_explicit_ghosts);
  fjw.SetGhostArea(1);  // set a very large ghost area to speed up jet finding
  fjw.SetR(fJetRadius);
  fjw.SetAlgorithm(jalgo);  
  fjw.SetMaxRap(1);
  fjw.Clear();

  if (tracks) {
    const Int_t Ntracks = tracks->GetEntries();
    for (Int_t iTracks = 0; iTracks < Ntracks; ++iTracks) {
      AliVParticle *t = static_cast<AliVParticle*>(tracks->At(iTracks));
      if (!t)
        continue;

      AliAODMCParticle *aodmcpart = dynamic_cast<AliAODMCParticle*>(t);
      if (aodmcpart && !aodmcpart->IsPhysicalPrimary())
	continue;
      
      if ((fJetType == 1 && t->Charge() == 0) ||
	  (fJetType == 2 && t->Charge() != 0))
	continue;

      if (t->Pt() < fJetConstituentMinPt)
	continue;

      fjw.AddInputVector(t->Px(), t->Py(), t->Pz(), t->P(), iTracks + 100);  
    }
  }

  if (clusters && fJetType != 1) {
    Double_t vert[3]={0};
    if (fAODVertex) 
      ((AliVVertex*)fAODVertex->At(0))->GetXYZ(vert);

    const Int_t Nclus = clusters->GetEntries();
    for (Int_t iClus = 0; iClus < Nclus; ++iClus) {
      AliVCluster *c = static_cast<AliVCluster*>(clusters->At(iClus));
      if (!c)
	continue;

      if (!c->IsEMCAL())
	continue;
      
      TLorentzVector nP;
      c->GetMomentum(nP, vert);

      if (nP.Pt() < fJetConstituentMinPt)
	continue;

      fjw.AddInputVector(nP.Px(), nP.Py(), nP.Pz(), nP.P(), -iClus - 100);
    }
  }
  
  // run jet finder
  fjw.Run();

  std::vector<fastjet::PseudoJet> jets_incl = fjw.GetInclusiveJets();
  AliDebug(1,Form("%d jets found", (Int_t)jets_incl.size()));

  TLorentzVector jet;

  Int_t njets = jets_incl.size();

  if (njets > 0) {
    //std::vector<fastjet::PseudoJet> jets_incl_sorted = fastjet::sorted_by_pt(jets_incl);
    for (Int_t i = 0; i < njets; i++) {
      if (jet.Pt() >= jets_incl[i].perp())
	continue;
      if (jets_incl[i].eta() < fJetMinEta || jets_incl[i].eta() > fJetMaxEta || jets_incl[i].phi() < fJetMinPhi || jets_incl[i].phi() > fJetMaxPhi)
	continue;
      jet.SetPxPyPzE(jets_incl[i].px(), jets_incl[i].py(), jets_incl[i].pz(), jets_incl[i].E());
    }
  }

  return jet;
}
