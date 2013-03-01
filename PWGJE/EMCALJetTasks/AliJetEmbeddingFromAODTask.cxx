// $Id: AliJetEmbeddingFromAODTask.cxx $
//
// Jet embedding from AOD task.
//
// Author: S.Aiola, C.Loizides

#include "AliJetEmbeddingFromAODTask.h"

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
  fIncludeNoITS(kTRUE),
  fUseNegativeLabels(kTRUE),
  fTrackEfficiency(1),
  fIsMC(kFALSE),
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
  fHistNotEmbedded(0)
{
  // Default constructor.
  SetSuffix("AODEmbedding");
  SetMarkMC(0);
  fAODfilterBits[0] = -1;
  fAODfilterBits[1] = -1;
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
  fIncludeNoITS(kTRUE),
  fUseNegativeLabels(kTRUE),
  fTrackEfficiency(1),
  fIsMC(kFALSE),
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
  fHistNotEmbedded(0)
{
  // Standard constructor.
  SetSuffix("AODEmbedding");
  SetMarkMC(0);
  fAODfilterBits[0] = -1;
  fAODfilterBits[1] = -1;
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
  TString fileName;

  do {
    if (i>0) {
      AliDebug(3,Form("Failed to open file %s...", fileName.Data()));
      if (fHistAODFileError)
	fHistAODFileError->Fill(fCurrentAODFileID);
    }

    fileName = GetNextFileName();

    if (fileName.IsNull())
      break;
    
    if (fileName.BeginsWith("alien://") && !gGrid) {
      AliInfo("Trying to connect to AliEn ...");
      TGrid::Connect("alien://");
    }

    AliDebug(3,Form("Trying to open file %s...", fileName.Data()));
    fCurrentAODFile = TFile::Open(fileName);

    i++;

  } while ((!fCurrentAODFile || fCurrentAODFile->IsZombie()) && i < fAttempts);

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

  AliDebug(2,Form("Will start embedding from entry %d", fCurrentAODEntry));
  
  if (fHistFileMatching)
    fHistFileMatching->Fill(fCurrentFileID, fCurrentAODFileID-1);
  
  return kTRUE;
}

//________________________________________________________________________
TString AliJetEmbeddingFromAODTask::GetNextFileName()
{
    if (fRandomAccess) 
      fCurrentAODFileID = TMath::Nint(gRandom->Rndm()*fFileList->GetEntriesFast());
    else
      fCurrentAODFileID++;

    if (fCurrentAODFileID >= fFileList->GetEntriesFast())
      return "";
    
    TObjString *objFileName = static_cast<TObjString*>(fFileList->At(fCurrentAODFileID));
    return objFileName->GetString();
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
  if (!fEsdTreeMode && fAODHeader) {
    AliAODHeader *aodHeader = static_cast<AliAODHeader*>(fAODHeader);

    if (fTriggerMask != AliVEvent::kAny) {
      UInt_t offlineTrigger = aodHeader->GetOfflineTrigger();
      
      if ((offlineTrigger & fTriggerMask) == 0) {
	AliDebug(2,Form("Event rejected due to physics selection. Event trigger mask: %d, trigger mask selection: %d.", offlineTrigger, fTriggerMask));
	return kFALSE;
      }
    }
    
    if (fMinCentrality >= 0) {
      AliCentrality *cent = aodHeader->GetCentralityP();
      Float_t centVal = cent->GetCentralityPercentile("V0M");
      if (centVal < fMinCentrality || centVal >= fMaxCentrality) {
	AliDebug(2,Form("Event rejected due to centrality selection. Event centrality: %f, centrality range selection: %f to %f", centVal, fMinCentrality, fMaxCentrality));
	return kFALSE;
      }
    }
  }

  if (fAODVertex) {
    Double_t vert[3]={0};
    ((AliVVertex*)fAODVertex->At(0))->GetXYZ(vert);
    if (TMath::Abs(vert[2]) > fZVertexCut) {
      AliDebug(2,Form("Event rejected due to Z vertex selection. Event Z vertex: %f, Z vertex cut: %f", vert[2], fZVertexCut));
      return kFALSE;
    }
  }

  return kTRUE;
}

//________________________________________________________________________
void AliJetEmbeddingFromAODTask::Run() 
{
  if (!GetNextEntry()) {
    if (fHistNotEmbedded)
      fHistNotEmbedded->Fill(fCurrentFileID);
    AliError("Unable to get the AOD event to embed. Nothing will be embedded.");
    return;
  }

  if (fOutMCParticles) {

    if (fCopyArray && fMCParticles)
      CopyMCParticles();

    if (fAODMCParticles) {
      for (Int_t i = 0; i < fAODMCParticles->GetEntriesFast(); i++) {
	AliAODMCParticle *part = static_cast<AliAODMCParticle*>(fAODMCParticles->At(i));
	if (!part) {
	  AliError(Form("Could not find MC particle %d in branch %s of tree %s!", i, fAODMCParticlesName.Data(), fAODTreeName.Data()));
	  continue;
	}
	
	if (!part->IsPhysicalPrimary()) 
	  continue;

	AliDebug(3, Form("Embedding MC particle with pT = %f, eta = %f, phi = %f", part->Pt(), part->Eta(), part->Phi()));
	AddMCParticle(part, i);
      }
    }
  }

  if (fOutTracks) {

    if (fCopyArray && fTracks)
      CopyTracks();

    if (fAODTracks) {
      AliDebug(2, Form("%d tracks will be processed for embedding.", fAODTracks->GetEntriesFast()));
      for (Int_t i = 0; i < fAODTracks->GetEntriesFast(); i++) {
	AliVTrack *track = static_cast<AliVTrack*>(fAODTracks->At(i));
	if (!track) {
	  AliError(Form("Could not find track %d in branch %s of tree %s!", i, fAODTrackName.Data(), fAODTreeName.Data()));
	  continue;
	}
	
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
	      else
		continue;
	    }
	    else {
	      type = 1;
	    }
	  }
	  else { /*not a good track*/
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
	
	if (fTrackEfficiency < 1) {
	  Double_t r = gRandom->Rndm();
	  if (fTrackEfficiency < r) 
	    continue;
	}
	
	Int_t label = 0;
	if (fIsMC) {
	  if (fUseNegativeLabels)
	    label = track->GetLabel();
	  else 
	    label = TMath::Abs(track->GetLabel());
	  
	  if (label == 0) {
	    AliWarning(Form("%s: Track %d with label==0", GetName(), i));
	    label = 99999;
	  }
	}

	AliDebug(3, Form("Embedding track with pT = %f, eta = %f, phi = %f, label = %d", track->Pt(), track->Eta(), track->Phi(), label));
	AddTrack(track->Pt(), track->Eta(), track->Phi(), type, track->GetTrackEtaOnEMCal(), track->GetTrackPhiOnEMCal(), isEmc, label);
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
	AddCluster(clus->E(), vect.Eta(), vect.Phi(), clus->GetLabel());
      }
    }
  }

  if (fOutCaloCells) {

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
	AliDebug(3,Form("Adding cell with amplitude %f, absolute ID %d, time %f", amp, cellNum, time));
	AddCell(amp, cellNum, time, mclabel);
      }
    }

    AliDebug(2,Form("Added cells = %d, total cells = %d", fAddedCells, totalCells));
  }
}
