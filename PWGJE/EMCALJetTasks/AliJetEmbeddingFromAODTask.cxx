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
  fTotalFiles(2200),
  fEsdTreeMode(kFALSE),
  fCurrentFileID(0),
  fCurrentAODFileID(-1),
  fCurrentAODFile(0),
  fPicoTrackVersion(0),
  fAODHeader(0),
  fAODVertex(0),
  fAODTracks(0),
  fAODClusters(0),
  fAODCaloCells(0),
  fAODMCParticles(0),
  fHistFileIDs(0)
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
  fTotalFiles(2200),
  fEsdTreeMode(kFALSE),
  fCurrentFileID(0),
  fCurrentAODFileID(-1),
  fCurrentAODFile(0),
  fPicoTrackVersion(0),
  fAODHeader(0),
  fAODVertex(0),
  fAODTracks(0),
  fAODClusters(0),
  fAODCaloCells(0),  
  fAODMCParticles(0),
  fHistFileIDs(0)
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
  
  fHistFileIDs = new TH2C("fHistFileIDs", "fHistFileIDs", fTotalFiles, -0.5, fTotalFiles - 0.5, fFileList->GetEntriesFast(), -0.5, fFileList->GetEntriesFast() -0.5);
  fHistFileIDs->GetXaxis()->SetTitle("File no. (PYTHIA)");
  fHistFileIDs->GetYaxis()->SetTitle("File no. (Embedded AOD)");
  fOutput->Add(fHistFileIDs);

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
  fCurrentAODFileID = fFileList->GetEntriesFast() * fCurrentFileID / fTotalFiles;
  AliInfo(Form("Start embedding from file ID %d", fCurrentAODFileID));
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

  if (fCurrentAODFileID >= fFileList->GetEntriesFast())
    return kFALSE;

  TObjString *objFileName = static_cast<TObjString*>(fFileList->At(fCurrentAODFileID));
  TString fileName = objFileName->GetString();
  
  if (fileName.BeginsWith("alien://") && !gGrid) {
      AliInfo("Trying to connect to AliEn ...");
      TGrid::Connect("alien://");
  }

  fCurrentAODFile = TFile::Open(fileName);

  if (!fCurrentAODFile || fCurrentAODFile->IsZombie()) {
    return kFALSE;
  }

  const TList *clist = fCurrentAODFile->GetStreamerInfoCache();
  if(clist) {
    TStreamerInfo *cinfo = static_cast<TStreamerInfo*>(clist->FindObject("AliPicoTrack"));
    if(cinfo) 
      fPicoTrackVersion = cinfo->GetClassVersion();
    else
      fPicoTrackVersion = 0;
  }

  if (fQAhistos)
    fHistFileIDs->Fill(fCurrentFileID, fCurrentAODFileID);
  
  fCurrentAODFileID++;
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliJetEmbeddingFromAODTask::GetNextEntry() 
{
  static TTree *tree = 0;
  static Int_t entry = 0;

  Int_t attempts = 0;

  do {
    
    if (!fCurrentAODFile || !tree || entry >= tree->GetEntries()) {
      if (!OpenNextFile())
	return kFALSE;
      
      tree = static_cast<TTree*>(fCurrentAODFile->Get(fAODTreeName));
      if (!tree)
	return kFALSE;

      if (!fAODHeaderName.IsNull()) 
	tree->SetBranchAddress(fAODHeaderName, &fAODHeader);

      if (!fAODVertexName.IsNull()) 
	tree->SetBranchAddress(fAODVertexName, &fAODVertex);
      
      if (!fAODTrackName.IsNull()) 
	tree->SetBranchAddress(fAODTrackName, &fAODTracks);
      
      if (!fAODClusName.IsNull()) 
	tree->SetBranchAddress(fAODClusName, &fAODClusters);
      
      if (!fAODCellsName.IsNull()) 
	tree->SetBranchAddress(fAODCellsName, &fAODCaloCells);

      if (!fAODMCParticlesName.IsNull()) 
	tree->SetBranchAddress(fAODMCParticlesName, &fAODMCParticles);
      
      entry = 0;
    }
    
    tree->GetEntry(entry);
    entry++;
    attempts++;

    if (attempts == 1000) 
      AliWarning("After 1000 attempts no event has been accepted by the event selection (trigger, centrality...)!");

  } while (!IsAODEventSelected() && tree);

  return (tree!=0);
}

//________________________________________________________________________
Bool_t AliJetEmbeddingFromAODTask::IsAODEventSelected()
{
  if (!fEsdTreeMode && fAODHeader) {
    AliAODHeader *aodHeader = static_cast<AliAODHeader*>(fAODHeader);

    UInt_t offlineTrigger = aodHeader->GetOfflineTrigger();
  
    if ((offlineTrigger & fTriggerMask) == 0) {
      AliDebug(2, Form("Event rejected due to physics selection. Event trigger mask: %d, trigger mask selection: %d.", offlineTrigger, fTriggerMask));
      return kFALSE;
    }
    
    AliCentrality *cent = aodHeader->GetCentralityP();
    Float_t centVal = cent->GetCentralityPercentile("V0M");
    if (centVal < fMinCentrality || centVal >= fMaxCentrality) {
      AliDebug(2, Form("Event rejected due to centrality selection. Event centrality: %f, centrality range selection: %f to %f", centVal, fMinCentrality, fMaxCentrality));
      return kFALSE;
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
    AliError("Unable to get AOD event to embed. Nothing will be embedded.");
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
	
	AliDebug(3, Form("Embedding track with pT = %f, eta = %f, phi = %f", track->Pt(), track->Eta(), track->Phi()));
	AddTrack(track->Pt(), track->Eta(), track->Phi(), type, track->GetTrackEtaOnEMCal(), track->GetTrackPhiOnEMCal(), isEmc, track->GetLabel());
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
	Short_t mclabel = 0;
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
