//-------------------------------------------------------------------------
//     Task for the Analysis Framework 
// Updates the an already created AOD with the PWG2 information taken from 
// the ESD.
//  - Puts the per-track information into the AliPWG2AODTrack container, 
//    together with the link to the original AliAODTrack
//
//     Author: Adam Kisiel, OSU, Adam.Kisiel@cern.ch
//-------------------------------------------------------------------------
#include <TChain.h>
#include <TFile.h>
#include <TList.h> 

#include "AliAnalysisTaskPWG2AODUpdate.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliESDtrack.h"
#include "AliESDVertex.h"
#include "AliLog.h"
#include "AliPWG2AODTrack.h"

ClassImp(AliAnalysisTaskPWG2AODUpdate)

////////////////////////////////////////////////////////////////////////

AliAnalysisTaskPWG2AODUpdate::AliAnalysisTaskPWG2AODUpdate():
  fESD(0x0),
  fAOD(0x0),
  fPWG2AODTracks(0x0)
{
  // Default constructor
}

AliAnalysisTaskPWG2AODUpdate::AliAnalysisTaskPWG2AODUpdate(const char* name):
  AliAnalysisTaskSE(name),
  fESD(0x0),
  fAOD(0x0),
  fPWG2AODTracks(0x0)
{
  // Default constructor
  DefineInput (0, TChain::Class());
  DefineOutput(0, TTree::Class());
}

AliAnalysisTaskPWG2AODUpdate::AliAnalysisTaskPWG2AODUpdate(const AliAnalysisTaskPWG2AODUpdate &task):
  AliAnalysisTaskSE(),
  fESD(0x0),
  fAOD(0x0),
  fPWG2AODTracks(0x0)
{
  // Copy
  fESD = task.fESD;
  fAOD =  task.fAOD;
  fPWG2AODTracks = task.fPWG2AODTracks;
}

AliAnalysisTaskPWG2AODUpdate& AliAnalysisTaskPWG2AODUpdate::operator=(const AliAnalysisTaskPWG2AODUpdate &task)
{
  // Assignment
  if (&task == this) return *this;
  TTask::operator=(task);

  fESD = task.fESD;
  fAOD =  task.fAOD;
  fPWG2AODTracks = task.fPWG2AODTracks;

  return *this;
}

void AliAnalysisTaskPWG2AODUpdate::UserCreateOutputObjects()
{
  // Get links to input and output containters  
  fAOD   = AODEvent();
  
  fPWG2AODTracks = new TClonesArray("AliPWG2AODTrack", 0);
  const char *name = "pwg2aodtracks";
  fPWG2AODTracks->SetName(name);
  
  AddAODBranch("TClonesArray", &fPWG2AODTracks);
}

void AliAnalysisTaskPWG2AODUpdate::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event
  //
  fESD = dynamic_cast<AliESDEvent *> (InputEvent());
  if (!fESD) { printf ("No input ESD !!! Not adding PWG2 information\n"); }
  fAOD = AODEvent();
  if (!fAOD) { printf ("No output AOD !!! Not adding PWG2 information\n"); }
 
  Double_t pos[3];
  Double_t p[3];
  Double_t tpcentrancepoint[3];
  Double_t tpcexitpoint[3];
  TBits    sharemap;
  TBits    clustermap;
  Int_t    pTracks = 0;
  
  for (Int_t i = 0; i < 3; i++)  { tpcentrancepoint[i] = 0.; tpcexitpoint[i] = 0; }

  // Get primary vertex
  const AliESDVertex *vtx = fESD->GetPrimaryVertex();
  
  vtx->GetXYZ(pos); // position
    
  // Tracks (primary and orphan)
  Int_t nTracks = fAOD->GetNTracks();

  printf("NUMBER OF AOD TRACKS %5d\n", nTracks);
  
  // Must start from scratch for each event
  fPWG2AODTracks->Delete();
  fPWG2AODTracks->Expand(nTracks);
  TClonesArray &tPWG2AODTracks = *fPWG2AODTracks;
    
  // Loop over AOD tracks
  for (Int_t nTrack = 0; nTrack < nTracks; ++nTrack) {
    AliAODTrack *aodTrack = fAOD->GetTrack(nTrack);
    Short_t trackId = aodTrack->GetID();

    // Get the corresponding ESD track
    AliESDtrack *esdTrack = fESD->GetTrack(trackId);

    esdTrack->GetPxPyPz(p);

    if (TMath::Abs(p[0] - aodTrack->Px())>0.000000001) {
      printf("Got different momenta !!! %f %f\n", p[0], aodTrack->Px());
    }

    sharemap = esdTrack->GetTPCSharedMap();
    clustermap = esdTrack->GetTPCClusterMap();
	
    esdTrack->GetInnerXYZ(tpcentrancepoint);
    tpcentrancepoint[2] -= pos[2];
	
    esdTrack->GetOuterXYZ(tpcexitpoint);
    tpcexitpoint[2] -= pos[2];

    // Add the PWG2 info into the AOD
    new (tPWG2AODTracks[pTracks++]) AliPWG2AODTrack(tpcentrancepoint,
						    tpcexitpoint,
						    sharemap,
						    clustermap,
						    aodTrack);
  }
  
  printf("PWG2=%d\n", fPWG2AODTracks->GetEntries());

  return;
}

