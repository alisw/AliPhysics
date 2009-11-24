/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

//*************************************************************************
// Class AliAnalysisTaskIPInfo
// AliAnalysisTask to extract from the ESD the IP position and sigma
// as well as to estimate the primary vertex and tracks DCA resolution.
// Uses external class AliIntSpotEstimator
//
// Author: ruben.shahoyan@cern.ch
//*************************************************************************

#include <TChain.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TBranch.h>
#include <TClonesArray.h>
#include <TObjArray.h>
#include <TH1F.h>
#include <TH2F.h>  
#include <TCanvas.h>

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDtrack.h"
#include "AliExternalTrackParam.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliVertexerTracks.h"
#include "AliESDInputHandler.h"
#include "AliAnalysisTaskIPInfo.h"
#include "AliIntSpotEstimator.h"


ClassImp(AliAnalysisTaskIPInfo)

//________________________________________________________________________
AliAnalysisTaskIPInfo::AliAnalysisTaskIPInfo(const char *name) 
: AliAnalysisTask(name, "IP analysis"),
  fESD(0),fOutput(0),fTracks(50)
{

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TList container
  DefineOutput(0, TList::Class());  //My private output
  //
  for (int i=0;i<kNEst;i++) { // default settings
    fIPEst[i] = 0;
    SetOptions(i);
  }
}

//________________________________________________________________________
AliAnalysisTaskIPInfo::~AliAnalysisTaskIPInfo()
{
  // Destructor
  if (fOutput) {
    delete fOutput;
    fOutput = 0;
  }
}

//________________________________________________________________________
void AliAnalysisTaskIPInfo::SetOptions(Int_t estID, Bool_t recoVtx,
				       Double_t outcut,Int_t nPhiBins,Int_t nestb,
				       Double_t estmin,Double_t estmax,
				       Int_t ntrBins,Int_t ntMn,Int_t ntMx,
				       Int_t nPBins,Double_t pmn,Double_t pmx)
{
  // set options for estimators
  if (estID<0 || estID>= kNEst) return;
  fRecoVtx[estID] = recoVtx;
  fNPhiBins[estID] = nPhiBins;
  fNEstb[estID]    = nestb;
  fNTrBins[estID]  = ntrBins;
  fNPBins[estID]   = nPBins;
  fNTrMin[estID]   = ntMn;
  fNTrMax[estID]   = ntMx;
  fOutCut[estID]   = outcut;
  fEstMin[estID]   = estmin;
  fEstMax[estID]   = estmax;
  fPMin[estID]     = pmn;
  fPMax[estID]     = pmx;
}

//________________________________________________________________________
void AliAnalysisTaskIPInfo::ConnectInputData(Option_t *) 
{
  // Connect ESD or AOD here
  // Called once
  //
  TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
  if (!tree) {
    Printf("ERROR: Could not read chain from input slot 0");
  } 
  else {
    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());    
    if (!esdH) Printf("ERROR: Could not get ESDInputHandler");
    else fESD = esdH->GetEvent();
  }
  //
  return;
}

//________________________________________________________________________
void AliAnalysisTaskIPInfo::CreateOutputObjects()
{
  // Create estimators
  // Several histograms are more conveniently managed in a TList
  fOutput = new TList;
  fOutput->SetOwner();
  //
  TString nm = GetName();
  nm += "_TPC";
  fIPEst[kTPC] = new AliIntSpotEstimator(nm.Data(),fOutCut[kTPC],fNPhiBins[kTPC],
					 fNEstb[kTPC],fEstMin[kTPC],fEstMax[kTPC],
					 fNTrBins[kTPC],fNTrMin[kTPC],fNTrMax[kTPC],
					 fNPBins[kTPC],fPMin[kTPC],fPMax[kTPC]);
  fIPEst[kTPC]->GetVertexer()->SetTPCMode();
  fIPEst[kTPC]->GetVertexer()->SetConstraintOff();
  //
  nm = GetName();
  nm += "_ITSTPC";
  fIPEst[kITS] = new AliIntSpotEstimator(nm.Data(),fOutCut[kITS],fNPhiBins[kITS],
					 fNEstb[kITS],fEstMin[kITS],fEstMax[kITS],
					 fNTrBins[kITS],fNTrMin[kITS],fNTrMax[kITS],
					 fNPBins[kITS],fPMin[kITS],fPMax[kITS]);
  //
  fIPEst[kITS]->GetVertexer()->SetITSMode();
  fIPEst[kITS]->GetVertexer()->SetConstraintOff();
  //
  fOutput->Add(fIPEst[kTPC]);
  fOutput->Add(fIPEst[kITS]);
  //
  return;
}

//________________________________________________________________________
void AliAnalysisTaskIPInfo::Exec(Option_t *) 
{
  // Main loop
  // Called for each event
  
  if (!fESD) {
    Printf("ERROR: fESD not available");
    return;
  }
  //
  fIPEst[kTPC]->GetVertexer()->SetFieldkG( fESD->GetMagneticField() );
  fIPEst[kITS]->GetVertexer()->SetFieldkG( fESD->GetMagneticField() );
  const AliESDVertex *vtx;
  UShort_t *trackID;
  Int_t ntracks;
  //
  // Process vertex made of TPC tracks only
  vtx = fRecoVtx[kTPC] ? fIPEst[kTPC]->GetVertexer()->FindPrimaryVertex(fESD) : fESD->GetPrimaryVertexTPC();
  if (vtx) {
    ntracks = vtx->GetNIndices();
    trackID = (UShort_t*)vtx->GetIndices();
    for (int i=ntracks;i--;) fTracks.Add((TObject*)fESD->GetTrack(trackID[i])->GetTPCInnerParam());
    fIPEst[kTPC]->ProcessEvent(&fTracks);
    fTracks.Clear();
  }
  //
  // Process vertex made of TPC+ITS tracks
  vtx = fRecoVtx[kITS] ? fIPEst[kITS]->GetVertexer()->FindPrimaryVertex(fESD) : fESD->GetPrimaryVertex();
  if (vtx) {
    ntracks = vtx->GetNIndices();
    trackID = (UShort_t*)vtx->GetIndices();
    for (int i=ntracks;i--;) fTracks.Add((TObject*)fESD->GetTrack(trackID[i]));
    fIPEst[kITS]->ProcessEvent(&fTracks);
    fTracks.Clear();
  }
  //  
  PostData(0, fOutput);
  //
  return;
}      

//________________________________________________________________________
void AliAnalysisTaskIPInfo::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query
  fOutput = dynamic_cast<TList*> (GetOutputData(0));
  if (!fOutput) {     
    Printf("ERROR: fOutput not available");
    return;
  }
  //
  return;
}
