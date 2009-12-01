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
#include <TNtuple.h>  
#include <TCanvas.h>

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDtrack.h"
#include "AliExternalTrackParam.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include "AliVertexerTracks.h"
#include "AliESDInputHandler.h"
#include "AliAnalysisTaskIPInfo.h"
#include "AliIntSpotEstimator.h"
#include "AliMultiplicity.h"


ClassImp(AliAnalysisTaskIPInfo)

const Char_t* AliAnalysisTaskIPInfo::fEstNames[kNEst] = {"ITSTPC","TPC","SPD"};

//________________________________________________________________________
AliAnalysisTaskIPInfo::AliAnalysisTaskIPInfo(const char *name) 
: AliAnalysisTask(name, "IP analysis"),
  fESD(0),fESDfriend(0),fOutput(0),fTracks(50)
{

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TList container
  DefineOutput(0, TList::Class());  //My private output
  //
  for (int i=0;i<kNEst;i++) fIPEst[i] = 0;
  //
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
void AliAnalysisTaskIPInfo::SetIPCenIni(Int_t estID, Double_t x,Double_t y,Double_t z)
{
  // set initial estimate of the IP center
  if (estID<0 || estID>= kNEst) return;
  fIPCenIni[estID][0] = x;
  fIPCenIni[estID][1] = y;
  fIPCenIni[estID][2] = z;
}

//________________________________________________________________________
void AliAnalysisTaskIPInfo::SetOptions(Int_t estID, Bool_t recoVtx,
				       Double_t outcut,Int_t ntrIP,Int_t nPhiBins,Int_t nestb,
				       Double_t estmin,Double_t estmax,
				       Int_t ntrBins,Int_t ntMn,Int_t ntMx,
				       Int_t nPBins,Double_t pmn,Double_t pmx,Bool_t fillNt)
{
  // set options for estimators
  if (estID<0 || estID>= kNEst) return;
  fNTrMinIP[estID] = ntrIP;
  fRecoVtx[estID]  = recoVtx;
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
  fFillNt[estID]   = fillNt;
}

//________________________________________________________________________
void AliAnalysisTaskIPInfo::ConnectInputData(Option_t *) 
{
  // Connect ESD or AOD here
  // Called once
  //
  AliInfo("HERE");
  TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
  if (!tree) {
    Printf("ERROR: Could not read chain from input slot 0");
  } 
  else {
    tree->SetBranchAddress("ESDfriend.",&fESDfriend);
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
  Double_t err[3]={0.0400,0.0400,7.5};
  //
  for (int i=0;i<kNEst;i++) if (fNEstb[i]>1) {
      TString nm   = GetName();
      TString nmes = fEstNames[i];
      nm += nmes;
      fIPEst[i] = new AliIntSpotEstimator(nm.Data(),fOutCut[i],fNTrMinIP[i],fNPhiBins[i],
					  fNEstb[i],fEstMin[i],fEstMax[i],
					 fNTrBins[i],fNTrMin[i],fNTrMax[i],
					 fNPBins[i],fPMin[i],fPMax[i],fFillNt[i]);
      AliESDVertex *initVertex = new AliESDVertex(fIPCenIni[i],err);
      fIPEst[i]->GetVertexer()->SetVtxStart(initVertex);
      delete initVertex;
      fIPEst[i]->GetVertexer()->SetConstraintOff();
      fIPEst[i]->GetVertexer()->SetMinClusters(2);
      fIPEst[i]->SetIPCenIni(fIPCenIni[i]);
      if (nmes == "TPC") fIPEst[i]->GetVertexer()->SetTPCMode();
      else               fIPEst[i]->GetVertexer()->SetITSMode();
      //
      fOutput->Add(fIPEst[i]);
      if (fIPEst[i]->GetNtuple()) fOutput->Add(fIPEst[i]->GetNtuple());
    }
  //
  return;
}

//________________________________________________________________________
void AliAnalysisTaskIPInfo::Exec(Option_t *) 
{
  static TClonesArray tracks("AliExternalTrackParam",50);
  //
  // Main loop
  // Called for each event
  if (!fESD) {
    Printf("ERROR: fESD not available");
    return;
  }
  fESD->SetESDfriend(fESDfriend);
  //
  const AliESDVertex *vtx;
  UShort_t *trackID;
  Int_t ntracks;
  //
  for (int ie=0;ie<kNEst;ie++) {
    //
    if (!fIPEst[ie]) continue;
    if (ie==kTPC) {
      fIPEst[kTPC]->GetVertexer()->SetFieldkG( fESD->GetMagneticField() );
      vtx = fRecoVtx[kTPC] ? fIPEst[kTPC]->GetVertexer()->FindPrimaryVertex(fESD) : fESD->GetPrimaryVertexTPC();
      if (vtx) {
	ntracks = vtx->GetNIndices();
	trackID = (UShort_t*)vtx->GetIndices();
	for (int i=ntracks;i--;) fTracks.Add((TObject*)fESD->GetTrack(trackID[i])->GetTPCInnerParam());
	fIPEst[kTPC]->ProcessEvent(&fTracks);
	fTracks.Clear();
      }
    }
    else if (ie==kITSTPC) {
      fIPEst[kITSTPC]->GetVertexer()->SetFieldkG( fESD->GetMagneticField() );
      vtx = fRecoVtx[kITSTPC] ? fIPEst[kITSTPC]->GetVertexer()->FindPrimaryVertex(fESD) : fESD->GetPrimaryVertex();
      if (vtx) {
	ntracks = vtx->GetNIndices();
	trackID = (UShort_t*)vtx->GetIndices();
	for (int i=ntracks;i--;) fTracks.Add((TObject*)fESD->GetTrack(trackID[i]));
	fIPEst[kITSTPC]->ProcessEvent(&fTracks);
	fTracks.Clear();
      }
    }
    else if (ie==kSPD) {
      fIPEst[kSPD]->GetVertexer()->SetFieldkG( fESD->GetMagneticField() );
      ntracks = CreateSPDTracklets(tracks);
      for (int i=ntracks;i--;) fTracks.Add((TObject*)tracks[i]);
      fIPEst[kSPD]->ProcessEvent(&fTracks);
      fTracks.Clear();
    }
  }
  //
  PostData(0, fOutput);
  //
  return;
}      

//________________________________________________________________________
Int_t AliAnalysisTaskIPInfo::CreateSPDTracklets(TClonesArray& tracks)
{
  // create traclets from multiplicity class
  double cv[21] = {
    25e-4,
    0   ,  25e-4,
    0   ,     0,  40e-2,
    0   ,     0,      0,   1e-2,
    0   ,     0,      0,      0,   1e-2,
    0   ,     0,      0,      0,      0,   1e-2
  };
  //
  double xyzt[3],pxyz[3];
  int nSPDtracks = 0;
  tracks.Delete();
  const AliMultiplicity *mult = fESD->GetMultiplicity();
  const AliESDVertex *spdVtx  = fESD->GetPrimaryVertexSPD();
  int nTracklets = 0;
  if (mult && spdVtx && (nTracklets=mult->GetNumberOfTracklets())>2 ) {
    const Double_t kRLay1=3.9, kRLay2=7.6;
    double xv = spdVtx->GetXv();
    double yv = spdVtx->GetYv();
    double zv = spdVtx->GetZv();
    for (int i=0;i<nTracklets;i++) { // get cluster coordinates from tracklet
      double phi1 = mult->GetPhi(i);
      double tht1 = mult->GetTheta(i);
      double phi2 = phi1 - mult->GetDeltaPhi(i);
      double tht2 = tht1 - mult->GetDeltaTheta(i);
      double cs = TMath::Cos(phi1);
      double sn = TMath::Sin(tht1);
      double det = xv*sn+yv*sn;
      det = det*det + kRLay1*kRLay1;
      double t  = -(xv*cs+yv*sn) + TMath::Sqrt(det);
      double x1 = cs*t;
      double y1 = sn*t;
      double z1 = zv + TMath::Sqrt(x1*x1+y1*y1)/TMath::Tan(tht1);
      x1 += xv;
      y1 += yv;
      //
      cs = TMath::Cos(phi2);
      sn = TMath::Sin(tht2);
      det = xv*sn+yv*sn;
      det = det*det + kRLay2*kRLay2;
      t  = -(xv*cs+yv*sn) + TMath::Sqrt(det);
      double dx = cs*t;
      double dy = sn*t;
      double dz = zv + TMath::Sqrt(dx*dx+dy*dy)/TMath::Tan(tht2);
      dx += xv-x1;
      dy += yv-y1;
      dz += -z1;
      double dr = TMath::Sqrt(dx*dx+dy*dy+dz*dz);
      pxyz[0] = dx/dr; // direction cosines
      pxyz[1] = dy/dr;
      pxyz[2] = dz/dr;
      t = (xv-x1)*pxyz[0] + (yv-y1)*pxyz[1] + (zv-z1)*pxyz[2];
      xyzt[0] = x1 + t*pxyz[0];  // PCA to vertex
      xyzt[1] = y1 + t*pxyz[1];
      xyzt[2] = z1 + t*pxyz[2];
      //
      new(tracks[nSPDtracks++]) AliExternalTrackParam(xyzt,pxyz,cv,0);
    }
  }
  return nSPDtracks;
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
