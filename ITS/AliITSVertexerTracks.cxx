/**************************************************************************
 * Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
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

//-----------------------------------------------------------------
//    Implementation of the vertexer from tracks
//    It accepts V2 and ESD tracks
//
// Origin: A.Dainese, Padova, 
//         andrea.dainese@pd.infn.it
//         M.Masera,  Torino, 
//         massimo.masera@to.infn.it
//-----------------------------------------------------------------

//---- standard headers ----
#include <Riostream.h>
//---- Root headers --------
#include <TFile.h>
#include <TTree.h>
#include <TMatrixD.h>
//---- AliRoot headers -----
#include "AliESDVertex.h"
#include "AliITSVertexerTracks.h"
#include "AliESD.h"
#include "AliESDtrack.h"
#include "AliVertexerTracks.h"


ClassImp(AliITSVertexerTracks)


//----------------------------------------------------------------------------
AliITSVertexerTracks::AliITSVertexerTracks():AliITSVertexer(),
fInFile(0),
fOutFile(0),
fMinTracks(0),
fMaxChi2PerTrack(0),
fTrkArray(0),
fTrksToSkip(0),
fNTrksToSkip(0) {
//
// Default constructor
//
  SetVtxStart();
  SetMinTracks();

}
//----------------------------------------------------------------------------
AliITSVertexerTracks::AliITSVertexerTracks(TFile *inFile,TFile *outFile,
                                           Int_t fEv,Int_t lEv,
                                           Double_t xStart,Double_t yStart):
fInFile(inFile),
fOutFile(outFile),
fMinTracks(0),
fMaxChi2PerTrack(0),
fTrkArray(0),
fTrksToSkip(0),
fNTrksToSkip(0) {
//
// Standard constructor
//
  fCurrentVertex = 0;
  SetFirstEvent(fEv);
  SetLastEvent(lEv);
  SetVtxStart(xStart,yStart);
  SetMinTracks();
  SetDebug();
}
//----------------------------------------------------------------------------
AliITSVertexerTracks::AliITSVertexerTracks(TString fn,
					   Double_t xStart,Double_t yStart)
                                          :AliITSVertexer(fn),
fInFile(0),
fOutFile(0),
fMinTracks(0),
fMaxChi2PerTrack(0),
fTrkArray(0),
fTrksToSkip(0),
fNTrksToSkip(0) {
//
// Alternative constructor
//
  SetVtxStart(xStart,yStart);
  SetMinTracks();
}

//______________________________________________________________________
AliITSVertexerTracks::AliITSVertexerTracks(const AliITSVertexerTracks &vtxr) : AliITSVertexer(vtxr),
fInFile(vtxr.fInFile),
fOutFile(vtxr.fOutFile),
fMinTracks(vtxr.fMinTracks),
fMaxChi2PerTrack(vtxr.fMaxChi2PerTrack),
fTrkArray(vtxr.fTrkArray),
fTrksToSkip(vtxr.fTrksToSkip),
fNTrksToSkip(vtxr.fNTrksToSkip) {
  // Copy constructor

}

//______________________________________________________________________
AliITSVertexerTracks& AliITSVertexerTracks::operator=(const AliITSVertexerTracks&  vtxr ){
  //Assignment operator
  this->~AliITSVertexerTracks();
  new(this) AliITSVertexerTracks(vtxr);
  return *this;
}

//-----------------------------------------------------------------------------
AliITSVertexerTracks::~AliITSVertexerTracks() {
  // Default Destructor
  // The objects poited by the following pointers are not owned
  // by this class and are not deleted

  fCurrentVertex = 0;
  fInFile        = 0;
  fOutFile       = 0;
}
//-----------------------------------------------------------------------------
Bool_t AliITSVertexerTracks::CheckField() const {
//
// Check if the conv. const. has been set
//
  Double_t field = AliTracker::GetBz(); // in kG

  if(field<1 || field>6) {
    printf("AliITSVertexerTracks::CheckField():\n ERROR: AliKalmanTrack::fConvConst not set\n Use AliKalmanTrack::SetConvConst() or AliITSVertexerTracks::SetField()\n");
    return kFALSE;
  }
  printf("AliITSVertexerTracks::CheckField():  Using B = %3.1f T\n",field);
  return kTRUE;
}
//---------------------------------------------------------------------------
void AliITSVertexerTracks::ComputeMaxChi2PerTrack(Int_t nTracks) {
//
// Max. contr. to the chi2 has been tuned as a function of multiplicity
//
  if(nTracks < 7) { fMaxChi2PerTrack = 1.e6;
  } else { fMaxChi2PerTrack = 100.; }

  return;
}
//---------------------------------------------------------------------------
void AliITSVertexerTracks::FindVertices() {
//
// Vertices for all events from fFirstEvent to fLastEvent
//

  // Check if the conv. const. has been set
  if(!CheckField()) return;

  TDirectory *curdir = 0;

  // loop over events
  for(Int_t ev=fFirstEvent; ev<=fLastEvent; ev++) {
    if(ev % 100 == 0 || fDebug) printf("--- Processing event %d of %d ---\n",ev,fLastEvent);

    FindPrimaryVertexForCurrentEvent(ev);

    if(!fCurrentVertex) {
      printf("AliITSVertexerTracks::FindVertices(): no tracks tree for event %d\n",ev);
      continue;
    }

    if(fDebug) fCurrentVertex->PrintStatus();

    // write vertex to file
    TString vtxName = "Vertex_";
    vtxName += ev;
    fCurrentVertex->SetName(vtxName.Data()); 
    fCurrentVertex->SetTitle("VertexerTracks");
    //WriteCurrentVertex();
    curdir = gDirectory;
    fOutFile->cd();
    fCurrentVertex->Write();
    curdir->cd();
    fCurrentVertex = 0;
  } // loop over events

  return;
}
//---------------------------------------------------------------------------
void AliITSVertexerTracks::FindVerticesESD() {
//
// Vertices for all events from fFirstEvent to fLastEvent
//

 // Check if the conv. const. has been set
 if(!CheckField()) return;

 TDirectory *curdir = 0;

 fInFile->cd();
 TTree *esdTree = (TTree*)fInFile->Get("esdTree");
 if(!esdTree) {
     printf("AliITSVertexerTracks::FindVerticesESD(): no tree in file!\n");
   return;
 }
 Int_t nev = (Int_t)esdTree->GetEntries();
 Int_t ev;
 // loop on events in tree
 for(Int_t i=0; i<nev; i++) {
   AliESD *esdEvent = new AliESD;
   esdTree->SetBranchAddress("ESD",&esdEvent);
   if(!esdTree->GetEvent(i)) {
     printf("AliITSVertexerTracks::FindVerticesESD(): not an ESD!\n");
     delete esdEvent;
     return;
   }
   ev = (Int_t)esdEvent->GetEventNumber();
   if(ev<fFirstEvent || ev>fLastEvent) { delete esdEvent; continue; }
   if(ev % 100 == 0 || fDebug)
     printf("--- Processing event %d of %d ---\n",ev,fLastEvent-fFirstEvent);

   FindPrimaryVertexForCurrentEvent(esdEvent);

   if(!fCurrentVertex) {
     printf("AliITSVertexerTracks::FindVertixesESD():\n no vertex for event %d\n",ev);
     continue;
   }

   if(fDebug) fCurrentVertex->PrintStatus();

   // write vertex to file
   TString vtxName = "Vertex_";
   vtxName += ev;
   fCurrentVertex->SetName(vtxName.Data());
   fCurrentVertex->SetTitle("VertexerTracks");
   //WriteCurrentVertex();
   curdir = gDirectory;
   fOutFile->cd();
   fCurrentVertex->Write();
   curdir->cd();
   fCurrentVertex = 0;
   esdEvent = 0;
   delete esdEvent;
 } // end loop over events

 return;
}
//----------------------------------------------------------------------------
Int_t AliITSVertexerTracks::PrepareTracks(TTree &trkTree) {
  //
  // Propagate tracks to initial vertex position and store them in a TObjArray
  //
  Double_t maxd0rphi = 3.;  
  Double_t alpha,xlStart,d0rphi;
  Int_t    nTrks    = 0;
  Bool_t   skipThis;
  Int_t    nEntries = (Int_t)trkTree.GetEntries();

  Double_t field=AliTracker::GetBz();

  if(!fTrkArray.IsEmpty()) fTrkArray.Clear();
  fTrkArray.Expand(nEntries);

  if(fDebug) {
    printf(" PrepareTracks()\n");
    //    trkTree.Print();
  }

  for(Int_t i=0; i<nEntries; i++) {
    // check tracks to skip
    skipThis = kFALSE;
    for(Int_t j=0; j<fNTrksToSkip; j++) { 
      if(i==fTrksToSkip[j]) {
	if(fDebug) printf("skipping track: %d\n",i);
	skipThis = kTRUE;
      }
    }
    if(skipThis) continue;

    AliESDtrack *track = new AliESDtrack; 
    trkTree.SetBranchAddress("tracks",&track);
    trkTree.GetEvent(i);

    // propagate track to vtxSeed
    alpha  = track->GetAlpha();
    xlStart = fNominalPos[0]*TMath::Cos(alpha)+fNominalPos[1]*TMath::Sin(alpha);
    track->AliExternalTrackParam::PropagateTo(xlStart,field);   // to vtxSeed

    // select tracks with d0rphi < maxd0rphi
   
    d0rphi = TMath::Abs(track->GetD(fNominalPos[0],fNominalPos[1],field));
    if(d0rphi > maxd0rphi) { delete track; continue; }
   

    fTrkArray.AddLast(track);
    
    nTrks++; 
    if(fDebug)cout<<" :-) nTrks, d0rphi "<<nTrks<<"  "<<d0rphi<<endl;
   
  }
  if(fTrksToSkip) delete [] fTrksToSkip;

  return nTrks;
} 
//----------------------------------------------------------------------------
void AliITSVertexerTracks::PrintStatus() const {
//
// Print status
//  printf(" Initial position (%f,%f)\n",fNominalPos[0],fNominalPos[1]);
  printf(" Number of tracks in array: %d\n",(Int_t)fTrkArray.GetEntriesFast());
  printf(" Minimum # tracks required in fit: %d\n",fMinTracks);

  return;
}
//----------------------------------------------------------------------------
AliVertex* AliITSVertexerTracks::VertexForSelectedTracks(AliESD *esdEvent,Int_t nofCand, Int_t *trkPos, Int_t opt){

  //
  // Computes the vertex for selected tracks 
  // trkPos=vector with track positions in ESD
  //
  Double_t vtx[3];
  esdEvent->GetVertex()->GetXYZ(vtx);
  TTree *trkTree = new TTree("TreeT","tracks");
  AliESDtrack *esdTrack = 0;
  trkTree->Branch("tracks","AliESDtrack",&esdTrack);
  for(Int_t i=0; i<nofCand;i++){
    esdTrack = (AliESDtrack*)esdEvent->GetTrack(trkPos[i]);

    if(!esdTrack->GetStatus()&AliESDtrack::kTPCin) continue; 
    if(!esdTrack->GetStatus()&AliESDtrack::kITSin) continue; 
    if(!esdTrack->GetStatus()&AliESDtrack::kITSrefit) continue;

    Int_t nclus=esdTrack->GetNcls(0); // check number of clusters in ITS
    if(nclus<6) continue;
    trkTree->Fill();
  }
  delete esdTrack;
  Int_t nTrks = PrepareTracks(*trkTree);
  //delete trkTree;//  :-)) 
  if(fDebug) printf(" tracks prepared: %d\n",nTrks);
  if(nTrks < fMinTracks) {
    printf("TooFewTracks\n");
    AliVertex *theVert=new AliVertex();
    theVert->SetDispersion(999);
    theVert->SetNContributors(-5);
    return theVert;
  }
  
  AliVertexerTracks *vertexer=new AliVertexerTracks(vtx[0],vtx[1]);
  vertexer->SetFinderAlgorithm(opt);
  AliVertex *theVert=(AliVertex*)vertexer->VertexForSelectedTracks(&fTrkArray);
// beware: newvt object should be deleted by the caller
  AliVertex *newvt = new AliVertex(*theVert); 
  delete vertexer;
  return newvt;
}
//----------------------------------------------------------------------------
AliESDVertex* AliITSVertexerTracks::FindPrimaryVertexForCurrentEvent(Int_t evnumb) {
//
// Vertex for current event
//
  fCurrentVertex = 0;

  // get tree with tracks from input file
  fInFile->cd();
  TTree *esdTree = (TTree*)fInFile->Get("esdTree");

 if(!esdTree) {
     printf("AliITSVertexerTracks::FindPrimaryVertexForCurrentEvent(): no tree in file!\n");
   return fCurrentVertex;
 }
 AliESD *esdEvent = new AliESD;
 esdTree->SetBranchAddress("ESD",&esdEvent);
 esdTree->GetEvent(evnumb);
 return FindPrimaryVertexForCurrentEvent(esdEvent);
}
//----------------------------------------------------------------------------
AliESDVertex* AliITSVertexerTracks::FindPrimaryVertexForCurrentEvent(AliESD *esdEvent)
{
//
// Vertex for current ESD event
//
  fCurrentVertex = 0;
  Double_t vtx[3],cvtx[6];

  Int_t entr = (Int_t)esdEvent->GetNumberOfTracks();
  TTree *trkTree = new TTree("TreeT","tracks");
  AliESDtrack *esdTrack = 0;
  trkTree->Branch("tracks","AliESDtrack",&esdTrack);

  for(Int_t i=0; i<entr; i++) {
    AliESDtrack *et = esdEvent->GetTrack(i);
    esdTrack = new AliESDtrack(*et);
    if(!esdTrack->GetStatus()&AliESDtrack::kITSin) continue;
    if(!esdTrack->GetStatus()&AliESDtrack::kITSrefit) continue;
    Int_t nclus=esdTrack->GetNcls(0); // check number of clusters in ITS
    if(nclus<5) continue;

    trkTree->Fill();
  }
  delete esdTrack;

  // preselect tracks and propagate them to initial vertex position
  Int_t nTrks = PrepareTracks(*trkTree);
  delete trkTree;
  if(fDebug) printf(" tracks prepared: %d\n",nTrks);
  if(nTrks < fMinTracks) { TooFewTracks(); return fCurrentVertex; }

  // Set initial vertex position from ESD
  esdEvent->GetVertex()->GetXYZ(vtx);
  SetVtxStart(vtx[0],vtx[1]);

  // VERTEX FITTER
  ComputeMaxChi2PerTrack(nTrks);
  VertexFitter();
  if(fDebug) printf(" vertex fit completed\n");

  // store vertex information in ESD
  fCurrentVertex->GetXYZ(vtx);
  fCurrentVertex->GetCovMatrix(cvtx);

  Double_t tp[3];
  esdEvent->GetVertex()->GetTruePos(tp);
  fCurrentVertex->SetTruePos(tp);

  return fCurrentVertex;
}
//---------------------------------------------------------------------------
void AliITSVertexerTracks::SetSkipTracks(Int_t n,Int_t *skipped) {
//
// Mark the tracks not ot be used in the vertex finding
//
  fNTrksToSkip = n;
  fTrksToSkip = new Int_t[n]; 
  for(Int_t i=0;i<n;i++) fTrksToSkip[i] = skipped[i]; 
  return; 
}
//---------------------------------------------------------------------------
void AliITSVertexerTracks::TooFewTracks() {
//
// When the number of tracks is < fMinTracks the vertex is set to (0,0,0)
// and the number of tracks to -1
//
  fCurrentVertex = new AliESDVertex(0.,0.,-1);
  return;
}
//---------------------------------------------------------------------------
void AliITSVertexerTracks::VertexFitter() {
//
// The optimal estimate of the vertex position is given by a "weighted 
// average of tracks positions"
// Original method: CMS Note 97/0051
//
  if(fDebug) { 
    printf(" VertexFitter(): start\n");
    PrintStatus();
  }
  AliVertexerTracks *vertexer=new AliVertexerTracks(fNominalPos[0],fNominalPos[1]);
  vertexer->SetFinderAlgorithm(1);
  AliVertex *thevert=(AliVertex*)vertexer->VertexForSelectedTracks(&fTrkArray);
  Double_t initPos[3];
  thevert->GetXYZ(initPos);
  //  cout<<"Finder: "<<initPos[0]<<"; "<<initPos[1]<<"; "<<initPos[2]<<endl;
  delete vertexer;


  Int_t i,j,k,step=0;
  TMatrixD rv(3,1);
  TMatrixD vV(3,3);
  rv(0,0) = initPos[0];
  rv(1,0) = initPos[1];
  rv(2,0) = 0.;
  Double_t xlStart,alpha;
  Double_t rotAngle;
  Double_t cosRot,sinRot;
  Double_t cc[15];
  Int_t nUsedTrks;
  Double_t chi2,chi2i;
  Int_t arrEntries = (Int_t)fTrkArray.GetEntries();
  AliESDtrack *t = 0;
  Int_t failed = 0;

  Int_t *skipTrack = new Int_t[arrEntries];
  for(i=0; i<arrEntries; i++) skipTrack[i]=0;

  // 3 steps:
  // 1st - first estimate of vtx using all tracks
  // 2nd - apply cut on chi2 max per track
  // 3rd - estimate of global chi2
  for(step=0; step<3; step++) {
    if(fDebug) printf(" step = %d\n",step);
    chi2 = 0.;
    nUsedTrks = 0;

    TMatrixD sumWiri(3,1);
    TMatrixD sumWi(3,3);
    for(i=0; i<3; i++) {
      sumWiri(i,0) = 0.;
      for(j=0; j<3; j++) sumWi(j,i) = 0.;
    }

    // loop on tracks  
    for(k=0; k<arrEntries; k++) {
      if(skipTrack[k]) continue;
      // get track from track array
      t = (AliESDtrack*)fTrkArray.At(k);
      alpha = t->GetAlpha();
      xlStart = initPos[0]*TMath::Cos(alpha)+initPos[1]*TMath::Sin(alpha);
      t->AliExternalTrackParam::PropagateTo(xlStart,AliTracker::GetBz());   // to vtxSeed
      rotAngle = alpha;
      if(alpha<0.) rotAngle += 2.*TMath::Pi();
      cosRot = TMath::Cos(rotAngle);
      sinRot = TMath::Sin(rotAngle);
      
      // vector of track global coordinates
      TMatrixD ri(3,1);
      ri(0,0) = t->GetX()*cosRot-t->GetY()*sinRot;
      ri(1,0) = t->GetX()*sinRot+t->GetY()*cosRot;
      ri(2,0) = t->GetZ();

      // matrix to go from global (x,y,z) to local (y,z);
      TMatrixD qQi(2,3);
      qQi(0,0) = -sinRot;
      qQi(0,1) = cosRot;
      qQi(0,2) = 0.;
      qQi(1,0) = 0.;
      qQi(1,1) = 0.;
      qQi(1,2) = 1.;

      // covariance matrix of local (y,z) - inverted
      TMatrixD uUi(2,2);
      t->GetExternalCovariance(cc);
      uUi(0,0) = cc[0];
      uUi(0,1) = cc[1];
      uUi(1,0) = cc[1];
      uUi(1,1) = cc[2];

      // weights matrix: wWi = qQiT * uUiInv * qQi
      if(uUi.Determinant() <= 0.) continue;
      TMatrixD uUiInv(TMatrixD::kInverted,uUi);
      TMatrixD uUiInvQi(uUiInv,TMatrixD::kMult,qQi);
      TMatrixD wWi(qQi,TMatrixD::kTransposeMult,uUiInvQi);

      // track chi2
      TMatrixD deltar = rv; deltar -= ri;
      TMatrixD wWideltar(wWi,TMatrixD::kMult,deltar);
      chi2i = deltar(0,0)*wWideltar(0,0)+
              deltar(1,0)*wWideltar(1,0)+
	      deltar(2,0)*wWideltar(2,0);


      if(step==1 && chi2i > fMaxChi2PerTrack) {
	skipTrack[k] = 1;
	continue;
      }

      // add to total chi2
      chi2 += chi2i;

      TMatrixD wWiri(wWi,TMatrixD::kMult,ri); 

      sumWiri += wWiri;
      sumWi   += wWi;

      nUsedTrks++;
    } // end loop on tracks

    if(nUsedTrks < fMinTracks) {
      failed=1;
      continue;
    }

    Double_t determinant = sumWi.Determinant();
    //cerr<<" determinant: "<<determinant<<endl;
    if(determinant < 100.)  { 
      printf("det(V) = 0\n");       
      failed=1;
      continue;
    }

    // inverted of weights matrix
    TMatrixD invsumWi(TMatrixD::kInverted,sumWi);
    vV = invsumWi;
     
    // position of primary vertex
    rv.Mult(vV,sumWiri);

  } // end loop on the 3 steps

  delete [] skipTrack;
  delete t;

  if(failed) { 
    TooFewTracks(); 
    return; 
  }

  Double_t position[3];
  position[0] = rv(0,0);
  position[1] = rv(1,0);
  position[2] = rv(2,0);
  Double_t covmatrix[6];
  covmatrix[0] = vV(0,0);
  covmatrix[1] = vV(0,1);
  covmatrix[2] = vV(1,1);
  covmatrix[3] = vV(0,2);
  covmatrix[4] = vV(1,2);
  covmatrix[5] = vV(2,2);
  
  // store data in the vertex object
  fCurrentVertex = new AliESDVertex(position,covmatrix,chi2,nUsedTrks);

  if(fDebug) {
    printf(" VertexFitter(): finish\n");
    printf(" rv = ( %f , %f , %f )\n\n",rv(0,0),rv(1,0),rv(2,0));
    fCurrentVertex->PrintStatus();
  }

  return;
}
//----------------------------------------------------------------------------
AliESDVertex *AliITSVertexerTracks::VertexOnTheFly(TTree &trkTree) {
//
// Return vertex from tracks in trkTree
//
  if(fCurrentVertex) fCurrentVertex = 0;

  // get tracks and propagate them to initial vertex position
  Int_t nTrks = PrepareTracks(*(&trkTree));
  if(fDebug) printf(" tracks prepared: %d\n",nTrks);
  if(nTrks < fMinTracks) { TooFewTracks(); return fCurrentVertex; }

  // VERTEX FITTER
  ComputeMaxChi2PerTrack(nTrks);
  VertexFitter();
  if(fDebug) printf(" vertex fit completed\n");

  return fCurrentVertex;
}
//----------------------------------------------------------------------------



