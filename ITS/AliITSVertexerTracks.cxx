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
#include <TKey.h>
#include <TFile.h>
#include <TTree.h>
#include <TMatrixD.h>
//---- AliRoot headers -----
#include "AliITSStrLine.h"
#include "AliITStrackV2.h"
#include "AliESDVertex.h"
#include "AliITSVertexerTracks.h"
#include "AliESD.h"
#include "AliESDtrack.h"


ClassImp(AliITSVertexerTracks)


//----------------------------------------------------------------------------
AliITSVertexerTracks::AliITSVertexerTracks():AliITSVertexer() {
//
// Default constructor
//
  fInFile  = 0;
  fOutFile = 0;
  SetVtxStart();
  SetMinTracks();
  fTrksToSkip = 0;
  fNTrksToSkip = 0;
  for(Int_t i=0; i<3; i++)fInitPos[i] = 0.;
}
//----------------------------------------------------------------------------
AliITSVertexerTracks::AliITSVertexerTracks(TFile *inFile,TFile *outFile,
                                           Double_t field,
                                           Int_t fEv,Int_t lEv,
                                           Double_t xStart,Double_t yStart) {
//
// Standard constructor
//
  fCurrentVertex = 0;
  fInFile  = inFile;
  fOutFile = outFile;
  SetFirstEvent(fEv);
  SetLastEvent(lEv);
  SetField(field);
  SetVtxStart(xStart,yStart);
  SetMinTracks();
  fTrksToSkip = 0;
  fNTrksToSkip = 0;
  for(Int_t i=0; i<3; i++) fInitPos[i] = 0.;
  SetDebug();
}
//----------------------------------------------------------------------------
AliITSVertexerTracks::AliITSVertexerTracks(Double_t field, TString fn,
					   Double_t xStart,Double_t yStart)
                                          :AliITSVertexer(fn) {
//
// Alternative constructor
//
  fInFile  = 0;
  fOutFile = 0;
  SetField(field);
  SetVtxStart(xStart,yStart);
  SetMinTracks();
  fTrksToSkip = 0;
  fNTrksToSkip = 0;
  for(Int_t i=0; i<3; i++) fInitPos[i] = 0.;
}
//______________________________________________________________________
AliITSVertexerTracks::AliITSVertexerTracks(const AliITSVertexerTracks &vtxr) : AliITSVertexer(vtxr) {
  // Copy constructor
  // Copies are not allowed. The method is protected to avoid misuse.
  Error("AliITSVertexerTracks","Copy constructor not allowed\n");
}

//______________________________________________________________________
AliITSVertexerTracks& AliITSVertexerTracks::operator=(const AliITSVertexerTracks& /* vtxr */){
  // Assignment operator
  // Assignment is not allowed. The method is protected to avoid misuse.
  Error("= operator","Assignment operator not allowed\n");
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
  AliITStrackV2 t;
  Double_t cc    = t.GetConvConst();
  Double_t field = 100./0.299792458/cc;

  if(field<0.1 || field>0.6) {
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

    FindVertexForCurrentEvent(ev);

    if(!fCurrentVertex) {
      printf("AliITSVertexerTracks::FindVertixes(): no tracks tree for event %d\n",ev);
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
  TKey *key=0;
  TIter next(fInFile->GetListOfKeys());
  // loop on events in file
  while ((key=(TKey*)next())!=0) {
    AliESD *esdEvent=(AliESD*)key->ReadObj();
    if(!esdEvent) { 
      printf("AliITSVertexerTracks::FindVerticesESD(): not an ESD!\n"); 
      return; 
    }
    Int_t ev = (Int_t)esdEvent->GetEventNumber();
    if(ev<fFirstEvent || ev>fLastEvent) { delete esdEvent; continue; }
    if(ev % 100 == 0 || fDebug) 
      printf("--- Processing event %d of %d ---\n",ev,fLastEvent-fFirstEvent);

    FindVertexForCurrentEvent(esdEvent);

    if(!fCurrentVertex) {
      printf("AliITSVertexerTracks::FindVertixesESD():\n no vertex for event %d\n",ev);
      continue;
    }

    cout<<"VERTICE TROVATO\n";
    if(fDebug) fCurrentVertex->PrintStatus();

    // write the ESD to file
    curdir = gDirectory;
    Char_t ename[100];
    sprintf(ename,"%d",ev);
    fOutFile->cd();
    esdEvent->Dump();
    esdEvent->Write(ename,TObject::kOverwrite);
    curdir->cd();
    fCurrentVertex = 0;
  } // loop over events

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

  if(!fTrkArray.IsEmpty()) fTrkArray.Clear();
  fTrkArray.Expand(nEntries);

  if(fDebug) {
    printf(" PrepareTracks()\n");
    trkTree.Print();
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

    AliITStrackV2 *itstrack = new AliITStrackV2; 
    trkTree.SetBranchAddress("tracks",&itstrack);
    trkTree.GetEvent(i);


    // propagate track to vtxSeed
    alpha  = itstrack->GetAlpha();
    xlStart = fNominalPos[0]*TMath::Cos(alpha)+fNominalPos[1]*TMath::Sin(alpha);
    itstrack->PropagateTo(3.,0.0023,65.19); // to beam pipe (0.8 mm of Be) 
    itstrack->PropagateTo(xlStart,0.,0.);   // to vtxSeed

    // select tracks with d0rphi < maxd0rphi
    d0rphi = TMath::Abs(itstrack->GetD(fNominalPos[0],fNominalPos[1]));
    if(d0rphi > maxd0rphi) { delete itstrack; continue; }
   
    fTrkArray.AddLast(itstrack);

    nTrks++; 
  }

  if(fTrksToSkip) delete [] fTrksToSkip;

  return nTrks;
} 
//----------------------------------------------------------------------------
void AliITSVertexerTracks::PrintStatus() const {
//
// Print status
//
  printf(" Initial position (%f,%f)\n",fNominalPos[0],fNominalPos[1]);
  printf(" Vertex position after vertex finder (%f, %f, %f)\n",fInitPos[0],fInitPos[1],fInitPos[2]);
  printf(" Number of tracks in array: %d\n",(Int_t)fTrkArray.GetEntriesFast());
  printf(" Minimum # tracks required in fit: %d\n",fMinTracks);

  return;
}
//----------------------------------------------------------------------------
AliESDVertex* AliITSVertexerTracks::FindVertexForCurrentEvent(Int_t evnumb) {
//
// Vertex for current event
//
  fCurrentVertex = 0;

  // get tree with tracks from input file
  TString treeName = "TreeT_ITS_";
  treeName += evnumb;
  TTree *trkTree=(TTree*)fInFile->Get(treeName.Data());
  if(!trkTree) return fCurrentVertex;


  // get tracks and propagate them to initial vertex position
  Int_t nTrks = PrepareTracks(*trkTree);
  delete trkTree;
  if(fDebug) printf(" tracks prepared: %d\n",nTrks);
  if(nTrks < fMinTracks) { TooFewTracks(); return fCurrentVertex; }

  // VERTEX FINDER
  VertexFinder();

  // VERTEX FITTER
  ComputeMaxChi2PerTrack(nTrks);
  VertexFitter();
  if(fDebug) printf(" vertex fit completed\n");

  return fCurrentVertex;
}
//----------------------------------------------------------------------------
AliESDVertex* AliITSVertexerTracks::FindVertexForCurrentEvent(AliESD *esdEvent)
{
//
// Vertex for current ESD event
//
  fCurrentVertex = 0;
  Double_t vtx[3],cvtx[6];

  // put tracks reco in ITS in a tree
  Int_t entr = (Int_t)esdEvent->GetNumberOfTracks();
  TTree *trkTree = new TTree("TreeT_ITS","its tracks");
  AliITStrackV2 *itstrack = 0;
  trkTree->Branch("tracks","AliITStrackV2",&itstrack,entr,0);

  for(Int_t i=0; i<entr; i++) {
    AliESDtrack *esdTrack = (AliESDtrack*)esdEvent->GetTrack(i);
    if(!esdTrack->GetStatus()&AliESDtrack::kITSin)
      { delete esdTrack; continue; }
    try {
      itstrack = new AliITStrackV2(*esdTrack);
    }
    catch (const Char_t *msg) {
        Warning("FindVertexForCurrentEvent",msg);
        delete esdTrack;
        continue;
    }

    trkTree->Fill();
    itstrack = 0;
    delete esdTrack; 
  }
  delete itstrack;

  // preselect tracks and propagate them to initial vertex position
  Int_t nTrks = PrepareTracks(*trkTree);
  delete trkTree;
  if(fDebug) printf(" tracks prepared: %d\n",nTrks);
  if(nTrks < fMinTracks) { TooFewTracks(); return fCurrentVertex; }

  // Set initial vertex position from ESD
  esdEvent->GetVertex()->GetXYZ(vtx);
  SetVtxStart(vtx[0],vtx[1]);

  // VERTEX FINDER
  VertexFinder();

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

  esdEvent->SetVertex(fCurrentVertex);

  cout<<"Vertex: "<<vtx[0]<<", "<<vtx[1]<<", "<<vtx[2]<<endl;
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
void AliITSVertexerTracks::VertexFinder() {

  // Get estimate of vertex position in (x,y) from tracks DCA
  // Then this estimate is stored to the data member fInitPos   
  // (previous values are overwritten)


 
  /*
******* TEMPORARY!!! FOR TEST ONLY!!! **********************************

fInitPos[0] = fNominalPos[0]+gRandom->Gaus(0.,0.0100); // 100 micron gaussian smearing
fInitPos[1] = fNominalPos[1]+gRandom->Gaus(0.,0.0100); // 100 micron gaussian smearing
  */

  fInitPos[2] = 0.;
  for(Int_t i=0;i<2;i++)fInitPos[i]=fNominalPos[i];

  Int_t nacc = (Int_t)fTrkArray.GetEntriesFast();

  Double_t aver[3]={0.,0.,0.};
  Int_t ncombi = 0;
  AliITStrackV2 *track1;
  AliITStrackV2 *track2;
  for(Int_t i=0; i<nacc; i++){
    track1 = (AliITStrackV2*)fTrkArray.At(i);
    if(fDebug>5){
      Double_t xv,par[5];
      track1->GetExternalParameters(xv,par);
      cout<<"Track in position "<<i<<" xr= "<<xv<<endl;
      for(Int_t ii=0;ii<5;ii++)cout<<par[ii]<<" ";
      cout<<endl;
    }
    Double_t mom1[3];
    Double_t alpha = track1->GetAlpha();
    Double_t azim = TMath::ASin(track1->GetSnp())+alpha;
    Double_t theta = TMath::Pi()/2. - TMath::ATan(track1->GetTgl());
    mom1[0] = TMath::Sin(theta)*TMath::Cos(azim);
    mom1[1] = TMath::Sin(theta)*TMath::Sin(azim);
    mom1[2] = TMath::Cos(theta);

    Double_t pos1[3];
    Double_t mindist = TMath::Cos(alpha)*fNominalPos[0]+TMath::Sin(alpha)*fNominalPos[1];
    track1->GetGlobalXYZat(mindist,pos1[0],pos1[1],pos1[2]);
    AliITSStrLine *line1 = new AliITSStrLine(pos1,mom1);
    for(Int_t j=i+1; j<nacc; j++){
      track2 = (AliITStrackV2*)fTrkArray.At(j);
      Double_t mom2[3];
      alpha = track2->GetAlpha();
      azim = TMath::ASin(track2->GetSnp())+alpha;
      theta = TMath::Pi()/2. - TMath::ATan(track2->GetTgl());
      mom2[0] = TMath::Sin(theta)*TMath::Cos(azim);
      mom2[1] = TMath::Sin(theta)*TMath::Sin(azim);
      mom2[2] = TMath::Cos(theta);
      Double_t pos2[3];
      mindist = TMath::Cos(alpha)*fNominalPos[0]+TMath::Sin(alpha)*fNominalPos[1];
      track2->GetGlobalXYZat(mindist,pos2[0],pos2[1],pos2[2]);
      AliITSStrLine *line2 = new AliITSStrLine(pos2,mom2);
      Double_t crosspoint[3];
      Int_t retcode = line2->Cross(line1,crosspoint);
      if(retcode<0){
	if(fDebug>10)cout<<" i= "<<i<<",   j= "<<j<<endl;
	if(fDebug>10)cout<<"bad intersection\n";
	line1->PrintStatus();
	line2->PrintStatus();
      }
      else {
	ncombi++;
	for(Int_t jj=0;jj<3;jj++)aver[jj]+=crosspoint[jj];
	if(fDebug>10)cout<<" i= "<<i<<",   j= "<<j<<endl;
	if(fDebug>10)cout<<"\n Cross point: ";
	if(fDebug>10)cout<<crosspoint[0]<<" "<<crosspoint[1]<<" "<<crosspoint[2]<<endl;
      }
      delete line2;
    }
    delete line1;
  }
  if(ncombi>0){
    for(Int_t jj=0;jj<3;jj++)fInitPos[jj] = aver[jj]/ncombi;
  }
  else {
    Warning("VertexFinder","Finder did not succed");
  }


  //************************************************************************
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


  Int_t i,j,k,step=0;
  TMatrixD rv(3,1);
  TMatrixD vV(3,3);
  rv(0,0) = fInitPos[0];
  rv(1,0) = fInitPos[1];
  rv(2,0) = 0.;
  Double_t xlStart,alpha;
  Double_t rotAngle;
  Double_t cosRot,sinRot;
  Double_t cc[15];
  Int_t nUsedTrks;
  Double_t chi2,chi2i;
  Int_t arrEntries = (Int_t)fTrkArray.GetEntries();
  AliITStrackV2 *t = 0;
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
      t = (AliITStrackV2*)fTrkArray.At(k);
      alpha = t->GetAlpha();
      xlStart = fInitPos[0]*TMath::Cos(alpha)+fInitPos[1]*TMath::Sin(alpha);
      t->PropagateTo(xlStart,0.,0.);   // to vtxSeed
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

  // VERTEX FINDER
  VertexFinder();

  // VERTEX FITTER
  ComputeMaxChi2PerTrack(nTrks);
  VertexFitter();
  if(fDebug) printf(" vertex fit completed\n");

  return fCurrentVertex;
}
//----------------------------------------------------------------------------



