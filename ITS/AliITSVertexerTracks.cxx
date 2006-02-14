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
  fDCAcut=0;
}
//----------------------------------------------------------------------------
AliITSVertexerTracks::AliITSVertexerTracks(TFile *inFile,TFile *outFile,
                                           const AliMagF *map,
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
  SetFieldMap(map);
  SetVtxStart(xStart,yStart);
  SetMinTracks();
  fTrksToSkip = 0;
  fNTrksToSkip = 0;
  fDCAcut=0;
  SetDebug();
}
//----------------------------------------------------------------------------
AliITSVertexerTracks::AliITSVertexerTracks(const AliMagF *map, TString fn,
					   Double_t xStart,Double_t yStart)
                                          :AliITSVertexer(fn) {
//
// Alternative constructor
//
  fInFile  = 0;
  fOutFile = 0;
  SetFieldMap(map);
  SetVtxStart(xStart,yStart);
  SetMinTracks();
  fTrksToSkip = 0;
  fNTrksToSkip = 0;
  fDCAcut=0;
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

    FindPrimaryVertexForCurrentEvent(ev);

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

    FindPrimaryVertexForCurrentEvent(esdEvent);

    if(!fCurrentVertex) {
      printf("AliITSVertexerTracks::FindVertixesESD():\n no vertex for event %d\n",ev);
      continue;
    }


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
  Int_t    nTrks    = 0;
  Bool_t   skipThis;
  Double_t d0rphi;
  Int_t    nEntries = (Int_t)trkTree.GetEntries();

  if(!fTrkArray.IsEmpty()) fTrkArray.Clear();
  fTrkArray.Expand(nEntries);

  if(fDebug) {
    printf(" PrepareTracks()\n");
    trkTree.Print();
  }
  cout<<" entr tree its tracks = "<<nEntries<<endl;
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
    d0rphi=Prepare(itstrack);
    if(d0rphi> maxd0rphi) { if(fDebug)cout<<"    !!!! d0rphi "<<d0rphi<<endl;continue; }
    fTrkArray.AddLast(itstrack);
    
    nTrks++; 
    if(fDebug)cout<<" :-) nTrks, d0rphi "<<nTrks<<"  "<<d0rphi<<endl;
   
  }
  if(fTrksToSkip) delete [] fTrksToSkip;

  return nTrks;
} 
//----------------------------------------------------------------------------
Int_t AliITSVertexerTracks::PrepareTracks(AliESD* esdEvent,Int_t nofCand, Int_t *trkPos) {
  //
  // Propagate tracks to initial vertex position and store them in a TObjArray
  //
  Int_t    nTrks    = 0;
  Double_t maxd0rphi = 3.; 
  Double_t d0rphi; 

  if(!fTrkArray.IsEmpty()) fTrkArray.Clear();
  fTrkArray.Expand(100);

  if(fDebug) {
    printf(" PrepareTracks()\n");
  }
  AliITStrackV2* itstrack;
  
  for(Int_t i=0; i<nofCand;i++){
    AliESDtrack *esdTrack = (AliESDtrack*)esdEvent->GetTrack(trkPos[i]);
    UInt_t status=esdTrack->GetStatus();
    if ((status&AliESDtrack::kTPCin)==0)continue;
    if ((status&AliESDtrack::kITSin)==0)continue;
    if ((status&AliESDtrack::kITSrefit)==0) continue;

    itstrack = new AliITStrackV2(*esdTrack);
    d0rphi=Prepare(itstrack);
    if(d0rphi> maxd0rphi) { if(fDebug)cout<<"    !!!! d0rphi "<<d0rphi<<endl;continue; }
    Int_t nclus=itstrack->GetNumberOfClusters();

    if(nclus<6){delete itstrack;continue;}
    fTrkArray.AddLast(itstrack);
    
    nTrks++; 
    if(fDebug)cout<<" :-) nTrks, d0rphi "<<nTrks<<"  "<<d0rphi<<endl;
    //delete itstrack;
  }
 

  if(fTrksToSkip) delete [] fTrksToSkip;
  return nTrks;
}
//----------------------------------------------------------------------------
Double_t AliITSVertexerTracks::Prepare(AliITStrackV2* itstrack){
 
  //
  Double_t alpha,xlStart,d0rphi; 
  // propagate track to vtxSeed
  alpha  = itstrack->GetAlpha();
  xlStart = fNominalPos[0]*TMath::Cos(alpha)+fNominalPos[1]*TMath::Sin(alpha);
  if(itstrack->GetX()>3.)itstrack->PropagateTo(3.,0.0023,65.19); // to beam pipe (0.8 mm of Be) 
  itstrack->PropagateTo(xlStart,0.,0.); 
  // select tracks with d0rphi < maxd0rphi
  d0rphi = TMath::Abs(itstrack->GetD(fNominalPos[0],fNominalPos[1]));
  return d0rphi;
          
}
//----------------------------------------------------------------------------
void AliITSVertexerTracks::PrintStatus() const {
//
// Print status
//
  printf(" Initial position (%f,%f)\n",fNominalPos[0],fNominalPos[1]);
  printf(" Vertex position after vertex finder:\n");
  fSimpVert.Print();
  printf(" Number of tracks in array: %d\n",(Int_t)fTrkArray.GetEntriesFast());
  printf(" Minimum # tracks required in fit: %d\n",fMinTracks);

  return;
}
//----------------------------------------------------------------------------
AliESDVertex* AliITSVertexerTracks::FindPrimaryVertexForCurrentEvent(Int_t evnumb) {
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
AliESDVertex* AliITSVertexerTracks::FindPrimaryVertexForCurrentEvent(AliESD *esdEvent)
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
        Warning("FindPrimaryVertexForCurrentEvent",msg);
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
//----------------------------------------------------------------------------
AliITSSimpleVertex* AliITSVertexerTracks::VertexForSelectedTracks(AliESD *esdEvent,Int_t nofCand, Int_t *trkPos,  Int_t opt){

  //
  // Computes the vertex for selected tracks 
  // trkPos=vector with track positions in ESD
  // values of opt -> see AliITSVertexerTracks.h
  //
  Double_t vtx[3]={0,0,0};

  Int_t nTrks = PrepareTracks(esdEvent,nofCand, trkPos);
  //delete trkTree;//  :-)) 
  if(fDebug) printf(" tracks prepared: %d\n",nTrks);
  if(nTrks < fMinTracks) {
    fSimpVert.SetXYZ(vtx);
    fSimpVert.SetDispersion(999);
    fSimpVert.SetNContributors(-5);
    return &fSimpVert;
  }
 
  // Set initial vertex position from ESD
  esdEvent->GetVertex()->GetXYZ(vtx);
  SetVtxStart(vtx[0],vtx[1]);
  if(opt==1)  StrLinVertexFinderMinDist(1);
  if(opt==2)  StrLinVertexFinderMinDist(0);
  if(opt==3)  HelixVertexFinder();
  if(opt==4)  VertexFinder(1);
  if(opt==5)  VertexFinder(0);
  return &fSimpVert;
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
void AliITSVertexerTracks::VertexFinder(Int_t OptUseWeights) {

  // Get estimate of vertex position in (x,y) from tracks DCA
  // Then this estimate is stored to the data member fSimpVert  
  // (previous values are overwritten)

 
  Double_t initPos[3];
  initPos[2] = 0.;
  for(Int_t i=0;i<2;i++)initPos[i]=fNominalPos[i];
  Int_t nacc = (Int_t)fTrkArray.GetEntriesFast();
  Double_t aver[3]={0.,0.,0.};
  Double_t aversq[3]={0.,0.,0.};
  Double_t sigmasq[3]={0.,0.,0.};
  Double_t sigma=0;
  Int_t ncombi = 0;
  AliITStrackV2 *track1;
  AliITStrackV2 *track2;
  Double_t alpha,mindist;

  for(Int_t i=0; i<nacc; i++){
    track1 = (AliITStrackV2*)fTrkArray.At(i);
    alpha=track1->GetAlpha();
    mindist = TMath::Cos(alpha)*fNominalPos[0]+TMath::Sin(alpha)*fNominalPos[1];
    AliITSStrLine *line1 = new AliITSStrLine();
    track1->ApproximateHelixWithLine(mindist,line1);
   
    if(fDebug>5){
      Double_t xv,par[5];
      track1->GetExternalParameters(xv,par);
      cout<<"Track in position "<<i<<" xr= "<<xv<<endl;
      for(Int_t ii=0;ii<5;ii++)cout<<par[ii]<<" ";
      cout<<endl;
    }

    for(Int_t j=i+1; j<nacc; j++){
      track2 = (AliITStrackV2*)fTrkArray.At(j);
      alpha=track2->GetAlpha();
      mindist = TMath::Cos(alpha)*fNominalPos[0]+TMath::Sin(alpha)*fNominalPos[1];
      AliITSStrLine *line2 = new AliITSStrLine();
      track2->ApproximateHelixWithLine(mindist,line2);
      Double_t distCA=line2->GetDCA(line1);
      if(fDCAcut<=0 || (fDCAcut>0&&distCA<fDCAcut)){
	Double_t pnt1[3],pnt2[3],crosspoint[3];

	if(OptUseWeights<=0){
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
	    for(Int_t jj=0;jj<3;jj++)aversq[jj]+=(crosspoint[jj]*crosspoint[jj]);
	    if(fDebug>10)cout<<" i= "<<i<<",   j= "<<j<<endl;
	    if(fDebug>10)cout<<"\n Cross point: ";
	    if(fDebug>10)cout<<crosspoint[0]<<" "<<crosspoint[1]<<" "<<crosspoint[2]<<endl;
	  }
	}
	if(OptUseWeights>0){
	  Int_t retcode = line1->CrossPoints(line2,pnt1,pnt2);
	  if(retcode>=0){
	    Double_t alpha, cs, sn;
	    alpha=track1->GetAlpha();
	    cs=TMath::Cos(alpha); sn=TMath::Sin(alpha);	  
	    Double_t sx1=sn*sn*track1->GetSigmaY2(), sy1=cs*cs*track1->GetSigmaY2();
	    alpha=track2->GetAlpha();
	    cs=TMath::Cos(alpha); sn=TMath::Sin(alpha);
	    Double_t sx2=sn*sn*track2->GetSigmaY2(), sy2=cs*cs*track2->GetSigmaY2();
	    Double_t sz1=track1->GetSigmaZ2(), sz2=track2->GetSigmaZ2();
	    Double_t wx1=sx2/(sx1+sx2), wx2=1.- wx1;
	    Double_t wy1=sy2/(sy1+sy2), wy2=1.- wy1;
	    Double_t wz1=sz2/(sz1+sz2), wz2=1.- wz1;
	    crosspoint[0]=wx1*pnt1[0] + wx2*pnt2[0]; 
	    crosspoint[1]=wy1*pnt1[1] + wy2*pnt2[1]; 
	    crosspoint[2]=wz1*pnt1[2] + wz2*pnt2[2];
	  
	    ncombi++;
	    for(Int_t jj=0;jj<3;jj++)aver[jj]+=crosspoint[jj];
	    for(Int_t jj=0;jj<3;jj++)aversq[jj]+=(crosspoint[jj]*crosspoint[jj]);
	  }
	}
      }
      delete line2;
    }
    delete line1;
  }
  if(ncombi>0){
    for(Int_t jj=0;jj<3;jj++){
      initPos[jj] = aver[jj]/ncombi;
      aversq[jj]/=ncombi;
      sigmasq[jj]=aversq[jj]-initPos[jj]*initPos[jj];
      sigma+=sigmasq[jj];
    }
    sigma=TMath::Sqrt(TMath::Abs(sigma));
  }
  else {
    Warning("VertexFinder","Finder did not succed");
    sigma=999;
  }
  fSimpVert.SetXYZ(initPos);
  fSimpVert.SetDispersion(sigma);
  fSimpVert.SetNContributors(ncombi);
}
//---------------------------------------------------------------------------
void AliITSVertexerTracks::HelixVertexFinder() {

  // Get estimate of vertex position in (x,y) from tracks DCA
  // Then this estimate is stored to the data member fSimpVert  
  // (previous values are overwritten)


  Double_t initPos[3];
  initPos[2] = 0.;
  for(Int_t i=0;i<2;i++)initPos[i]=fNominalPos[i];

  Int_t nacc = (Int_t)fTrkArray.GetEntriesFast();

  Double_t aver[3]={0.,0.,0.};
  Double_t averquad[3]={0.,0.,0.};
  Double_t sigmaquad[3]={0.,0.,0.};
  Double_t sigma=0;
  Int_t ncombi = 0;
  AliITStrackV2 *track1;
  AliITStrackV2 *track2;
  Double_t distCA;
  Double_t x, par[5];
  Double_t alpha, cs, sn;
  Double_t crosspoint[3];
  for(Int_t i=0; i<nacc; i++){
    track1 = (AliITStrackV2*)fTrkArray.At(i);
    

    for(Int_t j=i+1; j<nacc; j++){
      track2 = (AliITStrackV2*)fTrkArray.At(j);
      

      distCA=track2->PropagateToDCA(track1);

      if(fDCAcut<=0 ||(fDCAcut>0&&distCA<fDCAcut)){
	track1->GetExternalParameters(x,par);
	alpha=track1->GetAlpha();
	cs=TMath::Cos(alpha); sn=TMath::Sin(alpha);
	Double_t x1=x*cs - par[0]*sn;
	Double_t y1=x*sn + par[0]*cs;
	Double_t z1=par[1];
	Double_t sx1=sn*sn*track1->GetSigmaY2(), sy1=cs*cs*track1->GetSigmaY2(); 

	track2->GetExternalParameters(x,par);
	alpha=track2->GetAlpha();
	cs=TMath::Cos(alpha); sn=TMath::Sin(alpha);
	Double_t x2=x*cs - par[0]*sn;
	Double_t y2=x*sn + par[0]*cs;
	Double_t z2=par[1];
	Double_t sx2=sn*sn*track2->GetSigmaY2(), sy2=cs*cs*track2->GetSigmaY2(); 
	//	printf("Track %d pos=(%f,%f,%f) - dca=%f\n",i,x1,y1,z1,distCA);
	//printf("Track %d pos=(%f,%f,%f)\n",j,x2,y2,z2);

	Double_t sz1=track1->GetSigmaZ2(), sz2=track2->GetSigmaZ2();
	Double_t wx1=sx2/(sx1+sx2), wx2=1.- wx1;
	Double_t wy1=sy2/(sy1+sy2), wy2=1.- wy1;
	Double_t wz1=sz2/(sz1+sz2), wz2=1.- wz1;
	crosspoint[0]=wx1*x1 + wx2*x2; 
	crosspoint[1]=wy1*y1 + wy2*y2; 
	crosspoint[2]=wz1*z1 + wz2*z2;

	ncombi++;
	for(Int_t jj=0;jj<3;jj++)aver[jj]+=crosspoint[jj];
	for(Int_t jj=0;jj<3;jj++)averquad[jj]+=(crosspoint[jj]*crosspoint[jj]);
      }
    }
      
  }
  if(ncombi>0){
    for(Int_t jj=0;jj<3;jj++){
      initPos[jj] = aver[jj]/ncombi;
      averquad[jj]/=ncombi;
      sigmaquad[jj]=averquad[jj]-initPos[jj]*initPos[jj];
      sigma+=sigmaquad[jj];
    }
    sigma=TMath::Sqrt(TMath::Abs(sigma));
  }
  else {
    Warning("VertexFinder","Finder did not succed");
    sigma=999;
  }
  fSimpVert.SetXYZ(initPos);
  fSimpVert.SetDispersion(sigma);
  fSimpVert.SetNContributors(ncombi);
}
//---------------------------------------------------------------------------
void AliITSVertexerTracks::StrLinVertexFinderMinDist(Int_t OptUseWeights){

  // Calculate the point at minimum distance to prepared tracks 
  // Then this estimate is stored to the data member fSimpVert  
  // (previous values are overwritten)
  
  Double_t initPos[3];
  initPos[2] = 0.;
  Double_t sigma=0;
  for(Int_t i=0;i<2;i++)initPos[i]=fNominalPos[i];
  const Int_t knacc = (Int_t)fTrkArray.GetEntriesFast();

  AliITStrackV2 *track1;
  Double_t (*vectP0)[3]=new Double_t [knacc][3];
  Double_t (*vectP1)[3]=new Double_t [knacc][3];
  
  Double_t sum[3][3];
  Double_t dsum[3]={0,0,0};
  for(Int_t i=0;i<3;i++)
    for(Int_t j=0;j<3;j++)sum[i][j]=0;
  for(Int_t i=0; i<knacc; i++){
    track1 = (AliITStrackV2*)fTrkArray.At(i);
    Double_t alpha=track1->GetAlpha();
    Double_t mindist = TMath::Cos(alpha)*fNominalPos[0]+TMath::Sin(alpha)*fNominalPos[1];
    AliITSStrLine *line1 = new AliITSStrLine();
    track1->ApproximateHelixWithLine(mindist,line1);

    Double_t p0[3],cd[3];
    line1->GetP0(p0);
    line1->GetCd(cd);
    Double_t p1[3]={p0[0]+cd[0],p0[1]+cd[1],p0[2]+cd[2]};
    vectP0[i][0]=p0[0];
    vectP0[i][1]=p0[1];
    vectP0[i][2]=p0[2];
    vectP1[i][0]=p1[0];
    vectP1[i][1]=p1[1];
    vectP1[i][2]=p1[2];
    
    Double_t matr[3][3];
    Double_t dknow[3];
    if(OptUseWeights==0)GetStrLinDerivMatrix(p0,p1,matr,dknow);
    if(OptUseWeights==1){
      Double_t sigmasq[3];
      sigmasq[0]=track1->GetSigmaY2();
      sigmasq[1]=track1->GetSigmaY2();
      sigmasq[2]=track1->GetSigmaZ2();
      GetStrLinDerivMatrix(p0,p1,sigmasq,matr,dknow);
    }

    for(Int_t iii=0;iii<3;iii++){
      dsum[iii]+=dknow[iii]; 
      for(Int_t lj=0;lj<3;lj++) sum[iii][lj]+=matr[iii][lj];
    }
    delete line1;
  }
  
  Double_t vett[3][3];
  Double_t det=GetDeterminant3X3(sum);
  
   if(det!=0){
     for(Int_t zz=0;zz<3;zz++){
       for(Int_t ww=0;ww<3;ww++){
	 for(Int_t kk=0;kk<3;kk++) vett[ww][kk]=sum[ww][kk];
       }
       for(Int_t kk=0;kk<3;kk++) vett[kk][zz]=dsum[kk];
       initPos[zz]=GetDeterminant3X3(vett)/det;
     }


     for(Int_t i=0; i<knacc; i++){
       Double_t p0[3]={0,0,0},p1[3]={0,0,0};
       for(Int_t ii=0;ii<3;ii++){
	 p0[ii]=vectP0[i][ii];
	 p1[ii]=vectP1[i][ii];
       }
       sigma+=GetStrLinMinDist(p0,p1,initPos);
     }

     sigma=TMath::Sqrt(sigma);
   }else{
    Warning("VertexFinder","Finder did not succed");
    sigma=999;
  }
  delete vectP0;
  delete vectP1;
  fSimpVert.SetXYZ(initPos);
  fSimpVert.SetDispersion(sigma);
  fSimpVert.SetNContributors(knacc);
}
//_______________________________________________________________________
Double_t AliITSVertexerTracks::GetDeterminant3X3(Double_t matr[][3]){
  //
  Double_t det=matr[0][0]*matr[1][1]*matr[2][2]-matr[0][0]*matr[1][2]*matr[2][1]-matr[0][1]*matr[1][0]*matr[2][2]+matr[0][1]*matr[1][2]*matr[2][0]+matr[0][2]*matr[1][0]*matr[2][1]-matr[0][2]*matr[1][1]*matr[2][0];
 return det;
}
//____________________________________________________________________________
void AliITSVertexerTracks::GetStrLinDerivMatrix(Double_t *p0,Double_t *p1,Double_t m[][3],Double_t *d){

  //
  Double_t x12=p0[0]-p1[0];
  Double_t y12=p0[1]-p1[1];
  Double_t z12=p0[2]-p1[2];
  Double_t kk=x12*x12+y12*y12+z12*z12;
  m[0][0]=2-2/kk*x12*x12;
  m[0][1]=-2/kk*x12*y12;
  m[0][2]=-2/kk*x12*z12;
  m[1][0]=-2/kk*x12*y12;
  m[1][1]=2-2/kk*y12*y12;
  m[1][2]=-2/kk*y12*z12;
  m[2][0]=-2/kk*x12*z12;
  m[2][1]=-2*y12*z12;
  m[2][2]=2-2/kk*z12*z12;
  d[0]=2*p0[0]-2/kk*p0[0]*x12*x12-2/kk*p0[2]*x12*z12-2/kk*p0[1]*x12*y12;
  d[1]=2*p0[1]-2/kk*p0[1]*y12*y12-2/kk*p0[0]*x12*y12-2/kk*p0[2]*z12*y12;
  d[2]=2*p0[2]-2/kk*p0[2]*z12*z12-2/kk*p0[0]*x12*z12-2/kk*p0[1]*z12*y12;

}
//____________________________________________________________________________
void AliITSVertexerTracks::GetStrLinDerivMatrix(Double_t *p0,Double_t *p1,Double_t *sigmasq,Double_t m[][3],Double_t *d){
  //
  Double_t x12=p1[0]-p0[0];
  Double_t y12=p1[1]-p0[1];
  Double_t z12=p1[2]-p0[2];

  Double_t den= x12*x12*sigmasq[1]*sigmasq[2]+y12*y12*sigmasq[0]*sigmasq[2]+z12*z12*sigmasq[0]*sigmasq[1];

  Double_t kk= 2*(x12*x12/sigmasq[0]+y12*y12/sigmasq[1]+z12*z12/sigmasq[2]);

  Double_t cc[3];
  cc[0]=-x12/sigmasq[0];
  cc[1]=-y12/sigmasq[1];
  cc[2]=-z12/sigmasq[2];

  Double_t ww=(-p0[0]*x12*sigmasq[1]*sigmasq[2]-p0[1]*y12*sigmasq[0]*sigmasq[2]-p0[2]*z12*sigmasq[0]*sigmasq[1])/den;

  Double_t ss= -p0[0]*cc[0]-p0[1]*cc[1]-p0[2]*cc[2];

  Double_t aa[3];
  aa[0]=x12*sigmasq[1]*sigmasq[2]/den;
  aa[1]=y12*sigmasq[0]*sigmasq[2]/den;
  aa[2]=z12*sigmasq[0]*sigmasq[1]/den;

  m[0][0]=aa[0]*(aa[0]*kk+2*cc[0])+2*cc[0]*aa[0]+2/sigmasq[0];
  m[0][1]=aa[1]*(aa[0]*kk+2*cc[0])+2*cc[1]*aa[0];
  m[0][2]=aa[2]*(aa[0]*kk+2*cc[0])+2*cc[2]*aa[0];

  m[1][0]=aa[0]*(aa[1]*kk+2*cc[1])+2*cc[0]*aa[1];
  m[1][1]=aa[1]*(aa[1]*kk+2*cc[1])+2*cc[1]*aa[1]+2/sigmasq[1];
  m[1][2]=aa[2]*(aa[1]*kk+2*cc[1])+2*cc[2]*aa[1];

  m[2][0]=aa[0]*(aa[2]*kk+2*cc[2])+2*cc[0]*aa[2];
  m[2][1]=aa[1]*(aa[2]*kk+2*cc[2])+2*cc[1]*aa[2];
  m[2][2]=aa[2]*(aa[2]*kk+2*cc[2])+2*cc[2]*aa[2]+2/sigmasq[2];

  d[0]=-ww*(aa[0]*kk+2*cc[0])-2*ss*aa[0]+2*p0[0]/sigmasq[0];
  d[1]=-ww*(aa[1]*kk+2*cc[1])-2*ss*aa[1]+2*p0[1]/sigmasq[1];
  d[2]=-ww*(aa[2]*kk+2*cc[2])-2*ss*aa[2]+2*p0[2]/sigmasq[2];

}
//_____________________________________________________________________________
Double_t AliITSVertexerTracks::GetStrLinMinDist(Double_t *p0,Double_t *p1,Double_t *x0){
  //
  Double_t x12=p0[0]-p1[0];
  Double_t y12=p0[1]-p1[1];
  Double_t z12=p0[2]-p1[2];
  Double_t x10=p0[0]-x0[0];
  Double_t y10=p0[1]-x0[1];
  Double_t z10=p0[2]-x0[2];
  return ((x10*x10+y10*y10+z10*z10)*(x12*x12+y12*y12+z12*z12)-(x10*x12+y10*y12+z10*z12)*(x10*x12+y10*y12+z10*z12))/(x12*x12+y12*y12+z12*z12);

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
  Double_t initPos[3];
  fSimpVert.GetXYZ(initPos);
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
      xlStart = initPos[0]*TMath::Cos(alpha)+initPos[1]*TMath::Sin(alpha);
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



