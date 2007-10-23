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
//    Implementation of the vertexer from ESD tracks
//
// Origin: AliITSVertexerTracks
//         A.Dainese, Padova, 
//         andrea.dainese@pd.infn.it
//         M.Masera,  Torino, 
//         massimo.masera@to.infn.it
// Moved to STEER and adapted to ESD tracks: 
//          F.Prino,  Torino, prino@to.infn.it
//-----------------------------------------------------------------

//---- Root headers --------
#include <TSystem.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TTree.h>
#include <TMatrixD.h>
//---- AliRoot headers -----
#include "AliStrLine.h"
#include "AliVertexerTracks.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"

ClassImp(AliVertexerTracks)


//----------------------------------------------------------------------------
AliVertexerTracks::AliVertexerTracks():
TObject(),
fVert(),
fCurrentVertex(0),
fFieldkG(-999.),
fConstraint(kFALSE),
fOnlyFitter(kFALSE),
fMinTracks(1),
fMinITSClusters(5),
fTrkArray(),
fTrksToSkip(0),
fNTrksToSkip(0),
fDCAcut(0),
fAlgo(1),
fNSigma(3),
fMaxd0z0(0.5),
fITSin(kTRUE),
fITSrefit(kTRUE),
fnSigmaForUi00(1.5),
fDebug(0)
{
//
// Default constructor
//
  SetVtxStart();
  SetVtxStartSigma();
  SetMinTracks();
  SetMinITSClusters();
  SetNSigmad0(); 
  SetMaxd0z0(); 
}
//----------------------------------------------------------------------------
AliVertexerTracks::AliVertexerTracks(Double_t fieldkG):
TObject(),
fVert(),
fCurrentVertex(0),
fFieldkG(-999.),
fConstraint(kFALSE),
fOnlyFitter(kFALSE),
fMinTracks(1),
fMinITSClusters(5),
fTrkArray(),
fTrksToSkip(0),
fNTrksToSkip(0),
fDCAcut(0),
fAlgo(1),
fNSigma(3),
fMaxd0z0(0.5),
fITSin(kTRUE),
fITSrefit(kTRUE),
fnSigmaForUi00(1.5),
fDebug(0)
{
//
// Standard constructor
//
  SetVtxStart();
  SetVtxStartSigma();
  SetMinTracks();
  SetMinITSClusters();
  SetNSigmad0(); 
  SetMaxd0z0();
  SetFieldkG(fieldkG);
}
//-----------------------------------------------------------------------------
AliVertexerTracks::~AliVertexerTracks() 
{
  // Default Destructor
  // The objects pointed by the following pointer are not owned
  // by this class and are not deleted
  fCurrentVertex = 0;
  if(fTrksToSkip) { delete [] fTrksToSkip; fTrksToSkip=NULL; }
}
//----------------------------------------------------------------------------
AliESDVertex* AliVertexerTracks::FindPrimaryVertex(const AliESDEvent *esdEvent)
{
//
// Primary vertex for current ESD event
// (Two iterations: 
//  1st with 5*fNSigma*sigma cut w.r.t. to initial vertex
//      + cut on sqrt(d0d0+z0z0) if fConstraint=kFALSE  
//  2nd with fNSigma*sigma cut w.r.t. to vertex found in 1st iteration) 
// All ESD tracks with inside the beam pipe are then propagated to found vertex
//
  fCurrentVertex = 0;

  // accept 1-track case only if constraint is available
  if(!fConstraint && fMinTracks==1) fMinTracks=2;

  // read tracks from ESD
  Int_t nTrksTot = (Int_t)esdEvent->GetNumberOfTracks();
  if(nTrksTot<=0) {
    if(fDebug) printf("TooFewTracks\n");
    TooFewTracks(esdEvent);
    return fCurrentVertex;
  } 

  TDirectory * olddir = gDirectory;
  TFile f("VertexerTracks.root","recreate");
  TTree *trkTree = new TTree("TreeT","tracks");
  AliESDtrack *esdTrack = 0;
  trkTree->Branch("tracks","AliESDtrack",&esdTrack);

  Bool_t   skipThis;
  for(Int_t i=0; i<nTrksTot; i++) {
    AliESDtrack *et = esdEvent->GetTrack(i);
    // check tracks to skip
    skipThis = kFALSE;
    for(Int_t j=0; j<fNTrksToSkip; j++) { 
      if(et->GetID()==fTrksToSkip[j]) {
	if(fDebug) printf("skipping track: %d\n",i);
	skipThis = kTRUE;
      }
    }
    if(skipThis) continue;
    if(fITSin) {
      if(!(et->GetStatus()&AliESDtrack::kITSin)) continue;
      if(fITSrefit && !(et->GetStatus()&AliESDtrack::kITSrefit)) continue;
      Int_t nclus=et->GetNcls(0); // check number of clusters in ITS
      if(nclus<fMinITSClusters) continue;
    }
    esdTrack = new AliESDtrack(*et);
    trkTree->Fill();
    delete esdTrack;
  }
  
  // If fConstraint=kFALSE
  // run VertexFinder(1) to get rough estimate of initVertex (x,y)
  if(!fConstraint) {
    // fill fTrkArray, for VertexFinder()
    if(!fTrkArray.IsEmpty()) fTrkArray.Delete();
    PrepareTracks(*trkTree,0);
    Double_t cutsave = fDCAcut;  fDCAcut = 0.1; // 1 mm
    VertexFinder(1); // using weights, cutting dca < fDCAcut
    fDCAcut = cutsave;
    fTrkArray.Delete();
    if(fVert.GetNContributors()>0) {
      fVert.GetXYZ(fNominalPos);
      fNominalPos[0] = fVert.GetXv();
      fNominalPos[1] = fVert.GetYv();
      fNominalPos[2] = fVert.GetZv();
      if(fDebug) printf("No mean vertex: VertexFinder gives (%f, %f, %f)\n",fNominalPos[0],fNominalPos[1],fNominalPos[2]);
    } else {
      fNominalPos[0] = 0.;
      fNominalPos[1] = 0.;
      fNominalPos[2] = 0.;
      if(fDebug) printf("No mean vertex and VertexFinder failed\n");
    }
  }
  
  // TWO ITERATIONS:
  //
  // ITERATION 1
  // propagate tracks to fNominalPos vertex
  // preselect them:
  // if(constraint) reject for |d0|>5*fNSigma*sigma w.r.t. fNominal... vertex
  // else  reject for |d0|\oplus|z0| > 5 mm w.r.t. fNominal... vertex
  // ITERATION 2
  // propagate tracks to best between initVertex and fCurrentVertex
  // preselect tracks (reject for |d0|>fNSigma*sigma w.r.t. best 
  //                   between initVertex and fCurrentVertex) 
  for(Int_t iter=0; iter<2; iter++) {
    if(fOnlyFitter && iter==0) continue; 
    Int_t nTrksPrep = PrepareTracks(*trkTree,iter+1);
    if(fDebug) printf(" tracks prepared - iteration %d: %d\n",iter+1,nTrksPrep);
    if(nTrksPrep < fMinTracks) {
      if(fDebug) printf("TooFewTracks\n");
      TooFewTracks(esdEvent);
      if(fDebug) fCurrentVertex->PrintStatus();
      fTrkArray.Delete();
      delete trkTree;
      f.Close();
      gSystem->Unlink("VertexerTracks.root");
      olddir->cd();
      return fCurrentVertex; 
    }

    // vertex finder
    if(!fOnlyFitter) {
      if(nTrksPrep==1){
	if(fDebug) printf("Just one track\n");
	OneTrackVertFinder();
      }else{
	switch (fAlgo) {
        case 1: StrLinVertexFinderMinDist(1); break;
        case 2: StrLinVertexFinderMinDist(0); break;
        case 3: HelixVertexFinder();          break;
        case 4: VertexFinder(1);              break;
        case 5: VertexFinder(0);              break;
        default: printf("Wrong algorithm\n"); break;  
	}
      }
      if(fDebug) printf(" Vertex finding completed\n");
    }

    // vertex fitter
    VertexFitter(fConstraint);
    if(fDebug) printf(" Vertex fit completed\n");
    if(iter==0) fTrkArray.Delete();
  } // end loop on the two iterations


  if(fConstraint) {
    if(fOnlyFitter) {
      fCurrentVertex->SetTitle("VertexerTracksWithConstraintOnlyFitter");
    } else {
      fCurrentVertex->SetTitle("VertexerTracksWithConstraint");
    }
  } else {
    fCurrentVertex->SetTitle("VertexerTracksNoConstraint");
  }

    
  // set indices of used tracks
  UShort_t *indices = 0;
  AliESDtrack *ett = 0;
  if(fCurrentVertex->GetNContributors()>0) {
    indices = new UShort_t[fCurrentVertex->GetNContributors()];
    for(Int_t jj=0;jj<(Int_t)fTrkArray.GetEntriesFast();jj++) {
      ett = (AliESDtrack*)fTrkArray.At(jj);
      indices[jj] = (UShort_t)ett->GetID();
    }
    fCurrentVertex->SetIndices(fCurrentVertex->GetNContributors(),indices);
  }
  delete [] indices;

  delete trkTree;
  f.Close();
  gSystem->Unlink("VertexerTracks.root");
  olddir->cd();

  fTrkArray.Delete();

  if(fTrksToSkip) { delete [] fTrksToSkip; fTrksToSkip=NULL; }


  if(fDebug) fCurrentVertex->PrintStatus();
  if(fDebug) fCurrentVertex->PrintIndices();

  
  return fCurrentVertex;
}
//------------------------------------------------------------------------
Double_t AliVertexerTracks::GetDeterminant3X3(Double_t matr[][3])
{
  //
  Double_t det=matr[0][0]*matr[1][1]*matr[2][2]-matr[0][0]*matr[1][2]*matr[2][1]-matr[0][1]*matr[1][0]*matr[2][2]+matr[0][1]*matr[1][2]*matr[2][0]+matr[0][2]*matr[1][0]*matr[2][1]-matr[0][2]*matr[1][1]*matr[2][0];
 return det;
}
//-------------------------------------------------------------------------
void AliVertexerTracks::GetStrLinDerivMatrix(Double_t *p0,Double_t *p1,Double_t (*m)[3],Double_t *d)
{
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
//--------------------------------------------------------------------------  
void AliVertexerTracks::GetStrLinDerivMatrix(Double_t *p0,Double_t *p1,Double_t *sigmasq,Double_t (*m)[3],Double_t *d)
{
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
//--------------------------------------------------------------------------   
Double_t AliVertexerTracks::GetStrLinMinDist(Double_t *p0,Double_t *p1,Double_t *x0)
{
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
void AliVertexerTracks::OneTrackVertFinder() 
{
  // find vertex for events with 1 track, using DCA to nominal beam axis
  if(fDebug) printf("Number of prepared tracks =%d - Call OneTrackVertFinder",fTrkArray.GetEntries());
  AliESDtrack *track1;
  track1 = (AliESDtrack*)fTrkArray.At(0);
  Double_t field=GetFieldkG();
  Double_t alpha=track1->GetAlpha();
  Double_t mindist = TMath::Cos(alpha)*fNominalPos[0]+TMath::Sin(alpha)*fNominalPos[1];
  Double_t pos[3],dir[3]; 
  track1->GetXYZAt(mindist,field,pos);
  track1->GetPxPyPzAt(mindist,field,dir);
  AliStrLine *line1 = new AliStrLine(pos,dir);
  Double_t p1[3]={fNominalPos[0],fNominalPos[1],0.}; 
  Double_t p2[3]={fNominalPos[0],fNominalPos[1],10.}; 
  AliStrLine *zeta=new AliStrLine(p1,p2,kTRUE);
  Double_t crosspoint[3]={0.,0.,0.};
  Double_t sigma=999.;
  Int_t nContrib=-1;
  Int_t retcode = zeta->Cross(line1,crosspoint);
  if(retcode>=0){
    sigma=line1->GetDistFromPoint(crosspoint);
    nContrib=1;
  }
  delete zeta;
  delete line1;
  fVert.SetXYZ(crosspoint);
  fVert.SetDispersion(sigma);
  fVert.SetNContributors(nContrib);  
}
//---------------------------------------------------------------------------
void AliVertexerTracks::HelixVertexFinder() 
{
  // Get estimate of vertex position in (x,y) from tracks DCA


  Double_t initPos[3];
  initPos[2] = 0.;
  for(Int_t i=0;i<2;i++)initPos[i]=fNominalPos[i];
  Double_t field=GetFieldkG();

  Int_t nacc = (Int_t)fTrkArray.GetEntriesFast();

  Double_t aver[3]={0.,0.,0.};
  Double_t averquad[3]={0.,0.,0.};
  Double_t sigmaquad[3]={0.,0.,0.};
  Double_t sigma=0;
  Int_t ncombi = 0;
  AliESDtrack *track1;
  AliESDtrack *track2;
  Double_t distCA;
  Double_t x, par[5];
  Double_t alpha, cs, sn;
  Double_t crosspoint[3];
  for(Int_t i=0; i<nacc; i++){
    track1 = (AliESDtrack*)fTrkArray.At(i);
    

    for(Int_t j=i+1; j<nacc; j++){
      track2 = (AliESDtrack*)fTrkArray.At(j);

      distCA=track2->PropagateToDCA(track1,field);
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
    Warning("HelixVertexFinder","Finder did not succed");
    sigma=999;
  }
  fVert.SetXYZ(initPos);
  fVert.SetDispersion(sigma);
  fVert.SetNContributors(ncombi);
}
//----------------------------------------------------------------------------
Int_t AliVertexerTracks::PrepareTracks(TTree &trkTree,Int_t optImpParCut) 
{
//
// Propagate tracks to initial vertex position and store them in a TObjArray
//
  Double_t maxd0rphi; 
  Double_t maxd0z0 = fMaxd0z0; // default is 5 mm  
  Int_t    nTrks    = 0;
  Double_t sigmaCurr[3];
  Double_t normdistx,normdisty;
  Float_t  d0z0[2],covd0z0[3]; 
  Double_t sigma;
  Double_t field=GetFieldkG();

  AliESDVertex *initVertex = new AliESDVertex(fNominalPos,fNominalCov,1,1);

  Int_t    nEntries = (Int_t)trkTree.GetEntries();
  if(!fTrkArray.IsEmpty()) fTrkArray.Delete();

  if(fDebug) {
    printf(" PrepareTracks()\n");
  }

  for(Int_t i=0; i<nEntries; i++) {
    AliESDtrack *track = new AliESDtrack; 
    trkTree.SetBranchAddress("tracks",&track);
    trkTree.GetEvent(i);

    // propagate track to vertex
    if(optImpParCut<=1 || fOnlyFitter) { // optImpParCut==1 or 0
      track->RelateToVertex(initVertex,field,100.);
    } else {              // optImpParCut==2
      fCurrentVertex->GetSigmaXYZ(sigmaCurr);
      normdistx = TMath::Abs(fCurrentVertex->GetXv()-fNominalPos[0])/TMath::Sqrt(sigmaCurr[0]*sigmaCurr[0]+fNominalCov[0]);
      normdisty = TMath::Abs(fCurrentVertex->GetYv()-fNominalPos[1])/TMath::Sqrt(sigmaCurr[1]*sigmaCurr[1]+fNominalCov[2]);
      if(normdistx < 3. && normdisty < 3. &&
	 (sigmaCurr[0]+sigmaCurr[1])<(TMath::Sqrt(fNominalCov[0])+TMath::Sqrt(fNominalCov[2]))) {
	track->RelateToVertex(fCurrentVertex,field,100.);
      } else {
	track->RelateToVertex(initVertex,field,100.);
      }
    }

    track->GetImpactParameters(d0z0,covd0z0);
    sigma = TMath::Sqrt(covd0z0[0]);
    maxd0rphi = fNSigma*sigma;
    if(optImpParCut==1) maxd0rphi *= 5.;



    if(fDebug) printf("trk %d; lab %d; |d0| = %f;  d0 cut = %f; |z0| = %f; |d0|oplus|z0| = %f; d0z0 cut = %f\n",i,track->GetLabel(),TMath::Abs(d0z0[0]),maxd0rphi,TMath::Abs(d0z0[1]),TMath::Sqrt(d0z0[0]*d0z0[0]+d0z0[1]*d0z0[1]),maxd0z0);

    // during iterations 1 and 2, if fConstraint=kFALSE,
    // select tracks with d0oplusz0 < maxd0z0
    if(optImpParCut>=1 && !fConstraint && nEntries>=3 && 
       fVert.GetNContributors()>0) {
      if(TMath::Sqrt(d0z0[0]*d0z0[0]+d0z0[1]*d0z0[1]) > maxd0z0) { 
	if(fDebug) printf("     rejected\n");
	delete track; continue; 
      }
    }

    // select tracks with d0rphi < maxd0rphi
    if(optImpParCut>0 && TMath::Abs(d0z0[0]) > maxd0rphi) { 
      if(fDebug) printf("     rejected\n");
      delete track; continue; 
    }

    fTrkArray.AddLast(track);
    nTrks++; 
  }

  delete initVertex;

  return nTrks;
} 
//---------------------------------------------------------------------------
AliESDVertex* AliVertexerTracks::RemoveTracksFromVertex(AliESDVertex *inVtx,
							TTree *trksTree,
							Float_t *diamondxy) 
{
//
// Removes tracks in trksTree from fit of inVtx
//

  if(!strstr(inVtx->GetTitle(),"VertexerTracksWithConstraint")) {
    printf("ERROR: primary vertex has no constraint: cannot remove tracks\n");
    return 0x0;
  }
  if(!strstr(inVtx->GetTitle(),"VertexerTracksWithConstraintOnlyFitter"))
    printf("WARNING: result of tracks' removal will be only approximately correct\n");

  TMatrixD rv(3,1);
  rv(0,0) = inVtx->GetXv();
  rv(1,0) = inVtx->GetYv();
  rv(2,0) = inVtx->GetZv();
  TMatrixD vV(3,3);
  Double_t cov[6];
  inVtx->GetCovMatrix(cov);
  vV(0,0) = cov[0];
  vV(0,1) = cov[1]; vV(1,0) = cov[1];
  vV(1,1) = cov[2];
  vV(0,2) = cov[3]; vV(2,0) = cov[3];
  vV(1,2) = cov[4]; vV(2,1) = cov[4]; 
  vV(2,2) = cov[5];

  TMatrixD sumWi(TMatrixD::kInverted,vV);
  TMatrixD sumWiri(sumWi,TMatrixD::kMult,rv);

  Int_t nUsedTrks = inVtx->GetNContributors();
  Double_t chi2 = inVtx->GetChi2();

  AliESDtrack *track = 0;
  trksTree->SetBranchAddress("tracks",&track);
  Int_t ntrks = trksTree->GetEntries();
  for(Int_t i=0;i<ntrks;i++) {
    trksTree->GetEvent(i);
    if(!inVtx->UsesTrack(track->GetID())) {
      printf("track %d was not used in vertex fit\n",track->GetID());
      continue;
    }
    Double_t alpha = track->GetAlpha();
    Double_t xl = diamondxy[0]*TMath::Cos(alpha)+diamondxy[1]*TMath::Sin(alpha);
    track->AliExternalTrackParam::PropagateTo(xl,GetFieldkG()); 
    // vector of track global coordinates
    TMatrixD ri(3,1);
    // covariance matrix of ri
    TMatrixD wWi(3,3);
    
    // get space point from track
    if(!TrackToPoint(track,ri,wWi)) continue;

    TMatrixD wWiri(wWi,TMatrixD::kMult,ri); 

    sumWi -= wWi;
    sumWiri -= wWiri;

    // track chi2
    TMatrixD deltar = rv; deltar -= ri;
    TMatrixD wWideltar(wWi,TMatrixD::kMult,deltar);
    Double_t chi2i = deltar(0,0)*wWideltar(0,0)+
                     deltar(1,0)*wWideltar(1,0)+
	             deltar(2,0)*wWideltar(2,0);
    // remove from total chi2
    chi2 -= chi2i;

    nUsedTrks--;
    if(nUsedTrks<2) {
      printf("Trying to remove too many tracks!\n");
      return 0x0;
    }
  }

  TMatrixD rvnew(3,1);
  TMatrixD vVnew(3,3);

  // new inverted of weights matrix
  TMatrixD invsumWi(TMatrixD::kInverted,sumWi);
  vVnew = invsumWi;
  // new position of primary vertex
  rvnew.Mult(vVnew,sumWiri);

  Double_t position[3];
  position[0] = rvnew(0,0);
  position[1] = rvnew(1,0);
  position[2] = rvnew(2,0);
  cov[0] = vVnew(0,0);
  cov[1] = vVnew(0,1);
  cov[2] = vVnew(1,1);
  cov[3] = vVnew(0,2);
  cov[4] = vVnew(1,2);
  cov[5] = vVnew(2,2);
  
  // store data in the vertex object
  AliESDVertex *outVtx = new AliESDVertex(position,cov,chi2,nUsedTrks);
  outVtx->SetTitle(inVtx->GetTitle());
  UShort_t *inindices = inVtx->GetIndices();
  UShort_t *outindices = new UShort_t[outVtx->GetNContributors()];
  Int_t j=0;
  Bool_t copyindex;
  for(Int_t k=0; k<inVtx->GetNIndices(); k++) {
    copyindex=kTRUE;
    for(Int_t l=0; l<ntrks; l++) {
      trksTree->GetEvent(l);
      if(inindices[k]==track->GetID()) copyindex=kFALSE;
    }
    if(copyindex) {
      outindices[j] = inindices[k]; j++;
    }
  }
  outVtx->SetIndices(outVtx->GetNContributors(),outindices);
  delete [] outindices;

  if(fDebug) {
    printf("Vertex before removing tracks:\n");
    inVtx->PrintStatus();
    inVtx->PrintIndices();
    printf("Vertex after removing tracks:\n");
    outVtx->PrintStatus();
    outVtx->PrintIndices();
  }

  return outVtx;
}
//---------------------------------------------------------------------------
void AliVertexerTracks::SetSkipTracks(Int_t n,Int_t *skipped) 
{
//
// Mark the tracks not to be used in the vertex reconstruction.
// Tracks are identified by AliESDtrack::GetID()
//
  fNTrksToSkip = n;  fTrksToSkip = new Int_t[n]; 
  for(Int_t i=0;i<n;i++) fTrksToSkip[i] = skipped[i]; 
  return; 
}
//---------------------------------------------------------------------------
void  AliVertexerTracks::SetVtxStart(AliESDVertex *vtx) 
{ 
//
// Set initial vertex knowledge
//
  vtx->GetXYZ(fNominalPos);
  vtx->GetCovMatrix(fNominalCov);
  SetConstraintOn();
  return; 
}
//---------------------------------------------------------------------------
void AliVertexerTracks::StrLinVertexFinderMinDist(Int_t optUseWeights)
{
  AliESDtrack *track1;
  Double_t field=GetFieldkG();
  const Int_t knacc = (Int_t)fTrkArray.GetEntriesFast();
  TClonesArray *linarray = new TClonesArray("AliStrLine",1000);
  TClonesArray &lines = *linarray;
  for(Int_t i=0; i<knacc; i++){
    track1 = (AliESDtrack*)fTrkArray.At(i);
    Double_t alpha=track1->GetAlpha();
    Double_t mindist = TMath::Cos(alpha)*fNominalPos[0]+TMath::Sin(alpha)*fNominalPos[1];
    Double_t pos[3],dir[3],sigmasq[3]; 
    track1->GetXYZAt(mindist,field,pos);
    track1->GetPxPyPzAt(mindist,field,dir);
    sigmasq[0]=TMath::Sin(alpha)*TMath::Sin(alpha)*track1->GetSigmaY2();
    sigmasq[1]=TMath::Cos(alpha)*TMath::Cos(alpha)*track1->GetSigmaY2();
    sigmasq[2]=track1->GetSigmaZ2();
    TMatrixD ri(3,1);
    TMatrixD wWi(3,3);
    if(!TrackToPoint(track1,ri,wWi)) continue;
    Double_t wmat[9];
    Int_t iel=0;
    for(Int_t ia=0;ia<3;ia++){
      for(Int_t ib=0;ib<3;ib++){
	wmat[iel]=wWi(ia,ib);
	iel++;
      }    
    }
    new(lines[i]) AliStrLine(pos,sigmasq,wmat,dir);     
  }
  fVert=TrackletVertexFinder(linarray,optUseWeights);
  linarray->Delete();
  delete linarray;
}
//---------------------------------------------------------------------------
AliESDVertex AliVertexerTracks::TrackletVertexFinder(TClonesArray *lines, Int_t optUseWeights)
{
  // Calculate the point at minimum distance to prepared tracks 
  
  const Int_t knacc = (Int_t)lines->GetEntriesFast();
  Double_t initPos[3]={0.,0.,0.};

  Double_t (*vectP0)[3]=new Double_t [knacc][3];
  Double_t (*vectP1)[3]=new Double_t [knacc][3];
  
  Double_t sum[3][3];
  Double_t dsum[3]={0,0,0};
  TMatrixD sumWi(3,3);
  for(Int_t i=0;i<3;i++){
    for(Int_t j=0;j<3;j++){
      sum[i][j]=0;
      sumWi(i,j)=0.;
    }
  }

  for(Int_t i=0; i<knacc; i++){
    AliStrLine* line1 = (AliStrLine*)lines->At(i); 
    Double_t p0[3],cd[3],sigmasq[3];
    Double_t wmat[9];
    line1->GetP0(p0);
    line1->GetCd(cd);
    line1->GetSigma2P0(sigmasq);
    line1->GetWMatrix(wmat);
    TMatrixD wWi(3,3);
    Int_t iel=0;
    for(Int_t ia=0;ia<3;ia++){
      for(Int_t ib=0;ib<3;ib++){
	wWi(ia,ib)=wmat[iel];
	iel++;
      }    
    }

    sumWi+=wWi;

    Double_t p1[3]={p0[0]+cd[0],p0[1]+cd[1],p0[2]+cd[2]};
    vectP0[i][0]=p0[0];
    vectP0[i][1]=p0[1];
    vectP0[i][2]=p0[2];
    vectP1[i][0]=p1[0];
    vectP1[i][1]=p1[1];
    vectP1[i][2]=p1[2];
    
    Double_t matr[3][3];
    Double_t dknow[3];
    if(optUseWeights==0)GetStrLinDerivMatrix(p0,p1,matr,dknow);
    else GetStrLinDerivMatrix(p0,p1,sigmasq,matr,dknow);


    for(Int_t iii=0;iii<3;iii++){
      dsum[iii]+=dknow[iii]; 
      for(Int_t lj=0;lj<3;lj++) sum[iii][lj]+=matr[iii][lj];
    }
  }
 
  TMatrixD invsumWi(TMatrixD::kInverted,sumWi);
  Double_t covmatrix[6];
  covmatrix[0] = invsumWi(0,0);
  covmatrix[1] = invsumWi(0,1);
  covmatrix[2] = invsumWi(1,1);
  covmatrix[3] = invsumWi(0,2);
  covmatrix[4] = invsumWi(1,2);
  covmatrix[5] = invsumWi(2,2);

  Double_t vett[3][3];
  Double_t det=GetDeterminant3X3(sum);
  Double_t sigma=0;
  
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
    sigma=999;
  }
  AliESDVertex theVert(initPos,covmatrix,99999.,knacc);
  theVert.SetDispersion(sigma);
  delete vectP0;
  delete vectP1;
  return theVert;
}
//---------------------------------------------------------------------------
Bool_t AliVertexerTracks::TrackToPoint(AliESDtrack *t,
				       TMatrixD &ri,TMatrixD &wWi,
				       Bool_t uUi3by3) const 
{
//
// Extract from the AliESDtrack the global coordinates ri and covariance matrix
// wWi of the space point that it represents (to be used in VertexFitter())
//

  
  Double_t rotAngle = t->GetAlpha();
  if(rotAngle<0.) rotAngle += 2.*TMath::Pi();
  Double_t cosRot = TMath::Cos(rotAngle);
  Double_t sinRot = TMath::Sin(rotAngle);

  ri(0,0) = t->GetX()*cosRot-t->GetY()*sinRot;
  ri(1,0) = t->GetX()*sinRot+t->GetY()*cosRot;
  ri(2,0) = t->GetZ();



  if(!uUi3by3) {
    // matrix to go from global (x,y,z) to local (y,z);
    TMatrixD qQi(2,3);
    qQi(0,0) = -sinRot;
    qQi(0,1) = cosRot;
    qQi(0,2) = 0.;
    qQi(1,0) = 0.;
    qQi(1,1) = 0.;
    qQi(1,2) = 1.;

    // covariance matrix of local (y,z) - inverted
    Double_t cc[15];
    t->GetExternalCovariance(cc);
    TMatrixD uUi(2,2);
    uUi(0,0) = cc[0];
    uUi(0,1) = cc[1];
    uUi(1,0) = cc[1];
    uUi(1,1) = cc[2];
    //printf(" Ui :\n");
    //printf(" %f   %f\n",uUi(0,0),uUi(0,1));
    //printf(" %f   %f\n",uUi(1,0),uUi(1,1));

    if(uUi.Determinant() <= 0.) return kFALSE;
    TMatrixD uUiInv(TMatrixD::kInverted,uUi);
  
    // weights matrix: wWi = qQiT * uUiInv * qQi
    TMatrixD uUiInvQi(uUiInv,TMatrixD::kMult,qQi);
    TMatrixD m(qQi,TMatrixD::kTransposeMult,uUiInvQi);
    wWi = m;
  } else {
    if(fVert.GetNContributors()<1) AliFatal("Vertex from finder is empty");
    // matrix to go from global (x,y,z) to local (x,y,z);
    TMatrixD qQi(3,3);
    qQi(0,0) = cosRot;
    qQi(0,1) = sinRot;
    qQi(0,2) = 0.;
    qQi(1,0) = -sinRot;
    qQi(1,1) = cosRot;
    qQi(1,2) = 0.;
    qQi(2,0) = 0.;
    qQi(2,1) = 0.;
    qQi(2,2) = 1.;
   
    // covariance of fVert along the track  
    Double_t p[3],pt,ptot;
    t->GetPxPyPz(p);
    pt = TMath::Sqrt(p[0]*p[0]+p[1]*p[1]);
    ptot = TMath::Sqrt(pt*pt+p[2]*p[2]);
    Double_t cphi = p[0]/pt;               //cos(phi)=px/pt
    Double_t sphi = p[1]/pt;               //sin(phi)=py/pt
    Double_t clambda = pt/ptot;            //cos(lambda)=pt/ptot
    Double_t slambda = p[2]/ptot;            //sin(lambda)=pz/ptot
    Double_t covfVert[6];
    fVert.GetCovMatrix(covfVert);
    Double_t covfVertalongt = 
       covfVert[0]*cphi*cphi*clambda*clambda 
      +covfVert[1]*2.*cphi*sphi*clambda*clambda
      +covfVert[3]*2.*cphi*clambda*slambda 
      +covfVert[2]*sphi*sphi*clambda*clambda 
      +covfVert[4]*2.*sphi*clambda*slambda 
      +covfVert[5]*slambda*slambda; 
    Double_t cc[15];
    t->GetExternalCovariance(cc);
    // covariance matrix of local (x,y,z) - inverted
    TMatrixD uUi(3,3);
    uUi(0,0) = covfVertalongt * fnSigmaForUi00 * fnSigmaForUi00;
    if(fDebug) printf("=====> sqrtUi00 cm  %f\n",TMath::Sqrt(uUi(0,0)));
    uUi(0,1) = 0.;
    uUi(0,2) = 0.;
    uUi(1,0) = 0.;
    uUi(1,1) = cc[0];
    uUi(1,2) = cc[1];
    uUi(2,0) = 0.;
    uUi(2,1) = cc[1];
    uUi(2,2) = cc[2];
    //printf(" Ui :\n");
    //printf(" %f   %f\n",uUi(0,0),uUi(0,1));
    //printf(" %f   %f\n",uUi(1,0),uUi(1,1));
  
    if(uUi.Determinant() <= 0.) return kFALSE;
    TMatrixD uUiInv(TMatrixD::kInverted,uUi);
  
    // weights matrix: wWi = qQiT * uUiInv * qQi
    TMatrixD uUiInvQi(uUiInv,TMatrixD::kMult,qQi);
    TMatrixD m(qQi,TMatrixD::kTransposeMult,uUiInvQi);
    wWi = m;
  }


  return kTRUE;
} 
//---------------------------------------------------------------------------
void AliVertexerTracks::TooFewTracks(const AliESDEvent* esdEvent) 
{
//
// When the number of tracks is < fMinTracks
//

  // deal with vertices not found
  Double_t pos[3],err[3];
  Int_t    ncontr=0;
  pos[0] = fNominalPos[0];
  err[0] = TMath::Sqrt(fNominalCov[0]);
  pos[1] = fNominalPos[1];
  err[1] = TMath::Sqrt(fNominalCov[2]);
  pos[2] = esdEvent->GetVertex()->GetZv();
  err[2] = esdEvent->GetVertex()->GetZRes();
  if(err[0]>1. && esdEvent->GetVertex()->GetNContributors()<=0) 
    ncontr = -1; // (x,y,z) = (0,0,0) 
  if(err[0]>1. && esdEvent->GetVertex()->GetNContributors()>0) 
    ncontr = -2; // (x,y,z) = (0,0,z_fromSPD) 
  if(err[0]<1. && esdEvent->GetVertex()->GetNContributors()<=0) 
    ncontr = -3; // (x,y,z) = (x_mean,y_mean,0)
  if(err[0]<1. && esdEvent->GetVertex()->GetNContributors()>0) 
    ncontr = -4; // (x,y,z) = (x_mean,y_mean,z_fromSPD)
  fCurrentVertex = 0;
  fCurrentVertex = new AliESDVertex(pos,err);
  fCurrentVertex->SetNContributors(ncontr);

  if(fConstraint) {
    fCurrentVertex->SetTitle("VertexerTracksWithConstraint");
  } else {
    fCurrentVertex->SetTitle("VertexerTracksNoConstraint");
  }

  return;
}
//---------------------------------------------------------------------------
void AliVertexerTracks::VertexFinder(Int_t optUseWeights) 
{

  // Get estimate of vertex position in (x,y) from tracks DCA
 
  Double_t initPos[3];
  initPos[2] = 0.;
  for(Int_t i=0;i<2;i++)initPos[i]=fNominalPos[i];
  Int_t nacc = (Int_t)fTrkArray.GetEntriesFast();
  Double_t aver[3]={0.,0.,0.};
  Double_t aversq[3]={0.,0.,0.};
  Double_t sigmasq[3]={0.,0.,0.};
  Double_t sigma=0;
  Int_t ncombi = 0;
  AliESDtrack *track1;
  AliESDtrack *track2;
  Double_t pos[3],dir[3]; 
  Double_t alpha,mindist;
  Double_t field=GetFieldkG();

  for(Int_t i=0; i<nacc; i++){
    track1 = (AliESDtrack*)fTrkArray.At(i);
    alpha=track1->GetAlpha();
    mindist = TMath::Cos(alpha)*fNominalPos[0]+TMath::Sin(alpha)*fNominalPos[1];
    track1->GetXYZAt(mindist,field,pos);
    track1->GetPxPyPzAt(mindist,field,dir);
    AliStrLine *line1 = new AliStrLine(pos,dir); 

   //    AliStrLine *line1 = new AliStrLine();
   //    track1->ApproximateHelixWithLine(mindist,field,line1);
   
    for(Int_t j=i+1; j<nacc; j++){
      track2 = (AliESDtrack*)fTrkArray.At(j);
      alpha=track2->GetAlpha();
      mindist = TMath::Cos(alpha)*fNominalPos[0]+TMath::Sin(alpha)*fNominalPos[1];
      track2->GetXYZAt(mindist,field,pos);
      track2->GetPxPyPzAt(mindist,field,dir);
      AliStrLine *line2 = new AliStrLine(pos,dir); 
    //      AliStrLine *line2 = new AliStrLine();
    //  track2->ApproximateHelixWithLine(mindist,field,line2);
      Double_t distCA=line2->GetDCA(line1);
      //printf("%d   %d   %f\n",i,j,distCA);
       if(fDCAcut<=0 || (fDCAcut>0&&distCA<fDCAcut)){
	Double_t pnt1[3],pnt2[3],crosspoint[3];

	if(optUseWeights<=0){
	  Int_t retcode = line2->Cross(line1,crosspoint);
	  if(retcode>=0){
	    ncombi++;
	    for(Int_t jj=0;jj<3;jj++)aver[jj]+=crosspoint[jj];
	    for(Int_t jj=0;jj<3;jj++)aversq[jj]+=(crosspoint[jj]*crosspoint[jj]);
	  }
	}
	if(optUseWeights>0){
	  Int_t retcode = line1->CrossPoints(line2,pnt1,pnt2);
	  if(retcode>=0){
	    Double_t cs, sn;
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
      //printf("%f\n",initPos[jj]);
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
  fVert.SetXYZ(initPos);
  fVert.SetDispersion(sigma);
  fVert.SetNContributors(ncombi);
}
//---------------------------------------------------------------------------
void AliVertexerTracks::VertexFitter(Bool_t useConstraint) 
{
//
// The optimal estimate of the vertex position is given by a "weighted 
// average of tracks positions".
// Original method: V. Karimaki, CMS Note 97/0051
//
  Double_t initPos[3];
  if(!fOnlyFitter) {
    fVert.GetXYZ(initPos);
  } else {
    initPos[0]=fNominalPos[0];
    initPos[1]=fNominalPos[1];
    initPos[2]=fNominalPos[2];
  }

  Int_t nTrks = (Int_t)fTrkArray.GetEntries();
  if(nTrks==1) useConstraint=kTRUE;
  if(fDebug) { 
    printf("--- VertexFitter(): start\n");
    printf(" Number of tracks in array: %d\n",nTrks);
    printf(" Minimum # tracks required in fit: %d\n",fMinTracks);
    printf(" Vertex position after finder: %f,%f,%f\n",initPos[0],initPos[1],initPos[2]);
    if(useConstraint) printf(" This vertex will be used in fit: (%f+-%f,%f+-%f)\n",fNominalPos[0],TMath::Sqrt(fNominalCov[0]),fNominalPos[1],TMath::Sqrt(fNominalCov[2])); 
  }

  // special treatment for few-tracks fits (e.g. secondary vertices)
  Bool_t uUi3by3 = kFALSE; if(nTrks<5 && !useConstraint) uUi3by3 = kTRUE;

  Int_t i,j,k,step=0;
  TMatrixD rv(3,1);
  TMatrixD vV(3,3);
  rv(0,0) = initPos[0];
  rv(1,0) = initPos[1];
  rv(2,0) = 0.;
  Double_t xlStart,alpha;
  Int_t nUsedTrks;
  Double_t chi2,chi2i,chi2b;
  AliESDtrack *t = 0;
  Int_t failed = 0;

  // initial vertex covariance matrix
  TMatrixD vVb(3,3);
  vVb(0,0) = fNominalCov[0];
  vVb(0,1) = fNominalCov[1];
  vVb(0,2) = 0.;
  vVb(1,0) = fNominalCov[1];
  vVb(1,1) = fNominalCov[2];
  vVb(1,2) = 0.;
  vVb(2,0) = 0.;
  vVb(2,1) = 0.;
  vVb(2,2) = fNominalCov[5];
  TMatrixD vVbInv(TMatrixD::kInverted,vVb);
  TMatrixD rb(3,1);
  rb(0,0) = fNominalPos[0];
  rb(1,0) = fNominalPos[1];
  rb(2,0) = fNominalPos[2];
  TMatrixD vVbInvrb(vVbInv,TMatrixD::kMult,rb);


  // 2 steps:
  // 1st - estimate of vtx using all tracks
  // 2nd - estimate of global chi2
  for(step=0; step<2; step++) {
    if(fDebug) printf(" step = %d\n",step);
    chi2 = 0.;
    nUsedTrks = 0;

    if(step==1) { initPos[0]=rv(0,0); initPos[0]=rv(1,0); }

    TMatrixD sumWiri(3,1);
    TMatrixD sumWi(3,3);
    for(i=0; i<3; i++) {
      sumWiri(i,0) = 0.;
      for(j=0; j<3; j++) sumWi(j,i) = 0.;
    }

    // mean vertex constraint
    if(useConstraint) {
      for(i=0;i<3;i++) {
	sumWiri(i,0) += vVbInvrb(i,0);
	for(k=0;k<3;k++) sumWi(i,k) += vVbInv(i,k);
      }
      // chi2
      TMatrixD deltar = rv; deltar -= rb;
      TMatrixD vVbInvdeltar(vVbInv,TMatrixD::kMult,deltar);
      chi2b = deltar(0,0)*vVbInvdeltar(0,0)+
              deltar(1,0)*vVbInvdeltar(1,0)+
	      deltar(2,0)*vVbInvdeltar(2,0);
      chi2 += chi2b;
    }

    // loop on tracks  
    for(k=0; k<nTrks; k++) {
      // get track from track array
      t = (AliESDtrack*)fTrkArray.At(k);
      alpha = t->GetAlpha();
      xlStart = initPos[0]*TMath::Cos(alpha)+initPos[1]*TMath::Sin(alpha);
      // to vtxSeed (from finder)
      t->AliExternalTrackParam::PropagateTo(xlStart,GetFieldkG());   
 
      // vector of track global coordinates
      TMatrixD ri(3,1);
      // covariance matrix of ri
      TMatrixD wWi(3,3);

      // get space point from track
      if(!TrackToPoint(t,ri,wWi,uUi3by3)) continue;
      TMatrixD wWiri(wWi,TMatrixD::kMult,ri); 

      // track chi2
      TMatrixD deltar = rv; deltar -= ri;
      TMatrixD wWideltar(wWi,TMatrixD::kMult,deltar);
      chi2i = deltar(0,0)*wWideltar(0,0)+
              deltar(1,0)*wWideltar(1,0)+
	      deltar(2,0)*wWideltar(2,0);

      // add to total chi2
      chi2 += chi2i;

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

    if(step==0) { 
      // inverted of weights matrix
      TMatrixD invsumWi(TMatrixD::kInverted,sumWi);
      vV = invsumWi;
      // position of primary vertex
      rv.Mult(vV,sumWiri);
    }
  } // end loop on the 2 steps

  if(failed) { 
    if(fDebug) printf("TooFewTracks\n");
    fCurrentVertex = new AliESDVertex(0.,0.,-1);
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
  
  // for correct chi2/ndf, count constraint as additional "track"
  if(fConstraint) nUsedTrks++;
  // store data in the vertex object
  fCurrentVertex = new AliESDVertex(position,covmatrix,chi2,nUsedTrks);

  if(fDebug) {
    printf(" Vertex after fit:\n");
    fCurrentVertex->PrintStatus();
    printf("--- VertexFitter(): finish\n");
  }

  return;
}
//----------------------------------------------------------------------------
AliESDVertex* AliVertexerTracks::VertexForSelectedTracks(TTree *trkTree,Bool_t optUseFitter,Bool_t optPropagate) 
{
//
// Return vertex from tracks in trkTree
//

  SetConstraintOff();

  // get tracks and propagate them to initial vertex position
  Int_t nTrksPrep = PrepareTracks(*trkTree,0);
  if(nTrksPrep <  TMath::Max(2,fMinTracks) ) {
    if(fDebug) printf("TooFewTracks\n");
    fCurrentVertex = new AliESDVertex(0.,0.,-1);
    return fCurrentVertex;
  }
 
  switch (fAlgo) {
    case 1: StrLinVertexFinderMinDist(1); break;
    case 2: StrLinVertexFinderMinDist(0); break;
    case 3: HelixVertexFinder();          break;
    case 4: VertexFinder(1);              break;
    case 5: VertexFinder(0);              break;
    default: printf("Wrong algorithm\n"); break;  
  }
  
  if(fDebug) printf(" Vertex finding completed\n");

  // vertex fitter
  if(optUseFitter){
    //SetVtxStart(&fVert);
    VertexFitter(fConstraint);
    if(fDebug) printf(" Vertex fit completed\n");
  }else{
    Double_t position[3]={fVert.GetXv(),fVert.GetYv(),fVert.GetZv()};
    Double_t covmatrix[6];
    fVert.GetCovMatrix(covmatrix);
    Double_t chi2=99999.;
    Int_t    nUsedTrks=fVert.GetNContributors();
    fCurrentVertex = new AliESDVertex(position,covmatrix,chi2,nUsedTrks);    
  }
  fCurrentVertex->SetDispersion(fVert.GetDispersion());


  // set indices of used tracks and propagate track to found vertex
  UShort_t *indices = 0;
  AliESDtrack *eta = 0;
  if(fCurrentVertex->GetNContributors()>0) {
    indices = new UShort_t[fCurrentVertex->GetNContributors()];
    for(Int_t jj=0;jj<(Int_t)fTrkArray.GetEntriesFast();jj++) {
      eta = (AliESDtrack*)fTrkArray.At(jj);
      indices[jj] = (UShort_t)eta->GetID();
      if(optPropagate&&optUseFitter){
	if(TMath::Sqrt(fCurrentVertex->GetXv()*fCurrentVertex->GetXv()+fCurrentVertex->GetYv()*fCurrentVertex->GetYv())<3.) {
	  eta->RelateToVertex(fCurrentVertex,GetFieldkG(),100.);
	  if(fDebug) printf("Track %d propagated to found vertex\n",jj);
	}else{
	  AliWarning("Found vertex outside beam pipe!");
	}
      }
    }
    fCurrentVertex->SetIndices(fCurrentVertex->GetNContributors(),indices);
  }
  delete [] indices;
  fTrkArray.Delete();
  
  return fCurrentVertex;
}
//----------------------------------------------------------------------------
AliESDVertex* AliVertexerTracks::VertexForSelectedTracks(TObjArray *trkArray,Bool_t optUseFitter, Bool_t optPropagate) 
{
//
// Return vertex from array of tracks
//

  // get tracks and propagate them to initial vertex position
  Int_t nTrks = trkArray->GetEntriesFast();
  if(nTrks < TMath::Max(2,fMinTracks) ) {
    if(fDebug) printf("TooFewTracks\n");
    fCurrentVertex = new AliESDVertex(0.,0.,-1);
    return fCurrentVertex;
  }
  TDirectory * olddir = gDirectory;
  TFile f("VertexerTracks.root","recreate");
  TTree *trkTree = new TTree("TreeT","tracks");
  AliESDtrack *esdTrack = 0;
  trkTree->Branch("tracks","AliESDtrack",&esdTrack);
  for(Int_t i=0; i<nTrks; i++){
    esdTrack = (AliESDtrack*)trkArray->At(i);
    trkTree->Fill();
  }
    
  AliESDVertex *vtx =  VertexForSelectedTracks(trkTree,optUseFitter,optPropagate);
  delete trkTree;
  f.Close();
  gSystem->Unlink("VertexerTracks.root");
  olddir->cd();
  return vtx;
}
//--------------------------------------------------------------------------
