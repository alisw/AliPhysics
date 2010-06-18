/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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

#include <TClonesArray.h>
#include <TTree.h>
#include "AliLog.h"
#include "AliESDVertex.h"
#include "AliITSgeomTGeo.h"
#include "AliITSRecPoint.h"
#include "AliITSReconstructor.h"
#include "AliITSVertexerCosmics.h"
#include "AliStrLine.h"

//------------------------------------------------------------------------
// This class implements a method to construct a "fake" primary
// vertex for cosmic events in which the muon crosses one of 5 inner
// ITS layers. A fake primary vertex is needed for the reconstruction,
// with e.g. AliITStrackerSA, of the two tracks produced by the muon 
// in the ITS.
//   We build pairs of clusters on a given layer and define the fake vertex as
// the mid-point of the straight line joining the two clusters.
//   We use the innermost layer that has at least two clusters.
//   We reject the background by requiring at least one cluster on the outer
// layer, closer than fMaxDistOnOuterLayer to the tracklet prolongation.
//   We can reject (potentially pathological) events with the muon track
// tangential to the layer by the requiring the radial position of
// the vertex to be smaller than fMaxVtxRadius.
//   Due to background clusters, more than one vertex per event can 
// be found. We consider the first found.
//   The errors on x,y,z of the vertex are calculated as errors on the mean
// of clusters coordinates. Non-diag elements of vertex cov. mat. are set to 0.
//   The number of contributors set in the AliESDVertex object is the
// number of the layer on which the tracklet was built; if this number is -1, 
// the procedure could not find a vertex position and by default 
// the vertex coordinates are set to (0,0,0) with large errors (100,100,100)
//
// Origin: A.Dainese, andrea.dainese@lnl.infn.it
//-------------------------------------------------------------------------

ClassImp(AliITSVertexerCosmics)

//-------------------------------------------------------------------------
AliITSVertexerCosmics::AliITSVertexerCosmics():AliITSVertexer(),
fMaxDistOnOuterLayer(0),
fMinDist2Vtxs(0)
{
  // Default constructor
  SetFirstLastModules(0,0,79);
  SetFirstLastModules(1,80,239);
  SetFirstLastModules(2,240,323);
  SetFirstLastModules(3,324,499);
  SetFirstLastModules(4,500,1247);
  SetFirstLastModules(5,1248,2197);
  /*
  SetMaxVtxRadius(0,3.5);
  SetMaxVtxRadius(1,6.5);
  SetMaxVtxRadius(2,14.5);
  SetMaxVtxRadius(3,23.5);
  SetMaxVtxRadius(4,37.5);
  SetMaxVtxRadius(5,42.5);
  */  
  SetMaxVtxRadius(0,5.5);
  SetMaxVtxRadius(1,8.5);
  SetMaxVtxRadius(2,18.5);
  SetMaxVtxRadius(3,28.5);
  SetMaxVtxRadius(4,39.5);
  SetMaxVtxRadius(5,48.5);
  
  SetMaxDistOnOuterLayer();
  SetMinDist2Vtxs();
}
//--------------------------------------------------------------------------
AliESDVertex* AliITSVertexerCosmics::FindVertexForCurrentEvent(TTree *itsClusterTree) 
{
  // Defines the AliESDVertex for the current event

  fCurrentVertex = 0;

  TClonesArray *recpoints=new TClonesArray("AliITSRecPoint",10000);
  itsClusterTree->SetBranchAddress("ITSRecPoints",&recpoints);

  Int_t lay,lad,det; 

  Double_t pos[3]={0.,0.,0.};
  Double_t err[3]={100.,100.,100.};
  Int_t ncontributors = -1;

  // Search for innermost layer with at least two clusters 
  // on two different modules
  Int_t ilayer=0,ilayer2=0;
  Int_t nHitModulesSPDinner=0;
  while(ilayer<AliITSgeomTGeo::GetNLayers()) {
    if(AliITSReconstructor::GetRecoParam()->GetLayersToSkip(ilayer)) {
      ilayer++;
      continue;
    }
    Int_t nHitModules=0;
    for(Int_t imodule=fFirst[ilayer]; imodule<=fLast[ilayer]; imodule++) {
      itsClusterTree->GetEvent(imodule);
      AliITSgeomTGeo::GetModuleId(imodule,lay,lad,det);
      lay -= 1;  // AliITSgeomTGeo gives layer from 1 to 6, we want 0 to 5
      if(lay!=ilayer) AliFatal("Layer mismatch!");
      if(recpoints->GetEntriesFast()>0) nHitModules++;
    }
    if(ilayer==0) nHitModulesSPDinner=nHitModules;
    if(nHitModules>=2) break;
    ilayer++;
  }

  ilayer2=ilayer+1;
  while(ilayer2<6) {
    if(!AliITSReconstructor::GetRecoParam()->GetLayersToSkip(ilayer2)) break;
    ilayer2++;
  }

  // try tracklet on SPD2 and point on SPD1
  if(ilayer==1 && !AliITSReconstructor::GetRecoParam()->GetLayersToSkip(0) &&
     nHitModulesSPDinner>0) { ilayer=0; ilayer2=1; }

  if(ilayer>4 || ilayer2>5) {
    AliWarning("Not enough clusters");
    delete recpoints;
    fCurrentVertex = new AliESDVertex(pos,err,"cosmics");
    fCurrentVertex->SetTitle("cosmics fake vertex (failed)");
    fCurrentVertex->SetNContributors(ncontributors);
    return fCurrentVertex;
  }


  const Int_t arrSize = 1000;
  Float_t xclInnLay[arrSize],yclInnLay[arrSize],zclInnLay[arrSize],modclInnLay[arrSize];
  Float_t e2xclInnLay[arrSize],e2yclInnLay[arrSize],e2zclInnLay[arrSize];
  Float_t e2xclOutLay[arrSize],e2yclOutLay[arrSize],e2zclOutLay[arrSize];
  Int_t nclInnLayStored=0;
  Float_t xclOutLay[arrSize],yclOutLay[arrSize],zclOutLay[arrSize],modclOutLay[arrSize];
  Int_t nclOutLayStored=0;
  Int_t nRecPoints,nRecPointsInnLay=0;

  Float_t gc[3],gcov[5];

  Float_t matchOutLayValue;
  Float_t distxyInnLay,distxyInnLayBest=0.;
  Double_t p1[3],p2[3],p3[3];

  Float_t xvtx,yvtx,zvtx,rvtx;
  Int_t i1InnLayBest=-1,i2InnLayBest=-1;


  // Collect clusters in the selected layer and the outer one
  for(Int_t imodule=fFirst[ilayer]; imodule<=fLast[ilayer2]; imodule++) {
    itsClusterTree->GetEvent(imodule);
    AliITSgeomTGeo::GetModuleId(imodule,lay,lad,det);
    lay -= 1; // AliITSgeomTGeo gives layer from 1 to 6, we want 0 to 5
    nRecPoints=recpoints->GetEntriesFast();
    if(imodule<=fLast[ilayer]) nRecPointsInnLay += nRecPoints;
    //printf("cosmics: module %d clusters %d\n",imodule,nRecPoints);
    for(Int_t irp=0; irp<nRecPoints; irp++) {
      AliITSRecPoint *rp=(AliITSRecPoint*)recpoints->UncheckedAt(irp);
      // Local coordinates of this recpoint
      rp->GetGlobalXYZ(gc);
      if(lay==ilayer) { // store InnLay clusters
	xclInnLay[nclInnLayStored]=gc[0];
	yclInnLay[nclInnLayStored]=gc[1];
	zclInnLay[nclInnLayStored]=gc[2];
	rp->GetGlobalCov(gcov);
	e2xclInnLay[nclInnLayStored]=gcov[0];
	e2yclInnLay[nclInnLayStored]=gcov[3];
	e2zclInnLay[nclInnLayStored]=gcov[5];
	modclInnLay[nclInnLayStored]=imodule;
	nclInnLayStored++;
      }
      if(lay==ilayer2) { // store OutLay clusters
	xclOutLay[nclOutLayStored]=gc[0];
	yclOutLay[nclOutLayStored]=gc[1];
	zclOutLay[nclOutLayStored]=gc[2];
	rp->GetGlobalCov(gcov);
	e2xclOutLay[nclOutLayStored]=gcov[0];
	e2yclOutLay[nclOutLayStored]=gcov[3];
	e2zclOutLay[nclOutLayStored]=gcov[5];
	modclOutLay[nclOutLayStored]=imodule;
	nclOutLayStored++;
      }
      if(nclInnLayStored>arrSize || nclOutLayStored>arrSize) {
	//AliFatal("More than arrSize clusters per layer");
	AliWarning("Too many clusters per layer");
	delete recpoints;
	fCurrentVertex = new AliESDVertex(pos,err,"cosmics");
	fCurrentVertex->SetTitle("cosmics fake vertex (failed)");
	fCurrentVertex->SetNContributors(ncontributors);
	return fCurrentVertex;
      }
    }// end clusters in a module
  }// end modules


  Int_t i1InnLay,i2InnLay,iOutLay;

  // build fake vertices
  //printf("Building tracklets on layer %d\n",ilayer);

  // InnLay - first cluster
  for(i1InnLay=0; i1InnLay<nclInnLayStored; i1InnLay++) { 
    p1[0]=xclInnLay[i1InnLay]; 
    p1[1]=yclInnLay[i1InnLay]; 
    p1[2]=zclInnLay[i1InnLay];
    // InnLay - second cluster
    for(i2InnLay=i1InnLay+1; i2InnLay<nclInnLayStored; i2InnLay++) { 
      if(modclInnLay[i1InnLay]==modclInnLay[i2InnLay]) continue;
      p2[0]=xclInnLay[i2InnLay]; 
      p2[1]=yclInnLay[i2InnLay]; 
      p2[2]=zclInnLay[i2InnLay];
      // look for point on OutLay
      AliStrLine innLayline(p1,p2,kTRUE);
      for(iOutLay=0; iOutLay<nclOutLayStored; iOutLay++) {
	p3[0]=xclOutLay[iOutLay]; 
	p3[1]=yclOutLay[iOutLay]; 
	p3[2]=zclOutLay[iOutLay];
	//printf("(%f,%f) (%f,%f)     (%f,%f) %f\n",p1[0],p1[1],p2[0],p2[1],p3[0],p3[1],innLayline.GetDistFromPoint(p3));
	matchOutLayValue=innLayline.GetDistFromPoint(p3);
	distxyInnLay = (p1[0]-p2[0])*(p1[0]-p2[0])+(p1[1]-p2[1])*(p1[1]-p2[1]);
	if(matchOutLayValue<fMaxDistOnOuterLayer &&
	   distxyInnLay>distxyInnLayBest) { 
	  //printf("found\n");
	  distxyInnLayBest=distxyInnLay;
	  i1InnLayBest=i1InnLay;
	  i2InnLayBest=i2InnLay;
	}
      }
    } // InnLay - second cluster
  } // InnLay - first cluster

  if(i1InnLayBest>-1 && i2InnLayBest>-1) { 
    xvtx = 0.5*(xclInnLay[i1InnLayBest]+xclInnLay[i2InnLayBest]);
    yvtx = 0.5*(yclInnLay[i1InnLayBest]+yclInnLay[i2InnLayBest]);
    zvtx = 0.5*(zclInnLay[i1InnLayBest]+zclInnLay[i2InnLayBest]);
    rvtx = TMath::Sqrt(xvtx*xvtx+yvtx*yvtx);
    if(rvtx<fMaxVtxRadius[ilayer]) {
      ncontributors = ilayer;
      pos[0] = xvtx;
      pos[1] = yvtx;
      pos[2] = zvtx;
      err[0]=TMath::Sqrt(0.25*(e2xclInnLay[i1InnLayBest]+e2xclInnLay[i2InnLayBest])); 
      err[1]=TMath::Sqrt(0.25*(e2yclInnLay[i1InnLayBest]+e2yclInnLay[i2InnLayBest])); 
      err[2]=TMath::Sqrt(0.25*(e2zclInnLay[i1InnLayBest]+e2zclInnLay[i2InnLayBest]));
    }

  } else { // give it a try exchanging InnLay and OutLay

    // OutLay - first cluster
    for(i1InnLay=0; i1InnLay<nclOutLayStored; i1InnLay++) { 
      p1[0]=xclOutLay[i1InnLay]; 
      p1[1]=yclOutLay[i1InnLay]; 
      p1[2]=zclOutLay[i1InnLay];
      // OutLay - second cluster
      for(i2InnLay=i1InnLay+1; i2InnLay<nclOutLayStored; i2InnLay++) { 
	if(modclOutLay[i1InnLay]==modclOutLay[i2InnLay]) continue;
	p2[0]=xclOutLay[i2InnLay]; 
	p2[1]=yclOutLay[i2InnLay]; 
	p2[2]=zclOutLay[i2InnLay];
	// look for point on InnLay
	AliStrLine outLayline(p1,p2,kTRUE);
	for(iOutLay=0; iOutLay<nclInnLayStored; iOutLay++) {
	  p3[0]=xclInnLay[iOutLay]; 
	  p3[1]=yclInnLay[iOutLay]; 
	  p3[2]=zclInnLay[iOutLay];
	  //printf(" %f\n",InnLayline.GetDistFromPoint(p3));
	  matchOutLayValue=outLayline.GetDistFromPoint(p3);
	  distxyInnLay = (p1[0]-p2[0])*(p1[0]-p2[0])+(p1[1]-p2[1])*(p1[1]-p2[1]);
	  if(matchOutLayValue<fMaxDistOnOuterLayer &&
	     distxyInnLay>distxyInnLayBest) { 
	    distxyInnLayBest=distxyInnLay;
	    i1InnLayBest=i1InnLay;
	    i2InnLayBest=i2InnLay;
	  }
	}
      } // OutLay - second cluster
    } // OutLay - first cluster

    if(i1InnLayBest>-1 && i2InnLayBest>-1) { 
      xvtx = 0.5*(xclOutLay[i1InnLayBest]+xclOutLay[i2InnLayBest]);
      yvtx = 0.5*(yclOutLay[i1InnLayBest]+yclOutLay[i2InnLayBest]);
      zvtx = 0.5*(zclOutLay[i1InnLayBest]+zclOutLay[i2InnLayBest]);
      rvtx = TMath::Sqrt(xvtx*xvtx+yvtx*yvtx);
      if(rvtx<fMaxVtxRadius[ilayer]) {
	ncontributors = ilayer2;
	pos[0] = xvtx;
	pos[1] = yvtx;
	pos[2] = zvtx;
	err[0]=TMath::Sqrt(0.25*(e2xclOutLay[i1InnLayBest]+e2xclOutLay[i2InnLayBest])); 
	err[1]=TMath::Sqrt(0.25*(e2yclOutLay[i1InnLayBest]+e2yclOutLay[i2InnLayBest])); 
	err[2]=TMath::Sqrt(0.25*(e2zclOutLay[i1InnLayBest]+e2zclOutLay[i2InnLayBest]));
      }
    }
  } // give it a try exchanging InnLay and OutLay

  fCurrentVertex = new AliESDVertex(pos,err,"cosmics");
  fCurrentVertex->SetTitle("cosmics fake vertex");
  fCurrentVertex->SetNContributors(ncontributors);
  //fCurrentVertex->Print();
  if(fComputeMultiplicity) FindMultiplicity(itsClusterTree);

  delete recpoints;

  return fCurrentVertex;
}  

//-------------------------------------------------------------------------
void AliITSVertexerCosmics::PrintStatus() const 
{
  // Print current status
  printf("=======================================================\n");
  printf(" fMaxDistOnOuterLayer: %f\n",fMaxDistOnOuterLayer);
  printf(" fMaxVtxRadius[0]:  %f\n",fMaxVtxRadius[0]);
  printf(" fMinDist2Vtxs:  %f\n",fMinDist2Vtxs);
  printf("=======================================================\n");
}
//-------------------------------------------------------------------------
