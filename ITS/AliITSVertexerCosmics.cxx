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
#include "AliLog.h"
#include "AliESDVertex.h"
#include "AliRunLoader.h"
#include "AliITSLoader.h"
#include "AliITSgeom.h"
#include "AliITSRecPoint.h"
#include "AliITSVertexerCosmics.h"
#include "AliStrLine.h"

//------------------------------------------------------------------------
// This class implements a method to construct a "fake" primary
// vertex for cosmic events in which the muon crosses the SPD inner
// layer (SPD1). A fake primary vertex is needed for the reconstruction,
// with e.g. AliITStrackerSA, of the two tracks produced by the muon 
// in the ITS.
//   We build pairs of clusters on SPD1 and define the fake vertex as
// the mid-point of the straight line joining the two clusters.
//   We reject the backgroung by requiring at least one clusters on SPD2
// closer than fMaxDistOnSPD2 to the tracklet prolongation.
//   We can reject (potentially pathological) events with the muon track
// tangential to the SPD1 layer by the requiring the radial position of
// the vertex to be smaller than fMaxVtxRadius.
//   Due to background clusters, more than one vertex per event can 
// be found. We consider the first found.
//   The number of contributors set in the AliESDVertex object is the
// number of vertices found in the event; if this number is <1, 
// the procedure could not find a vertex position and by default 
// the vertex coordinates are set to (0,0,0) with large errors (100,100,100)
// Number of contributors = 0  --> No SPD1 tracklets matching criteria 
// Number of contributors = -1 --> No SPD1 recpoints
//
// Origin: A.Dainese, andrea.dainese@lnl.infn.it
//-------------------------------------------------------------------------

ClassImp(AliITSVertexerCosmics)

//-------------------------------------------------------------------------
AliITSVertexerCosmics::AliITSVertexerCosmics():AliITSVertexer(),
fFirstSPD1(0),
fLastSPD1(0),
fFirstSPD2(0),
fLastSPD2(0),
fMaxDistOnSPD2(0),
fMaxVtxRadius(0),
fMinDist2Vtxs(0)
{
  // Default constructor
  SetSPD1Modules();
  SetSPD2Modules();
  SetMaxDistOnSPD2();
  SetMaxVtxRadius();
  SetMinDist2Vtxs();
}
//--------------------------------------------------------------------------
AliESDVertex* AliITSVertexerCosmics::FindVertexForCurrentEvent(Int_t evnumber) 
{
  // Defines the AliESDVertex for the current event

  fCurrentVertex = 0;
  AliRunLoader *rl =AliRunLoader::GetRunLoader();
  AliITSLoader* itsLoader = (AliITSLoader*)rl->GetLoader("ITSLoader");
  AliITSgeom* geom = itsLoader->GetITSgeom();
  itsLoader->LoadRecPoints();
  rl->GetEvent(evnumber);

  TTree *rpTree = itsLoader->TreeR();

  TClonesArray *recpoints=new TClonesArray("AliITSRecPoint",10000);
  rpTree->SetBranchAddress("ITSRecPoints",&recpoints);

  Double_t xclspd1[100],yclspd1[100],zclspd1[100],modclspd1[100];
  Int_t nclspd1stored=0;
  Double_t xclspd2[100],yclspd2[100],zclspd2[100],modclspd2[100];
  Int_t nclspd2stored=0;
  Int_t nrecpoints,nrecpointsSPD1=0;

  Double_t gc[3]={0.,0.,0.};
  Double_t lc[3]={0.,0.,0.};
  Int_t lay,lad,det; 

  Double_t x[100],y[100],z[100],p1[3],p2[3],p3[3];
  Int_t nvtxs;
  Bool_t good,matchtospd2;
  Double_t xvtx,yvtx,zvtx,rvtx;

  for(Int_t imodule=fFirstSPD1; imodule<fLastSPD2; imodule++) { // SPD
    rpTree->GetEvent(imodule);
    geom->GetModuleId(imodule,lay,lad,det);
    nrecpoints=recpoints->GetEntriesFast();
    if(imodule<fLastSPD1) nrecpointsSPD1 += nrecpoints;
    //printf("cosmics: module %d clusters %d\n",imodule,nrecpoints);
    for(Int_t irp=0; irp<nrecpoints; irp++) {
      AliITSRecPoint *rp=(AliITSRecPoint*)recpoints->UncheckedAt(irp);
      // Local coordinates of this recpoint
      lc[0]=rp->GetDetLocalX();
      lc[2]=rp->GetDetLocalZ();
      geom->LtoG(imodule,lc,gc); // global coordinates
      if(lay==1) { // store SPD1 clusters
	xclspd1[nclspd1stored]=gc[0];
	yclspd1[nclspd1stored]=gc[1];
	zclspd1[nclspd1stored]=gc[2];
	modclspd1[nclspd1stored]=imodule;
	nclspd1stored++;
      }
      if(lay==2) { // store SPD2 clusters
	xclspd2[nclspd2stored]=gc[0];
	yclspd2[nclspd2stored]=gc[1];
	zclspd2[nclspd2stored]=gc[2];
	modclspd2[nclspd2stored]=imodule;
	nclspd2stored++;
      }
      if(nclspd1stored>100 || nclspd2stored>100) 
	AliFatal("More than 100 clusters per layer in SPD");
    }// end clusters in a module
  }// end SPD modules for a given event

  // build fake vertices
  nvtxs=0;
  // SPD1 - first cluster
  for(Int_t i1spd1=0; i1spd1<nclspd1stored; i1spd1++) { 
    p1[0]=xclspd1[i1spd1]; p1[1]=yclspd1[i1spd1]; p1[2]=zclspd1[i1spd1];
    // SPD1 - second cluster
    for(Int_t i2spd1=i1spd1+1; i2spd1<nclspd1stored; i2spd1++) { 
      if(modclspd1[i1spd1]==modclspd1[i2spd1]) continue;
      p2[0]=xclspd1[i2spd1]; p2[1]=yclspd1[i2spd1]; p2[2]=zclspd1[i2spd1];
      // look for point on SPD2
      AliStrLine spd1line(p1,p2,kTRUE);
      matchtospd2 = kFALSE;
      for(Int_t ispd2=0; ispd2<nclspd2stored; ispd2++) {
	p3[0]=xclspd2[ispd2]; p3[1]=yclspd2[ispd2]; p3[2]=zclspd2[ispd2];
	//printf(" %f\n",spd1line.GetDistFromPoint(p3));
	if(spd1line.GetDistFromPoint(p3)<fMaxDistOnSPD2) 
	  { matchtospd2 = kTRUE; break; }
      }
      if(!matchtospd2) continue;
      xvtx = 0.5*(xclspd1[i1spd1]+xclspd1[i2spd1]);
      yvtx = 0.5*(yclspd1[i1spd1]+yclspd1[i2spd1]);
      zvtx = 0.5*(zclspd1[i1spd1]+zclspd1[i2spd1]);
      rvtx = TMath::Sqrt(xvtx*xvtx+yvtx*yvtx);
      if(rvtx>fMaxVtxRadius) continue;
      good = kTRUE;
      for(Int_t iv=0; iv<nvtxs; iv++) {
	if(TMath::Sqrt((xvtx- x[iv])*(xvtx- x[iv])+
		       (yvtx- y[iv])*(yvtx- y[iv])+
		       (zvtx- z[iv])*(zvtx- z[iv])) < fMinDist2Vtxs) 
	  good = kFALSE;
      }
      if(good) {
	x[nvtxs]=xvtx;
	y[nvtxs]=yvtx;
	z[nvtxs]=zvtx;
	nvtxs++;
      }
    } // SPD1 - second cluster
  } // SPD1 - first cluster


  Double_t pos[3]={0.,0.,0.};
  Double_t err[3]={100.,100.,100.};
  if(nvtxs) { 
    pos[0]=x[0]; 
    pos[1]=y[0]; 
    pos[2]=z[0];
    err[0]=0.1; 
    err[1]=0.1; 
    err[2]=0.1;
  }
  if(!nrecpointsSPD1) nvtxs=-1;
  fCurrentVertex = new AliESDVertex(pos,err,"cosmics");
  fCurrentVertex->SetNContributors(nvtxs);
  fCurrentVertex->SetTitle("cosmics fake vertex");

  //if(nvtxs>0) fCurrentVertex->Print();

  delete recpoints;
  itsLoader->UnloadRecPoints();

  return fCurrentVertex;
}  
//-------------------------------------------------------------------------
void AliITSVertexerCosmics::FindVertices()
{
  // computes the vertices of the events in the range FirstEvent - LastEvent
  AliRunLoader *rl = AliRunLoader::GetRunLoader();
  AliITSLoader* itsLoader =  (AliITSLoader*) rl->GetLoader("ITSLoader");
  itsLoader->ReloadRecPoints();
  for(Int_t i=fFirstEvent;i<=fLastEvent;i++){
    //  printf("Processing event %d\n",i);
    rl->GetEvent(i);
    FindVertexForCurrentEvent(i);
    if(fCurrentVertex){
      WriteCurrentVertex();
    }
  }
}
//-------------------------------------------------------------------------
void AliITSVertexerCosmics::PrintStatus() const 
{
  // Print current status
  printf("=======================================================\n");
  printf(" fMaxDistOnSPD2: %f\n",fMaxDistOnSPD2);
  printf(" fMaxVtxRadius:  %f\n",fMaxVtxRadius);
  printf(" fMinDist2Vtxs:  %f\n",fMinDist2Vtxs);
  printf(" First layer first and last modules: %d, %d\n",fFirstSPD1,fLastSPD1);
  printf(" Second layer first and last modules: %d, %d\n",fFirstSPD2,fLastSPD2);
  printf("=======================================================\n");
}
//-------------------------------------------------------------------------
