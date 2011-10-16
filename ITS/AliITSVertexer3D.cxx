/**************************************************************************
 * Copyright(c) 2006-2008, ALICE Experiment at CERN, All rights reserved. *
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
#include <TTree.h>
#include "AliRunLoader.h"
#include "AliESDVertex.h"
#include "AliLog.h"
#include "AliStrLine.h"
#include "AliTracker.h"
#include "AliITSDetTypeRec.h"
#include "AliITSRecPoint.h"
#include "AliITSRecPointContainer.h"
#include "AliITSgeomTGeo.h"
#include "AliVertexerTracks.h"
#include "AliITSVertexer3D.h"
#include "AliITSVertexerZ.h"
#include "AliITSSortTrkl.h"
/////////////////////////////////////////////////////////////////
// this class implements a method to determine
// the 3 coordinates of the primary vertex
// optimized for 
// p-p collisions
////////////////////////////////////////////////////////////////

const Int_t    AliITSVertexer3D::fgkMaxNumOfClDefault = 300;
const Int_t    AliITSVertexer3D::fgkMaxNumOfClRebinDefault = 500;
const Int_t    AliITSVertexer3D::fgkMaxNumOfClDownscaleDefault = 1000;
const Float_t  AliITSVertexer3D::fgk3DBinSizeDefault = 0.1;

ClassImp(AliITSVertexer3D)

/* $Id$ */

//______________________________________________________________________
AliITSVertexer3D::AliITSVertexer3D():
  AliITSVertexer(),
  fLines("AliStrLine",1000),
  fVert3D(),
  fCoarseDiffPhiCut(0.),
  fFineDiffPhiCut(0.),
  fCutOnPairs(0.),
  fCoarseMaxRCut(0.),
  fMaxRCut(0.),
  fMaxRCut2(0.),
  fZCutDiamond(0.),
  fMaxZCut(0.),
  fDCAcut(0.),
  fDiffPhiMax(0.),
  fMeanPSelTrk(0.),
  fMeanPtSelTrk(0.),
  fUsedCluster(kMaxCluPerMod*kNSPDMod),
  fZHisto(0),
  fDCAforPileup(0.),
  fDiffPhiforPileup(0.),
  fBinSizeR(0.),
  fBinSizeZ(0.),
  fPileupAlgo(0),
  fMaxNumOfCl(fgkMaxNumOfClDefault),
  fMaxNumOfClForRebin(fgkMaxNumOfClRebinDefault),
  fMaxNumOfClForDownScale(fgkMaxNumOfClDownscaleDefault),
  fNRecPLay1(0),
  fNRecPLay2(0),
  f3DBinSize(fgk3DBinSizeDefault),
  fDoDownScale(kFALSE),
  fGenerForDownScale(0),
  f3DPeak(),
  fHighMultAlgo(1),
  fSwitchAlgorithm(kFALSE)
{
  // Default constructor
  SetCoarseDiffPhiCut();
  SetFineDiffPhiCut();
  SetCutOnPairs();
  SetCoarseMaxRCut();
  SetMaxRCut();
  SetMaxRCutAlgo2();
  SetZCutDiamond();
  SetMaxZCut();
  SetDCACut();
  SetDiffPhiMax();
  SetMeanPSelTracks();
  SetMeanPtSelTracks();
  SetMinDCAforPileup();
  SetDeltaPhiforPileup();
  SetPileupAlgo();
  SetBinSizeR();
  SetBinSizeZ();
  Double_t binsize=0.02; // default 200 micron
  Int_t nbins=static_cast<Int_t>(1+2*fZCutDiamond/binsize);
  fZHisto=new TH1F("hz","",nbins,-fZCutDiamond,-fZCutDiamond+binsize*nbins);
  fGenerForDownScale=new TRandom3(987654321);
}

//______________________________________________________________________
AliITSVertexer3D::~AliITSVertexer3D() {
  // Destructor
  fLines.Clear("C");
  if(fZHisto) delete fZHisto;
  if(fGenerForDownScale) delete fGenerForDownScale;
}

//______________________________________________________________________
void AliITSVertexer3D::ResetVert3D(){
  //
  ResetVertex();
  fVert3D.SetXv(0.);
  fVert3D.SetYv(0.);
  fVert3D.SetZv(0.);
  fVert3D.SetDispersion(0.);
  fVert3D.SetNContributors(0);
  fUsedCluster.ResetAllBits(0);
}
//______________________________________________________________________
AliESDVertex* AliITSVertexer3D::FindVertexForCurrentEvent(TTree *itsClusterTree){
  // Defines the AliESDVertex for the current event
  ResetVert3D();
  AliDebug(1,"FindVertexForCurrentEvent - 3D - PROCESSING NEXT EVENT");
  fLines.Clear("C");
  fCurrentVertex = NULL;

  Int_t nolines = FindTracklets(itsClusterTree,0);
  Int_t rc;
  if(nolines>=2){
    if(fSwitchAlgorithm) {
      rc = Prepare3DVertexPbPb();
      FindVertex3D(itsClusterTree);
    } else {
      rc=Prepare3DVertex(0);
      if(fVert3D.GetNContributors()>0){
	fLines.Clear("C");
	nolines = FindTracklets(itsClusterTree,1);
	if(nolines>=2){
	  rc=Prepare3DVertex(1);
	  if(fPileupAlgo == 2 && rc == 0) FindVertex3DIterative();
	  else if(fPileupAlgo!=2 && rc == 0) FindVertex3D(itsClusterTree);
	  if(rc!=0) fVert3D.SetNContributors(0); // exclude this vertex      
	}
      }
    }
  }

  if(!fCurrentVertex){
    AliITSVertexerZ vertz(GetNominalPos()[0],GetNominalPos()[1]);
    vertz.SetDetTypeRec(GetDetTypeRec());
    AliDebug(1,"Call Vertexer Z\n");
    vertz.SetLowLimit(-fZCutDiamond);
    vertz.SetHighLimit(fZCutDiamond);
    AliESDVertex* vtxz = vertz.FindVertexForCurrentEvent(itsClusterTree);
    if(vtxz){
      Double_t position[3]={GetNominalPos()[0],GetNominalPos()[1],vtxz->GetZv()};
      Double_t covmatrix[6];
      vtxz->GetCovMatrix(covmatrix);
      Double_t chi2=99999.;
      Int_t    nContr=vtxz->GetNContributors();
      fCurrentVertex = new AliESDVertex(position,covmatrix,chi2,nContr);    
      fCurrentVertex->SetDispersion(vtxz->GetDispersion());
      fCurrentVertex->SetTitle("vertexer: Z");
      fCurrentVertex->SetName("SPDVertexZ");
      delete vtxz;
    }

  }  
  if(fComputeMultiplicity) FindMultiplicity(itsClusterTree);
  return fCurrentVertex;
}  

//______________________________________________________________________
void AliITSVertexer3D::FindVertex3D(TTree *itsClusterTree){
 
  Double_t vRadius=TMath::Sqrt(fVert3D.GetXv()*fVert3D.GetXv()+fVert3D.GetYv()*fVert3D.GetYv());
  if(vRadius<GetPipeRadius() && fVert3D.GetNContributors()>0){
    Double_t position[3]={fVert3D.GetXv(),fVert3D.GetYv(),fVert3D.GetZv()};
    Double_t covmatrix[6];
    fVert3D.GetCovMatrix(covmatrix);
    Double_t chi2=99999.;
    Int_t    nContr=fVert3D.GetNContributors();
    fCurrentVertex = new AliESDVertex(position,covmatrix,chi2,nContr);    
    fCurrentVertex->SetTitle("vertexer: 3D");
    fCurrentVertex->SetName("SPDVertex3D");
    fCurrentVertex->SetDispersion(fVert3D.GetDispersion());
    fNoVertices=1;
    
    switch(fPileupAlgo){
    case 0: PileupFromZ(); break;
    case 1: FindOther3DVertices(itsClusterTree); break;
    case 3: break; // no pileup algo  
    default: AliError("Wrong pileup algorithm"); break;
    }
    if(fNoVertices==1){
      fVertArray = new AliESDVertex[1];
      fVertArray[0]=(*fCurrentVertex);	  
    }
  }
}

//______________________________________________________________________
void AliITSVertexer3D::FindVertex3DIterative(){
  //

  Int_t nLines=fLines.GetEntriesFast();
  Int_t maxPoints=nLines*(nLines-1)/2;
  Double_t* xP=new Double_t[maxPoints];
  Double_t* yP=new Double_t[maxPoints];
  Double_t* zP=new Double_t[maxPoints];
  Int_t* index1=new Int_t[maxPoints];
  Int_t* index2=new Int_t[maxPoints];
  Double_t xbeam=fVert3D.GetXv();
  Double_t ybeam=fVert3D.GetYv();

  Int_t iPoint=0;
  for(Int_t ilin1=0; ilin1<nLines; ilin1++){
    AliStrLine *l1 = (AliStrLine*)fLines.At(ilin1);
    for(Int_t ilin2=ilin1+1; ilin2<nLines; ilin2++){
      AliStrLine *l2 = (AliStrLine*)fLines.At(ilin2);
      Double_t dca=l1->GetDCA(l2);
      if(dca > fDCAcut || dca<0.00001) continue;
      Double_t point[3];
      Int_t retc = l1->Cross(l2,point);
      if(retc<0)continue;
      Double_t rad=TMath::Sqrt(point[0]*point[0]+point[1]*point[1]);
      if(rad>fCoarseMaxRCut)continue;
      Double_t distFromBeam=TMath::Sqrt((point[0]-xbeam)*(point[0]-xbeam)+(point[1]-ybeam)*(point[1]-ybeam));
      if(distFromBeam>fMaxRCut2) continue;
      xP[iPoint]=point[0];
      yP[iPoint]=point[1];
      zP[iPoint]=point[2];
      index1[iPoint]=ilin1;
      index2[iPoint]=ilin2;
      iPoint++;
    }
  }
  Int_t npoints=iPoint++;
  Int_t index=0;
  Short_t* mask=new Short_t[npoints];
  for(Int_t ip=0;ip<npoints;ip++) mask[ip]=-1;
 
  for(Int_t ip1=0;ip1<npoints;ip1++){
    if(mask[ip1]==-1) mask[ip1]=index++;
    for(Int_t ip2=ip1+1; ip2<npoints; ip2++){
      if(mask[ip2]==mask[ip1] && mask[ip2]!=-1) continue;
      Double_t dist2=(xP[ip1]-xP[ip2])*(xP[ip1]-xP[ip2]);
      dist2+=(yP[ip1]-yP[ip2])*(yP[ip1]-yP[ip2]);
      dist2+=(zP[ip1]-zP[ip2])*(zP[ip1]-zP[ip2]);
      if(dist2<fCutOnPairs*fCutOnPairs){ 
	if(mask[ip2]==-1) mask[ip2]=mask[ip1];
	else{
	  for(Int_t ip=0; ip<npoints;ip++){
	    if(mask[ip]==mask[ip2]) mask[ip]=mask[ip1];
	  }
	}
      }
    }
  }


  // Count multiplicity of trackelts in clusters
  UInt_t* isIndUsed=new UInt_t[index+1];
  for(Int_t ind=0;ind<index+1;ind++) isIndUsed[ind]=0;
  for(Int_t ip=0; ip<npoints;ip++){
    Int_t ind=mask[ip];
    isIndUsed[ind]++;
  }

  // Count clusters/vertices and sort according to multiplicity
  Int_t nClusters=0;
  Int_t* sortedIndex=new Int_t[index+1];
  for(Int_t ind1=0;ind1<index+1;ind1++){
    if(isIndUsed[ind1]<=1) isIndUsed[ind1]=0;
    else nClusters++;
    UInt_t cap=9999999;
    if(ind1>0) cap=isIndUsed[sortedIndex[ind1-1]];
    UInt_t bigger=0;
    Int_t biggerindex=-1;
    for(Int_t ind2=0;ind2<index+1;ind2++){
      Bool_t use=kTRUE;
      for(Int_t ind3=0; ind3<ind1; ind3++)
	if(ind2==sortedIndex[ind3]) use=kFALSE;
      if(use && isIndUsed[ind2]>bigger && isIndUsed[ind2]<=cap){
	bigger=isIndUsed[ind2];
	biggerindex=ind2;
      }
    }
    sortedIndex[ind1]=biggerindex;    
  }
  AliDebug(3,Form("Number of clusters before merging = %d\n",nClusters));

  // Assign lines to clusters/vertices and merge clusters which share 1 line
  Int_t nClustersAfterMerge=nClusters;
  Int_t* belongsTo=new Int_t[nLines];
  for(Int_t ilin=0; ilin<nLines; ilin++) belongsTo[ilin]=-1;
  for(Int_t iclu=0;iclu<nClusters;iclu++){
    Int_t actualCluIndex=iclu;
    for(Int_t ip=0; ip<npoints;ip++){
      if(mask[ip]==sortedIndex[iclu]){
	Int_t ind1=index1[ip];
	if(belongsTo[ind1]==-1) belongsTo[ind1]=actualCluIndex;
	else if(belongsTo[ind1]<actualCluIndex){
	  Int_t newCluIndex=belongsTo[ind1];
	  for(Int_t ilin=0; ilin<nLines; ilin++){
	    if(belongsTo[ilin]==actualCluIndex) belongsTo[ilin]=newCluIndex;
	  }
	  AliDebug(10,Form("Merged cluster %d with %d\n",actualCluIndex,newCluIndex));
	  actualCluIndex=newCluIndex;
	  nClustersAfterMerge--;
	}
	Int_t ind2=index2[ip];      
	if(belongsTo[ind2]==-1) belongsTo[ind2]=actualCluIndex;
	else if(belongsTo[ind2]<actualCluIndex){
	  Int_t newCluIndex=belongsTo[ind2];
	  for(Int_t ilin=0; ilin<nLines; ilin++){
	    if(belongsTo[ilin]==actualCluIndex) belongsTo[ilin]=newCluIndex;
	  }
	  AliDebug(10,Form("Merged cluster %d with %d\n",actualCluIndex,newCluIndex));
	  actualCluIndex=newCluIndex;
	  nClustersAfterMerge--;
	}
      }
    }
  }
  AliDebug(3,Form("Number of clusters after merging = %d\n",nClustersAfterMerge));
  
  // Count lines associated to each cluster/vertex
  UInt_t *cluSize=new UInt_t[nClusters];
  for(Int_t iclu=0;iclu<nClusters;iclu++){ 
    cluSize[iclu]=0;
    for(Int_t ilin=0; ilin<nLines; ilin++){
      if(belongsTo[ilin]==iclu) cluSize[iclu]++;
    }
  }

  // Count good vertices (>1 associated tracklet)
  UInt_t nGoodVert=0;
  for(Int_t iclu=0;iclu<nClusters;iclu++){ 
    AliDebug(3,Form("Vertex %d Size=%d\n",iclu,cluSize[iclu]));
    if(cluSize[iclu]>1) nGoodVert++;
  }
    
  AliDebug(1,Form("Number of good vertices = %d\n",nGoodVert));
  // Calculate vertex coordinates for each cluster
  if(nGoodVert>0){
    fVertArray = new AliESDVertex[nGoodVert];
    Int_t iVert=0;
    for(Int_t iclu=0;iclu<nClusters;iclu++){
      Int_t size=cluSize[iclu];
      if(size>1){
	AliStrLine **arrlin = new AliStrLine*[size];
	Int_t nFilled=0;
	for(Int_t ilin=0; ilin<nLines; ilin++){
	  if(belongsTo[ilin]==iclu){
	    arrlin[nFilled++] = dynamic_cast<AliStrLine*>(fLines[ilin]);
	  }
	}      
	AliDebug(3,Form("Vertex %d  N associated tracklets = %d out of %d\n",iVert,size,nFilled));

	fVertArray[iVert]=AliVertexerTracks::TrackletVertexFinder(arrlin,nFilled);
	Double_t peak[3];
	fVertArray[iVert].GetXYZ(peak);
	AliStrLine **arrlin2 = new AliStrLine*[size];
	Int_t nFilled2=0;	
	for(Int_t i=0; i<nFilled;i++){
	  AliStrLine *l1 = arrlin[i];	  
	  if(l1->GetDistFromPoint(peak)< fDCAcut)
	    arrlin2[nFilled2++] = dynamic_cast<AliStrLine*>(l1);
	}
	if(nFilled2>1){
	  AliDebug(3,Form("Vertex %d  recalculated with %d tracklets\n",iVert,nFilled2));
	  fVertArray[iVert]=AliVertexerTracks::TrackletVertexFinder(arrlin2,nFilled2);
	}
 	delete [] arrlin;
 	delete [] arrlin2;
	++iVert;
      }
    }
    
    if(nGoodVert > 1){
      fIsPileup = kTRUE;
      fNTrpuv = fVertArray[1].GetNContributors();
      fZpuv = fVertArray[1].GetZv();
    }
    
    Double_t vRadius=TMath::Sqrt(fVertArray[0].GetXv()*fVertArray[0].GetXv()+fVertArray[0].GetYv()*fVertArray[0].GetYv());
    if(vRadius<GetPipeRadius() && fVertArray[0].GetNContributors()>0){
      Double_t position[3]={fVertArray[0].GetXv(),fVertArray[0].GetYv(),fVertArray[0].GetZv()};
      Double_t covmatrix[6];
      fVertArray[0].GetCovMatrix(covmatrix);
      Double_t chi2=99999.;
      Int_t    nContr=fVertArray[0].GetNContributors();
      fCurrentVertex = new AliESDVertex(position,covmatrix,chi2,nContr);    
      fCurrentVertex->SetTitle("vertexer: 3D");
      fCurrentVertex->SetName("SPDVertex3D");
      fCurrentVertex->SetDispersion(fVertArray[0].GetDispersion());  
    }
  }

  delete [] index1;
  delete [] index2;
  delete [] mask;
  delete [] isIndUsed;
  delete [] sortedIndex;
  delete [] belongsTo;
  delete [] cluSize;
  delete [] xP;
  delete [] yP;
  delete [] zP;
}
//______________________________________________________________________
void AliITSVertexer3D::FindVertex3DIterativeMM(){
  // Defines the AliESDVertex for the current event
  Int_t numsor=fLines.GetEntriesFast()*(fLines.GetEntriesFast()-1)/2;
  //cout<<"AliITSVertexer3D::FindVertexForCurentEvent: Number of tracklets selected for vertexing "<<fLines.GetEntriesFast()<<"; Number of pairs: "<<numsor<<endl;
  AliITSSortTrkl srt(fLines,numsor,fCutOnPairs,fCoarseMaxRCut); 
  srt.FindClusters();	
  AliInfo(Form("Number of vertices: %d",srt.GetNumberOfClusters()));
      
  fNoVertices = srt.GetNumberOfClusters();
  //printf("fNoVertices = %d \n",fNoVertices);
  if(fNoVertices>0){
    fVertArray = new AliESDVertex[fNoVertices];
    for(Int_t kk=0; kk<srt.GetNumberOfClusters(); kk++){
      Int_t size = 0;
      Int_t *labels = srt.GetTrackletsLab(kk,size);
      /*
	Int_t *pairs = srt.GetClusters(kk);
	Int_t nopai = srt.GetSizeOfCluster(kk);
	cout<<"***** Vertex number "<<kk<<".  Pairs: \n";
	for(Int_t jj=0;jj<nopai;jj++){
	cout<<pairs[jj]<<" - ";
	if(jj>0 & jj%8==0)cout<<endl;
	}
	cout<<endl;
	cout<<"***** Vertex number "<<kk<<".  Labels: \n";
      */
      AliStrLine **tclo = new AliStrLine* [size];
      for(Int_t jj=0;jj<size;jj++){
	//	    cout<<labels[jj]<<" - ";
	//	    if(jj>0 & jj%8==0)cout<<endl;
	tclo[jj] = dynamic_cast<AliStrLine*>(fLines[labels[jj]]);
      }
      //	  cout<<endl;
      delete []labels;
      fVertArray[kk]=AliVertexerTracks::TrackletVertexFinder(tclo,size);
      delete [] tclo;
      //	  fVertArray[kk].PrintStatus();
      if(kk == 1){
	// at least one second vertex is present
	fIsPileup = kTRUE;
	fNTrpuv = fVertArray[kk].GetNContributors();
	fZpuv = fVertArray[kk].GetZv();
      }
    }
    Double_t vRadius=TMath::Sqrt(fVertArray[0].GetXv()*fVertArray[0].GetXv()+fVertArray[0].GetYv()*fVertArray[0].GetYv());
    if(vRadius<GetPipeRadius() && fVertArray[0].GetNContributors()>0){
      Double_t position[3]={fVertArray[0].GetXv(),fVertArray[0].GetYv(),fVertArray[0].GetZv()};
      Double_t covmatrix[6];
      fVertArray[0].GetCovMatrix(covmatrix);
      Double_t chi2=99999.;
      Int_t    nContr=fVertArray[0].GetNContributors();
      fCurrentVertex = new AliESDVertex(position,covmatrix,chi2,nContr);    
      fCurrentVertex->SetTitle("vertexer: 3D");
      fCurrentVertex->SetName("SPDVertex3D");
      fCurrentVertex->SetDispersion(fVertArray[0].GetDispersion());  
    }
  }

}  

//______________________________________________________________________
Bool_t AliITSVertexer3D::DistBetweenVertices(AliESDVertex &a, AliESDVertex &b, Double_t test, Double_t &dist){
  // method to compare the distance between vertices a and b with "test"
  //it returns kTRUE is the distance is less or equal to test
  dist = (a.GetX()-b.GetX()) * (a.GetX()-b.GetX());
  dist +=  (a.GetY()-b.GetY()) * (a.GetY()-b.GetY());
  dist +=  (a.GetZ()-b.GetZ()) * (a.GetZ()-b.GetZ());
  dist = TMath::Sqrt(dist);
  if(dist <= test)return kTRUE;
  return kFALSE;
}


//______________________________________________________________________
Int_t AliITSVertexer3D::FindTracklets(TTree *itsClusterTree, Int_t optCuts){
  // All the possible combinations between recpoints on layer 1and 2 are
  // considered. Straight lines (=tracklets)are formed. 
  // The tracklets are processed in Prepare3DVertex

  TClonesArray *itsRec  = 0;
  if(optCuts==0) fZHisto->Reset();
  // gc1 are local and global coordinates for layer 1
  Float_t gc1f[3]={0.,0.,0.};
  Double_t gc1[3]={0.,0.,0.};
  // gc2 are local and global coordinates for layer 2
  Float_t gc2f[3]={0.,0.,0.};
  Double_t gc2[3]={0.,0.,0.};
  AliITSRecPointContainer* rpcont=AliITSRecPointContainer::Instance();
  rpcont->FetchClusters(0,itsClusterTree);
  if(!rpcont->IsSPDActive()){
    AliWarning("No SPD rec points found, 3D vertex not calculated");
    return -1;
  }

  // Set values for cuts
  Double_t xbeam=GetNominalPos()[0]; 
  Double_t ybeam=GetNominalPos()[1];
  Double_t zvert=0.;
  Double_t deltaPhi=fCoarseDiffPhiCut;
  Double_t deltaR=fCoarseMaxRCut;
  Double_t dZmax=fZCutDiamond;
  if(optCuts==1){
    xbeam=fVert3D.GetXv();
    ybeam=fVert3D.GetYv();
    zvert=fVert3D.GetZv();
    deltaPhi = fDiffPhiMax; 
    deltaR=fMaxRCut;
    dZmax=fMaxZCut;
    if(fPileupAlgo == 2){
      dZmax=fZCutDiamond;
      deltaR=fMaxRCut2;
    }
  } else if(optCuts==2){
    xbeam=fVert3D.GetXv();
    ybeam=fVert3D.GetYv();
    deltaPhi = fDiffPhiforPileup;
    deltaR=fMaxRCut;
  }

  fNRecPLay1=rpcont->GetNClustersInLayerFast(1);
  fNRecPLay2=rpcont->GetNClustersInLayerFast(2);
  if(fNRecPLay1 == 0 || fNRecPLay2 == 0){
    AliDebug(1,Form("No RecPoints in at least one SPD layer (%d %d)",fNRecPLay1,fNRecPLay2));
    return -1;
  }
  AliDebug(1,Form("RecPoints on Layer 1,2 = %d, %d\n",fNRecPLay1,fNRecPLay2));
  fDoDownScale=kFALSE;
  fSwitchAlgorithm=kFALSE;

  Float_t factDownScal=1.;
  Int_t origLaddersOnLayer2=fLadOnLay2;

  switch(fHighMultAlgo){
  case 0: 
    if(fNRecPLay1>fMaxNumOfClForDownScale || fNRecPLay2>fMaxNumOfClForDownScale){
      if(optCuts==2) return -1; // do not try to search for pileup
      SetLaddersOnLayer2(2);
      fDoDownScale=kTRUE;
      factDownScal=(Float_t)fMaxNumOfClForDownScale*(Float_t)fMaxNumOfClForDownScale/(Float_t)fNRecPLay1/(Float_t)fNRecPLay2;
      if(optCuts==1){
	factDownScal*=(fCoarseDiffPhiCut/fDiffPhiMax)*10;
	if(factDownScal>1.){
	  fDoDownScale=kFALSE;
	  SetLaddersOnLayer2(origLaddersOnLayer2);
	}
      }
      if(fDoDownScale)AliDebug(1,Form("Too many recpoints on SPD(%d %d ), downscale by %f",fNRecPLay1,fNRecPLay2,factDownScal));
    }
    break;
  case 1: 
    if(fNRecPLay1>fMaxNumOfCl || fNRecPLay2>fMaxNumOfCl) {
      if(optCuts==2) return -1; // do not try to search for pileup
      fSwitchAlgorithm=kTRUE;
    }
    break;
  default: break; // no pileup algo  
  }

  if(!fDoDownScale && !fSwitchAlgorithm){
    if(fNRecPLay1>fMaxNumOfClForRebin || fNRecPLay2>fMaxNumOfClForRebin){
      SetLaddersOnLayer2(2);
    }
  }

  Double_t a[3]={xbeam,ybeam,0.}; 
  Double_t b[3]={xbeam,ybeam,10.};
  AliStrLine zeta(a,b,kTRUE);
  static Double_t bField=TMath::Abs(AliTracker::GetBz()/10.); //T
  SetMeanPPtSelTracks(bField);

  Int_t nolines = 0;
  // Loop on modules of layer 1
  Int_t firstL1 = TMath::Max(0,AliITSgeomTGeo::GetModuleIndex(1,1,1));
  Int_t lastL1 = AliITSgeomTGeo::GetModuleIndex(2,1,1)-1;
  for(Int_t modul1= firstL1; modul1<=lastL1;modul1++){   // Loop on modules of layer 1
    if(!fUseModule[modul1]) continue;
    
    UShort_t ladder=modul1/4+1; // ladders are numbered starting from 1
    TClonesArray *prpl1=rpcont->UncheckedGetClusters(modul1);
    Int_t nrecp1 = prpl1->GetEntries();
    for(Int_t j=0;j<nrecp1;j++){
      if(j>kMaxCluPerMod) continue;
      UShort_t idClu1=modul1*kMaxCluPerMod+j;
      if(fUsedCluster.TestBitNumber(idClu1)) continue;
      if(fDoDownScale && !fSwitchAlgorithm){
	if(fGenerForDownScale->Rndm()>factDownScal) continue;
      }
      AliITSRecPoint *recp1 = (AliITSRecPoint*)prpl1->At(j);
      recp1->GetGlobalXYZ(gc1f);
      for(Int_t ico=0;ico<3;ico++)gc1[ico]=gc1f[ico];

      Double_t phi1 = TMath::ATan2(gc1[1]-ybeam,gc1[0]-xbeam);
      if(phi1<0)phi1=2*TMath::Pi()+phi1;
      for(Int_t ladl2=0 ; ladl2<fLadOnLay2*2+1;ladl2++){
	for(Int_t k=0;k<4;k++){
	  Int_t ladmod=fLadders[ladder-1]+ladl2;
 	  if(ladmod>AliITSgeomTGeo::GetNLadders(2)) ladmod=ladmod-AliITSgeomTGeo::GetNLadders(2);
	  Int_t modul2=AliITSgeomTGeo::GetModuleIndex(2,ladmod,k+1);
	  if(modul2<0)continue;
	  if(!fUseModule[modul2]) continue;
	  itsRec=rpcont->UncheckedGetClusters(modul2);
	  Int_t nrecp2 = itsRec->GetEntries();
	  for(Int_t j2=0;j2<nrecp2;j2++){
	    if(j2>kMaxCluPerMod) continue;
	    UShort_t idClu2=modul2*kMaxCluPerMod+j2;
	    if(fUsedCluster.TestBitNumber(idClu2)) continue;

	    AliITSRecPoint *recp2 = (AliITSRecPoint*)itsRec->At(j2);
	    recp2->GetGlobalXYZ(gc2f);
	    for(Int_t ico=0;ico<3;ico++)gc2[ico]=gc2f[ico];
	    Double_t phi2 = TMath::ATan2(gc2[1]-ybeam,gc2[0]-xbeam);
	    if(phi2<0)phi2=2*TMath::Pi()+phi2;
	    Double_t diff = TMath::Abs(phi2-phi1); 
	    if(diff>TMath::Pi())diff=2.*TMath::Pi()-diff; 
	    if(optCuts==0 && diff<fDiffPhiforPileup){
	      Double_t r1=TMath::Sqrt(gc1[0]*gc1[0]+gc1[1]*gc1[1]);
	      Double_t zc1=gc1[2];
	      Double_t r2=TMath::Sqrt(gc2[0]*gc2[0]+gc2[1]*gc2[1]);
	      Double_t zc2=gc2[2];
	      Double_t zr0=(r2*zc1-r1*zc2)/(r2-r1); //Z @ null radius
	      fZHisto->Fill(zr0);
	    }
	    if(diff>deltaPhi)continue;
	    AliStrLine line(gc1,gc2,kTRUE);
	    Double_t cp[3];
	    Int_t retcode = line.Cross(&zeta,cp);
	    if(retcode<0)continue;
	    Double_t dca = line.GetDCA(&zeta);
	    if(dca<0.) continue;
	    if(dca>deltaR)continue;
	    Double_t deltaZ=cp[2]-zvert;
	    if(TMath::Abs(deltaZ)>dZmax)continue;


	    if(nolines == 0){
	      if(fLines.GetEntriesFast()>0)fLines.Clear("C");
	    }
	    Float_t cov[6];
	    recp2->GetGlobalCov(cov);


	    Double_t rad1=TMath::Sqrt(gc1[0]*gc1[0]+gc1[1]*gc1[1]);
	    Double_t rad2=TMath::Sqrt(gc2[0]*gc2[0]+gc2[1]*gc2[1]);
	    Double_t factor=(rad1+rad2)/(rad2-rad1); //factor to account for error on tracklet direction 

	    Double_t curvErr=0;
	    if(bField>0.00001){
	      Double_t curvRadius=fMeanPtSelTrk/(0.3*bField)*100; //cm 
	      Double_t dRad=TMath::Sqrt((gc1[0]-gc2[0])*(gc1[0]-gc2[0])+(gc1[1]-gc2[1])*(gc1[1]-gc2[1]));
	      Double_t aux=dRad/2.+rad1;
	      curvErr=TMath::Sqrt(curvRadius*curvRadius-dRad*dRad/4.)-TMath::Sqrt(curvRadius*curvRadius-aux*aux); //cm
	    }
	    Double_t sigmasq[3];
	    sigmasq[0]=(cov[0]+curvErr*curvErr/2.)*factor*factor;
	    sigmasq[1]=(cov[3]+curvErr*curvErr/2.)*factor*factor;
	    sigmasq[2]=cov[5]*factor*factor;

	    // Multiple scattering
	    Double_t pOverMass=fMeanPSelTrk/0.140;
 	    Double_t beta2=pOverMass*pOverMass/(1+pOverMass*pOverMass);
 	    Double_t p2=fMeanPSelTrk*fMeanPSelTrk;
 	    Double_t rBP=GetPipeRadius();
 	    Double_t dBP=0.08/35.3; // 800 um of Be
 	    Double_t dL1=0.01; //approx. 1% of radiation length  
 	    Double_t theta2BP=14.1*14.1/(beta2*p2*1e6)*dBP;
 	    Double_t theta2L1=14.1*14.1/(beta2*p2*1e6)*dL1;
	    Double_t rtantheta1=(rad2-rad1)*TMath::Tan(TMath::Sqrt(theta2L1));
	    Double_t rtanthetaBP=(rad1-rBP)*TMath::Tan(TMath::Sqrt(theta2BP));
 	    for(Int_t ico=0; ico<3;ico++){    
 	      sigmasq[ico]+=rtantheta1*rtantheta1*factor*factor/3.;
 	      sigmasq[ico]+=rtanthetaBP*rtanthetaBP*factor*factor/3.;
 	    }
	    Double_t wmat[9]={1.,0.,0.,0.,1.,0.,0.,0.,1.};
	    if(sigmasq[0]!=0.) wmat[0]=1./sigmasq[0];
	    if(sigmasq[1]!=0.) wmat[4]=1./sigmasq[1];
	    if(sigmasq[2]!=0.) wmat[8]=1./sigmasq[2];
	    new(fLines[nolines++])AliStrLine(gc1,sigmasq,wmat,gc2,kTRUE,idClu1,idClu2);

	  }
	}
      }
    }
  }

  SetLaddersOnLayer2(origLaddersOnLayer2);

  if(nolines == 0)return -2;
  return nolines;
}

//______________________________________________________________________
Int_t  AliITSVertexer3D::Prepare3DVertex(Int_t optCuts){
  // Finds the 3D vertex information using tracklets
  Int_t retcode = -1;
  Double_t xbeam=GetNominalPos()[0];
  Double_t ybeam=GetNominalPos()[1];
  Double_t zvert=0.;
  Double_t deltaR=fCoarseMaxRCut;
  Double_t dZmax=fZCutDiamond;
  if(optCuts==1){
    xbeam=fVert3D.GetXv();
    ybeam=fVert3D.GetYv();
    zvert=fVert3D.GetZv();
    deltaR=fMaxRCut;
    dZmax=fMaxZCut;
    if(fPileupAlgo == 2){ 
      dZmax=fZCutDiamond;
      deltaR=fMaxRCut2;
    }
  }else if(optCuts==2){
    xbeam=fVert3D.GetXv();
    ybeam=fVert3D.GetYv();
    deltaR=fMaxRCut;
  }

  Double_t origBinSizeR=fBinSizeR;
  Double_t origBinSizeZ=fBinSizeZ;
  Bool_t rebinned=kFALSE;
  if(fDoDownScale){
    SetBinSizeR(0.05);
    SetBinSizeZ(0.05);
    rebinned=kTRUE;
  }else{
    if(optCuts==0 && (fNRecPLay1>fMaxNumOfClForRebin || fNRecPLay2>fMaxNumOfClForRebin)){
      SetBinSizeR(0.1);
      SetBinSizeZ(0.2);
      rebinned=kTRUE;
    }
  }
  Double_t rl=-fCoarseMaxRCut;
  Double_t rh=fCoarseMaxRCut;
  Double_t zl=-fZCutDiamond;
  Double_t zh=fZCutDiamond;
  Int_t nbr=(Int_t)((rh-rl)/fBinSizeR+0.0001);
  Int_t nbz=(Int_t)((zh-zl)/fBinSizeZ+0.0001);
  Int_t nbrcs=(Int_t)((rh-rl)/(fBinSizeR*2.)+0.0001);
  Int_t nbzcs=(Int_t)((zh-zl)/(fBinSizeZ*2.)+0.0001);

  TH3F *h3d = new TH3F("h3d","xyz distribution",nbr,rl,rh,nbr,rl,rh,nbz,zl,zh);
  TH3F *h3dcs = new TH3F("h3dcs","xyz distribution",nbrcs,rl,rh,nbrcs,rl,rh,nbzcs,zl,zh);

  // cleanup of the TCLonesArray of tracklets (i.e. fakes are removed)
  Int_t vsiz = fLines.GetEntriesFast();
  Int_t *validate = new Int_t [vsiz];
  for(Int_t i=0; i<vsiz;i++)validate[i]=0;
  for(Int_t i=0; i<vsiz-1;i++){
    AliStrLine *l1 = (AliStrLine*)fLines.At(i);
    for(Int_t j=i+1;j<fLines.GetEntriesFast();j++){
      AliStrLine *l2 = (AliStrLine*)fLines.At(j);
      Double_t dca=l1->GetDCA(l2);
      if(dca > fDCAcut || dca<0.00001) continue;
      Double_t point[3];
      Int_t retc = l1->Cross(l2,point);
      if(retc<0)continue;
      Double_t deltaZ=point[2]-zvert;
      if(TMath::Abs(deltaZ)>dZmax)continue;
      Double_t rad=TMath::Sqrt(point[0]*point[0]+point[1]*point[1]);
      if(rad>fCoarseMaxRCut)continue;
      Double_t deltaX=point[0]-xbeam;
      Double_t deltaY=point[1]-ybeam;
      Double_t raddist=TMath::Sqrt(deltaX*deltaX+deltaY*deltaY);
      if(raddist>deltaR)continue;
      validate[i]=1;
      validate[j]=1;
      h3d->Fill(point[0],point[1],point[2]);
      h3dcs->Fill(point[0],point[1],point[2]);
    }
  }

  Int_t numbtracklets=0;
  for(Int_t i=0; i<vsiz;i++)if(validate[i]>=1)numbtracklets++;
  if(numbtracklets<2){
    delete [] validate; 
    delete h3d; 
    delete h3dcs; 
    SetBinSizeR(origBinSizeR);
    SetBinSizeZ(origBinSizeZ);
    return retcode; 
  }

  for(Int_t i=0; i<fLines.GetEntriesFast();i++){
    if(validate[i]<1)fLines.RemoveAt(i);
  }
  fLines.Compress();
  AliDebug(1,Form("Number of tracklets (after compress)%d ",fLines.GetEntriesFast()));
  delete [] validate;

  // Exit here if Pileup Algorithm 2 has been chosen during second loop
  if(fPileupAlgo == 2 && optCuts==1){
    delete h3d; 
    delete h3dcs;     
    SetBinSizeR(origBinSizeR);
    SetBinSizeZ(origBinSizeZ);
    return 0;
  }

  //        Find peaks in histos

  Double_t peak[3]={0.,0.,0.};
  Int_t ntrkl,ntimes;
  FindPeaks(h3d,peak,ntrkl,ntimes);  
  delete h3d;
  Double_t binsizer=(rh-rl)/nbr;
  Double_t binsizez=(zh-zl)/nbz;
  if(optCuts==0 && (ntrkl<=2 || ntimes>1)){
    ntrkl=0;
    ntimes=0;
    FindPeaks(h3dcs,peak,ntrkl,ntimes);  
    binsizer=(rh-rl)/nbrcs;
    binsizez=(zh-zl)/nbzcs;
    if(ntrkl==1 || ntimes>1){
      delete h3dcs; 
      SetBinSizeR(origBinSizeR);
      SetBinSizeZ(origBinSizeZ);
      return retcode;
    }
  }
  delete h3dcs;

  Double_t bs=(binsizer+binsizez)/2.;
  for(Int_t i=0; i<fLines.GetEntriesFast();i++){
    AliStrLine *l1 = (AliStrLine*)fLines.At(i);
    if(l1->GetDistFromPoint(peak)>2.5*bs)fLines.RemoveAt(i);
  }
  fLines.Compress();
  AliDebug(1,Form("Number of tracklets (after 2nd compression) %d",fLines.GetEntriesFast()));

  // Finer Histo in limited range in case of high mult.
  if(rebinned){
    SetBinSizeR(0.01);
    SetBinSizeZ(0.01);
    Double_t xl=peak[0]-0.3;
    Double_t xh=peak[0]+0.3;
    Double_t yl=peak[1]-0.3;
    Double_t yh=peak[1]+0.3;
    zl=peak[2]-0.5;
    zh=peak[2]+0.5;
    Int_t nbxfs=(Int_t)((xh-xl)/fBinSizeR+0.0001);
    Int_t nbyfs=(Int_t)((yh-yl)/fBinSizeR+0.0001);
    Int_t nbzfs=(Int_t)((zh-zl)/fBinSizeZ+0.0001);

    TH3F *h3dfs = new TH3F("h3dfs","xyz distribution",nbxfs,xl,xh,nbyfs,yl,yh,nbzfs,zl,zh);
    for(Int_t i=0; i<fLines.GetEntriesFast()-1;i++){
      AliStrLine *l1 = (AliStrLine*)fLines.At(i);
      for(Int_t j=i+1;j<fLines.GetEntriesFast();j++){
	AliStrLine *l2 = (AliStrLine*)fLines.At(j);
	Double_t dca=l1->GetDCA(l2);
	if(dca > fDCAcut || dca<0.00001) continue;
	Double_t point[3];
	Int_t retc = l1->Cross(l2,point);
	if(retc<0)continue;
	Double_t deltaZ=point[2]-zvert;
	if(TMath::Abs(deltaZ)>dZmax)continue;
	Double_t rad=TMath::Sqrt(point[0]*point[0]+point[1]*point[1]);
	if(rad>fCoarseMaxRCut)continue;
	Double_t deltaX=point[0]-xbeam;
	Double_t deltaY=point[1]-ybeam;
	Double_t raddist=TMath::Sqrt(deltaX*deltaX+deltaY*deltaY);
	if(raddist>deltaR)continue;
	h3dfs->Fill(point[0],point[1],point[2]);
      }
    }
    ntrkl=0;
    ntimes=0;

    Double_t newpeak[3]={0.,0.,0.};
    FindPeaks(h3dfs,newpeak,ntrkl,ntimes);  
    if(ntimes==1){
      for(Int_t iCoo=0; iCoo<3; iCoo++) peak[iCoo]=newpeak[iCoo];
      binsizer=fBinSizeR;
      binsizez=fBinSizeZ;
    }
    delete h3dfs;
    bs=(binsizer+binsizez)/2.;
    for(Int_t i=0; i<fLines.GetEntriesFast();i++){
      AliStrLine *l1 = (AliStrLine*)fLines.At(i);
      if(l1->GetDistFromPoint(peak)>2.5*bs)fLines.RemoveAt(i);
    }
    fLines.Compress();
    AliDebug(1,Form("Number of tracklets (after 3rd compression) %d",fLines.GetEntriesFast()));
  }
  SetBinSizeR(origBinSizeR);
  SetBinSizeZ(origBinSizeZ);


  //         Second selection loop


  if(fLines.GetEntriesFast()>1){
    retcode=0;
    //  find a first candidate for the primary vertex
    fVert3D=AliVertexerTracks::TrackletVertexFinder(&fLines,0); 
    // make a further selection on tracklets based on this first candidate
    fVert3D.GetXYZ(peak);
    AliDebug(1,Form("FIRST V candidate: %f ; %f ; %f",peak[0],peak[1],peak[2]));
    Int_t *validate2 = new Int_t [fLines.GetEntriesFast()];
    for(Int_t i=0; i<fLines.GetEntriesFast();i++) validate2[i]=1; 
    for(Int_t i=0; i<fLines.GetEntriesFast();i++){
      if(validate2[i]==0) continue; 
      AliStrLine *l1 = (AliStrLine*)fLines.At(i);
      if(l1->GetDistFromPoint(peak)> fDCAcut)fLines.RemoveAt(i);
      if(optCuts==2){ // temporarily only for pileup
	for(Int_t j=i+1; j<fLines.GetEntriesFast();j++){
	  AliStrLine *l2 = (AliStrLine*)fLines.At(j);
	  if(l1->GetDCA(l2)<0.00001){ 
	    Int_t delta1=(Int_t)l1->GetIdPoint(0)-(Int_t)l2->GetIdPoint(0);
	    Int_t delta2=(Int_t)l1->GetIdPoint(1)-(Int_t)l2->GetIdPoint(1);
	    Int_t deltamod1=(Int_t)l1->GetIdPoint(0)/kMaxCluPerMod
	      -(Int_t)l2->GetIdPoint(0)/kMaxCluPerMod;
	    Int_t deltamod2=(Int_t)l1->GetIdPoint(1)/kMaxCluPerMod
	      -(Int_t)l2->GetIdPoint(1)/kMaxCluPerMod;
	    // remove tracklets sharing a point
	    if( (delta1==0 && deltamod2==0)  || 
		(delta2==0 && deltamod1==0)  ) validate2[j]=0; 
	  }
	}
      }
    }
    for(Int_t i=0; i<fLines.GetEntriesFast();i++){
      if(validate2[i]==0)  fLines.RemoveAt(i);
    }
    delete [] validate2;
    fLines.Compress();
    AliDebug(1,Form("Number of tracklets (after 3rd compression) %d",fLines.GetEntriesFast()));
    if(fLines.GetEntriesFast()>1){// this new tracklet selection is used
      fVert3D=AliVertexerTracks::TrackletVertexFinder(&fLines,0);
    }
  }
  return retcode;  
}

//________________________________________________________
Int_t  AliITSVertexer3D::Prepare3DVertexPbPb(){
  // Finds the 3D vertex information in Pb-Pb events using tracklets
  AliDebug(1,"High multiplicity event.\n");

  Int_t nxy=(Int_t)(2.*fCoarseMaxRCut/f3DBinSize);
  Double_t xymi= -nxy*f3DBinSize/2.;
  Double_t xyma= nxy*f3DBinSize/2.;
  Int_t nz=(Int_t)(2.*fZCutDiamond/f3DBinSize);
  Double_t zmi=-nz*f3DBinSize/2.;
  Double_t zma=nz*f3DBinSize/2.;
  Int_t nolines=fLines.GetEntriesFast();
  TH3F *h3dv = new TH3F("h3dv","3d tracklets",nxy,xymi,xyma,nxy,xymi,xyma,nz,zmi,zma);
  
  for(Int_t itra=0; itra<nolines; itra++){
    Double_t wei = GetFraction(itra);
    //printf("tracklet %d ) - weight %f \n",itra,wei);
    if(wei>1.e-6){
      AliStrLine *str=(AliStrLine*)fLines.At(itra);
      Double_t t1,t2;
      if(str->GetParamAtRadius(fCoarseMaxRCut,t1,t2)){
	do{
	  Double_t punt[3];
	  str->ComputePointAtT(t1,punt);
	  h3dv->Fill(punt[0],punt[1],punt[2],wei);
	  t1+=f3DBinSize/3.;
	} while(t1<t2);
      }
    }
  }
  Int_t noftrk,noftim;
  FindPeaks(h3dv,f3DPeak,noftrk,noftim); // arg: histo3d, peak, # of contrib., # of other peak with same magnitude
  
  
  // Remove all the tracklets which are not passing near peak
  
  while(nolines--){
    AliStrLine *str=(AliStrLine*)fLines.At(nolines);
    Double_t dist = str->GetDistFromPoint(f3DPeak);
    if(dist>(2.*f3DBinSize)) fLines.RemoveAt(nolines);
    }
  fLines.Compress();
  nolines=fLines.GetEntriesFast();

  delete h3dv;

  Int_t *validate2 = new Int_t [fLines.GetEntriesFast()];
  for(Int_t i=0; i<fLines.GetEntriesFast();i++) validate2[i]=1; 
  for(Int_t i=0; i<fLines.GetEntriesFast();i++){
    if(validate2[i]==0) continue; 
    AliStrLine *l1 = (AliStrLine*)fLines.At(i);
    if(l1->GetDistFromPoint(f3DPeak)> fDCAcut)fLines.RemoveAt(i);
    for(Int_t j=i+1; j<fLines.GetEntriesFast();j++){
      AliStrLine *l2 = (AliStrLine*)fLines.At(j);
      if(l1->GetDCA(l2)<0.00001){ 
	Int_t delta1=(Int_t)l1->GetIdPoint(0)-(Int_t)l2->GetIdPoint(0);
	Int_t delta2=(Int_t)l1->GetIdPoint(1)-(Int_t)l2->GetIdPoint(1);
	Int_t deltamod1=(Int_t)l1->GetIdPoint(0)/kMaxCluPerMod
	  -(Int_t)l2->GetIdPoint(0)/kMaxCluPerMod;
	Int_t deltamod2=(Int_t)l1->GetIdPoint(1)/kMaxCluPerMod
	  -(Int_t)l2->GetIdPoint(1)/kMaxCluPerMod;
	// remove tracklets sharing a point
	if( (delta1==0 && deltamod2==0)  || 
	    (delta2==0 && deltamod1==0)  ) validate2[j]=0; 
	
      }
    }
  }
  for(Int_t i=0; i<fLines.GetEntriesFast();i++){
    if(validate2[i]==0)  fLines.RemoveAt(i);
  }
  
  delete [] validate2;
  fLines.Compress();

  
  AliDebug(1,Form("Number of tracklets (after 3rd compression) %d",fLines.GetEntriesFast()));

  fVert3D=AliVertexerTracks::TrackletVertexFinder(&fLines,0); 
  fVert3D.GetXYZ(f3DPeak);
  
  return 0;  
}

//________________________________________________________
void AliITSVertexer3D::SetMeanPPtSelTracks(Double_t fieldTesla){
  // Sets mean values of Pt based on the field
  // for P (used in multiple scattering) the most probable value is used
  if(TMath::Abs(fieldTesla-0.5)<0.01){
    SetMeanPSelTracks(0.375);
    SetMeanPtSelTracks(0.630);
  }else if(TMath::Abs(fieldTesla-0.4)<0.01){
    SetMeanPSelTracks(0.375);
    SetMeanPtSelTracks(0.580);
  }else if(TMath::Abs(fieldTesla-0.2)<0.01){
    SetMeanPSelTracks(0.375);
    SetMeanPtSelTracks(0.530);
  }else if(fieldTesla<0.00001){
    SetMeanPSelTracks(0.375);
    SetMeanPtSelTracks(0.230);
  }else{
    SetMeanPSelTracks();
    SetMeanPtSelTracks();
  }
}

//________________________________________________________
void AliITSVertexer3D::FindPeaks(TH3F* histo, Double_t *peak, Int_t &nOfTracklets, Int_t &nOfTimes){
  // Finds bin with max contents in 3D histo of tracket intersections
  TAxis *xax = histo->GetXaxis();  
  TAxis *yax = histo->GetYaxis();
  TAxis *zax = histo->GetZaxis();
  peak[0]=0.;
  peak[1]=0.;
  peak[2]=0.;
  nOfTracklets = 0;
  nOfTimes=0;
  Int_t peakbin[3]={0,0,0};
  Int_t peak2bin[3]={-1,-1,-1};
  Int_t bc2=-1;
  for(Int_t i=xax->GetFirst();i<=xax->GetLast();i++){
    Double_t xval = xax->GetBinCenter(i);
    for(Int_t j=yax->GetFirst();j<=yax->GetLast();j++){
      Double_t yval = yax->GetBinCenter(j);
      for(Int_t k=zax->GetFirst();k<=zax->GetLast();k++){
	Double_t zval = zax->GetBinCenter(k);
	Int_t bc =(Int_t)histo->GetBinContent(i,j,k);
	if(bc==0) continue;
	if(bc>nOfTracklets){
	  nOfTracklets=bc;
	  peak[2] = zval;
	  peak[1] = yval;
	  peak[0] = xval;
	  peakbin[2] = k;
	  peakbin[1] = j;
	  peakbin[0] = i;
	  peak2bin[2] = -1;
	  peak2bin[1] = -1;
	  peak2bin[0] = -1;
	  bc2=-1;
	  nOfTimes = 1;
	}else if(bc==nOfTracklets){
	  if(TMath::Abs(i-peakbin[0])<=1 && TMath::Abs(j-peakbin[1])<=1 && TMath::Abs(k-peakbin[2])<=1){
	    peak2bin[2] = k;
	    peak2bin[1] = j;
	    peak2bin[0] = i;
	    bc2=bc;
	    nOfTimes = 1;
	  }else{
	    nOfTimes++;
	  }
	}
      }
    }
  }
  if(peak2bin[0]>=-1 && bc2!=-1){ // two contiguous peak-cells with same contents
    peak[0]=0.5*(xax->GetBinCenter(peakbin[0])+xax->GetBinCenter(peak2bin[0]));
    peak[1]=0.5*(yax->GetBinCenter(peakbin[1])+yax->GetBinCenter(peak2bin[1]));
    peak[2]=0.5*(zax->GetBinCenter(peakbin[2])+zax->GetBinCenter(peak2bin[2]));
    nOfTracklets+=bc2;
    nOfTimes=1;
  }
}
//________________________________________________________
void AliITSVertexer3D::MarkUsedClusters(){
  // Mark clusters of tracklets used in vertex claulation
  for(Int_t i=0; i<fLines.GetEntriesFast();i++){
    AliStrLine *lin = (AliStrLine*)fLines.At(i);
    Int_t idClu1=lin->GetIdPoint(0);
    Int_t idClu2=lin->GetIdPoint(1);
    fUsedCluster.SetBitNumber(idClu1);
    fUsedCluster.SetBitNumber(idClu2);
  }
}
//________________________________________________________
Int_t AliITSVertexer3D::RemoveTracklets(){
  // Remove trackelts close to first found vertex
  Double_t vert[3]={fVert3D.GetXv(),fVert3D.GetYv(),fVert3D.GetZv()};
  Int_t nRemoved=0;
  for(Int_t i=0; i<fLines.GetEntriesFast();i++){
    AliStrLine *lin = (AliStrLine*)fLines.At(i);
    if(lin->GetDistFromPoint(vert)<fDCAforPileup){
      Int_t idClu1=lin->GetIdPoint(0);
      Int_t idClu2=lin->GetIdPoint(1);
      fUsedCluster.SetBitNumber(idClu1);
      fUsedCluster.SetBitNumber(idClu2);
      fLines.RemoveAt(i);
      ++nRemoved;
    }
  }
  fLines.Compress();
  return nRemoved;
}
//________________________________________________________
void AliITSVertexer3D::FindOther3DVertices(TTree *itsClusterTree){
  // pileup identification based on 3D vertexing with not used clusters

  fVertArray = new AliESDVertex[kMaxPileupVertices+1];
  fVertArray[0]=(*fCurrentVertex);
  Int_t nFoundVert=1;
  for(Int_t iPilV=1; iPilV<=kMaxPileupVertices; iPilV++){
    MarkUsedClusters();
    fLines.Clear("C");
    Int_t nolines = FindTracklets(itsClusterTree,2);
    if(nolines>=2){
      Int_t nr=RemoveTracklets();
      nolines-=nr;
      if(nolines>=2){
	Int_t rc=Prepare3DVertex(2);
	if(rc==0){ 
	  fVert3D=AliVertexerTracks::TrackletVertexFinder(&fLines,0);
	  if(fVert3D.GetNContributors()>=fMinTrackletsForPilup){
	    fIsPileup=kTRUE;
	    fVertArray[nFoundVert]=fVert3D;
	    nFoundVert++;
	    if(nFoundVert==2){
	      fZpuv=fVert3D.GetZv();
	      fNTrpuv=fVert3D.GetNContributors();
	    }
	  }
	}
      }
    }
  }
  fNoVertices=nFoundVert;
}
//______________________________________________________________________
void AliITSVertexer3D::PileupFromZ(){
  // Calls the pileup algorithm of ALiITSVertexerZ
  Int_t binmin, binmax;
  Int_t nPeaks=AliITSVertexerZ::GetPeakRegion(fZHisto,binmin,binmax);   
  if(nPeaks==2)AliWarning("2 peaks found");
  Int_t firstPeakCont=0;
  Double_t firstPeakPos=0.;
  for(Int_t i=binmin-1;i<=binmax+1;i++){
    firstPeakCont+=static_cast<Int_t>(fZHisto->GetBinContent(i));
    firstPeakPos+=fZHisto->GetBinContent(i)*fZHisto->GetBinCenter(i);
  }
  if(firstPeakCont>0){ 
    firstPeakPos/=firstPeakCont;
    Int_t ncontr2=0;
    if(firstPeakCont>fMinTrackletsForPilup){     
      Float_t secPeakPos;
      ncontr2=AliITSVertexerZ::FindSecondPeak(fZHisto,binmin,binmax,secPeakPos);
      if(ncontr2>=fMinTrackletsForPilup){ 
	fIsPileup=kTRUE;
	fNoVertices=2;
	AliESDVertex secondVert(secPeakPos,0.1,ncontr2);
	fVertArray = new AliESDVertex[2];
	fVertArray[0]=(*fCurrentVertex);
	fVertArray[1]=secondVert;
	fZpuv=secPeakPos;
	fNTrpuv=ncontr2;
      }
    }
  }
}

//________________________________________________________
Double_t AliITSVertexer3D::GetFraction(Int_t itr) const {
  // this method is used to fill a 3D histogram representing
  // the trajectories of the candidate tracklets
  // The computed fraction is used as a weight at filling time
  AliStrLine *str = (AliStrLine*)fLines.At(itr);
  Double_t spigolo=10.;
  Double_t cd[3];
  str->GetCd(cd);
  Double_t par=0.;
  Double_t maxl=TMath::Sqrt(3.)*spigolo;
 // intersection with a plane normal to the X axis 
  if(TMath::AreEqualAbs(cd[0],0.,1.e-9)){
    par=1000000.;
  }
  else {
    par=spigolo/cd[0];
  }
  Double_t zc=cd[2]*par;
  Double_t yc=cd[1]*par;
  if((-spigolo<=yc && yc<=spigolo) && (-spigolo<=zc && zc<=spigolo))return TMath::Abs(par/maxl);
 // intersection with a plane normal to the Y axis
  if(TMath::AreEqualAbs(cd[1],0.,1.e-9)){
    par=1000000.;
  }
  else {
    par=spigolo/cd[1];
  }
  zc=cd[2]*par;
  Double_t xc=cd[0]*par;
  if((-spigolo<=xc && xc<=spigolo) && (-spigolo<=zc && zc<=spigolo))return TMath::Abs(par/maxl);
 // intersection with a plane normal to the Z axis
  if(TMath::AreEqualAbs(cd[2],0.,1.e-9)){
    par=1000000.;
  }
  else {
    par=spigolo/cd[2];
  }
  yc=cd[1]*par;
  xc=cd[0]*par;
  if((-spigolo<=xc && xc<=spigolo) && (-spigolo<=yc && yc<=spigolo))return TMath::Abs(par/maxl);
  // control should never reach the following lines
  AliError(Form("anomalous tracklet direction for tracklet %d in fLines\n",itr));
  str->PrintStatus();
  return 0.;
}

//________________________________________________________
void AliITSVertexer3D::PrintStatus() const {
  // Print current status
  printf("========= First step selections =====================\n");
  printf("Cut on diamond (Z) %f\n",fZCutDiamond);
  printf("Loose cut on Delta Phi %f\n",fCoarseDiffPhiCut);
  printf("Loose cut on tracklet DCA to Z axis %f\n",fCoarseMaxRCut);
  printf("Cut on DCA - tracklet to tracklet and to vertex %f\n",fDCAcut);
  printf("========= Second step selections ====================\n");
  printf("Cut on tracklet-to-first-vertex Z distance %f\n",fMaxZCut);
  printf("Max Phi difference: %f\n",fDiffPhiMax);
  printf("Cut on tracklet DCA to beam axis %f\n",fMaxRCut);
  printf("Cut on tracklet DCA to beam axis (algo2) %f\n",fMaxRCut2);
  printf("========= Pileup selections =========================\n");
  printf("Pileup algo: %d\n",fPileupAlgo);
  printf("Min DCA to 1st vertex for pileup (algo 0 and 1): %f\n",fDCAforPileup);
  printf("Cut on distance between pair-vertices  (algo 2): %f\n",fCutOnPairs);
  printf("Maximum number of clusters on L1 or L2 for downscale: %d\n",fMaxNumOfClForDownScale);
  printf("Maximum number of clusters on L1 or L2 for histo rebin: %d\n",fMaxNumOfClForRebin);
  printf("=======================================================\n");
}

