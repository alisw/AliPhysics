/**************************************************************************
 * Copyright(c) 2009-2010, ALICE Experiment at CERN, All rights reserved. *
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
#include<Riostream.h>
#include<TClonesArray.h>
#include<TMath.h>
#include "AliStrLine.h"
#include "AliITSSortTrkl.h"

/* $Id$ */

////////////////////////////////////////////////////////////////////////
//           Helper class for finding multiple primary vertices       //
//           To be used by AliITSVertexer3D                           //
//           It is based on the association of pairs of tracklets     //
//           obtained by matching reconstructed points onthe first    //
//           2 layers of SPD                                          //
//           Origin M. Masera masera@to.infn.it                       //
////////////////////////////////////////////////////////////////////////


ClassImp(AliITSSortTrkl)

//______________________________________________________________________
AliITSSortTrkl::AliITSSortTrkl():TObject(),
fkSize(0),
fIndex(0),
fPairs(NULL),
fClustersTmp(NULL),
fClusters(NULL),
fNoClus(0),
fSize(NULL),
fMark(),
fCut(0.),
fCoarseMaxRCut(0.)
{
  // Default constructor
}

//______________________________________________________________________
AliITSSortTrkl::AliITSSortTrkl(Int_t n, Double_t cut):TObject(),
fkSize(n),
fIndex(0),
fPairs(),
fClustersTmp(),
fClusters(NULL),
fNoClus(0),
fSize(NULL),
fMark(n),
fCut(cut),
fCoarseMaxRCut(0.) {
  // Standard constructor
  fPairs = new AliITSTracklPairs* [fkSize];
  fClustersTmp = new Int_t* [fkSize-1]; 
}

//______________________________________________________________________
AliITSSortTrkl::AliITSSortTrkl(TClonesArray &tclo, Int_t n, Double_t cut, Double_t rcut):TObject(),
fkSize(n),
fIndex(0),
fPairs(),
fClustersTmp(),
fClusters(NULL),
fNoClus(0),
fSize(NULL),
fMark(n),
fCut(cut),
fCoarseMaxRCut(rcut){
  // Constructor based on a TClonesArray of AliStrLine
  Int_t numtrack=tclo.GetEntriesFast();
  if(numtrack<2){
    AliFatal(Form("Insufficient number of tracks (%d )",numtrack));
    return;
  }
  TObject *member = tclo[0];
  TString str(member->ClassName());
  if(!(str.Contains("AliStrLine"))){
    AliFatal(Form("Wrong type of class in input TClonesArray (%s )",str.Data()));
    return;
  }
  Int_t siz = numtrack*(numtrack-1)/2;
  if(fkSize < siz){
    AliError(Form("fkSize is too small. It is %d and it should be %d",fkSize,siz)); 
  }
  fPairs = new AliITSTracklPairs* [fkSize];
  fClustersTmp = new Int_t* [fkSize-1];
  for(Int_t i=0;i<numtrack-1;i++){
    AliStrLine *one = (AliStrLine*)tclo[i];
    for(Int_t j=i+1;j<numtrack;j++){
      AliStrLine *two = (AliStrLine*)tclo[j];
      Double_t point[3];
      one->Cross(two,point);
      Double_t dca = one->GetDCA(two);
      if(dca>fCut)continue;
      Double_t rad=TMath::Sqrt(point[0]*point[0]+point[1]*point[1]);
      if(rad>fCoarseMaxRCut)continue;
      AddPairs(i,j,dca,point);
    }
  } 
}

//______________________________________________________________________
AliITSSortTrkl::~AliITSSortTrkl(){
  // Destructor
  if(fPairs){
    for(Int_t i=0;i<fIndex;i++)delete fPairs[i];
    delete [] fPairs;
  }
  DeleteClustersTmp();
  if(fClusters){
    for(Int_t i=0;i<fNoClus;i++)delete []fClusters[i];
    delete fClusters;
    delete []fSize;
  } 
}

//______________________________________________________________________
void AliITSSortTrkl::DeleteClustersTmp(){
  // fClustersTmp is deleted
  if(fClustersTmp){
    for(Int_t i=0;i<fIndex-1;i++)delete []fClustersTmp[i];
    delete fClustersTmp;
    fClustersTmp = NULL;
  }
}

//______________________________________________________________________
Int_t AliITSSortTrkl::AddPairs(Int_t t1, Int_t t2, Double_t dca, Double_t *coo){
  // Add a tracklet pair at current position
  fPairs[fIndex] = new AliITSTracklPairs(t1,t2,dca,coo);
  //  cout<<"Pair "<<fIndex<<" Tracklets "<<t1<<" and "<<t2<<". DCA= "<<dca<<" Crossing "<<coo[0]<<" "<<coo[1]<<" "<<coo[2]<<endl;
  fIndex++;
  return fIndex-1;
}

//______________________________________________________________________
Int_t AliITSSortTrkl::FindClusters(){
  // find clusters
  fMark.ResetAllBits();
  PrepareClustersTmp();
  for(Int_t i=0; i<fIndex-1;i++){
    Int_t *v=fClustersTmp[i];
    Clustering(i,v);
    if(v[0]>0){
      AliInfo(Form("Clusters starting from pair %d : %d",i,v[0]));
      Int_t dim;
      v[v[0]+1]=i;
      Int_t *arr=FindLabels(v+1,2*(v[0]+1),dim);
      cout<<"Pairs involved \n";
      for(Int_t j=1;j<=v[0];j++){
	cout<<v[j]<<" ";
	if(j%10==0)cout<<endl;
      }
      cout<<endl;
      cout<<"Tracklets involved\n";
      for(Int_t j=0;j<dim;j++){
	cout<<arr[j]<<" ";
	if(j%10==0 && j>0)cout<<endl;
      }
      cout<<endl;
      for(Int_t j=0;j<fIndex;j++){
	if(fMark.TestBitNumber(j))continue;
	for(Int_t k=0;k<dim;k++){
	  if(fPairs[j]->HasTrack(arr[k])){
	    fMark.SetBitNumber(j);
	    //	    cout<<"Marked pair "<<j<<" because has tracklet "<<arr[k]<<endl;
	    k=dim;
	  }
	}
      }
      delete []arr;
    }
  }
  Cleanup();
  return fNoClus;
}

//______________________________________________________________________
void AliITSSortTrkl::Cleanup(){
  // empty arrays are eliminated, the others are sorted according
  // to cluster multiplicity
  Int_t siz[fIndex-1];
  Int_t index[fIndex-1];
  fNoClus=0;
  for(Int_t i=0; i<fIndex-1;i++){
    Int_t *v=fClustersTmp[i];
    if(v[0]>0)fNoClus++;
    siz[i]=v[0];
  }
  if(fNoClus == 0)return;
  TMath::Sort(fIndex-1,siz,index);
  fClusters = new Int_t* [fNoClus];
  fSize = new Int_t [fNoClus];
  for(Int_t i=0;i<fNoClus;i++){
    Int_t curind = index[i];
    Int_t *v=fClustersTmp[curind];
    fClusters[i] = new Int_t[v[0]+1];
    fSize[i]=v[0]+1;
    Int_t *vext=fClusters[i];
    vext[0]=curind;
    for(Int_t j=1;j<fSize[i];j++)vext[j]=v[j];
  }
}

//______________________________________________________________________
Int_t* AliITSSortTrkl::GetTrackletsLab(Int_t index, Int_t& dim) const {
  // Returns the tracklet labels corresponding to cluster index
  // Calling code must take care of memory deallocation
  if(fNoClus <=0){
    dim = 0;
    return NULL;
  }
  Int_t dimmax = 2*GetSizeOfCluster(index);
  if(dimmax<=0){
    dim = 0;
    return NULL;
  } 
  Int_t *v = GetClusters(index);
  return FindLabels(v,dimmax,dim);
}

//______________________________________________________________________
Int_t* AliITSSortTrkl::FindLabels(Int_t *v, Int_t dimmax, Int_t& dim) const {
  // Returns the tracklet labels corresponding to te list of pairs 
  // contained in v.
  // Calling code must take care of memory deallocation
  Int_t *arr = new Int_t [dimmax];
  Int_t j=0;
  for(Int_t i=0; i<dimmax/2; i++){
    AliITSTracklPairs *pai=GetPairsAt(v[i]);
    arr[j++]=pai->GetTrack1();
    arr[j++]=pai->GetTrack2();
  }
  SortAndClean(dimmax,arr,dim);
  return arr;
}
//______________________________________________________________________
void AliITSSortTrkl::SortAndClean(Int_t numb, Int_t *arr, Int_t& numb2){
  // Array arr (with numb elements) is sorted in ascending order. 
  // Then possible reoccurrences
  // of elements are eliminated. numb2 is the number of remaining elements
  // after cleanup.
  if(numb<=0)return;
  Int_t index[numb];
  TMath::Sort(numb,arr,index,kFALSE);
  Int_t tmp[numb];
  numb2 = 0;
  tmp[0] = arr[index[0]];
  for(Int_t i=1;i<numb;i++){
    if(arr[index[i]] != tmp[numb2]){
      ++numb2;
      tmp[numb2]=arr[index[i]];
    }
  }
  ++numb2;
  for(Int_t i=0;i<numb;i++){
    if(i<numb2){
      arr[i]=tmp[i];
    }
    else {
      arr[i]=0;
    }
  }
}

//______________________________________________________________________
void AliITSSortTrkl::PrepareClustersTmp(){
  // prepare arrays of clusters
  for(Int_t i=0; i<fIndex-1;i++){
    fClustersTmp[i] = new Int_t [fIndex-i];
    Int_t *v = fClustersTmp[i];
    v[0]=0;
    for(Int_t j=1;j<fIndex-i;j++)v[j]=-1;
  }
}

//______________________________________________________________________
void AliITSSortTrkl::Clustering(Int_t i,Int_t *v){
  // recursive method to build up clusters starting from point i
  if(fMark.TestBitNumber(i))return;
  AliITSTracklPairs *p1 = fPairs[i];
  for(Int_t j=i+1;j<fIndex;j++){
    if(fMark.TestBitNumber(j))continue;
    AliITSTracklPairs *p2 = fPairs[j];
    Double_t dist = p1->GetDistance(*p2);
    //    AliInfo(Form("  ******* i %d , j %d . Distance %g ",i,j,dist));
    if(dist<=fCut){
      Int_t dimclu=v[0];
      Bool_t already = kFALSE;
      for(Int_t k=1;k<=dimclu;k++){
	if(v[k]==j)already=kTRUE;
      }
      if(!already){
	++dimclu;
	v[0]=dimclu;
	fMark.SetBitNumber(j);
	fMark.SetBitNumber(i);
	v[dimclu]=j;
	Clustering(j,v);
      }
    }
  }
}



