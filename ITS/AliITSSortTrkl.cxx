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
fClustersTmp(NULL),
fClusters(NULL),
fNoClus(0),
fSize(NULL),
fMark(n),
fCut(cut),
fCoarseMaxRCut(0.) {
  // Standard constructor
  fPairs = new AliITSTracklPairs* [fkSize];
  fClustersTmp = new Int_t* [fkSize-1];
  for(Int_t i=0; i<fkSize-1; i++) fClustersTmp[i]=NULL;
}

//______________________________________________________________________
AliITSSortTrkl::AliITSSortTrkl(TClonesArray &tclo, Int_t n, Double_t cut, Double_t rcut):TObject(),
fkSize(n),
fIndex(0),
fPairs(),
fClustersTmp(NULL),
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
  for(Int_t i=0; i<fkSize-1; i++) fClustersTmp[i]=NULL;
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
    for(Int_t i=0;i<fIndex-1;i++){
      if(fClustersTmp[i]){
	delete []fClustersTmp[i];
      }
    }
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
  if(fIndex<2){
    AliWarning(Form("fIndex = %d",fIndex));
    fNoClus = 0;
    return fNoClus;
  }
  fMark.ResetAllBits();
  PrepareClustersTmp();
  //  cout<<"AliITSSortTrkl::FindClusters fkSize "<<fkSize<<"; fIndex "<<fIndex<<endl;
  for(Int_t i=0; i<fIndex-1;i++){
    Int_t *v=fClustersTmp[i];
    AliDebug(25,Form("Starting clustering for pair number %d",i));
    Clustering(i,v);  
    if(v[0]>0){
      AliDebug(25,Form("Clusters starting from pair %d : %d",i,v[0]));
      Int_t dim;
      v[v[0]+1]=i;
      // arr will contain the labels of the tracks associated to
      // the cluster of pairs starting from pair i.
      Int_t *arr=FindLabels(v+1,2*(v[0]+1),dim);

      /*
      cout<<"AliITSSortTrkl::FindClusters: Pairs involved \n";
      for(Int_t j=1;j<=v[0]+1;j++){
	cout<<v[j]<<" ";
	if(j%10==0)cout<<endl;
      }
      cout<<endl;
      cout<<"AliITSSortTrkl::FindClusters: Tracklets involved\n";
      for(Int_t j=0;j<dim;j++){
	cout<<arr[j]<<" ";
	if(j%10==0 && j>0)cout<<endl;
      }
      cout<<endl;
      */

      //In the following loop, all the pairs having 
      // one tracklet already associated to the cluster starting
      // at pair i are marked, in order to be excluded from further
      // searches    
      for(Int_t j=0;j<fIndex;j++){
	if(fMark.TestBitNumber(j))continue;
	for(Int_t k=0;k<dim;k++){
	  if(fPairs[j]->HasTrack(arr[k])){
	    fMark.SetBitNumber(j);
	    // cout<<"Marked pair "<<j<<" because has tracklet "<<arr[k]<<endl;
	    k=dim;
	  }
	}
      }
    
      delete []arr;
    }
  }
  // The following method builds the array fClusters
  Cleanup();
  return fNoClus;
}

//______________________________________________________________________
void AliITSSortTrkl::Cleanup(){
  // empty arrays are eliminated, the others are sorted according
  // to cluster multiplicity
  AliDebug(25,Form("fIndex = %d",fIndex));
  Int_t *siz = new Int_t[fIndex-1];
  Int_t *index = new Int_t[fIndex-1];
  fNoClus=0;
  for(Int_t i=0; i<fIndex-1;i++){
    Int_t *v=fClustersTmp[i];
    if(v[0]>0)fNoClus++;
    siz[i]=v[0];
  }
  if(fNoClus == 0){
    delete []siz;
    delete [] index;
    return;
  }
  AliDebug(25,Form("fNoClus = %d",fNoClus));
  TMath::Sort(fIndex-1,siz,index);
  fClusters = new Int_t* [fNoClus];
  fSize = new Int_t [fNoClus];
  for(Int_t i=0;i<fNoClus;i++){
    //    cout<<"AliITSSortTrkl::Cleanup: Cluster number "<<i<<"; Index= "<<index[i]<<endl;
    Int_t curind = index[i];
    Int_t *v=fClustersTmp[curind];
    fClusters[i] = new Int_t[v[0]+1];
    fSize[i]=v[0]+1;
    //    cout<<"AliITSSortTrkl::Cleanup: Size = "<<fSize[i]<<endl;
    Int_t *vext=fClusters[i];
    vext[0]=curind;
    for(Int_t j=1;j<fSize[i];j++)vext[j]=v[j];
    /*
    for(Int_t j=0;j<fSize[i];j++){
      cout<<vext[j]<<" ";
      if(j%10 == 0 && j!=0)cout<<endl;
    }
    cout<<endl;
    */
  }
  delete []siz;
  delete [] index;
}

//______________________________________________________________________
Int_t* AliITSSortTrkl::GetTrackletsLab(Int_t index, Int_t& dim) const {
  // Returns the tracklet labels corresponding to cluster index
  // Calling code must take care of memory deallocation
  AliDebug(25,Form("called with parameters %d %d",index,dim));
  if(fNoClus <=0){
    dim = 0;
    return NULL;
  }
  Int_t dimmax = 2*GetSizeOfCluster(index);
  AliDebug(25,Form("dimmax = %d",dimmax));
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

  //  cout<<"AliITSSortTrkl::Findlabels parameters "<<v<<" dimmax = "<<dimmax<<" dim "<<dim<<endl;
  Int_t *arr = new Int_t [dimmax];
  Int_t j=0;
  for(Int_t i=0; i<dimmax/2; i++){
    AliITSTracklPairs *pai=GetPairsAt(v[i]);
    arr[j++]=pai->GetTrack1();
    //    cout<<"AliITSSortTrkl::FindLabels - i="<<i<<" arr["<<j-1<<"]= "<<arr[j-1]<<endl;
    arr[j++]=pai->GetTrack2();
    //    cout<<"AliITSSortTrkl::FindLabels arr["<<j-1<<"]= "<<arr[j-1]<<endl;
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
   
  //  cout<<"AliITSSortTrkl::SortAndClean - Parameters: numb= "<<numb<<" numb2= "<<numb2<<endl;
  if(numb<=0)return;
  Int_t* index = new Int_t[numb];
  TMath::Sort(numb,arr,index,kFALSE);
  Int_t* tmp = new Int_t[numb];
  numb2 = 0;
  tmp[0] = arr[index[0]];
  for(Int_t i=1;i<numb;i++){
    if(arr[index[i]] != tmp[numb2]){
      ++numb2;
      tmp[numb2]=arr[index[i]];
    }
  }
  ++numb2;
  /*
  cout<<"AliITSSortTrkl::SortAndClean - numb2 = "<<numb2<<endl;
  for(Int_t i=0;i<numb;i++){
    if(i<numb2){
      arr[i]=tmp[i];
      cout<<"arr["<<i<<"]= "<<arr[i]<<endl;
    }
    else {
      arr[i]=0;
    }
  }
  */
  delete [] index;
  delete [] tmp;
}

//______________________________________________________________________
void AliITSSortTrkl::PrepareClustersTmp(){
  // prepare arrays of clusters
  for(Int_t i=0; i<fIndex-1;i++){
    fClustersTmp[i] = new Int_t [fIndex+1];
    Int_t *v = fClustersTmp[i];
    v[0]=0;
    for(Int_t j=1;j<fIndex+1;j++)v[j]=-1;
  }
}


//______________________________________________________________________
void AliITSSortTrkl::Clustering(Int_t i,Int_t *v){
  // recursive method to build up clusters starting from point i
  AliDebug(25,Form("Clustering called for pair %d",i));
  if(fMark.TestBitNumber(i)){
    AliDebug(25,Form("Leaving Clustering for pair %d - nothing done",i));
  return;
  }
  AliITSTracklPairs *p1 = fPairs[i];
  for(Int_t j=0;j<fIndex;j++){
    if(j == i) continue;
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
	fMark.SetBitNumber(i);
	AliDebug(25,Form("Marked pair %d - and call Clustering for pair %d",i,j));
	v[dimclu]=j;
	Clustering(j,v);
	fMark.SetBitNumber(j);
	AliDebug(25,Form("Marked pair %d",j));
      }
    }
  }
  AliDebug(25,Form("Leaving Clustering for pair %d ",i));
}





