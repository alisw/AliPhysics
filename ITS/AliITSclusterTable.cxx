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

//////////////////////////////////////////////////////////////////////////
// Class used to simplify some operations with clusters.                 //
// -Function FillArray fills an array wich contains, for each            //
//  ITS module, an array with the indices of all the clusters detected   //
//  by the module. The indices correspond to the cluster indices in class// 
//  AliITSlayer of AliITStrackerV2.                                      //
//  This function is used in AliITStrackerSA::FindTracks.                // 
// -Function FillArrayLabel fills an array wich contains, for each       //
//  particle label, and for each layer, the information on clusters:     //
//  0 if there is no cluster, 1 if there is a cluster with this label.   //
//  This function is used to define trackable tracks.                    //   
///////////////////////////////////////////////////////////////////////////


#include "AliITSclusterTable.h"
#include "AliITSclusterV2.h"
#include "AliITSgeom.h"
#include "AliITStrackerSA.h"
#include<TBranch.h>
#include<TClonesArray.h>
#include<TTree.h>
ClassImp(AliITSclusterTable)

//__________________________________________________________
AliITSclusterTable::AliITSclusterTable(){
// Default constructor
  fDet=0;
  fNCl = 0;
  fLbl=0;
  fGeom = 0;
  fTracker = 0;
}

//__________________________________________________________
AliITSclusterTable::AliITSclusterTable(AliITSgeom* geom, AliITStrackerSA* tracker) {
// Standard constructor
  fDet=0;
  Int_t nm = geom->GetIndexMax();
  fNCl = new Int_t[nm];
  for(Int_t i=0;i<nm;i++){fNCl[i]=0;}
  fLbl=0;
  fGeom = geom;
  fTracker = tracker;
}

//______________________________________________________________________
AliITSclusterTable::AliITSclusterTable(const AliITSclusterTable &tab) : 
                    TObject(tab) {
  // Copy constructor
  // Copies are not allowed. The method is protected to avoid misuse.
  Error("AliITSclusterTable","Copy constructor not allowed\n");
}

//______________________________________________________________________
AliITSclusterTable& AliITSclusterTable::operator=(const 
                    AliITSclusterTable& /* tab */){
  // Assignment operator
  // Assignment is not allowed. The method is protected to avoid misuse.
  Error("= operator","Assignment operator not allowed\n");
  return *this;
}

//__________________________________________________________
AliITSclusterTable::~AliITSclusterTable(){
// Destructor
  Int_t nm = fGeom->GetIndexMax();
  if(fDet){
    for (Int_t i=0; i<nm; i++){
      delete fDet[i];
    }
    delete fDet;
  }
  /*  if(fLbl){
    for (Int_t i=0; i<nm; i++){
      delete fLbl[i];
    }
    delete fLbl;
    }*/
  if (fLbl)delete [] fLbl;
  if (fNCl)delete [] fNCl;
}

//__________________________________________________________

void AliITSclusterTable::FillArray(TTree* clusterTree,Int_t evnumber){
  
  //
  Int_t nm = fGeom->GetIndexMax();
  fDet = new TArrayI*[nm];
  fTracker->SetEventNumber(evnumber);
  fTracker->LoadClusters(clusterTree);

  TArrayI** vect = new TArrayI*[fGeom->GetNlayers()];
  
  Int_t * firstmod = new Int_t[fGeom->GetNlayers()+1];
  firstmod[fGeom->GetNlayers()]=fGeom->GetIndexMax();  // upper limit
  for(Int_t nlayer=0;nlayer<fGeom->GetNlayers();nlayer++){
    firstmod[nlayer] = fGeom->GetModuleIndex(nlayer+1,1,1);
    Int_t ncl = fTracker->GetNumberOfClustersLayer(nlayer);
    vect[nlayer]=new TArrayI(ncl); 
    for(Int_t j=0;j<ncl;j++){
      AliITSclusterV2* cl = fTracker->GetClusterLayer(nlayer,j);
      vect[nlayer]->AddAt(cl->GetDetectorIndex()+firstmod[nlayer],j);
    }
  }
    
  TBranch *brancht=(TBranch*)clusterTree->GetBranch("Clusters");
  if(!brancht) Warning("FillArray","No cluster branch");
  TClonesArray* clus = new TClonesArray("AliITSclusterV2",10000);
  brancht->SetAddress(&clus);
 
  
  for(Int_t mod=0;mod<nm;mod++){
    Int_t nc=0;
    clusterTree->GetEvent(mod);
    Int_t ncl = clus->GetEntries();
    fDet[mod] = new TArrayI(ncl);
    fNCl[mod]= ncl;
    Int_t nlr = FindIndex(fGeom->GetNlayers(),firstmod,mod);
    if(nlr<0){
      Fatal("FillArray","Wrong module number %d . Limits: %d , %d",mod,firstmod[0],firstmod[fGeom->GetNlayers()+1]);
      return;
    }
    else {
      for(Int_t n=0;n<vect[nlr]->GetSize();n++){
	Int_t mm=vect[nlr]->At(n);
	if(mm==mod) {fDet[mod]->AddAt(n,nc); nc+=1; }
      }
    }
  }


  fTracker->UnloadClusters();

  for(Int_t n=0;n<fGeom->GetNlayers();n++)delete vect[n];
  delete vect;
  delete [] firstmod;
}

//_________________________________________________________________
void AliITSclusterTable::FillArrayLabel(Int_t numberofparticles,TTree* clusterTree,Int_t evnumber){
  //


  fLbl = new TArrayI*[numberofparticles];
  fTracker->SetEventNumber(evnumber);
  fTracker->LoadClusters(clusterTree);
  const Int_t knm =fGeom->GetNlayers();
  for(Int_t nlab=0;nlab<numberofparticles;nlab++){
    fLbl[nlab] = new TArrayI(knm);
    Int_t * nn = new Int_t[knm]; 
    for(Int_t i=0;i<knm;i++)nn[i]=0;
    for(Int_t nlayer=0;nlayer<knm;nlayer++){
      Int_t ncl = fTracker->GetNumberOfClustersLayer(nlayer);
      while(ncl--){
	AliITSclusterV2* cl = fTracker->GetClusterLayer(nlayer,ncl);
	if(cl->IsUsed()==1) continue;
	if(cl->GetLabel(0)==nlab || cl->GetLabel(1)==nlab || cl->GetLabel(2)==nlab){
	  nn[nlayer]+=1;
	  cl->Use();
	  break;
	}
      }     
      fLbl[nlab]->AddAt(nn[nlayer],nlayer);
    }
    delete [] nn;
  }

  fTracker->UnloadClusters();
}

//_______________________________________________________________
Int_t AliITSclusterTable::ThisParticleIsTrackable(Int_t label,Int_t numberofpoints){

  //Returns 1 if particle with label "label" is trackable. 
   
  Int_t nb=0;
  for(Int_t i=0;i<fGeom->GetNlayers();i++){
    Int_t ncl = fLbl[label]->At(i);
    if(ncl>0) nb++;
  }
  if(nb>=numberofpoints) return 1;
  else return 0;

}

//_______________________________________________________________
Int_t AliITSclusterTable::FindIndex(Int_t ndim, Int_t *ptr, Int_t value){
// ptr[ndim+1] is an array of integers
// ndim is its dimension minus 1
// value is a number such as:  ptr[0]<=value < ptr[ndim]
// if ptr[i]<=value<ptr[i+1] the index i is returned;  -1 if out of bounds
  Int_t retval = -1;
  for(Int_t k=0;k<ndim;k++){
    if(value>=ptr[k] && value <ptr[k+1]){
      retval = k;
      break;
    }
  }
  return retval;
}
