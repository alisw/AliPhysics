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
// -Function FillArrayCoorAngles fills arrays wich contains, for each    //
//  layer, the global coordinates, errors on x,y,z and angles lambda     //
//  and phi for each cluster                                             //
///////////////////////////////////////////////////////////////////////////


#include <stdlib.h>
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
  fPhiList = 0; 
  fLambdaList = 0;
  fXList = 0;
  fYList =0;
  fZList =0;
  fSxList = 0;
  fSyList =0;
  fSzList =0;
  for(Int_t i=0;i<3;i++){fPrimaryVertex[i]=0;}
}

//__________________________________________________________
AliITSclusterTable::AliITSclusterTable(AliITSgeom* geom, AliITStrackerSA* tracker,Double_t* primaryVertex) {
// Standard constructor
  fDet=0;
  Int_t nm = geom->GetIndexMax();
  fNCl = new Int_t[nm];
  for(Int_t i=0;i<nm;i++){fNCl[i]=0;}
  fLbl=0;
  fGeom = geom;
  fTracker = tracker;
  for(Int_t i=0;i<3;i++){fPrimaryVertex[i]=primaryVertex[i];}
  fPhiList = 0;
  fLambdaList = 0;
  fXList = 0;
  fYList = 0;
  fZList =0;
  fSxList = 0;
  fSyList =0;
  fSzList =0;
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
    delete[] fDet;
  }
  if (fLbl)delete [] fLbl; // memory leak!
  if (fNCl)delete [] fNCl;
  for (Int_t i = 0; i < fGeom->GetNlayers(); i++) {
    if (fPhiList) delete fPhiList[i];
    if (fLambdaList) delete fLambdaList[i];
    if (fXList) delete fXList[i];
    if (fYList) delete fYList[i];
    if (fZList) delete fZList[i];
    if (fSxList) delete fSxList[i];
    if (fSyList) delete fSyList[i];
    if (fSzList) delete fSzList[i];
  }
  if (fPhiList) delete[] fPhiList;
  if (fLambdaList) delete[] fLambdaList;
  if (fXList) delete[] fXList;
  if (fYList) delete[] fYList;
  if (fZList) delete[] fZList;
  if (fSxList) delete[] fSxList;
  if (fSyList) delete[] fSyList;
  if (fSzList) delete[] fSzList;
}

//__________________________________________________________

void AliITSclusterTable::FillArray(TTree* clusterTree){
  
  //
  Int_t nm = fGeom->GetIndexMax();
  fDet = new TArrayI*[nm];

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
      exit(1);
    }
    else {
      for(Int_t n=0;n<vect[nlr]->GetSize();n++){
	Int_t mm=vect[nlr]->At(n);
	if(nc>=fDet[mod]->GetSize()) fDet[mod]->Set(nc*2+10);
	if(mm==mod) {(*fDet[mod])[nc]=n; nc+=1; }
      }
    }
  }

  clus->Delete();
  delete clus;
  for(Int_t n=0;n<fGeom->GetNlayers();n++)delete vect[n];
  delete [] vect;
  delete [] firstmod;
}

//_________________________________________________________________
void AliITSclusterTable::FillArrayLabel(Int_t numberofparticles){
  //

  fLbl = new TArrayI*[numberofparticles];
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

void AliITSclusterTable::FillArrayCoorAngles(){
  //Fill arrays with phi,lambda and indices of clusters for each layer
 
  fPhiList = new TArrayD*[fGeom->GetNlayers()];
  fLambdaList = new TArrayD*[fGeom->GetNlayers()];
  fXList = new TArrayF*[fGeom->GetNlayers()];
  fYList = new TArrayF*[fGeom->GetNlayers()];
  fZList = new TArrayF*[fGeom->GetNlayers()];
  fSxList = new TArrayF*[fGeom->GetNlayers()];
  fSyList = new TArrayF*[fGeom->GetNlayers()];
  fSzList = new TArrayF*[fGeom->GetNlayers()];

  Int_t * firstmod = new Int_t[fGeom->GetNlayers()+1];
  firstmod[fGeom->GetNlayers()]=fGeom->GetIndexMax();  // upper limit

  for(Int_t nlay=0;nlay<fGeom->GetNlayers();nlay++){
    firstmod[nlay] = fGeom->GetModuleIndex(nlay+1,1,1);
    Int_t ncl = fTracker->GetNumberOfClustersLayer(nlay);
    fPhiList[nlay] = new TArrayD(ncl);
    fLambdaList[nlay]=new TArrayD(ncl);
    fXList[nlay]=new TArrayF(ncl);
    fYList[nlay]=new TArrayF(ncl);
    fZList[nlay]=new TArrayF(ncl);
    fSxList[nlay]=new TArrayF(ncl);
    fSyList[nlay]=new TArrayF(ncl);
    fSzList[nlay]=new TArrayF(ncl);

    for(Int_t j=0;j<ncl;j++){
      AliITSclusterV2* cl = fTracker->GetClusterLayer(nlay,j);
      Double_t phi=0;Double_t lambda=0;
      Float_t x=0;Float_t y=0;Float_t z=0;
      Float_t sx=0;Float_t sy=0;Float_t sz=0;
      Int_t module = cl->GetDetectorIndex()+firstmod[nlay];
      GetCoorAngles(cl,module,phi,lambda,x,y,z);
      GetCoorErrors(cl,module,sx,sy,sz);
      fPhiList[nlay]->AddAt(phi,j);
      fLambdaList[nlay]->AddAt(lambda,j);
      fXList[nlay]->AddAt(x,j);
      fYList[nlay]->AddAt(y,j);
      fZList[nlay]->AddAt(z,j);
      fSxList[nlay]->AddAt(sx,j);
      fSyList[nlay]->AddAt(sy,j);
      fSzList[nlay]->AddAt(sz,j);
      
    }
    
  }
  
  delete [] firstmod;
}
void AliITSclusterTable::GetCoorAngles(AliITSclusterV2* cl,Int_t module,Double_t &phi,Double_t &lambda, Float_t &x, Float_t &y,Float_t &z){
  //Returns values of phi (azimuthal) and lambda angles for a given cluster
  
  Double_t rot[9];     fGeom->GetRotMatrix(module,rot);
  Int_t lay,lad,det; fGeom->GetModuleId(module,lay,lad,det);
  Float_t tx,ty,tz;  fGeom->GetTrans(lay,lad,det,tx,ty,tz);     

  Double_t alpha=TMath::ATan2(rot[1],rot[0])+TMath::Pi();
  Double_t phi1=TMath::Pi()/2+alpha;
  if (lay==1) phi1+=TMath::Pi();

  Float_t cp=TMath::Cos(phi1), sp=TMath::Sin(phi1);
  Float_t r=tx*cp+ty*sp;

  x= r*cp - cl->GetY()*sp;
  y= r*sp + cl->GetY()*cp;
  z=cl->GetZ();
  
  phi=TMath::ATan2(y,x);
  lambda=TMath::ATan2(z-fPrimaryVertex[2],TMath::Sqrt((x-fPrimaryVertex[0])*(x-fPrimaryVertex[0])+(y-fPrimaryVertex[1])*(y-fPrimaryVertex[1])));
}

void AliITSclusterTable::GetCoorErrors(AliITSclusterV2* cl, Int_t module,Float_t &sx,Float_t &sy, Float_t &sz){

  //returns x,y,z of cluster in global coordinates

  Double_t rot[9];     fGeom->GetRotMatrix(module,rot);
  Int_t lay,lad,det; fGeom->GetModuleId(module,lay,lad,det);
 
  Double_t alpha=TMath::ATan2(rot[1],rot[0])+TMath::Pi();
  Double_t phi=TMath::Pi()/2+alpha;
  if (lay==1) phi+=TMath::Pi();

  Float_t cp=TMath::Cos(phi), sp=TMath::Sin(phi);

  sx = TMath::Sqrt(sp*sp*cl->GetSigmaY2());
  sy = TMath::Sqrt(cp*cp*cl->GetSigmaY2());
  sz = TMath::Sqrt(cl->GetSigmaZ2());

}
