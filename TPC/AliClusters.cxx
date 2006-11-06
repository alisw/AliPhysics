/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Time Projection Chamber clusters objects                                //
//
//  Origin: Marian Ivanov , GSI Darmstadt
//                                                                           //
//                                                                           //
//                                                                          //
///////////////////////////////////////////////////////////////////////////////
#include "TError.h"
#include "TClass.h"
#include  <TROOT.h>
#include "AliComplexCluster.h"
#include "AliClusters.h"
#include "TMarker.h"


const Int_t kDefSize = 1;  //defalut size


ClassImp(AliClusters) 

//*****************************************************************************
//
//_____________________________________________________________________________
AliClusters::AliClusters()
            :AliSegmentID(), 
	     fClusters(0),
             fNclusters(0),
             fClass(0)
{  
  //
  //default constructor
  //

}
//________________________________________________________________________
AliClusters::AliClusters(const AliClusters &param)
            :AliSegmentID(), 
	     fClusters(0),
             fNclusters(0),
             fClass(0)
{
  //
  //  copy constructor - dummy
  //
  fNclusters = param.fNclusters;
}
AliClusters & AliClusters::operator =(const AliClusters & param)
{
  //
  // assignment operator - dummy
  //
  fNclusters=param.fNclusters;
  return (*this);
}
//________________________________________________________________________
AliClusters::~AliClusters()
{
   //
   //default destructor
  //
   if (fClusters !=0) fClusters->Clear();
   delete fClusters;
}

//_________________________________________________________________________

Bool_t AliClusters::SetClass(const Text_t *classname)
{
  //
  //set class of stored object
  if ( fClass !=0 ) {
    //    delete fClass;
    fClass = 0;
  }
 
  if (!gROOT)
      ::Fatal("AliClusters::SetClass", "ROOT system not initialized");
   
   fClass = gROOT->GetClass(classname);
   if (!fClass) {
      Error("AliClusters", "%s is not a valid class name", classname);
      return kFALSE;
   }
   if (!fClass->InheritsFrom(TObject::Class())) {
      Error("AliClusters", "%s does not inherit from TObject", classname);
      return kFALSE; 
   } 
   return kTRUE;
}

//_____________________________________________________________________________
void AliClusters::SetArray(Int_t length)
{
  //
  // construct Clones array of object
  //
  if (fClusters!=0) delete fClusters;
  if (fClass==0){
     Error("AliClusters", "cluster type not initialised \n SetClass before!");
     return;
  }
  fClusters = new TClonesArray(fClass->GetName(),length);
}



//_____________________________________________________________________________
const  TObject* AliClusters::operator[](Int_t i)
{
  //
  // return cluster at internal position i
  //
  if (fClusters==0) return 0;
  return fClusters->UncheckedAt(i);
}
//_____________________________________________________________________________
void  AliClusters::Sort()
{
  // sort cluster 
  if (fClusters!=0) fClusters->Sort();
}

//_____________________________________________________________________________
TObject * AliClusters::InsertCluster( const TObject * c) 
{ 
  //
  // Add a simulated cluster copy to the list
  //
  if (fClass==0) { 
    Error("AliClusters", "class type not specified");
    return 0; 
  }
  if(!fClusters) fClusters=new TClonesArray(fClass->GetName(),1000);
  TClonesArray &lclusters = *fClusters;
  return new(lclusters[fNclusters++]) AliComplexCluster(*((AliComplexCluster*)c));
}

//_____________________________________________________________________________
Int_t AliClusters::Find(Double_t y) const 
{
  //
  // return index of cluster nearest to given y position
  //
  AliComplexCluster* cl;
  cl=(AliComplexCluster*)fClusters->UncheckedAt(0);
  if (y <= cl->GetY()) return 0;  
  cl=(AliComplexCluster*)fClusters->UncheckedAt(fNclusters-1);
  if (y > cl->GetY()) return fNclusters; 
  Int_t b=0, e=fNclusters-1, m=(b+e)/2;
  for (; b<e; m=(b+e)/2) {
    cl = (AliComplexCluster*)fClusters->UncheckedAt(m);
    if (y > cl->GetY()) b=m+1;
    else e=m; 
  }
  return m;
}


//_____________________________________________________________________________

void AliClusters::DrawClusters(Float_t shiftx, Float_t shifty, 
				  Int_t color, Int_t size, Int_t style)
{

  if (fClusters==0) return;  
  //draw marker for each of cluster
  Int_t ncl=fClusters->GetEntriesFast();
  for (Int_t i=0;i<ncl;i++){
    AliComplexCluster *cl = (AliComplexCluster*)fClusters->UncheckedAt(i);
    TMarker * marker = new TMarker;
    marker->SetX(cl->GetX()+shiftx);
    marker->SetY(cl->GetY()+shifty);
    marker->SetMarkerSize(size);
    marker->SetMarkerStyle(style);
    marker->SetMarkerColor(color);
    marker->Draw();    
  }
}

