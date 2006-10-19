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

/* $Id$ */

////////////////////////////////////////////////////
//  Stand alone track class                       //
//  Origin:  Elisabetta Crescio                   //
//  e-mail:  crescio@to.infn.it                   //
//  it is a V2 track with a possible number       //
//  of cluster equal to fgkMaxNumberOfClusters    //
////////////////////////////////////////////////////

#include "AliITSgeom.h"
#include "AliITStrackSA.h"


ClassImp(AliITStrackSA)

//_____________________________________
AliITStrackSA:: AliITStrackSA() : AliITStrackMI(),
fNSA(0)
{
// Default constructor  
  SetNumberOfClusters(0);
  SetNumberOfClustersSA(0);
  ResetIndexSA();
  for(Int_t nlay=0;nlay<fgkLayers;nlay++){ 
    SetNumberOfMarked(nlay,0);
  }
  ResetMarked();
}


//___________________________________________________
AliITStrackSA::AliITStrackSA(const AliITStrackMI& t) : 
AliITStrackMI(t),
fNSA(0){
// Copy a V2 track into a SA track
  SetNumberOfClustersSA(0);
  ResetIndexSA();
  for(Int_t nlay=0;nlay<fgkLayers;nlay++){ 
    SetNumberOfMarked(nlay,0);
  }
  ResetMarked();

}
//___________________________________________________
AliITStrackSA::AliITStrackSA(const AliITStrackSA& t) : 
AliITStrackMI(t),
fNSA(t.fNSA){
// Copy constructor


  ResetIndexSA();
  ResetMarked();
  Int_t number = t.GetNumberOfClustersSA();
  SetNumberOfClustersSA(number);
  for(Int_t nlay=0;nlay<fgkLayers;nlay++){
    SetNumberOfMarked(nlay,t.GetNumberOfMarked(nlay));
  }
  for(Int_t i=0;i<number;i++){
    fSain[i]=t.fSain[i];
  }
  for(Int_t nlay=0;nlay<fgkLayers;nlay++){
    for(Int_t i=0;i<t.GetNumberOfMarked(nlay);i++){
      fCluMark[nlay][i]=t.fCluMark[nlay][i];
    }
  }
}
//____________________________________________________
AliITStrackSA::AliITStrackSA(AliITSgeom* geom,Int_t layer, Int_t ladder, Int_t detector, Double_t Ycoor, Double_t Zcoor, Double_t phi, Double_t tanlambda, Double_t curv, Int_t lab ):
fNSA(0) 
{
  // standard constructor. Used for ITS standalone tracking

  if(!geom){
    Fatal("AliITStrackSA","ITS geometry not found - Abort\n");
    return;
  }
  // get the azimuthal angle of the detector containing the innermost
  // cluster of this track (data member fAlpha)
  Float_t rotmatr[9];
  geom->GetRotMatrix(layer,ladder,detector,rotmatr);
  Double_t sAlpha=TMath::ATan2(rotmatr[1],rotmatr[0])+TMath::Pi();
  sAlpha+=TMath::Pi()/2.;
  if(layer==1) sAlpha+=TMath::Pi();


  // get the radius of this detector. Procedure taken from the 
  // AliITStrackerV2 constructor
  Float_t x=0,y=0,z=0;
  geom->GetTrans(layer,ladder,detector,x,y,z);

  Double_t fi=TMath::ATan2(rotmatr[1],rotmatr[0])+TMath::Pi();
  fi+=TMath::Pi()/2;
  if (layer==1) fi+=TMath::Pi();
  Double_t cp=TMath::Cos(fi), sp=TMath::Sin(fi);
  Double_t sX=x*cp+y*sp;


  fdEdx = 0;

  Double_t conv=GetBz()*kB2C;
  Double_t sC[] = {0.000009, // 0.000009
                   0.,
                   0.000003, //0.000030
                   0.,
	           0.,
	           0.000001, //0.000001
	           0.,
	           0.,
	           0.,
	           0.000002, //0.000002
	           0.,
	           0.,
	           0.,
	           0.,
		   0.000001/(conv*conv)}; //0.0000001

  Double_t sP[] = {Ycoor,
		   Zcoor,
                   TMath::Sin(phi-sAlpha),
		   tanlambda,
		   curv/conv};

  Set(sX,sAlpha,sP,sC);

  for(Int_t i=0; i<kMaxLayer; i++) fIndex[i] = 0;  // to be set explicitely

  for(Int_t i=0; i<4; i++) fdEdxSample[i] = 0; 

  SetNumberOfClusters(0);
  SetNumberOfClustersSA(0);
  for(Int_t nlay=0;nlay<fgkLayers;nlay++) SetNumberOfMarked(nlay,0);
  ResetIndexSA();
  ResetMarked();
  SetChi2(0);
  SetMass(0.139);    // pion mass
  SetLabel(lab); 
  
}

//____________________________________________________________
void AliITStrackSA::AddClusterSA(Int_t layer, Int_t clnumb) {
  // add one clusters to the list (maximum number=fgkMaxNumberOfClusters)
  Int_t presnum = GetNumberOfClustersSA();
  if(presnum>=fgkMaxNumberOfClusters){
    Warning("AddClusterSA","Maximum number of clusters already reached. Nothing is done\n");
    return;
  }

  fSain[presnum] = (layer<<28)+clnumb;  
  presnum++;
  SetNumberOfClustersSA(presnum);
}

//____________________________________________________________
void AliITStrackSA::AddClusterMark(Int_t layer, Int_t clnumb) {
  // add one clusters to the list (maximum number=fgkMaxNumberOfClusters)
  Int_t presnum = GetNumberOfMarked(layer);
  //  printf("presnum=%d\n",presnum);
  if(presnum>=fgkMaxNumberOfClustersL){
    Warning("AddClusterMark","Maximum number of clusters already reached. Nothing is done\n");
    return;
  }

  fCluMark[layer][presnum] = clnumb;  
  presnum++;
  SetNumberOfMarked(layer,presnum);
}

//____________________________________________________________
void AliITStrackSA::AddClusterV2(Int_t layer,Int_t clnumb) {
  // add one clusters to the list (maximum number=6)
  Int_t presnum = GetNumberOfClusters();
  if(presnum>=kMaxLayer){
    Warning("AddClusterV2","Maximum number of clusters already reached. Nothing is done\n");
    return;
   }    

  fIndex[presnum] = (layer<<28)+clnumb;  
  presnum++;
  SetNumberOfClusters(presnum);
}

//_____________________________________________________________
void AliITStrackSA::ResetMarked(){

  //Reset array of marked clusters
  for(Int_t nlay=0;nlay<fgkLayers;nlay++){
    for(Int_t k=0; k<fgkMaxNumberOfClustersL; k++) fCluMark[nlay][k]=0;
  }
}














