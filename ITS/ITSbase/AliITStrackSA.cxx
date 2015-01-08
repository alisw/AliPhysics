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
//  of cluster equal to kMaxNumberOfClusters    //
////////////////////////////////////////////////////

#include "AliITSgeomTGeo.h"
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
  for(Int_t nlay=0;nlay<AliITSgeomTGeo::GetNLayers();nlay++){ 
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
  for(Int_t nlay=0;nlay<AliITSgeomTGeo::GetNLayers();nlay++){ 
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
  for(Int_t nlay=0;nlay<AliITSgeomTGeo::GetNLayers();nlay++){
    SetNumberOfMarked(nlay,t.GetNumberOfMarked(nlay));
  }
  for(Int_t i=0;i<number;i++){
    fSain[i]=t.fSain[i];
  }
  for(Int_t nlay=0;nlay<AliITSgeomTGeo::GetNLayers();nlay++){
    for(Int_t i=0;i<t.GetNumberOfMarked(nlay);i++){
      fCluMark[nlay][i]=t.fCluMark[nlay][i];
    }
  }
}
//____________________________________________________
AliITStrackSA::AliITStrackSA(Int_t layer, Int_t ladder, Int_t detector, Double_t Ycoor, Double_t Zcoor, Double_t phi, Double_t tanlambda, Double_t curv, Int_t lab ):
fNSA(0) 
{
  // standard constructor. Used for ITS standalone tracking

  // get the azimuthal angle of the detector containing the innermost
  // cluster of this track (data member fAlpha)

  TGeoHMatrix m; AliITSgeomTGeo::GetOrigMatrix(layer,ladder,detector,m);
  const TGeoHMatrix *tm=AliITSgeomTGeo::GetTracking2LocalMatrix(layer,ladder,detector);
  m.Multiply(tm);
  Double_t txyz[3]={0.}, xyz[3]={0.};
  m.LocalToMaster(txyz,xyz);
  Double_t sAlpha=TMath::ATan2(xyz[1],xyz[0]);

  if (sAlpha<0) sAlpha+=TMath::TwoPi();
  else if (sAlpha>=TMath::TwoPi()) sAlpha-=TMath::TwoPi();

  Double_t sX=TMath::Sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]);

  Init(sAlpha,sX,Ycoor,Zcoor,phi,tanlambda,curv,lab);

}
//____________________________________________________
AliITStrackSA::AliITStrackSA(Double_t alpha, Double_t radius, Double_t Ycoor, Double_t Zcoor, Double_t phi, Double_t tanlambda, Double_t curv, Int_t lab ):
fNSA(0) 
{
  // standard constructor. Used for ITS standalone tracking

  // get the azimuthal angle of the detector containing the innermost
  // cluster of this track (data member fAlpha)

  if (alpha<0) alpha+=TMath::TwoPi();
  else if (alpha>=TMath::TwoPi()) alpha-=TMath::TwoPi();
  Init(alpha,radius,Ycoor,Zcoor,phi,tanlambda,curv,lab);
}
//____________________________________________________
  void AliITStrackSA::Init(Double_t alpha, Double_t radius, Double_t Ycoor, Double_t Zcoor, Double_t phi, Double_t tanlambda, Double_t curv, Int_t lab ){
    // initialize parameters

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
                   TMath::Sin(phi-alpha),
		   tanlambda,
		   curv/conv};


  // dealing with the case B=0 (taken from AliTPCtrack.cxx)
  Double_t mostProbablePt=AliExternalTrackParam::GetMostProbablePt();
  Double_t p0=TMath::Sign(1/mostProbablePt,sP[4]);
  Double_t w0=sC[14]/(sC[14] + p0*p0), w1=p0*p0/(sC[14] + p0*p0);
  sP[4] = w0*p0 + w1*sP[4];
  sC[14]*=w1;
                                                                              
  Set(radius,alpha,sP,sC);

  for(Int_t i=0; i<AliITSgeomTGeo::GetNLayers(); i++) fIndex[i] = 0;  // to be set explicitely

  for(Int_t i=0; i<4; i++) fdEdxSample[i] = 0; 

  SetNumberOfClusters(0);
  SetNumberOfClustersSA(0);
  for(Int_t nlay=0;nlay<AliITSgeomTGeo::GetNLayers();nlay++) SetNumberOfMarked(nlay,0);
  ResetIndexSA();
  ResetMarked();
  SetChi2(0);
  SetMass(0.139);    // pion mass
  SetLabel(lab); 
  
}

//____________________________________________________________
void AliITStrackSA::AddClusterSA(Int_t layer, Int_t clnumb) {
  // add one clusters to the list (maximum number=kMaxNumberOfClusters)
  Int_t presnum = GetNumberOfClustersSA();
  if(presnum>=kMaxNumberOfClusters){
    Warning("AddClusterSA","Maximum number of clusters already reached. Nothing is done\n");
    return;
  }

  fSain[presnum] = (layer<<28)+clnumb;  
  presnum++;
  SetNumberOfClustersSA(presnum);
}

//____________________________________________________________
void AliITStrackSA::AddClusterMark(Int_t layer, Int_t clnumb) {
  // add one clusters to the list (maximum number=kMaxNumberOfClusters)
  Int_t presnum = GetNumberOfMarked(layer);
  //  printf("presnum=%d\n",presnum);
  if(presnum>=kMaxNumberOfClustersL){
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
  if(presnum>=AliITSgeomTGeo::GetNLayers()){
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
  for(Int_t nlay=0;nlay<AliITSgeomTGeo::GetNLayers();nlay++){
    for(Int_t k=0; k<kMaxNumberOfClustersL; k++) fCluMark[nlay][k]=0;
  }
}

