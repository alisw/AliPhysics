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


////////////////////////////////////////////////////////
//  Stand alone track class UPGRADE                   //
//  Authors: A.Mastroserio                            //
//           C.Terrevoli                              //
//           annalisa.mastroserio@cern.ch             //      
//           cristina.terrevoli@ba.infn.it            //                                                    
////////////////////////////////////////////////////////

#include "AliESDtrack.h"
#include "AliITStrackU.h"
#include "AliITSsegmentationUpgrade.cxx"

ClassImp(AliITStrackU)

//_____________________________________
  AliITStrackU:: AliITStrackU() : 
    AliITStrackV2(),
    fNLayers(0),
    fNU(0),
    fExpQ(40)
{
  // Default constructor  fgMaxNLayer
  SetNumberOfClusters(0);
  SetNumberOfClustersU(0);
  ResetIndexU();
  for(Int_t nlay=0;nlay<fgMaxNLayer;nlay++){ 
    fDy[nlay]=0; fDz[nlay]=0; fSigmaY[nlay]=0; fSigmaZ[nlay]=0; fSigmaYZ[nlay]=0;
    fClIndex[nlay]=-1; 
    fNy[nlay]=0; fNz[nlay]=0; fNormQ[nlay]=0; fNormChi2[nlay]=1000;
    SetNumberOfMarked(nlay,0);
  }
  ResetMarked();
}
//_____________________________________
AliITStrackU:: AliITStrackU(Int_t nlay) : 
  AliITStrackV2(),
  fNLayers(nlay),
  fNU(0),
  fExpQ(40)
{
  // Constructor
  SetNumberOfClusters(0);
  SetNumberOfClustersU(0);
  ResetIndexU();
  for(Int_t nl=0;nl<fgMaxNLayer;nl++){
    fDy[nl]=0; fDz[nl]=0; fSigmaY[nl]=0; fSigmaZ[nl]=0; fSigmaYZ[nl]=0;
    fClIndex[nl]=-1;
    fNy[nl]=0; fNz[nl]=0; fNormQ[nl]=0; fNormChi2[nl]=1000;
    SetNumberOfMarked(nl,0);
  }
  ResetMarked();
}

//_____________________________________________________
AliITStrackU::AliITStrackU(AliESDtrack& t,Bool_t c):
  AliITStrackV2(t,c),
  fNLayers(0),
  fNU(0),
  fExpQ(40)
{
  //------------------------------------------------------------------
  // Copy a V2 track into a U track
  // -> to be checked
  //------------------------------------------------------------------

  for(Int_t nlay=0;nlay<fgMaxNLayer;nlay++){
    fDy[nlay]=0; fDz[nlay]=0; fSigmaY[nlay]=0; fSigmaZ[nlay]=0; fSigmaYZ[nlay]=0;
    fClIndex[nlay]=-1; fNy[nlay]=0; fNz[nlay]=0; fNormQ[nlay]=0; fNormChi2[nlay]=1000;  
  }
}
//___________________________________________________

AliITStrackU::AliITStrackU(const AliITStrackU& t, Bool_t trackMI) : 
  AliITStrackV2(t),
  fNLayers(t.fNLayers),
  fNU(t.fNU),
  fExpQ(t.fExpQ)
{
  // Copy constructor

  ResetIndexU();
  ResetMarked();
  Int_t number = t.GetNumberOfClustersU();
  SetNumberOfClustersU(number);
  for(Int_t lay=0;lay<fgMaxNLayer;lay++){
    SetNumberOfMarked(lay,t.GetNumberOfMarked(lay));
  }
  for(Int_t i=0;i<number;i++){
    fSain[i]=t.fSain[i];
  }
  for(Int_t nlay=0;nlay<fNLayers;nlay++){
    for(Int_t i=0;i<t.GetNumberOfMarked(nlay);i++){
      fCluMark[nlay][i]=t.fCluMark[nlay][i];
    }
    fClIndex[nlay]= t.fClIndex[nlay]; fNy[nlay]=t.fNy[nlay]; fNz[nlay]=t.fNz[nlay]; fNormQ[nlay]=t.fNormQ[nlay]; fNormChi2[nlay] = t.fNormChi2[nlay];
  } 
  
  if(trackMI){
    fLab = t.fLab;
    fFakeRatio = t.fFakeRatio;
    for(Int_t i=0; i<fgMaxNLayer; i++) {fDy[i]=t.fDy[i]; fDz[i]=t.fDz[i];
      fSigmaY[i]=t.fSigmaY[i]; fSigmaZ[i]=t.fSigmaZ[i]; fSigmaYZ[i]=t.fSigmaYZ[i]; }
  }

}

//____________________________________________________
AliITStrackU::AliITStrackU(Double_t alpha, Double_t radius, Double_t Ycoor, Double_t Zcoor, Double_t phi, Double_t tanlambda, Double_t curv, Int_t lab, Int_t nlay ):
  fNLayers(nlay),
  fNU(0), 
  fExpQ(40)
{
  
  for(Int_t i=0; i<fNLayers; i++) { fClIndex[i]=-1; fNy[i]=0; fNz[i]=0; fNormQ[i]=0; fNormChi2[i]=1000; }
  // standard constructor. Used for ITSUpgrade standalone tracking

  // get the azimuthal angle of the detector containing the innermost
  // cluster of this track (data member fAlpha)
  for(Int_t i=0; i<fgMaxNLayer; i++) {fDy[i]=0; fDz[i]=0; fSigmaY[i]=0; fSigmaZ[i]=0; fSigmaYZ[i]=0;}

  if (alpha<0) alpha+=TMath::TwoPi();
  else if (alpha>=TMath::TwoPi()) alpha-=TMath::TwoPi();
  Init(alpha,radius,Ycoor,Zcoor,phi,tanlambda,curv,lab/*,nlay*/);
}
//____________________________________________________
void AliITStrackU::Init(Double_t alpha, Double_t radius, Double_t Ycoor, Double_t Zcoor, Double_t phi, Double_t tanlambda, Double_t curv, Int_t lab/*, Int_t nlay*/){
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

  for(Int_t i=0; i<fNLayers; i++) fIndex[i] = 0;  // to be set explicitely

  for(Int_t i=0; i<4; i++) fdEdxSample[i] = 0; 

  SetNumberOfClusters(0);
  SetNumberOfClustersU(0);
  for(Int_t nl=0;nl<fNLayers;nl++) SetNumberOfMarked(nl,0);
  ResetIndexU();
  ResetMarked();
  SetChi2(0);
  SetMass(0.139);    // pion mass
  SetLabel(lab); 
  
}

//____________________________________________________________
void AliITStrackU::AddClusterU(Int_t layer, Int_t clnumb) {
  // add one clusters to the list (maximum number=kMaxNumberOfClusters)
  Int_t presnum = GetNumberOfClustersU();
  if(presnum>=kMaxNumberOfClusters){
    Warning("AddClusterU","Maximum number of clusters already reached. Nothing is done\n");
    return;
  }

  fSain[presnum] = (layer<<28)+clnumb;  
  presnum++;
  SetNumberOfClustersU(presnum);
}

//____________________________________________________________
void AliITStrackU::AddClusterMark(Int_t layer, Int_t clnumb) {
  // add one clusters to the list (maximum number=kMaxNumberOfClusters)
  Int_t presnum = GetNumberOfMarked(layer);
  //printf("presnum=%d\n",presnum);
  if(presnum>=kMaxNumberOfClustersL){
    Warning("AddClusterMark","Maximum number of clusters already reached. Nothing is done\n");
    return;
  }

  fCluMark[layer][presnum] = clnumb;  
  presnum++;
  SetNumberOfMarked(layer,presnum);
}

//____________________________________________________________
void AliITStrackU::AddClusterV2(Int_t layer,Int_t clnumb) {
  // add one clusters to the list (maximum number=6)
  Int_t presnum = GetNumberOfClusters();
  if(presnum>=fNLayers){
    Warning("AddClusterV2","Maximum number of clusters already reached. Nothing is done\n");
    return;
  }    

  fIndex[presnum] = (layer<<28)+clnumb;  
  presnum++;
  SetNumberOfClusters(presnum);
}

//_____________________________________________________________
void AliITStrackU::ResetMarked(){

  //Reset array of marked clusters
  for(Int_t nlay=0;nlay<fNLayers;nlay++){
    for(Int_t k=0; k<kMaxNumberOfClustersL; k++) fCluMark[nlay][k]=0;
  }
}
//_____________________________________________________________
Double_t AliITStrackU::GetPredictedChi2MI(Double_t cy, Double_t cz, Double_t cerry, Double_t cerrz, Double_t covyz) const
{
  //-----------------------------------------------------------------
  // This function calculates a predicted chi2 increment.
  //-----------------------------------------------------------------
  Double_t p[2]={cy, cz};
  Double_t cov[3]={cerry*cerry, covyz, cerrz*cerrz};
  return AliExternalTrackParam::GetPredictedChi2(p,cov);
}

//____________________________________________________________________________
Bool_t AliITStrackU::UpdateMI(const AliCluster *c, Double_t chi2, Int_t index) {
  //------------------------------------------------------------------
  //This function updates track parameters
  //------------------------------------------------------------------
  Double_t dy=c->GetY() - GetY(), dz=c->GetZ() - GetZ();
  Int_t layer = (index & 0xf0000000) >> 28;
  fDy[layer] = dy;
  fDz[layer] = dz;
  fSigmaY[layer] = TMath::Sqrt(c->GetSigmaY2()+GetSigmaY2());
  fSigmaZ[layer] = TMath::Sqrt(c->GetSigmaZ2()+GetSigmaZ2());
  fSigmaYZ[layer] = c->GetSigmaYZ()+GetSigmaZY();


  return Update(c,chi2,index);
}
 
