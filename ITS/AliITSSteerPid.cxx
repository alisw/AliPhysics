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

/////////////////////////////////////////////////////////////////////////
//Steering Class for PID in the ITS                                      //
//The PID is based on the likelihood of all the four ITS' layers,        //
//without using the truncated mean for the dE/dx. The response           //
//functions for each layer are convoluted Landau-Gaussian functions.     // 
// Origin: Elena Bruna bruna@to.infn.it, Massimo Masera masera@to.infn.it//
/////////////////////////////////////////////////////////////////////////

#include <TMath.h>

#include "AliITSSteerPid.h"

ClassImp(AliITSSteerPid)

  //______________________________________________________________
AliITSSteerPid::AliITSSteerPid():
fClonarr2(0),
fVect2(0),
fVect2lay1(0),
fVect2lay2(0),
fVect2lay3(0),
fVect2lay4(0),
fFitTree(0),
fItem(0),
fPCenter(0),
fPWidth(0)
{
  // default constructor
 }
//______________________________________________________________
AliITSSteerPid::~AliITSSteerPid(){
  // destructor
  delete fClonarr2;
  delete fFitTree;
}

//______________________________________________________________________
AliITSSteerPid::AliITSSteerPid(const AliITSSteerPid &ob) :TObject(ob),
fClonarr2(ob.fClonarr2),
fVect2(ob.fVect2),
fVect2lay1(ob.fVect2lay1),
fVect2lay2(ob.fVect2lay2),
fVect2lay3(ob.fVect2lay3),
fVect2lay4(ob.fVect2lay4),
fFitTree(ob.fFitTree),
fItem(ob.fItem),
fPCenter(ob.fPCenter),
fPWidth(ob.fPWidth) {
  // Copy constructor
}

//______________________________________________________________________
AliITSSteerPid& AliITSSteerPid::operator=(const AliITSSteerPid& ob){
  // Assignment operator
  this->~AliITSSteerPid();
  new(this) AliITSSteerPid(ob);
  return *this;
}

//______________________________________________________________
void AliITSSteerPid::InitLayer(TString fileITS,TString fileFitPar){
  // it opens the files useful for the PID 
  TFile *fClonarr2=new TFile (fileITS,"r");
  fVect2=(TClonesArray*)fClonarr2->Get("vectfitits_0");//truncated mean
  fVect2lay1=(TClonesArray*)fClonarr2->Get("vectfitits_1");//lay 1
  fVect2lay2=(TClonesArray*)fClonarr2->Get("vectfitits_2");//lay 2
  fVect2lay3=(TClonesArray*)fClonarr2->Get("vectfitits_3");//lay 3
  fVect2lay4=(TClonesArray*)fClonarr2->Get("vectfitits_4");//lay 4
 
  TFile *fFitPar=new TFile (fileFitPar);
  fFitTree=(TTree*)fFitPar->Get("tree");
 
}

//______________________________________________________________
AliITSPidParItem* AliITSSteerPid::GetItemLayer(Int_t nolay,Float_t mom){
    // it gives an AliITSPidParItem object for a given momentum and ITS layer
  if(nolay==1) return Item(fVect2lay1,mom);
  if(nolay==2) return Item(fVect2lay2,mom);
  if(nolay==3) return Item(fVect2lay3,mom);
  if(nolay==4) return Item(fVect2lay4,mom);
  if(nolay!=1&&nolay!=2&&nolay!=3&&nolay!=4) {
    fItem=new AliITSPidParItem();
    return fItem;

  }
  return 0;
}

//______________________________________________________________
void AliITSSteerPid::GetParFitLayer(Int_t nolay,Float_t mom,Double_t *parp,Double_t *park,Double_t *parpi){
  //it gives the parameters of the convoluted functions (WL, MP, WG) for 
  //protons, kaons and pions for a given momentum and ITS layer

  Double_t parfit0pro[3]={0,0,0};
  Double_t parfit1pro[3]={0,0,0};
  Double_t parfit3pro[3]={0,0,0};
  Double_t parfit0kao[3]={0,0,0};
  Double_t parfit1kao[3]={0,0,0};
  Double_t parfit3kao[3]={0,0,0};
  Double_t parfit0pi[3]={0,0,0};
  Double_t parfit1pi[3]={0,0,0};
  Double_t parfit3pi[3]={0,0,0};
 
  fFitTree->SetBranchAddress("par0pro",parfit0pro);
  fFitTree->SetBranchAddress("par1pro",parfit1pro);
  fFitTree->SetBranchAddress("par3pro",parfit3pro);

  fFitTree->SetBranchAddress("par0kao",parfit0kao);
  fFitTree->SetBranchAddress("par1kao",parfit1kao);
  fFitTree->SetBranchAddress("par3kao",parfit3kao);

  fFitTree->SetBranchAddress("par0pi",parfit0pi);
  fFitTree->SetBranchAddress("par1pi",parfit1pi);
  fFitTree->SetBranchAddress("par3pi",parfit3pi);
  fFitTree->GetEvent(nolay);
 
  GetLangausProPars(mom,parfit0pro,parfit1pro,parfit3pro,parp);
  GetLangausKaoPars(mom,parfit0kao,parfit1kao,parfit3kao,park);
  GetLangausPiPars(mom,parfit0pi,parfit1pi,parfit3pi,parpi);


}//______________________________________________________________
void AliITSSteerPid::GetLangausProPars(Float_t mom,Double_t *parfit0,Double_t *parfit1,Double_t *parfit3,Double_t *par){
 
  //It finds the parameters of the convoluted Landau-Gaussian response 
  //function for protons (Width Landau, Most Probable, Width Gaussian)
  par[0]=parfit0[0]+parfit0[1]/mom;
  par[1]=parfit1[0]/(mom*mom)+parfit1[1]/(mom*mom)*TMath::Log(mom*mom)+parfit1[2];
  par[2]=parfit3[0]/(mom*mom)+parfit3[1]/(mom*mom)*TMath::Log(mom*mom)+parfit3[2];
}
//______________________________________________________________
void AliITSSteerPid::GetLangausKaoPars(Float_t mom,Double_t *parfit0,Double_t *parfit1,Double_t *parfit3,Double_t *par){
  // It finds the parameters of the convoluted Landau-Gaussian response 
  //function for kaons (Width Landau, Most Probable, Width Gaussian)

  par[0]=parfit0[0]+parfit0[1]/(mom*mom);
  par[1]=parfit1[0]/(mom*mom)+parfit1[1]/(mom*mom)*TMath::Log(mom*mom)+parfit1[2];
  par[2]=parfit3[0]/(mom*mom)+parfit3[1]/(mom*mom)*TMath::Log(mom*mom)+parfit3[2];
}

//______________________________________________________________
void AliITSSteerPid::GetLangausPiPars(Float_t mom,Double_t *parfit0,Double_t *parfit1,Double_t *parfit3,Double_t *par){
  //It finds the parameters of the convoluted Landau-Gaussian response 
  //function for pions (Width Landau, Most Probable, Width Gaussian)

  par[0]=parfit0[0]/(mom*mom)+parfit0[1]/(mom*mom)*TMath::Log(mom*mom)+parfit0[2];
  par[1]=parfit1[0]/(mom)+parfit1[1]/(mom)*TMath::Log(mom*mom)+parfit1[2];
  par[2]=parfit3[0]/(mom*mom)+parfit3[1]/(mom*mom)*TMath::Log(mom*mom)+parfit3[2];
}



//______________________________________________________________
AliITSPidParItem* AliITSSteerPid::Item(TClonesArray *Vect,Float_t mom){

  //it gives an AliITSPidParItem object taken from the TClonesArray. 
  Int_t mybin=-1;

  AliITSPidParItem* punt;
  
  for (Int_t a=0;a<50;a++){
    punt=(AliITSPidParItem*)Vect->At(a);
    Float_t centerp=punt->GetMomentumCenter(); 
    Float_t widthp=punt->GetWidthMom();
    if (mom>(centerp-widthp/2) && mom<=(centerp+widthp/2)) mybin=a; 
  }
  if (mybin!=-1) fItem=(AliITSPidParItem*)Vect->At(mybin);
  else {
    fPCenter=0;
    fPWidth=0;
    for (Int_t ii=0;ii<52;ii++) fBuff[ii]=0;
    fItem = new AliITSPidParItem(fPCenter,fPWidth,fBuff);
  }
  
  return fItem;



}
