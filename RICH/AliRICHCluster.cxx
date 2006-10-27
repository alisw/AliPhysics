//  **************************************************************************
//  * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
//  *                                                                        *
//  * Author: The ALICE Off-line Project.                                    *
//  * Contributors are mentioned in the code where appropriate.              *
//  *                                                                        *
//  * Permission to use, copy, modify and distribute this software and its   *
//  * documentation strictly for non-commercial purposes is hereby granted   *
//  * without fee, provided that the above copyright notice appears in all   *
//  * copies and that both the copyright notice and this permission notice   *
//  * appear in the supporting documentation. The authors make no claims     *
//  * about the suitability of this software for any purpose. It is          *
//  * provided "as is" without express or implied warranty.                  *
//  **************************************************************************

#include "AliRICHCluster.h"  //class header
#include <TMinuit.h>         //Solve()
#include <TClonesArray.h>    //Solve()

ClassImp(AliRICHCluster)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliRICHCluster::CoG()
{
// Calculates naive cluster position as a center of gravity of its digits.
// Arguments: none 
//   Returns: shape of the cluster i.e. the box which fully contains the cluster      
  if(fDigs==0) return;                                  //no digits in this cluster
  fX=fY=0;                                              //set cluster position to (0,0) to start to collect contributions
  for(Int_t iDig=0;iDig<fDigs->GetEntriesFast();iDig++){//digits loop
    AliRICHDigit *pDig=(AliRICHDigit*)fDigs->At(iDig);  //get pointer to next digit
    Float_t q=pDig->Q();                                //get QDC 
    fX += pDig->LorsX()*q;fY +=pDig->LorsY()*q;         //add digit center weighted by QDC
  }//digits loop
  fX/=fQ;fY/=fQ;                                        //final center of gravity

  CorrSin();
  
  fSt=kCoG;
}//CoG()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliRICHCluster::CorrSin() 
{
// Correction of cluster x position due to sinoid, see HMPID TDR  page 30
// Arguments: none
//   Returns: none
  AliRICHDigit dig(Ch(),100,1,fX,fY);                                               //tmp digit to get it center
  Float_t x=fX-dig.LorsX();  
  fX+=3.31267e-2*TMath::Sin(2*TMath::Pi()/0.8*x)-2.66575e-3*TMath::Sin(4*TMath::Pi()/0.8*x)+2.80553e-3*TMath::Sin(6*TMath::Pi()/0.8*x)+0.0070;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliRICHCluster::FitFunc(Int_t &iNpars, Double_t *, Double_t &chi2, Double_t *par, Int_t )
{
// Cluster fit function 
// par[0]=x par[1]=y par[2]=q for the first Mathieson shape
// par[3]=x par[4]=y par[5]=q for the second Mathieson shape and so on up to iNpars/3 Mathieson shapes
// For each pad of the cluster calculates the difference between actual pad charge and the charge induced to this pad by all Mathieson distributions 
// Then the chi2 is calculated as the sum of this value squared for all pad in the cluster.  
// Arguments: iNpars - number of parameters which is number of local maxima of cluster * 3
//            chi2   - function result to be minimised 
//            par   - parameters array of size iNpars            
//   Returns: none  
  AliRICHCluster *pClu=(AliRICHCluster*)gMinuit->GetObjectFit();
  Int_t iNshape = iNpars/3;
    
  chi2 = 0;
  for(Int_t i=0;i<pClu->Size();i++){                                       //loop on all pads of the cluster
    Double_t dQpadMath = 0;                                                //pad charge collector  
    for(Int_t j=0;j<iNshape;j++){                                          //Mathiesons loop as all of them may contribute to this pad
      dQpadMath+=par[3*j+2]*pClu->Dig(i)->Mathieson(par[3*j],par[3*j+1]);  // par[3*j+2] is charge par[3*j] is x par[3*j+1] is y of current Mathieson
    }
    chi2 +=TMath::Power((pClu->Dig(i)->Q()-dQpadMath),2);                  //
  }                                                                             //loop on all pads of the cluster     
}//FitFunction()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliRICHCluster::Print(Option_t* opt)const
{
//Print current cluster  
  const char *status=0;
  switch(fSt){
    case      kFor: status="formed"     ;break;
    case      kUnf: status="unfolded"   ;break;
    case      kCoG: status="coged"      ;break;
    case      kEmp: status="empty"      ;break;
  }
  Printf("%s cs=%2i, Size=%2i (x=%7.3f cm,y=%7.3f cm,Q=%4i qdc), %s",
         opt,Ch(),Size(),X(),Y(),Q(),status);
  for(Int_t i=0;i<Size();i++) Dig(i)->Print();    
}//Print()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t AliRICHCluster::Solve(TClonesArray *pCluLst,Bool_t isTryUnfold)
{
//This methode is invoked when the cluster is formed to solve it. Solve the cluster means to try to unfold the cluster
//into the local maxima number of clusters. This methode is invoked by AliRICHRconstructor::Dig2Clu() on cluster by cluster basis.  
//At this point, cluster contains a list of digits, cluster charge and size is precalculated in AddDigit(), position is preset to (-1,-1) in ctor,
//status is preset to kFormed in AddDigit(), chamber-sector info is preseted to actual values in AddDigit()
//Method first finds number of local maxima and if it's more then one tries to unfold this cluster into local maxima number of clusters
//Arguments: pCluLst     - cluster list pointer where to add new cluster(s)
//           isTryUnfold - flag to switch on/off unfolding   
//  Returns: number of local maxima of original cluster

//Phase 0. Initialise TMinuit  
  const Int_t kMaxLocMax=6;                                                            //max allowed number of loc max for fitting
  TMinuit *pMinuit = new TMinuit(3*kMaxLocMax);                                        //init MINUIT with this number of parameters (3 params per mathieson)
  pMinuit->SetObjectFit((TObject*)this);  pMinuit->SetFCN(AliRICHCluster::FitFunc);    //set fit function
  Double_t aArg=-1,parStart,parStep,parLow,parHigh;     Int_t iErrFlg;                 //tmp vars for TMinuit
  pMinuit->mnexcm("SET PRI",&aArg,1,iErrFlg);                                          //suspend all printout from TMinuit 
  pMinuit->mnexcm("SET NOW",&aArg,0,iErrFlg);                                          //suspend all warning printout from TMinuit
//Phase 1. Find number of local maxima. Strategy is to check if the current pad has QDC more then all neigbours   
  Int_t iLocMaxCnt=0;
  for(Int_t iDig1=0;iDig1<Size();iDig1++) {                                             //first digits loop
    AliRICHDigit *pDig1 = Dig(iDig1);                                                   //take next digit
    Int_t iHowManyMoreCnt = 0;                                                          //counts how many neighbouring pads has QDC more then current one
    for(Int_t iDig2=0;iDig2<Size();iDig2++) {                                           //loop on all digits again
      if(iDig1==iDig2) continue;                                                        //the same digit, no need to compare 
      AliRICHDigit *pDig2 = Dig(iDig2);                                                 //take second digit to compare with the first one
      Int_t dist = TMath::Sign(Int_t(pDig1->PadX()-pDig2->PadX()),1)+TMath::Sign(Int_t(pDig1->PadY()-pDig2->PadY()),1);//distance between pads
      if(dist==1)                                                                       //means dig2 is a neighbour of dig1
         if(pDig2->Q()>=pDig1->Q()) iHowManyMoreCnt++;                                 //count number of pads with Q more then Q of current pad
    }//second digits loop
    if(iHowManyMoreCnt==0&&iLocMaxCnt<kMaxLocMax){                                     //this pad has Q more then any neighbour so it's local maximum
        pMinuit->mnparm(3*iLocMaxCnt  ,Form("x%i",iLocMaxCnt),parStart=pDig1->LorsX(),parStep=0.01,parLow=0,parHigh=0,iErrFlg);
        pMinuit->mnparm(3*iLocMaxCnt+1,Form("y%i",iLocMaxCnt),parStart=pDig1->LorsY(),parStep=0.01,parLow=0,parHigh=0,iErrFlg);
        pMinuit->mnparm(3*iLocMaxCnt+2,Form("q%i",iLocMaxCnt),parStart=pDig1->Q()    ,parStep=0.01,parLow=0,parHigh=0,iErrFlg);
        iLocMaxCnt++;
    }//if this pad is local maximum
  }//first digits loop
//Phase 2. Fit loc max number of Mathiesons or add this current cluster to the list
  Int_t iCluCnt=pCluLst->GetEntriesFast();                                          //get current number of clusters already stored in the list by previous operations
  if(isTryUnfold==kTRUE && iLocMaxCnt<kMaxLocMax){                                        //resonable number of local maxima to fit and user requested it
    pMinuit->mnexcm("MIGRAD" ,&aArg,0,iErrFlg);                                     //start fitting
    if (!iErrFlg) { // Only if MIGRAD converged normally
      Double_t fitX,fitY,fitQ,d1,d2,d3; TString sName;                                //vars to get results from TMinuit
      for(Int_t i=0;i<iLocMaxCnt;i++){//local maxima loop
	pMinuit->mnpout(3*i   ,sName,  fitX, d1 , d2, d3, iErrFlg);
	pMinuit->mnpout(3*i+1 ,sName,  fitY, d1 , d2, d3, iErrFlg);
	pMinuit->mnpout(3*i+2 ,sName,  fitQ, d1 , d2, d3, iErrFlg);
	if (TMath::Abs(fitQ)>2147483647.0) fitQ = TMath::Sign((Double_t)2147483647,fitQ);//???????????????
	new ((*pCluLst)[iCluCnt++]) AliRICHCluster(Ch(),fitX,fitY,(Int_t)fitQ);	    //add new unfolded clusters
      }//local maxima loop
    }
  }else{//do not unfold since number of loc max is unresonably high or user's baned unfolding 
    CoG();
    new ((*pCluLst)[iCluCnt++]) AliRICHCluster(*this);  //add this raw cluster 
  }
  delete pMinuit;
  return iLocMaxCnt;
}//Solve()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
