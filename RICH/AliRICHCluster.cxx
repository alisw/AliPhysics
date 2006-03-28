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
#include <AliStack.h>        //FindCfm(), Solve() 
#include <TParticle.h>       //FindCfm()
#include <TClonesArray.h>    //Solve() Test()

ClassImp(AliRICHCluster)
//__________________________________________________________________________________________________
void AliRICHCluster::CoG()
{
// Calculates naive cluster position as a center of gravity of its digits.
// Arguments: none 
//   Returns: shape of the cluster i.e. the box which fully contains the cluster      
  if(fDigs==0) return;                                  //no digits in this cluster
  fX=fY=0;                                              //set cluster position to (0,0) to start to collect contributions
  for(Int_t iDig=0;iDig<fDigs->GetEntriesFast();iDig++){//digits loop
    AliRICHDigit *pDig=(AliRICHDigit*)fDigs->At(iDig);  //get pointer to next digit
    TVector pad=pDig->Pad(); Double_t q=pDig->Qdc();    //get pad adn QDC of this digit 
    TVector2 x2=AliRICHParam::Pad2Loc(pad);             //calculate center of the pad in LORS
    fX += x2.X()*q;fY +=x2.Y()*q;                       //sum up digit centers weighted by QDC
  }//digits loop
  fX/=fQdc;fY/=fQdc;                                    //final center of gravity

  TVector2 center = AliRICHParam::Pad2Loc(AliRICHParam::Loc2Pad(TVector2(fX,fY)));//center of the pad containing calculated cluster position
  fX += AliRICHParam::CogCorr(fX-center.X());                                     //correct cluster position for sinoid

  fStatus=kCoG;
}//CoG()
//__________________________________________________________________________________________________
void AliRICHCluster::FitFunc(Int_t &iNpars, Double_t *, Double_t &chi2, Double_t *par, Int_t )
{
// Cluster fit function 
// par[0]=x par[1]=y par[2]=q for the first Mathieson shape
// par[3]=x par[4]=y par[5]=q for the second Mathieson shape and so on up to iNpars/3 Mathieson shapes
// We need to calculate QpadExp - QpadMathieson summup over all pads of the cluster
// Here QpadExp is a actual charge of the pad, QpadMathieson is calculated charge of the pad induced by all Mathiesons
// Arguments: iNpars - number of parameters which is number of local maxima of cluster * 3
//            chi2   - function result to be minimised 
//            par   - parameters array of size iNpars            
//   Returns: none  
  AliRICHCluster *pClu=(AliRICHCluster*)gMinuit->GetObjectFit();
  Int_t iNmathiesons = iNpars/3;
    
  TVector2 curMathiesonPos;
  chi2 = 0;
  for(Int_t i=0;i<pClu->Size();i++){                                            //loop on all pads of the cluster
    TVector    pad          = pClu->Dig(i)->Pad();
    Double_t dQpadExp       = pClu->Dig(i)->Qdc();
    Double_t dQpadMathieson = 0;
    for(Int_t j=0;j<iNmathiesons;j++){                                          //Mathiesons loop as all of them may contribute to this pad
      curMathiesonPos.Set(par[3*j],par[3*j+1]);                                 //get position of current Mathieson
      dQpadMathieson += par[3*j+2]*AliRICHParam::FracQdc(curMathiesonPos,pad);  //sums up contributions to the current pad from all Mathiesons
    }
    chi2 += TMath::Power((dQpadMathieson-dQpadExp),2);                          //
  }                                                                             //loop on all pads of the cluster     
}//FitFunction()
//__________________________________________________________________________________________________
void AliRICHCluster::Print(Option_t* opt)const
{
//Print current cluster  
  const char *status=0;
  switch(fStatus){
    case      kFormed: status="formed"     ;break;
    case    kUnfolded: status="unfolded"   ;break;
    case         kCoG: status="coged"      ;break;
    case       kEmpty: status="empty"      ;break;
  }
  Int_t iNdigs=0;  if(fDigs) iNdigs=fDigs->GetEntriesFast();
    
  Printf("%s cs=%2i, Size=%2i (x=%7.3f cm,y=%7.3f cm,Q=%4i qdc), %s",
         opt,fCham,iNdigs,fX,fY,fQdc,status);
  for(Int_t i=0;i<iNdigs;i++) Dig(i)->Print();    
}//Print()
//__________________________________________________________________________________________________
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
      AliRICHDigit *pDig2 = Dig(iDig2);                                                 //take second digit to compare with the first one
      if(iDig1==iDig2) continue;                                                        //the same digit, no need to compare 
      Int_t dist = TMath::Sign(Int_t(pDig1->PadX()-pDig2->PadX()),1)+TMath::Sign(Int_t(pDig1->PadY()-pDig2->PadY()),1);//distance between pads
      if(dist==1)                                                                       //means dig2 is a neighbour of dig1
         if(pDig2->Qdc()>=pDig1->Qdc()) iHowManyMoreCnt++;                              //count number of pads with Q more then Q of current pad
    }//second digits loop
    if(iHowManyMoreCnt==0&&iLocMaxCnt<=kMaxLocMax){                                     //this pad has Q more then any neighbour so it's local maximum
        TVector2 x2=AliRICHParam::Pad2Loc(pDig1->Pad());                                //take pad center position and use it as parameter for current Mathienson shape
        pMinuit->mnparm(3*iLocMaxCnt  ,Form("x%i",iLocMaxCnt),parStart=x2.X()      ,parStep=0.01,parLow=0,parHigh=0,iErrFlg);
        pMinuit->mnparm(3*iLocMaxCnt+1,Form("y%i",iLocMaxCnt),parStart=x2.Y()      ,parStep=0.01,parLow=0,parHigh=0,iErrFlg);
        pMinuit->mnparm(3*iLocMaxCnt+2,Form("q%i",iLocMaxCnt),parStart=pDig1->Qdc(),parStep=0.01,parLow=0,parHigh=0,iErrFlg);
        iLocMaxCnt++;
    }//if this pad is local maximum
  }//first digits loop
//Phase 2. Fit loc max number of Mathiesons or add this current cluster to the list
  Int_t iCluCnt=pCluLst->GetEntriesFast();                                          //get current number of clusters already stored in the list by previous operations
  if(isTryUnfold==kTRUE && iLocMaxCnt<=kMaxLocMax){                                        //resonable number of local maxima to fit and user requested it
    pMinuit->mnexcm("MIGRAD" ,&aArg,0,iErrFlg);                                     //start fitting
    Double_t fitX,fitY,fitQ,d1,d2,d3; TString sName;                                //vars to get results from TMinuit
    for(Int_t i=0;i<iLocMaxCnt;i++){//local maxima loop
      pMinuit->mnpout(3*i   ,sName,  fitX, d1 , d2, d3, iErrFlg);
      pMinuit->mnpout(3*i+1 ,sName,  fitY, d1 , d2, d3, iErrFlg);
      pMinuit->mnpout(3*i+2 ,sName,  fitQ, d1 , d2, d3, iErrFlg);
      new ((*pCluLst)[iCluCnt++]) AliRICHCluster(C(),fitX,fitY,(Int_t)fitQ);        //add new unfolded clusters
    }//local maxima loop
  }else{//do not unfold since number of loc max is unresonably high or user's baned unfolding 
    CoG();
    new ((*pCluLst)[iCluCnt++]) AliRICHCluster(*this);  //add this raw cluster 
  }
  delete pMinuit;
  return iLocMaxCnt;
}//Solve()
//__________________________________________________________________________________________________
void AliRICHCluster::Test(Double_t x,Double_t y,Double_t e,Bool_t isTryUnfold)
{
//This is to test all cluster functionality
//Uses AddDigit() to add a predifined pad structure and then calls Solve   
  TVector2 hitX2(x,y);
  Int_t iQtot=AliRICHParam::TotQdc(hitX2,e);
  if(iQtot==0){
    Printf("Provided hit position out of sensitive area");
    return;
  }
  TVector area=AliRICHParam::Loc2Area(hitX2);
  TVector pad(2);
  AliRICHCluster clu;
  for(pad[1]=area[1];pad[1]<=area[3];pad[1]++){//affected pads loop first y
    for(pad[0]=area[0];pad[0]<=area[2];pad[0]++){//then x               
      Double_t dQpad=iQtot*AliRICHParam::FracQdc(hitX2,pad);//charge fraction from Mathieson centered at x to pad
      clu.DigAdd(new AliRICHDigit(pad,dQpad));
    }//affected pads loop 
  }
                 Printf("Initial hit    :  (%.2f,%.2f) Qtot=%i E=%.2f eV",x,y,iQtot,e*1e9);
  clu.CoG();  clu.Print("Initial cluster:");
  TClonesArray *pCluLst=new TClonesArray("AliRICHCluster",1);
  clu.Solve(pCluLst,isTryUnfold);  
  ((AliRICHCluster *)pCluLst->At(0))->Print("Solved cluster:");
  
  delete pCluLst; clu.Reset();
}//Test()
//__________________________________________________________________________________________________
void AliRICHCluster::Test()
{
//Test cluster builder by a number of predefined digit patterns
//Arguments: none
//  Returns: none
  AliRICHCluster clu; Int_t ch,padx,pady,qdc; TClonesArray *pCluLst=new TClonesArray("AliRICHCluster",10);
  Printf("2 digits vertical cluster");
  clu.DigAdd(new AliRICHDigit(ch=1,padx=3,pady=3,qdc=101));
  clu.DigAdd(new AliRICHDigit(ch=1,padx=3,pady=4,qdc=202)); clu.Print("Formed cluster:");
  clu.Solve(pCluLst,kTRUE);  pCluLst->Print();
  delete pCluLst;
}
