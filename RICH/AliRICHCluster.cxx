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

#include "AliRICHCluster.h"
#include <TMinuit.h>  //Solve()
 
ClassImp(AliRICHCluster)
//__________________________________________________________________________________________________
void AliRICHCluster::Print(Option_t*)const
{
//Print current cluster  
  const char *status=0;
  switch(fStatus){
    case      kFormed: status="formed"     ;break;
    case    kUnfolded: status="unfolded"   ;break;
    case         kCoG: status="CoGed"      ;break;
    case       kEmpty: status="empty"      ;break;
  }
  Int_t iNdigs=0;  if(fDigits) iNdigs=fDigits->GetEntriesFast();
    
  Printf("cfm=%10i, cs=%2i, Size=%2i Maxima=%2i, Shape=%5i, pos=(%7.3f,%7.3f) Q=%6i, %s",
            fCFM,fChamber,Size(),Nlocmax(),fShape,fX,fY,fQdc,status,iNdigs);
  for(Int_t i=0;i<iNdigs;i++) Digit(i)->Print();    
}//Print()


//__________________________________________________________________________________________________
TMinuit* AliRICHCluster::Solve()
{
//At this point, cluster contains a list of digits, cluster charge is precalculated as a sum of digits charges (in AddDigit()),
//position is preset to (-1,-1) (in ctor), status is preset to kFormed in (AddDigit()), chamber-sector info is preseted to actual value (in AddDigit())
//Here we decide what to do with this cluster: unfold or just calculate center of gravity
//Arguments: none
//  Returns: pointer to fitter or 0 if no unfolding decided 
  TMinuit *pMinuit=0;
  if(Size()>=2 && AliRICHParam::IsResolveClusters())
    pMinuit=Unfold();
  else
    CoG(0);
  return pMinuit;
}//Solve()
//__________________________________________________________________________________________________
void AliRICHCluster::FitFunc(Int_t &iNpars, Double_t *, Double_t &chi2, Double_t *aPar, Int_t )
{
//Cluster fit function 
//par[0]=x par[1]=y par[2]=q for the first Mathieson shape
//par[3]=x par[4]=y par[5]=q for the second Mathieson shape and so on up to iNpars/3 Mathieson shapes
//We need to calculate Qpad - Qpadmath summup over all pads of the cluster
//Here Qpad is a actual charge of the pad, Qpadmath is calculated charge of the pad induced by all Mathiesons
//Arguments: iNpars - number of parameters which is number of local maxima of cluster * 3
//           chi2   - function result to be minimised 
//           aPar   - parametrs array of size iNpars            
//  Returns: none  
  AliRICHCluster *pClu=(AliRICHCluster*)gMinuit->GetObjectFit();
  Int_t iNmathiesons = iNpars/3;
    
  TVector2 curMathiesonPos;
  chi2 = 0;
  for(Int_t i=0;i<pClu->Size();i++){//digits loop
    TVector    pad     = pClu->Digit(i)->Pad();
    Double_t dQpad     = pClu->Digit(i)->Qdc();
    Double_t dQpadmath = 0;
    for(Int_t j=0;j<iNmathiesons;j++){//all Mathiesons may contribute to this pad
      curMathiesonPos.Set(aPar[3*j],aPar[3*j+1]);//get current position for current Mathieson
      dQpadmath += aPar[3*j+2]*AliRICHParam::FracQdc(curMathiesonPos,pad);//sums up contributions to the current pad from all Mathiesons
    }
    chi2 += TMath::Power((dQpadmath-dQpad),2);
  }//digits loop     
}//RichClusterFitFunction()
//__________________________________________________________________________________________________
TMinuit* AliRICHCluster::Unfold()
{
//This methode is invoked from Solve() when decided to unfold this cluster
//Method first finds number of local maxima and if it's more then one tries to unfold this cluster into local maxima number of clusters
//Arguments: none
//  Returns: pointer to fitter for retriving parameters
    
  TMinuit *pMinuit = new TMinuit(15); //init MINUIT with max 15 parameters (maxim 5 mathiesons, 3 params per matheson )
  pMinuit->SetObjectFit((TObject*)this);
  pMinuit->SetFCN(AliRICHCluster::FitFunc);//set fit function
  Double_t aArg=-1,parStart,parStep,parLow,parHigh; Int_t iErrFlg; //tmp for MINUIT parameters definitions
  pMinuit->mnexcm("SET PRI" ,&aArg,1,iErrFlg); //suspend all printout from TMinuit 
  
  Int_t iLocMaxCnt=0;
//Strategy is to check if the current pad has QDC more then all neigbours 
  for(Int_t iDig1=0;iDig1<Size();iDig1++) {//first digits loop
    AliRICHDigit *pDig1 = Digit(iDig1);//take the current digit
    Int_t iHowManyMoreCnt = 0;//counts how many neighbouring pads has QDC more then current one
    for(Int_t iDig2=0;iDig2<Size();iDig2++) {//loop on all digits again
      AliRICHDigit *pDig2 = Digit(iDig2);
      if(iDig1==iDig2) continue;             //no need to compare 
      Int_t dist = TMath::Sign(Int_t(pDig1->PadX()-pDig2->PadX()),1)+TMath::Sign(Int_t(pDig1->PadY()-pDig2->PadY()),1);//distance between pads
      if(dist==1)//means pads are neighbours
         if(pDig2->Qdc()>=pDig1->Qdc()) iHowManyMoreCnt++;//count number of pads with Q more then Q of current pad
    }//second digits loop
    if(iHowManyMoreCnt==0&&iLocMaxCnt<6){//this pad has Q more then any neighbour so it's local maximum
        TVector2 x2=AliRICHParam::Pad2Loc(pDig1->Pad());//take pad center position and use it as parameter for current Mathienson shape
        pMinuit->mnparm(3*iLocMaxCnt  ,Form("x%i",iLocMaxCnt),parStart=x2.X()      ,parStep=0.01,parLow=0,parHigh=0,iErrFlg);
        pMinuit->mnparm(3*iLocMaxCnt+1,Form("y%i",iLocMaxCnt),parStart=x2.Y()      ,parStep=0.01,parLow=0,parHigh=0,iErrFlg);
        pMinuit->mnparm(3*iLocMaxCnt+2,Form("q%i",iLocMaxCnt),parStart=pDig1->Qdc(),parStep=0.01,parLow=0,parHigh=0,iErrFlg);//
        iLocMaxCnt++;
    }//if this pad is local maximum
  }//first digits loop
  
  fSize+=iLocMaxCnt;
  if(iLocMaxCnt>0&&iLocMaxCnt<6){ //resonable number of local maxima to fit
    Double_t aArg=0;
    pMinuit->mnexcm("MIGRAD",&aArg,0,iErrFlg);//start fitting
    fStatus=kUnfolded;
  }else{
    delete pMinuit;
    pMinuit=0;
    CoG(0);
  }
  return pMinuit;  
}//Unfold()
//__________________________________________________________________________________________________
void AliRICHCluster::CoG(Int_t nLocals)
{
//Calculates naive cluster position as a center of gravity of its digits.
//Also determines the box fully contaning this cluster
//Arguments:     
  Float_t xmin=999,ymin=999,xmax=0,ymax=0;
  fX=fY=0;
  for(Int_t iDig=0;iDig<Size();iDig++) {
    AliRICHDigit *pDig=Digit(iDig);
    TVector pad=pDig->Pad(); Double_t q=pDig->Qdc();
    TVector2 x2=AliRICHParam::Pad2Loc(pad);
    fX += x2.X()*q;fY +=x2.Y()*q;
    if(pad[0]<xmin)xmin=pad[0];if(pad[0]>xmax)xmax=pad[0];if(pad[1]<ymin)ymin=pad[1];if(pad[1]>ymax)ymax=pad[1];
   }
   fX/=fQdc;fY/=fQdc;//Center of Gravity

   TVector2 center = AliRICHParam::Pad2Loc(AliRICHParam::Loc2Pad(TVector2(fX,fY)));
   fX += AliRICHParam::CogCorr(fX-center.X());//correct cluster position for sinoid

   fShape=Int_t(100*(xmax-xmin+1)+ymax-ymin+1);//find box containing cluster
   fSize+=nLocals;
   fStatus=kCoG;
}//CoG()
//__________________________________________________________________________________________________
void AliRICHCluster::Test(const TVector2 &hitX2,Double_t dEloss)
{
//This is to test all cluster functionality
//Method uses AddDigit() to add a predifined pad structure and then calls Solve   
  Int_t iQtot=AliRICHParam::TotQdc(hitX2,dEloss);
  if(iQtot==0){
    Printf("Provided hit position out of sensitive area");
    return;
  }
  TVector area=AliRICHParam::Loc2Area(hitX2);
  TVector pad(2);
  for(pad[1]=area[1];pad[1]<=area[3];pad[1]++){//affected pads loop first y
    for(pad[0]=area[0];pad[0]<=area[2];pad[0]++){//then x               
      Double_t dQpad=iQtot*AliRICHParam::FracQdc(hitX2,pad);//charge fraction from Mathieson centered at x to pad
      AddDigit(new AliRICHDigit(3,(Int_t)pad[0],(Int_t)pad[1],dQpad));
    }//affected pads loop 
  }
  TMinuit *pMinuit=Solve();
  Print();
  Printf("Initial hit (%.2f,%.2f) Qtot=%i Eloss=%.2f",hitX2.X(),hitX2.Y(),iQtot,dEloss);
  Double_t d1,d2,d3; Int_t iErrFlg;TString sName; //tmp vars for TMinuit
  Double_t x,y,q;
  for(Int_t i=0;i<Nlocmax();i++){//retrive fitting results
    pMinuit->mnpout(3*i   ,sName,  x, d1 , d2, d3, iErrFlg);
    pMinuit->mnpout(3*i+1 ,sName,  y, d1 , d2, d3, iErrFlg);
    pMinuit->mnpout(3*i+2 ,sName,  q, d1 , d2, d3, iErrFlg);
  }
  Printf(" Fitted hit (%.2f,%.2f) Qfit=%.0f",x,y,q);
  delete pMinuit; pMinuit=0;  Reset();
}//Test()
