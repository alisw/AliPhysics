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

#include <TVirtualFitter.h>  //Solve()
#include <TMinuit.h>         //Solve()
#include <TClonesArray.h>    //Solve()
#include <TMarker.h>         //Draw()

#include "AliLog.h"          //FitFunc()

#include "AliHMPIDCluster.h"  //class header

Bool_t AliHMPIDCluster::fgDoCorrSin=kTRUE;

ClassImp(AliHMPIDCluster)
    
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDCluster::CoG()
{
// Calculates naive cluster position as a center of gravity of its digits.
// Arguments: none 
//   Returns: none
  Int_t minPadX=999,minPadY=999,maxPadX=-1,maxPadY=-1;      //for box finding  
  if(fDigs==0) return;                                      //no digits in this cluster
  fX=fY=fQRaw=0;                                            //init summable parameters
  Int_t maxQpad=-1,maxQ=-1;                                 //to calculate the pad with the highest charge
  AliHMPIDDigit *pDig=0x0;
  for(Int_t iDig=0;iDig<fDigs->GetEntriesFast();iDig++){    //digits loop
    pDig=(AliHMPIDDigit*)fDigs->At(iDig);                   //get pointer to next digit

    if(pDig->PadPcX() > maxPadX) maxPadX = pDig->PadPcX();  // find the minimum box that contain the cluster  MaxX                            
    if(pDig->PadPcY() > maxPadY) maxPadY = pDig->PadPcY();  //                                                MaxY
    if(pDig->PadPcX() < minPadX) minPadX = pDig->PadPcX();  //                                                MinX   
    if(pDig->PadPcY() < minPadY) minPadY = pDig->PadPcY();  //                                                MinY   
    
    Float_t q=pDig->Q();                                    //get QDC 
    fX += pDig->LorsX()*q;fY +=pDig->LorsY()*q;             //add digit center weighted by QDC
    fQRaw+=q;                                               //increment total charge 
    if(q>maxQ) {maxQpad = pDig->Pad();maxQ=(Int_t)q;}       // to find pad with highest charge
  }//digits loop
  
  fBox=(maxPadX-minPadX+1)*100+maxPadY-minPadY+1;           // dimension of the box: format Xdim*100+Ydim
  
  if ( fQRaw != 0 )   fX/=fQRaw;fY/=fQRaw;                  //final center of gravity
   
  if(fDigs->GetEntriesFast()>1&&fgDoCorrSin)CorrSin();       //correct it by sinoid   
  
  fQ  = fQRaw;                                              // Before starting fit procedure, Q and QRaw must be equal
  fCh=pDig->Ch();                                           //initialize chamber number
  fMaxQpad = maxQpad; fMaxQ=maxQ;                           //store max charge pad to the field
  fChi2=0;                                                  // no Chi2 to find
  fNlocMax=0;                                               // proper status from this method
  fSt=kCoG;
}//CoG()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDCluster::CorrSin() 
{
// Correction of cluster x position due to sinoid, see HMPID TDR  page 30
// Arguments: none
//   Returns: none
  Int_t pc,px,py;
  AliHMPIDParam::Lors2Pad(fX,fY,pc,px,py);             //tmp digit to get it center
  Float_t x=fX-AliHMPIDParam::LorsX(pc,px);                    //diff between cluster x and center of the pad contaning this cluster   
  fX+=3.31267e-2*TMath::Sin(2*TMath::Pi()/0.8*x)-2.66575e-3*TMath::Sin(4*TMath::Pi()/0.8*x)+2.80553e-3*TMath::Sin(6*TMath::Pi()/0.8*x)+0.0070;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDCluster::Draw(Option_t*)
{
  TMarker *pMark=new TMarker(X(),Y(),5); pMark->SetUniqueID(fSt);pMark->SetMarkerColor(kBlue); pMark->Draw();
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDCluster::FitFunc(Int_t &iNpars, Double_t* deriv, Double_t &chi2, Double_t *par, Int_t iflag)
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
  
  AliHMPIDCluster *pClu=(AliHMPIDCluster*)TVirtualFitter::GetFitter()->GetObjectFit();

  Int_t nPads = pClu->Size();  
  
  chi2 = 0;
  
  Int_t iNshape = iNpars/3;
  
  for(Int_t i=0;i<nPads;i++){                                                          //loop on all pads of the cluster
    Double_t dQpadMath = 0;
    for(Int_t j=0;j<iNshape;j++){                                                      //Mathiesons loop as all of them may contribute to this pad
      Double_t fracMathi = pClu->Dig(i)->IntMathieson(par[3*j],par[3*j+1]);
      dQpadMath+=par[3*j+2]*fracMathi;                                                 // par[3*j+2] is charge par[3*j] is x par[3*j+1] is y of current Mathieson
    }
    if(dQpadMath>0 && pClu->Dig(i)->Q()>0) {
      chi2 +=TMath::Power((pClu->Dig(i)->Q()-dQpadMath),2)/pClu->Dig(i)->Q();          //chi2 function to be minimized
    }
  }
//---calculate gradients...  
  if(iflag==2) {
    Double_t **derivPart;

    derivPart = new Double_t*[iNpars];

    for(Int_t j=0;j<iNpars;j++){                                                      
      deriv[j] = 0;
      derivPart[j] = new Double_t[nPads];
      for(Int_t i=0;i<nPads;i++){                                                          
        derivPart[j][i] = 0;
      }
    }

    for(Int_t i=0;i<nPads;i++){                                                          //loop on all pads of the cluster
      for(Int_t j=0;j<iNshape;j++){                                                      //Mathiesons loop as all of them may contribute to this pad
        Double_t fracMathi = pClu->Dig(i)->IntMathieson(par[3*j],par[3*j+1]);
        derivPart[3*j  ][i] += par[3*j+2]*(pClu->Dig(i)->Mathieson(par[3*j]-pClu->Dig(i)->LorsX()-0.5*AliHMPIDParam::SizePadX())-
                                           pClu->Dig(i)->Mathieson(par[3*j]-pClu->Dig(i)->LorsX()+0.5*AliHMPIDParam::SizePadX()))*
                                           pClu->Dig(i)->IntPartMathi(par[3*j+1],2);
        derivPart[3*j+1][i] += par[3*j+2]*(pClu->Dig(i)->Mathieson(par[3*j+1]-pClu->Dig(i)->LorsY()-0.5*AliHMPIDParam::SizePadY())-
                                           pClu->Dig(i)->Mathieson(par[3*j+1]-pClu->Dig(i)->LorsY()+0.5*AliHMPIDParam::SizePadY()))*
                                           pClu->Dig(i)->IntPartMathi(par[3*j],1);
        derivPart[3*j+2][i] += fracMathi;
      }
    }
                                                                                         //loop on all pads of the cluster     
    for(Int_t i=0;i<nPads;i++){                                                          //loop on all pads of the cluster
      Double_t dQpadMath = 0;                                                            //pad charge collector  
      for(Int_t j=0;j<iNshape;j++){                                                      //Mathiesons loop as all of them may contribute to this pad
        Double_t fracMathi = pClu->Dig(i)->IntMathieson(par[3*j],par[3*j+1]);
        dQpadMath+=par[3*j+2]*fracMathi;                                                 
        if(dQpadMath>0 && pClu->Dig(i)->Q()>0) {
          deriv[3*j]   += 2/pClu->Dig(i)->Q()*(pClu->Dig(i)->Q()-dQpadMath)*derivPart[3*j  ][i];
          deriv[3*j+1] += 2/pClu->Dig(i)->Q()*(pClu->Dig(i)->Q()-dQpadMath)*derivPart[3*j+1][i];
          deriv[3*j+2] += 2/pClu->Dig(i)->Q()*(pClu->Dig(i)->Q()-dQpadMath)*derivPart[3*j+2][i];
        }
      }
    }
    //delete array...
    for(Int_t i=0;i<iNpars;i++) delete [] derivPart[i]; delete [] derivPart;
  }
//---gradient calculations ended
  
}//FitFunction()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDCluster::Print(Option_t* opt)const
{
//Print current cluster  
  const char *status=0;
  switch(fSt){
    case        kFrm  : status="formed        "   ;break;
    case        kUnf  : status="unfolded (fit)"   ;break;
    case        kCoG  : status="coged         "   ;break;
    case        kLo1  : status="locmax 1 (fit)"   ;break;
    case        kMax  : status="exceeded (cog)"   ;break;
    case        kNot  : status="not done (cog)"   ;break;
    case        kEmp  : status="empty         "   ;break;
    case        kEdg  : status="edge     (fit)"   ;break;
    case 	kSi1  : status="size 1   (cog)"   ;break;
    case 	kNoLoc: status="no LocMax(fit)"   ;break;
    case 	kAbn  : status="Abnormal fit  "   ;break;
    
    default:            status="??????"          ;break;   
  }
  Double_t ratio=0;
  if(Q()>0&&QRaw()>0) ratio = Q()/QRaw()*100;
  Printf("%sCLU: ch=%i                 (%7.3f,%7.3f) Q=%8.3f Qraw=%8.3f(%3.0f%%) Size=%2i DimBox=%i LocMax=%i Chi2=%7.3f   %s",
         opt,Ch(),X(),Y(),Q(),QRaw(),ratio,Size(),fBox,fNlocMax,fChi2,status);
  if(fDigs) fDigs->Print();    
}//Print()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t AliHMPIDCluster::Solve(TClonesArray *pCluLst,Bool_t isTryUnfold)
{
//This methode is invoked when the cluster is formed to solve it. Solve the cluster means to try to unfold the cluster
//into the local maxima number of clusters. This methode is invoked by AliHMPIDRconstructor::Dig2Clu() on cluster by cluster basis.  
//At this point, cluster contains a list of digits, cluster charge and size is precalculated in AddDigit(), position is preset to (-1,-1) in ctor,
//status is preset to kFormed in AddDigit(), chamber-sector info is preseted to actual values in AddDigit()
//Method first finds number of local maxima and if it's more then one tries to unfold this cluster into local maxima number of clusters
//Arguments: pCluLst     - cluster list pointer where to add new cluster(s)
//           isTryUnfold - flag to switch on/off unfolding   
//  Returns: number of local maxima of original cluster
  const Int_t kMaxLocMax=6;                                                              //max allowed number of loc max for fitting
//  
  CoG();                                                                                 //First calculate CoG for the given cluster
  Int_t iCluCnt=pCluLst->GetEntriesFast();                                               //get current number of clusters already stored in the list by previous operations
  if(isTryUnfold==kFALSE || Size()==1) {                                                 //if cluster contains single pad there is no way to improve the knowledge 
    (isTryUnfold)?fSt=kSi1:fSt=kNot;
    new ((*pCluLst)[iCluCnt++]) AliHMPIDCluster(*this);  //add this raw cluster 
    return 1;
  } 
  
//Phase 0. Initialise Fitter  
  Double_t arglist[10];
  Int_t ierflg = 0;
  TVirtualFitter *fitter = TVirtualFitter::Fitter(this,3*6);                            //initialize Fitter

  arglist[0] = -1;
  ierflg = fitter->ExecuteCommand("SET PRI", arglist, 1);                               // no printout
  ierflg = fitter->ExecuteCommand("SET NOW", arglist, 0);                               //no warning messages
  arglist[0] =  1;
  ierflg = fitter->ExecuteCommand("SET GRA", arglist, 1);                               //force Fitter to use my gradient

  fitter->SetFCN(AliHMPIDCluster::FitFunc);

//  arglist[0] = 1;
//  ierflg = fitter->ExecuteCommand("SET ERR", arglist ,1);
  
// Set starting values and step sizes for parameters
    
//Phase 1. Find number of local maxima. Strategy is to check if the current pad has QDC more then all neigbours. Also find the box contaning the cluster   
  fNlocMax=0;

  for(Int_t iDig1=0;iDig1<Size();iDig1++) {                                               //first digits loop
    
    AliHMPIDDigit *pDig1 = Dig(iDig1);                                                   //take next digit    
    Int_t iCnt = 0;                                                                      //counts how many neighbouring pads has QDC more then current one
    
    for(Int_t iDig2=0;iDig2<Size();iDig2++) {                                            //loop on all digits again
      
      if(iDig1==iDig2) continue;                                                         //the same digit, no need to compare 
      AliHMPIDDigit *pDig2 = Dig(iDig2);                                                 //take second digit to compare with the first one
      Int_t dist = TMath::Sign(Int_t(pDig1->PadChX()-pDig2->PadChX()),1)+TMath::Sign(Int_t(pDig1->PadChY()-pDig2->PadChY()),1);//distance between pads
      if(dist==1)                                                                        //means dig2 is a neighbour of dig1
         if(pDig2->Q()>=pDig1->Q()) iCnt++;                                              //count number of pads with Q more then Q of current pad
      
    }//second digits loop
    
    if(iCnt==0&&fNlocMax<kMaxLocMax){                                                    //this pad has Q more then any neighbour so it's local maximum
      
      Double_t xStart=pDig1->LorsX();Double_t yStart=pDig1->LorsY();
      Double_t xMin=xStart-AliHMPIDParam::SizePadX();
      Double_t xMax=xStart+AliHMPIDParam::SizePadX();
      Double_t yMin=yStart-AliHMPIDParam::SizePadY();
      Double_t yMax=yStart+AliHMPIDParam::SizePadY();
      
      ierflg = fitter->SetParameter(3*fNlocMax  ,Form("x%i",fNlocMax),xStart,0.1,xMin,xMax);    // X,Y,Q initial values of the loc max pad
      ierflg = fitter->SetParameter(3*fNlocMax+1,Form("y%i",fNlocMax),yStart,0.1,yMin,yMax);    // X, Y constrained to be near the loc max
      ierflg = fitter->SetParameter(3*fNlocMax+2,Form("q%i",fNlocMax),pDig1->Q(),0.1,0,100000); // Q constrained to be positive
      
      fNlocMax++;
      
    }//if this pad is local maximum
  }//first digits loop
  
//Phase 2. Fit loc max number of Mathiesons or add this current cluster to the list
// case 1 -> no loc max found
 if ( fNlocMax == 0) {                                                                       // case of no local maxima found: pads with same charge...
   
   ierflg = fitter->SetParameter(3*fNlocMax  ,Form("x%i",fNlocMax),fX,0.1,0,0);              // Init values taken from CoG() -> fX,fY,fQRaw
   ierflg = fitter->SetParameter(3*fNlocMax+1,Form("y%i",fNlocMax),fY,0.1,0,0);              //
   ierflg = fitter->SetParameter(3*fNlocMax+2,Form("q%i",fNlocMax),fQRaw,0.1,0,100000);      //
   
   fNlocMax = 1;
   fSt=kNoLoc;
 }

// case 2 -> loc max found. Check # of loc maxima 
 if ( fNlocMax >= kMaxLocMax)  {                                                          // if # of local maxima exceeds kMaxLocMax...
   fSt = kMax;   new ((*pCluLst)[iCluCnt++]) AliHMPIDCluster(*this);                      //...add this raw cluster  
   } else {                                                                               //or resonable number of local maxima to fit and user requested it
  // Now ready for minimization step
   arglist[0] = 500;                                                                      //number of steps and sigma on pads charges
   arglist[1] = 1.;                                                                       //

   ierflg = fitter->ExecuteCommand("SIMPLEX",arglist,2);                                  //start fitting with Simplex
   if (!ierflg)
     fitter->ExecuteCommand("MIGRAD" ,arglist,2);                                         //fitting improved by Migrad
   if(ierflg) {
     Double_t strategy=2;
     ierflg = fitter->ExecuteCommand("SET STR",&strategy,1);                              //change level of strategy 
     if(!ierflg) {
       ierflg = fitter->ExecuteCommand("SIMPLEX",arglist,2);                              //start fitting with Simplex
       if (!ierflg)
         fitter->ExecuteCommand("MIGRAD" ,arglist,2);                                     //fitting improved by Migrad
     }
   }        
   if(ierflg) fSt=kAbn;                                                                   //no convergence of the fit...
   Double_t dummy; char sName[80];                                                        //vars to get results from Minuit
   Double_t edm, errdef;
   Int_t nvpar, nparx;
  
   for(Int_t i=0;i<fNlocMax;i++){                                                        //store the local maxima parameters
     fitter->GetParameter(3*i   ,sName,  fX, fErrX , dummy, dummy);                      // X
     fitter->GetParameter(3*i+1 ,sName,  fY, fErrY , dummy, dummy);                      // Y
     fitter->GetParameter(3*i+2 ,sName,  fQ, fErrQ , dummy, dummy);                      // Q
     fitter->GetStats(fChi2, edm, errdef, nvpar, nparx);                                 //get fit infos
      if(fSt!=kAbn) {         
        if(fNlocMax!=1)fSt=kUnf;                                                         // if unfolded
        if(fNlocMax==1&&fSt!=kNoLoc) fSt=kLo1;                                           // if only 1 loc max
        if ( !IsInPc()) fSt = kEdg;                                                      // if Out of Pc
        if(fSt==kNoLoc) fNlocMax=0;                                                      // if with no loc max (pads with same charge..)
      }
      new ((*pCluLst)[iCluCnt++]) AliHMPIDCluster(*this);	                         //add new unfolded cluster
   }
 }

 return fNlocMax;
 
}//Solve()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
