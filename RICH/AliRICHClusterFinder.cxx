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


#include "AliRICHClusterFinder.h"
#include "AliRICHMap.h"
#include <TMinuit.h>
#include <TParticle.h>
#include <TVector3.h>
#include <AliLoader.h>
#include <AliStack.h>
#include <AliRun.h>

void RICHMinMathieson(Int_t &npar, Double_t *gin, Double_t &chi2, Double_t *par, Int_t iflag);

ClassImp(AliRICHClusterFinder)
//__________________________________________________________________________________________________
AliRICHClusterFinder::AliRICHClusterFinder(AliRICH *pRICH)   
{//main ctor
  fRICH = pRICH;
  AliDebug(1,"main ctor Start.");
  
  fDigitMap = 0;
  fRawCluster.Reset();
  fResolvedCluster.Reset();
  AliDebug(1,"main ctor Stop.");
}//main ctor
//__________________________________________________________________________________________________
void AliRICHClusterFinder::Exec()
{
//Main method of cluster finder. Loops on  events and chambers, everything else is done in FindClusters()  
  AliDebug(1,"Exec Start.");
    
  R()->GetLoader()                ->LoadDigits();   
//  R()->GetLoader()->GetRunLoader()->LoadHeader(); 
  R()->GetLoader()->GetRunLoader()->LoadKinematics(); //header is already loaded

  for(Int_t iEventN=0;iEventN<gAlice->GetEventsPerRun();iEventN++){//events loop
    AliDebug(1,Form("Processing event %i...",iEventN));
    R()->GetLoader()->GetRunLoader()->GetEvent(iEventN);
    
    R()->GetLoader()->MakeTree("R");  R()->MakeBranch("R");
    R()->ResetDigits();               R()->ResetClusters();
    
    R()->GetLoader()->TreeD()->GetEntry(0);
    for(Int_t iChamber=1;iChamber<=kNchambers;iChamber++){//chambers loop
      FindClusters(iChamber);
    }//chambers loop
    R()->GetLoader()->TreeR()->Fill();  R()->GetLoader()->WriteRecPoints("OVERWRITE");//write out clusters for current event
  }//events loop  
  
  R()->ResetDigits();//reset and unload everything
  R()->ResetClusters();
  R()->GetLoader()                ->UnloadDigits(); 
  R()->GetLoader()                ->UnloadRecPoints();  
//  R()->GetLoader()->GetRunLoader()->UnloadHeader();
  R()->GetLoader()->GetRunLoader()->UnloadKinematics();

  AliDebug(1,"Stop.");      
}//Exec()
//__________________________________________________________________________________________________
void AliRICHClusterFinder::FindClusters(Int_t iChamber)
{
//Loops on digits for a given chamber, forms raw clusters, then tries to resolve them if requested
  Int_t iNdigits=R()->Digits(iChamber)->GetEntriesFast();
  AliDebug(1,Form("Start for chamber %i with %i digits.",iChamber,iNdigits));  
  
  if(iNdigits==0)return;//no digits for a given chamber, nothing to do

  fDigitMap=new AliRICHMap(R()->Digits(iChamber));//create digit map for the given chamber

  for(Int_t iDigN=0;iDigN<iNdigits;iDigN++){//digits loop for a given chamber    
    AliRICHdigit *dig=(AliRICHdigit*)R()->Digits(iChamber)->At(iDigN);
    Int_t i=dig->X();   Int_t j=dig->Y();
    if(fDigitMap->TestHit(i,j)==kUsed) continue;//this digit is already taken, go after next digit
	
    FormRawCluster(i,j);//form raw cluster starting from (i,j) pad 
    AliDebug(1,"After FormRawCluster:");ToAliDebug(1,fRawCluster.Print());  
    FindLocalMaxima();  //find number of local maxima and initial center of gravity
    AliDebug(1,"After FindLocalMaxima:");ToAliDebug(1,fRawCluster.Print());  
    
    if(AliRICHParam::IsResolveClusters()&&fRawCluster.Size()>1&&fRawCluster.Size()<6){
      FitCoG(); //serialization of resolved clusters will happen inside
    }else{//cluster size=1 or resolving is switched off
      WriteRawCluster();//simply output the formed raw cluster without deconvolution
    }
    fRawCluster.Reset(); fResolvedCluster.Reset();
  }//digits loop for a given chamber

  delete fDigitMap;
  
  AliDebug(1,"Stop.");  
}//FindClusters()
//__________________________________________________________________________________________________
void AliRICHClusterFinder::FindClusterContribs(AliRICHcluster *pCluster)
{
//Finds cerenkov-feedback-mip mixture for a given cluster
  AliDebug(1,"Start.");ToAliDebug(1,pCluster->Print());

  R()->GetLoader()->GetRunLoader()->LoadHeader();
  AliStack *pStack = R()->GetLoader()->GetRunLoader()->Stack();
  if(!pStack)
  {AliInfo("No Stack found!!! No contrib to cluster found.");return;}
  
  TObjArray *pDigits = pCluster->Digits();
  if(!pDigits) return; //??????????
  Int_t iNmips=0,iNckovs=0,iNfeeds=0;
  TArrayI contribs(3*pCluster->Size());
  Int_t *pindex = new Int_t[3*pCluster->Size()];
  for(Int_t iDigN=0;iDigN<pCluster->Size();iDigN++) {//loop on digits of a given cluster
    contribs[3*iDigN]  =((AliRICHdigit*)pDigits->At(iDigN))->GetTrack(0);
    if (contribs[3*iDigN] >= 10000000) contribs[3*iDigN] = 0;
    contribs[3*iDigN+1]=((AliRICHdigit*)pDigits->At(iDigN))->GetTrack(1);
    if (contribs[3*iDigN+1] >= 10000000) contribs[3*iDigN+1] = 0;
    contribs[3*iDigN+2]=((AliRICHdigit*)pDigits->At(iDigN))->GetTrack(2);
    if (contribs[3*iDigN+2] >= 10000000) contribs[3*iDigN+2] = 0;
  }//loop on digits of a given cluster
  TMath::Sort(contribs.GetSize(),contribs.GetArray(),pindex);
  for(Int_t iDigN=0;iDigN<3*pCluster->Size()-1;iDigN++) {//loop on digits to sort tids
    AliDebug(1,Form("%4i for digit n. %4i",contribs[pindex[iDigN]],iDigN));
    if(contribs[pindex[iDigN]]!=contribs[pindex[iDigN+1]]) {
      TParticle* particle = pStack->Particle(contribs[pindex[iDigN]]);
      Int_t code   = particle->GetPdgCode();
      Double_t charge = 0;
      if(particle->GetPDG()) charge=particle->GetPDG()->Charge();
      AliDebug(1,Form(" charge of particle %f",charge));

      if(code==50000050) iNckovs++;
      if(code==50000051) iNfeeds++;
      if(charge!=0) iNmips++;
    }
  }//loop on digits to sort Tid
  
  if (contribs[pindex[3*pCluster->Size()-1]]!=kBad) {

     TParticle* particle = pStack->Particle(contribs[pindex[3*pCluster->Size()-1]]);
     Int_t code   = particle->GetPdgCode();
     Double_t charge = 0;
     if(particle->GetPDG()) charge=particle->GetPDG()->Charge();
     AliDebug(1,Form(" charge of particle %f",charge));
     if(code==50000050) iNckovs++;
     if(code==50000051) iNfeeds++;
     if(charge!=0) iNmips++;
  }
    
  pCluster->CFM(iNckovs,iNfeeds,iNmips);
//  
  delete [] pindex; 
  ToAliDebug(1,pCluster->Print());
  AliDebug(1,"Stop.");
}//FindClusterContribs()
//__________________________________________________________________________________________________
void  AliRICHClusterFinder::FormRawCluster(Int_t i, Int_t j)
{
//Builds the raw cluster (before deconvolution). Starts from the first pad (i,j) then calls itself recursevly  for all neighbours.
  AliDebug(1,Form("Start with digit(%i,%i) Q=%f",i,j,((AliRICHdigit*)fDigitMap->GetHit(i,j))->Q()));
  
  fRawCluster.AddDigit((AliRICHdigit*) fDigitMap->GetHit(i,j));//take this pad in cluster
  fDigitMap->FlagHit(i,j);//flag this pad as taken  

  Int_t listX[4], listY[4];    //  Now look recursively for all neighbours
  for (Int_t iNei=0;iNei<R()->P()->PadNeighbours(i,j,listX,listY);iNei++)
    if(fDigitMap->TestHit(listX[iNei],listY[iNei])==kUnused) FormRawCluster(listX[iNei],listY[iNei]);    
}//FormRawCluster()
//__________________________________________________________________________________________________
void AliRICHClusterFinder::FindLocalMaxima()
{
//find number of local maxima in the current raw cluster and then calculates initial center of gravity
  fNlocals=0;
  for(Int_t iDig1=0;iDig1<fRawCluster.Size();iDig1++) {
    Int_t iNotMax = 0;
    AliRICHdigit *pDig1 = (AliRICHdigit *)fRawCluster.Digits()->At(iDig1);
    TVector pad1 = pDig1->Pad();
    Int_t padQ1 = (Int_t)(pDig1->Q()+0.1);
    Int_t padC1 = pDig1->ChFbMi();
    for(Int_t iDig2=0;iDig2<fRawCluster.Size();iDig2++) {
      AliRICHdigit *pDig2 = (AliRICHdigit *)fRawCluster.Digits()->At(iDig2);
      TVector pad2 = pDig2->Pad();
      Int_t padQ2 = (Int_t)(pDig2->Q()+0.1);
      if(iDig1==iDig2) continue;
      Int_t diffx = TMath::Sign(Int_t(pad1[0]-pad2[0]),1);
      Int_t diffy = TMath::Sign(Int_t(pad1[1]-pad2[1]),1);
      if((diffx+diffy)<=1) {
         if(padQ2>=padQ1) iNotMax++;
      }
    }
    if(iNotMax==0) {
      TVector2 x2=AliRICHParam::Pad2Loc(pad1);
      fLocalX[fNlocals]=x2.X();fLocalY[fNlocals]=x2.Y();
      fLocalQ[fNlocals] = (Double_t)padQ1;
      fLocalC[fNlocals] = padC1;
      fNlocals++;
    }
  }
  fRawCluster.CoG(fNlocals); //first initial approximation of the CoG...to start minimization.
}//FindLocalMaxima()
//__________________________________________________________________________________________________
void AliRICHClusterFinder::WriteRawCluster()
{
//Add the current raw cluster to the list of clusters
  AliDebug(1,"Start.");
  
  FindClusterContribs(&fRawCluster);  
  R()->AddCluster(fRawCluster);
  
  ToAliDebug(1,fRawCluster.Print()); AliDebug(1,"Stop."); 
}//WriteRawCluster()
//__________________________________________________________________________________________________
void AliRICHClusterFinder::WriteResolvedCluster()
{
//Add the current resolved cluster to the list of clusters
  AliDebug(1,"Start.");
  
  FindClusterContribs(&fResolvedCluster);  
  R()->AddCluster(fResolvedCluster);
  
  ToAliDebug(1,fResolvedCluster.Print()); AliDebug(1,"Stop.");  
}//WriteResolvedCluster()
//__________________________________________________________________________________________________
void AliRICHClusterFinder::FitCoG()
{
//Fits cluster of size  by the corresponding number of Mathieson shapes.
//This methode is only invoked in case everything is ok to start deconvolution  
  AliDebug(1,"Start with:"); ToAliDebug(1,fRawCluster.Print());
  
  Double_t arglist;
  Int_t ierflag = 0;

  TMinuit *pMinuit = new TMinuit(3*fNlocals-1);
  pMinuit->mninit(5,10,7);
  
  arglist = -1;
  pMinuit->mnexcm("SET PRI",&arglist, 1, ierflag);
  pMinuit->mnexcm("SET NOW",&arglist, 0, ierflag);
  
  TString chname;
  Int_t ierflg;
  
  pMinuit->SetObjectFit((TObject*)this);
  pMinuit->SetFCN(RICHMinMathieson);
  
  Double_t vstart,lower, upper;
  Double_t stepX= 0.001;
  Double_t stepY= 0.001;
  Double_t stepQ= 0.0001;
  
  for(Int_t i=0;i<fNlocals;i++) {
    vstart   = fLocalX[i];
    lower    = vstart - 2*AliRICHParam::PadSizeX();
    upper    = vstart + 2*AliRICHParam::PadSizeX();
    pMinuit->mnparm(3*i  ,Form("xCoG  %i",i),vstart,stepX,lower,upper,ierflag);
    vstart   = fLocalY[i];
    lower    = vstart - 2*AliRICHParam::PadSizeY();
    upper    = vstart + 2*AliRICHParam::PadSizeY();
    pMinuit->mnparm(3*i+1,Form("yCoG  %i",i),vstart,stepY,lower,upper,ierflag);
    if(i==fNlocals-1) break;                    // last parameter is constrained
    vstart = fLocalQ[i]/fRawCluster.Q();
    lower  = 0;
    upper  = 1;
    pMinuit->mnparm(3*i+2,Form("qfrac %i",i),vstart,stepQ,lower,upper,ierflag);
  }
  
  arglist = -1;
  pMinuit->mnexcm("SET NOGR",&arglist, 1, ierflag);
  arglist = 1;
  pMinuit->mnexcm("SET ERR", &arglist, 1,ierflg);
  arglist = -1;
  pMinuit->mnexcm("SIMPLEX",&arglist, 0, ierflag);
  pMinuit->mnexcm("MIGRAD",&arglist, 0, ierflag);
  pMinuit->mnexcm("EXIT" ,&arglist, 0, ierflag);

  Double_t xCoG[50],yCoG[50],qfracCoG[50];
  Double_t eps, b1, b2;

  Double_t qfraclast=0;  
  for(Int_t i=0;i<fNlocals;i++) {
    pMinuit->mnpout(3*i  ,chname,     xCoG[i], eps , b1, b2, ierflg);
    pMinuit->mnpout(3*i+1,chname,     yCoG[i], eps , b1, b2, ierflg);
    if(i==fNlocals-1) break;
    pMinuit->mnpout(3*i+2,chname, qfracCoG[i], eps , b1, b2, ierflg);
    qfraclast+=qfracCoG[i];
   }
  qfracCoG[fNlocals-1] = 1 - qfraclast;
  delete pMinuit;

  for(Int_t i=0;i<fNlocals;i++){//resolved positions loop
    fResolvedCluster.Fill(&fRawCluster,xCoG[i],yCoG[i],qfracCoG[i],fLocalC[i]);
    WriteResolvedCluster();
  }
  AliDebug(1,"Stop.");
}//FitCoG()
//__________________________________________________________________________________________________
void RICHMinMathieson(Int_t &npar, Double_t *, Double_t &chi2, Double_t *par, Int_t )
{
//Mathieson minimization function 
  
  AliRICHcluster *pRawCluster = ((AliRICHClusterFinder*)gMinuit->GetObjectFit())->GetRawCluster();

  TVector2 centroid[50];
  Double_t q[50];
  Int_t nFunctions = (npar+1)/3;
  Double_t qfract = 0;
  for(Int_t i=0;i<nFunctions;i++) {
    centroid[i].Set(par[3*i],par[3*i+1]);
    if(i==nFunctions-1) break;
    q[i]=par[3*i+2];
    qfract+=q[i];
  }
  q[nFunctions-1] = 1 - qfract;
    
  chi2 = 0;
  Int_t qtot = pRawCluster->Q();
  for(Int_t i=0;i<pRawCluster->Size();i++) {
    TVector  pad=((AliRICHdigit *)pRawCluster->Digits()->At(i))->Pad();
    Double_t padQ = ((AliRICHdigit *)pRawCluster->Digits()->At(i))->Q();
    Double_t qfracpar=0;
    for(Int_t j=0;j<nFunctions;j++) {
      qfracpar += q[j]*AliRICHParam::FracQdc(centroid[j],pad);
    }
    chi2 += TMath::Power((qtot*qfracpar-padQ),2)/padQ;
  }     
}//RICHMinMathieson()
//__________________________________________________________________________________________________
