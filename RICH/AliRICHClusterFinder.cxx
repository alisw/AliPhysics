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
//  Info("main ctor","Start.");
  
  fRICH = pRICH;
  
  fHitMap = 0;
  fRawCluster.Reset();
  fResolvedCluster.Reset();
  
}//main ctor
//__________________________________________________________________________________________________
void AliRICHClusterFinder::FindLocalMaxima()
{
//  Find number of local maxima in the cluster
//  Info("FindLocalMaxima","Start.");
  
  fNlocals = 0;
  for(Int_t iDig1=0;iDig1<fRawCluster.Size();iDig1++) {
    Int_t iNotMax = 0;
    AliRICHdigit *pDig1 = (AliRICHdigit *)fRawCluster.Digits()->At(iDig1);
    Int_t padX1 = pDig1->X();
    Int_t padY1 = pDig1->Y();
    Int_t padQ1 = (Int_t)(pDig1->Q()+0.1);
    Int_t padC1 = pDig1->CombiPid();
    for(Int_t iDig2=0;iDig2<fRawCluster.Size();iDig2++) {
      AliRICHdigit *pDig2 = (AliRICHdigit *)fRawCluster.Digits()->At(iDig2);
      Int_t padX2 = pDig2->X();
      Int_t padY2 = pDig2->Y();
      Int_t padQ2 = (Int_t)(pDig2->Q()+0.1);
      if(iDig1==iDig2) continue;
      Int_t diffx = TMath::Sign(padX1-padX2,1);
      Int_t diffy = TMath::Sign(padY1-padY2,1);
      if((diffx+diffy)<=1) {
         if(padQ2>=padQ1) iNotMax++;
      }
    }
    if(iNotMax==0) {
      TVector2 x2=AliRICHParam::Pad2Loc(padX1,padY1);
      fLocalX[fNlocals]=x2.X();fLocalY[fNlocals]=x2.Y();
      fLocalQ[fNlocals] = (Double_t)padQ1;
      fLocalC[fNlocals] = padC1;
      fNlocals++;
    }
  }
}//FindLocalMaxima()
//__________________________________________________________________________________________________
void AliRICHClusterFinder::Exec()
{
  Info("Exec","Start.");
  
  
  Rich()->GetLoader()->LoadDigits(); 
  
  Rich()->GetLoader()->GetRunLoader()->LoadHeader();
  Rich()->GetLoader()->GetRunLoader()->LoadKinematics();

  for(Int_t iEventN=0;iEventN<gAlice->GetEventsPerRun();iEventN++){//events loop
    Info("Exec","Event %i processed.",iEventN+1);
//    gAlice->GetRunLoader()->GetEvent(iEventN);
    Rich()->GetLoader()->GetRunLoader()->GetEvent(iEventN);
    
    Rich()->GetLoader()->MakeTree("R");  Rich()->MakeBranch("R");
    Rich()->ResetDigits();  Rich()->ResetClusters();
    
    Rich()->GetLoader()->TreeD()->GetEntry(0);
    for(Int_t iChamber=1;iChamber<=kNCH;iChamber++){//chambers loop
      FindClusters(iChamber);
    }//chambers loop
    Rich()->GetLoader()->TreeR()->Fill();
    Rich()->GetLoader()->WriteRecPoints("OVERWRITE");
  }//events loop  
  Rich()->GetLoader()->UnloadDigits(); Rich()->GetLoader()->UnloadRecPoints();  
  Rich()->ResetDigits();  Rich()->ResetClusters();
  
  Rich()->GetLoader()->GetRunLoader()->UnloadHeader();
  Rich()->GetLoader()->GetRunLoader()->UnloadKinematics();

  Info("Exec","Stop.");      
}//Exec()
//__________________________________________________________________________________________________
void AliRICHClusterFinder::FindClusters(Int_t iChamber)
{
  //finds neighbours and fill the tree with raw clusters
  Int_t nDigits=Rich()->Digits(iChamber)->GetEntriesFast();
//  Info("FindClusters","Start for Chamber %i with %i digits.",iChamber,nDigits);  
  if(nDigits==0)return;

  fHitMap=new AliRICHMap(Rich()->Digits(iChamber));//create digit map for the given chamber

  for(Int_t iDig=0;iDig<nDigits;iDig++){    
    AliRICHdigit *dig=(AliRICHdigit*)Rich()->Digits(iChamber)->At(iDig);
    Int_t i=dig->X();   Int_t j=dig->Y();
    if(fHitMap->TestHit(i,j)==kUsed) continue;
	
    FormRawCluster(i,j);
    FindLocalMaxima();
//    cout << " fNlocals in FindCluster " << fNlocals << endl;
    fRawCluster.CoG(fNlocals); // first initial approxmation of the CoG...to start minimization.
    
    if(AliRICHParam::IsResolveClusters()) {
      ResolveCluster(); // ResolveCluster serialization will happen inside
    } else {
      WriteRawCluster(); // simply output of the RawCluster found without deconvolution
    }
    fRawCluster.Reset();
    fResolvedCluster.Reset();
  }//digits loop

  delete fHitMap;
//  Info("FindClusters","Stop.");
  
}//FindClusters()
//__________________________________________________________________________________________________
void AliRICHClusterFinder::FindClusterContribs(AliRICHcluster *pCluster)
{
// finds CombiPid for a given cluster
// Info("FindClusterContribs","Start");
  
  TObjArray *pDigits = pCluster->Digits();
  Int_t iNmips=0,iNckovs=0,iNfeeds=0;
  TArrayI contribs(3*pCluster->Size());
  Int_t *pindex = new Int_t[3*pCluster->Size()];
  for(Int_t iDigN=0;iDigN<pCluster->Size();iDigN++) {//loop on digits of a given cluster
    contribs[3*iDigN]  =((AliRICHdigit*)pDigits->At(iDigN))->Tid(0);
    contribs[3*iDigN+1]=((AliRICHdigit*)pDigits->At(iDigN))->Tid(1);
    contribs[3*iDigN+2]=((AliRICHdigit*)pDigits->At(iDigN))->Tid(2);
  }//loop on digits of a given cluster
  TMath::Sort(contribs.GetSize(),contribs.GetArray(),pindex);
  for(Int_t iDigN=0;iDigN<3*pCluster->Size()-1;iDigN++) {//loop on digits to sort Tid
    if(contribs[pindex[iDigN]]!=contribs[pindex[iDigN+1]]) {
      Int_t code   = Rich()->GetLoader()->GetRunLoader()->Stack()->Particle(contribs[pindex[iDigN]])->GetPdgCode();
      Double_t charge = Rich()->GetLoader()->GetRunLoader()->Stack()->Particle(contribs[pindex[iDigN]])->GetPDG()->Charge();

      if(code==50000050) iNckovs++;
      else if(code==50000051) iNfeeds++;
      else if(charge!=0) iNmips++;
      if (contribs[pindex[iDigN+1]]==kBad) break;
    }
  }//loop on digits to sort Tid
  pCluster->SetCombiPid(iNckovs,iNfeeds,iNmips);
//  pCluster->Print();
  delete [] pindex; 
}//FindClusterContribs()
//__________________________________________________________________________________________________
void  AliRICHClusterFinder::FormRawCluster(Int_t i, Int_t j)
{
// Builder of the final Raw Cluster (before deconvolution)  
  if(GetDebug()) Info("FormRawCluster","Start with digit(%i,%i)",i,j);
  
  fRawCluster.AddDigit((AliRICHdigit*) fHitMap->GetHit(i,j));
  fHitMap->FlagHit(i,j);// Flag hit as taken  

  Int_t listX[4], listY[4];    //  Now look recursively for all neighbours
  for (Int_t iNeighbour=0;iNeighbour<Rich()->Param()->PadNeighbours(i,j,listX,listY);iNeighbour++)
    if(fHitMap->TestHit(listX[iNeighbour],listY[iNeighbour])==kUnused) 
                      FormRawCluster(listX[iNeighbour],listY[iNeighbour]);    
}//FormRawCluster()
//__________________________________________________________________________________________________
void AliRICHClusterFinder::ResolveCluster()
{// Decluster algorithm
  if(GetDebug()) {Info("ResolveCluster","Start."); fRawCluster.Print();}
  
  switch (fRawCluster.Size()) {
  
  case 1:                     // nothing to decluster: cluster size = 1
    WriteRawCluster();break; 
  default:                     // cluster size > 1: if=2 FitCoG; if>2 Resolve and FitCoG
    FitCoG();break;
  }     
}//ResolveCluster()
//__________________________________________________________________________________________________
void AliRICHClusterFinder::WriteRawCluster()
{
// out the current raw cluster
  Info("WriteRawCluster","Start.");
  
  FindClusterContribs(&fRawCluster);
  fRawCluster.Dump();
  Rich()->AddCluster(fRawCluster);
//  fRawCluster.Print();
}//WriteRawCluster()
//__________________________________________________________________________________________________
void AliRICHClusterFinder::WriteResolvedCluster()
{
// out the current resolved cluster
  Info("WriteResolvedCluster","Start.");
  
//  FindClusterContribs(&fResolvedCluster);
  Rich()->AddCluster(fResolvedCluster);
  
}//WriteResolvedCluster()
//__________________________________________________________________________________________________
void AliRICHClusterFinder::FitCoG()
{// Fit cluster size 2 by single Mathieson
  Info("FitCoG","Start with size %3i and local maxima %3i",fRawCluster.Size(),fNlocals);
  
  Double_t arglist;
  Int_t ierflag = 0;

//  fRawCluster.Print();
  if(fNlocals==0) fRawCluster.Print();
  if(fNlocals==0||fNlocals>3) {WriteRawCluster();return;}
  
  TMinuit *pMinuit = new TMinuit(3*fNlocals-1);
  pMinuit->mninit(5,10,7);
  
  arglist = -1;
  pMinuit->mnexcm("SET PRI",&arglist, 1, ierflag);
  
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
//    cout << Form("xCoG  %i",i) << "start " << vstart << "lower " << lower << "upper " << upper << endl;
    vstart   = fLocalY[i];
    lower    = vstart - 2*AliRICHParam::PadSizeY();
    upper    = vstart + 2*AliRICHParam::PadSizeY();
    pMinuit->mnparm(3*i+1,Form("yCoG  %i",i),vstart,stepY,lower,upper,ierflag);
//    cout << Form("yCoG  %i",i) << "start " << vstart << "lower " << lower << "upper " << upper << endl;
    if(i==fNlocals-1) break;                    // last parameter is constrained
    vstart = fLocalQ[i]/fRawCluster.Q();
    lower  = 0;
    upper  = 1;
    pMinuit->mnparm(3*i+2,Form("qfrac %i",i),vstart,stepQ,lower,upper,ierflag);
//    cout << Form("qCoG  %i",i) << "start " << vstart << "lower " << lower << "upper " << upper << endl;

  }
  
  arglist = -1;
  
  pMinuit->mnexcm("SET NOGR",&arglist, 1, ierflag);
  pMinuit->mnexcm("SET NOW",&arglist, 1, ierflag);
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

  for(Int_t i=0;i<fNlocals;i++) {

    if(fNlocals==5) {
    cout << " Start writing  deconvolved cluster n." << i << endl;
    cout << " fRawCluster " << &fRawCluster << " xCoG " << xCoG[i] << " yCoG " << yCoG[i] << " qfrac " << qfracCoG[i] << endl;
    cout << " Combipid " << fLocalC[i] << " for maximum n. " << i << endl; 
  }
    fResolvedCluster.Fill(&fRawCluster,xCoG[i],yCoG[i],qfracCoG[i],fLocalC[i]);
//    fResolvedCluster.Print();
    WriteResolvedCluster();
  }
if(fNlocals==5)  Info("CoG","Stop.");
}//FitCoG()
//__________________________________________________________________________________________________
void RICHMinMathieson(Int_t &npar, Double_t *, Double_t &chi2, Double_t *par, Int_t )
{
// Mathieson minimization function 
  
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
    Int_t    padX = ((AliRICHdigit *)pRawCluster->Digits()->At(i))->X();
    Int_t    padY = ((AliRICHdigit *)pRawCluster->Digits()->At(i))->Y();
    Double_t padQ = ((AliRICHdigit *)pRawCluster->Digits()->At(i))->Q();
    Double_t qfracpar=0;
    for(Int_t j=0;j<nFunctions;j++) {
      qfracpar += q[j]*AliRICHParam::FracQdc(centroid[j],padX,padY);
//      cout << " function n. " << j+1 << " q " << q[j] << " fracMat " << AliRICHParam::FracQdc(centroid[j],padX,padY) 
//           << " xpar " << centroid[j].X() << " ypar " << centroid[j].Y() << endl;
    }
    chi2 += TMath::Power((qtot*qfracpar-padQ),2)/padQ;
//    cout << " chi2 " << chi2 << " pad n. " << i << " qfracpar " << qfracpar << endl;
  }     
//  if(iflag==3) {
//    cout << " chi2 final " << chi2 << endl;
//  }
}//RICHMinMathieson()
