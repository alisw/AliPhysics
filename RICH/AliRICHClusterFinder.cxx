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
#include <TVector3.h>
#include <AliLoader.h>
#include <AliRun.h>

void RICHMinMathieson(Int_t &npar, Double_t *gin, Double_t &chi2, Double_t *par, Int_t iflag);

ClassImp(AliRICHClusterFinder)
//__________________________________________________________________________________________________
AliRICHClusterFinder::AliRICHClusterFinder(AliRICH *pRICH)   
{//main ctor
  Info("main ctor","Start.");
  
  fRICH = pRICH;
  
  fHitMap = 0;  
  
}//main ctor
//__________________________________________________________________________________________________
void AliRICHClusterFinder::FindLocalMaxima()
{// Split the cluster according to the number of maxima inside
  Info("FindLocalMaxima","Start.");
  Int_t Nlocal = 0;
  Int_t localX[100],localY[100];
  for(Int_t iDig1=0;iDig1<fCurrentCluster.Size();iDig1++) {
    Int_t iNotMax = 0;
    AliRICHdigit *pDig1 = (AliRICHdigit *)fCurrentCluster.Digits()->At(iDig1);
    Int_t padX1 = pDig1->X();
    Int_t padY1 = pDig1->Y();
    Double_t padQ1 = pDig1->Q();
    for(Int_t iDig2=0;iDig2<fCurrentCluster.Size();iDig2++) {
      AliRICHdigit *pDig2 = (AliRICHdigit *)fCurrentCluster.Digits()->At(iDig2);
      Int_t padX2 = pDig2->X();
      Int_t padY2 = pDig2->Y();
      Double_t padQ2 = pDig2->Q();
      if(iDig1==iDig2) continue;
      Int_t diffx = TMath::Sign(padX1-padX2,1);
      Int_t diffy = TMath::Sign(padY1-padY2,1);
      if((diffx+diffy)<=1) {
         if(padQ2>padQ1) iNotMax++;
      }
    }
    if(iNotMax==0) {
    localX[Nlocal] = padX1;
    localY[Nlocal] = padY1;
    Nlocal++;
    }
  }   
}//FindLocalMaxima()
//__________________________________________________________________________________________________
void AliRICHClusterFinder::Exec()
{
  Info("Exec","Start.");
  
  
  Rich()->GetLoader()->LoadDigits(); 
  
  for(Int_t iEventN=0;iEventN<gAlice->GetEventsPerRun();iEventN++){//events loop
    gAlice->GetRunLoader()->GetEvent(iEventN);
    
    Rich()->GetLoader()->MakeTree("R");  Rich()->MakeBranch("R");
    Rich()->ResetDigits();  Rich()->ResetClusters();
    
    Rich()->GetLoader()->TreeD()->GetEntry(0);
    for(Int_t iChamber=1;iChamber<=kNCH;iChamber++){//chambers loop
     FindRawClusters(iChamber);
        
    }//chambers loop
    
    Rich()->GetLoader()->TreeR()->Fill();
    Rich()->GetLoader()->WriteRecPoints("OVERWRITE");
  }//events loop  
  Rich()->GetLoader()->UnloadDigits(); Rich()->GetLoader()->UnloadRecPoints();  
  Rich()->ResetDigits();  Rich()->ResetClusters();
  Info("Exec","Stop.");      
}//Exec()
//__________________________________________________________________________________________________
void AliRICHClusterFinder::FindRawClusters(Int_t iChamber)
{//finds neighbours and fill the tree with raw clusters
  Int_t nDigits=Rich()->Digits(iChamber)->GetEntriesFast();
  Info("FindRawClusters","Start for Chamber %i with %i digits.",iChamber,nDigits);  
  if(nDigits==0)return;

  fHitMap=new AliRICHMap(Rich()->Digits(iChamber));//create digit map for the given chamber

  for(Int_t iDig=0;iDig<nDigits;iDig++){    
    AliRICHdigit *dig=(AliRICHdigit*)Rich()->Digits(iChamber)->At(iDig);
    Int_t i=dig->X();   Int_t j=dig->Y();
    if(fHitMap->TestHit(i,j)==kUsed) continue;
	
    FormRawCluster(i,j);
	
    if(AliRICHParam::IsResolveClusters()) {
      ResolveCluster(); // ResolveCluster serialization will happen inside
    } else {
      WriteRawCluster(); // simply output of the RawCluster found without deconvolution
    }    
    fCurrentCluster.Reset();
  }//digits loop

  delete fHitMap;
  Info("FindRawClusters","Stop.");
  
}//FindRawClusters()
//__________________________________________________________________________________________________
void  AliRICHClusterFinder::FormRawCluster(Int_t i, Int_t j)
{// Builder of the final Raw Cluster (before deconvolution)  
  Info("FormRawCluster","Start with digit(%i,%i)",i,j);
  
  fCurrentCluster.AddDigit((AliRICHdigit*) fHitMap->GetHit(i,j));
  fHitMap->FlagHit(i,j);// Flag hit as taken  

  Int_t listX[4], listY[4];    //  Now look recursively for all neighbours
  for (Int_t iNeighbour=0;iNeighbour<Rich()->Param()->PadNeighbours(i,j,listX,listY);iNeighbour++)
    if(fHitMap->TestHit(listX[iNeighbour],listY[iNeighbour])==kUnused) 
                      FormRawCluster(listX[iNeighbour],listY[iNeighbour]);    
}//AddDigit2Cluster()
//__________________________________________________________________________________________________
void AliRICHClusterFinder::ResolveCluster()
{// Decluster algorithm
  Info("ResolveCluster","Start.");    

  fCurrentCluster.CoG(); // first initial approxmation of the CoG...to start minimization.
  fCurrentCluster.Print();
  switch (fCurrentCluster.Size()) {
  
  case 1:
    WriteRawCluster(); break;
  case 2:
    FitCoG();
    WriteRawCluster(); break;
    
  default:
    WriteRawCluster(); break;
  }     
}//ResolveCluster()
//__________________________________________________________________________________________________
void AliRICHClusterFinder::WriteRawCluster()
{// out the current RawCluster
  Info("WriteRawCluster","Start.");
  
  Rich()->AddCluster(fCurrentCluster);
  
}//WriteRawCluster()
//__________________________________________________________________________________________________
void AliRICHClusterFinder::FitCoG()
{// Fit cluster size 2 by single Mathieson
  Info("FitCoG","Start.");
  
  TMinuit *pMinuit = new TMinuit(2);
  
  Double_t arglist;
  Int_t ierflag = 0;

  static Double_t vstart[2];
  static Double_t lower[2], upper[2];
  static Double_t step[2]={0.001,0.001};

  TString chname;
  Int_t ierflg;
  
  pMinuit->SetObjectFit((TObject*)this);
  pMinuit->SetFCN(RICHMinMathieson);
  pMinuit->mninit(5,10,7);

  vstart[0] = fCurrentCluster.X();
  vstart[1] = fCurrentCluster.Y();

  lower[0] = vstart[0] - 2*AliRICHParam::PadSizeX();
  upper[0] = vstart[0] + 2*AliRICHParam::PadSizeX();
  lower[1] = vstart[1] - 2*AliRICHParam::PadSizeY();
  upper[1] = vstart[1] + 2*AliRICHParam::PadSizeY();


  pMinuit->mnparm(0," x position ",vstart[0],step[0],lower[0],upper[0],ierflag);
  pMinuit->mnparm(1," y position ",vstart[1],step[1],lower[1],upper[1],ierflag);

  arglist = -1;

  pMinuit->SetPrintLevel(-1);
  pMinuit->mnexcm("SET NOGR",&arglist, 1, ierflag);
  pMinuit->mnexcm("SET NOW",&arglist, 1, ierflag);
  arglist = 1;
  pMinuit->mnexcm("SET ERR", &arglist, 1,ierflg);
  arglist = -1;
  pMinuit->mnexcm("SIMPLEX",&arglist, 0, ierflag);
  pMinuit->mnexcm("MIGRAD",&arglist, 0, ierflag);
  pMinuit->mnexcm("EXIT" ,&arglist, 0, ierflag);  
  Double_t xCoG,yCoG;
  Double_t eps, b1, b2;
  pMinuit->mnpout(0,chname, xCoG, eps , b1, b2, ierflg);
  pMinuit->mnpout(1,chname, yCoG, eps , b1, b2, ierflg);
  delete pMinuit;
}
//__________________________________________________________________________________________________
void RICHMinMathieson(Int_t &, Double_t *, Double_t &chi2, Double_t *par, Int_t iflag)
{// Minimization function of Mathieson
  
  AliRICHcluster *pRawCluster = ((AliRICHClusterFinder*)gMinuit->GetObjectFit())->GetCurrentCluster();

  TVector3 centroid(par[0],par[1],0);
  
  chi2 = 0;
  Int_t qtot = pRawCluster->Q();
  for(Int_t i=0;i<pRawCluster->Size();i++) {
    Int_t    padX = ((AliRICHdigit *)pRawCluster->Digits()->At(i))->X();
    Int_t    padY = ((AliRICHdigit *)pRawCluster->Digits()->At(i))->Y();
    Double_t padQ = ((AliRICHdigit *)pRawCluster->Digits()->At(i))->Q();
     chi2 += TMath::Power((qtot*AliRICHParam::Loc2PadFrac(centroid,padX,padY)-padQ),2)/padQ;
  }   
  if(iflag == 3)
    {
            cout << " --- end convergence...summary --- " << endl;
            cout << " x position " << par[0] << endl;
            cout << " y position " << par[1] << endl;
            cout << " chi2       " << chi2 << endl;
    }
}//RICHMinMathieson()
