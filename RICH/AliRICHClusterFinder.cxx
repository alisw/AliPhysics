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
#include "AliRICH.h"
#include "AliRICHMap.h"
#include "AliRICHParam.h"
#include <AliLoader.h>
#include <AliRun.h>


ClassImp(AliRICHClusterFinder)
//__________________________________________________________________________________________________
AliRICHClusterFinder::AliRICHClusterFinder(AliRICH *pRICH)   
{//main ctor
  Info("main ctor","Start.");
  
  fRICH = pRICH;
  
  fHitMap = 0;  
  
}//main ctor
//__________________________________________________________________________________________________
void AliRICHClusterFinder::FindLocalMaxima(AliRICHcluster *pRawCluster)
{// Split the cluster according to the number of maxima inside
  Info("SplitbyLocalMaxima","Start.");
  Int_t Nlocal = 0;
  Int_t localX[100],localY[100];
  for(Int_t iDig1=0;iDig1<pRawCluster->Size();iDig1++) {
    Int_t iNotMax = 0;
    AliRICHdigit *pDig1 = (AliRICHdigit *)pRawCluster->Digits()->At(iDig1);
    Int_t padX1 = pDig1->X();
    Int_t padY1 = pDig1->Y();
    Double_t padQ1 = pDig1->Q();
    for(Int_t iDig2=0;iDig2<pRawCluster->Size();iDig2++) {
      AliRICHdigit *pDig2 = (AliRICHdigit *)pRawCluster->Digits()->At(iDig2);
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

  fHitMap=new AliRICHMap(Rich()->Digits(iChamber));

  AliRICHcluster *pRawCluster;
    
  for(Int_t iDig=0;iDig<nDigits;iDig++){    
    AliRICHdigit *dig=(AliRICHdigit*)Rich()->Digits(iChamber)->At(iDig);
    Int_t i=dig->X();   Int_t j=dig->Y();
    if(fHitMap->TestHit(i,j)==kUsed) continue;
	
    pRawCluster = new AliRICHcluster;
    FormRawCluster(i,j,pRawCluster);
	
    if(AliRICHParam::IsResolveClusters()) {
      ResolveCluster(pRawCluster); // ResolveCluster serialization will happen inside
    } else {
      WriteRawCluster(pRawCluster); // simply output of the RawCluster found without deconvolution
    }    
    delete pRawCluster;
    
  }//digits loop

  delete fHitMap;
  Info("FindRawClusters","Stop.");
  
}//FindRawClusters()
//__________________________________________________________________________________________________
void  AliRICHClusterFinder::FormRawCluster(Int_t i, Int_t j, AliRICHcluster *pCluster)
{// Builder of the final Raw Cluster (before deconvolution)  
  Info("FormRawCluster","Start with digit(%i,%i)",i,j);
  
//  Int_t idx = fHitMap->GetHitIndex(i,j);
  AliRICHdigit* pDigit = (AliRICHdigit*) fHitMap->GetHit(i,j);
  pCluster->AddDigit(pDigit);
  
  fHitMap->FlagHit(i,j);// Flag hit as taken  

  Int_t listX[4], listY[4];    //  Now look recursively for all neighbours
  for (Int_t iNeighbour=0;iNeighbour<Rich()->Param()->PadNeighbours(i,j,listX,listY);iNeighbour++)
    if(fHitMap->TestHit(listX[iNeighbour],listY[iNeighbour])==kUnused) 
                      FormRawCluster(listX[iNeighbour],listY[iNeighbour],pCluster);    
}//AddDigit2Cluster()
//__________________________________________________________________________________________________
void AliRICHClusterFinder::ResolveCluster(AliRICHcluster *pRawCluster)
{// Decluster algorithm
  Info("ResolveCluster","Start.");    
  
  pRawCluster->SetStatus(kRaw);// just dummy to compile...
     
}//ResolveCluster()
//__________________________________________________________________________________________________
void AliRICHClusterFinder::WriteRawCluster(AliRICHcluster *pRawCluster)
{// out the current RawCluster
  Info("ResolveCluster","Start.");
  
  pRawCluster->SetStatus(kRaw);
  
  
}//WriteRawCluster()
//__________________________________________________________________________________________________
