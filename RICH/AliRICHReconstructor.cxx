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

#include "AliRICHReconstructor.h" //class header
#include <AliRawReader.h>         //Reconstruct(...)
#include "AliRICHv1.h"         //Reconstruct(...) 
#include <AliRun.h>               //ConvertDigits uses gAlice
#include <TMinuit.h>              //Dig2Clu()
ClassImp(AliRICHReconstructor)

void AliRICHReconstructor::Reconstruct(AliRunLoader *pAL,AliRawReader* pRR)const
{
 
  AliLoader *pRL=pAL->GetDetectorLoader("RICH");
  AliRICH *pRich=(AliRICH*)pAL->GetAliRun()->GetDetector("RICH");
  
  AliRICHDigit dig; //tmp digit, raw digit will be converted to it
  TClonesArray *pDigList=new TClonesArray("AliRICHDigit"); Int_t iDigCnt=0; //tmp list of digits for single chamber only
  
  Int_t iEvtN=0;
  while(pRR->NextEvent()){//events loop
    pAL->GetEvent(iEvtN++);
    pRL->MakeTree("R");  pRich->MakeBranch("R");
    
    for(Int_t iChN=1;iChN<=7;iChN++){//chambers loop
      pRR->Select(AliRICHDigit::kRichRawId,2*iChN-2,2*iChN-1);//select only DDL files for the current chamber      
      UInt_t w32=0;
      while(pRR->ReadNextInt(w32)){//raw records loop (in selected DDL files)
        UInt_t ddl=pRR->GetDDLID(); //returns 0,1,2 ... 13
        dig.Raw2Dig(ddl,w32);  
        AliDebug(1,Form("Ch=%i DDL=%i raw=0x%x digit=(%3i,%3i,%3i,%3i) Q=%5.2f",iChN,ddl,w32,dig.Chamber(),dig.Sector(),dig.PadX(),dig.PadY(),dig.Qdc()));
        new((*pDigList)[iDigCnt++]) AliRICHDigit(dig); //add this digit to the tmp list
      }//raw records loop
      if(iDigCnt) Dig2Clu(pDigList,pRich->Clusters(iChN));//cluster finder for the current chamber if any digits present
      pRR->Reset();        
      pDigList->Delete();  iDigCnt=0;//clean up list of digits for the current chamber
    }//chambers loop
    pRL->TreeR()->Fill();            //fill tree for current event
    pRL->WriteRecPoints("OVERWRITE");//write out clusters for current event
    pRich->ClustersReset();
  }//events loop  
  pRL->UnloadRecPoints();  
}//Reconstruct raw data
//__________________________________________________________________________________________________
void AliRICHReconstructor::Dig2Clu(TClonesArray *pDigList,TClonesArray *pCluList)const
{
//Finds all clusters for a given digits list provided not empty. Currently digits list provided is a list of all digits for a single chamber.
//Puts all found clusters in the given clusters list. 
//Arguments: pDigList - list of digits provided not empty
//  Returns: none    
  TMatrixF digMap(1,AliRICHParam::NpadsX(),1,AliRICHParam::NpadsY());  digMap=(Float_t)-1; //digit map for one chamber reseted to -1
  for(Int_t iDigN=0 ; iDigN < pDigList->GetEntriesFast() ; iDigN++){ //digits loop to fill digits map
    AliRICHDigit *pDig= (AliRICHDigit*)pDigList->At(iDigN);
    digMap( pDig->PadX(), pDig->PadY() )=iDigN;                     //(padx,pady) cell takes digit number
  }
  
  AliRICHCluster clu; Int_t iCluCnt=0; //tmp cluster and cluster counter
  
  for(Int_t iDigN=0;iDigN<pDigList->GetEntriesFast();iDigN++){//digits loop    
    AliRICHDigit *pDig=(AliRICHDigit*)pDigList->At(iDigN);
    if(!(pDig=UseDig(pDig->PadX(),pDig->PadY(),pDigList,&digMap))) continue;  //this digit is already taken in FormCluster(), go after next digit
    FormCluster(&clu,pDig,pDigList,&digMap);                            //form cluster starting from this digit
    TMinuit *pMinuit=clu.Solve();                              //solve this cluster
    
    if(pMinuit){//means cluster is solved into local maxima number of clusters, so add all of them in loop
      Double_t x=-1,y=-1,q=-1;TString str; Double_t b1,b2,b3; Int_t iErrFlg;//tmp to withdraw resulting parameters
      for(Int_t i=0;i<clu.Nlocmax();i++){//take resulting parameters 
        pMinuit->mnpout(3*i  ,str, x, b1, b2, b3, iErrFlg);
        pMinuit->mnpout(3*i+1,str, y, b1, b2, b3, iErrFlg);
        pMinuit->mnpout(3*i+2,str, q, b1, b2, b3, iErrFlg);
        clu.Set(x,y,(Int_t)q);
        new((*pCluList)[iCluCnt++]) AliRICHCluster(clu);
      }
      delete pMinuit;
    }else//means cluster is solved as simple center of gravity cluster, add it
      new((*pCluList)[iCluCnt++]) AliRICHCluster(clu); 
    clu.Reset();//make current cluster empty
  }//digits loop
}//Dig2Clu()
//__________________________________________________________________________________________________
void  AliRICHReconstructor::FormCluster(AliRICHCluster *pClu,AliRICHDigit *pDig,TClonesArray *pDigList,TMatrixF *pDigMap)const
{
//Forms the initial cluster as a sum of all adjascent digits. Starts from the given digit
//then calls itself recursevly  for all neighbours.
//Arguments: pClu - pointer to cluster being formed
//  Returns: none
  pClu->AddDigit(pDig);//take this digit in cluster and mark it as taken

  Int_t x[4],y[4];
  
  Int_t iNnei=AliRICHParam::PadNeighbours(pDig->PadX(),pDig->PadY(),x,y);//returns in x,y all possible neighbours of the given one
  for (Int_t i=0;i<iNnei;i++)
    if((pDig=UseDig(x[i],y[i],pDigList,pDigMap))) FormCluster(pClu,pDig,pDigList,pDigMap);    
}//FormCluster()
