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

#include "AliHMPIDReconstructor.h" //class header
#include "AliHMPID.h"              //Reconstruct() 
#include "AliHMPIDCluster.h"       //CluQA() 
#include "AliHMPIDParam.h"         //FillEsd() 
#include <AliESD.h>               //FillEsd()
#include <AliRunLoader.h>         //Reconstruct() for simulated digits
#include <AliRawReader.h>         //Reconstruct() for raw digits
#include <AliRun.h>               //Reconstruct()
#include <TH1F.h>                 //CluQA() 
#include <TH2F.h>                 //CluQA() 
#include <TCanvas.h>              //CluQA() 
#include <TNtupleD.h>             //CheckPR()
ClassImp(AliHMPIDReconstructor)

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDReconstructor::CluQA(AliRunLoader *pAL)
{
// Quality assesment plots for clusters. 
// This methode takes list of digits and form list of clusters again in order to 
// calculate cluster shape and cluster particle mixture    
  AliLoader *pRL=pAL->GetDetectorLoader("HMPID");  AliHMPID *pRich=(AliHMPID*)pAL->GetAliRun()->GetDetector("HMPID");//get pointers for HMPID and HMPID loader
  Int_t iNevt=pAL->GetNumberOfEvents();  if(iNevt==0)             {AliInfoClass("No events");return;}   
                                         if(pRL->LoadDigits())    {AliInfoClass("No digits file");return;}
                                            pAL->LoadHeader();
                                            pAL->LoadKinematics(); 
//                                            AliStack *pStack=pAL->Stack();
  TH1::AddDirectory(kFALSE);
  
        
  TH1F*    pQ=new TH1F("RiAllQ"  ,"Charge All"           ,4000 ,0  ,4000);// Q hists
  TH1F* pCerQ=new TH1F("RiCerQ"  ,"Charge Ckov"          ,4000 ,0  ,4000);
  TH1F* pMipQ=new TH1F("RiMipQ"  ,"Charge MIP"           ,4000 ,0  ,4000);
  
  TH1F*    pS=new TH1F("RichCluSize"    ,"Cluster size;size"         ,100  ,0  ,100 );// size hists
  TH1F* pCerS=new TH1F("RichCluCerSize" ,"Ckov size;size"            ,100  ,0  ,100 );
  TH1F* pMipS=new TH1F("RichCluMipSize" ,"MIP size;size"             ,100  ,0  ,100 );
  
  TH2F*    pM=new TH2F("RichCluMap"     ,"Cluster map;x [cm];y [cm]" ,1000 ,0  ,AliHMPIDDigit::SizePcX(),1000,0,AliHMPIDDigit::SizePcY()); // maps
  TH2F* pMipM=new TH2F("RichCluMipMap"  ,"MIP map;x [cm];y [cm]"     ,1000 ,0  ,AliHMPIDDigit::SizePcX(),1000,0,AliHMPIDDigit::SizePcY());
  TH2F* pCerM=new TH2F("RichCluCerMap"  ,"Ckov map;x [cm];y [cm]"    ,1000 ,0  ,AliHMPIDDigit::SizePcX(),1000,0,AliHMPIDDigit::SizePcY());
 
  
  
  for(Int_t iEvt=0;iEvt<iNevt;iEvt++){
    pAL->GetEvent(iEvt);               
    pRL->TreeD()->GetEntry(0); 
    TClonesArray *pCluLst=new TClonesArray("AliHMPIDCluster");//tmp list of clusters for this event
    
    for(Int_t iCh=0;iCh<7;iCh++) Dig2Clu(pRich->DigLst(iCh),pCluLst,kFALSE);//cluster finder for all chamber if any digits present
    
    for(Int_t iClu=0;iClu<pCluLst->GetEntriesFast();iClu++){
      AliHMPIDCluster *pClu = (AliHMPIDCluster*)pCluLst->At(iClu);
      Int_t cfm=0; for(Int_t iDig=0;iDig<pClu->Size();iDig++)  cfm+=pClu->Dig(iDig)->Ch(); //collect ckov-fee-mip structure of current cluster ?????
      Int_t iNckov=cfm/1000000;      Int_t iNfee =cfm%1000000/1000;      Int_t iNmip =cfm%1000000%1000; 

                                             pQ   ->Fill(pClu->Q()) ; pS   ->Fill(pClu->Size()) ; pM    ->Fill(pClu->X(),pClu->Y()); //all clusters                                      
      if(iNckov!=0 && iNfee==0 && iNmip==0) {pCerQ->Fill(pClu->Q()) ; pCerS->Fill(pClu->Size()) ; pCerM ->Fill(pClu->X(),pClu->Y());}//ckov only cluster
      if(iNckov==0 && iNfee==0 && iNmip!=0) {pMipQ->Fill(pClu->Q()) ; pMipS->Fill(pClu->Size()) ; pMipM ->Fill(pClu->X(),pClu->Y());}//mip only cluster
                                       
    }//clusters loop   
    pCluLst->Clear();delete pCluLst;
  }//events loop
  
  pRL->UnloadDigits(); pAL->UnloadKinematics(); pAL->UnloadHeader();
  TCanvas *pC=new TCanvas("RichCluQA",Form("QA for cluster from %i events",iNevt),1000,900); pC->Divide(3,3);
  pC->cd(1);    pM->Draw();          pC->cd(2);    pQ->Draw();       pC->cd(3);    pS->Draw();        
  pC->cd(4); pMipM->Draw();          pC->cd(5); pMipQ->Draw();       pC->cd(6); pMipS->Draw();        
  pC->cd(7); pCerM->Draw();          pC->cd(8); pCerQ->Draw();       pC->cd(9); pCerS->Draw();        
}//CluQA()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDReconstructor::Dig2Clu(TClonesArray *pDigLst,TClonesArray *pCluLst,Bool_t isTryUnfold)
{
// Finds all clusters for a given digits list provided not empty. Currently digits list is a list of all digits for a single chamber.
// Puts all found clusters in separate lists, one per clusters. 
// Arguments: pDigLst     - list of digits provided not empty
//            pCluLst     - list of clusters, provided empty     
//            isTryUnfold - flag to choose between CoG and Mathieson fitting  
//  Returns: none    
  TMatrixF digMap(AliHMPIDDigit::kPadAllX,AliHMPIDDigit::kPadAllY);  digMap=(Float_t)-1;   //digit map for single chamber reseted to -1
  for(Int_t iDig=0;iDig<pDigLst->GetEntriesFast();iDig++){                                 //digits loop to fill digits map
    AliHMPIDDigit *pDig= (AliHMPIDDigit*)pDigLst->At(iDig);                                //get current digit
    digMap( pDig->PadChX(), pDig->PadChY() )=iDig;                                             //fill the map, (padx,pady) cell takes digit index
  }                                                                                        //digits loop to fill digits map 
  
  AliHMPIDCluster clu;                                                                      //tmp cluster to be used as current
  
  for(Int_t iDig=0;iDig<pDigLst->GetEntriesFast();iDig++){                                 //digits loop to form clusters list
    AliHMPIDDigit *pDig=(AliHMPIDDigit*)pDigLst->At(iDig);                                 //take current digit
    if(!(pDig=UseDig(pDig->PadChX(),pDig->PadChY(),pDigLst,&digMap))) continue;            //this digit is already taken in FormClu(), go after next digit
    FormClu(&clu,pDig,pDigLst,&digMap);                                                    //form cluster starting from this digit by recursion
    clu.Solve(pCluLst,isTryUnfold);                                                        //solve this cluster and add all unfolded clusters to provided list  
    clu.Reset();                                                                           //empty current cluster
  }                                                                                        //digits loop to form clusters list
}//Dig2Clu()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void  AliHMPIDReconstructor::FormClu(AliHMPIDCluster *pClu,AliHMPIDDigit *pDig,TClonesArray *pDigLst,TMatrixF *pDigMap)
{
// Forms the initial cluster as a combination of all adjascent digits. Starts from the given digit
// then calls itself recursevly  for all possible neighbours.
// Arguments: pClu - pointer to cluster being formed
//  Returns: none   ???????????????????
  pClu->DigAdd(pDig);//take this digit in cluster

  Int_t x[4],y[4];
  
  Int_t cnt=0;  Int_t iPadX=pDig->PadPcX(); Int_t iPadY=pDig->PadPcY();
  if(iPadX != 0)                        {x[cnt]=iPadX-1; y[cnt]=iPadY;   cnt++;}       //left
  if(iPadX != AliHMPIDDigit::kPadPcX-1) {x[cnt]=iPadX+1; y[cnt]=iPadY;   cnt++;}       //right
  if(iPadY != 0)                        {x[cnt]=iPadX;   y[cnt]=iPadY-1; cnt++;}       //down
  if(iPadY != AliHMPIDDigit::kPadPcY-1) {x[cnt]=iPadX;   y[cnt]=iPadY+1; cnt++;}       //up
  
  for (Int_t i=0;i<cnt;i++)
    if((pDig=UseDig(x[i],y[i],pDigLst,pDigMap))) FormClu(pClu,pDig,pDigLst,pDigMap);   //check if this neighbour pad fired and mark it as taken  
  
}//FormClu()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDReconstructor::Reconstruct(AliRunLoader *pAL)const
{
//Invoked  by AliReconstruction to convert digits to clusters i.e. reconstruct simulated data
//Arguments: pAL - ALICE run loader pointer
//  Returns: none    
  AliDebug(1,"Start.");
  AliLoader *pRL=pAL->GetDetectorLoader("HMPID");
  AliHMPID *pRich=(AliHMPID*)pAL->GetAliRun()->GetDetector("HMPID");//get pointers for HMPID and HMPID loader
  pRL->LoadDigits();   
  pRL->LoadRecPoints("recreate");
  
  for(Int_t iEvtN=0;iEvtN<pAL->GetNumberOfEvents();iEvtN++){//events loop
    pAL->GetEvent(iEvtN); AliDebug(1,Form("Processing event %i...",iEvtN)); //switch current directory to next event    
    pRL->TreeD()->GetEntry(0);  pRL->MakeTree("R");  pRich->MakeBranch("R");  //load digits to memory  and create branches for clusters              
    
    for(Int_t iCh=0;iCh<7;iCh++) Dig2Clu(pRich->DigLst(iCh),pRich->CluLst(iCh));//cluster finder 
      
    pRL->TreeR()->Fill();            //fill tree for current event
    pRL->WriteRecPoints("OVERWRITE");//write out clusters for current event
    pRich->DigReset(); pRich->CluReset();
  }//events loop  

  pRL->UnloadDigits(); 
  pRL->UnloadRecPoints();  
    
  AliDebug(1,"Stop.");      
}//Reconstruct(for simulated digits)
//__________________________________________________________________________________________________
void AliHMPIDReconstructor::Reconstruct(AliRunLoader *pAL,AliRawReader* pRR)const
{
//Invoked  by AliReconstruction to convert raw digits from DDL files to clusters
//Arguments: pAL - ALICE run loader pointer
//           pRR - ALICE raw reader pointer  
//  Returns: none    
  AliLoader *pRL=pAL->GetDetectorLoader("HMPID");  AliHMPID *pRich=(AliHMPID*)pAL->GetAliRun()->GetDetector("HMPID");//get pointers for HMPID and HMPID loader
  
  AliHMPIDDigit dig; //tmp digit, raw digit will be converted to it
  
  Int_t iEvtN=0;
  while(pRR->NextEvent()){//events loop
    pAL->GetEvent(iEvtN++);
    pRL->MakeTree("R");  pRich->MakeBranch("R");
    
    for(Int_t iCh=0;iCh<7;iCh++){//chambers loop
      TClonesArray *pDigLst=new TClonesArray("AliHMPIDDigit"); Int_t iDigCnt=0; //tmp list of digits for single chamber
      pRR->Select("HMPID",2*iCh,2*iCh+1);//select only DDL files for the current chamber      
      UInt_t w32=0;
      while(pRR->ReadNextInt(w32)){//raw records loop (in selected DDL files)
        UInt_t ddl=pRR->GetDDLID(); //returns 0,1,2 ... 13
        dig.Raw(ddl,w32);  
        AliDebug(1,Form("Ch=%i DDL=%i raw=0x%x digit=(%3i,%3i,%3i,%3i) Q=%5.2f",iCh,ddl,w32,dig.Ch(),dig.Pc(),dig.PadPcX(),dig.PadPcY(),dig.Q()));
        new((*pDigLst)[iDigCnt++]) AliHMPIDDigit(dig); //add this digit to the tmp list
      }//raw records loop
      if(iDigCnt) Dig2Clu(pDigLst,pRich->CluLst(iCh));//cluster finder for the current chamber if any digits present
      pRR->Reset();        
      pDigLst->Delete();  iDigCnt=0;//clean up list of digits for the current chamber
    }//chambers loop
    pRL->TreeR()->Fill();            //fill tree for current event
    pRL->WriteRecPoints("OVERWRITE");//write out clusters for current event
    pRich->CluReset();
  }//events loop  
  pRL->UnloadRecPoints();  
}//Reconstruct raw data
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void AliHMPIDReconstructor::FillESD(AliRunLoader *, AliESD *pESD) const
{
// Calculates probability to be a electron-muon-pion-kaon-proton
// from the given Cerenkov angle and momentum assuming no initial particle composition
// (i.e. apriory probability to be the particle of the given sort is the same for all sorts)

  AliPID ppp; //needed
  Double_t pid[AliPID::kSPECIES],h[AliPID::kSPECIES];
   
  for(Int_t iTrk=0;iTrk<pESD->GetNumberOfTracks();iTrk++){//ESD tracks loop
    AliESDtrack *pTrk = pESD->GetTrack(iTrk);// get next reconstructed track
    if(pTrk->GetHMPIDsignal()<=0){//HMPID does not find anything reasonable for this track, assign 0.2 for all species
      for(Int_t iPart=0;iPart<AliPID::kSPECIES;iPart++) pid[iPart]=1.0/AliPID::kSPECIES;
      pTrk->SetHMPIDpid(pid);
      continue;
    } 
    Double_t pmod = pTrk->GetP();
    Double_t hTot=0;
    for(Int_t iPart=0;iPart<AliPID::kSPECIES;iPart++){
      Double_t mass = AliPID::ParticleMass(iPart);
      Double_t cosThetaTh = TMath::Sqrt(mass*mass+pmod*pmod)/(AliHMPIDParam::Instance()->MeanIdxRad()*pmod);
      if(cosThetaTh<1) //calculate the height of theortical theta ckov on the gaus of experimental one
        h[iPart] =TMath::Gaus(TMath::ACos(cosThetaTh),pTrk->GetHMPIDsignal(),TMath::Sqrt(pTrk->GetHMPIDchi2()),kTRUE);
      
      else             //beta < 1/ref. idx. => no light at all  
        h[iPart] =0 ;       
      hTot    +=h[iPart]; //total height of all theoretical heights for normalization
    }//species loop
     
    Double_t hMin=TMath::Gaus(pTrk->GetHMPIDsignal()-4*TMath::Sqrt(pTrk->GetHMPIDchi2()),pTrk->GetHMPIDsignal(),TMath::Sqrt(pTrk->GetHMPIDchi2()),kTRUE);//5 sigma protection
    
    for(Int_t iPart=0;iPart<AliPID::kSPECIES;iPart++)//species loop to assign probabilities
      if(hTot>hMin)  
        pid[iPart]=h[iPart]/hTot;
      else                               //all theoretical values are far away from experemental one
        pid[iPart]=1.0/AliPID::kSPECIES; 
    pTrk->SetHMPIDpid(pid);
  }//ESD tracks loop
  //last line is to check if the nearest thetacerenkov to the teorethical one is within 5 sigma, otherwise no response (equal prob to every particle
}//FillESD
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
