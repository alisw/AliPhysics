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
#include "AliRICH.h"              //Reconstruct(...) 
#include <AliRawReader.h>         //Reconstruct(...)
#include <AliRun.h>               //ConvertDigits uses gAlice
#include <AliESD.h>               //RichAna()
#include <AliStack.h>             //RichAna()
#include <TFile.h>                //RichAna()
#include <TFile.h>                //RichAna()
#include <TParticle.h>            //RichAna()
#include <TH1F.h>                 //CluQA() 
#include <TH2F.h>                 //CluQA() 
#include <TCanvas.h>              //CluQA() 
#include <TNtupleD.h>             //CheckPR()
ClassImp(AliRICHReconstructor)

//__________________________________________________________________________________________________
void AliRICHReconstructor::CluQA(AliRunLoader *pAL)
{
// Quality assesment plots for clusters. 
// This methode takes list of digits and form list of clusters again in order to 
// calculate cluster shape and cluster particle mixture    
  AliLoader *pRL=pAL->GetDetectorLoader("RICH");  AliRICH *pRich=(AliRICH*)pAL->GetAliRun()->GetDetector("RICH");//get pointers for RICH and RICH loader
  Int_t iNevt=pAL->GetNumberOfEvents();  if(iNevt==0)             {AliInfoClass("No events");return;}   
                                         if(pRL->LoadDigits())    {AliInfoClass("No digits file");return;}
                                            pAL->LoadHeader();
                                            pAL->LoadKinematics(); 
//                                            AliStack *pStack=pAL->Stack();
  TH1::AddDirectory(kFALSE);
  
  TH1F *pN   =new TH1F("RichCluPerEvt"  ,"Number of clusters per event;number" ,100                    ,0  ,99                                                        );
  TH1F* pQ   =new TH1F("RichCluQdc"     ,"Cluster QDC;q [QDC]"                 ,4000 ,0  ,4000                                    );
  TH1F* pS   =new TH1F("RichCluSize"    ,"Cluster size;size"                   ,100                    ,0  ,100                                                       );
  TH2F* pM   =new TH2F("RichCluMap"     ,"Cluster map;x [cm];y [cm]"           ,1000                   ,0  ,AliRICHParam::PcSizeX()   ,1000,0,AliRICHParam::PcSizeY() );
        
  TH1F* pMipQ=new TH1F("RichCluMipQdc"  ,"MIP QDC;q [QDC]"                     ,4000 ,0  ,4000                                    );
  TH1F* pMipS=new TH1F("RichCluMipSize" ,"MIP size;size"                       ,100                    ,0  ,100                                                       );
  TH2F* pMipM=new TH2F("RichCluMipMap"  ,"MIP map;x [cm];y [cm]"               ,1000                   ,0  ,AliRICHParam::PcSizeX()   ,1000,0,AliRICHParam::PcSizeY() );
  
  TH1F* pCerQ=new TH1F("RichCluCerQdc"  ,"Ckov QDC;q [QDC]"                    ,4000 ,0  ,4000                                    );
  TH1F* pCerS=new TH1F("RichCluCerSize" ,"Ckov size;size"                      ,100                    ,0  ,100                                                       );
  TH2F* pCerM=new TH2F("RichCluCerMap"  ,"Ckov map;x [cm];y [cm]"              ,1000                   ,0  ,AliRICHParam::PcSizeX()   ,1000,0,AliRICHParam::PcSizeY() );
    
  TH1F* pFeeQ=new TH1F("RichCluFeeQdc"  ,"Fee QDC;q [QDC]"                     ,4000 ,0  ,4000);
  TH1F* pFeeS=new TH1F("RichCluFeeSize" ,"Fee size;size"                       ,100                    ,0  ,100                                                       );
  TH2F* pFeeM=new TH2F("RichCluFeeMap"  ,"Fee map;x [cm];y [cm]"               ,1000                   ,0  ,AliRICHParam::PcSizeX()   ,1000,0,AliRICHParam::PcSizeY() );
 
  
  TClonesArray *pCluLst=new TClonesArray("AliRICHCluster");//tmp list of clusters
  
  for(Int_t iEvtN=0; iEvtN<iNevt; iEvtN++){
    pAL->GetEvent(iEvtN);               
    pRL->TreeD()->GetEntry(0); 
    for(Int_t iChN=1;iChN<=7;iChN++){//chambers loop
      if(pRich->Digs(iChN)->GetEntriesFast()>0) Dig2Clu(pRich->Digs(iChN),pCluLst,kFALSE);//cluster finder for the current chamber if any digits present
    }//chambers loop
    
    pCluLst->Print();
    for(Int_t iCluN=0 ; iCluN < pCluLst->GetEntriesFast() ; iCluN++){
      AliRICHCluster *pClu = (AliRICHCluster*)pCluLst->At(iCluN);
      Int_t iNckov=0,iNfee=0,iNmip=0; 
                                             pQ   ->Fill(pClu->Q()) ; pS   ->Fill(pClu->Size()) ; pM    ->Fill(pClu->X(),pClu->Y()); //all clusters                                      
      if(iNckov!=0 && iNfee==0 && iNmip==0) {pCerQ->Fill(pClu->Q()) ; pCerS->Fill(pClu->Size()) ; pCerM ->Fill(pClu->X(),pClu->Y());}//ckov only cluster
      if(iNckov==0 && iNfee!=0 && iNmip==0) {pFeeQ->Fill(pClu->Q()) ; pFeeS->Fill(pClu->Size()) ; pFeeM ->Fill(pClu->X(),pClu->Y());}//feed only cluster
      if(iNckov==0 && iNfee==0 && iNmip!=0) {pMipQ->Fill(pClu->Q()) ; pMipS->Fill(pClu->Size()) ; pMipM ->Fill(pClu->X(),pClu->Y());}//mip only cluster
                                       
    }//clusters loop   
  }//events loop
  
  pRL->UnloadDigits(); pAL->UnloadKinematics(); pAL->UnloadHeader();
  TCanvas *pC=new TCanvas("RichCluQA",Form("QA for cluster from %i events",iNevt),1000,900); pC->Divide(3,3);
  pC->cd(1); pN->Draw();          pC->cd(2); pQ->Draw();          pC->cd(3); pS->Draw();        
  pC->cd(4); pN->Draw();          pC->cd(5); pMipQ->Draw();       pC->cd(6); pMipS->Draw();        
  pC->cd(7); pN->Draw();          pC->cd(8); pCerQ->Draw();       pC->cd(9); pCerS->Draw();        
}//CluQA()
//__________________________________________________________________________________________________
void AliRICHReconstructor::CheckPR()
{
//Pattern recognition with stack particles
  TFile *pFile = new TFile("$(HOME)/RPR.root","RECREATE","RICH Pattern Recognition");
  TNtupleD *hn = new TNtupleD("hn","ntuple","Pmod:Charge:TrackTheta:TrackPhi:TrackX:TrackY:MinX:MinY:ChargeMIP:ThetaCerenkov:NPhotons:MipIndex:Chamber:Particle");
//  printf("\n\n");
//  printf("Pattern Recognition done for event %5i",0);
  AliRICH *pRich=((AliRICH*)gAlice->GetDetector("RICH"));
  AliMagF * magf = gAlice->Field();
  AliTracker::SetFieldMap(magf,kTRUE);
  for(Int_t iEvtN=0;iEvtN<pRich->GetLoader()->GetRunLoader()->GetNumberOfEvents();iEvtN++) {
    pRich->GetLoader()->GetRunLoader()->GetEvent(iEvtN);
    AliRICHTracker *tr = new AliRICHTracker();
    tr->RecWithStack(hn);
    AliInfoClass(Form("Pattern Recognition done for event %i \b",iEvtN));
//    printf("\b\b\b\b\b%5i",iEvtN+1);
  }
  printf("\n\n");
  pFile->Write();pFile->Close();
}
//__________________________________________________________________________________________________
void AliRICHReconstructor::Dig2Clu(TClonesArray *pDigLst,TClonesArray *pCluLst,Bool_t isTryUnfold)
{
//Finds all clusters for a given digits list provided not empty. Currently digits list is a list of all digits for a single chamber.
//If pStack not 0 then also finds Ckov-Fee-Mip composition for formed clusters.  
//Puts all found clusters in the given clusters list. 
//Arguments: pDigLst     - list of digits provided not empty
//           pCluLst     - list of clusters, provided empty     
//           isTryUnfold - flag to choose between CoG and Mathieson fitting  
//  Returns: none    
  TMatrixF digMap(1,AliRICHParam::NpadsX(),1,AliRICHParam::NpadsY());  digMap=(Float_t)-1; //digit map for one chamber reseted to -1
  for(Int_t iDigN=0 ; iDigN < pDigLst->GetEntriesFast() ; iDigN++){                        //digits loop to fill digits map
    AliRICHDigit *pDig= (AliRICHDigit*)pDigLst->At(iDigN);                                   //get current digit
    digMap( pDig->PadX(), pDig->PadY() )=iDigN;                                              //fill the map, (padx,pady) cell takes digit number
  }                                                                                        //digits loop to fill digits map 
  
  AliRICHCluster clu;                                                                      //tmp cluster to be used as current
  
  for(Int_t iDigN=0;iDigN<pDigLst->GetEntriesFast();iDigN++){                              //digits loop to form clusters list
    AliRICHDigit *pDig=(AliRICHDigit*)pDigLst->At(iDigN);                                  //take current digit
    if(!(pDig=UseDig(pDig->PadX(),pDig->PadY(),pDigLst,&digMap))) continue;                //this digit is already taken in FormClu(), go after next digit
    FormClu(&clu,pDig,pDigLst,&digMap);                                                    //form cluster starting from this digit by recursion
    clu.Solve(pCluLst,isTryUnfold);                                                        //solve this cluster and add all unfolded clusters to provided list  
    clu.Reset();                                                                           //empty current cluster
  }                                                                                        //digits loop to form clusters list
}//Dig2Clu()
//__________________________________________________________________________________________________
void  AliRICHReconstructor::FormClu(AliRICHCluster *pClu,AliRICHDigit *pDig,TClonesArray *pDigLst,TMatrixF *pDigMap)
{
//Forms the initial cluster as a sum of all adjascent digits. Starts from the given digit
//then calls itself recursevly  for all neighbours.
//Arguments: pClu - pointer to cluster being formed
//  Returns: none
  pClu->DigAdd(pDig);//take this digit in cluster

  Int_t x[4],y[4];
  
  Int_t iNnei=AliRICHParam::PadNeighbours(pDig->PadX(),pDig->PadY(),x,y);//returns in x,y all possible neighbours of the given padx,pady
  for (Int_t i=0;i<iNnei;i++)
    if((pDig=UseDig(x[i],y[i],pDigLst,pDigMap))) FormClu(pClu,pDig,pDigLst,pDigMap);   //check if this neighbour hit and mark it as taken  
}//FormClu()
//__________________________________________________________________________________________________
void AliRICHReconstructor::Reconstruct(AliRunLoader *pAL)const
{
//Invoked  by AliReconstruction to convert digits to clusters i.e. reconstruct simulated data
//Arguments: pAL - ALICE run loader pointer
//  Returns: none    
  AliDebug(1,"Start.");
  AliLoader *pRL=pAL->GetDetectorLoader("RICH");  AliRICH *pRich=(AliRICH*)pAL->GetAliRun()->GetDetector("RICH");//get pointers for RICH and RICH loader
  pRL->LoadDigits();   
  
  for(Int_t iEvtN=0;iEvtN<pAL->GetNumberOfEvents();iEvtN++){//events loop
    pAL->GetEvent(iEvtN++); AliDebug(1,Form("Processing event %i...",iEvtN)); //switch current directory to next event    
    pRL->TreeD()->GetEntry(0);  pRL->MakeTree("R");  pRich->MakeBranch("R");  //load digits to memory  and create branches for clusters              
    for(Int_t iChN=1;iChN<=7;iChN++){//chambers loop
      if(pRich->Digs(iChN)->GetEntriesFast()>0) Dig2Clu(pRich->Digs(iChN),pRich->Clus(iChN));//cluster finder for the current chamber if any digits present
    }//chambers loop
    pRL->TreeR()->Fill();            //fill tree for current event
    pRL->WriteRecPoints("OVERWRITE");//write out clusters for current event
    pRich->DigReset(); pRich->CluReset();
  }//events loop  

  pRL->UnloadDigits(); 
  pRL->UnloadRecPoints();  
    
  AliDebug(1,"Stop.");      
}//Reconstruct(for simulated digits)
//__________________________________________________________________________________________________
void AliRICHReconstructor::Reconstruct(AliRunLoader *pAL,AliRawReader* pRR)const
{
//Invoked  by AliReconstruction to convert raw digits from DDL files to clusters
//Arguments: pAL - ALICE run loader pointer
//           pRR - ALICE raw reader pointer  
//  Returns: none    
  AliLoader *pRL=pAL->GetDetectorLoader("RICH");  AliRICH *pRich=(AliRICH*)pAL->GetAliRun()->GetDetector("RICH");//get pointers for RICH and RICH loader
  
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
        AliDebug(1,Form("Ch=%i DDL=%i raw=0x%x digit=(%3i,%3i,%3i,%3i) Q=%5.2f",iChN,ddl,w32,dig.C(),dig.S(),dig.PadX(),dig.PadY(),dig.Qdc()));
        new((*pDigList)[iDigCnt++]) AliRICHDigit(dig); //add this digit to the tmp list
      }//raw records loop
      if(iDigCnt) Dig2Clu(pDigList,pRich->Clus(iChN));//cluster finder for the current chamber if any digits present
      pRR->Reset();        
      pDigList->Delete();  iDigCnt=0;//clean up list of digits for the current chamber
    }//chambers loop
    pRL->TreeR()->Fill();            //fill tree for current event
    pRL->WriteRecPoints("OVERWRITE");//write out clusters for current event
    pRich->CluReset();
  }//events loop  
  pRL->UnloadRecPoints();  
}//Reconstruct raw data
//__________________________________________________________________________________________________
void AliRICHReconstructor::RichAna(Int_t iNevMin,Int_t iNevMax,Bool_t askPatRec)
{
  TFile *pFile=TFile::Open("AliESDs.root","read");
  if(!pFile || !pFile->IsOpen()) {AliInfoClass("ESD file not open.");return;}      //open AliESDs.root                                                                    
  TTree *pTree = (TTree*) pFile->Get("esdTree");
  if(!pTree){AliInfoClass("ESD not found.");return;}                               //get ESD tree  
  AliInfoClass("ESD found. Go ahead!");
  
  AliRICH *pRich=((AliRICH*)gAlice->GetDetector("RICH"));

  AliMagF * magf = gAlice->Field();
  AliTracker::SetFieldMap(magf,kTRUE);
  pRich->GetLoader()->GetRunLoader()->LoadHeader();
  pRich->GetLoader()->GetRunLoader()->LoadKinematics();
  TString var1 = "Pmod:Charge:TrackTheta:TrackPhi:MinX:MinY:ThetaCerenkov:NPhotons:";
  TString var2 = "ChargeMIP:Chamber:TOF:LengthTOF:prob1:prob2:prob3:";
  TString var3 = "ErrPar1:ErrPar2:ErrPar3:Th1:Th2:Th3:nPhotBKG:pdgCode";
  TString varList = var1+var2+var3;
  
  Double_t dx,dy;
  Double_t hnvec[30];

//  TFile *pFileRA = new TFile("$(HOME)/RichAna.root","RECREATE","RICH Pattern Recognition");
  TFile *pFileRA = new TFile("./RichAna.root","RECREATE","RICH Pattern Recognition");
  TNtupleD *hn = new TNtupleD("hn","ntuple",varList);
  if(iNevMin<0) iNevMin=0;
  if(iNevMin>iNevMax) {iNevMin=0;iNevMax=0;}  
  if(iNevMax==0) iNevMax=999999;
  if(pRich->GetLoader()->GetRunLoader()->GetNumberOfEvents()<iNevMax) iNevMax = pRich->GetLoader()->GetRunLoader()->GetNumberOfEvents();
  AliESD *pESD=new AliESD;  pTree->SetBranchAddress("ESD", &pESD);
  for(Int_t iEvtN=iNevMin;iEvtN<iNevMax;iEvtN++) {
    pTree->GetEvent(iEvtN);
    pRich->GetLoader()->GetRunLoader()->GetEvent(iEvtN);
    pRich->GetLoader()->LoadRecPoints();
    pRich->GetLoader()->TreeR()->GetEntry(0);
    AliStack *pStack = pRich->GetLoader()->GetRunLoader()->Stack();
//Pattern recognition started
    if(pESD->GetNumberOfTracks()) {
      Int_t iNtracks=pESD->GetNumberOfTracks();
      AliInfoClass(Form("Start with %i tracks",iNtracks));
      for(Int_t iTrackN=0;iTrackN<iNtracks;iTrackN++){//ESD tracks loop
        if(iTrackN%100==0)AliInfoClass(Form("Track %i to be processed",iTrackN));
        AliRICHTracker *pTrRich = new AliRICHTracker();
        if(askPatRec==kTRUE) pTrRich->PropagateBack(pESD);
        AliESDtrack *pTrack = pESD->GetTrack(iTrackN);// get next reconstructed track
        
        Int_t lab=TMath::Abs(pTrack->GetLabel());
        TParticle *pPart=pStack->Particle(lab);
        Int_t code=pPart->GetPdgCode();

        pTrack->GetRICHdxdy(dx,dy);
        hnvec[0]=pTrack->GetP();
        hnvec[1]=pTrack->GetSign();
//  cout << " Track momentum " << pTrack->GetP() << " charge " << pTrack->GetSign() << endl;
  
        pTrack->GetRICHthetaPhi(hnvec[2],hnvec[3]);
        pTrack->GetRICHdxdy(hnvec[4],hnvec[5]);
        hnvec[6]=pTrack->GetRICHsignal();
        hnvec[7]=pTrack->GetRICHnclusters();
        hnvec[9]=pTrack->GetRICHcluster()/1000000;
        hnvec[8]=pTrack->GetRICHcluster()-hnvec[9]*1000000;
        hnvec[10]=pTrack->GetTOFsignal();
        hnvec[11]=pTrack->GetIntegratedLength();
        Double_t prob[5];
        pTrack->GetRICHpid(prob);
        hnvec[12]=prob[0]+prob[1]+prob[2];
        hnvec[13]=prob[3];
        hnvec[14]=prob[4];
        hnvec[15]=pTrRich->fErrPar[2];
        hnvec[16]=pTrRich->fErrPar[3];
        hnvec[17]=pTrRich->fErrPar[4];
        for(Int_t i=0;i<3;i++) {
          Double_t mass = AliRICHParam::fgMass[i+2];
          Double_t refIndex=AliRICHParam::Instance()->IdxC6F14(AliRICHParam::EckovMean());
          Double_t cosThetaTh = TMath::Sqrt(mass*mass+pTrack->GetP()*pTrack->GetP())/(refIndex*pTrack->GetP());
          hnvec[18+i]=0;
          if(cosThetaTh>=1) continue;
          hnvec[18+i]= TMath::ACos(cosThetaTh);
        }
//        if(askPatRec==kTRUE) hnvec[21]=pTrRich->fnPhotBKG; else hnvec[21]=0;
        hnvec[22]=code;
        hn->Fill(hnvec);
      }
    }
    AliInfoClass(Form("Pattern Recognition done for event %i",iEvtN));
  }
  pFileRA->Write();pFileRA->Close();// close RichAna.root
  delete pESD;  pFile->Close();//close AliESDs.root
}//RichAna()
//__________________________________________________________________________________________________
void AliRICHReconstructor::Test(Bool_t isTryUnfold)
{
// Test the cluster finding algorithm by providing predifined set of digits
// Arguments: none
//   Returns: none  
  TClonesArray *pDigTst=new TClonesArray("AliRICHDigit");       TClonesArray *pCluTst=new TClonesArray("AliRICHCluster");     
  Int_t iDigCnt=0;
  Int_t c,padx,pady,qdc;
//ckov cluster  
  new((*pDigTst)[iDigCnt++]) AliRICHDigit(c=1,padx= 89,pady=13,qdc= 10);
  new((*pDigTst)[iDigCnt++]) AliRICHDigit(c=1,padx= 90,pady=13,qdc=  7);
  new((*pDigTst)[iDigCnt++]) AliRICHDigit(c=1,padx= 90,pady=12,qdc=  6);
  new((*pDigTst)[iDigCnt++]) AliRICHDigit(c=1,padx= 91,pady=12,qdc=  7);
//mip cluster  
  new((*pDigTst)[iDigCnt++]) AliRICHDigit(c=1,padx= 99,pady=21,qdc=  9);
  new((*pDigTst)[iDigCnt++]) AliRICHDigit(c=1,padx= 99,pady=22,qdc= 26);
  new((*pDigTst)[iDigCnt++]) AliRICHDigit(c=1,padx=100,pady=21,qdc= 39);
  new((*pDigTst)[iDigCnt++]) AliRICHDigit(c=1,padx=100,pady=22,qdc=109);
  new((*pDigTst)[iDigCnt++]) AliRICHDigit(c=1,padx=100,pady=23,qdc=  7);
  new((*pDigTst)[iDigCnt++]) AliRICHDigit(c=1,padx=101,pady=22,qdc= 11);
  
  pDigTst->Print("","Initial digit:");
  Dig2Clu(pDigTst,pCluTst,isTryUnfold);   
  pCluTst->Print("","Solved cluster:"); 
  delete pDigTst; delete pCluTst;
}//Test()
//__________________________________________________________________________________________________
