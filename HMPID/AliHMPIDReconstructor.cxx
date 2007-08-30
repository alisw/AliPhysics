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
//.
// HMPID base class to reconstruct an event
//.
//.
//.
#include "AliHMPIDReconstructor.h" //class header
#include "AliHMPID.h"              //Reconstruct() 
#include "AliHMPIDCluster.h"       //Dig2Clu()
#include "AliHMPIDParam.h"         //FillEsd() 
#include <AliCDBEntry.h>           //ctor
#include <AliCDBManager.h>         //ctor
#include <AliESDEvent.h>           //FillEsd()
#include <AliRawReader.h>          //Reconstruct() for raw digits
ClassImp(AliHMPIDReconstructor)

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliHMPIDReconstructor::AliHMPIDReconstructor():AliReconstructor(),fUserCut(0),fDaqSig(0),fDig(0),fClu(0)
{
//
//ctor
//
  AliHMPIDParam::Instance();                                                        //geometry loaded for reconstruction
  fUserCut = new Int_t[7];
  fClu=new TObjArray(AliHMPIDParam::kMaxCh+1); fClu->SetOwner(kTRUE);
  fDig=new TObjArray(AliHMPIDParam::kMaxCh+1); fDig->SetOwner(kTRUE);
  
  for(int i=AliHMPIDParam::kMinCh;i<=AliHMPIDParam::kMaxCh;i++){ 
    fDig->AddAt(new TClonesArray("AliHMPIDDigit"),i);
    TClonesArray *pClus = new TClonesArray("AliHMPIDCluster");
    pClus->SetUniqueID(i);
    fClu->AddAt(pClus,i);
  }
  
  AliCDBEntry *pUserCutEnt =AliCDBManager::Instance()->Get("HMPID/Calib/UserCut");    //contains TObjArray of 14 TObject with n. of sigmas to cut charge 
  if(pUserCutEnt) {
    TObjArray *pUserCut = (TObjArray*)pUserCutEnt->GetObject(); 
    for(Int_t iCh=AliHMPIDParam::kMinCh;iCh<=AliHMPIDParam::kMaxCh;iCh++){                  //chambers loop 
      fUserCut[iCh] = pUserCut->At(iCh)->GetUniqueID();
      Printf("HMPID: UserCut successfully loaded for chamber %i -> %i ",iCh,fUserCut[iCh]);
    }
  }

  
  AliCDBEntry *pDaqSigEnt =AliCDBManager::Instance()->Get("HMPID/Calib/DaqSig");  //contains TObjArray of TObjArray 14 TMatrixF sigmas values for pads 
  if(!pDaqSigEnt) AliFatal("No pedestals from DAQ!");
  fDaqSig = (TObjArray*)pDaqSigEnt->GetObject();
  for(Int_t iCh=AliHMPIDParam::kMinCh;iCh<=AliHMPIDParam::kMaxCh;iCh++){                  //chambers loop 
    Printf(" HMPID: DaqSigCut successfully loaded for chamber %i -> %i ",iCh,(Int_t)fDaqSig->At(iCh)->GetUniqueID());
  }
}//AliHMPIDReconstructor
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDReconstructor::Dig2Clu(TObjArray *pDigAll,TObjArray *pCluAll,Bool_t isTryUnfold)
{
// Finds all clusters for a given digits list provided not empty. Currently digits list is a list of all digits for a single chamber.
// Puts all found clusters in separate lists, one per clusters. 
// Arguments: pDigAll     - list of digits for all chambers 
//            pCluAll     - list of clusters for all chambers
//            isTryUnfold - flag to choose between CoG and Mathieson fitting  
//  Returns: none    
  TMatrixF padMap(AliHMPIDParam::kMinPx,AliHMPIDParam::kMaxPcx,AliHMPIDParam::kMinPy,AliHMPIDParam::kMaxPcy);  //pads map for single chamber 0..159 x 0..143 
  
  for(Int_t iCh=AliHMPIDParam::kMinCh;iCh<=AliHMPIDParam::kMaxCh;iCh++){                  //chambers loop 
    TClonesArray *pDigCur=(TClonesArray*)pDigAll->At(iCh);                                //get list of digits for current chamber
    if(pDigCur->GetEntriesFast()==0) continue;                                            //no digits for current chamber
  
    padMap=(Float_t)-1;                                                                   //reset map to -1 (means no digit for this pad)  
    TClonesArray *pCluCur=(TClonesArray*)pCluAll->At(iCh);                                //get list of clusters for current chamber
    
    for(Int_t iDig=0;iDig<pDigCur->GetEntriesFast();iDig++){                              //digits loop to fill pads map
      AliHMPIDDigit *pDig= (AliHMPIDDigit*)pDigCur->At(iDig);                             //get current digit
      padMap( pDig->PadChX(), pDig->PadChY() )=iDig;                                      //fill the map, (padx,pady) cell takes digit index
    }//digits loop to fill digits map 
    
    AliHMPIDCluster clu;                                                                  //tmp cluster to be used as current
    for(Int_t iDig=0;iDig<pDigCur->GetEntriesFast();iDig++){                              //digits loop for current chamber
      AliHMPIDDigit *pDig=(AliHMPIDDigit*)pDigCur->At(iDig);                              //take current digit
      if(!(pDig=UseDig(pDig->PadChX(),pDig->PadChY(),pDigCur,&padMap))) continue;         //this digit is already taken in FormClu(), go after next digit
      FormClu(&clu,pDig,pDigCur,&padMap);                                                 //form cluster starting from this digit by recursion
      
      clu.Solve(pCluCur,isTryUnfold);                                                     //solve this cluster and add all unfolded clusters to provided list  
      clu.Reset();                                                                        //empty current cluster
    }//digits loop for current chamber
  }//chambers loop
}//Dig2Clu()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void  AliHMPIDReconstructor::FormClu(AliHMPIDCluster *pClu,AliHMPIDDigit *pDig,TClonesArray *pDigLst,TMatrixF *pDigMap)
{
// Forms the initial cluster as a combination of all adjascent digits. Starts from the given digit then calls itself recursevly  for all neighbours.
// Arguments: pClu - pointer to cluster being formed
//   Returns: none   
  pClu->DigAdd(pDig);//take this digit in cluster

  Int_t cnt=0,cx[4],cy[4];
  
  if(pDig->PadPcX() != AliHMPIDParam::kMinPx){cx[cnt]=pDig->PadChX()-1; cy[cnt]=pDig->PadChY()  ;cnt++;}       //left
  if(pDig->PadPcX() != AliHMPIDParam::kMaxPx){cx[cnt]=pDig->PadChX()+1; cy[cnt]=pDig->PadChY()  ;cnt++;}       //right
  if(pDig->PadPcY() != AliHMPIDParam::kMinPy){cx[cnt]=pDig->PadChX()  ; cy[cnt]=pDig->PadChY()-1;cnt++;}       //down
  if(pDig->PadPcY() != AliHMPIDParam::kMaxPy){cx[cnt]=pDig->PadChX()  ; cy[cnt]=pDig->PadChY()+1;cnt++;}       //up
  
  for (Int_t i=0;i<cnt;i++){//neighbours loop
    if((pDig=UseDig(cx[i],cy[i],pDigLst,pDigMap))) FormClu(pClu,pDig,pDigLst,pDigMap);   //check if this neighbour pad fired and mark it as taken  
  }//neighbours loop  
}//FormClu()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDReconstructor::Reconstruct(TTree *pDigTree,TTree *pCluTree)const
{
//Invoked  by AliReconstruction to convert digits to clusters i.e. reconstruct simulated data
//Arguments: pDigTree - pointer to Digit tree
//           pCluTree - poitner to Cluster tree
//  Returns: none    
  AliDebug(1,"Start.");
  for(Int_t iCh=AliHMPIDParam::kMinCh;iCh<=AliHMPIDParam::kMaxCh;iCh++) {
    pCluTree->Branch(Form("HMPID%d",iCh),&((*fClu)[iCh]),4000,0);
    pDigTree->SetBranchAddress(Form("HMPID%d",iCh),&((*fDig)[iCh]));
  }   
  pDigTree->GetEntry(0);
  Dig2Clu(fDig,fClu);              //cluster finder 
  pCluTree->Fill();                //fill tree for current event
  
  for(Int_t iCh=AliHMPIDParam::kMinCh;iCh<=AliHMPIDParam::kMaxCh;iCh++){
    fDig->At(iCh)->Clear();
    fClu->At(iCh)->Clear();
  }
  
  AliDebug(1,"Stop.");      
}//Reconstruct(for simulated digits)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDReconstructor::Reconstruct(AliRunLoader *pAL,AliRawReader* pRR)const
{
//Invoked  by AliReconstruction to convert raw digits from DDL files to clusters
//Arguments: pAL - ALICE run loader pointer
//           pRR - ALICE raw reader pointer
//  Returns: none
  AliLoader *pRL=pAL->GetDetectorLoader("HMPID");  AliHMPID *pRich=(AliHMPID*)pAL->GetAliRun()->GetDetector("HMPID");//get pointers for HMPID and HMPID loader

  AliHMPIDDigit dig; //tmp digit, raw digit will be converted to it

  TObjArray digLst; Int_t iDigCnt[7]; for(Int_t i=0;i<7;i++){digLst.AddAt(new TClonesArray("AliHMPIDDigit"),i); iDigCnt[i]=0;} //tmp list of digits for allchambers

  Int_t iEvtN=0;
  while(pRR->NextEvent()){//events loop
    pAL->GetEvent(iEvtN++);
    pRL->MakeTree("R");  pRich->MakeBranch("R");

    for(Int_t iCh=0;iCh<7;iCh++){//chambers loop
      pRR->Select("HMPID",2*iCh,2*iCh+1);//select only DDL files for the current chamber
      UInt_t w32=0;
      while(pRR->ReadNextInt(w32)){//raw records loop (in selected DDL files)
        UInt_t ddl=pRR->GetDDLID(); //returns 0,1,2 ... 13
        dig.Raw(w32,ddl,pRR);
        AliDebug(1,Form("Ch=%i DDL=%i raw=0x%x digit=(%3i,%3i,%3i,%3i) Q=%5.2f",iCh,ddl,w32,dig.Ch(),dig.Pc(),dig.PadPcX(),dig.PadPcY(),dig.Q()));
        if(!IsDigSurvive(&dig)) continue;                                             //sigma cut test
        new((*((TClonesArray*)digLst.At(iCh)))[iDigCnt[iCh]++]) AliHMPIDDigit(dig); //add this digit to the tmp list
      }//raw records loop
      pRR->Reset();
    }//chambers loop

    Dig2Clu(&digLst,pRich->CluLst());//cluster finder for all chambers
    for(Int_t i=0;i<7;i++){digLst.At(i)->Delete(); iDigCnt[i]=0;}                    //clean up list of digits for all chambers

    pRL->TreeR()->Fill();            //fill tree for current event
    pRL->WriteRecPoints("OVERWRITE");//write out clusters for current event
    pRich->CluReset();
  }//events loop
  pRL->UnloadRecPoints();
}//Reconstruct raw data
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDReconstructor::ConvertDigits(AliRawReader *pRR,TTree *pDigTree)const
{
//Invoked  by AliReconstruction to convert raw digits from DDL files to digits
//Arguments: pRR - ALICE raw reader pointer
//           pDigTree - pointer to Digit tree
//  Returns: none    
  AliDebug(1,"Start.");
  AliHMPIDDigit dig; //tmp digit, raw digit will be converted to it
  for(Int_t iCh=AliHMPIDParam::kMinCh;iCh<=AliHMPIDParam::kMaxCh;iCh++) {
    pDigTree->Branch(Form("HMPID%d",iCh),&((*fDig)[iCh]),4000,0);
    pRR->Select("HMPID",2*iCh,2*iCh+1);//select only DDL files for the current chamber      
    Int_t iDigCnt=0;
    UInt_t w32=0;
                while(pRR->ReadNextInt(w32)){//raw records loop (in selected DDL files)
      UInt_t ddl=pRR->GetDDLID(); //returns 0,1,2 ... 13
      if(!dig.Raw(w32,ddl,pRR)) continue;  
      AliDebug(1,Form("Ch=%i DDL=%i raw=0x%x digit=(%3i,%3i,%3i,%3i) Q=%5.2f",iCh,ddl,w32,dig.Ch(),dig.Pc(),dig.PadPcX(),dig.PadPcY(),dig.Q()));
      if(!IsDigSurvive(&dig)) continue;                                       //sigma cut test
      new((*((TClonesArray*)fDig->At(iCh)))[iDigCnt++]) AliHMPIDDigit(dig); //add this digit to the tmp list
    }//raw records loop
    pRR->Reset();        
  }//chambers loop
  pDigTree->Fill();
  for(Int_t iCh=AliHMPIDParam::kMinCh;iCh<=AliHMPIDParam::kMaxCh;iCh++)fDig->At(iCh)->Clear();
  AliDebug(1,"Stop.");
}//Reconstruct digits from raw digits
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDReconstructor::FillESD(AliRunLoader *, AliESDEvent *pESD) const
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
      if(cosThetaTh<1) //calculate the height of theoretical theta ckov on the gaus of experimental one
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
}//FillESD()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
