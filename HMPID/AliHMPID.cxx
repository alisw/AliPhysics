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

#include "AliHMPID.h"
#include "AliHMPIDHit.h"   //OccupancyPrint(), HitQa()
#include "AliHMPIDDigit.h" //
#include <TParticle.h>  //SummaryOfEvent(), HitQa()
#include <TBenchmark.h>  //HitQA()
#include <TPDGCode.h>    //HitQA()
#include <AliStack.h>   //SummaryOfEvent(), HitQa()
#include <AliRun.h>     //HitQa() 
#include <AliMC.h>       //ctor
#include <AliHeader.h>
#include <TH1F.h>        //HitQA()
#include <AliLog.h>      //in many methods to print AliInfo 
ClassImp(AliHMPID)    
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliHMPID::AliHMPID(const char *name, const char *title):AliDetector(name,title),fSdi(0),fDig(0),fClu(0)
{
//Named ctor
  AliDebug(1,"Start.");
//AliDetector ctor deals with Hits and Digits (reset them to 0, does not create them)
  HitCreate();          gAlice->GetMCApp()->AddHitList(fHits);
  AliDebug(1,"Stop.");
}//AliHMPID::AliHMPID(const char *name, const char *title)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliHMPID::~AliHMPID()
{
//dtor
  AliDebug(1,"Start.");

  
  if(fHits)      delete fHits;
  if(fDigits)    delete fDigits;
  if(fSdi)       delete fSdi;
  if(fDig)      {fDig->Delete();   delete fDig;}
  if(fClu)      {fClu->Delete();   delete fClu;}
  AliDebug(1,"Stop.");    
}//AliHMPID::~AliHMPID()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPID::MakeBranch(Option_t* option)
{
//Create Tree branches for the HMPID.
  AliDebug(1,Form("Start with option= %s.",option));
    
  const Int_t kBufSize = 4000;
      
  const char *cH = strstr(option,"H");
  const char *cD = strstr(option,"D");
  const char *cR = strstr(option,"R");
  const char *cS = strstr(option,"S");

  if(cH&&         TreeH()){HitCreate();                       MakeBranchInTree(         TreeH(),     "HMPID"     ,&fHits       ,kBufSize,0);}
  if(cS&&fLoader->TreeS()){SdiCreate();                       MakeBranchInTree(fLoader->TreeS(),     "HMPID"     ,&fSdi        ,kBufSize,0);}
  if(cD&&fLoader->TreeD()){DigCreate();for(Int_t i=0;i<7;i++) MakeBranchInTree(fLoader->TreeD(),Form("HMPID%d",i),&((*fDig)[i]),kBufSize,0);}
  if(cR&&fLoader->TreeR()){CluCreate();for(Int_t i=0;i<7;i++) MakeBranchInTree(fLoader->TreeR(),Form("HMPID%d",i),&((*fClu)[i]),kBufSize,0);}   
  
  AliDebug(1,"Stop.");   
}//void AliHMPID::MakeBranch(Option_t* option)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPID::SetTreeAddress()
{
//Set branch address for the Hits and Digits Tree.
  AliDebug(1,"Start.");
  if(fLoader->TreeH() && fLoader->TreeH()->GetBranch("HMPID" )){HitCreate();                      fLoader->TreeH()->SetBranchAddress(     "HMPID"     ,&fHits       );}
  if(fLoader->TreeS() && fLoader->TreeS()->GetBranch("HMPID" )){SdiCreate();                      fLoader->TreeS()->SetBranchAddress(     "HMPID"     ,&fSdi        );}
  if(fLoader->TreeD() && fLoader->TreeD()->GetBranch("HMPID0")){DigCreate(); for(int i=0;i<7;i++) fLoader->TreeD()->SetBranchAddress(Form("HMPID%d",i),&((*fDig)[i]));}
  if(fLoader->TreeR() && fLoader->TreeR()->GetBranch("HMPID0")){CluCreate(); for(int i=0;i<7;i++) fLoader->TreeR()->SetBranchAddress(Form("HMPID%d",i),&((*fClu)[i]));}
  AliDebug(1,"Stop.");
}//void AliHMPID::SetTreeAddress()
//__________________________________________________________________________________________________
// AliHMPIDHit* AliHMPID::Hit(Int_t tid)const
// {
// // Search for the first HMPID hit belonging to the given tid
//   GetLoader()->LoadHits();
//   for(Int_t iPrimN=0;iPrimN<GetLoader()->TreeH()->GetEntries();iPrimN++){//prims loop      
//     GetLoader()->TreeH()->GetEntry(iPrimN);
//     for(Int_t iHitN=0;iHitN<Hits()->GetEntries();iHitN++){
//       AliHMPIDHit *pHit=(AliHMPIDHit*)Hits()->At(iHitN);
//       if(tid==pHit->Track()) {GetLoader()->UnloadHits();return pHit;}
//     }//hits
//   }//prims loop
//   GetLoader()->UnloadHits();
//   return 0;
// }
//__________________________________________________________________________________________________
void AliHMPID::HitPrint(Int_t iEvtN)const
{
//Prints a list of HMPID hits for a given event. Default is event number 0.
  if(GetLoader()->GetRunLoader()->GetEvent(iEvtN)) return;    
  AliInfo(Form("List of HMPID hits for event %i",iEvtN));
  if(GetLoader()->LoadHits()) return;
  
  Int_t iTotalHits=0;
  for(Int_t iPrimN=0;iPrimN<GetLoader()->TreeH()->GetEntries();iPrimN++){//prims loop
    GetLoader()->TreeH()->GetEntry(iPrimN);      
    Hits()->Print();
    iTotalHits+=Hits()->GetEntries();
  }
  GetLoader()->UnloadHits();
  AliInfo(Form("totally %i hits",iTotalHits));
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPID::SdiPrint(Int_t iEvt)const
{
//prints a list of HMPID sdigits  for a given event
  if(GetLoader()->GetRunLoader()->GetEvent(iEvt)) return;    
  Info("PrintSDigits","List of HMPID sdigits for event %i",iEvt);
  if(GetLoader()->LoadSDigits()) return;
  
  GetLoader()->TreeS()->GetEntry(0);
  SdiLst()->Print();
  GetLoader()->UnloadSDigits();
  Printf("totally %i sdigits",SdiLst()->GetEntries());
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPID::DigPrint(Int_t iEvt)const
{
//prints a list of HMPID digits  for a given event
  if(GetLoader()->GetRunLoader()->GetEvent(iEvt)) return;    
  Printf("List of HMPID digits for event %i",iEvt);
  if(GetLoader()->LoadDigits()) return;
  
  GetLoader()->TreeD()->GetEntry(0);
  DigLst()->Print();
  GetLoader()->UnloadDigits();
  Printf("totally %i Digits",DigLst()->GetEntries());
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPID::CluPrint(Int_t iEvtN)const
{
//prints a list of HMPID clusters  for a given event
  Printf("List of HMPID clusters for event %i",iEvtN);
  GetLoader()->GetRunLoader()->GetEvent(iEvtN);    
  if(GetLoader()->LoadRecPoints()) return;
  
  Int_t iCluCnt=0;
  GetLoader()->TreeR()->GetEntry(0);
  for(Int_t iCh=0;iCh<7;iCh++){
    TClonesArray *pCluLst=(TClonesArray*)fClu->At(iCh);    iCluCnt+=pCluLst->GetEntries();    pCluLst->Print();
  }
  GetLoader()->UnloadRecPoints();
  Printf("totally %i clusters for event %i",iCluCnt,iEvtN);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
