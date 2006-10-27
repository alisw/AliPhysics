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

#include "AliRICH.h"
#include "AliRICHHit.h"   //OccupancyPrint(), HitQa()
#include "AliRICHDigit.h" //OccupancyPrint()
#include <TParticle.h>  //SummaryOfEvent(), HitQa()
#include <TBenchmark.h>  //HitQA()
#include <TPDGCode.h>    //HitQA()
#include <AliStack.h>   //OccupancyPrint(), SummaryOfEvent(), HitQa()
#include <AliRun.h>     //HitQa() 
#include <AliMC.h>       //ctor
#include <AliHeader.h>
#include <AliGenEventHeader.h>
#include <AliGenHijingEventHeader.h>
#include <TH1F.h>        //HitQA()
#include <AliLog.h>      //in many methods to print AliInfo 
ClassImp(AliRICH)    
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliRICH::AliRICH(const char *name, const char *title):AliDetector(name,title),fSdi(0),fDig(0),fClu(0)
{
//Named ctor
  AliDebug(1,"Start.");
//AliDetector ctor deals with Hits and Digits (reset them to 0, does not create them)
  HitCreate();          gAlice->GetMCApp()->AddHitList(fHits);
  AliDebug(1,"Stop.");
}//AliRICH::AliRICH(const char *name, const char *title)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliRICH::~AliRICH()
{
//dtor
  AliDebug(1,"Start.");

  
  if(fHits)      delete fHits;
  if(fDigits)    delete fDigits;
  if(fSdi)       delete fSdi;
  if(fDig)      {fDig->Delete();   delete fDig;}
  if(fClu)      {fClu->Delete();   delete fClu;}
  AliDebug(1,"Stop.");    
}//AliRICH::~AliRICH()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliRICH::MakeBranch(Option_t* option)
{
//Create Tree branches for the RICH.
  AliDebug(1,Form("Start with option= %s.",option));
    
  const Int_t kBufSize = 4000;
      
  const char *cH = strstr(option,"H");
  const char *cD = strstr(option,"D");
  const char *cR = strstr(option,"R");
  const char *cS = strstr(option,"S");

  if(cH&&         TreeH()){HitCreate();                       MakeBranchInTree(         TreeH(),     "RICH"     ,&fHits       ,kBufSize,0);}
  if(cS&&fLoader->TreeS()){SdiCreate();                       MakeBranchInTree(fLoader->TreeS(),     "RICH"     ,&fSdi        ,kBufSize,0);}
  if(cD&&fLoader->TreeD()){DigCreate();for(Int_t i=0;i<7;i++) MakeBranchInTree(fLoader->TreeD(),Form("RICH%d",i),&((*fDig)[i]),kBufSize,0);}
  if(cR&&fLoader->TreeR()){CluCreate();for(Int_t i=0;i<7;i++) MakeBranchInTree(fLoader->TreeR(),Form("RICH%d",i),&((*fClu)[i]),kBufSize,0);}   
  
  AliDebug(1,"Stop.");   
}//void AliRICH::MakeBranch(Option_t* option)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliRICH::SetTreeAddress()
{
//Set branch address for the Hits and Digits Tree.
  AliDebug(1,"Start.");
  if(fLoader->TreeH() && fLoader->TreeH()->GetBranch("RICH" )){HitCreate();                      fLoader->TreeH()->SetBranchAddress(     "RICH"     ,&fHits       );}
  if(fLoader->TreeS() && fLoader->TreeS()->GetBranch("RICH" )){SdiCreate();                      fLoader->TreeS()->SetBranchAddress(     "RICH"     ,&fSdi        );}
  if(fLoader->TreeD() && fLoader->TreeD()->GetBranch("RICH0")){DigCreate(); for(int i=0;i<7;i++) fLoader->TreeD()->SetBranchAddress(Form("RICH%d",i),&((*fDig)[i]));}
  if(fLoader->TreeR() && fLoader->TreeR()->GetBranch("RICH0")){CluCreate(); for(int i=0;i<7;i++) fLoader->TreeR()->SetBranchAddress(Form("RICH%d",i),&((*fClu)[i]));}
  AliDebug(1,"Stop.");
}//void AliRICH::SetTreeAddress()
//__________________________________________________________________________________________________
// AliRICHHit* AliRICH::Hit(Int_t tid)const
// {
// // Search for the first RICH hit belonging to the given tid
//   GetLoader()->LoadHits();
//   for(Int_t iPrimN=0;iPrimN<GetLoader()->TreeH()->GetEntries();iPrimN++){//prims loop      
//     GetLoader()->TreeH()->GetEntry(iPrimN);
//     for(Int_t iHitN=0;iHitN<Hits()->GetEntries();iHitN++){
//       AliRICHHit *pHit=(AliRICHHit*)Hits()->At(iHitN);
//       if(tid==pHit->Track()) {GetLoader()->UnloadHits();return pHit;}
//     }//hits
//   }//prims loop
//   GetLoader()->UnloadHits();
//   return 0;
// }
//__________________________________________________________________________________________________
void AliRICH::HitPrint(Int_t iEvtN)const
{
//Prints a list of RICH hits for a given event. Default is event number 0.
  if(GetLoader()->GetRunLoader()->GetEvent(iEvtN)) return;    
  AliInfo(Form("List of RICH hits for event %i",iEvtN));
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
void AliRICH::SdiPrint(Int_t iEvt)const
{
//prints a list of RICH sdigits  for a given event
  if(GetLoader()->GetRunLoader()->GetEvent(iEvt)) return;    
  Info("PrintSDigits","List of RICH sdigits for event %i",iEvt);
  if(GetLoader()->LoadSDigits()) return;
  
  GetLoader()->TreeS()->GetEntry(0);
  SdiLst()->Print();
  GetLoader()->UnloadSDigits();
  Printf("totally %i sdigits",SdiLst()->GetEntries());
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliRICH::DigPrint(Int_t iEvt)const
{
//prints a list of RICH digits  for a given event
  if(GetLoader()->GetRunLoader()->GetEvent(iEvt)) return;    
  Printf("List of RICH digits for event %i",iEvt);
  if(GetLoader()->LoadDigits()) return;
  
  GetLoader()->TreeD()->GetEntry(0);
  DigLst()->Print();
  GetLoader()->UnloadDigits();
  Printf("totally %i Digits",DigLst()->GetEntries());
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliRICH::OccupancyPrint(Int_t iEvtNreq)
{
//prints occupancy for each chamber in a given event
  Int_t iEvtNmin,iEvtNmax;
  if(iEvtNreq==-1){
    iEvtNmin=0;
    iEvtNmax=gAlice->GetEventsPerRun();
  } else { 
    iEvtNmin=iEvtNreq;iEvtNmax=iEvtNreq+1;
  }
    
  if(GetLoader()->GetRunLoader()->LoadHeader()) return;    
  if(GetLoader()->GetRunLoader()->LoadKinematics()) return;    
  
//  Info("Occupancy","for event %i",iEvtN);
  if(GetLoader()->LoadHits()) return;
  if(GetLoader()->LoadDigits()) return;

  
  for(Int_t iEvtN=iEvtNmin;iEvtN<iEvtNmax;iEvtN++){    
    Int_t nDigCh[7]={0,0,0,0,0,0,0};  
    Int_t iChHits[7]={0,0,0,0,0,0,0};
    Int_t nPrim[7]={0,0,0,0,0,0,0};
    Int_t nSec[7]={0,0,0,0,0,0,0};
    AliInfo(Form("events processed %i",iEvtN));
    if(GetLoader()->GetRunLoader()->GetEvent(iEvtN)) return;    
    AliStack *pStack = GetLoader()->GetRunLoader()->Stack();
    for(Int_t iPrimN=0;iPrimN<GetLoader()->TreeH()->GetEntries();iPrimN++){//prims loop
      GetLoader()->TreeH()->GetEntry(iPrimN);      
      for(Int_t iHitN=0;iHitN<Hits()->GetEntries();iHitN++){
        AliRICHHit *pHit = (AliRICHHit*)Hits()->At(iHitN);
        if(pHit->E()>0){
          iChHits[pHit->Ch()]++;
          if(pStack->Particle(pHit->GetTrack())->Rho()<0.01) nPrim[pHit->Ch()]++;else nSec[pHit->Ch()]++;
        }
      }
    }
    
    GetLoader()->TreeD()->GetEntry(0);
    for(Int_t iCh=0;iCh<7;iCh++){
      for(Int_t iDig=0;iDig<DigLst(iCh)->GetEntries();iDig++){
        AliRICHDigit *pDig=(AliRICHDigit*)DigLst(iCh)->At(iDig);
        nDigCh[pDig->Ch()]++;
      }
      Printf("Occupancy for chamber %i = %4.2f %% and charged prim tracks %i and sec. tracks %i with total %i",
        iCh,Float_t(nDigCh[iCh])*100/AliRICHDigit::kPadAll,nPrim[iCh],nSec[iCh],iChHits[iCh]);
    }      
    
    
  }//events loop
  GetLoader()->UnloadHits();
  GetLoader()->UnloadDigits();
  GetLoader()->GetRunLoader()->UnloadHeader();    
  GetLoader()->GetRunLoader()->UnloadKinematics();    
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliRICH::CluPrint(Int_t iEvtN)const
{
//prints a list of RICH clusters  for a given event
  Printf("List of RICH clusters for event %i",iEvtN);
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
void AliRICH::SummaryOfEvent(Int_t iEvtN) const
{
//prints a summary for a given event
  AliInfo(Form("Summary of event %i",iEvtN));
  GetLoader()->GetRunLoader()->GetEvent(iEvtN);    
  if(GetLoader()->GetRunLoader()->LoadHeader()) return;
  if(GetLoader()->GetRunLoader()->LoadKinematics()) return;
  AliStack *pStack=GetLoader()->GetRunLoader()->Stack();
  
  AliGenEventHeader* pGenHeader =  gAlice->GetHeader()->GenEventHeader();
  if(pGenHeader->InheritsFrom("AliGenHijingEventHeader")) {
    AliInfo(Form(" Hijing event with impact parameter b = %.2f (fm)",((AliGenHijingEventHeader*) pGenHeader)->ImpactParameter()));
  }
  Int_t nChargedPrimaries=0;
  for(Int_t i=0;i<pStack->GetNtrack();i++) {
    TParticle *pParticle = pStack->Particle(i);
    if(pParticle->IsPrimary()&&pParticle->GetPDG()->Charge()!=0) nChargedPrimaries++;
    }
  AliInfo(Form("Total number of         primaries %i",pStack->GetNprimary()));
  AliInfo(Form("Total number of charged primaries %i",nChargedPrimaries));
  AliInfo(Form("Total n. of tracks in stack(+sec) %i",pStack->GetNtrack()));
  GetLoader()->GetRunLoader()->UnloadHeader();
  GetLoader()->GetRunLoader()->UnloadKinematics();
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
