// **************************************
// Task used for the correction of determiantion of reconstructed jet spectra
// Compares input (gen) and output (rec) jets   
// *******************************************


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

 
#include <TROOT.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THnSparse.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TRefArray.h>
#include  "TDatabasePDG.h"

#include "AliAnalysisTaskJetBackgroundSubtract.h"
#include "AliAnalysisManager.h"
#include "AliAODHandler.h"
#include "AliAODTrack.h"
#include "AliAODJet.h"
#include "AliAODEvent.h"
#include "AliInputEventHandler.h"
#include "AliAODJetEventBackground.h"


ClassImp(AliAnalysisTaskJetBackgroundSubtract)

AliAnalysisTaskJetBackgroundSubtract::~AliAnalysisTaskJetBackgroundSubtract(){
  delete fJBArray;
  delete fOutJetArrayList;
  delete fInJetArrayList;
}

AliAnalysisTaskJetBackgroundSubtract::AliAnalysisTaskJetBackgroundSubtract(): 
  AliAnalysisTaskSE(),
  fAODOut(0x0),
  fAODIn(0x0),  
  fAODExtension(0x0),
  fJBArray(new TObjArray()),
  fBackgroundBranch(""),
  fNonStdFile(""),
  fReplaceString1("B0"),
  fReplaceString2("B%d"),
  fSubtraction(k4Area),
  fInJetArrayList(0x0),
  fOutJetArrayList(0x0),
  fh2CentvsRho(0x0),
  fh2CentvsSigma(0x0),
  fh2ShiftEta(0x0),
  fh2ShiftPhi(0x0),
  fh2ShiftEtaLeading(0x0),
  fh2ShiftPhiLeading(0x0),
  fHistList(0x0)  
{

}

AliAnalysisTaskJetBackgroundSubtract::AliAnalysisTaskJetBackgroundSubtract(const char* name):

  AliAnalysisTaskSE(name),
  fAODOut(0x0),
  fAODIn(0x0),  
  fAODExtension(0x0),
  fJBArray(new TObjArray()),
  fBackgroundBranch(""),
  fNonStdFile(""),
  fReplaceString1("B0"),
  fReplaceString2("B%d"),
  fSubtraction(k4Area),
  fInJetArrayList(0x0),
  fOutJetArrayList(0x0),
  fh2CentvsRho(0x0),
  fh2CentvsSigma(0x0),
  fh2ShiftEta(0x0),
  fh2ShiftPhi(0x0),
  fh2ShiftEtaLeading(0x0),
  fh2ShiftPhiLeading(0x0),
  fHistList(0x0)  
{
 DefineOutput(1, TList::Class());  
}



Bool_t AliAnalysisTaskJetBackgroundSubtract::Notify()
{
  //
  fAODIn = dynamic_cast<AliAODEvent*>(InputEvent());

  ResetOutJets();

  // Now we also have the Input Event Available! Fillvthe pointers in the list

  fInJetArrayList->Clear();
  fOutJetArrayList->Clear();

  for(int iJB = 0;iJB<fJBArray->GetEntries();iJB++){
    TObjString *ostr = (TObjString*)fJBArray->At(iJB);
 
  
    TClonesArray* jarray = 0;      
    if(!jarray&&fAODOut){
      jarray = (TClonesArray*)(fAODOut->FindListObject(ostr->GetString().Data()));
    }
    if(!jarray&&fAODExtension){
      jarray = (TClonesArray*)(fAODExtension->GetAOD()->FindListObject(ostr->GetString().Data()));
    }
    if(!jarray&&fAODIn){
      jarray = (TClonesArray*)(fAODIn->FindListObject(ostr->GetString().Data()));
    }

    if(!jarray){
      if(fDebug){
	Printf("%s:%d Input jet branch %s not found",(char*)__FILE__,__LINE__,ostr->GetString().Data());
      }
      continue;
    }
    
    TString newName(ostr->GetString().Data());
    newName.ReplaceAll(fReplaceString1.Data(),Form(fReplaceString2.Data(),fSubtraction));
    TClonesArray* jarrayOut = 0;      
    if(!jarrayOut&&fAODOut){
      jarrayOut = (TClonesArray*)(fAODOut->FindListObject(newName.Data()));
    }
    if(!jarrayOut&&fAODExtension){
      jarrayOut = (TClonesArray*)(fAODExtension->GetAOD()->FindListObject(newName.Data()));
    }

    if(!jarrayOut){
      if(fDebug){
	Printf("%s:%d Output jet branch %s not found",(char*)__FILE__,__LINE__,newName.Data());
	PrintAODContents();
      }
      continue;
    }
    if(jarrayOut&&jarray){
      fOutJetArrayList->Add(jarrayOut);
      fInJetArrayList->Add(jarray);
    }
  }
  return kTRUE;
}

void AliAnalysisTaskJetBackgroundSubtract::UserCreateOutputObjects()
{

  //
  // Create the output container
  //
  // Connect the AOD

  if (fDebug > 1) printf("AnalysisTaskJetBackgroundSubtract::UserCreateOutputObjects() \n");
  if(fNonStdFile.Length()!=0){
    
    // case that we have an AOD extension we need to fetch the jets from the extended output
    // we identifay the extension aod event by looking for the branchname
    AliAODHandler *aodH = dynamic_cast<AliAODHandler*>(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler());

    fAODExtension = (aodH?aodH->GetExtension(fNonStdFile.Data()):0);
    
    if(!fAODExtension){
      if(fDebug>1)Printf("AODExtension found for %s",fNonStdFile.Data());
    }
  }
  fAODOut = AODEvent();

  // usually we do not have the input already here

  if(!fInJetArrayList)fInJetArrayList =new TList();
  if(!fOutJetArrayList)fOutJetArrayList =new TList();

  for(int iJB = 0;iJB<fJBArray->GetEntries();iJB++){
    TObjString *ostr = (TObjString*)fJBArray->At(iJB);
    TString newName(ostr->GetString().Data());
    if(!newName.Contains(fReplaceString1.Data())){
      Printf("%s:%d cannot replace string %s in %s",(char*)__FILE__,__LINE__,fReplaceString1.Data(),
	     newName.Data());
      continue;
    }

    // add a new branch to the output for the background subtracted jets take the names from
    // the input jets and replace the background flag names
    TClonesArray *tca = new TClonesArray("AliAODJet", 0);
    newName.ReplaceAll(fReplaceString1.Data(),Form(fReplaceString2.Data(),fSubtraction));
    if(fDebug){
      Printf("%s:%d created branch \n %s from \n %s",(char*)__FILE__,__LINE__,newName.Data(),
	     ostr->GetString().Data());
    }
    tca->SetName(newName.Data());
    AddAODBranch("TClonesArray",&tca,fNonStdFile.Data());
  }
  

  if(!fHistList)fHistList = new TList();
  fHistList->SetOwner();

  for(int iJB = 0;iJB<fJBArray->GetEntries();iJB++){
    TObjString *ostr = (TObjString*)fJBArray->At(iJB);
    TString oldName(ostr->GetString().Data()); 
    TString newName(ostr->GetString().Data()); 
    if(!newName.Contains(fReplaceString1.Data())){
      Printf("%s:%d cannot replace string %s in %s",(char*)__FILE__,__LINE__,fReplaceString1.Data(),
	     newName.Data());
      continue;
    }
    newName.ReplaceAll(fReplaceString1.Data(),Form(fReplaceString2.Data(),fSubtraction));
    TH2F *hTmp = new TH2F(Form("h2PtInPtOut_%d",iJB),Form(";%s p_{T}; %s p_{T}",oldName.Data(),newName.Data()),200,0,200.,200,0.,200.);
    fHistList->Add(hTmp);
  }

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  //
  //  Histogram booking, add som control histograms here
  //    

 
    fh2CentvsRho = new TH2F("fh2CentvsRho","centrality vs background density",100,0.,100.,100,0.,150.);
    fh2CentvsSigma = new TH2F("fh2CentvsSigma","centrality vs backgroun sigma",100,0.,100.,100,0.,150.);
    fHistList->Add(fh2CentvsRho);
    fHistList->Add(fh2CentvsSigma);

   if(fSubtraction==k4Area){
    fh2ShiftEta = new TH2F("fh2ShiftEta","extended correction Eta",100,-0.9,0.9,100,-0.9,0.9);
    fh2ShiftPhi = new TH2F("fh2ShiftPhi","extended correction Phi",100,0.,6.5,100,0.,6.5);
    fh2ShiftEtaLeading = new TH2F("fh2ShiftEtaLeading","extended correction Eta",100,-0.9,0.9,100,-0.9,0.9);
    fh2ShiftPhiLeading = new TH2F("fh2ShiftPhiLeading","extended correction Phi",100,0.,6.5,100,0.,6.5);
     fHistList->Add(fh2ShiftEta);
     fHistList->Add(fh2ShiftPhi);
     fHistList->Add(fh2ShiftEtaLeading);
     fHistList->Add(fh2ShiftPhiLeading);
   }
    

  // =========== Switch on Sumw2 for all histos ===========
  for (Int_t i=0; i<fHistList->GetEntries(); ++i) {
    TH1 *h1 = dynamic_cast<TH1*>(fHistList->At(i));
    if (h1){
      h1->Sumw2();
      continue;
    }
    THnSparse *hn = dynamic_cast<THnSparse*>(fHistList->At(i));
    if(hn)hn->Sumw2();
  }
  TH1::AddDirectory(oldStatus);

  if(fBackgroundBranch.Length()==0)
    AliError(Form("%s:%d No BackgroundBranch defined",(char*)__FILE__,__LINE__));
  if(fJBArray->GetEntries()==0)
    AliError(Form("%s:%d No Jet Branches defined defined",(char*)__FILE__,__LINE__));
}

void AliAnalysisTaskJetBackgroundSubtract::Init()
{
  //
  // Initialization
  //
  if (fDebug > 1) printf("AnalysisTaskJetBackgroundSubtract::Init() \n");
}

void AliAnalysisTaskJetBackgroundSubtract::UserExec(Option_t */*option*/)
{

  if (fDebug > 1) printf("AnalysisTaskJetBackgroundSubtract::UserExec() \n");
  ResetOutJets();
  if(fBackgroundBranch.Length()==0||fJBArray->GetEntries()==0){
    if(fDebug)Printf("%s:%d No background subtraction done",(char*)__FILE__,__LINE__);
    PostData(1,fHistList);
  }
  if(fJBArray->GetEntries()!=fInJetArrayList->GetEntries()){
    if(fDebug)Printf("%s:%d different Array  sizes %d %d %d",(char*)__FILE__,__LINE__,fJBArray->GetEntries(),fInJetArrayList->GetEntries(),fOutJetArrayList->GetEntries());
    PostData(1,fHistList);
  }



  AliAODJetEventBackground*  evBkg = 0;
  TClonesArray*              bkgClusters = 0;
  TClonesArray*              bkgClustersRC = 0;
  TString bkgClusterName(fBackgroundBranch.Data());
  bkgClusterName.ReplaceAll(Form("%s_",AliAODJetEventBackground::StdBranchName()),"");
  TString bkgClusterRCName(Form("%s%s",bkgClusterName.Data(),"RandomCone")); 

  if(!evBkg&&!bkgClusters&&!bkgClustersRC&&fAODOut){
    evBkg = (AliAODJetEventBackground*)(fAODOut->FindListObject(fBackgroundBranch.Data()));
    bkgClusters = (TClonesArray*)(fAODOut->FindListObject(bkgClusterName.Data()));
    bkgClustersRC = (TClonesArray*)(fAODOut->FindListObject(bkgClusterRCName.Data()));

    if(fDebug&&bkgClusters)Printf("%s:%d Background cluster branch %s found",(char*)__FILE__,__LINE__,bkgClusterName.Data());
    if(fDebug&&bkgClustersRC)Printf("%s:%d Background cluster RC branch %s found",(char*)__FILE__,__LINE__,bkgClusterRCName.Data());
    if(fDebug&&evBkg)Printf("%s:%d Backgroundbranch %s found",(char*)__FILE__,__LINE__,fBackgroundBranch.Data());
  }
  if(!evBkg&&!bkgClusters&&!bkgClustersRC&&fAODExtension){
    evBkg = (AliAODJetEventBackground*)(fAODExtension->GetAOD()->FindListObject(fBackgroundBranch.Data()));
    bkgClusters = (TClonesArray*)(fAODExtension->GetAOD()->FindListObject(bkgClusterName.Data()));
    bkgClustersRC = (TClonesArray*)(fAODExtension->GetAOD()->FindListObject(bkgClusterRCName.Data()));
    if(fDebug&&bkgClusters)Printf("%s:%d Background cluster branch %s found",(char*)__FILE__,__LINE__,bkgClusterName.Data());
    if(fDebug&&bkgClustersRC)Printf("%s:%d Background cluster RC branch %s found",(char*)__FILE__,__LINE__,bkgClusterRCName.Data());

    if(fDebug&&evBkg)Printf("%s:%d Backgroundbranch %s found",(char*)__FILE__,__LINE__,fBackgroundBranch.Data());
  }

  if(!evBkg&&!bkgClusters&&!bkgClustersRC&&fAODIn){
    evBkg = (AliAODJetEventBackground*)(fAODIn->FindListObject(fBackgroundBranch.Data()));
    bkgClusters = (TClonesArray*)(fAODIn->FindListObject(bkgClusterName.Data()));
    bkgClustersRC = (TClonesArray*)(fAODIn->FindListObject(bkgClusterRCName.Data()));

    if(fDebug&&bkgClusters)Printf("%s:%d Background cluster branch %s found",(char*)__FILE__,__LINE__,bkgClusterName.Data());
    if(fDebug&&bkgClustersRC)Printf("%s:%d Background cluster RC branch %s found",(char*)__FILE__,__LINE__,bkgClusterRCName.Data());
    if(fDebug&&evBkg)Printf("%s:%d Backgroundbranch %s found",(char*)__FILE__,__LINE__,fBackgroundBranch.Data());
  }

  if(!evBkg&&(fSubtraction==kArea||fSubtraction==kRhoRecalc||fSubtraction==k4Area)){
    if(fDebug){
      Printf("%s:%d Backgroundbranch %s not found",(char*)__FILE__,__LINE__,fBackgroundBranch.Data());
      PrintAODContents();
    }
    PostData(1,fHistList);
    return;
  }

  if(!bkgClusters&&(fSubtraction==kRhoRecalc)){
    if(fDebug){
      Printf("%s:%d Background cluster branch %s not found",(char*)__FILE__,__LINE__,bkgClusterName.Data());
      PrintAODContents();
    }
    PostData(1,fHistList);
    return;
  }

  if(!bkgClustersRC&&(fSubtraction==kRhoRC)){
    if(fDebug){
      Printf("%s:%d Background cluster RC branch %s not found",(char*)__FILE__,__LINE__,bkgClusterRCName.Data());
      PrintAODContents();
    }
    PostData(1,fHistList);
    return;
  }
  // LOOP over all jet branches and subtract the background

   Float_t rho = 0;
   Float_t sigma=0.;
   Double_t meanarea = 0;
   TLorentzVector backgroundv;
   Float_t cent=0.;
   if(fAODOut)cent = fAODOut->GetHeader()->GetCentrality();
   if(evBkg)sigma=evBkg->GetSigma(1); 

   if(fSubtraction==kArea) rho = evBkg->GetBackground(1);
   if(fSubtraction==k4Area){
     rho = evBkg->GetBackground(0);
     sigma=evBkg->GetSigma(0);
   }

   if(fSubtraction==kRhoRecalc){
     meanarea=evBkg->GetMeanarea(1);
     rho =RecalcRho(bkgClusters,meanarea);
   }
   if(fSubtraction==kRhoRC) rho=RhoRC(bkgClustersRC);
   fh2CentvsRho->Fill(cent,rho);
   fh2CentvsSigma->Fill(cent,sigma);
   
   for(int iJB = 0;iJB<fInJetArrayList->GetEntries();iJB++){
    TClonesArray* jarray = (TClonesArray*)fInJetArrayList->At(iJB);
    TClonesArray* jarrayOut = (TClonesArray*)fOutJetArrayList->At(iJB);
    
    if(!jarray||!jarrayOut){
      Printf("%s:%d Array not found %d: %p %p",(char*)__FILE__,__LINE__,iJB,jarray,jarrayOut);
      continue;
    }
    TH2F* h2PtInOut = (TH2F*)fHistList->FindObject(Form("h2PtInPtOut_%d",iJB));
    // loop over all jets
    Int_t nOut = 0;
      

    for(int i = 0;i < jarray->GetEntriesFast();i++){
      AliAODJet *jet = (AliAODJet*)jarray->At(i);
      AliAODJet tmpNewJet(*jet);
      Bool_t bAdd = false;
     

      if(fSubtraction==kArea){	
	Double_t background = rho * jet->EffectiveAreaCharged();
	Float_t ptSub = jet->Pt() - background;	
	if(fDebug>2){
	  Printf("%s:%d Jet %d %3.3f %3.3f",(char*)__FILE__,__LINE__,i,jet->Pt(),ptSub);
	}
	if(ptSub<0){
	  // optionally rescale it and keep??
	  bAdd = RescaleJetMomentum(&tmpNewJet,0.1);
	  if(h2PtInOut)h2PtInOut->Fill(jet->Pt(),0.1);
	}
	else{
	  bAdd = RescaleJetMomentum(&tmpNewJet,ptSub);
	  if(h2PtInOut)h2PtInOut->Fill(jet->Pt(),ptSub);
	}
	// add background estimates to the new jet object
	// allows to recover old p_T and rho...
	tmpNewJet.SetBgEnergy(background,0);
      }// kAREA
      else if(fSubtraction==kRhoRecalc){
 	Double_t background = rho * jet->EffectiveAreaCharged();
	Float_t ptSub = jet->Pt() - background;	
        if(fDebug>2){
	  Printf("%s:%d Jet %d %3.3f %3.3f %3.3f %3.3f",(char*)__FILE__,__LINE__,i,jet->Pt(),ptSub,background,rho);}
	if(ptSub<0){
	  // optionally rescale it and keep??
	  bAdd = RescaleJetMomentum(&tmpNewJet,0.1);
	  if(h2PtInOut)h2PtInOut->Fill(jet->Pt(),0.1);
	}
	else{
	  bAdd = RescaleJetMomentum(&tmpNewJet,ptSub);
	  if(h2PtInOut)h2PtInOut->Fill(jet->Pt(),ptSub);
	}
	// add background estimates to the new jet object
	// allows to recover old p_T and rho...
	tmpNewJet.SetBgEnergy(background,0);

      }//kRhoRecalc
       else if(fSubtraction==kRhoRC){
	Double_t background = rho * jet->EffectiveAreaCharged();
	Float_t ptSub = jet->Pt() - background;	
	if(fDebug>2){	Printf("%s:%d Jet %d %3.3f %3.3f %3.3f %3.3f",(char*)__FILE__,__LINE__,i,jet->Pt(),ptSub,background,rho);}
	if(ptSub<0){
	  // optionally rescale it and keep??
	  bAdd = RescaleJetMomentum(&tmpNewJet,0.1);
	  if(h2PtInOut)h2PtInOut->Fill(jet->Pt(),0.1);
	}
	else{
	  bAdd = RescaleJetMomentum(&tmpNewJet,ptSub);
	  if(h2PtInOut)h2PtInOut->Fill(jet->Pt(),ptSub);
	}
	// add background estimates to the new jet object
	// allows to recover old p_T and rho...
	tmpNewJet.SetBgEnergy(background,0);

       }//kRhoRC

       else if(fSubtraction==k4Area&&jet->VectorAreaCharged()){

	 
	 backgroundv.SetPxPyPzE(rho*(jet->VectorAreaCharged())->Px(),rho*(jet->VectorAreaCharged())->Py(),rho*(jet->VectorAreaCharged())->Pz(),rho*(jet->VectorAreaCharged())->E());
	 if((backgroundv.E()>jet->E())&&(backgroundv.Pt()>jet->Pt())){
       
	   // optionally rescale it and keep??
	   bAdd = RescaleJetMomentum(&tmpNewJet,0.1);
	   if(h2PtInOut)h2PtInOut->Fill(jet->Pt(),0.1);
	 }
	 else{
	   bAdd = RescaleJet4vector(&tmpNewJet,backgroundv);
	 }
	 // add background estimates to the new jet object
	 // allows to recover old p_T and rho...
	 tmpNewJet.SetBgEnergy(backgroundv.P(),0);
	 
       }//kArea4vector  






      if(bAdd){
        AliAODJet *newJet = new ((*jarrayOut)[nOut++]) AliAODJet(tmpNewJet);
	// what about track references, clear for now...
	if(fSubtraction==k4Area){
         if(h2PtInOut)h2PtInOut->Fill(jet->Pt(),jet->Pt()-newJet->Pt());
         fh2ShiftEta->Fill(jet->Eta(),newJet->Eta());
         fh2ShiftPhi->Fill(jet->Phi(),newJet->Phi());
         if(i==0){fh2ShiftEtaLeading->Fill(jet->Eta(),newJet->Eta());
	   fh2ShiftPhiLeading->Fill(jet->Phi(),newJet->Phi());}}


	newJet->GetRefTracks()->Clear();
      }






    }



    // subtract the background
    

    // remove jets??

    // sort jets...

  }
  PostData(1, fHistList);
}

void AliAnalysisTaskJetBackgroundSubtract::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
  if (fDebug > 1) printf("AnalysisJetBackgroundSubtract: Terminate() \n");
}

Bool_t AliAnalysisTaskJetBackgroundSubtract::RescaleJetMomentum(AliAODJet *jet,Float_t pT){
  // keep the direction and the jet mass
  if(pT<=0)return kFALSE;
  Double_t pTold = jet->Pt();
  Double_t scale  = pT/pTold;
  Double_t mass  = jet->M();
  Double_t pNew = jet->P() * scale;
  jet->SetPxPyPzE(scale*jet->Px(),scale*jet->Py(),scale*jet->Pz(),TMath::Sqrt(mass*mass+pNew*pNew));
 


  return kTRUE;
}

Bool_t AliAnalysisTaskJetBackgroundSubtract::RescaleJet4vector(AliAODJet *jet,TLorentzVector backgroundv){
  
  if(backgroundv.Pt()<0.) return kFALSE;
  jet->SetPxPyPzE(jet->Px()-backgroundv.Px(),jet->Py()-backgroundv.Py(),jet->Pz()-backgroundv.Pz(),jet->E()-backgroundv.E());
 
 return kTRUE;
}








Double_t AliAnalysisTaskJetBackgroundSubtract::RecalcRho(TClonesArray* bkgClusters,Double_t meanarea){
  
  Double_t ptarea=0.;
  Int_t count=0;
  Double_t rho=0.; 
  const Double_t Rlimit2=0.8*0.8;  //2*jet radius.
  TClonesArray* jarray=0;
  
  for(int iJB = 0;iJB<fInJetArrayList->GetEntries();iJB++){
    TObjString *ostr = (TObjString*)fInJetArrayList->At(iJB);
    TString jetref=ostr->GetString().Data();
    if(jetref.Contains("ANTIKT04")){ 
      jarray = (TClonesArray*)fInJetArrayList->At(iJB);
    }
  }
  if(!jarray)return rho;
  if(jarray->GetEntries()>=2){ 
    AliAODJet *first = (AliAODJet*)(jarray->At(0)); 
    AliAODJet *second= (AliAODJet*)(jarray->At(1)); 
    for(Int_t k=0;k<bkgClusters->GetEntriesFast();k++){
      AliAODJet *clus = (AliAODJet*)(bkgClusters->At(k));
      if(TMath::Abs(clus->Eta())>0.5) continue;
      if((clus->EffectiveAreaCharged())<0.1*meanarea) continue; 
      Double_t distance1=(first->Eta()-clus->Eta())*(first->Eta()-clus->Eta())+
	(first->Phi()-clus->Phi())*(first->Phi()-clus->Phi());
      Double_t distance2= (second->Eta()-clus->Eta())*(second->Eta()-clus->Eta())+
	(second->Phi()-clus->Phi())*(second->Phi()-clus->Phi());
      if((distance1<Rlimit2)||(distance2<Rlimit2)) continue;    
      ptarea=ptarea+clus->Pt()/clus->EffectiveAreaCharged(); 
      count=count+1;}
    if(count!=0) rho=ptarea/count; 
  }        
  return rho;
}

   Double_t AliAnalysisTaskJetBackgroundSubtract::RhoRC(TClonesArray* bkgClustersRC){
  
       Double_t ptarea=0.;
       Int_t count=0;
       Double_t rho=0.; 
       const Double_t Rlimit2=0.8*0.8;  //2*jet radius.
       TClonesArray* jarray=0;
        for(int iJB = 0;iJB<fInJetArrayList->GetEntries();iJB++){
	  TObjString *ostr = (TObjString*)fInJetArrayList->At(iJB);
	  TString jetref=ostr->GetString().Data();
	  if(jetref.Contains("ANTIKT04")){ 
	    jarray = (TClonesArray*)fInJetArrayList->At(iJB);
	  }
	}
	if(!jarray)return rho;

         if(jarray->GetEntries()>=2){ 
	   AliAODJet *first = (AliAODJet*)(jarray->At(0)); 
	   AliAODJet *second=(AliAODJet*)(jarray->At(1)); 
         for(Int_t k=0;k<bkgClustersRC->GetEntriesFast();k++){
	   AliAODJet *clus = (AliAODJet*)(bkgClustersRC->At(k));
	   if(TMath::Abs(clus->Eta())>0.5) continue;
	   Double_t distance1=(first->Eta()-clus->Eta())*(first->Eta()-clus->Eta())+
	     (first->Phi()-clus->Phi())*(first->Phi()-clus->Phi());
	   Double_t distance2= (second->Eta()-clus->Eta())*(second->Eta()-clus->Eta())+
	     (second->Phi()-clus->Phi())*(second->Phi()-clus->Phi());
	   if((distance1<Rlimit2)||(distance2<Rlimit2)) continue;    
	   ptarea=ptarea+clus->Pt()/clus->EffectiveAreaCharged(); 
	   count=count+1;}
         if(count!=0) rho=ptarea/count;  }
         return rho;
}









void AliAnalysisTaskJetBackgroundSubtract::ResetOutJets(){
  if(!fOutJetArrayList)return;
  for(int iJB = 0;iJB<fOutJetArrayList->GetEntries();iJB++){
    TClonesArray* jarray = (TClonesArray*)fOutJetArrayList->At(iJB);
    if(jarray)jarray->Delete();
  }
}


void AliAnalysisTaskJetBackgroundSubtract::PrintAODContents(){
  if(fAODIn){
    Printf("%s:%d >>>>>> Input",(char*)__FILE__,__LINE__);
    fAODIn->Print();
  }
  if(fAODExtension){
    Printf("%s:%d >>>>>> Extenstion",(char*)__FILE__,__LINE__);
    fAODExtension->GetAOD()->Print();
  }
  if(fAODOut){
    Printf("%s:%d >>>>>> Output",(char*)__FILE__,__LINE__);
    fAODOut->Print();
  }
}
