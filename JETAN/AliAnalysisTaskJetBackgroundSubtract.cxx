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
#include <THnSparse.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <TList.h>
#include <TLorentzVector.h>
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
  fSubtraction(kArea),
  fInJetArrayList(0x0),
  fOutJetArrayList(0x0),
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
  fSubtraction(kArea),
  fInJetArrayList(0x0),
  fOutJetArrayList(0x0),
  fHistList(0x0)  
{
 DefineOutput(1, TList::Class());  
}



Bool_t AliAnalysisTaskJetBackgroundSubtract::Notify()
{
  //
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
    // 
    // case that we have an AOD extension we need to fetch the jets from the extended output
    // we identifay the extension aod event by looking for the branchname
    AliAODHandler *aodH = dynamic_cast<AliAODHandler*>(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler());
    fAODExtension = aodH->GetExtension(fNonStdFile.Data());
    
    if(!fAODExtension){
      if(fDebug>1)Printf("AODExtension found for %s",fNonStdFile.Data());
    }
  }
  fAODOut = AODEvent();
  fAODIn = dynamic_cast<AliAODEvent*>(InputEvent());


  // collect the jet branches 

  if(!fInJetArrayList)fInJetArrayList =new TList();
  if(!fOutJetArrayList)fOutJetArrayList =new TList();

  for(int iJB = 0;iJB<fJBArray->GetEntries();iJB++){
    TClonesArray* jarray = 0;      
    TObjString *ostr = (TObjString*)fJBArray->At(iJB);
    if(!jarray&&fAODOut){
      jarray = (TClonesArray*)(fAODOut->FindListObject(ostr->GetString().Data()));
    }
    if(!jarray&&fAODExtension){
      jarray = (TClonesArray*)(fAODExtension->GetAOD()->FindListObject(ostr->GetString().Data()));
    }
    if(!jarray){
      if(fDebug)Printf("%s:%d jet branch %s not found",(char*)__FILE__,__LINE__,ostr->GetString().Data());
      continue;
    }

    TString newName(jarray->GetName());
    if(!newName.Contains(fReplaceString1.Data())){
      Printf("%s:%d cannot replace string %s in %s",(char*)__FILE__,__LINE__,fReplaceString1.Data(),
	     jarray->GetName());
      continue;
    }

    fInJetArrayList->Add(jarray);

    // add a new branch to the output for the background subtracted jets take the names from
    // the input jets and replace the background flag names
    TClonesArray *tca = new TClonesArray("AliAODJet", 0);

    newName.ReplaceAll(fReplaceString1.Data(),Form(fReplaceString2.Data(),fSubtraction));
    if(fDebug){
      Printf("%s:%d created branch \n %s from \n %s",(char*)__FILE__,__LINE__,newName.Data(),
	     jarray->GetName());
    }
    tca->SetName(newName.Data());
    AddAODBranch("TClonesArray",&tca,fNonStdFile.Data());
    fOutJetArrayList->Add(tca);
  }
  



  if(!fHistList)fHistList = new TList();
  fHistList->SetOwner();

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  //
  //  Histogram booking, add som control histograms here
  //    

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

  ResetOutJets();
  if(fBackgroundBranch.Length()==0||fJBArray->GetEntries()==0){
    PostData(1,fHistList);
  }
    if (fDebug > 1) printf("AnalysisTaskJetBackgroundSubtract::UserExec() \n");


  static AliAODJetEventBackground*  evBkg = 0;
  if(!evBkg&&fAODOut){
    evBkg = (AliAODJetEventBackground*)(fAODOut->FindListObject(fBackgroundBranch.Data()));
  }
  if(!evBkg&&fAODExtension){
    evBkg = (AliAODJetEventBackground*)(fAODExtension->GetAOD()->FindListObject(fBackgroundBranch.Data()));
  }
  if(!evBkg&&fAODIn){
    evBkg = (AliAODJetEventBackground*)(fAODIn->FindListObject(fBackgroundBranch.Data()));
  }

  if(!evBkg){
    if(fDebug)Printf("%s:%d Backroundbranch %s not found",(char*)__FILE__,__LINE__,fBackgroundBranch.Data());
    PostData(1,fHistList);
  }


  // LOOP over all jet branches and subtract the background

  Float_t rho = 0;
  if(fSubtraction==kArea)rho =  evBkg->GetBackground(2);
  
  for(int iJB = 0;iJB<fInJetArrayList->GetEntries();iJB++){
    TClonesArray* jarray = (TClonesArray*)fInJetArrayList->At(iJB);
    TClonesArray* jarrayOut = (TClonesArray*)fOutJetArrayList->At(iJB);
    if(!jarray||!jarrayOut)continue;
    
    // loop over all jets
    Int_t nOut = 0;
    for(int i = 0;i < jarray->GetEntriesFast();i++){
      AliAODJet *jet = (AliAODJet*)jarray->At(i);
      AliAODJet tmpNewJet(*jet);
      Bool_t bAdd = false;
      if(fSubtraction==kArea){	
	Float_t ptSub = jet->Pt() -  rho * jet->EffectiveAreaCharged();
	if(fDebug>2){
	  Printf("%s:%d Jet %d %3.3f %3.3f",(char*)__FILE__,__LINE__,i,jet->Pt(),ptSub);
	}
	if(ptSub<0){
	  // optionally rescale it and keep??
	  bAdd = RescaleJetMomentum(&tmpNewJet,0.1);
	}
	else{
	  bAdd = RescaleJetMomentum(&tmpNewJet,ptSub);
	}
      }// kAREA

      if(bAdd){
        new ((*jarrayOut)[nOut++]) AliAODJet(tmpNewJet);
	// optionally add background estimates to the jet object? CKB
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

void AliAnalysisTaskJetBackgroundSubtract::ResetOutJets(){
  if(!fOutJetArrayList)return;
  for(int iJB = 0;iJB<fOutJetArrayList->GetEntries();iJB++){
    TClonesArray* jarray = (TClonesArray*)fOutJetArrayList->At(iJB);
    if(jarray)jarray->Delete();
  }
}
