#include <TChain.h>
#include <TList.h>

#include <TTree.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TCanvas.h>

#include "AliVEvent.h"
#include "AliVParticle.h"

#include "AliESDEvent.h"
#include "AliESDtrack.h"

#include "AliAODEvent.h"
#include "AliAODTrack.h"

#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliAnalysisManager.h"
#include "AliMCEvent.h"
#include "TParticle.h"
#include "AliStack.h"

#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODHandler.h"
#include "AliAODTrack.h"
#include "AliAODJet.h"
#include "AliAODMCParticle.h"

#include "AliAnalysisTaskQGSep.h"
#include "AliAnalysisHelperJetTasks.h"
ClassImp(AliAnalysisTaskQGSep)

//________________________________________________________________________
AliAnalysisTaskQGSep::AliAnalysisTaskQGSep(const char *name)
  : AliAnalysisTaskSE(name),
    fBranchRec("jets"),
    fUseMC(kFALSE),
    fUseAOD(kFALSE),
    fXsection(0),
    fWeight(0),
    fMyAODEvent(0),
    fOutputList(0),
    fpHistPtAvEQ(0),
    fpHistPtAvEG(0),
    fpHistDrEQ(0),
    fpHistDrEG(0),
    fpHistDrE(0),
    fpHistPtAvE(0),
    fpHistDrE3(0),
    fpHistPtAvE3(0)

{
  // Constructor

  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
void AliAnalysisTaskQGSep::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

  fOutputList = new TList();

  if(fUseMC){

    //histos for quarks
    fpHistPtAvEQ = new TProfile("fpHistPtAvEQ", "", 100, 0, 100, 0, 10);
    fOutputList->Add(fpHistPtAvEQ);
    fpHistDrEQ = new TProfile("fpHistDrEQ", "", 100, 0, 100, 0, 0.5);
    fOutputList->Add(fpHistDrEQ);

    //histos for gluons
    fpHistPtAvEG = new TProfile("fpHistPtAvEG", "", 100, 0, 100, 0., 10.);
    fOutputList->Add(fpHistPtAvEG);
    fpHistDrEG = new TProfile("fpHistDrEG", "", 100, 0, 100, 0, 0.5);
    fOutputList->Add(fpHistDrEG);
  }
  
  //histos for full spetra, no separation
  fpHistDrE = new TProfile("fpHistDrE", "", 100, 0, 100, 0, 0.5);
  fOutputList->Add(fpHistDrE);
  fpHistPtAvE = new TProfile("fpHistPtAvE", "", 100, 0, 100, 0., 10.);
  fOutputList->Add(fpHistPtAvE);
  fpHistDrE3 = new TProfile("fpHistDrE3", "", 100, 0, 100, 0, 0.5);
  fOutputList->Add(fpHistDrE3);
  fpHistPtAvE3 = new TProfile("fpHistPtAvE3", "", 100, 0, 100, 0., 10.);
  fOutputList->Add(fpHistPtAvE3);

  for(Int_t i = 0; i < fOutputList->GetEntries(); i++){
    TH1 * h = dynamic_cast<TH1*>(fOutputList->At(i));
    if(h) h->Sumw2();
  }

  if(fDebug)Printf("~# QGSep User objects created");
}
  
//__________________________________________________________________
Bool_t AliAnalysisTaskQGSep::Notify()
{
  
  if(fUseMC){
    //
    // Implemented Notify() to read the cross sections
    // and number of trials from pyxsec.root
    // 
    
    
    TTree *tree = AliAnalysisManager::GetAnalysisManager()->GetTree();
    Float_t xsection = 0;
    Float_t ftrials  = 1;
    
    Float_t fAvgTrials = 1;
    if(tree){
      TFile *curfile = tree->GetCurrentFile();
      if (!curfile) {
	Error("Notify","No current file");
	return kFALSE;
      }
      AliAnalysisHelperJetTasks::PythiaInfoFromFile(curfile->GetName(),xsection,ftrials);
      // construct a poor man average trials 
      Float_t nEntries = (Float_t)tree->GetTree()->GetEntries();
      if(ftrials>=nEntries)fAvgTrials = ftrials/nEntries; // CKB take this into account for normalisation
    }  
  
    if(xsection>0){
      fXsection  = xsection;
      fWeight = fXsection/fAvgTrials;
    }
    else fWeight = 1;
    return kTRUE;
  }
  
  else{
    fWeight = 1;
    return kTRUE;
  }
  
}

//________________________________________________________________________
void AliAnalysisTaskQGSep::UserExec(Option_t *)
{
  // Main loop
  // Called for each event
  
  if(fUseAOD){
    AliVEvent *event = InputEvent();
    if (!event) {
      Error("UserExec", "Could not retrieve event");
      return;
    }
  
  fMyAODEvent = dynamic_cast<AliAODEvent*> (InputEvent());
  }

  else{
        //  assume that the AOD is in the general output...
    fMyAODEvent  = AODEvent();
    if(!fMyAODEvent){
      Printf("%s:%d AODEvent not found in the Output",(char*)__FILE__,__LINE__);
      return;
    }    
  }

  if (fMyAODEvent){
    if(fUseMC)
      LoopAODMC();
    else
      LoopAOD();
  }
	
  // Post output data.
  PostData(1, fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskQGSep::Terminate(Option_t*)
{
  // Draw result to the screen
  // Called once at the end of the query
  
  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputList) {
    Error("Terminate","fOutputList not available");
    return;
  }
  
}

//________________________________________________________________________
void AliAnalysisTaskQGSep::LoopAOD(){
  AliAODJet recJets[4];
  AliAODJet jets[4];
  AliAODJet rJets[4];
  Int_t nRecJets = 0;
  
  //array of reconstructed jets from the AOD input
  TClonesArray *aodRecJets = dynamic_cast<TClonesArray*>(fMyAODEvent->FindListObject(fBranchRec.Data()));
  if(!aodRecJets){
    PostData(1, fOutputList);
    return;
  }
   
  // reconstructed jets
  nRecJets = aodRecJets->GetEntries(); 
  if(fDebug)Printf("--- Jets found in bRec: %d", nRecJets);
  nRecJets = TMath::Min(nRecJets, 4);
  for(int ir = 0;ir < nRecJets;++ir)
    {
      AliAODJet *tmp = dynamic_cast<AliAODJet*>(aodRecJets->At(ir));
      if(!tmp)continue;
      jets[ir] = *tmp;
    }
  
  Int_t counter = 0;
  Int_t tag = 0;
  Int_t nGenSel = 0;
  
  TLorentzVector v[4];
  Double_t eSum = 0.;
  Double_t pxSum = 0.;
  Double_t pySum = 0.;
  Double_t pzSum = 0.;
  
  for(Int_t i = 0; i < nRecJets; i++)
    {
      if(nRecJets == 1)
	{
	  rJets[nGenSel] = jets[i];
	  v[nGenSel].SetPxPyPzE(jets[i].Px(), jets[i].Py(), jets[i].Pz(), jets[i].E());
	  eSum += jets[i].E();
	  pxSum += jets[i].Px();
	  pySum += jets[i].Py();
	  pzSum += jets[i].Pz();		    
	  nGenSel++;
	}
      else
	{
	  counter = 0;
	  tag = 0;
	  for(Int_t j = 0; j < nRecJets; j++)
	    {
	      if(i!=j)
		{
		  Double_t dRij = jets[i].DeltaR(&jets[j]);
		  counter++;
		  if(dRij > 2*0.4) tag++;
		}
	    }
	  if(counter!=0)
	    {
	      if(tag/counter == 1)
		{
		  rJets[nGenSel] = jets[i];
		  v[nGenSel].SetPxPyPzE(jets[i].Px(), jets[i].Py(), jets[i].Pz(), jets[i].E());
		  eSum += jets[i].E();
		  pxSum += jets[i].Px();
		  pySum += jets[i].Py();
		  pzSum += jets[i].Pz();		    
		  nGenSel++;
		}
	    }
	}
    }
   
  nRecJets = nGenSel;
  
  if(nRecJets == 0){
    PostData(1, fOutputList);
    return;
  }

  TLorentzVector vB;
  vB.SetPxPyPzE(pxSum, pySum, pzSum, eSum);
  Double_t e[4];
  for(Int_t i = 0; i < nRecJets; i++){
    v[i].Boost(-vB.Px()/vB.E(),-vB.Py()/vB.E(),-vB.Pz()/vB.E());
    e[i] = v[i].E();
  }

  Int_t idxj[4];
  TMath::Sort(nRecJets, e, idxj);
  for(Int_t i = 0; i < nRecJets; i++){
    recJets[i] = rJets[idxj[i]];
  }
  
   
  for(Int_t iJ = 0; iJ < nRecJets; iJ++){
    Double_t pTsum = 0.;
    TRefArray * tra = dynamic_cast<TRefArray*>(recJets[iJ].GetRefTracks());
    if(!tra) continue;
    Int_t nAODtracks = TMath::Min(1000, tra->GetEntries());
    Double_t dR[1000];
    for(Int_t iT = 0; iT < nAODtracks; iT++){
      AliAODTrack * jetTrack = dynamic_cast<AliAODTrack*>(tra->At(iT));
      if(!jetTrack) continue;
      pTsum += jetTrack->Pt();
      dR[iT] = recJets[iJ].DeltaR(jetTrack);
    }
     
 
    fpHistPtAvE->Fill(recJets[iJ].E(), (Double_t)pTsum/(Double_t)nAODtracks, fWeight);
           
    if(iJ > 1) //fill mulit-jet histo
      fpHistPtAvE3->Fill(recJets[iJ].E(), (Double_t)pTsum/(Double_t)nAODtracks, fWeight);
      
    Int_t idxAOD[1000];
    TMath::Sort(nAODtracks, dR, idxAOD, kFALSE);
     
    Double_t pTsum90Inv=0.;
    for(Int_t iT = 0; iT < nAODtracks; iT++){
      AliAODTrack * track = dynamic_cast<AliAODTrack*>(tra->At(idxAOD[iT]));
      if(!track) continue;
      pTsum90Inv += track->Pt();
      if(pTsum90Inv >= 0.9*pTsum){
	Double_t deltaR = recJets[iJ].DeltaR(track);

	fpHistDrE->Fill(recJets[iJ].E(), deltaR, fWeight);

	if(iJ > 1)
	  fpHistDrE3->Fill(recJets[iJ].E(), deltaR, fWeight);

	break;
      }
    }    
  }
}
 
//__________________________________________________________________
void AliAnalysisTaskQGSep::LoopAODMC(){
  
  AliAODJet recJets[4];
  AliAODJet jets[4];
  AliAODJet rJets[4];
  Int_t nRecJets = 0;
  
  //array of reconstructed jets from the AOD input
  TClonesArray *aodRecJets = dynamic_cast<TClonesArray*>(fMyAODEvent->FindListObject(fBranchRec.Data()));
  if(!aodRecJets){
    PostData(1, fOutputList);
    return;
  }
  
  // reconstructed jets
  nRecJets = aodRecJets->GetEntries(); 
  if(fDebug)Printf("--- Jets found in bRec: %d", nRecJets);
  nRecJets = TMath::Min(nRecJets, 4);
  for(int ir = 0;ir < nRecJets;++ir)
    {
      AliAODJet *tmp = dynamic_cast<AliAODJet*>(aodRecJets->At(ir));
      if(!tmp)continue;
      jets[ir] = *tmp;
    }
  
  Int_t counter = 0;
  Int_t tag = 0;
  Int_t nGenSel = 0;
  
  TLorentzVector v[4];
  Double_t eSum = 0.;
  Double_t pxSum = 0.;
  Double_t pySum = 0.;
  Double_t pzSum = 0.;
  
  for(Int_t i = 0; i < nRecJets; i++)
    {
      if(nRecJets == 1)
	{
	  rJets[nGenSel] = jets[i];
	  v[nGenSel].SetPxPyPzE(jets[i].Px(), jets[i].Py(), jets[i].Pz(), jets[i].E());
	  eSum += jets[i].E();
	  pxSum += jets[i].Px();
	  pySum += jets[i].Py();
	  pzSum += jets[i].Pz();		    
	  nGenSel++;
	}
      else
	{
	  counter = 0;
	  tag = 0;
	  for(Int_t j = 0; j < nRecJets; j++)
	    {
	      if(i!=j)
		{
		  Double_t dRij = jets[i].DeltaR(&jets[j]);
		  counter++;
		  if(dRij > 2*0.4) tag++;
		}
	    }
	  if(counter!=0)
	    {
	      if(tag/counter == 1)
		{
		  rJets[nGenSel] = jets[i];
		  v[nGenSel].SetPxPyPzE(jets[i].Px(), jets[i].Py(), jets[i].Pz(), jets[i].E());
		  eSum += jets[i].E();
		  pxSum += jets[i].Px();
		  pySum += jets[i].Py();
		  pzSum += jets[i].Pz();		    
		  nGenSel++;
		}
	    }
	}
    }
  
  nRecJets = nGenSel;
  
  if(nRecJets == 0){
    PostData(1, fOutputList);
    return;
   }
   
  TLorentzVector vB;
  vB.SetPxPyPzE(pxSum, pySum, pzSum, eSum);
  Double_t e[4];
  for(Int_t i = 0; i < nRecJets; i++){
    v[i].Boost(-vB.Px()/vB.E(),-vB.Py()/vB.E(),-vB.Pz()/vB.E());
    e[i] = v[i].E();
  }

  Int_t idxj[4];
  TMath::Sort(nRecJets, e, idxj);
  for(Int_t i = 0; i < nRecJets; i++){
    recJets[i] = rJets[idxj[i]];
  }
  
  Int_t nMCtracks = 0;
  TClonesArray *tca = dynamic_cast<TClonesArray*>(fMyAODEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if(!tca) return;
  nMCtracks = TMath::Min(tca->GetEntries(), 10000);
  
  //sort AliAODMCParticles in pT
  Double_t pTMC[10000];
  for(Int_t iMC = 0; iMC < nMCtracks; iMC++){
    AliAODMCParticle *partMC = dynamic_cast<AliAODMCParticle*>(tca->At(iMC));
    if(!partMC) continue;
    pTMC[iMC]=partMC->Pt();
  }
  
  Int_t idxMC[10000];
  TMath::Sort(nMCtracks, pTMC, idxMC);
  
  
  Int_t flagQ[4], flagG[4];   
  for(Int_t iJ = 0; iJ < nRecJets; iJ++){
    //flag jet as q/g
    
    //look for highest momentum parton in the jet cone
    for(Int_t iMC = 0; iMC < nMCtracks; iMC++){
      AliAODMCParticle *partMC = dynamic_cast<AliAODMCParticle*>(tca->At(idxMC[iMC]));
      if(!partMC) continue;
      Double_t r = recJets[iJ].DeltaR(partMC);
      if(r < 0.4){
	if(TMath::Abs(partMC->GetPdgCode()) < 9)
	  flagQ[iJ] = 1;
	if(TMath::Abs(partMC->GetPdgCode()) == 21)
	  flagG[iJ] = 1;
	break;
      }
    }

    Double_t pTsum = 0.;
    TRefArray * tra = dynamic_cast<TRefArray*>(recJets[iJ].GetRefTracks());
    if(!tra) continue;
    Int_t nAODtracks = TMath::Min(1000, tra->GetEntries());
    Double_t dR[1000];
    for(Int_t iT = 0; iT < nAODtracks; iT++){
      AliAODTrack * jetTrack = dynamic_cast<AliAODTrack*>(tra->At(iT));
      if(!jetTrack) continue;
      pTsum += jetTrack->Pt();
      dR[iT] = recJets[iJ].DeltaR(jetTrack);
    }

    fpHistPtAvE->Fill(recJets[iJ].E(), (Double_t)pTsum/(Double_t)nAODtracks, fWeight);
    if(flagQ[iJ] == 1) fpHistPtAvEQ->Fill(recJets[iJ].E(), (Double_t)pTsum/(Double_t)nAODtracks, fWeight);
    if(flagG[iJ] == 1) fpHistPtAvEG->Fill(recJets[iJ].E(), (Double_t)pTsum/(Double_t)nAODtracks, fWeight);
        
    if(iJ > 1) //fill mulit-jet histo
      fpHistPtAvE3->Fill(recJets[iJ].E(), (Double_t)pTsum/(Double_t)nAODtracks, fWeight);
      
    Int_t idxAOD[1000];
    TMath::Sort(nAODtracks, dR, idxAOD, kFALSE);

    Double_t pTsum90Inv=0.;
    for(Int_t iT = 0; iT < nAODtracks; iT++){
      AliAODTrack * track = dynamic_cast<AliAODTrack*>(tra->At(idxAOD[iT]));
      if(!track) continue;
      pTsum90Inv += track->Pt();
      if(pTsum90Inv >= 0.9*pTsum){
	Double_t deltaR = recJets[iJ].DeltaR(track);

	fpHistDrE->Fill(recJets[iJ].E(), deltaR, fWeight);
	if(flagQ[iJ] == 1) fpHistDrEQ->Fill(recJets[iJ].E(), deltaR, fWeight);
	if(flagG[iJ] == 1) fpHistDrEG->Fill(recJets[iJ].E(), deltaR, fWeight);

	if(iJ > 1)
	  fpHistDrE3->Fill(recJets[iJ].E(), deltaR, fWeight);

	break;
      }
    }    
  }
}
