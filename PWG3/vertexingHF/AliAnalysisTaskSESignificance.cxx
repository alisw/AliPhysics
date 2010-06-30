/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

/////////////////////////////////////////////////////////////
//
// AliAnalysisTaskSESignificane to calculate effects on 
// significance of D mesons  cut 
// Authors: G. Ortona, ortona@to.infn.it
// F. Prino, prino@to.infn.it
// Renu Bala, bala@to.infn.it
// Chiara Bianchin, cbianchi@pd.infn.it
/////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TClonesArray.h>
#include <TCanvas.h>
#include <TList.h>
#include <TFile.h>
#include <TString.h>
#include <TH1F.h>
#include <TDatabasePDG.h>

#include <AliLog.h>
#include "AliAnalysisManager.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoDecayHF4Prong.h"
#include "AliAODRecoCascadeHF.h"

#include "AliAnalysisTaskSE.h"
#include "AliRDHFCutsDplustoKpipi.h"
#include "AliRDHFCutsD0toKpipipi.h"
#include "AliRDHFCutsDstoKKpi.h"
#include "AliRDHFCutsDStartoKpipi.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliRDHFCutsLctopKpi.h"
#include "AliMultiDimVector.h"

#include "AliAnalysisTaskSESignificance.h"

ClassImp(AliAnalysisTaskSESignificance)


//________________________________________________________________________
AliAnalysisTaskSESignificance::AliAnalysisTaskSESignificance():
  AliAnalysisTaskSE(),
  fOutput(0),
  fCutList(0),
  fHistNEvents(0),
  fUpmasslimit(1.965),
  fLowmasslimit(1.765),
  fRDCuts(0),
  fNPtBins(0),
  fReadMC(kFALSE),
  fDecChannel(0),
  fSelectionlevel(0)
{
  // Default constructor
}

//________________________________________________________________________
AliAnalysisTaskSESignificance::AliAnalysisTaskSESignificance(const char *name, TList* listMDV,AliRDHFCuts *rdCuts,Int_t decaychannel,Int_t selectionlevel):
  AliAnalysisTaskSE(name),
  fOutput(0),
  fCutList(listMDV),
  fHistNEvents(0),
  fUpmasslimit(0),
  fLowmasslimit(0),
  fRDCuts(rdCuts),
  fNPtBins(0),
  fReadMC(kFALSE),
  fDecChannel(decaychannel),
  fSelectionlevel(selectionlevel)
{

  Int_t pdg=421;
  switch(fDecChannel){
  case 0:
    pdg=411;
    break;
  case 1:
    pdg=421;
    break;
  case 2:
    pdg=413;
    break;
  case 3:
    pdg=431;
    break;
  case 4:
    pdg=421;
    break;
  case 5:
    pdg=4122;
    break;
  }

  SetMassLimits(0.1,pdg); //check range
  fNPtBins=fRDCuts->GetNPtBins();

  if(fDebug>1)fRDCuts->PrintAll();
   // Output slot #1 writes into a TList container
  DefineOutput(1,TList::Class());  //My private output
  DefineOutput(2,TList::Class());
  CheckConsistency();
}

 //________________________________________________________________________
AliAnalysisTaskSESignificance::~AliAnalysisTaskSESignificance()
{
  // Destructor
  if (fOutput) {
    delete fOutput;
    fOutput = 0;
  }
  if (fCutList) {
    delete fCutList;
    fCutList = 0;
  }
  if(fHistNEvents){
    delete fHistNEvents;
    fHistNEvents=0;
  }

/*
  if(fRDCuts) {
    delete fRDCuts;
    fRDCuts= 0;
  } 
*/
  
}
//_________________________________________________________________
Bool_t AliAnalysisTaskSESignificance::CheckConsistency(){

  Bool_t result = kTRUE;

  const Int_t nvars=fRDCuts->GetNVars();//ForOpt();
  //Float_t *vars = new Float_t[nvars];
  Bool_t *varsforopt = fRDCuts->GetVarsForOpt();
  Bool_t *uppervars = fRDCuts->GetIsUpperCut();
  TString *names = fRDCuts->GetVarNames();

  for(Int_t i=0;i<fNPtBins;i++){
    TString mdvname=Form("multiDimVectorPtBin%d",i);
    Int_t ic=0;
 
    for(Int_t ivar=0;ivar<nvars;ivar++){
      if(varsforopt[ivar]){
	Float_t min = ((AliMultiDimVector*)fCutList->FindObject(mdvname.Data()))->GetMinLimit(ic);
	Float_t max = ((AliMultiDimVector*)fCutList->FindObject(mdvname.Data()))->GetMaxLimit(ic);
	if(min==max){
	  printf("AliAnalysisTaskSESignificance::CheckConsistency: ERROR! \n tight and loose cut for optimization variable number %d are the same in ptbin %d\n",ic,i);
	  result = kFALSE;
	}
	Bool_t lowermdv = ((AliMultiDimVector*)fCutList->FindObject(mdvname.Data()))->GetGreaterThan(ic);
	if(uppervars[ivar]&&lowermdv){
	  printf("AliAnalysisTaskSESignificance::CheckConsistency: WARNING! %s is declared as uppercut, but as been given tighter cut larger then loose cut in ptbin %d \n ---Task will use swapped Tight/Loose cuts \n ",names[ivar].Data(),i);
	  ((AliMultiDimVector*)fCutList->FindObject(mdvname.Data()))->SwapLimits(ic);
	  result = kTRUE;
	}
	if(!uppervars[ivar]&&!lowermdv){
	  printf("AliAnalysisTaskSESignificance::CheckConsistency: WARNING! %s is declared as lower cut, but as been given tighter cut smaller then loose cut in ptbin %d \n ---Task will use swapped Tight/Loose cuts \n",names[ivar].Data(),i);
	  ((AliMultiDimVector*)fCutList->FindObject(mdvname.Data()))->SwapLimits(ic);
	  result = kTRUE;
	}
	ic++;
      }
    }
  }
  return result;
}
//_________________________________________________________________
void  AliAnalysisTaskSESignificance::SetMassLimits(Float_t range, Int_t pdg){
  Float_t mass=0;
  Int_t abspdg=TMath::Abs(pdg);
  mass=TDatabasePDG::Instance()->GetParticle(abspdg)->Mass();
  fUpmasslimit = mass+range;
  fLowmasslimit = mass-range;
}
//_________________________________________________________________
void  AliAnalysisTaskSESignificance::SetMassLimits(Float_t lowlimit, Float_t uplimit){
  if(uplimit>lowlimit)
    {
      fUpmasslimit = lowlimit;
      fLowmasslimit = uplimit;
    }
}



//________________________________________________________________________
void AliAnalysisTaskSESignificance::LocalInit()
{
  // Initialization

  if(fDebug > 1) printf("AnalysisTaskSESignificance::Init() \n");

  TList *mdvList =  new TList();
  mdvList->SetOwner();
  mdvList = fCutList;
  
  PostData(2,mdvList);

  return;
}
//________________________________________________________________________
void AliAnalysisTaskSESignificance::UserCreateOutputObjects()
{
  // Create the output container
 
  if(fDebug > 1) printf("AnalysisTaskSESignificance::UserCreateOutputObjects() \n");

  // Several histograms are more conveniently managed in a TList
  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("OutputHistos");

  //same number of steps in each multiDimVectorPtBin%d !
  Int_t nHist=((AliMultiDimVector*)fCutList->FindObject("multiDimVectorPtBin0"))->GetNTotCells();
  cout<<"ncells = "<<nHist<<" n ptbins = "<<fNPtBins<<endl;
  nHist=nHist*fNPtBins;
  cout<<"Total = "<<nHist<<endl;
  for(Int_t i=0;i<nHist;i++){

    TString hisname;
    TString signame;
    TString bkgname;
    TString rflname;
    TString title;
    
    hisname.Form("hMass_%d",i);
    signame.Form("hSig_%d",i);
    bkgname.Form("hBkg_%d",i);
    rflname.Form("hRfl_%d",i);
    
    title.Form("Invariant mass;M[GeV/c^{2}];Entries");

    fMassHist[i]=new TH1F(hisname.Data(),title.Data(),100,fLowmasslimit,fUpmasslimit);
    fMassHist[i]->Sumw2();
    fOutput->Add(fMassHist[i]);

    if(fReadMC){
      fSigHist[i]=new TH1F(signame.Data(),title.Data(),100,fLowmasslimit,fUpmasslimit);
      fSigHist[i]->Sumw2();
      fOutput->Add(fSigHist[i]);
      
      fBkgHist[i]=new TH1F(bkgname.Data(),title.Data(),100,fLowmasslimit,fUpmasslimit);
      fBkgHist[i]->Sumw2();
      fOutput->Add(fBkgHist[i]);

      fRflHist[i]=new TH1F(rflname.Data(),title.Data(),100,fLowmasslimit,fUpmasslimit);
      fRflHist[i]->Sumw2();
      fOutput->Add(fRflHist[i]);

    }
  }

  fHistNEvents=new TH1F("fHistNEvents","Number of AODs scanned",4,0,4.);
  fHistNEvents->GetXaxis()->SetBinLabel(1,"nEventsAnal");
  fHistNEvents->GetXaxis()->SetBinLabel(2,"nEvNotSelected");
  fHistNEvents->GetXaxis()->SetBinLabel(3,"nCandidatesSelected");
  fHistNEvents->GetXaxis()->SetBinLabel(4,"nTotEntries Mass hists");
  fHistNEvents->GetXaxis()->SetNdivisions(1,kFALSE);
  fOutput->Add(fHistNEvents);

 
  return;
}

//________________________________________________________________________
void AliAnalysisTaskSESignificance::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event:
  // heavy flavor candidates association to MC truth
   
  AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());
  if(fDebug>2) printf("Analysing decay %d\n",fDecChannel);
  // Post the data already here
  PostData(1,fOutput);
  TClonesArray *arrayProng =0;

  if(!aod && AODEvent() && IsStandardAOD()) { 
    // In case there is an AOD handler writing a standard AOD, use the AOD 
    // event in memory rather than the input (ESD) event.    
    aod = dynamic_cast<AliAODEvent*> (AODEvent());
    // in this case the braches in the deltaAOD (AliAOD.VertexingHF.root)
    // have to taken from the AOD event hold by the AliAODExtension
    AliAODHandler* aodHandler = (AliAODHandler*) 
      ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
    if(aodHandler->GetExtensions()) {
      
      AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
      AliAODEvent *aodFromExt = ext->GetAOD();
   
      
      switch(fDecChannel){
      case 0:
	arrayProng=(TClonesArray*)aodFromExt->GetList()->FindObject("Charm3Prong");
	break; 
      case 1:
	arrayProng=(TClonesArray*)aodFromExt->GetList()->FindObject("D0toKpi");
	break; 
      case 2:
	arrayProng=(TClonesArray*)aodFromExt->GetList()->FindObject("Dstar");
	break; 
      case 3:
	arrayProng=(TClonesArray*)aodFromExt->GetList()->FindObject("Charm3Prong");
	break; 
      case 4:
	arrayProng=(TClonesArray*)aodFromExt->GetList()->FindObject("Charm4Prong");
	break; 
      case 5:
	arrayProng=(TClonesArray*)aodFromExt->GetList()->FindObject("Charm3Prong");
	break; 
      }
    }
  } else {
    switch(fDecChannel){
    case 0:
      arrayProng=(TClonesArray*)aod->GetList()->FindObject("Charm3Prong");
      break; 
    case 1:
      arrayProng=(TClonesArray*)aod->GetList()->FindObject("D0toKpi");
      break; 
    case 2:
      arrayProng=(TClonesArray*)aod->GetList()->FindObject("Dstar");
      break; 
    case 3:
      arrayProng=(TClonesArray*)aod->GetList()->FindObject("Charm3Prong");
      break; 
    case 4:
      arrayProng=(TClonesArray*)aod->GetList()->FindObject("Charm4Prong");
      break; 
    case 5:
      arrayProng=(TClonesArray*)aod->GetList()->FindObject("Charm3Prong");
      break; 
    }
  }
  if(!arrayProng) {
    printf("AliAnalysisTaskSESignificance::UserExec: branch not found!\n");
    return;
  }
  
  TClonesArray *arrayMC=0;
  AliAODMCHeader *mcHeader=0;

  // load MC particles
  if(fReadMC){
    
    arrayMC =  (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
    if(!arrayMC) {
      printf("AliAnalysisTaskSESignificance::UserExec: MC particles branch not found!\n");
      //    return;
    }
    
    // load MC header
    mcHeader =  (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if(!mcHeader) {
      printf("AliAnalysisTaskSESignificance::UserExec: MC header branch not found!\n");
      return;
    }
  }

  fHistNEvents->Fill(0); // count event

  AliAODRecoDecayHF *d=0;
  Int_t nprongs=1;
  Int_t *pdgdaughters=0x0;
  Int_t prongpdg=211;
  

  switch(fDecChannel){
  case 0:
    //D+
    pdgdaughters =new Int_t[3];
    pdgdaughters[0]=211;//pi
    pdgdaughters[1]=321;//K
    pdgdaughters[2]=211;//pi
    nprongs=3;
    prongpdg=411;
    break;
  case 1:
    //D0
    pdgdaughters =new Int_t[2];
    pdgdaughters[0]=211;//pi 
    pdgdaughters[1]=321;//K
    nprongs=2;
    prongpdg=421;
    break;
  case 2:
    //D*
    pdgdaughters =new Int_t[3];
    pdgdaughters[1]=211;//pi
    pdgdaughters[0]=321;//K
    pdgdaughters[2]=211;//pi (soft?)
    nprongs=3;
    prongpdg=413;
    break;
  case 3:
    //Ds
    pdgdaughters =new Int_t[3];
    pdgdaughters[0]=211;//pi
    pdgdaughters[1]=321;//K
    pdgdaughters[2]=321;//K
    nprongs=3;
    prongpdg=431;
    break;
  case 4:
    //D0
    pdgdaughters =new Int_t[4];
    pdgdaughters[0]=321;
    pdgdaughters[1]=211;
    pdgdaughters[2]=211;
    pdgdaughters[3]=211;
    nprongs=4;
    prongpdg=421;
    break;
  case 5:
    //Lambda_c
    pdgdaughters =new Int_t[3];
    //check order
    pdgdaughters[0]=2212;//p
    pdgdaughters[1]=321;//K
    pdgdaughters[2]=211;//pi
    nprongs=3;
    prongpdg=4122;
    break;
  }

  Int_t nHistpermv=((AliMultiDimVector*)fCutList->FindObject("multiDimVectorPtBin0"))->GetNTotCells();
  Int_t nProng = arrayProng->GetEntriesFast();
  if(fDebug>1) printf("Number of D2H: %d\n",nProng);

  if(!fRDCuts->IsEventSelected(aod)){
    fHistNEvents->Fill(1);
    return;
  }

  for (Int_t iProng = 0; iProng < nProng; iProng++) {
    d=(AliAODRecoDecayHF*)arrayProng->UncheckedAt(iProng);
    Double_t invMass=d->InvMass(nprongs,(UInt_t*)pdgdaughters);
    Double_t invMassC=0;
    
    Int_t isSelected=fRDCuts->IsSelected(d,fSelectionlevel);
    
    if(isSelected) {
      fHistNEvents->Fill(2); // count selected with loosest cuts
      if(fDebug>1) printf("+++++++Is Selected\n");
    
      const Int_t nvars=fRDCuts->GetNVarsForOpt();
      Float_t *vars = new Float_t[nvars];
      Int_t nVals=0;
      
      fRDCuts->GetCutVarsForOpt(d,vars,nvars,pdgdaughters);
      Int_t ptbin=fRDCuts->PtBin(d->Pt());
      if(ptbin==-1) continue;
      TString mdvname=Form("multiDimVectorPtBin%d",ptbin);
      ULong64_t *addresses = ((AliMultiDimVector*)fCutList->FindObject(mdvname.Data()))->GetGlobalAddressesAboveCuts(vars,(Float_t)d->Pt(),nVals);
      if(fDebug>1)printf("nvals = %d\n",nVals);
      for(Int_t ivals=0;ivals<nVals;ivals++){
	if(addresses[ivals]>=((AliMultiDimVector*)fCutList->FindObject(mdvname.Data()))->GetNTotCells()){
	  if (fDebug>1) printf("Overflow!!\n");
	  return;
	}
	fMassHist[ptbin*nHistpermv+addresses[ivals]]->Fill(invMass);
	fHistNEvents->Fill(3);
	if(fReadMC){
	  Int_t lab=-1;
	  lab = d->MatchToMC(prongpdg,arrayMC,nprongs,pdgdaughters);
	  if(lab>=0){ //signal
	    AliAODMCParticle *dMC = (AliAODMCParticle*)arrayMC->At(lab);
	    Int_t pdgMC = dMC->GetPdgCode();
	    
	    if(pdgMC==+prongpdg) fSigHist[ptbin*nHistpermv+addresses[ivals]]->Fill(invMass);
	    else fRflHist[ptbin*nHistpermv+addresses[ivals]]->Fill(invMass);
	  }
	  else{ //background
	    fBkgHist[ptbin*nHistpermv+addresses[ivals]]->Fill(invMass);
	  }
	
	}
      }


      if(fDecChannel==AliAnalysisTaskSESignificance::kD0toKpi || fDecChannel==AliAnalysisTaskSESignificance::kLambdactopKpi){ //repeat the cycle to checks cuts passed as D0bar
	if(fDecChannel==AliAnalysisTaskSESignificance::kD0toKpi){
	  pdgdaughters[0]=321;
	  pdgdaughters[1]=211; 
	}
	if(fDecChannel==AliAnalysisTaskSESignificance::kLambdactopKpi){
	  pdgdaughters[0]=211;
	  pdgdaughters[1]=321;
 	  pdgdaughters[2]=2212;
	}

	invMass=d->InvMass(nprongs,(UInt_t*)pdgdaughters);
	fRDCuts->GetCutVarsForOpt(d,vars,nvars,pdgdaughters);
	addresses = ((AliMultiDimVector*)fCutList->FindObject(mdvname.Data()))->GetGlobalAddressesAboveCuts(vars,(Float_t)d->Pt(),nVals);
	for(Int_t ivals=0;ivals<nVals;ivals++){
	  fMassHist[ptbin*nHistpermv+addresses[ivals]]->Fill(invMassC);
	  if(fReadMC){
	    Int_t lab=-1;
	    lab = d->MatchToMC(-prongpdg,arrayMC,nprongs,pdgdaughters);
	    if(lab>=0){ //signal
	      AliAODMCParticle *dMC = (AliAODMCParticle*)arrayMC->At(lab);
	      Int_t pdgMC = dMC->GetPdgCode();
	      if(pdgMC==-prongpdg) fSigHist[ptbin*nHistpermv+addresses[ivals]]->Fill(invMassC);
	      else fRflHist[ptbin*nHistpermv+addresses[ivals]]->Fill(invMassC);
	    }
	    else{ //background
	      fBkgHist[ptbin*nHistpermv+addresses[ivals]]->Fill(invMassC);
	    }
	    
	  }
	}
      }
  
    }// end if selected

    
  }

  PostData(1,fOutput);    
  return;
}


//________________________________________________________________________
void AliAnalysisTaskSESignificance::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
  if(fDebug > 1) printf("AnalysisTaskSESignificance: Terminate() \n");


  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {     
    printf("ERROR: fOutput not available\n");
    return;
  }
 
  fCutList = dynamic_cast<TList*> (GetOutputData(2));
  if (!fCutList) {     
    printf("ERROR: fCutList not available\n");
    return;
  }

  AliMultiDimVector* mdvtmp=(AliMultiDimVector*)fCutList->FindObject("multiDimVectorPtBin0");
  if (!mdvtmp){
    cout<<"multidimvec not found in TList"<<endl;
    fCutList->ls();
    return;
  }
  Int_t nHist=mdvtmp->GetNTotCells();
  TCanvas *c1=new TCanvas("c1","Invariant mass distribution - loose cuts",500,500);
  Bool_t drawn=kFALSE;
  for(Int_t i=0;i<nHist;i++){

    TString hisname;
    TString signame;
    TString bkgname;
    TString rflname;
    
    hisname.Form("hMass_%d",i);
    signame.Form("hSig_%d",i);
    bkgname.Form("hBkg_%d",i);
    rflname.Form("hRfl_%d",i);

    fMassHist[i]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));

    if (!drawn && fMassHist[i]->GetEntries() > 0){
      c1->cd();
      fMassHist[i]->Draw();
      drawn=kTRUE;
    }
    
    if(fReadMC){
      fSigHist[i]=dynamic_cast<TH1F*>(fOutput->FindObject(signame.Data()));
      fBkgHist[i]=dynamic_cast<TH1F*>(fOutput->FindObject(bkgname.Data()));
      fRflHist[i]=dynamic_cast<TH1F*>(fOutput->FindObject(bkgname.Data()));
    }
    
  }
  

 
  
  return;
}

