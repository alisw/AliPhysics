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

/* $Id$ */

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
  fUseSelBit(kFALSE),
  fBFeedDown(kBoth),
  fDecChannel(0),
  fPDGmother(0),
  fNProngs(0),
  fBranchName(""),
  fSelectionlevel(0),
  fNVars(0),
  fNBins(100),
  fPartOrAndAntiPart(0),
  fDsChannel(0)
{
  // Default constructor
  SetPDGCodes();
  SetDsChannel(kPhi);

  for(Int_t i=0;i<4;i++) fPDGdaughters[i]=0;
  for(Int_t i=0;i<2;i++) {fPDGD0ToKpi[i]=0;fPDGDStarToD0pi[i]=0;}
  for(Int_t i=0;i<kMaxCutVar;i++) fVars[i]=0.;
  for(Int_t i=0;i<kMaxNHist;i++) {
    fMassHist[i]=0;
    fSigHist[i]=0;
    fBkgHist[i]=0;
    fRflHist[i]=0;
  }

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
  fUseSelBit(kFALSE),
  fBFeedDown(kBoth),
  fDecChannel(decaychannel),
  fPDGmother(0),
  fNProngs(0),
  fBranchName(""),
  fSelectionlevel(selectionlevel),
  fNVars(0),
  fNBins(100),
  fPartOrAndAntiPart(0),
  fDsChannel(0)
{

  SetPDGCodes();
  SetDsChannel(kPhi);
  if (fDecChannel!=2) SetMassLimits(0.15,fPDGmother); //check range
  else {
    Float_t min = 0.13;
    Float_t max = 0.19;
    SetMassLimits(min, max);
  }
  fNPtBins=fRDCuts->GetNPtBins();

  fNVars=fRDCuts->GetNVarsForOpt();
  if(fNVars>kMaxCutVar) AliFatal(Form("Too large number of cut variables, maximum is %d",kMaxCutVar));
  
  if(fDebug>1)fRDCuts->PrintAll();
   // Output slot #1 writes into a TList container
  DefineOutput(1,TList::Class());  //My private output
  DefineOutput(2,TList::Class());
  DefineOutput(3,AliRDHFCuts::Class()); //class of the cuts
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
void AliAnalysisTaskSESignificance::SetPDGCodes(){
  // sets channel dependent quantities
 
  switch(fDecChannel){
  case 0:
    //D+
    fPDGmother=411;
    fNProngs=3;
    fPDGdaughters[0]=211;//pi
    fPDGdaughters[1]=321;//K
    fPDGdaughters[2]=211;//pi
    fPDGdaughters[3]=0; //empty
    fBranchName="Charm3Prong";
    break;
  case 1:
    //D0
    fPDGmother=421;
    fNProngs=2;
    fPDGdaughters[0]=211;//pi 
    fPDGdaughters[1]=321;//K
    fPDGdaughters[2]=0; //empty
    fPDGdaughters[3]=0; //empty
    fBranchName="D0toKpi";
    break;
  case 2:
    //D*
    fPDGmother=413;
    fNProngs=3;
    fPDGdaughters[1]=211;//pi
    fPDGdaughters[0]=321;//K
    fPDGdaughters[2]=211;//pi (soft?)
    fPDGdaughters[3]=0; //empty
    fBranchName="Dstar";
    break;
  case 3:
    //Ds
    fPDGmother=431;
    fNProngs=3;
    fPDGdaughters[0]=321;//K
    fPDGdaughters[1]=321;//K
    fPDGdaughters[2]=211;//pi
    fPDGdaughters[3]=0; //empty
    fBranchName="Charm3Prong";
    break;
  case 4:
    //D0 in 4 prongs
    fPDGmother=421;
    fNProngs=4;
    fPDGdaughters[0]=321;
    fPDGdaughters[1]=211;
    fPDGdaughters[2]=211;
    fPDGdaughters[3]=211;
    fBranchName="Charm4Prong";
    break;
  case 5:
    //Lambda_c
    fPDGmother=4122;
    fNProngs=3;
    fPDGdaughters[0]=2212;//p
    fPDGdaughters[1]=321;//K
    fPDGdaughters[2]=211;//pi
    fPDGdaughters[3]=0; //empty
    fBranchName="Charm3Prong";
    break;
  }
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
	  AliFatal(Form("tight and loose cut for optimization variable number %d are the same in ptbin %d\n",ic,i));
	  return kFALSE;
	}
	Bool_t lowermdv = ((AliMultiDimVector*)fCutList->FindObject(mdvname.Data()))->GetGreaterThan(ic);
	if(uppervars[ivar]&&lowermdv){
	  AliFatal(Form("%s is declared as uppercut, but as been given tighter cut larger then loose cut in ptbin %d \n ---please check your cuts \n ",names[ivar].Data(),i));
	  return kFALSE;
	}
	if(!uppervars[ivar]&&!lowermdv){
	  AliFatal(Form("%s is declared as lower cut, but as been given tighter cut smaller then loose cut in ptbin %d \n ---please check your cuts \n",names[ivar].Data(),i));
	  return kFALSE;
	}
	ic++;
      }
    }
  }
  return result;
}
//_________________________________________________________________
void AliAnalysisTaskSESignificance::SetBFeedDown(FeedDownEnum flagB){
  if(fReadMC)fBFeedDown=flagB;
  else {
    AliInfo("B feed down not allowed without MC info\n");
    fBFeedDown=kBoth;
  }
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
      fUpmasslimit = uplimit;
      fLowmasslimit = lowlimit;
    }
}



//________________________________________________________________________
void AliAnalysisTaskSESignificance::LocalInit()
{
  // Initialization

  if(fDebug > 1) printf("AnalysisTaskSESignificance::Init() \n");

  switch(fDecChannel){
  case 0:
    {
      AliRDHFCutsDplustoKpipi* copycut=new AliRDHFCutsDplustoKpipi(*(static_cast<AliRDHFCutsDplustoKpipi*>(fRDCuts)));
      // Post the data
      PostData(3,copycut);
    }
    break;
  case 1:
    {
      AliRDHFCutsD0toKpi* copycut=new AliRDHFCutsD0toKpi(*(static_cast<AliRDHFCutsD0toKpi*>(fRDCuts)));
      // Post the data
      PostData(3,copycut);
    }
    break;
  case 2:
    {
      AliRDHFCutsDStartoKpipi* copycut=new AliRDHFCutsDStartoKpipi(*(static_cast<AliRDHFCutsDStartoKpipi*>(fRDCuts)));
      // Post the data
      PostData(3,copycut);
    }
    break;
  case 3:
    {
      AliRDHFCutsDstoKKpi* copycut=new AliRDHFCutsDstoKKpi(*(static_cast<AliRDHFCutsDstoKKpi*>(fRDCuts)));
      // Post the data
      PostData(3,copycut);
    }
    break;
  case 4:
    {
      AliRDHFCutsD0toKpipipi* copycut=new AliRDHFCutsD0toKpipipi(*(static_cast<AliRDHFCutsD0toKpipipi*>(fRDCuts)));
      // Post the data
      PostData(3,copycut);
    }
    break;
  case 5:
    {
      AliRDHFCutsLctopKpi* copycut=new AliRDHFCutsLctopKpi(*(static_cast<AliRDHFCutsLctopKpi*>(fRDCuts)));
      // Post the data
      PostData(3,copycut);
    }
    break;

  default:
    return;
  }

  TList *mdvList =  new TList();
  mdvList->SetOwner();
  mdvList = fCutList;
  
  PostData(2,mdvList);


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

    fMassHist[i]=new TH1F(hisname.Data(),title.Data(),fNBins,fLowmasslimit,fUpmasslimit);
    fMassHist[i]->Sumw2();
    fOutput->Add(fMassHist[i]);

    if(fReadMC){
      fSigHist[i]=new TH1F(signame.Data(),title.Data(),fNBins,fLowmasslimit,fUpmasslimit);
      fSigHist[i]->Sumw2();
      fOutput->Add(fSigHist[i]);
      
      fBkgHist[i]=new TH1F(bkgname.Data(),title.Data(),fNBins,fLowmasslimit,fUpmasslimit);
      fBkgHist[i]->Sumw2();
      fOutput->Add(fBkgHist[i]);

      if(fDecChannel != AliAnalysisTaskSESignificance::kDplustoKpipi){
	fRflHist[i]=new TH1F(rflname.Data(),title.Data(),fNBins,fLowmasslimit,fUpmasslimit);
	fRflHist[i]->Sumw2();
	fOutput->Add(fRflHist[i]);
      }
    }
  }

  fHistNEvents=new TH1F("fHistNEvents","Number of AODs scanned",8,-0.5,7.5);
  fHistNEvents->GetXaxis()->SetBinLabel(1,"nEventsAnal");
  fHistNEvents->GetXaxis()->SetBinLabel(2,"nEvSelected (vtx)");
  fHistNEvents->GetXaxis()->SetBinLabel(3,"nCandidatesSelected");
  fHistNEvents->GetXaxis()->SetBinLabel(4,"nTotEntries Mass hists");
  fHistNEvents->GetXaxis()->SetBinLabel(5,"Pile-up Rej");
  fHistNEvents->GetXaxis()->SetBinLabel(6,"N. of 0SMH");
  if(fReadMC){
    fHistNEvents->GetXaxis()->SetBinLabel(7,"MC Cand from c");
    fHistNEvents->GetXaxis()->SetBinLabel(8,"MC Cand from b");
  } else{
    fHistNEvents->GetXaxis()->SetBinLabel(7,"N candidates");
  }
  fHistNEvents->GetXaxis()->SetNdivisions(1,kFALSE);
  fOutput->Add(fHistNEvents);

 
  PostData(1,fOutput);    
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
      arrayProng=(TClonesArray*)aodFromExt->GetList()->FindObject(fBranchName.Data());
      
   }
  } else if(aod) {
    arrayProng=(TClonesArray*)aod->GetList()->FindObject(fBranchName.Data());
  }
  if(!aod || !arrayProng) {
    AliError("AliAnalysisTaskSESignificance::UserExec:Branch not found!\n");
    return;
  }
  
  // fix for temporary bug in ESDfilter 
  // the AODs with null vertex pointer didn't pass the PhysSel
  if(!aod->GetPrimaryVertex() || TMath::Abs(aod->GetMagneticField())<0.001) return;
  TClonesArray *arrayMC=0;
  AliAODMCHeader *mcHeader=0;

  // load MC particles
  if(fReadMC){
    
    arrayMC =  (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
    if(!arrayMC) {
      AliError("AliAnalysisTaskSESignificance::UserExec:MC particles branch not found!\n");
      return;
    }
    
    // load MC header
    mcHeader =  (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if(!mcHeader) {
      AliError("AliAnalysisTaskSESignificance::UserExec:MC header branch not found!\n");
      return;
    }
  }

  fHistNEvents->Fill(0); // count event

  AliAODRecoDecayHF *d=0;
  


  Int_t nHistpermv=((AliMultiDimVector*)fCutList->FindObject("multiDimVectorPtBin0"))->GetNTotCells();
  Int_t nProng = arrayProng->GetEntriesFast();
  if(fDebug>1) printf("Number of D2H: %d\n",nProng);

  // trigger class for PbPb C0SMH-B-NOPF-ALLNOTRD, C0SMH-B-NOPF-ALL
  TString trigclass=aod->GetFiredTriggerClasses();
  if(trigclass.Contains("C0SMH-B-NOPF-ALLNOTRD") || trigclass.Contains("C0SMH-B-NOPF-ALL")) fHistNEvents->Fill(5);

  if(fRDCuts->IsEventSelected(aod)) {
    fHistNEvents->Fill(1);
  }else{
    if(fRDCuts->GetWhyRejection()==1) // rejected for pileup
      fHistNEvents->Fill(4);
    return;
  }
  
  for (Int_t iProng = 0; iProng < nProng; iProng++) {
    fHistNEvents->Fill(6);
    d=(AliAODRecoDecayHF*)arrayProng->UncheckedAt(iProng);

    AliAODRecoCascadeHF* DStarToD0pi = NULL;
    AliAODRecoDecayHF2Prong* D0Particle = NULL;
    fPDGDStarToD0pi[0] = 421; fPDGDStarToD0pi[1] = 211;
    fPDGD0ToKpi[0] = 321; fPDGD0ToKpi[1] = 211;

    Bool_t isSelBit=kTRUE;
    if(fUseSelBit){
      if(fDecChannel==0) {
	isSelBit=d->HasSelectionBit(AliRDHFCuts::kDplusCuts);
      }else{ 
	if(fDecChannel==1) {
	  isSelBit=d->HasSelectionBit(AliRDHFCuts::kD0toKpiCuts);
	}else{
	  if(fDecChannel==2) {
	    isSelBit=d->HasSelectionBit(AliRDHFCuts::kDstarCuts);
	  }else{
	    if(fDecChannel==3) {
	      isSelBit=d->HasSelectionBit(AliRDHFCuts::kDsCuts);
	    }else{
	      if(fDecChannel==5) isSelBit=d->HasSelectionBit(AliRDHFCuts::kLcCuts);
	    }
	  }
	}
      }
    }
    if(!isSelBit) continue; 

    if (fDecChannel==2) {
      DStarToD0pi = (AliAODRecoCascadeHF*)arrayProng->At(iProng);
      if (!DStarToD0pi->GetSecondaryVtx()) continue;
      D0Particle = (AliAODRecoDecayHF2Prong*)DStarToD0pi->Get2Prong();
      if (!D0Particle) continue;
    }
    
    Bool_t isFidAcc = fRDCuts->IsInFiducialAcceptance(d->Pt(),d->Y(fPDGmother));
    Int_t isSelected=fRDCuts->IsSelected(d,fSelectionlevel,aod);

    if(fReadMC && fBFeedDown!=kBoth && isSelected){
      Int_t labD = d->MatchToMC(fPDGmother,arrayMC,fNProngs,fPDGdaughters);
      if(labD>=0){
	AliAODMCParticle *partD = (AliAODMCParticle*)arrayMC->At(labD);
	Int_t pdgGranma = CheckOrigin(partD, arrayMC);
	Int_t abspdgGranma = TMath::Abs(pdgGranma);
	if ((abspdgGranma > 500 && abspdgGranma < 600) || (abspdgGranma > 5000 && abspdgGranma < 6000)) {
	  //feed down particle
	  AliDebug(2,Form("Particle has a b-meson, or b-baryon mother (pdg code mother = %d )--> not coming from a c-quark, skipping...", pdgGranma));
	  fHistNEvents->Fill(7);
	  if(fBFeedDown==kCharmOnly) isSelected=kFALSE; //from beauty
	}
	else { 
	  //prompt particle
	  fHistNEvents->Fill(6);
	  if(fBFeedDown==kBeautyOnly)isSelected=kFALSE;
	} 
      }
    }
    
    if(isSelected&&isFidAcc) {
      fHistNEvents->Fill(2); // count selected with loosest cuts
      if(fDebug>1) printf("+++++++Is Selected\n");
    
      Int_t nVals=0;
      if(fDecChannel==3) SetPDGdaughterDstoKKpi();
      fRDCuts->GetCutVarsForOpt(d,fVars,fNVars,fPDGdaughters,aod);
      Int_t ptbin=fRDCuts->PtBin(d->Pt());
      if(ptbin==-1) continue;
      TString mdvname=Form("multiDimVectorPtBin%d",ptbin);
      AliMultiDimVector* muvec=(AliMultiDimVector*)fCutList->FindObject(mdvname.Data());

      ULong64_t *addresses = muvec->GetGlobalAddressesAboveCuts(fVars,(Float_t)d->Pt(),nVals);
      if(fDebug>1)printf("nvals = %d\n",nVals);
      for(Int_t ivals=0;ivals<nVals;ivals++){
	if(addresses[ivals]>=muvec->GetNTotCells()){
	  if (fDebug>1) printf("Overflow!!\n");
	  delete [] addresses;
	  return;
	}
	
	fHistNEvents->Fill(3);
	
	//fill the histograms with the appropriate method
	switch (fDecChannel){
	case 0:
	  FillDplus(d,arrayMC,(Int_t)(ptbin*nHistpermv+addresses[ivals]),isSelected);
	  break;
	case 1:
	  FillD02p(d,arrayMC,(Int_t)(ptbin*nHistpermv+addresses[ivals]),isSelected);
	  break;
	case 2:
	  FillDstar(DStarToD0pi,arrayMC,(Int_t)(ptbin*nHistpermv+addresses[ivals]),isSelected);
	  break;
	case 3:
	  if(isSelected&1){
	    FillDs(d,arrayMC,(Int_t)(ptbin*nHistpermv+addresses[ivals]),isSelected,1);
	  }
	  break;
	case 4:
	  FillD04p(d,arrayMC,(Int_t)(ptbin*nHistpermv+addresses[ivals]),isSelected);
	  break;
	case 5:
	  FillLambdac(d,arrayMC,(Int_t)(ptbin*nHistpermv+addresses[ivals]),isSelected);
	  break;
	default:
	  break;
	}
	
      }
      
      if (fDecChannel==3 && isSelected&2){
	SetPDGdaughterDstopiKK();
	nVals=0;
	fRDCuts->GetCutVarsForOpt(d,fVars,fNVars,fPDGdaughters,aod);
	delete [] addresses;
	addresses = muvec->GetGlobalAddressesAboveCuts(fVars,(Float_t)d->Pt(),nVals);
	if(fDebug>1)printf("nvals = %d\n",nVals);
	for(Int_t ivals=0;ivals<nVals;ivals++){
	  if(addresses[ivals]>=muvec->GetNTotCells()){
	    if (fDebug>1) printf("Overflow!!\n");
	    delete [] addresses;	    
	    return;
	  }
	  FillDs(d,arrayMC,(Int_t)(ptbin*nHistpermv+addresses[ivals]),isSelected,0);	  
	  
	}

      }
      
      delete [] addresses;
    }// end if selected
    
  }
  

  PostData(1,fOutput);    
  return;
}

//***************************************************************************

// Methods used in the UserExec


//********************************************************************************************

//Methods to fill the istograms with MC information, one for each candidate
//NB: the implementation for each candidate is responsibility of the corresponding developer

void AliAnalysisTaskSESignificance::FillDplus(AliAODRecoDecayHF* d,TClonesArray *arrayMC,Int_t index,Int_t isSel){
    //D+ channel
  if(!isSel){
    AliError("Candidate not selected\n");
    return;
  }

  if(fPartOrAndAntiPart*d->GetCharge()<0)return;
  Int_t pdgdaughters[3] = {211,321,211};
  Double_t mass=d->InvMass(3,(UInt_t*)pdgdaughters);

  fMassHist[index]->Fill(mass);


  if(fReadMC){
    Int_t lab=-1;
    lab = d->MatchToMC(411,arrayMC,3,pdgdaughters);
    if(lab>=0){ //signal
      fSigHist[index]->Fill(mass);
    } else{ //background
      fBkgHist[index]->Fill(mass);
    } 
  }   
}


void AliAnalysisTaskSESignificance::FillD02p(AliAODRecoDecayHF* d,TClonesArray *arrayMC,Int_t index,Int_t isSel){

  //D0->Kpi channel

  //mass histograms
  Int_t pdgdaughtersD0[2]={211,321};//pi,K 
  Int_t pdgdaughtersD0bar[2]={321,211};//K,pi 
  Double_t masses[2];
  masses[0]=d->InvMass(fNProngs,(UInt_t*)pdgdaughtersD0); //D0
  masses[1]=d->InvMass(fNProngs,(UInt_t*)pdgdaughtersD0bar); //D0bar

  if((isSel==1 || isSel==3) && fPartOrAndAntiPart>=0) fMassHist[index]->Fill(masses[0]);
  if(isSel>=2 && fPartOrAndAntiPart<=0) fMassHist[index]->Fill(masses[1]);



  //MC histograms
  if(fReadMC){

    Int_t matchtoMC=-1;

    //D0
    Int_t pdgdaughters[2];
    pdgdaughters[0]=211;//pi 
    pdgdaughters[1]=321;//K

    matchtoMC = d->MatchToMC(fPDGmother,arrayMC,fNProngs,pdgdaughters);

    Int_t prongPdgPlus=421,prongPdgMinus=(-1)*421;
    if((isSel==1 || isSel==3) && fPartOrAndAntiPart>=0){ //D0
      if(matchtoMC>=0){
	AliAODMCParticle *dMC = (AliAODMCParticle*)arrayMC->At(matchtoMC);
	Int_t pdgMC = dMC->GetPdgCode();
	
	if(pdgMC==prongPdgPlus) fSigHist[index]->Fill(masses[0]);
	else fRflHist[index]->Fill(masses[0]);
	
      } else fBkgHist[index]->Fill(masses[0]);
      
    }
    if(isSel>=2 && fPartOrAndAntiPart<=0){ //D0bar
      if(matchtoMC>=0){
	AliAODMCParticle *dMC = (AliAODMCParticle*)arrayMC->At(matchtoMC);
	Int_t pdgMC = dMC->GetPdgCode();
	
	if(pdgMC==prongPdgMinus) fSigHist[index]->Fill(masses[1]);
	else fRflHist[index]->Fill(masses[1]);
      } else fBkgHist[index]->Fill(masses[1]);
    }
  }
}

void AliAnalysisTaskSESignificance::FillDstar(AliAODRecoCascadeHF* dstarD0pi,TClonesArray *arrayMC,Int_t index,Int_t isSel){
    //D* channel

    AliInfo("Dstar selected\n");
    
    Double_t mass = dstarD0pi->DeltaInvMass();

    if((isSel==1 || isSel==3) && fPartOrAndAntiPart>=0) fMassHist[index]->Fill(mass);
    if(isSel>=2 && fPartOrAndAntiPart<=0) fMassHist[index]->Fill(mass);
	
    if(fReadMC) {
       Int_t matchtoMC = -1; 
       matchtoMC = dstarD0pi->MatchToMC(413,421,fPDGDStarToD0pi, fPDGD0ToKpi, arrayMC);

       Int_t prongPdgDStarPlus=413,prongPdgDStarMinus=(-1)*413;
	
       if ((isSel==1 || isSel==3) && fPartOrAndAntiPart>=0) { 
          //D*+
      	  if(matchtoMC>=0) {
             AliAODMCParticle *dMC = (AliAODMCParticle*)arrayMC->At(matchtoMC);
	     Int_t pdgMC = dMC->GetPdgCode();
	
	     if (pdgMC==prongPdgDStarPlus) fSigHist[index]->Fill(mass);
	     else {
	        dstarD0pi->SetCharge(-1*dstarD0pi->GetCharge());
		mass =	dstarD0pi->DeltaInvMass();
		fRflHist[index]->Fill(mass);
		dstarD0pi->SetCharge(-1*dstarD0pi->GetCharge());
	      }	
      	    } 
	  else fBkgHist[index]->Fill(mass);
        }
    
	if (isSel>=2 && fPartOrAndAntiPart<=0) { 
           //D*-
      	   if (matchtoMC>=0) {
	      AliAODMCParticle *dMC = (AliAODMCParticle*)arrayMC->At(matchtoMC);
	      Int_t pdgMC = dMC->GetPdgCode();
	
	      if (pdgMC==prongPdgDStarMinus) fSigHist[index]->Fill(mass);
	      else {
		 dstarD0pi->SetCharge(-1*dstarD0pi->GetCharge());
		 mass =	dstarD0pi->DeltaInvMass();
		 fRflHist[index]->Fill(mass);
		 dstarD0pi->SetCharge(-1*dstarD0pi->GetCharge());
	      }
      	    } 
	    else fBkgHist[index]->Fill(mass);
    	  }
       }

}


void AliAnalysisTaskSESignificance::FillDs(AliAODRecoDecayHF* d,TClonesArray *arrayMC,Int_t index,Int_t isSel,Int_t optDecay){

  // Fill Ds histos


  Int_t pdgDsKKpi[3]={321,321,211};//K,K,pi 
  Int_t pdgDspiKK[3]={211,321,321};//pi,K,K 
  Double_t masses[2];  
  masses[0]=d->InvMass(fNProngs,(UInt_t*)pdgDsKKpi); //Ds
  masses[1]=d->InvMass(fNProngs,(UInt_t*)pdgDspiKK); //Dsbar 

  Int_t labDs=-1;
  if(fReadMC){
    labDs = d->MatchToMC(431,arrayMC,3,pdgDsKKpi);
  }

  Int_t isKKpi=isSel&1;
  Int_t ispiKK=isSel&2;
  Int_t isPhiKKpi=isSel&4;
  Int_t isPhipiKK=isSel&8;
  Int_t isK0starKKpi=isSel&16;
  Int_t isK0starpiKK=isSel&32;


  if(fDsChannel==kPhi && (isPhiKKpi==0 && isPhipiKK==0)) return;
  if(fDsChannel==kK0star && (isK0starKKpi==0 && isK0starpiKK==0)) return;
   
  if (optDecay==1){ 
    if(isKKpi  && fPartOrAndAntiPart*d->GetCharge()>=0) {
      if(fDsChannel==kPhi && isPhiKKpi==0) return;
      if(fDsChannel==kK0star && isK0starKKpi==0) return;
      
      fMassHist[index]->Fill(masses[0]); 
      
      if(fReadMC){
	if(labDs>=0){
	  Int_t labDau0=((AliAODTrack*)d->GetDaughter(0))->GetLabel();
	  AliAODMCParticle* p=(AliAODMCParticle*)arrayMC->UncheckedAt(labDau0);
	  Int_t pdgCode0=TMath::Abs(p->GetPdgCode());
	  if(pdgCode0==321) {
	    fSigHist[index]->Fill(masses[0]); //signal
	  }else{
	    fRflHist[index]->Fill(masses[0]); //Reflected signal
	  }
	}else{
	  fBkgHist[index]->Fill(masses[0]); // Background
	}
      }
    }
  }
  
  if (optDecay==0){ 
    if(ispiKK && fPartOrAndAntiPart*d->GetCharge()>=0){
      if(fDsChannel==kPhi && isPhipiKK==0) return;
      if(fDsChannel==kK0star && isK0starpiKK==0) return;
      
      fMassHist[index]->Fill(masses[1]);
      
      if(fReadMC){
	if(labDs>=0){
	  Int_t labDau0=((AliAODTrack*)d->GetDaughter(0))->GetLabel();
	  AliAODMCParticle* p=(AliAODMCParticle*)arrayMC->UncheckedAt(labDau0);
	  Int_t pdgCode0=TMath::Abs(p->GetPdgCode());	
	  if(pdgCode0==211) {	  	  
	    fSigHist[index]->Fill(masses[1]);
	  }else{
	    fRflHist[index]->Fill(masses[1]);
	  }
	}else{
	  fBkgHist[index]->Fill(masses[1]);
	}
      }
    }
  }
}

void AliAnalysisTaskSESignificance::FillD04p(AliAODRecoDecayHF* /*d*/,TClonesArray */*arrayMC*/,Int_t /*index*/,Int_t /*matchtoMC*/){
  //D0->Kpipipi channel
  AliInfo("D0 in 4 prongs channel not implemented\n");

}

void AliAnalysisTaskSESignificance::FillLambdac(AliAODRecoDecayHF* d,TClonesArray* arrayMC,Int_t index,Int_t isSel){

  // Mass hypothesis
  Double_t masses[2];
  Int_t pdgdaughtersLc[3]={2212,321,211}; //p,K,pi
  masses[0]=d->InvMass(fNProngs,(UInt_t*)pdgdaughtersLc);
  Int_t pdgdaughtersLc2[3]={211,321,2212}; //pi,K,p
  masses[1]=d->InvMass(fNProngs,(UInt_t*)pdgdaughtersLc2);


  if(fPartOrAndAntiPart==0 || fPartOrAndAntiPart==d->GetCharge()) {
    
    // isSel=1 : p K pi ; isSel=2 : pi K p ;
    if(isSel==1 || isSel==3) fMassHist[index]->Fill(masses[0]);
    if(isSel>=2) fMassHist[index]->Fill(masses[1]);
    
    // Check the MC truth
    if(fReadMC){

      Int_t ispKpi = 0;
      Int_t ispiKp = 0;
      ispKpi = isSel&1;
      ispiKp = isSel&2;
      Int_t matchtoMC = -1;
      //
      Int_t pPDG = 2212; // p 
      Int_t kPDG = 321;  // K
      Int_t piPDG = 211; // pi
      Int_t absPdgMom = 4122;
      matchtoMC = d->MatchToMC(absPdgMom,arrayMC,fNProngs,pdgdaughtersLc);

      if(matchtoMC>=0){

	AliAODMCParticle *dMC = (AliAODMCParticle*)arrayMC->At(matchtoMC);
	Int_t pdgMC = dMC->GetPdgCode();
	if (TMath::Abs(pdgMC)!=absPdgMom) AliInfo("What's up, isn't it a lambdac ?!");
	Int_t labDau0 = ((AliAODTrack*)d->GetDaughter(0))->GetLabel();
	AliAODMCParticle* p0 = (AliAODMCParticle*)arrayMC->UncheckedAt(labDau0);
	Int_t pdgCode0 = TMath::Abs(p0->GetPdgCode());
	Int_t labDau1 = ((AliAODTrack*)d->GetDaughter(1))->GetLabel();
	AliAODMCParticle* p1 = (AliAODMCParticle*)arrayMC->UncheckedAt(labDau1);
	Int_t pdgCode1 = TMath::Abs(p1->GetPdgCode());
	Int_t labDau2 = ((AliAODTrack*)d->GetDaughter(2))->GetLabel();
	AliAODMCParticle* p2 = (AliAODMCParticle*)arrayMC->UncheckedAt(labDau2);
	Int_t pdgCode2 = TMath::Abs(p2->GetPdgCode());

	// Fill in the histograms in case of p K pi decays
	if(ispKpi==1){
	  if(pdgCode0==pPDG && pdgCode1==kPDG && pdgCode2==piPDG){
	    fSigHist[index]->Fill(masses[0]);
	  } else {
	    fRflHist[index]->Fill(masses[0]);
	  }
	}
	// Fill in the histograms in case of pi K p decays
	if(ispiKp==2){
	  if(pdgCode0==piPDG && pdgCode1==kPDG && pdgCode2==pPDG){
	    fSigHist[index]->Fill(masses[1]);
	  } else {
	    fRflHist[index]->Fill(masses[1]);
	  }
	}
      } else {
	if(ispKpi==1) fBkgHist[index]->Fill(masses[0]);
	if(ispiKp==2) fBkgHist[index]->Fill(masses[1]);
      }
    }
  }


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
      if(fDecChannel != AliAnalysisTaskSESignificance::kDplustoKpipi) fRflHist[i]=dynamic_cast<TH1F*>(fOutput->FindObject(rflname.Data()));
    }
    
  }
  

 
  
  return;
}
//_________________________________________________________________________________________________
Int_t AliAnalysisTaskSESignificance::CheckOrigin(const AliAODMCParticle* mcPart, const TClonesArray* mcArray)const{

	//
	// checking whether the very mother of the D0 is a charm or a bottom quark
	//

	Int_t pdgGranma = 0;
	Int_t mother = 0;
	mother = mcPart->GetMother();
	Int_t istep = 0;
	while (mother >0 ){
		istep++;
		AliDebug(2,Form("mother at step %d = %d", istep, mother));
		AliAODMCParticle* mcGranma = dynamic_cast<AliAODMCParticle*>(mcArray->At(mother));
		if(!mcGranma) break;
		pdgGranma = mcGranma->GetPdgCode();
		AliDebug(2,Form("Pdg mother at step %d = %d", istep, pdgGranma));
		Int_t abspdgGranma = TMath::Abs(pdgGranma);
		if ((abspdgGranma > 500 && abspdgGranma < 600) || (abspdgGranma > 5000 && abspdgGranma < 6000)) {
			break;
		}
		mother = mcGranma->GetMother();
	}
	return pdgGranma;
}
//_________________________________________________________________________________________________
