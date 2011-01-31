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
  fBFeedDown(kBoth),
  fDecChannel(0),
  fSelectionlevel(0),
  fNBins(100),
  fPartOrAndAntiPart(0)
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
  fBFeedDown(kBoth),
  fDecChannel(decaychannel),
  fSelectionlevel(selectionlevel),
  fNBins(100),
  fPartOrAndAntiPart(0)
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

  SetMassLimits(0.15,pdg); //check range
  fNPtBins=fRDCuts->GetNPtBins();

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
  fHistNEvents->GetXaxis()->SetBinLabel(7,"MC Cand from c");
  fHistNEvents->GetXaxis()->SetBinLabel(8,"MC Cand from b");
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
      AliWarning("AliAnalysisTaskSESignificance::UserExec:MC particles branch not found!\n");
      //    return;
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
  Int_t nprongs=1;
  Int_t *pdgdaughters=0x0;
  Int_t absPdgMom=411;
  

  switch(fDecChannel){
  case 0:
    //D+
    pdgdaughters =new Int_t[3];
    pdgdaughters[0]=211;//pi
    pdgdaughters[1]=321;//K
    pdgdaughters[2]=211;//pi
    nprongs=3;
    absPdgMom=411;
    break;
  case 1:
    //D0
    pdgdaughters =new Int_t[2];
    pdgdaughters[0]=211;//pi 
    pdgdaughters[1]=321;//K
    nprongs=2;
    absPdgMom=421;
    break;
  case 2:
    //D*
    pdgdaughters =new Int_t[3];
    pdgdaughters[1]=211;//pi
    pdgdaughters[0]=321;//K
    pdgdaughters[2]=211;//pi (soft?)
    nprongs=3;
    absPdgMom=413;
    break;
  case 3:
    //Ds
    pdgdaughters =new Int_t[3];
    pdgdaughters[0]=321;//K
    pdgdaughters[1]=321;//K
    pdgdaughters[2]=211;//pi
    nprongs=3;
    absPdgMom=431;
    break;
  case 4:
    //D0 in 4 prongs
    pdgdaughters =new Int_t[4];
    pdgdaughters[0]=321;
    pdgdaughters[1]=211;
    pdgdaughters[2]=211;
    pdgdaughters[3]=211;
    nprongs=4;
    absPdgMom=421;
    break;
  case 5:
    //Lambda_c
    pdgdaughters =new Int_t[3];
    pdgdaughters[0]=2212;//p
    pdgdaughters[1]=321;//K
    pdgdaughters[2]=211;//pi
    nprongs=3;
    absPdgMom=4122;
    break;
  }

  Int_t nHistpermv=((AliMultiDimVector*)fCutList->FindObject("multiDimVectorPtBin0"))->GetNTotCells();
  Int_t nProng = arrayProng->GetEntriesFast();
  if(fDebug>1) printf("Number of D2H: %d\n",nProng);

  // trigger class for PbPb C0SMH-B-NOPF-ALLNOTRD, C0SMH-B-NOPF-ALL
  TString trigclass=aod->GetFiredTriggerClasses();
  if(trigclass.Contains("C0SMH-B-NOPF-ALLNOTRD") || trigclass.Contains("C0SMH-B-NOPF-ALL")) fHistNEvents->Fill(5);

  if(fRDCuts->IsEventSelected(aod)) fHistNEvents->Fill(1);
  else{
    if(fRDCuts->GetWhyRejection()==1) // rejected for pileup
      fHistNEvents->Fill(4);
    return;
  }
  
  for (Int_t iProng = 0; iProng < nProng; iProng++) {
    
    d=(AliAODRecoDecayHF*)arrayProng->UncheckedAt(iProng);
    
    Bool_t isFidAcc = fRDCuts->IsInFiducialAcceptance(d->Pt(),d->Y(absPdgMom));
    Int_t isSelected=fRDCuts->IsSelected(d,fSelectionlevel,aod);

    if(fReadMC && fBFeedDown!=kBoth && isSelected){
      Int_t labD = d->MatchToMC(absPdgMom,arrayMC,nprongs,pdgdaughters);
      if(labD>=0){
	AliAODMCParticle *partD = (AliAODMCParticle*)arrayMC->At(labD);
	Int_t label=partD->GetMother();
	AliAODMCParticle *mot = (AliAODMCParticle*)arrayMC->At(label);
	while(label>=0){//get first mother
	  mot = (AliAODMCParticle*)arrayMC->At(label);
	  label=mot->GetMother();
	}
	Int_t pdgMotCode = mot->GetPdgCode();
	
	if(TMath::Abs(pdgMotCode)<=4){
	  fHistNEvents->Fill(6);
	  if(fBFeedDown==kBeautyOnly)isSelected=kFALSE; //from primary charm
	}else{
	  fHistNEvents->Fill(7);
	  if(fBFeedDown==kCharmOnly) isSelected=kFALSE; //from beauty
	}
	
	/*
	if(TMath::Abs(pdgMotCode)==4 && fBFeedDown==kBeautyOnly) isSelected=kFALSE; //from primary charm
	if(TMath::Abs(pdgMotCode)==5 && fBFeedDown==kCharmOnly) isSelected=kFALSE; //from beauty
	*/
      }
    }
    

    if(isSelected&&isFidAcc) {
      fHistNEvents->Fill(2); // count selected with loosest cuts
      if(fDebug>1) printf("+++++++Is Selected\n");
    
      Double_t* invMass=0x0;
      Int_t nmasses;
      CalculateInvMasses(d,invMass,nmasses);

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
	  delete addresses;
	  return;
	}

	fHistNEvents->Fill(3);

	//fill the histograms with the appropriate method
	switch (fDecChannel){
	case 0:
	  FillDplus(d,arrayMC,(Int_t)(ptbin*nHistpermv+addresses[ivals]),invMass,isSelected);
	  break;
	case 1:
	  FillD02p(d,arrayMC,(Int_t)(ptbin*nHistpermv+addresses[ivals]),invMass,isSelected);
	  break;
	case 2:
	  FillDstar(d,arrayMC,(Int_t)(ptbin*nHistpermv+addresses[ivals]),invMass,isSelected);
	  break;
	case 3:
	  FillDs(d,arrayMC,(Int_t)(ptbin*nHistpermv+addresses[ivals]),invMass,isSelected);
	  break;
	case 4:
	  FillD04p(d,arrayMC,(Int_t)(ptbin*nHistpermv+addresses[ivals]),invMass,isSelected);
	  break;
	case 5:
	  FillLambdac(d,arrayMC,(Int_t)(ptbin*nHistpermv+addresses[ivals]),invMass,isSelected);
	  break;
	default:
	  break;
	}
	
      }
    }// end if selected
    
  }
  
  if(pdgdaughters) {delete [] pdgdaughters; pdgdaughters=NULL;}

  PostData(1,fOutput);    
  return;
}

//***************************************************************************

// Methods used in the UserExec

void AliAnalysisTaskSESignificance::CalculateInvMasses(AliAODRecoDecayHF* d,Double_t*& masses,Int_t& nmasses){
  //Calculates all the possible invariant masses for each candidate
  //NB: the implementation for each candidate is responsibility of the corresponding developer

  switch(fDecChannel){
  case 0:
    //D+ -- Giacomo, Renu ?
    {
      nmasses=1;
      masses=new Double_t[nmasses];
      Int_t pdgdaughters[3] = {211,321,211};
      masses[0]=d->InvMass(3,(UInt_t*)pdgdaughters);
    }
    break;
  case 1:
    //D0 (Kpi)  -- Chiara
    {
      const Int_t ndght=2;
      nmasses=2;
      masses=new Double_t[nmasses];
      Int_t pdgdaughtersD0[ndght]={211,321};//pi,K 
      masses[0]=d->InvMass(ndght,(UInt_t*)pdgdaughtersD0); //D0
      Int_t pdgdaughtersD0bar[ndght]={321,211};//K,pi 
      masses[1]=d->InvMass(ndght,(UInt_t*)pdgdaughtersD0bar); //D0bar
    }
    break;
  case 2:
    //D* -- ? Yifei, Alessandro
    {
      //.....
    }
    break;
  case 3:
    //Ds  -- Sergey, Sahdana
    {
      const Int_t ndght=3;
      nmasses=2;
      masses=new Double_t[nmasses];
      Int_t pdgdaughtersDsKKpi[ndght]={321,321,211};//K,K,pi 
      masses[0]=d->InvMass(ndght,(UInt_t*)pdgdaughtersDsKKpi); //Ds
      Int_t pdgdaughtersDspiKK[ndght]={211,321,321};//pi,K,K 
      masses[1]=d->InvMass(ndght,(UInt_t*)pdgdaughtersDspiKK); //Dsbar 

      //.....
    }
    break;
  case 4:
    //D0 (Kpipipi) -- ? Rossella, Fabio ???
    {
      //.....
    }
    break;
  case 5:
    //Lambda_c -- Rossella
    {
      //.....
    }
    break;
  default:
    break;
  }
}

//********************************************************************************************

//Methods to fill the istograms with MC information, one for each candidate
//NB: the implementation for each candidate is responsibility of the corresponding developer

void AliAnalysisTaskSESignificance::FillDplus(AliAODRecoDecayHF* d,TClonesArray *arrayMC,Int_t index,Double_t* masses,Int_t isSel){
    //D+ channel
  if(!isSel){
    AliError("Candidate not selected\n");
    return;
  }

  if(fPartOrAndAntiPart*d->GetCharge()<0)return;

  fMassHist[index]->Fill(masses[0]);

  Int_t pdgdaughters[3] = {211,321,211};

  if(fReadMC){
    Int_t lab=-1;
    lab = d->MatchToMC(411,arrayMC,3,pdgdaughters);
    if(lab>=0){ //signal
      fSigHist[index]->Fill(masses[0]);
    } else{ //background
      fBkgHist[index]->Fill(masses[0]);
    } 
  }   
}


void AliAnalysisTaskSESignificance::FillD02p(AliAODRecoDecayHF* d,TClonesArray *arrayMC,Int_t index,Double_t* masses,Int_t isSel){

  //D0->Kpi channel

  //mass histograms
  if(!masses){
    AliError("Masses not calculated\n");
    return;
  }
  if((isSel==1 || isSel==3) && fPartOrAndAntiPart>=0) fMassHist[index]->Fill(masses[0]);
  if(isSel>=2 && fPartOrAndAntiPart<=0) fMassHist[index]->Fill(masses[1]);



  //MC histograms
  if(fReadMC){

    Int_t matchtoMC=-1;

    //D0
    Int_t pdgdaughters[2];
    pdgdaughters[0]=211;//pi 
    pdgdaughters[1]=321;//K
    Int_t nprongs=2;
    Int_t absPdgMom=421;

    matchtoMC = d->MatchToMC(absPdgMom,arrayMC,nprongs,pdgdaughters);

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

void AliAnalysisTaskSESignificance::FillDstar(AliAODRecoDecayHF* /*d*/,TClonesArray */*arrayMC*/,Int_t /*index*/,Double_t* /*masses*/,Int_t /*matchtoMC*/){
    //D* channel

  AliInfo("Dstar channel not implemented\n");

}


void AliAnalysisTaskSESignificance::FillDs(AliAODRecoDecayHF* d,TClonesArray *arrayMC,Int_t index,Double_t* masses,Int_t isSel){

  //AliInfo("Ds channel not implemented\n");


  if(!masses){
    AliError("Masses not calculated\n");
    return;
  }
  
  Int_t pdgDstoKKpi[3]={321,321,211};
  Int_t labDs=-1;
  if(fReadMC){
    labDs = d->MatchToMC(431,arrayMC,3,pdgDstoKKpi);
  }
  
  
  if(isSel&1  && fPartOrAndAntiPart*d->GetCharge()>=0) {
    
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
  
  if(isSel&2 && fPartOrAndAntiPart*d->GetCharge()>=0){
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
  
     
  //MC histograms
 
  //Ds channel
	// Int_t isKKpi=0;
	// Int_t ispiKK=0;
        //   isKKpi=isSelected&1;
	//   ispiKK=isSelected&2;

	// Int_t labDau0=((AliAODTrack*)d->GetDaughter(0))->GetLabel();
	// AliAODMCParticle* p=(AliAODMCParticle*)arrayMC->UncheckedAt(labDau0);
	// Int_t pdgCode0=TMath::Abs(p->GetPdgCode());

	//       if(isKKpi){
	//        if(pdgCode0==321){
        //         fSigHist[index]->Fill(masses[0]);
	//        }else{
	//         fRflHist[index]->Fill(masses[0]);
	//        }
	//       }
	//       if(ispiKK){
	//        if(pdgCode0==211){
        //         fSigHist[index]->Fill(masses[1]);
	//        }else{
	//         fRflHist[index]->Fill(masses[1]);
	//        }
	//       }

}

void AliAnalysisTaskSESignificance::FillD04p(AliAODRecoDecayHF* /*d*/,TClonesArray */*arrayMC*/,Int_t /*index*/,Double_t* /*masses*/,Int_t /*matchtoMC*/){
  //D0->Kpipipi channel
  AliInfo("D0 in 4 prongs channel not implemented\n");

}

void AliAnalysisTaskSESignificance::FillLambdac(AliAODRecoDecayHF* /*d*/,TClonesArray */*arrayMC*/,Int_t /*index*/,Double_t* /*masses*/,Int_t /*matchtoMC*/){
  AliInfo("Lambdac channel not implemented\n");

  //Lambdac channel
  // Int_t ispKpi=0;
  // Int_t ispiKp=0;

  // ispKpi=isSelected&1;
  // ispiKp=isSelected&2;
  // if(matchtoMC>=0){	
  // Int_t labDau0=((AliAODTrack*)d->GetDaughter(0))->GetLabel();
  // AliAODMCParticle* p=(AliAODMCParticle*)arrayMC->UncheckedAt(labDau0);
  // Int_t pdgCode0=TMath::Abs(p->GetPdgCode());
  // if(ispKpi){
  //   if(pdgCode0==2212){
  //     fSigHist[index]->Fill(invMass[0]);
  //   }else{
  //     fRflHist[index]->Fill(invMass[0]);
  //   }
  // }
  // if(ispiKp){
  //   if(pdgCode0==211){
  //     fSigHist[index]->Fill(invMass[1]);
  //   }else{
  //     fRflHist[index]->Fill(invMass[1]);
  //   }
  // }
  // }else{
  //   fBkgHist[index]
  // }


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
//-------------------------------------------

