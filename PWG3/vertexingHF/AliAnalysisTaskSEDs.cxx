/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
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

///////////////////////////////////////////////////////////////////
//                                                               //
//  Analysis task to produce Ds candidates mass spectra          //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////

#include <TClonesArray.h>
#include <TNtuple.h>
#include <TList.h>
#include <TString.h>
#include <TH1F.h>
#include <TMath.h>
#include <TDatabasePDG.h>

#include "AliAnalysisManager.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliAnalysisVertexingHF.h"
#include "AliRDHFCutsDstoKKpi.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskSEDs.h"
#include "AliNormalizationCounter.h"

ClassImp(AliAnalysisTaskSEDs)


//________________________________________________________________________
AliAnalysisTaskSEDs::AliAnalysisTaskSEDs():
  AliAnalysisTaskSE(),
  fOutput(0), 
  fHistNEvents(0),
  fPtVsMass(0),
  fYVsPt(0),
  fYVsPtSig(0),
  fNtupleDs(0),
  fFillNtuple(0),
  fReadMC(kFALSE),
  fDoCutVarHistos(kTRUE),
  fUseSelectionBit(kFALSE),
  fNPtBins(0),
  fListCuts(0),
  fMassRange(0.8),
  fMassBinSize(0.002),
  fCounter(0),
  fProdCuts(0),
  fAnalysisCuts(0)
{
  // Default constructor
}

//________________________________________________________________________
AliAnalysisTaskSEDs::AliAnalysisTaskSEDs(const char *name, AliRDHFCutsDstoKKpi* productioncuts, AliRDHFCutsDstoKKpi* analysiscuts,Int_t fillNtuple):
  AliAnalysisTaskSE(name),
  fOutput(0),
  fHistNEvents(0),
  fPtVsMass(0),
  fYVsPt(0),
  fYVsPtSig(0),
  fNtupleDs(0),
  fFillNtuple(fillNtuple),
  fReadMC(kFALSE),
  fDoCutVarHistos(kTRUE),
  fUseSelectionBit(kFALSE),
  fNPtBins(0),
  fListCuts(0),
  fMassRange(0.8),
  fMassBinSize(0.002),
  fCounter(0),
  fProdCuts(productioncuts),
  fAnalysisCuts(analysiscuts)
{
  // Default constructor
  // Output slot #1 writes into a TList container
  Int_t nptbins=fAnalysisCuts->GetNPtBins();
  Float_t *ptlim=fAnalysisCuts->GetPtBinLimits();
  SetPtBins(nptbins,ptlim);
 
  DefineOutput(1,TList::Class());  //My private output

  DefineOutput(2,TList::Class());

  DefineOutput(3,AliNormalizationCounter::Class());
  
  if(fFillNtuple>0){
    // Output slot #4 writes into a TNtuple container
    DefineOutput(4,TNtuple::Class());  //My private output
  }
  
}

//________________________________________________________________________
void AliAnalysisTaskSEDs::SetPtBins(Int_t n, Float_t* lim){
  // define pt bins for analysis
  if(n>kMaxPtBins){
    printf("Max. number of Pt bins = %d\n",kMaxPtBins);
    fNPtBins=kMaxPtBins;
    fPtLimits[0]=0.;
    fPtLimits[1]=1.;
    fPtLimits[2]=3.;
    fPtLimits[3]=5.;
    fPtLimits[4]=10.;
    for(Int_t i=5; i<kMaxPtBins+1; i++) fPtLimits[i]=99999999.;
  }else{
    fNPtBins=n;
    for(Int_t i=0; i<fNPtBins+1; i++) fPtLimits[i]=lim[i];
    for(Int_t i=fNPtBins+1; i<kMaxPtBins+1; i++) fPtLimits[i]=99999999.;
  }
  if(fDebug > 1){
    printf("Number of Pt bins = %d\n",fNPtBins);
    for(Int_t i=0; i<fNPtBins; i++) printf(" Bin%d = %8.2f-%8.2f\n",i,fPtLimits[i],fPtLimits[i+1]);    
  }
}
//________________________________________________________________________
AliAnalysisTaskSEDs::~AliAnalysisTaskSEDs()
{
  // Destructor
  if (fOutput) {
    delete fOutput;
    fOutput = 0;
  }

  if(fHistNEvents){
    delete fHistNEvents;
    fHistNEvents=0;
  } 
 
  if (fListCuts) {
    delete fListCuts;
    fListCuts = 0;
  }

  for(Int_t i=0;i<4*fNPtBins;i++){
    
    if(fMassHist[i]){ delete fMassHist[i]; fMassHist[i]=0;}
    if(fMassHistPhi[i]){ delete fMassHistPhi[i]; fMassHistPhi[i]=0;}
    if(fMassHistK0st[i]){ delete fMassHistK0st[i]; fMassHistK0st[i]=0;}
    if(fCosPHist[i]){ delete fCosPHist[i]; fCosPHist[i]=0;}
    if(fDLenHist[i]){ delete fDLenHist[i]; fDLenHist[i]=0;}
    if(fSumd02Hist[i]){ delete fSumd02Hist[i]; fSumd02Hist[i]=0;}
    if(fSigVertHist[i]){ delete fSigVertHist[i]; fSigVertHist[i]=0;}
    if(fPtMaxHist[i]){ delete fPtMaxHist[i]; fPtMaxHist[i]=0;}
    if(fDCAHist[i]){ delete fDCAHist[i]; fDCAHist[i]=0;}
    if(fPtProng0Hist[i]){ delete fPtProng0Hist[i]; fPtProng0Hist[i]=0;}
    if(fPtProng1Hist[i]){ delete fPtProng1Hist[i]; fPtProng1Hist[i]=0;}
    if(fPtProng2Hist[i]){ delete fPtProng2Hist[i]; fPtProng2Hist[i]=0;}
    if(fDalitz[i]){ delete fDalitz[i]; fDalitz[i]=0;}
    if(fDalitzPhi[i]){ delete fDalitzPhi[i]; fDalitzPhi[i]=0;}
    if(fDalitzK0st[i]){ delete fDalitzK0st[i]; fDalitzK0st[i]=0;}

  }

  if(fPtVsMass){
    delete fPtVsMass;
    fPtVsMass=0;
  }
  if(fYVsPt){
    delete fYVsPt;
    fYVsPt=0;
  }
  if(fYVsPtSig){
    delete fYVsPtSig;
    fYVsPtSig=0;
  }
  if(fNtupleDs){
    delete fNtupleDs;
    fNtupleDs=0;
  }
  if(fCounter){
    delete fCounter;
    fCounter = 0;
  }
  
  if (fProdCuts) {
    delete fProdCuts;
    fProdCuts = 0;
  }
  if (fAnalysisCuts) {
    delete fAnalysisCuts;
    fAnalysisCuts = 0;
  }
}  

//________________________________________________________________________
void AliAnalysisTaskSEDs::Init()
{
  // Initialization

  if(fDebug > 1) printf("AnalysisTaskSEDs::Init() \n");

  fListCuts=new TList();
  fListCuts->SetOwner();
  fListCuts->SetName("CutObjects");

  AliRDHFCutsDstoKKpi *production = new AliRDHFCutsDstoKKpi(*fProdCuts);
  production->SetName("ProductionCuts");
  AliRDHFCutsDstoKKpi *analysis = new AliRDHFCutsDstoKKpi(*fAnalysisCuts);
  analysis->SetName("AnalysisCuts");
  
  fListCuts->Add(production);
  fListCuts->Add(analysis);
  PostData(2,fListCuts);
  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEDs::UserCreateOutputObjects()
{
  // Create the output container
  //
  if(fDebug > 1) printf("AnalysisTaskSEDs::UserCreateOutputObjects() \n");

  // Several histograms are more conveniently managed in a TList
  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("OutputHistos");

  Double_t massDs=TDatabasePDG::Instance()->GetParticle(431)->Mass();
  
  Int_t nInvMassBins=(Int_t)(fMassRange/fMassBinSize+0.5);
  if(nInvMassBins%2==1) nInvMassBins++;
  // Double_t minMass=massDs-1.0*nInvMassBins*fMassBinSize;
  Double_t minMass=massDs-0.5*nInvMassBins*fMassBinSize;
  //  Double_t maxMass=massDs+1.0*nInvMassBins*fMassBinSize;
  Double_t maxMass=massDs+0.5*nInvMassBins*fMassBinSize;

  TString hisname;
  TString htype;
  Int_t index;
  for(Int_t iType=0; iType<4; iType++){
    for(Int_t i=0;i<fNPtBins;i++){
      if(iType==0){
	htype="All";
	index=GetHistoIndex(i);
      }else if(iType==1){ 
	htype="Sig";
	index=GetSignalHistoIndex(i);
      }else if(iType==2){ 
	htype="Bkg";
	index=GetBackgroundHistoIndex(i);
      }else{ 
	htype="ReflSig";
	index=GetReflSignalHistoIndex(i);
      }
      hisname.Form("hMass%sPt%d",htype.Data(),i);
      fMassHist[index]=new TH1F(hisname.Data(),hisname.Data(),nInvMassBins,minMass,maxMass);
      fMassHist[index]->Sumw2();
      hisname.Form("hMass%sPt%dphi",htype.Data(),i);
      fMassHistPhi[index]=new TH1F(hisname.Data(),hisname.Data(),nInvMassBins,minMass,maxMass);
      fMassHistPhi[index]->Sumw2();
      hisname.Form("hMass%sPt%dk0st",htype.Data(),i);
      fMassHistK0st[index]=new TH1F(hisname.Data(),hisname.Data(),nInvMassBins,minMass,maxMass);
      fMassHistK0st[index]->Sumw2();
      hisname.Form("hCosP%sPt%d",htype.Data(),i);
      fCosPHist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.5,1.);
      fCosPHist[index]->Sumw2();
      hisname.Form("hDLen%sPt%d",htype.Data(),i);
      fDLenHist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.,0.5);
      fDLenHist[index]->Sumw2();
      hisname.Form("hSumd02%sPt%d",htype.Data(),i);
      fSumd02Hist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.,1.);
      fSumd02Hist[index]->Sumw2();
      hisname.Form("hSigVert%sPt%d",htype.Data(),i);
      fSigVertHist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.,0.1);
      fSigVertHist[index]->Sumw2();
      hisname.Form("hPtMax%sPt%d",htype.Data(),i);
      fPtMaxHist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.5,20.);
      fPtMaxHist[index]->Sumw2();
      hisname.Form("hPtCand%sPt%d",htype.Data(),i);
      fPtCandHist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.5,20.);
      fPtCandHist[index]->Sumw2();
      hisname.Form("hDCA%sPt%d",htype.Data(),i);
      fDCAHist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.,0.1);
      fDCAHist[index]->Sumw2();
      hisname.Form("hPtProng0%sPt%d",htype.Data(),i);
      fPtProng0Hist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.0,20.);
      fPtProng0Hist[index]->Sumw2();
      hisname.Form("hPtProng1%sPt%d",htype.Data(),i);
      fPtProng1Hist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.0,20.);
      fPtProng1Hist[index]->Sumw2();
      hisname.Form("hPtProng2%sPt%d",htype.Data(),i);
      fPtProng2Hist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.0,20.);
      fPtProng2Hist[index]->Sumw2();
      hisname.Form("hDalitz%sPt%d",htype.Data(),i);
      fDalitz[index]=new TH2F(hisname.Data(),hisname.Data(),100,0.,2.,100,0.,2.);
      fDalitz[index]->Sumw2();
      hisname.Form("hDalitz%sPt%dphi",htype.Data(),i);
      fDalitzPhi[index]=new TH2F(hisname.Data(),hisname.Data(),100,0.,2.,100,0.,2.);
      fDalitzPhi[index]->Sumw2();
      hisname.Form("hDalitz%sPt%dk0st",htype.Data(),i);
      fDalitzK0st[index]=new TH2F(hisname.Data(),hisname.Data(),100,0.,2.,100,0.,2.);
      fDalitzK0st[index]->Sumw2();
    }
  }

  for(Int_t i=0; i<4*fNPtBins; i++){
    fOutput->Add(fMassHist[i]);
    fOutput->Add(fMassHistPhi[i]);
    fOutput->Add(fMassHistK0st[i]);
    fOutput->Add(fCosPHist[i]);
    fOutput->Add(fDLenHist[i]);
    fOutput->Add(fSumd02Hist[i]);
    fOutput->Add(fSigVertHist[i]);
    fOutput->Add(fPtMaxHist[i]);
    fOutput->Add(fPtCandHist[i]);
    fOutput->Add(fDCAHist[i]);
    fOutput->Add(fPtProng0Hist[i]);
    fOutput->Add(fPtProng1Hist[i]);
    fOutput->Add(fPtProng2Hist[i]);
    fOutput->Add(fDalitz[i]);
    fOutput->Add(fDalitzPhi[i]);
    fOutput->Add(fDalitzK0st[i]);
  }

  fChanHist[0] = new TH1F("hChanAll", "KKpi and piKK candidates",64,-0.5,63.5);
  fChanHist[1] = new TH1F("hChanSig", "KKpi and piKK candidates",64,-0.5,63.5);
  fChanHist[2] = new TH1F("hChanBkg", "KKpi and piKK candidates",64,-0.5,63.5);
  fChanHist[3] = new TH1F("hChanReflSig", "KKpi and piKK candidates",64,-0.5,63.5);
  for(Int_t i=0;i<4;i++){
    fChanHist[i]->Sumw2();
    fChanHist[i]->SetMinimum(0);
    fOutput->Add(fChanHist[i]);
  }

  fHistNEvents = new TH1F("hNEvents", "number of events ",12,-0.5,11.5);
  fHistNEvents->GetXaxis()->SetBinLabel(1,"nEventsAnal");
  fHistNEvents->GetXaxis()->SetBinLabel(2,"n. passing IsEvSelected");
  fHistNEvents->GetXaxis()->SetBinLabel(3,"n. rejected due to trigger");
  fHistNEvents->GetXaxis()->SetBinLabel(4,"n. rejected due to not reco vertex");
  fHistNEvents->GetXaxis()->SetBinLabel(5,"n. rejected for contr vertex");
  fHistNEvents->GetXaxis()->SetBinLabel(6,"n. rejected for vertex out of accept");
  fHistNEvents->GetXaxis()->SetBinLabel(7,"n. rejected for pileup events");
  fHistNEvents->GetXaxis()->SetBinLabel(8,"no. of out centrality events");
  fHistNEvents->GetXaxis()->SetBinLabel(9,"no. of candidate");
  fHistNEvents->GetXaxis()->SetBinLabel(10,"no. of Ds after loose cuts");
  fHistNEvents->GetXaxis()->SetBinLabel(11,"no. of Ds after tight cuts");
  fHistNEvents->GetXaxis()->SetBinLabel(12,"no. of cand wo bitmask");

  fHistNEvents->GetXaxis()->SetNdivisions(1,kFALSE);

  fHistNEvents->Sumw2();
  fHistNEvents->SetMinimum(0);
  fOutput->Add(fHistNEvents);

  fPtVsMass=new TH2F("hPtVsMass","PtVsMass (prod. cuts)",nInvMassBins,minMass,maxMass,40,0.,20.);
  fYVsPt=new TH2F("hYVsPt","YvsPt (prod. cuts)",40,0.,20.,80,-2.,2.);
  fYVsPtSig=new TH2F("hYVsPtSig","YvsPt (MC, only sig., prod. cuts)",40,0.,20.,80,-2.,2.);

  fOutput->Add(fPtVsMass);
  fOutput->Add(fYVsPt);
  fOutput->Add(fYVsPtSig);

  //Counter for Normalization
  fCounter = new AliNormalizationCounter("NormalizationCounter");
  fCounter->Init();

  PostData(1,fOutput); 
  PostData(3,fCounter);   
  
  if(fFillNtuple>0){
    OpenFile(4); // 4 is the slot number of the ntuple
    
    fNtupleDs = new TNtuple("fNtupleDs","Ds","labDs:retcode:pdgcode0:Pt0:Pt1:Pt2:PtRec:PointingAngle:DecLeng:VxRec:VyRec:VzRec:InvMassKKpi:InvMasspiKK:sigvert:d00:d01:d02:dca:d0square:InvMassPhiKKpi:InvMassPhipiKK:cosinePiDsFrameKKpi:cosinePiDsFramepiKK:cosineKPhiFrameKKpi:cosineKPhiFramepiKK"); 
    
  }
  

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEDs::UserExec(Option_t */*option*/)
{
  // Ds selection for current event, fill mass histos and selecetion variable histo
  // separate signal and backgound if fReadMC is activated

  AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());

  TClonesArray *array3Prong = 0;
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
      array3Prong=(TClonesArray*)aodFromExt->GetList()->FindObject("Charm3Prong");
    }
  } else if(aod) {
    array3Prong=(TClonesArray*)aod->GetList()->FindObject("Charm3Prong");
  }

  if(!aod || !array3Prong) {
    printf("AliAnalysisTaskSEDs::UserExec: Charm3Prong branch not found!\n");
    return;
  }


  // fix for temporary bug in ESDfilter 
  // the AODs with null vertex pointer didn't pass the PhysSel
  if(!aod->GetPrimaryVertex() || TMath::Abs(aod->GetMagneticField())<0.001) return;


  fHistNEvents->Fill(0); // count event
  // Post the data already here
  PostData(1,fOutput);
  
  fCounter->StoreEvent(aod,fProdCuts,fReadMC);
  

  Bool_t isEvSel=fAnalysisCuts->IsEventSelected(aod);
  if(fAnalysisCuts->IsEventRejectedDueToTrigger())fHistNEvents->Fill(2);
  if(fAnalysisCuts->IsEventRejectedDueToNotRecoVertex())fHistNEvents->Fill(3);
  if(fAnalysisCuts->IsEventRejectedDueToVertexContributors())fHistNEvents->Fill(4);
  if(fAnalysisCuts->IsEventRejectedDueToZVertexOutsideFiducialRegion())fHistNEvents->Fill(5);
  if(fAnalysisCuts->IsEventRejectedDueToPileupSPD())fHistNEvents->Fill(6);
  if(fAnalysisCuts->IsEventRejectedDueToCentrality())fHistNEvents->Fill(7);
  
  
  
  if(!isEvSel)return;
  
  fHistNEvents->Fill(1);

  TClonesArray *arrayMC=0;
  AliAODMCHeader *mcHeader=0;

  // AOD primary vertex
  AliAODVertex *vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();
  //    vtx1->Print();
  
  // load MC particles
  if(fReadMC){
    
    arrayMC =  (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
    if(!arrayMC) {
      printf("AliAnalysisTaskSEDs::UserExec: MC particles branch not found!\n");
      return;
    }
    
  // load MC header
    mcHeader =  (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if(!mcHeader) {
      printf("AliAnalysisTaskSEDs::UserExec: MC header branch not found!\n");
      return;
    }
  }
  
  Int_t n3Prong = array3Prong->GetEntriesFast();
  if(fDebug>1) printf("Number of Ds->KKpi: %d\n",n3Prong);
  
  
  Int_t pdgDstoKKpi[3]={321,321,211};
  Int_t nSelectedloose=0;
  Int_t nSelectedtight=0;

  for (Int_t i3Prong = 0; i3Prong < n3Prong; i3Prong++) {
  
    AliAODRecoDecayHF3Prong *d = (AliAODRecoDecayHF3Prong*)array3Prong->UncheckedAt(i3Prong);
    fHistNEvents->Fill(8);
    
    if(fUseSelectionBit && !(d->HasSelectionBit(AliRDHFCuts::kDsCuts))){
      fHistNEvents->Fill(11);
      continue;
    }
    
    Bool_t unsetvtx=kFALSE;
    if(!d->GetOwnPrimaryVtx()){
      d->SetOwnPrimaryVtx(vtx1);
      unsetvtx=kTRUE;
    }
    
    Bool_t recVtx=kFALSE;
    AliAODVertex *origownvtx=0x0;
    Int_t retCodeProdCuts=fProdCuts->IsSelected(d,AliRDHFCuts::kCandidate,aod);
   
    if(retCodeProdCuts) {
      if(fProdCuts->GetIsPrimaryWithoutDaughters()){
   	    if(d->GetOwnPrimaryVtx()) origownvtx=new AliAODVertex(*d->GetOwnPrimaryVtx());	
   	    if(fProdCuts->RecalcOwnPrimaryVtx(d,aod))recVtx=kTRUE;
   	    else fProdCuts->CleanOwnPrimaryVtx(d,aod,origownvtx);
      }
    }  
    
    Double_t ptCand = d->Pt();
    Int_t iPtBin=TMath::BinarySearch(fNPtBins,fPtLimits,(Float_t)ptCand);
    Int_t retCodeAnalysisCuts=fAnalysisCuts->IsSelected(d,AliRDHFCuts::kCandidate,aod);
    Double_t rapid=d->YDs(); 
    fYVsPt->Fill(ptCand,rapid);

    Bool_t isFidAcc=fAnalysisCuts->IsInFiducialAcceptance(ptCand,rapid);
    
    if(retCodeProdCuts>0){
      if(isFidAcc){
        nSelectedloose++;
        fHistNEvents->Fill(9);
        if(retCodeAnalysisCuts>0)nSelectedtight++;
      }
    }
  
    if(retCodeAnalysisCuts<=0) continue;
    if(!isFidAcc) continue;
    fHistNEvents->Fill(10);
    
    Int_t index=GetHistoIndex(iPtBin);
    fPtCandHist[index]->Fill(ptCand);

    Int_t isKKpi=retCodeAnalysisCuts&1;
    Int_t ispiKK=retCodeAnalysisCuts&2;
    Int_t isPhiKKpi=retCodeAnalysisCuts&4;
    Int_t isPhipiKK=retCodeAnalysisCuts&8;
    Int_t isK0starKKpi=retCodeAnalysisCuts&16;    
    Int_t isK0starpiKK=retCodeAnalysisCuts&32;

    fChanHist[0]->Fill(retCodeAnalysisCuts);
 
    Int_t indexMCKKpi=-1;
    Int_t indexMCpiKK=-1;
    Int_t labDs=-1;
    Int_t pdgCode0=-999;
    
    if(fReadMC){
      labDs = d->MatchToMC(431,arrayMC,3,pdgDstoKKpi);
      if(labDs>=0){
	Int_t labDau0=((AliAODTrack*)d->GetDaughter(0))->GetLabel();
	AliAODMCParticle* p=(AliAODMCParticle*)arrayMC->UncheckedAt(TMath::Abs(labDau0));
	pdgCode0=TMath::Abs(p->GetPdgCode());

	if(isKKpi){
	  if(pdgCode0==321) {	  
	    indexMCKKpi=GetSignalHistoIndex(iPtBin);
	    fYVsPtSig->Fill(ptCand,rapid);
	    fChanHist[1]->Fill(retCodeAnalysisCuts);
	  }else{
	    indexMCKKpi=GetReflSignalHistoIndex(iPtBin);
	    fChanHist[3]->Fill(retCodeAnalysisCuts);
	  }
	}
	if(ispiKK){
	  if(pdgCode0==211) {	  
	    indexMCpiKK=GetSignalHistoIndex(iPtBin);
	    fYVsPtSig->Fill(ptCand,rapid);
	    fChanHist[1]->Fill(retCodeAnalysisCuts);
	  }else{
	    indexMCpiKK=GetReflSignalHistoIndex(iPtBin);
	    fChanHist[3]->Fill(retCodeAnalysisCuts);
	  }
	}
      }else{
	indexMCpiKK=GetBackgroundHistoIndex(iPtBin);
	indexMCKKpi=GetBackgroundHistoIndex(iPtBin);
	fChanHist[2]->Fill(retCodeAnalysisCuts);
      }
    }

    if(isKKpi){
      Double_t invMass=d->InvMassDsKKpi();
      fMassHist[index]->Fill(invMass);
      fPtVsMass->Fill(invMass,ptCand);
      if(isPhiKKpi) fMassHistPhi[index]->Fill(invMass); 
      if(isK0starKKpi) fMassHistK0st[index]->Fill(invMass);
      if(fReadMC  && indexMCKKpi!=-1){
	fMassHist[indexMCKKpi]->Fill(invMass);
	if(isPhiKKpi) fMassHistPhi[indexMCKKpi]->Fill(invMass);
	if(isK0starKKpi) fMassHistK0st[indexMCKKpi]->Fill(invMass);	  
      }
    }
    if(ispiKK){
      Double_t invMass=d->InvMassDspiKK();
      fMassHist[index]->Fill(invMass);
      fPtVsMass->Fill(invMass,ptCand);
      if(isPhipiKK) fMassHistPhi[index]->Fill(invMass);
      if(isK0starpiKK) fMassHistK0st[index]->Fill(invMass);
      if(fReadMC  && indexMCpiKK!=-1){
	fMassHist[indexMCpiKK]->Fill(invMass);
	if(isPhipiKK) fMassHistPhi[indexMCpiKK]->Fill(invMass);
	if(isK0starpiKK) fMassHistK0st[indexMCpiKK]->Fill(invMass);      
      }
    }

    if(fDoCutVarHistos){
      Double_t dlen=d->DecayLength();
      Double_t cosp=d->CosPointingAngle();
      Double_t pt0=d->PtProng(0);
      Double_t pt1=d->PtProng(1);
      Double_t pt2=d->PtProng(2);
      Double_t sigvert=d->GetSigmaVert();
      Double_t sumD02=d->Getd0Prong(0)*d->Getd0Prong(0)+d->Getd0Prong(1)*d->Getd0Prong(1)+d->Getd0Prong(2)*d->Getd0Prong(2);
      Double_t dca=d->GetDCA();      
      Double_t ptmax=0;
      for(Int_t i=0;i<3;i++){
	if(d->PtProng(i)>ptmax)ptmax=d->PtProng(i);
      }
      fCosPHist[index]->Fill(cosp);
      fDLenHist[index]->Fill(dlen);
      fSigVertHist[index]->Fill(sigvert);
      fSumd02Hist[index]->Fill(sumD02);
      fPtMaxHist[index]->Fill(ptmax);
      fDCAHist[index]->Fill(dca);
      fPtProng0Hist[index]->Fill(pt0);
      fPtProng1Hist[index]->Fill(pt1);
      fPtProng2Hist[index]->Fill(pt2);
      if(isKKpi){
	Double_t massKK=d->InvMass2Prongs(0,1,321,321);
	Double_t massKp=d->InvMass2Prongs(1,2,321,211);
	fDalitz[index]->Fill(massKK,massKp);
	if(isPhiKKpi) fDalitzPhi[index]->Fill(massKK,massKp);
	if(isK0starKKpi) fDalitzK0st[index]->Fill(massKK,massKp);
	if(fReadMC && indexMCKKpi!=-1){
	  fDalitz[indexMCKKpi]->Fill(massKK,massKp);
	  if(isPhiKKpi) fDalitzPhi[indexMCKKpi]->Fill(massKK,massKp);
	  if(isK0starKKpi) fDalitzK0st[indexMCKKpi]->Fill(massKK,massKp);
	  fCosPHist[indexMCKKpi]->Fill(cosp);
	  fDLenHist[indexMCKKpi]->Fill(dlen);
	  fSigVertHist[indexMCKKpi]->Fill(sigvert);
	  fSumd02Hist[indexMCKKpi]->Fill(sumD02);
	  fPtMaxHist[indexMCKKpi]->Fill(ptmax);
	  fPtCandHist[indexMCKKpi]->Fill(ptCand);
	  fDCAHist[indexMCKKpi]->Fill(dca);
	  fPtProng0Hist[indexMCKKpi]->Fill(pt0);
	  fPtProng1Hist[indexMCKKpi]->Fill(pt1);
	  fPtProng2Hist[indexMCKKpi]->Fill(pt2);
	}	
      }
      if(ispiKK){
	Double_t massKK=d->InvMass2Prongs(1,2,321,321);
	Double_t massKp=d->InvMass2Prongs(0,1,211,321);
	fDalitz[index]->Fill(massKK,massKp);
	if(isPhipiKK) fDalitzPhi[index]->Fill(massKK,massKp);
	if(isK0starpiKK) fDalitzK0st[index]->Fill(massKK,massKp);

	
	if(fReadMC && indexMCpiKK!=-1){
	  fDalitz[indexMCpiKK]->Fill(massKK,massKp);
	  if(isPhipiKK) fDalitzPhi[indexMCpiKK]->Fill(massKK,massKp);
	  if(isK0starpiKK) fDalitzK0st[indexMCpiKK]->Fill(massKK,massKp);
	  fCosPHist[indexMCpiKK]->Fill(cosp);
	  fDLenHist[indexMCpiKK]->Fill(dlen);
	  fSigVertHist[indexMCpiKK]->Fill(sigvert);
	  fSumd02Hist[indexMCpiKK]->Fill(sumD02);
	  fPtMaxHist[indexMCpiKK]->Fill(ptmax);
	  fPtCandHist[indexMCpiKK]->Fill(ptCand);
	  fDCAHist[indexMCpiKK]->Fill(dca);
	  fPtProng0Hist[indexMCpiKK]->Fill(pt0);
	  fPtProng1Hist[indexMCpiKK]->Fill(pt1);
	  fPtProng2Hist[indexMCpiKK]->Fill(pt2);
	}
      }
     
   
    }
    
    Float_t tmp[26];
    if(fFillNtuple>0){
      
      if ((fFillNtuple==1 && (isPhiKKpi || isPhipiKK)) || (fFillNtuple==2 && (isK0starKKpi || isK0starpiKK)) || (fFillNtuple==3 && (isKKpi || ispiKK))){
	
	tmp[0]=Float_t(labDs);
	tmp[1]=Float_t(retCodeAnalysisCuts);
	tmp[2]=Float_t(pdgCode0);  
	tmp[3]=d->PtProng(0);
	tmp[4]=d->PtProng(1);
	tmp[5]=d->PtProng(2);
	tmp[6]=d->Pt();
	tmp[7]=d->CosPointingAngle();
	tmp[8]=d->DecayLength();
	tmp[9]=d->Xv();
	tmp[10]=d->Yv();
	tmp[11]=d->Zv();
	tmp[12]=d->InvMassDsKKpi();
	tmp[13]=d->InvMassDspiKK();
	tmp[14]=d->GetSigmaVert();
	tmp[15]=d->Getd0Prong(0);
	tmp[16]=d->Getd0Prong(1);
	tmp[17]=d->Getd0Prong(2);
	tmp[18]=d->GetDCA();
	tmp[19]=d->Getd0Prong(0)*d->Getd0Prong(0)+d->Getd0Prong(1)*d->Getd0Prong(1)+d->Getd0Prong(2)*d->Getd0Prong(2);
	tmp[20]=d->InvMass2Prongs(0,1,321,321);
	tmp[21]=d->InvMass2Prongs(1,2,321,321);
	tmp[22]=d->CosPiDsLabFrameKKpi();	   
	tmp[23]=d->CosPiDsLabFramepiKK();	
	tmp[24]=d->CosPiKPhiRFrameKKpi();	   
	tmp[25]=d->CosPiKPhiRFramepiKK();	
	
 	
	fNtupleDs->Fill(tmp);
	PostData(4,fNtupleDs);
      }  
    }
    
    if(unsetvtx) d->UnsetOwnPrimaryVtx();
    if(recVtx)fProdCuts->CleanOwnPrimaryVtx(d,aod,origownvtx);
  }
 
  fCounter->StoreCandidates(aod,nSelectedloose,kTRUE);
  fCounter->StoreCandidates(aod,nSelectedtight,kFALSE);

  PostData(1,fOutput); 
  PostData(3,fCounter);    

  return;
}

//_________________________________________________________________

void AliAnalysisTaskSEDs::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
  if(fDebug > 1) printf("AnalysisTaskSEDs: Terminate() \n");
  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {     
    printf("ERROR: fOutput not available\n");
    return;
  }
  fHistNEvents = dynamic_cast<TH1F*>(fOutput->FindObject("hNEvents"));
  if(fHistNEvents){
    printf("Number of analyzed events = %d\n",(Int_t)fHistNEvents->GetBinContent(2));
  }else{
    printf("ERROR: fHistNEvents not available\n");
    return;
  }
  return;
}

