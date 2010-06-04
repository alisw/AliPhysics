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

/* $Id: $ */

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

ClassImp(AliAnalysisTaskSEDs)


//________________________________________________________________________
AliAnalysisTaskSEDs::AliAnalysisTaskSEDs():
  AliAnalysisTaskSE(),
  fOutput(0), 
  fHistNEvents(0),
  fReadMC(kFALSE),
  fNPtBins(0),
  fListCuts(0),
  fMassRange(0.2),
  fMassBinSize(0.002),
  fProdCuts(0),
  fAnalysisCuts(0)
{
  // Default constructor
}

//________________________________________________________________________
AliAnalysisTaskSEDs::AliAnalysisTaskSEDs(const char *name, AliRDHFCutsDstoKKpi* productioncuts, AliRDHFCutsDstoKKpi* analysiscuts):
  AliAnalysisTaskSE(name),
  fOutput(0),
  fHistNEvents(0),
  fReadMC(kFALSE),
  fNPtBins(0),
  fListCuts(0),
  fMassRange(0.2),
  fMassBinSize(0.002),
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
    for(Int_t i=0; i<fNPtBins+1; i++) printf(" Bin%d = %8.2f-%8.2f\n",i,fPtLimits[i],fPtLimits[i+1]);    
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
  if (fListCuts) {
    delete fListCuts;
    fListCuts = 0;
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
  AliRDHFCutsDstoKKpi *production = new AliRDHFCutsDstoKKpi();
  production=fProdCuts;
  AliRDHFCutsDstoKKpi *analysis = new AliRDHFCutsDstoKKpi();
  analysis=fAnalysisCuts;
  
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
  Double_t minMass=massDs-0.5*nInvMassBins*fMassBinSize;
  Double_t maxMass=massDs+0.5*nInvMassBins*fMassBinSize;

  TString hisname;
  Int_t index;
  for(Int_t i=0;i<fNPtBins;i++){
    index=GetHistoIndex(i);
    hisname.Form("hMassAllPt%d",i);
    fMassHist[index]=new TH1F(hisname.Data(),hisname.Data(),nInvMassBins,minMass,maxMass);
    fMassHist[index]->Sumw2();
    hisname.Form("hMassAllPt%dCuts",i);
    fMassHistCuts[index]=new TH1F(hisname.Data(),hisname.Data(),nInvMassBins,minMass,maxMass);
    fMassHistCuts[index]->Sumw2();
    hisname.Form("hCosPAllPt%d",i);
    fCosPHist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.5,1.);
    fCosPHist[index]->Sumw2();
    hisname.Form("hDLenAllPt%d",i);
    fDLenHist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.,0.5);
    fDLenHist[index]->Sumw2();
    hisname.Form("hSumd02AllPt%d",i);
    fSumd02Hist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.,1.);
    fSumd02Hist[index]->Sumw2();
     hisname.Form("hSigVertAllPt%d",i);
    fSigVertHist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.,0.1);
    fSigVertHist[index]->Sumw2();
    hisname.Form("hPtMaxAllPt%d",i);
    fPtMaxHist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.5,20.);
    fPtMaxHist[index]->Sumw2();
    hisname.Form("hPtCandAllPt%d",i);
    fPtCandHist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.5,20.);
    fPtCandHist[index]->Sumw2();
     hisname.Form("hDCAAllPt%d",i);
    fDCAHist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.,0.1);
    fDCAHist[index]->Sumw2();
    hisname.Form("hPtProng0AllPt%d",i);
    fPtProng0Hist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.0,20.);
    fPtProng0Hist[index]->Sumw2();
    hisname.Form("hPtProng1AllPt%d",i);
    fPtProng1Hist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.0,20.);
    fPtProng1Hist[index]->Sumw2();
    hisname.Form("hPtProng2AllPt%d",i);
    fPtProng2Hist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.0,20.);
    fPtProng2Hist[index]->Sumw2();

    hisname.Form("hDalitzAllPt%d",i);
    fDalitz[index]=new TH2F(hisname.Data(),hisname.Data(),100,0.,2.,100,0.,2.);
    fDalitz[index]->Sumw2();

    index=GetSignalHistoIndex(i);    
    hisname.Form("hMassSigPt%d",i);
    fMassHist[index]=new TH1F(hisname.Data(),hisname.Data(),nInvMassBins,minMass,maxMass);
    fMassHist[index]->Sumw2();
    hisname.Form("hMassSigPt%dCuts",i);
    fMassHistCuts[index]=new TH1F(hisname.Data(),hisname.Data(),nInvMassBins,minMass,maxMass);
    fMassHistCuts[index]->Sumw2();
    hisname.Form("hCosPSigPt%d",i);
    fCosPHist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.5,1.);
    fCosPHist[index]->Sumw2();
    hisname.Form("hDLenSigPt%d",i);
    fDLenHist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.,0.5);
    fDLenHist[index]->Sumw2();
    hisname.Form("hSumd02SigPt%d",i);
    fSumd02Hist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.,1.);
    fSumd02Hist[index]->Sumw2();
    hisname.Form("hSigVertSigPt%d",i);
    fSigVertHist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.,0.1);
    fSigVertHist[index]->Sumw2();
    hisname.Form("hPtMaxSigPt%d",i);
    fPtMaxHist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.0,20.);
    fPtMaxHist[index]->Sumw2();
    hisname.Form("hPtCandSigPt%d",i);
    fPtCandHist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.0,20.);
    fPtCandHist[index]->Sumw2();
    hisname.Form("hDCASigPt%d",i);
    fDCAHist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.,0.1);
    fDCAHist[index]->Sumw2();
    hisname.Form("hPtProng0SigPt%d",i);
    fPtProng0Hist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.0,20.);
    fPtProng0Hist[index]->Sumw2();
    hisname.Form("hPtProng1SigPt%d",i);
    fPtProng1Hist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.0,20.);
    fPtProng1Hist[index]->Sumw2();
    hisname.Form("hPtProng2SigPt%d",i);
    fPtProng2Hist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.0,20.);
    fPtProng2Hist[index]->Sumw2();
    hisname.Form("hDalitzSigPt%d",i);
    fDalitz[index]=new TH2F(hisname.Data(),hisname.Data(),100,0.,2.,100,0.,2.);
    fDalitz[index]->Sumw2();

    index=GetBackgroundHistoIndex(i);    
    hisname.Form("hMassBkgPt%d",i);
    fMassHist[index]=new TH1F(hisname.Data(),hisname.Data(),nInvMassBins,minMass,maxMass);
    fMassHist[index]->Sumw2();
    hisname.Form("hMassBkgPt%dCuts",i);
    fMassHistCuts[index]=new TH1F(hisname.Data(),hisname.Data(),nInvMassBins,minMass,maxMass);
    fMassHistCuts[index]->Sumw2();
    hisname.Form("hCosPBkgPt%d",i);
    fCosPHist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.5,1.);
    fCosPHist[index]->Sumw2();
    hisname.Form("hDLenBkgPt%d",i);
    fDLenHist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.,0.5);
    fDLenHist[index]->Sumw2();
    hisname.Form("hSumd02BkgPt%d",i);
    fSumd02Hist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.,1.);
    fSumd02Hist[index]->Sumw2();
    hisname.Form("hSigVertBkgPt%d",i);
    fSigVertHist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.,0.1);
    fSigVertHist[index]->Sumw2();
    hisname.Form("hPtMaxSigBkg%d",i);
    fPtMaxHist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.0,20.);
    fPtMaxHist[index]->Sumw2();
    hisname.Form("hPtCandBkgPt%d",i);
    fPtCandHist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.0,20.);
    fPtCandHist[index]->Sumw2();
    hisname.Form("hDCABkgPt%d",i);
    fDCAHist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.,0.1);
    fDCAHist[index]->Sumw2();
    hisname.Form("hPtProng0BkgPt%d",i);
    fPtProng0Hist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.0,20.);
    fPtProng0Hist[index]->Sumw2();
    hisname.Form("hPtProng1BkgPt%d",i);
    fPtProng1Hist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.0,20.);
    fPtProng1Hist[index]->Sumw2();
    hisname.Form("hPtProng2BkgPt%d",i);
    fPtProng2Hist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.0,20.);
    fPtProng2Hist[index]->Sumw2();
    hisname.Form("hDalitzBkgPt%d",i);
    fDalitz[index]=new TH2F(hisname.Data(),hisname.Data(),100,0.,2.,100,0.,2.);
    fDalitz[index]->Sumw2();

    index=GetReflSignalHistoIndex(i);    
    hisname.Form("hMassReflSigPt%d",i);
    fMassHist[index]=new TH1F(hisname.Data(),hisname.Data(),100,minMass,maxMass);
    fMassHist[index]->Sumw2();
    hisname.Form("hMassReflSigPt%dCuts",i);
    fMassHistCuts[index]=new TH1F(hisname.Data(),hisname.Data(),100,minMass,maxMass);
    fMassHistCuts[index]->Sumw2();
    hisname.Form("hCosPReflSigPt%d",i);
    fCosPHist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.5,1.);
    fCosPHist[index]->Sumw2();
    hisname.Form("hDLenReflSigPt%d",i);
    fDLenHist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.,0.5);
    fDLenHist[index]->Sumw2();
    hisname.Form("hSumd02ReflSigPt%d",i);
    fSumd02Hist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.,1.);
    fSumd02Hist[index]->Sumw2();
    hisname.Form("hSigVertReflSigPt%d",i);
    fSigVertHist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.,0.1);
    fSigVertHist[index]->Sumw2();
    hisname.Form("hPtMaxReflSigPt%d",i);
    fPtMaxHist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.0,20.);
    fPtMaxHist[index]->Sumw2();
    hisname.Form("hPtCandReflSigPt%d",i);
    fPtCandHist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.0,20.);
    fPtCandHist[index]->Sumw2();
    hisname.Form("hDCAReflSigPt%d",i);
    fDCAHist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.,0.1);
    fDCAHist[index]->Sumw2();
    hisname.Form("hPtProng0ReflSigPt%d",i);
    fPtProng0Hist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.0,20.);
    fPtProng0Hist[index]->Sumw2();
    hisname.Form("hPtProng1ReflSigPt%d",i);
    fPtProng1Hist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.0,20.);
    fPtProng1Hist[index]->Sumw2();
    hisname.Form("hPtProng2ReflSigPt%d",i);
    fPtProng2Hist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.0,20.);
    fPtProng2Hist[index]->Sumw2();
    hisname.Form("hDalitzReflSigPt%d",i);
    fDalitz[index]=new TH2F(hisname.Data(),hisname.Data(),100,0.,2.,100,0.,2.);
    fDalitz[index]->Sumw2();
  }

  for(Int_t i=0; i<4*fNPtBins; i++){
    fOutput->Add(fMassHist[i]);
    fOutput->Add(fMassHistCuts[i]);
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
  }

  fChanHist[0] = new TH1F("hChanAll", "KKpi and piKK candidates",4,-0.5,3.5);
  fChanHist[1] = new TH1F("hChanSig", "KKpi and piKK candidates",4,-0.5,3.5);
  fChanHist[2] = new TH1F("hChanBkg", "KKpi and piKK candidates",4,-0.5,3.5);
  fChanHist[3] = new TH1F("hChanReflSig", "KKpi and piKK candidates",4,-0.5,3.5);
  fChanHistCuts[0] = new TH1F("hChanAllCuts", "KKpi and piKK candidates",4,-0.5,3.5);
  fChanHistCuts[1] = new TH1F("hChanSigCuts", "KKpi and piKK candidates",4,-0.5,3.5);
  fChanHistCuts[2] = new TH1F("hChanBkgCuts", "KKpi and piKK candidates",4,-0.5,3.5);
  fChanHistCuts[3] = new TH1F("hChanReflSigCuts", "KKpi and piKK candidates",4,-0.5,3.5);
  for(Int_t i=0;i<4;i++){
    fChanHist[i]->Sumw2();
    fChanHist[i]->SetMinimum(0);
    fChanHistCuts[i]->Sumw2();
    fChanHistCuts[i]->SetMinimum(0);
    fOutput->Add(fChanHist[i]);
    fOutput->Add(fChanHistCuts[i]);
  }

  fHistNEvents = new TH1F("hNEvents", "Number of processed events",3,-1.5,1.5);
  fHistNEvents->Sumw2();
  fHistNEvents->SetMinimum(0);
  fOutput->Add(fHistNEvents);


  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEDs::UserExec(Option_t */*option*/)
{
  // Ds selection for current event, fill mass histos and selecetion variable histo
  // separate signal and backgound if fReadMC is activated

  AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());
  fHistNEvents->Fill(0); // count event
  // Post the data already here
  PostData(1,fOutput);
  

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
  } else {
    array3Prong=(TClonesArray*)aod->GetList()->FindObject("Charm3Prong");
  }

  if(!array3Prong) {
    printf("AliAnalysisTaskSEDs::UserExec: Charm3Prong branch not found!\n");
    return;
  }

 
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
  for (Int_t i3Prong = 0; i3Prong < n3Prong; i3Prong++) {
    AliAODRecoDecayHF3Prong *d = (AliAODRecoDecayHF3Prong*)array3Prong->UncheckedAt(i3Prong);
    
    
    Bool_t unsetvtx=kFALSE;
    if(!d->GetOwnPrimaryVtx()){
      d->SetOwnPrimaryVtx(vtx1);
      unsetvtx=kTRUE;
    }

    Int_t retCodeProductionCuts=fProdCuts->IsSelected(d,AliRDHFCuts::kCandidate);
    if(retCodeProductionCuts>0){
      Int_t isKKpi=retCodeProductionCuts&1;
      Int_t ispiKK=retCodeProductionCuts&2;
//       Int_t isPhi=retCodeProductionCuts&4;
//       Int_t isK0star=retCodeProductionCuts&8;
      Double_t ptCand = d->Pt();
      Int_t iPtBin=TMath::BinarySearch(fNPtBins,fPtLimits,(Float_t)ptCand);
      Int_t retCodeAnalysisCuts=fAnalysisCuts->IsSelected(d,AliRDHFCuts::kCandidate);
      Int_t isKKpiAC=retCodeAnalysisCuts&1;
      Int_t ispiKKAC=retCodeAnalysisCuts&2;
//       Int_t isPhiAC=retCodeAnalysisCuts&4;
//       Int_t isK0starAC=retCodeAnalysisCuts&8;

      Int_t labDs=-1;
      if(fReadMC){
	labDs = d->MatchToMC(431,arrayMC,3,pdgDstoKKpi);
      }

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

      Int_t index=GetHistoIndex(iPtBin);
      Int_t type=0;
      if(isKKpi) type+=1;
      if(ispiKK) type+=2;
      Int_t typeAC=0;
      if(isKKpiAC) typeAC+=1;
      if(ispiKKAC) typeAC+=2;
      fCosPHist[index]->Fill(cosp);
      fDLenHist[index]->Fill(dlen);
      fSigVertHist[index]->Fill(sigvert);
      fSumd02Hist[index]->Fill(sumD02);
      fPtMaxHist[index]->Fill(ptmax);
      fPtCandHist[index]->Fill(ptCand);
      fDCAHist[index]->Fill(dca);
      fChanHist[0]->Fill(type);
      fPtProng0Hist[index]->Fill(pt0);
      fPtProng1Hist[index]->Fill(pt1);
      fPtProng2Hist[index]->Fill(pt2);

      if(retCodeAnalysisCuts>0) fChanHistCuts[0]->Fill(typeAC);
      if(fReadMC){
	if(labDs>=0) {	  
	  index=GetSignalHistoIndex(iPtBin);
	  fChanHist[1]->Fill(type);	  
	  if(retCodeAnalysisCuts>0) fChanHistCuts[1]->Fill(typeAC);
	}else{
	  index=GetBackgroundHistoIndex(iPtBin);
	  fChanHist[2]->Fill(type);	  
	  if(retCodeAnalysisCuts>0) fChanHistCuts[2]->Fill(typeAC);
	}
      }
      if(isKKpi){
	index=GetHistoIndex(iPtBin);
	Double_t invMass=d->InvMassDsKKpi();
	fMassHist[index]->Fill(invMass);
	Double_t mass01=d->InvMass2Prongs(0,1,321,321);
	Double_t mass12=d->InvMass2Prongs(1,2,321,211);
	fDalitz[index]->Fill(mass01,mass12);
	if(retCodeAnalysisCuts>0 && isKKpiAC) fMassHistCuts[index]->Fill(invMass);
	if(fReadMC){
	  Int_t labDau0=((AliAODTrack*)d->GetDaughter(0))->GetLabel();
	  AliAODMCParticle* p=(AliAODMCParticle*)arrayMC->UncheckedAt(labDau0);
	  Int_t pdgCode0=TMath::Abs(p->GetPdgCode());
	  if(labDs>=0){
	    if(pdgCode0==321) {	  
	      index=GetSignalHistoIndex(iPtBin);
	    }else{
	      index=GetReflSignalHistoIndex(iPtBin);
	    }
	  }else{
	    index=GetBackgroundHistoIndex(iPtBin);
	  }
	  fMassHist[index]->Fill(invMass);
	  fCosPHist[index]->Fill(cosp);
	  fDLenHist[index]->Fill(dlen);
	  fSigVertHist[index]->Fill(sigvert);
	  fSumd02Hist[index]->Fill(sumD02);
	  fPtMaxHist[index]->Fill(ptmax);
	  fPtCandHist[index]->Fill(ptCand);
	  fDCAHist[index]->Fill(dca);
	  fPtProng0Hist[index]->Fill(pt0);
	  fPtProng1Hist[index]->Fill(pt1);
	  fPtProng2Hist[index]->Fill(pt2);
	  fDalitz[index]->Fill(mass01,mass12);
	  if(retCodeAnalysisCuts>0 && isKKpiAC) fMassHistCuts[index]->Fill(invMass);	  
	}	
      }
      if(ispiKK){
	index=GetHistoIndex(iPtBin);
	Double_t invMass=d->InvMassDspiKK();
	fMassHist[index]->Fill(invMass);
	if(retCodeAnalysisCuts>0 && ispiKKAC) fMassHistCuts[index]->Fill(invMass);
	if(fReadMC){
	  Int_t labDau0=((AliAODTrack*)d->GetDaughter(0))->GetLabel();
	  AliAODMCParticle* p=(AliAODMCParticle*)arrayMC->UncheckedAt(labDau0);
	  Int_t pdgCode0=TMath::Abs(p->GetPdgCode());
	  if(labDs>=0) {	  
	    if(pdgCode0==211) {	  
	      index=GetSignalHistoIndex(iPtBin);
	    }else{
	      index=GetReflSignalHistoIndex(iPtBin);
	    }
	  }else{
	    index=GetBackgroundHistoIndex(iPtBin);
	  }
	  fMassHist[index]->Fill(invMass);
	  fCosPHist[index]->Fill(cosp);
	  fDLenHist[index]->Fill(dlen);
	  fSigVertHist[index]->Fill(sigvert);
	  fSumd02Hist[index]->Fill(sumD02);
	  fPtMaxHist[index]->Fill(ptmax);
	  fPtCandHist[index]->Fill(ptCand);
	  fDCAHist[index]->Fill(dca);
	  fPtProng0Hist[index]->Fill(pt0);
	  fPtProng1Hist[index]->Fill(pt1);
	  fPtProng2Hist[index]->Fill(pt2);
	  if(retCodeAnalysisCuts>0 && ispiKKAC) fMassHistCuts[index]->Fill(invMass);
	}
      }
    }
    if(unsetvtx) d->UnsetOwnPrimaryVtx();
  }
 
   
  PostData(1,fOutput);    
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
  fChanHist[0] = dynamic_cast<TH1F*>(fOutput->FindObject("hChanAll"));
  fChanHist[1] = dynamic_cast<TH1F*>(fOutput->FindObject("hChanSig"));
  fChanHist[2] = dynamic_cast<TH1F*>(fOutput->FindObject("hChanBkg"));
  fChanHist[3] = dynamic_cast<TH1F*>(fOutput->FindObject("hChanReflSig"));
  fChanHistCuts[0] = dynamic_cast<TH1F*>(fOutput->FindObject("hChanAllCuts"));
  fChanHistCuts[1] = dynamic_cast<TH1F*>(fOutput->FindObject("hChanSigCuts"));
  fChanHistCuts[2] = dynamic_cast<TH1F*>(fOutput->FindObject("hChanBkgCuts"));
  fChanHistCuts[3] = dynamic_cast<TH1F*>(fOutput->FindObject("hChanReflSigCuts"));


  TString hisname;
  Int_t index;
  for(Int_t i=0;i<fNPtBins;i++){

    index=GetHistoIndex(i);
    hisname.Form("hMassAllPt%d",i);
    fMassHist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hMassAllPt%dCuts",i);
    fMassHistCuts[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hCosPAllPt%d",i);
    fCosPHist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hDLenAllPt%d",i);
    fDLenHist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hSumd02AllPt%d",i);
    fSumd02Hist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hSigVertAllPt%d",i);
    fSigVertHist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hPtMaxAllPt%d",i);
    fPtMaxHist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hPtCandAllPt%d",i);
    fPtCandHist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hDCAAllPt%d",i);
    fDCAHist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hPtProng0AllPt%d",i);
    fPtProng0Hist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hPtProng1AllPt%d",i);
    fPtProng1Hist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hPtProng2AllPt%d",i);
    fPtProng2Hist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hDalitzAllPt%d",i);
    fDalitz[index]=dynamic_cast<TH2F*>(fOutput->FindObject(hisname.Data()));

    index=GetSignalHistoIndex(i);    
    hisname.Form("hMassSigPt%d",i);
    fMassHist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hMassSigPt%dCuts",i);
    fMassHistCuts[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hCosPSigPt%d",i);
    fCosPHist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hDLenSigPt%d",i);
    fDLenHist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hSumd02SigPt%d",i);
    fSumd02Hist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hSigVertSigPt%d",i);
    fSigVertHist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hPtMaxSigPt%d",i);
    fPtMaxHist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hPtCandSigPt%d",i);
    fPtCandHist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hDCASigPt%d",i);
    fDCAHist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hPtProng0SigPt%d",i);
    fPtProng0Hist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hPtProng1SigPt%d",i);
    fPtProng1Hist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hPtProng2SigPt%d",i);
    fPtProng2Hist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hDalitzSigPt%d",i);
    fDalitz[index]=dynamic_cast<TH2F*>(fOutput->FindObject(hisname.Data()));

    index=GetBackgroundHistoIndex(i);    
    hisname.Form("hMassBkgPt%d",i);
    fMassHist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hMassBkgPt%dCuts",i);
    fMassHistCuts[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hCosPBkgPt%d",i);
    fCosPHist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hDLenBkgPt%d",i);
    fDLenHist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hSumd02BkgPt%d",i);
    fSumd02Hist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hSigVertBkgPt%d",i);
    fSigVertHist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hPtMaxBkgPt%d",i);
    fPtMaxHist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hPtCandBkgPt%d",i);
    fPtCandHist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hDCABkgPt%d",i);
    fDCAHist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hPtProng0BkgPt%d",i);
    fPtProng0Hist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hPtProng1BkgPt%d",i);
    fPtProng1Hist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hPtProng2BkgPt%d",i);
    fPtProng2Hist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hDalitzBkgPt%d",i);
    fDalitz[index]=dynamic_cast<TH2F*>(fOutput->FindObject(hisname.Data()));

    index=GetReflSignalHistoIndex(i);    
    hisname.Form("hMassReflSigPt%d",i);
    fMassHist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hMassReflSigPt%dCuts",i);
    fMassHistCuts[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hCosPReflSigPt%d",i);
    fCosPHist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hDLenReflSigPt%d",i);
    fDLenHist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hSumd02ReflSigPt%d",i);
    fSumd02Hist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hSigVertReflSigPt%d",i);
    fSigVertHist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hPtMaxReflSigPt%d",i);
    fPtMaxHist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hPtCandReflSigPt%d",i);
    fPtCandHist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hDCAReflSigPt%d",i);
    fDCAHist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hPtProng0ReflSigPt%d",i);
    fPtProng0Hist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hPtProng1ReflSigPt%d",i);
    fPtProng1Hist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hPtProng2ReflSigPt%d",i);
    fPtProng2Hist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hDalitzReflSigPt%d",i);
    fDalitz[index]=dynamic_cast<TH2F*>(fOutput->FindObject(hisname.Data()));

  }
  return;
}
