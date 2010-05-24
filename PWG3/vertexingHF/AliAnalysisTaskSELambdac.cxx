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
// AliAnalysisTaskSE for the extraction of signal(e.g D+) of heavy flavor
// decay candidates with the MC truth.
// Authors: Renu Bala, bala@to.infn.it
// F. Prino, prino@to.infn.it
// G. Ortona, ortona@to.infn.it
/////////////////////////////////////////////////////////////

#include <TClonesArray.h>
#include <TNtuple.h>
#include <TCanvas.h>
#include <TList.h>
#include <TString.h>
#include <TH1F.h>
#include <TH2F.h>
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
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskSELambdac.h"
#include "AliKFParticle.h"
#include "AliAODPidHF.h"
#include "AliRDHFCutsLctopKpi.h"
#include "AliRDHFCuts.h"

ClassImp(AliAnalysisTaskSELambdac)


//________________________________________________________________________
AliAnalysisTaskSELambdac::AliAnalysisTaskSELambdac():
AliAnalysisTaskSE(),
fOutput(0), 
fHistNEvents(0),
fTPCSignal3Sigma(0),
fTPCSignal3SigmaReK(0),
fTPCSignal3SigmaRep(0),
fTPCSignal2Sigma(0),
fTPCSignal2SigmaReK(0),
fTPCSignal2SigmaRep(0),
fNtupleLambdac(0),
fUpmasslimit(2.486),
fLowmasslimit(2.086),
fNPtBins(0),
fRDCutsAnalysis(0),
fRDCutsProduction(0),
fListCuts(0),
fFillNtuple(kFALSE),
fReadMC(kFALSE),
fMCPid(kFALSE),
fRealPid(kFALSE),
fResPid(kTRUE),
fUseKF(kFALSE),
fVHF(0)
{
   // Default constructor
}

//________________________________________________________________________
AliAnalysisTaskSELambdac::AliAnalysisTaskSELambdac(const char *name,Bool_t fillNtuple,AliRDHFCutsLctopKpi *lccutsana,AliRDHFCutsLctopKpi *lccutsprod):
AliAnalysisTaskSE(name),
fOutput(0),
fHistNEvents(0),
fTPCSignal3Sigma(0),
fTPCSignal3SigmaReK(0),
fTPCSignal3SigmaRep(0),
fTPCSignal2Sigma(0),
fTPCSignal2SigmaReK(0),
fTPCSignal2SigmaRep(0),
fNtupleLambdac(0),
fUpmasslimit(2.486),
fLowmasslimit(2.086),
fNPtBins(0),
fRDCutsAnalysis(lccutsana),
fRDCutsProduction(lccutsprod),
fListCuts(0),
fFillNtuple(fillNtuple),
fReadMC(kFALSE),
fMCPid(kFALSE),
fRealPid(kTRUE),
fResPid(kFALSE),
fUseKF(kFALSE),
fVHF(0)
{
  //Double_t ptlim[5]={0.,2.,3.,5,9999999.};
  //SetPtBinLimit(5, ptlim);
   SetPtBinLimit(fRDCutsAnalysis->GetNPtBins()+1,fRDCutsAnalysis->GetPtBinLimits());
  // Default constructor
   // Output slot #1 writes into a TList container
  DefineOutput(1,TList::Class());  //My private output
  DefineOutput(2,TList::Class());

  if(fFillNtuple){
    // Output slot #2 writes into a TNtuple container
    DefineOutput(3,TNtuple::Class());  //My private output
  }
}

//________________________________________________________________________
AliAnalysisTaskSELambdac::~AliAnalysisTaskSELambdac()
{
  // Destructor
  if (fOutput) {
    delete fOutput;
    fOutput = 0;
  }
  if (fVHF) {
    delete fVHF;
    fVHF = 0;
  }
  
 if(fRDCutsAnalysis){
    delete fRDCutsAnalysis;
    fRDCutsAnalysis = 0;
 }
 if(fRDCutsProduction){
    delete fRDCutsProduction;
    fRDCutsProduction = 0;
 }

 if (fListCuts) {
   delete fListCuts;
   fListCuts = 0;
 }
}  
//_________________________________________________________________
void  AliAnalysisTaskSELambdac::SetMassLimits(Float_t range){
  fUpmasslimit = 2.286+range;
  fLowmasslimit = 2.286-range;
}
//_________________________________________________________________
void  AliAnalysisTaskSELambdac::SetMassLimits(Float_t lowlimit, Float_t uplimit){
  if(uplimit>lowlimit)
    {
      fUpmasslimit = lowlimit;
      fLowmasslimit = uplimit;
    }
}


//________________________________________________________________________
void AliAnalysisTaskSELambdac::SetPtBinLimit(Int_t n, Float_t* lim){
  // define pt bins for analysis
  if(n>kMaxPtBins){
    printf("Max. number of Pt bins = %d\n",kMaxPtBins);
    fNPtBins=kMaxPtBins;
    fArrayBinLimits[0]=0.;
    fArrayBinLimits[1]=2.;
    fArrayBinLimits[2]=3.;
    fArrayBinLimits[3]=5.;
    for(Int_t i=4; i<kMaxPtBins+1; i++) fArrayBinLimits[i]=99999999.;
  }else{
    fNPtBins=n-1;
    fArrayBinLimits[0]=lim[0];
    for(Int_t i=1; i<fNPtBins+1; i++) 
      if(lim[i]>fArrayBinLimits[i-1]){
	fArrayBinLimits[i]=lim[i];
      }
      else {
	fArrayBinLimits[i]=fArrayBinLimits[i-1];
      }
    for(Int_t i=fNPtBins; i<kMaxPtBins+1; i++) fArrayBinLimits[i]=99999999.;
  }
  if(fDebug > 1){
    printf("Number of Pt bins = %d\n",fNPtBins);
    for(Int_t i=0; i<fNPtBins+1; i++) printf(" Bin%d = %8.2f-%8.2f\n",i,fArrayBinLimits[i],fArrayBinLimits[i+1]);    
  }
}
//_________________________________________________________________
Double_t  AliAnalysisTaskSELambdac::GetPtBinLimit(Int_t ibin){
  if(ibin>fNPtBins)return -1;
  return fArrayBinLimits[ibin];
} 

//_________________________________________________________________
void AliAnalysisTaskSELambdac::Init()
{
  // Initialization

  if(fDebug > 1) printf("AnalysisTaskSELambdac::Init() \n");

  fListCuts=new TList();
  AliRDHFCutsLctopKpi *production = new AliRDHFCutsLctopKpi();
    production=fRDCutsProduction;

  AliRDHFCutsLctopKpi *analysis = new AliRDHFCutsLctopKpi();
  analysis=fRDCutsAnalysis;

   fListCuts->Add(analysis);
   fListCuts->Add(production);
   PostData(2,fListCuts);
  return;
}

//________________________________________________________________________
void AliAnalysisTaskSELambdac::UserCreateOutputObjects()
{
  // Create the output container
  //
  if(fDebug > 1) printf("AnalysisTaskSELambdac::UserCreateOutputObjects() \n");

  // Several histograms are more conveniently managed in a TList
  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("OutputHistos");

  TString hisname;
  Int_t index=0;
  Int_t indexLS=0;
  for(Int_t i=0;i<fNPtBins;i++){

    index=GetHistoIndex(i);
    indexLS=GetLSHistoIndex(i);

    hisname.Form("hMassPt%d",i);
    fMassHist[index]=new TH1F(hisname.Data(),hisname.Data(),100,fLowmasslimit,fUpmasslimit);
    fMassHist[index]->Sumw2();
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
    fPtMaxHist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.5,5.);
    fPtMaxHist[index]->Sumw2();

    hisname.Form("hDCAAllPt%d",i);
    fDCAHist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.,0.1);
    fDCAHist[index]->Sumw2();



    hisname.Form("hMassPt%dTC",i);
    fMassHistTC[index]=new TH1F(hisname.Data(),hisname.Data(),100,fLowmasslimit,fUpmasslimit);
    fMassHistTC[index]->Sumw2();




    
    hisname.Form("hCosPAllPt%dLS",i);
    fCosPHistLS[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.5,1.);
    fCosPHistLS[index]->Sumw2();
    hisname.Form("hDLenAllPt%dLS",i);
    fDLenHistLS[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.,0.5);
    fDLenHistLS[index]->Sumw2();
    hisname.Form("hSumd02AllPt%dLS",i);
    fSumd02HistLS[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.,1.);
    fSumd02HistLS[index]->Sumw2();
    hisname.Form("hSigVertAllPt%dLS",i);
    fSigVertHistLS[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.,0.1);
    fSigVertHistLS[index]->Sumw2();
    hisname.Form("hPtMaxAllPt%dLS",i);
    fPtMaxHistLS[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.5,5.);
    fPtMaxHistLS[index]->Sumw2();
    
    hisname.Form("hDCAAllPt%dLS",i);
    fDCAHistLS[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.,0.1);
    fDCAHistLS[index]->Sumw2();
    
    hisname.Form("hLSPt%dLC",i);
    fMassHistLS[indexLS] = new TH1F(hisname.Data(),hisname.Data(),100,fLowmasslimit,fUpmasslimit);
    fMassHistLS[indexLS]->Sumw2();
    
    hisname.Form("hLSPt%dTC",i);
    fMassHistLSTC[indexLS] = new TH1F(hisname.Data(),hisname.Data(),100,fLowmasslimit,fUpmasslimit);
    fMassHistLSTC[indexLS]->Sumw2();


    
    index=GetSignalHistoIndex(i);    
    indexLS++;
    hisname.Form("hSigPt%d",i);
    fMassHist[index]=new TH1F(hisname.Data(),hisname.Data(),100,fLowmasslimit,fUpmasslimit);
    fMassHist[index]->Sumw2();
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
    fPtMaxHist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.5,5.);
    fPtMaxHist[index]->Sumw2();    

    hisname.Form("hDCASigPt%d",i);
    fDCAHist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.,0.1);
    fDCAHist[index]->Sumw2();    


    hisname.Form("hSigPt%dTC",i);
    fMassHistTC[index]=new TH1F(hisname.Data(),hisname.Data(),100,fLowmasslimit,fUpmasslimit);
    fMassHistTC[index]->Sumw2();

    hisname.Form("hLSPt%dLCnw",i);
    fMassHistLS[indexLS]=new TH1F(hisname.Data(),hisname.Data(),100,fLowmasslimit,fUpmasslimit);
    fMassHistLS[indexLS]->Sumw2();
    hisname.Form("hLSPt%dTCnw",i);
    fMassHistLSTC[indexLS]=new TH1F(hisname.Data(),hisname.Data(),100,fLowmasslimit,fUpmasslimit);
    fMassHistLSTC[indexLS]->Sumw2();


    
    hisname.Form("hCosPSigPt%dLS",i);
    fCosPHistLS[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.5,1.);
    fCosPHistLS[index]->Sumw2();
    hisname.Form("hDLenSigPt%dLS",i);
    fDLenHistLS[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.,0.5);
    fDLenHistLS[index]->Sumw2();
    hisname.Form("hSumd02SigPt%dLS",i);
    fSumd02HistLS[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.,1.);
    fSumd02HistLS[index]->Sumw2();
    hisname.Form("hSigVertSigPt%dLS",i);
    fSigVertHistLS[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.,0.1);
    fSigVertHistLS[index]->Sumw2();
    hisname.Form("hPtMaxSigPt%dLS",i);
    fPtMaxHistLS[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.5,5.);
    fPtMaxHistLS[index]->Sumw2();

    hisname.Form("hDCASigPt%dLS",i);
    fDCAHistLS[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.,0.1);
    fDCAHistLS[index]->Sumw2();
    


    index=GetBackgroundHistoIndex(i); 
    indexLS++;
    hisname.Form("hBkgPt%d",i);
    fMassHist[index]=new TH1F(hisname.Data(),hisname.Data(),100,fLowmasslimit,fUpmasslimit);
    fMassHist[index]->Sumw2();
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
    hisname.Form("hPtMaxBkgPt%d",i);
    fPtMaxHist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.5,5.);
    fPtMaxHist[index]->Sumw2();

    hisname.Form("hDCABkgPt%d",i);
    fDCAHist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.,0.1);
    fDCAHist[index]->Sumw2();


    hisname.Form("hBkgPt%dTC",i);
    fMassHistTC[index]=new TH1F(hisname.Data(),hisname.Data(),100,fLowmasslimit,fUpmasslimit);
    fMassHistTC[index]->Sumw2();

    hisname.Form("hLSPt%dLCntrip",i);
    fMassHistLS[indexLS]=new TH1F(hisname.Data(),hisname.Data(),100,fLowmasslimit,fUpmasslimit);
    fMassHistLS[indexLS]->Sumw2();
    hisname.Form("hLSPt%dTCntrip",i);
    fMassHistLSTC[indexLS]=new TH1F(hisname.Data(),hisname.Data(),100,fLowmasslimit,fUpmasslimit);
    fMassHistLSTC[indexLS]->Sumw2();

    
    hisname.Form("hCosPBkgPt%dLS",i);
    fCosPHistLS[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.5,1.);
    fCosPHistLS[index]->Sumw2();
    hisname.Form("hDLenBkgPt%dLS",i);
    fDLenHistLS[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.,0.5);
    fDLenHistLS[index]->Sumw2();
    hisname.Form("hSumd02BkgPt%dLS",i);
    fSumd02HistLS[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.,1.);
    fSumd02HistLS[index]->Sumw2();
    hisname.Form("hSigVertBkgPt%dLS",i);
    fSigVertHistLS[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.,0.1);
    fSigVertHistLS[index]->Sumw2();
    hisname.Form("hPtMaxBkgPt%dLS",i);
    fPtMaxHistLS[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.5,5.);
    fPtMaxHistLS[index]->Sumw2();
    hisname.Form("hDCABkgPt%dLS",i);
    fDCAHistLS[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.,0.1);
    fDCAHistLS[index]->Sumw2();
    

    indexLS++;
    hisname.Form("hLSPt%dLCntripsinglecut",i);
    fMassHistLS[indexLS]=new TH1F(hisname.Data(),hisname.Data(),100,fLowmasslimit,fUpmasslimit);
    fMassHistLS[indexLS]->Sumw2();
    hisname.Form("hLSPt%dTCntripsinglecut",i);
    fMassHistLSTC[indexLS]=new TH1F(hisname.Data(),hisname.Data(),100,fLowmasslimit,fUpmasslimit);
    fMassHistLSTC[indexLS]->Sumw2();

    indexLS++;
    hisname.Form("hLSPt%dLCspc",i);
    fMassHistLS[indexLS]=new TH1F(hisname.Data(),hisname.Data(),100,fLowmasslimit,fUpmasslimit);
    fMassHistLS[indexLS]->Sumw2();
    hisname.Form("hLSPt%dTCspc",i);
    fMassHistLSTC[indexLS]=new TH1F(hisname.Data(),hisname.Data(),100,fLowmasslimit,fUpmasslimit);
    fMassHistLSTC[indexLS]->Sumw2();
  }

  fTPCSignal3Sigma=new TH2F("fTPCSignal3Sigma","fTPCSignal3Sigma",100,0, 10, 100, 0, 100);
  fTPCSignal3SigmaReK=new TH2F("fTPCSignal3SigmaReK","fTPCSignal3SigmaRe",100,0, 10, 100, 0, 100);
  fTPCSignal3SigmaRep=new TH2F("fTPCSignal3SigmaRep","fTPCSignal3SigmaRe",100,0, 10, 100, 0, 100);
  fTPCSignal2Sigma=new TH2F("fTPCSignal2Sigma","fTPCSignal2Sigma",100,0, 10, 100, 0, 100);
  fTPCSignal2SigmaReK=new TH2F("fTPCSignal2SigmaReK","fTPCSignal2SigmaRe",100,0, 10, 100, 0, 100);
  fTPCSignal2SigmaRep=new TH2F("fTPCSignal2SigmaRep","fTPCSignal2SigmaRe",100,0, 10, 100, 0, 100);
  
  for(Int_t i=0; i<3*fNPtBins; i++){
    fOutput->Add(fMassHist[i]);
    fOutput->Add(fMassHistTC[i]);
    fOutput->Add(fCosPHist[i]);
    fOutput->Add(fDLenHist[i]);
    fOutput->Add(fSumd02Hist[i]);
    fOutput->Add(fSigVertHist[i]);
    fOutput->Add(fPtMaxHist[i]);
    fOutput->Add(fDCAHist[i]);
  }

    fOutput->Add(fTPCSignal3Sigma);
    fOutput->Add(fTPCSignal3SigmaReK);
    fOutput->Add(fTPCSignal3SigmaRep);
    fOutput->Add(fTPCSignal2Sigma);
    fOutput->Add(fTPCSignal2SigmaReK);
    fOutput->Add(fTPCSignal2SigmaRep);
  

  fHistNEvents = new TH1F("fHistNEvents", "Number of processed events; ; Events",3,-1.5,1.5);
  fHistNEvents->Sumw2();
  fHistNEvents->SetMinimum(0);
  fOutput->Add(fHistNEvents);
  
  

  if(fFillNtuple){
    //OpenFile(3); // 2 is the slot number of the ntuple
   
    fNtupleLambdac = new TNtuple("fNtupleLambdac","D +","pdg:Px:Py:Pz:PtTrue:VxTrue:VyTrue:VzTrue:Ptpi:PtK:Ptpi2:PtRec:PointingAngle:DecLeng:VxRec:VyRec:VzRec:InvMass:sigvert:d0Pi:d0K:d0Pi2:dca:d0square");  
    
  }
  
  return;
}

//________________________________________________________________________
void AliAnalysisTaskSELambdac::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event:
  // heavy flavor candidates association to MC truth

   AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());
   //tmp!
 fHistNEvents->Fill(0); // count event
  // Post the data already here
  PostData(1,fOutput);

  TClonesArray *array3Prong = 0;
  TClonesArray *arrayLikeSign =0;
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
      arrayLikeSign=(TClonesArray*)aodFromExt->GetList()->FindObject("LikeSign3Prong");
    }
  } else {
    array3Prong=(TClonesArray*)aod->GetList()->FindObject("Charm3Prong");
    arrayLikeSign=(TClonesArray*)aod->GetList()->FindObject("LikeSign3Prong");
  }

  if(!array3Prong) {
    printf("AliAnalysisTaskSELambdac::UserExec: Charm3Prong branch not found!\n");
    return;
  }
  if(!arrayLikeSign) {
    printf("AliAnalysisTaskSELambdac::UserExec: LikeSign3Prong branch not found!\n");
  //  return;
  }

 
  TClonesArray *arrayMC=0;
  AliAODMCHeader *mcHeader=0;

  // AOD primary vertex
  AliAODVertex *vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();
  
  // load MC particles
  if(fReadMC){
    
    arrayMC =  (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
    if(!arrayMC) {
      printf("AliAnalysisTaskSELambdac::UserExec: MC particles branch not found!\n");
      return;
    }
    
  // load MC header
    mcHeader =  (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if(!mcHeader) {
    printf("AliAnalysisTaskSELambdac::UserExec: MC header branch not found!\n");
    return;
    }
  }
  
  Int_t n3Prong = array3Prong->GetEntriesFast();
  
  
  Int_t nOS=0;
  Int_t index;
  for (Int_t i3Prong = 0; i3Prong < n3Prong; i3Prong++) {
    AliAODRecoDecayHF3Prong *d = (AliAODRecoDecayHF3Prong*)array3Prong->UncheckedAt(i3Prong);
    
    
    Bool_t unsetvtx=kFALSE;
    if(!d->GetOwnPrimaryVtx()){
      d->SetOwnPrimaryVtx(vtx1);
      unsetvtx=kTRUE;
    }

    if(fRDCutsProduction->IsSelected(d,AliRDHFCuts::kCandidate)) {
      Int_t iPtBin = -1;
      Double_t ptCand = d->Pt();
      
      for(Int_t ibin=0;ibin<fNPtBins&&iPtBin<0&&ptCand>fArrayBinLimits[0]&&ptCand<fArrayBinLimits[fNPtBins];ibin++){
	if(ptCand<fArrayBinLimits[ibin+1])iPtBin=ibin;
      }

      Bool_t passTightCuts=fRDCutsAnalysis->IsSelected(d,AliRDHFCuts::kCandidate);
      Int_t labDp=-1;
      Float_t deltaPx=0.;
      Float_t deltaPy=0.;
      Float_t deltaPz=0.;
      Float_t truePt=0.;
      Float_t xDecay=0.;
      Float_t yDecay=0.;
      Float_t zDecay=0.;
      Float_t pdgCode=-2;
      if(fReadMC){
	labDp = MatchToMCLambdac(d,arrayMC);   
	//
	if(labDp>0){
	  AliAODMCParticle *partDp = (AliAODMCParticle*)arrayMC->At(labDp);
	  AliAODMCParticle *dg0 = (AliAODMCParticle*)arrayMC->At(partDp->GetDaughter(0));  
	  deltaPx=partDp->Px()-d->Px();
	  deltaPy=partDp->Py()-d->Py();
	  deltaPz=partDp->Pz()-d->Pz();
	  truePt=partDp->Pt();
	  xDecay=dg0->Xv();	  
	  yDecay=dg0->Yv();	  
	  zDecay=dg0->Zv();
	  pdgCode=TMath::Abs(partDp->GetPdgCode());
	}else{
	  pdgCode=-1;
	}
      }

      Double_t invMasspKpi=-1.;
      Double_t invMasspiKp=-1.;
      Int_t pdgs[3]={0,0,0};
      Double_t field=aod->GetMagneticField();
      if(fReadMC && fMCPid){
       
       if(IspKpiMC(d,arrayMC)) {
        invMasspKpi=d->InvMassLcpKpi();
	if(fUseKF){
	 pdgs[0]=2212;pdgs[1]=321;pdgs[2]=211;
	 if(!VertexingKF(d,pdgs,field)) invMasspKpi=-1.;
	}
       }
       IspiKpReal(d);
       IspiKpResonant(d,field);
       if(IspiKpMC(d,arrayMC)) {
        invMasspiKp=d->InvMassLcpiKp();
	if(fUseKF){
	 pdgs[0]=211;pdgs[1]=321;pdgs[2]=2212;
	 if(!VertexingKF(d,pdgs,field)) invMasspiKp=-1.;
	}
       }
      }
      if(fRealPid){
       if(IspKpiReal(d)) {
        invMasspKpi=d->InvMassLcpKpi();
	if(fUseKF){
	 pdgs[0]=2212;pdgs[1]=321;pdgs[2]=211;
	 if(!VertexingKF(d,pdgs,field)) invMasspKpi=-1.;
	}
       }
       if(IspiKpReal(d)) {
        invMasspiKp=d->InvMassLcpiKp();
	if(fUseKF){
	 pdgs[0]=211;pdgs[1]=321;pdgs[2]=2212;
	 if(!VertexingKF(d,pdgs,field)) invMasspiKp=-1.;
	}
       }
      }
      if(fResPid){
       if(IspKpiResonant(d,field)) {
        invMasspKpi=d->InvMassLcpKpi();
	if(fUseKF){
	 pdgs[0]=2212;pdgs[1]=321;pdgs[2]=211;
	 if(!VertexingKF(d,pdgs,field)) invMasspKpi=-1.;
	}
       }
       if(IspiKpResonant(d,field)) {
        invMasspiKp=d->InvMassLcpiKp();
	if(fUseKF){
	 pdgs[0]=2212;pdgs[1]=321;pdgs[2]=211;
	 if(!VertexingKF(d,pdgs,field)) invMasspiKp=-1.;
	}
       }
      }
      if(!fResPid && !fRealPid && !fMCPid){
       invMasspiKp=d->InvMassLcpiKp(); 
       if(fUseKF){
         pdgs[0]=211;pdgs[1]=321;pdgs[2]=2212;
         if(!VertexingKF(d,pdgs,field)) invMasspiKp=-1.;
        }
       invMasspKpi=d->InvMassLcpKpi();
       if(fUseKF){
        pdgs[0]=2212;pdgs[1]=321;pdgs[2]=211;
        if(!VertexingKF(d,pdgs,field)) invMasspKpi=-1.;
       }
      }

      if(invMasspiKp<0. && invMasspKpi<0.) continue;


      Float_t tmp[24];
      if(fFillNtuple){  	  
	tmp[0]=pdgCode;
	tmp[1]=deltaPx;
	tmp[2]=deltaPy;
	tmp[3]=deltaPz;
	tmp[4]=truePt;
	tmp[5]=xDecay;	  
	tmp[6]=yDecay;	  
	tmp[7]=zDecay;	  
	tmp[8]=d->PtProng(0);
	tmp[9]=d->PtProng(1);
	tmp[10]=d->PtProng(2);
	tmp[11]=d->Pt();
	tmp[12]=d->CosPointingAngle();
	tmp[13]=d->DecayLength();
	tmp[14]=d->Xv();
	tmp[15]=d->Yv();
	tmp[16]=d->Zv();
	if(invMasspiKp>0.) tmp[17]=invMasspiKp;
	if(invMasspKpi>0.) tmp[17]=invMasspKpi;
	tmp[18]=d->GetSigmaVert();
	tmp[19]=d->Getd0Prong(0);
	tmp[20]=d->Getd0Prong(1);
	tmp[21]=d->Getd0Prong(2);
	tmp[22]=d->GetDCA();
	tmp[23]=d->Prodd0d0(); 
	fNtupleLambdac->Fill(tmp);
	PostData(3,fNtupleLambdac);
      }
      Double_t dlen=d->DecayLength();
      Double_t cosp=d->CosPointingAngle();
      Double_t sumD02=d->Getd0Prong(0)*d->Getd0Prong(0)+d->Getd0Prong(1)*d->Getd0Prong(1)+d->Getd0Prong(2)*d->Getd0Prong(2);
      Double_t dca=d->GetDCA();      
Double_t ptmax=0;
      for(Int_t i=0;i<3;i++){
	if(d->PtProng(i)>ptmax)ptmax=d->PtProng(i);
      }
      if(iPtBin>=0){
      
	index=GetHistoIndex(iPtBin);
	if(invMasspiKp>0. && invMasspKpi>0.){
	if(invMasspiKp>0.) fMassHist[index]->Fill(invMasspiKp,0.5);
	if(invMasspKpi>0.) fMassHist[index]->Fill(invMasspKpi,0.5);
	}else{
         if(invMasspiKp>0.) fMassHist[index]->Fill(invMasspiKp);
         if(invMasspKpi>0.) fMassHist[index]->Fill(invMasspKpi);
        }
	fCosPHist[index]->Fill(cosp);
	fDLenHist[index]->Fill(dlen);
	fSumd02Hist[index]->Fill(sumD02);
	fPtMaxHist[index]->Fill(ptmax);
	fDCAHist[index]->Fill(dca);
	
	if(passTightCuts){
	 if(invMasspiKp>0. && invMasspKpi>0.){
	  if(invMasspiKp>0.) fMassHistTC[index]->Fill(invMasspiKp,0.5);
	  if(invMasspKpi>0.) fMassHistTC[index]->Fill(invMasspKpi,0.5);
	 }else{
          if(invMasspiKp>0.) fMassHistTC[index]->Fill(invMasspiKp);
          if(invMasspKpi>0.) fMassHistTC[index]->Fill(invMasspKpi);
         }
	}
	
	if(fReadMC){
	  if(labDp>0) {
	    index=GetSignalHistoIndex(iPtBin);
	    if(invMasspiKp>0. && invMasspKpi>0.){
	     if(invMasspiKp>0.) fMassHist[index]->Fill(invMasspiKp,0.5);
	     if(invMasspKpi>0.) fMassHist[index]->Fill(invMasspKpi,0.5);
            }else{
             if(invMasspiKp>0.) fMassHist[index]->Fill(invMasspiKp);
             if(invMasspKpi>0.) fMassHist[index]->Fill(invMasspKpi);
            }
	    fCosPHist[index]->Fill(cosp);
	    fDLenHist[index]->Fill(dlen);
	    fSumd02Hist[index]->Fill(sumD02);
	    fPtMaxHist[index]->Fill(ptmax);
	    fDCAHist[index]->Fill(dca);
	    if(passTightCuts){
	     if(invMasspiKp>0. && invMasspKpi>0.){
	      if(invMasspiKp>0.) fMassHistTC[index]->Fill(invMasspiKp,0.5);
	      if(invMasspKpi>0.) fMassHistTC[index]->Fill(invMasspKpi,0.5);
             }else{
              if(invMasspiKp>0.) fMassHistTC[index]->Fill(invMasspiKp);
              if(invMasspKpi>0.) fMassHistTC[index]->Fill(invMasspKpi);
             }
	    }
	    
	  }else{
	    index=GetBackgroundHistoIndex(iPtBin);
	    if(invMasspiKp>0. && invMasspKpi>0.){
	     fMassHist[index]->Fill(invMasspiKp,0.5);
	     fMassHist[index]->Fill(invMasspKpi,0.5);
	    }else{
	     if(invMasspiKp>0.) fMassHist[index]->Fill(invMasspiKp);
	     if(invMasspKpi>0.) fMassHist[index]->Fill(invMasspKpi);
	   }
	    fCosPHist[index]->Fill(cosp);
	    fDLenHist[index]->Fill(dlen);
	    fSumd02Hist[index]->Fill(sumD02);
	    fPtMaxHist[index]->Fill(ptmax);
	    fDCAHist[index]->Fill(dca);
	    if(passTightCuts){
	    if(invMasspiKp>0. && invMasspKpi>0.){
	       fMassHistTC[index]->Fill(invMasspiKp,0.5);
	       fMassHistTC[index]->Fill(invMasspKpi,0.5);
	      }else{
	       if(invMasspiKp>0.) fMassHistTC[index]->Fill(invMasspiKp);
	       if(invMasspKpi>0.) fMassHistTC[index]->Fill(invMasspKpi);
	      }
	    }	
	  }
	}
      }
      /*
      //start OS analysis
      if(labDp<0)fHistOSbkg->Fill(d->InvMassDplus());
      fHistOS->Fill(d->InvMassDplus());
      */
      nOS++;
    }
    if(unsetvtx) d->UnsetOwnPrimaryVtx();
  }
 
  PostData(1,fOutput);    
  return;
}



//________________________________________________________________________
void AliAnalysisTaskSELambdac::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
  if(fDebug > 1) printf("AnalysisTaskSELambdac: Terminate() \n");

  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {     
    printf("ERROR: fOutput not available\n");
    return;
  }
 fHistNEvents = dynamic_cast<TH1F*>(fOutput->FindObject("fHistNEvents"));

 TString hisname;
 Int_t index=0;
 

 //Int_t indexLS=0;
 for(Int_t i=0;i<fNPtBins;i++){
    index=GetHistoIndex(i);
    hisname.Form("hMassPt%d",i);
    fMassHist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
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
     hisname.Form("hDCAAllPt%d",i);
     fDCAHist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
     hisname.Form("hMassPt%dTC",i);
    fMassHistTC[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    
    index=GetSignalHistoIndex(i);    
    hisname.Form("hSigPt%d",i);
    fMassHist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
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
    hisname.Form("hDCASigPt%d",i);
    fDCAHist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    
    hisname.Form("hSigPt%dTC",i);
    fMassHistTC[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    
    index=GetBackgroundHistoIndex(i); 
    hisname.Form("hBkgPt%d",i);
    fMassHist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
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
    hisname.Form("hDCABkgPt%d",i);
    fDCAHist[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
    hisname.Form("hBkgPt%dTC",i);
    fMassHistTC[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
 
 }
    fTPCSignal3Sigma=dynamic_cast<TH2F*>(fOutput->FindObject(fTPCSignal3Sigma));
    fTPCSignal3SigmaReK=dynamic_cast<TH2F*>(fOutput->FindObject(fTPCSignal3SigmaReK));
    fTPCSignal3SigmaRep=dynamic_cast<TH2F*>(fOutput->FindObject(fTPCSignal3SigmaRep));
    fTPCSignal2Sigma=dynamic_cast<TH2F*>(fOutput->FindObject(fTPCSignal2Sigma));
    fTPCSignal2SigmaReK=dynamic_cast<TH2F*>(fOutput->FindObject(fTPCSignal2SigmaReK));
    fTPCSignal2SigmaRep=dynamic_cast<TH2F*>(fOutput->FindObject(fTPCSignal2SigmaRep));

  if(fFillNtuple){
    fNtupleLambdac = dynamic_cast<TNtuple*>(GetOutputData(3));
  }

  //TCanvas *c1=new TCanvas("c1","D+ invariant mass distribution",500,500);
  //c1->cd();
  //TH1F *hMassPt=(TH1F*)fOutput->FindObject("hMassPt3TC");
  //hMassPt->Draw();
 
 return;
}

//________________________________________________________________________
Int_t AliAnalysisTaskSELambdac::MatchToMCLambdac(AliAODRecoDecayHF3Prong *d,TClonesArray *arrayMC){

  Int_t lambdacLab[3]={0,0,0};
  Int_t pdgs[3]={0,0,0};
  for(Int_t i=0;i<3;i++){
   AliAODTrack *daugh=(AliAODTrack*)d->GetDaughter(i);
   Int_t lab=daugh->GetLabel();
   if(lab<0) return 0;
   AliAODMCParticle *part= (AliAODMCParticle*)arrayMC->At(lab);
   if(!part) continue;
   pdgs[i]=part->GetPdgCode();
   Int_t partPdgcode = TMath::Abs(part->GetPdgCode());
   if(partPdgcode==211 || partPdgcode==321 || partPdgcode==2212){
        Int_t motherLabel=part->GetMother();
	     if(motherLabel<0) return 0;
   AliAODMCParticle *motherPart = (AliAODMCParticle*)arrayMC->At(motherLabel);
   if(!motherPart) continue;
        Int_t motherPdg = TMath::Abs(motherPart->GetPdgCode());
   if(motherPdg==4122) {
        if(GetLambdacDaugh(motherPart,arrayMC)){lambdacLab[i]=motherLabel;continue;}
	      }
   if(motherPdg==313 || motherPdg==2224 || motherPdg==3124){
         Int_t GmotherLabel=motherPart->GetMother();
	       if(GmotherLabel<0) return 0;
	 AliAODMCParticle *GmotherPart = (AliAODMCParticle*)arrayMC->At(GmotherLabel);
	 if(!GmotherPart) continue;
	       Int_t GmotherPdg = TMath::Abs(GmotherPart->GetPdgCode());
	 if(GmotherPdg==4122) {
	        if(GetLambdacDaugh(GmotherPart,arrayMC)) {lambdacLab[i]=GmotherLabel;continue;}
	      }
	}
     }
  }

 if(lambdacLab[0]==lambdacLab[1] && lambdacLab[1]==lambdacLab[2]) {return lambdacLab[0];}
 return 0;

}
//------------------------
Bool_t AliAnalysisTaskSELambdac::GetLambdacDaugh(AliAODMCParticle *part,TClonesArray *arrayMC){

 Int_t numberOfLambdac=0;
 if(TMath::Abs(part->GetPdgCode())!=4122) return kFALSE;
 Int_t daugh_tmp[2];
 daugh_tmp[0]=part->GetDaughter(0);
 daugh_tmp[1]=part->GetDaughter(1);
 Int_t nDaugh = (Int_t)part->GetNDaughters();
 if(nDaugh<2) return kFALSE;
 if(nDaugh>3) return kFALSE;
 AliAODMCParticle* pdaugh1 = (AliAODMCParticle*)arrayMC->At(part->GetDaughter(0));
   if(!pdaugh1) {return kFALSE;}
 Int_t number1 = TMath::Abs(pdaugh1->GetPdgCode());
 AliAODMCParticle* pdaugh2 = (AliAODMCParticle*)arrayMC->At(part->GetDaughter(1));
   if(!pdaugh2) {return kFALSE;}
  Int_t number2 = TMath::Abs(pdaugh2->GetPdgCode());

  if(nDaugh==3){
   Int_t thirdDaugh=part->GetDaughter(1)-1;
            AliAODMCParticle* pdaugh3 = (AliAODMCParticle*)arrayMC->At(thirdDaugh);
  Int_t number3 = TMath::Abs(pdaugh3->GetPdgCode());
   if((number1==321 && number2==211 && number3==2212) || (number1==211 && number2==321 && number3==2212) || (number1==211 && number2==2212 && number3==321) || (number1==321 && number2==2212 && number3==211) || (number1==2212 && number2==321 && number3==211) || (number1==2212 && number2==211 && number3==321)) {
   numberOfLambdac++;
   } 
  }

 if(nDaugh==2){

         //Lambda in Lambda pi
//  if((number1==211 && number2==3122)|| (number1==3122 && number2==211))return kFALSE;
  
  //Lambda resonant
  
  //Lambda -> p K*0
  //
  Int_t nfiglieK=0;

   if((number1==2212 && number2==313)){
     nfiglieK=pdaugh2->GetNDaughters();
     if(nfiglieK!=2) return kFALSE;
     AliAODMCParticle* pdaughK1 = (AliAODMCParticle*)arrayMC->At(pdaugh2->GetDaughter(0));
     AliAODMCParticle* pdaughK2 = (AliAODMCParticle*)arrayMC->At(pdaugh2->GetDaughter(1));
     if(!pdaughK1) return kFALSE;
     if(!pdaughK2) return kFALSE;
     if((TMath::Abs(pdaughK1->GetPdgCode())==211 && TMath::Abs(pdaughK2->GetPdgCode())==321) || (TMath::Abs(pdaughK1->GetPdgCode())==321 && TMath::Abs(pdaughK2->GetPdgCode())==211)) numberOfLambdac++;
    }

   if((number1==313 && number2==2212)){
    nfiglieK=pdaugh1->GetNDaughters();
    if(nfiglieK!=2) return kFALSE;
    AliAODMCParticle* pdaughK1 = (AliAODMCParticle*)arrayMC->At(pdaugh1->GetDaughter(0));
    AliAODMCParticle* pdaughK2 = (AliAODMCParticle*)arrayMC->At(pdaugh1->GetDaughter(1));
    if(!pdaughK1) return kFALSE;
    if(!pdaughK2) return kFALSE;
     if((TMath::Abs(pdaughK1->GetPdgCode())==211 && TMath::Abs(pdaughK2->GetPdgCode())==321) || (TMath::Abs(pdaughK1->GetPdgCode())==321 && TMath::Abs(pdaughK2->GetPdgCode())==211)) numberOfLambdac++;
   }

   //Lambda -> Delta++ k
   Int_t nfiglieDelta=0;
   if(number1==321 && number2==2224){
    nfiglieDelta=pdaugh2->GetNDaughters();
    if(nfiglieDelta!=2) return kFALSE;
    AliAODMCParticle *pdaughD1=(AliAODMCParticle*)arrayMC->At(pdaugh2->GetDaughter(0));
    AliAODMCParticle *pdaughD2=(AliAODMCParticle*)arrayMC->At(pdaugh2->GetDaughter(1));
    if(!pdaughD1) return kFALSE;
    if(!pdaughD2) return kFALSE;
    if((TMath::Abs(pdaughD1->GetPdgCode())==211 && TMath::Abs(pdaughD2->GetPdgCode())==2212) || (TMath::Abs(pdaughD1->GetPdgCode())==2212 && TMath::Abs(pdaughD2->GetPdgCode())==211)) numberOfLambdac++;
   }
   if(number1==2224 && number2==321){
    nfiglieDelta=pdaugh1->GetNDaughters();
    if(nfiglieDelta!=2) return kFALSE;
    AliAODMCParticle* pdaughD1 = (AliAODMCParticle*)arrayMC->At(pdaugh1->GetDaughter(0));
    AliAODMCParticle* pdaughD2 = (AliAODMCParticle*)arrayMC->At(pdaugh1->GetDaughter(1));
    if(!pdaughD1) return kFALSE;
    if(!pdaughD2) return kFALSE;
    if((TMath::Abs(pdaughD1->GetPdgCode())==211 && TMath::Abs(pdaughD2->GetPdgCode())==2212) || (TMath::Abs(pdaughD1->GetPdgCode())==2212 && TMath::Abs(pdaughD2->GetPdgCode())==211)) numberOfLambdac++;
   }
    

   //Lambdac -> Lambda(1520) pi
   Int_t nfiglieLa=0;
   if(number1==3124 && number2==211){
    nfiglieLa=pdaugh1->GetNDaughters();
    if(nfiglieLa!=2) return kFALSE;
    AliAODMCParticle *pdaughL1=(AliAODMCParticle*)arrayMC->At(pdaugh1->GetDaughter(0));
    AliAODMCParticle *pdaughL2=(AliAODMCParticle*)arrayMC->At(pdaugh1->GetDaughter(1));
    if(!pdaughL1) return kFALSE;
    if(!pdaughL2) return kFALSE;
    if((TMath::Abs(pdaughL1->GetPdgCode())==321 && TMath::Abs(pdaughL2->GetPdgCode())==2212) || (TMath::Abs(pdaughL1->GetPdgCode())==2212 && TMath::Abs(pdaughL2->GetPdgCode())==321)) numberOfLambdac++;
   }
   if(number1==211 && number2==3124){
    nfiglieLa=pdaugh2->GetNDaughters();
    if(nfiglieLa!=2) return kFALSE;
    AliAODMCParticle *pdaughL1=(AliAODMCParticle*)arrayMC->At(pdaugh2->GetDaughter(0));
    AliAODMCParticle *pdaughL2=(AliAODMCParticle*)arrayMC->At(pdaugh2->GetDaughter(1));
    if(!pdaughL1) return kFALSE;
    if(!pdaughL2) return kFALSE;
    if((TMath::Abs(pdaughL1->GetPdgCode())==321 && TMath::Abs(pdaughL2->GetPdgCode())==2212) || (TMath::Abs(pdaughL1->GetPdgCode())==2212 && TMath::Abs(pdaughL2->GetPdgCode())==321)) numberOfLambdac++;
   
   }
 }

 if(numberOfLambdac>0) {return kTRUE;}
         return kFALSE;
}
//-----------------------------
Bool_t AliAnalysisTaskSELambdac::IspKpiMC(AliAODRecoDecayHF3Prong *d,TClonesArray *arrayMC){

   Int_t lab[3]={0,0,0},pdgs[3]={0,0,0};
   for(Int_t i=0;i<3;i++){
    AliAODTrack *daugh=(AliAODTrack*)d->GetDaughter(i);
    lab[i]=daugh->GetLabel();
    if(lab[i]<0) return kFALSE;
    AliAODMCParticle *part= (AliAODMCParticle*)arrayMC->At(lab[i]);
    if(!part) return kFALSE;
    pdgs[i]=TMath::Abs(part->GetPdgCode());
   }

   if(pdgs[0]==2212 && pdgs[1]==321 && pdgs[2]==211) return kTRUE;

   return kFALSE;
}
//-----------------------------
Bool_t AliAnalysisTaskSELambdac::IspiKpMC(AliAODRecoDecayHF3Prong *d,TClonesArray *arrayMC){

   Int_t lab[3]={0,0,0},pdgs[3]={0,0,0};
   for(Int_t i=0;i<3;i++){
    AliAODTrack *daugh=(AliAODTrack*)d->GetDaughter(i);
    lab[i]=daugh->GetLabel();
    if(lab[i]<0) return kFALSE;
    AliAODMCParticle *part= (AliAODMCParticle*)arrayMC->At(lab[i]);
    if(!part) return kFALSE;
    pdgs[i]=TMath::Abs(part->GetPdgCode());
   }

   if(pdgs[2]==2212 && pdgs[1]==321 && pdgs[0]==211) {return kTRUE;}

   return kFALSE;
}
//--------------------------------------
Bool_t AliAnalysisTaskSELambdac::IspiKpReal(AliAODRecoDecayHF3Prong *d){

  AliAODPidHF* pidObj=new AliAODPidHF();
  Bool_t type[2]={kFALSE,kFALSE};
  for(Int_t i=0;i<3;i++){
//   Bool_t pid[3]={kFALSE,kFALSE,kFALSE};
   AliAODTrack *track=(AliAODTrack*)d->GetDaughter(i);
   AliAODPid *pidObjtrk=track->GetDetPid();
//   pidObj->CombinedProbability(track,pid);

//   Int_t i2=i-1;
//   type[i2]=pid[i];
   Double_t dedx=pidObjtrk->GetTPCsignal();
   Double_t mom = pidObjtrk->GetTPCmomentum();
   //printf("TPC mom= %f\n",mom);
    pidObj->SetSigma(3.);
   if(pidObj->IsProtonRaw(track,"TPC")){ 
    fTPCSignal3Sigma->Fill(mom,dedx);
 //   if(i==2) type[1]=kTRUE;
    pidObj->SetSigma(1.);
    //if(!pidObj->IsPionRaw(track,"TPC") && !pidObj->IsKaonRaw(track,"TPC"))
    if(!pidObj->IsKaonRaw(track,"TPC")){
     fTPCSignal3SigmaRep->Fill(mom,dedx);
    }
   }
    pidObj->SetSigma(3.);
   if(pidObj->IsKaonRaw(track,"TPC")) {
    fTPCSignal3Sigma->Fill(mom,dedx);
 //   if(i==1) type[0]=kTRUE;
    pidObj->SetSigma(1.);
    //if(!pidObj->IsPionRaw(track,"TPC") && !pidObj->IsProtonRaw(track,"TPC"))
    if(!pidObj->IsPionRaw(track,"TPC")){
     fTPCSignal3SigmaReK->Fill(mom,dedx);
    }
   }

    pidObj->SetSigma(2.);
   if(pidObj->IsProtonRaw(track,"TPC")){ 
    if(i==2) type[1]=kTRUE;
    fTPCSignal2Sigma->Fill(mom,dedx);
    pidObj->SetSigma(1.);
    //if(!pidObj->IsPionRaw(track,"TPC") && !pidObj->IsKaonRaw(track,"TPC"))
    if(!pidObj->IsKaonRaw(track,"TPC")){
     fTPCSignal2SigmaRep->Fill(mom,dedx);
    }
   }
   if(pidObj->IsKaonRaw(track,"TPC")) {
    if(i==1) type[0]=kTRUE;
    fTPCSignal2Sigma->Fill(mom,dedx);
    pidObj->SetSigma(1.);
    //if(!pidObj->IsPionRaw(track,"TPC") && !pidObj->IsProtonRaw(track,"TPC"))
    if(!pidObj->IsPionRaw(track,"TPC")){
     fTPCSignal2SigmaReK->Fill(mom,dedx);
    }
   }
  }

  if(type[0] && type[1]) {delete pidObj;return kTRUE;}
  delete pidObj;
 return kFALSE;
}
//--------------------------------------
Bool_t AliAnalysisTaskSELambdac::IspKpiReal(AliAODRecoDecayHF3Prong *d){

  AliAODPidHF* pidObj=new AliAODPidHF();
  Bool_t type[2]={kFALSE,kFALSE};
    pidObj->SetSigma(2.);
  for(Int_t i=0;i<2;i++){
   //Bool_t pid[3]={kFALSE,kFALSE,kFALSE};
   AliAODTrack *track=(AliAODTrack*)d->GetDaughter(i);
  // pidObj->CombinedProbability(track,pid);
  // if(i==0) type[i]=pid[2];
  // if(i==1) type[i]=pid[i];
   if(pidObj->IsProtonRaw(track,"TPC") && i==0) type[1]=kTRUE;
   if(pidObj->IsKaonRaw(track,"TPC") && i==1) type[0]=kTRUE;
  }
  if(type[0] && type[1]) {delete pidObj;return kTRUE;}
  delete pidObj;
 return kFALSE;
}
//------------------------------------
Bool_t AliAnalysisTaskSELambdac::VertexingKF(AliAODRecoDecayHF3Prong *d,Int_t *pdgs,Double_t field){
  
   Int_t iprongs[3]={0,1,2};
   Double_t mass[2]={0.,0.};
 //topological constr
   AliKFParticle *Lambdac=d->ApplyVertexingKF(iprongs,3,pdgs,kTRUE,field,mass);
   if(!Lambdac) return kFALSE;
  Double_t probTot=TMath::Prob(Lambdac->GetChi2(),Lambdac->GetNDF());
  if(probTot<fCutsKF[0]) return kFALSE;

  if(d->PtProng(0)<1.5 && d->Getd0Prong(0)>0.02) return kFALSE;
    if(d->PtProng(1)<1.5 && d->Getd0Prong(1)>0.02) return kFALSE;
      if(d->PtProng(2)<1.5 && d->Getd0Prong(2)>0.02) return kFALSE;
 //mass constr for K*
   Int_t ipRes[2];
   Int_t pdgres[2];
   mass[0]=0.8961;mass[1]=0.03;
   if(TMath::Abs(pdgs[0])==211){
    ipRes[0]=0;ipRes[1]=1;
    pdgres[0]=pdgs[0];pdgres[1]=321;
   }
   if(TMath::Abs(pdgs[2])==211){
    ipRes[0]=2;ipRes[1]=1;
    pdgres[0]=pdgs[2];pdgres[1]=321;
   }
   AliKFParticle *KappaStar=d->ApplyVertexingKF(ipRes,2,pdgres,kFALSE,field,mass);

  Double_t probKstar=TMath::Prob(KappaStar->GetChi2(),KappaStar->GetNDF());
  if(probKstar>fCutsKF[1]) {
    AliAODTrack *esdProng1=(AliAODTrack*)d->GetDaughter(ipRes[0]);
    AliAODTrack *esdProng2=(AliAODTrack*)d->GetDaughter(ipRes[1]);
    AliKFParticle prong1(*esdProng1,pdgres[0]);
    AliKFParticle prong2(*esdProng2,pdgres[1]);
    if(KappaStar->GetPt()<fCutsKF[2] && prong1.GetAngle(prong2)>fCutsKF[3]) return kFALSE;
  } 
 //mass constr for Lambda
   mass[0]=1.520;mass[1]=0.005;
   if(TMath::Abs(pdgs[0])==2212){
    ipRes[0]=0;ipRes[1]=1;
    pdgres[0]=pdgs[0];pdgres[1]=pdgs[1];
   }
   if(TMath::Abs(pdgs[2])==2212){
    ipRes[0]=2;ipRes[1]=1;
    pdgres[0]=pdgs[2];pdgres[1]=pdgs[1];
   }
   AliKFParticle *Lambda1520=d->ApplyVertexingKF(ipRes,2,pdgres,kFALSE,field,mass);
  Double_t probLa=TMath::Prob(Lambda1520->GetChi2(),Lambda1520->GetNDF());
  if(probLa>fCutsKF[4]) {
    AliAODTrack *esdProng1=(AliAODTrack*)d->GetDaughter(ipRes[0]);
    AliAODTrack *esdProng2=(AliAODTrack*)d->GetDaughter(ipRes[1]);
    AliKFParticle prong1(*esdProng1,pdgres[0]);
    AliKFParticle prong2(*esdProng2,pdgres[1]);
    if(Lambda1520->GetPt()<fCutsKF[5] && prong1.GetAngle(prong2)>fCutsKF[6]) return kFALSE;
  } 
 //mass constr for Delta
   mass[0]=1.2;mass[1]=0.15;
   ipRes[0]=0;ipRes[1]=2;
   pdgres[0]=pdgs[0];pdgres[2]=pdgs[2];
   AliKFParticle *Delta=d->ApplyVertexingKF(ipRes,2,pdgres,kFALSE,field,mass);
  Double_t probDelta=TMath::Prob(Delta->GetChi2(),Delta->GetNDF());
  if(probDelta>fCutsKF[7]) {
    AliAODTrack *esdProng1=(AliAODTrack*)d->GetDaughter(ipRes[0]);
    AliAODTrack *esdProng2=(AliAODTrack*)d->GetDaughter(ipRes[1]);
    AliKFParticle prong1(*esdProng1,pdgres[0]);
    AliKFParticle prong2(*esdProng2,pdgres[1]);
    if(Delta->GetPt()<fCutsKF[8] && prong1.GetAngle(prong2)>fCutsKF[9]) return kFALSE;
  } 
 return kTRUE;
}
//-------------------------------------
Bool_t AliAnalysisTaskSELambdac::IspiKpResonant(AliAODRecoDecayHF3Prong *d,Double_t field){
  
 //if lambda* -> pk
        Double_t mass[2]={1.520,0.005};
        Int_t ipRes[2]={1,2};
        Int_t pdgres[2]={321,2212};
        AliKFParticle *Lambda1520=d->ApplyVertexingKF(ipRes,2,pdgres,kFALSE,field,mass);
	Double_t probLa=TMath::Prob(Lambda1520->GetChi2(),Lambda1520->GetNDF());
	if(probLa>0.1) return kTRUE;
 //K* -> kpi
//        mass[0]=0.8961;mass[1]=0.03;
//        ipRes[0]=0;ipRes[1]=1;
//        pdgres[0]=211;pdgres[1]=321;
//        AliKFParticle *Kstar=d->ApplyVertexingKF(ipRes,2,pdgres,kFALSE,field,mass);
//	Double_t probKa=TMath::Prob(Kstar->GetChi2(),Kstar->GetNDF());
//	if(probKa>0.1) return kTRUE;

 return kFALSE;

}
//-------------------------------------
Bool_t AliAnalysisTaskSELambdac::IspKpiResonant(AliAODRecoDecayHF3Prong *d,Double_t field){
  
 //if lambda* -> pk
        Double_t mass[2]={1.520,0.005};
        Int_t ipRes[2]={0,1};
        Int_t pdgres[2]={2212,321};
        AliKFParticle *Lambda1520=d->ApplyVertexingKF(ipRes,2,pdgres,kFALSE,field,mass);
	Double_t probLa=TMath::Prob(Lambda1520->GetChi2(),Lambda1520->GetNDF());
	if(probLa>0.1) return kTRUE;
 //K* -> kpi
//        mass[0]=0.8961;mass[1]=0.03;
//        ipRes[0]=1;ipRes[1]=2;
//        pdgres[1]=211;pdgres[0]=321;
//        AliKFParticle *Kstar=d->ApplyVertexingKF(ipRes,2,pdgres,kFALSE,field,mass);
//	Double_t probKa=TMath::Prob(Kstar->GetChi2(),Kstar->GetNDF());
//	if(probKa>0.1) return kTRUE;

 return kFALSE;

}
