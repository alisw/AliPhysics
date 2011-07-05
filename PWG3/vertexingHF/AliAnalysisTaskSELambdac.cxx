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
// AliAnalysisTaskSE for the extraction of signal(e.g Lambdac) of heavy flavor
// decay candidates with the MC truth.
// Authors: r.romita@gsi.de
/////////////////////////////////////////////////////////////

#include <TClonesArray.h>
#include <TNtuple.h>
#include <TCanvas.h>
#include <TList.h>
#include <TString.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TDatabasePDG.h>

#include <AliAnalysisDataSlot.h>
#include <AliAnalysisDataContainer.h>
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
#include "AliKFVertex.h"
#include "AliESDVertex.h"
#include "AliTOFPIDResponse.h"
#include "AliAODpidUtil.h"
#include "AliAODPid.h"
#include "AliInputEventHandler.h"

ClassImp(AliAnalysisTaskSELambdac)


//________________________________________________________________________
AliAnalysisTaskSELambdac::AliAnalysisTaskSELambdac():
AliAnalysisTaskSE(),
fOutput(0), 
fHistNEvents(0),
fhChi2(0),
fhMassPtGreater3(0),
fhMassPtGreater3TC(0),
fhMassPtGreater3Kp(0),
fhMassPtGreater3KpTC(0),
fhMassPtGreater3Lpi(0),
fhMassPtGreater3LpiTC(0),
fhMassPtGreater2(0),
fhMassPtGreater2TC(0),
fhMassPtGreater2Kp(0),
fhMassPtGreater2KpTC(0),
fhMassPtGreater2Lpi(0),
fhMassPtGreater2LpiTC(0),
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
fAnalysis(kFALSE),
fVHF(0),
fFillVarHists(kTRUE),
fNentries(0),
fOutputMC(0),
fUtilPid(0)
{
   // Default constructor
}

//________________________________________________________________________
AliAnalysisTaskSELambdac::AliAnalysisTaskSELambdac(const char *name,Bool_t fillNtuple,AliRDHFCutsLctopKpi *lccutsana,AliRDHFCutsLctopKpi *lccutsprod):
AliAnalysisTaskSE(name),
fOutput(0),
fHistNEvents(0),
fhChi2(0),
fhMassPtGreater3(0),
fhMassPtGreater3TC(0),
fhMassPtGreater3Kp(0),
fhMassPtGreater3KpTC(0),
fhMassPtGreater3Lpi(0),
fhMassPtGreater3LpiTC(0),
fhMassPtGreater2(0),
fhMassPtGreater2TC(0),
fhMassPtGreater2Kp(0),
fhMassPtGreater2KpTC(0),
fhMassPtGreater2Lpi(0),
fhMassPtGreater2LpiTC(0),
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
fAnalysis(kFALSE),
fVHF(0),
fFillVarHists(kTRUE),
fNentries(0),
fOutputMC(0),
fUtilPid(0)
{
   SetPtBinLimit(fRDCutsAnalysis->GetNPtBins()+1,fRDCutsAnalysis->GetPtBinLimits());
  // Default constructor
   // Output slot #1 writes into a TList container
  DefineOutput(1,TList::Class());  //My private output
  DefineOutput(2,TList::Class());
  DefineOutput(3,TList::Class());
  DefineOutput(4,TH1F::Class());
  if(fFillNtuple){
    // Output slot #2 writes into a TNtuple container
    DefineOutput(5,TNtuple::Class());  //My private output
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
  if (fOutputMC) {
    delete fOutputMC;
    fOutputMC = 0;
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
 if (fNentries){
   delete fNentries;
   fNentries = 0;
   }
 if (fUtilPid){
   delete fUtilPid;
   fUtilPid = 0;
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
    fArrayBinLimits[3]=4.;
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
    for(Int_t i=0; i<fNPtBins; i++) printf(" Bin%d = %8.2f-%8.2f\n",i,fArrayBinLimits[i],fArrayBinLimits[i+1]);    
  }
}
//_________________________________________________________________
Double_t  AliAnalysisTaskSELambdac::GetPtBinLimit(Int_t ibin) const{
  if(ibin>fNPtBins)return -1;
  return fArrayBinLimits[ibin];
} 

//_________________________________________________________________
void AliAnalysisTaskSELambdac::Init()
{
  // Initialization

  if(fDebug > 1) printf("AnalysisTaskSELambdac::Init() \n");

  fListCuts=new TList();
  fListCuts->SetOwner();

  fListCuts->Add(new AliRDHFCutsLctopKpi(*fRDCutsAnalysis));
  fListCuts->Add(new AliRDHFCutsLctopKpi(*fRDCutsProduction));
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
    hisname.Form("hMassPt%dTC",i);
    fMassHistTC[index]=new TH1F(hisname.Data(),hisname.Data(),100,fLowmasslimit,fUpmasslimit);
    fMassHistTC[index]->Sumw2();

    hisname.Form("hMassPtLpi%d",i);
    fMassHistLpi[index]=new TH1F(hisname.Data(),hisname.Data(),100,fLowmasslimit,fUpmasslimit);
    fMassHistLpi[index]->Sumw2();
    hisname.Form("hMassPtLpi%dTC",i);
    fMassHistLpiTC[index]=new TH1F(hisname.Data(),hisname.Data(),100,fLowmasslimit,fUpmasslimit);
    fMassHistLpiTC[index]->Sumw2();

    hisname.Form("hMassPtKp%d",i);
    fMassHistKp[index]=new TH1F(hisname.Data(),hisname.Data(),100,fLowmasslimit,fUpmasslimit);
    fMassHistKp[index]->Sumw2();
    hisname.Form("hMassPtKp%dTC",i);
    fMassHistKpTC[index]=new TH1F(hisname.Data(),hisname.Data(),100,fLowmasslimit,fUpmasslimit);
    fMassHistKpTC[index]->Sumw2();
//signal
    index=GetSignalHistoIndex(i);    
    hisname.Form("hSigPt%d",i);
    fMassHist[index]=new TH1F(hisname.Data(),hisname.Data(),100,fLowmasslimit,fUpmasslimit);
    fMassHist[index]->Sumw2();
    hisname.Form("hSigPt%dTC",i);
    fMassHistTC[index]=new TH1F(hisname.Data(),hisname.Data(),100,fLowmasslimit,fUpmasslimit);
    fMassHistTC[index]->Sumw2();
    hisname.Form("hSigPtLpi%d",i);
    fMassHistLpi[index]=new TH1F(hisname.Data(),hisname.Data(),100,fLowmasslimit,fUpmasslimit);
    fMassHistLpi[index]->Sumw2();
    hisname.Form("hSigPtLpi%dTC",i);
    fMassHistLpiTC[index]=new TH1F(hisname.Data(),hisname.Data(),100,fLowmasslimit,fUpmasslimit);
    fMassHistLpiTC[index]->Sumw2();

    hisname.Form("hSigPtKp%d",i);
    fMassHistKp[index]=new TH1F(hisname.Data(),hisname.Data(),100,fLowmasslimit,fUpmasslimit);
    fMassHistKp[index]->Sumw2();
    hisname.Form("hSigPtKp%dTC",i);
    fMassHistKpTC[index]=new TH1F(hisname.Data(),hisname.Data(),100,fLowmasslimit,fUpmasslimit);
    fMassHistKpTC[index]->Sumw2();

    index=GetBackgroundHistoIndex(i); 
    hisname.Form("hBkgPt%d",i);
    fMassHist[index]=new TH1F(hisname.Data(),hisname.Data(),100,fLowmasslimit,fUpmasslimit);
    fMassHist[index]->Sumw2();
    hisname.Form("hBkgPt%dTC",i);
    fMassHistTC[index]=new TH1F(hisname.Data(),hisname.Data(),100,fLowmasslimit,fUpmasslimit);
    fMassHistTC[index]->Sumw2();
    hisname.Form("hBkgPtLpi%d",i);
    fMassHistLpi[index]=new TH1F(hisname.Data(),hisname.Data(),100,fLowmasslimit,fUpmasslimit);
    fMassHistLpi[index]->Sumw2();
    hisname.Form("hBkgPtLpi%dTC",i);
    fMassHistLpiTC[index]=new TH1F(hisname.Data(),hisname.Data(),100,fLowmasslimit,fUpmasslimit);
    fMassHistLpiTC[index]->Sumw2();

    hisname.Form("hBkgPtKp%d",i);
    fMassHistKp[index]=new TH1F(hisname.Data(),hisname.Data(),100,fLowmasslimit,fUpmasslimit);
    fMassHistKp[index]->Sumw2();
    hisname.Form("hBkgPtKp%dTC",i);
    fMassHistKpTC[index]=new TH1F(hisname.Data(),hisname.Data(),100,fLowmasslimit,fUpmasslimit);
    fMassHistKpTC[index]->Sumw2();
}
  
  for(Int_t i=0; i<3*fNPtBins; i++){
    fOutput->Add(fMassHist[i]);
    fOutput->Add(fMassHistTC[i]);
    fOutput->Add(fMassHistLpi[i]);
    fOutput->Add(fMassHistLpiTC[i]);
    fOutput->Add(fMassHistKp[i]);
    fOutput->Add(fMassHistKpTC[i]);
  }

  fHistNEvents = new TH1F("fHistNEvents", "Number of processed events; ; Events",3,-1.5,1.5);
  fHistNEvents->Sumw2();
  fHistNEvents->SetMinimum(0);
  fOutput->Add(fHistNEvents);

  fhChi2 = new TH1F("fhChi2", "Chi2",100,0.,10.);
  fhChi2->Sumw2();
  fOutput->Add(fhChi2);

  fhMassPtGreater3=new TH1F("fhMassPtGreater3","Pt > 3 GeV/c",100,fLowmasslimit,fUpmasslimit);
  fhMassPtGreater3->Sumw2();
  fOutput->Add(fhMassPtGreater3);
  fhMassPtGreater3TC=new TH1F("fhMassPtGreater3TC","Pt > 3 GeV/c",100,fLowmasslimit,fUpmasslimit);
  fhMassPtGreater3TC->Sumw2();
  fOutput->Add(fhMassPtGreater3TC);
  fhMassPtGreater3Kp=new TH1F("fhMassPtGreater3Kp","Pt > 3 GeV/c",100,fLowmasslimit,fUpmasslimit);
  fhMassPtGreater3Kp->Sumw2();
  fOutput->Add(fhMassPtGreater3Kp);
  fhMassPtGreater3KpTC=new TH1F("fhMassPtGreater3KpTC","Pt > 3 GeV/c",100,fLowmasslimit,fUpmasslimit);
  fhMassPtGreater3KpTC->Sumw2();
  fOutput->Add(fhMassPtGreater3KpTC);
  fhMassPtGreater3Lpi=new TH1F("fhMassPtGreater3Lpi","Pt > 3 GeV/c",100,fLowmasslimit,fUpmasslimit);
  fhMassPtGreater3Lpi->Sumw2();
  fOutput->Add(fhMassPtGreater3Lpi);
  fhMassPtGreater3LpiTC=new TH1F("fhMassPtGreater3LpiTC","Pt > 3 GeV/c",100,fLowmasslimit,fUpmasslimit);
  fhMassPtGreater3LpiTC->Sumw2();
  fOutput->Add(fhMassPtGreater3LpiTC);
  fhMassPtGreater2=new TH1F("fhMassPtGreater2","Pt > 2 GeV/c",100,fLowmasslimit,fUpmasslimit);
  fhMassPtGreater2->Sumw2();
  fOutput->Add(fhMassPtGreater2);
  fhMassPtGreater2TC=new TH1F("fhMassPtGreater2TC","Pt > 2 GeV/c",100,fLowmasslimit,fUpmasslimit);
  fhMassPtGreater2TC->Sumw2();
  fOutput->Add(fhMassPtGreater2TC);
  fhMassPtGreater2Kp=new TH1F("fhMassPtGreater2Kp","Pt > 2 GeV/c",100,fLowmasslimit,fUpmasslimit);
  fhMassPtGreater2Kp->Sumw2();
  fOutput->Add(fhMassPtGreater2Kp);
  fhMassPtGreater2KpTC=new TH1F("fhMassPtGreater2KpTC","Pt > 2 GeV/c",100,fLowmasslimit,fUpmasslimit);
  fhMassPtGreater2KpTC->Sumw2();
  fOutput->Add(fhMassPtGreater2KpTC);
  fhMassPtGreater2Lpi=new TH1F("fhMassPtGreater2Lpi","Pt > 2 GeV/c",100,fLowmasslimit,fUpmasslimit);
  fhMassPtGreater2Lpi->Sumw2();
  fOutput->Add(fhMassPtGreater2Lpi);
  fhMassPtGreater2LpiTC=new TH1F("fhMassPtGreater2LpiTC","Pt > 2 GeV/c",100,fLowmasslimit,fUpmasslimit);
  fhMassPtGreater2LpiTC->Sumw2();
  fOutput->Add(fhMassPtGreater2LpiTC);
  
  fOutputMC = new TList();
  fOutputMC->SetOwner();
  fOutputMC->SetName("QAMCHistos");

//  const char* nameoutput=GetOutputSlot(4)->GetContainer()->GetName();

  fNentries=new TH1F("fNentries", "Integral(1,2) = number of AODs *** Integral(2,3) = number of candidates selected with cuts *** Integral(3,4) = number of Lc selected with cuts *** Integral(4,5) = events with good vertex ***  Integral(5,6) = pt out of bounds", 11,-0.5,10.5);

  //ROS: qui il bin assignment e' modellato su D0 ma sicuramente iv arie
  fNentries->GetXaxis()->SetBinLabel(1,"nEventsAnal");
  fNentries->GetXaxis()->SetBinLabel(2,"nCandSel(Cuts)");
  fNentries->GetXaxis()->SetBinLabel(3,"nLcSelected");
  fNentries->GetXaxis()->SetBinLabel(4,"nEventsGoodVtxS");
  fNentries->GetXaxis()->SetBinLabel(5,"ptbin = -1");
  fNentries->GetXaxis()->SetBinLabel(6,"no daughter");
  fNentries->GetXaxis()->SetBinLabel(7,"nCandSel(Tr)");
  fNentries->GetXaxis()->SetBinLabel(8,"PID=0");
  fNentries->GetXaxis()->SetBinLabel(9,"PID=1");
  fNentries->GetXaxis()->SetBinLabel(10,"PID=2");
  fNentries->GetXaxis()->SetBinLabel(11,"PID=3");
  fNentries->GetXaxis()->SetNdivisions(1,kFALSE);

  hisname.Form("hMass");
  TH1F *hMassInv=new TH1F(hisname.Data(),hisname.Data(),100,fLowmasslimit,fUpmasslimit);
  fOutputMC->Add(hMassInv);
  hisname.Form("hbMass");
  TH1F *hBMassInv=new TH1F(hisname.Data(),hisname.Data(),100,fLowmasslimit,fUpmasslimit);
  fOutputMC->Add(hBMassInv);

  // proton specific
  hisname.Form("hpTOFSignal");
  TH1F *hProtonTOFSignal=new TH1F(hisname.Data(),hisname.Data(),100,12000.,50000.0);
  fOutputMC->Add(hProtonTOFSignal);
  hisname.Form("hbpTOFSignal");
  TH1F *hBProtonTOFSignal=new TH1F(hisname.Data(),hisname.Data(),100,12000.,50000.0);
  fOutputMC->Add(hBProtonTOFSignal);

  hisname.Form("hpTPCSignal");
  TH1F *hProtonTPCSignal=new TH1F(hisname.Data(),hisname.Data(),150,0.,150.0);
  fOutputMC->Add(hProtonTPCSignal);
  hisname.Form("hbpTPCSignal");
  TH1F *hBProtonTPCSignal=new TH1F(hisname.Data(),hisname.Data(),150,0.,150.0);
  fOutputMC->Add(hBProtonTPCSignal);
  
  hisname.Form("hpptProng");
  TH1F *hProtonPtProng=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.0);
  fOutputMC->Add(hProtonPtProng);
  hisname.Form("hbpptProng");
  TH1F *hBProtonPtProng=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.0);
  fOutputMC->Add(hBProtonPtProng);

  hisname.Form("hpRealTot");
  TH1F *hProtonRealTot=new TH1F(hisname.Data(),hisname.Data(),100,0.,10.0);
  fOutputMC->Add(hProtonRealTot);
  hisname.Form("hbpRealTot");
  TH1F *hBProtonRealTot=new TH1F(hisname.Data(),hisname.Data(),100,0.,10.0);
  fOutputMC->Add(hBProtonRealTot);
  hisname.Form("hpIDTot");
  TH1F *hProtonIDTot=new TH1F(hisname.Data(),hisname.Data(),100,0.,10.0);
  fOutputMC->Add(hProtonIDTot);
  hisname.Form("hpIDGood");
  TH1F *hProtonIDGood=new TH1F(hisname.Data(),hisname.Data(),100,0.,10.0);
  fOutputMC->Add(hProtonIDGood);
  hisname.Form("hbpIDGood");
  TH1F *hBProtonIDGood=new TH1F(hisname.Data(),hisname.Data(),100,0.,10.0);
  fOutputMC->Add(hBProtonIDGood);
  hisname.Form("hbpIDTot");
  TH1F *hBProtonIDTot=new TH1F(hisname.Data(),hisname.Data(),100,0.,10.0);
  fOutputMC->Add(hBProtonIDTot);
  hisname.Form("hpnonIDTot");
  TH1F *hProtonnonIDTot=new TH1F(hisname.Data(),hisname.Data(),100,0.,10.0);
  fOutputMC->Add(hProtonnonIDTot);
  hisname.Form("hbpnonIDTot");
  TH1F *hBProtonnonIDTot=new TH1F(hisname.Data(),hisname.Data(),100,0.,10.0);
  fOutputMC->Add(hBProtonnonIDTot);

  hisname.Form("hpd0Prong");
  TH1F *hProtond0Prong=new TH1F(hisname.Data(),hisname.Data(),100,-0.1,0.1);
  fOutputMC->Add(hProtond0Prong);
  hisname.Form("hbpd0Prong");
  TH1F *hBProtond0Prong=new TH1F(hisname.Data(),hisname.Data(),100,-0.1,0.1);
  fOutputMC->Add(hBProtond0Prong);
  hisname.Form("hbpSignalVspTOF");
  TH2F *hBpSignalVspTOF=new TH2F(hisname.Data(),hisname.Data(),100,0.,5.0,100,12000.,50000.0);
  fOutputMC->Add(hBpSignalVspTOF);
  hisname.Form("hbpSignalVspTPC");
  TH2F *hBpSignalVspTPC=new TH2F(hisname.Data(),hisname.Data(),100,0.,5.0,150,0.,150.0);
  fOutputMC->Add(hBpSignalVspTPC);
  hisname.Form("hpSignalVspTOF");
  TH2F *hpSignalVspTOF=new TH2F(hisname.Data(),hisname.Data(),100,0.,5.0,100,12000.,50000.0);
  fOutputMC->Add(hpSignalVspTOF);
  hisname.Form("hpSignalVspTPC");
  TH2F *hpSignalVspTPC=new TH2F(hisname.Data(),hisname.Data(),100,0.,5.0,150,0.,150.0);
  fOutputMC->Add(hpSignalVspTPC);

  hisname.Form("hpSigmaVspTOF");
  TH2F *hpSigmaVspTOF=new TH2F(hisname.Data(),hisname.Data(),100,0.,5.0,100,-10.0,10.0);
  fOutputMC->Add(hpSigmaVspTOF);
  hisname.Form("hpSigmaVspTPC");
  TH2F *hpSigmaVspTPC=new TH2F(hisname.Data(),hisname.Data(),100,0.,5.0,100,-10.0,10.0);
  fOutputMC->Add(hpSigmaVspTPC);
  hisname.Form("hbpSigmaVspTOF");
  TH2F *hBpSigmaVspTOF=new TH2F(hisname.Data(),hisname.Data(),100,0.,5.0,100,-10.0,10.0);
  fOutputMC->Add(hBpSigmaVspTOF);
  hisname.Form("hbpSigmaVspTPC");
  TH2F *hBpSigmaVspTPC=new TH2F(hisname.Data(),hisname.Data(),100,0.,5.0,100,-10.0,10.0);
  fOutputMC->Add(hBpSigmaVspTPC);

  //kaon specific

  hisname.Form("hKTOFSignal");
  TH1F *hKaonTOFSignal=new TH1F(hisname.Data(),hisname.Data(),100,12000.,50000.0);
  fOutputMC->Add(hKaonTOFSignal);
  hisname.Form("hbKTOFSignal");
  TH1F *hBKaonTOFSignal=new TH1F(hisname.Data(),hisname.Data(),100,12000.,50000.0);
  fOutputMC->Add(hBKaonTOFSignal);

  hisname.Form("hKTPCSignal");
  TH1F *hKaonTPCSignal=new TH1F(hisname.Data(),hisname.Data(),150,0.,150.0);
  fOutputMC->Add(hKaonTPCSignal);
  hisname.Form("hbKTPCSignal");
  TH1F *hBKaonTPCSignal=new TH1F(hisname.Data(),hisname.Data(),150,0.,150.0);
  fOutputMC->Add(hBKaonTPCSignal);

  hisname.Form("hKptProng");
  TH1F *hKaonPtProng=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.0);
  fOutputMC->Add(hKaonPtProng);
  hisname.Form("hbKptProng");
  TH1F *hBKaonPtProng=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.0);
  fOutputMC->Add(hBKaonPtProng);

  hisname.Form("hKRealTot");
  TH1F *hKaonRealTot=new TH1F(hisname.Data(),hisname.Data(),100,0.,10.0);
  fOutputMC->Add(hKaonRealTot);
  hisname.Form("hbKRealTot");
  TH1F *hBKaonRealTot=new TH1F(hisname.Data(),hisname.Data(),100,0.,10.0);
  fOutputMC->Add(hBKaonRealTot);
  hisname.Form("hKIDGood");
  TH1F *hKaonIDGood=new TH1F(hisname.Data(),hisname.Data(),100,0.,10.0);
  fOutputMC->Add(hKaonIDGood);
  hisname.Form("hKIDTot");
  TH1F *hKaonIDTot=new TH1F(hisname.Data(),hisname.Data(),100,0.,10.0);
  fOutputMC->Add(hKaonIDTot);
  hisname.Form("hbKIDGood");
  TH1F *hBKaonIDGood=new TH1F(hisname.Data(),hisname.Data(),100,0.,10.0);
  fOutputMC->Add(hBKaonIDGood);
  hisname.Form("hbKIDTot");
  TH1F *hBKaonIDTot=new TH1F(hisname.Data(),hisname.Data(),100,0.,10.0);
  fOutputMC->Add(hBKaonIDTot);
  hisname.Form("hKnonIDTot");
  TH1F *hKaonnonIDTot=new TH1F(hisname.Data(),hisname.Data(),100,0.,10.0);
  fOutputMC->Add(hKaonnonIDTot);
  hisname.Form("hbKnonIDTot");
  TH1F *hBKaonnonIDTot=new TH1F(hisname.Data(),hisname.Data(),100,0.,10.0);
  fOutputMC->Add(hBKaonnonIDTot);


  hisname.Form("hKd0Prong");
  TH1F *hKaond0Prong=new TH1F(hisname.Data(),hisname.Data(),100,-0.1,0.1);
  fOutputMC->Add(hKaond0Prong);
  hisname.Form("hbKd0Prong");
  TH1F *hBKaond0Prong=new TH1F(hisname.Data(),hisname.Data(),100,-0.1,0.1);
  fOutputMC->Add(hBKaond0Prong);
  hisname.Form("hbKSignalVspTOF");
  TH2F *hbKSignalVspTOF=new TH2F(hisname.Data(),hisname.Data(),100,0.,5.0,100,12000.,50000.0);
  fOutputMC->Add(hbKSignalVspTOF);
  hisname.Form("hbKSignalVspTPC");
  TH2F *hbKSignalVspTPC=new TH2F(hisname.Data(),hisname.Data(),100,0.,5.0,150,0.,150.0);
  fOutputMC->Add(hbKSignalVspTPC);
  hisname.Form("hKSignalVspTOF");
  TH2F *hKSignalVspTOF=new TH2F(hisname.Data(),hisname.Data(),100,0.,5.0,100,12000.,50000.0);
  fOutputMC->Add(hKSignalVspTOF);
  hisname.Form("hKSignalVspTPC");
  TH2F *hKSignalVspTPC=new TH2F(hisname.Data(),hisname.Data(),100,0.,5.0,150,0.,150.0);
  fOutputMC->Add(hKSignalVspTPC);

  hisname.Form("hKSigmaVspTOF");
  TH2F *hKSigmaVspTOF=new TH2F(hisname.Data(),hisname.Data(),100,0.,5.0,100,-10.0,10.0);
  fOutputMC->Add(hKSigmaVspTOF);

  hisname.Form("hKSigmaVspTPC");
  TH2F *hKSigmaVspTPC=new TH2F(hisname.Data(),hisname.Data(),100,0.,5.0,100,-10.0,10.0);
  fOutputMC->Add(hKSigmaVspTPC);

  hisname.Form("hbKSigmaVspTOF");
  TH2F *hBKSigmaVspTOF=new TH2F(hisname.Data(),hisname.Data(),100,0.,5.0,100,-10.0,10.0);
  fOutputMC->Add(hBKSigmaVspTOF);
  hisname.Form("hbKSigmaVspTPC");
  TH2F *hBKSigmaVspTPC=new TH2F(hisname.Data(),hisname.Data(),100,0.,5.0,100,-10.0,10.0);
  fOutputMC->Add(hBKSigmaVspTPC);

  // pion specific
  hisname.Form("hpiTOFSignal");
  TH1F *hPionTOFSignal=new TH1F(hisname.Data(),hisname.Data(),100,12000.,50000.0);
  fOutputMC->Add(hPionTOFSignal);
  hisname.Form("hbpiTOFSignal");
  TH1F *hBPionTOFSignal=new TH1F(hisname.Data(),hisname.Data(),100,12000.,50000.0);
  fOutputMC->Add(hBPionTOFSignal);

  hisname.Form("hpiTPCSignal");
  TH1F *hPionTPCSignal=new TH1F(hisname.Data(),hisname.Data(),100,30.,100.0);
  fOutputMC->Add(hPionTPCSignal);
  hisname.Form("hbpiTPCSignal");
  TH1F *hBPionTPCSignal=new TH1F(hisname.Data(),hisname.Data(),100,30.,100.0);
  fOutputMC->Add(hBPionTPCSignal);
  
  hisname.Form("hpiptProng");
  TH1F *hPionPtProng=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.0);
  fOutputMC->Add(hPionPtProng);
  hisname.Form("hbpiptProng");
  TH1F *hBPionPtProng=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.0);
  fOutputMC->Add(hBPionPtProng);

  hisname.Form("hpiRealTot");
  TH1F *hPionRealTot=new TH1F(hisname.Data(),hisname.Data(),100,0.,10.0);
  fOutputMC->Add(hPionRealTot);
  hisname.Form("hbpiRealTot");
  TH1F *hBPionRealTot=new TH1F(hisname.Data(),hisname.Data(),100,0.,10.0);
  fOutputMC->Add(hBPionRealTot);
  hisname.Form("hpiIDGood");
  TH1F *hPionIDGood=new TH1F(hisname.Data(),hisname.Data(),100,0.,10.0);
  fOutputMC->Add(hPionIDGood);
  hisname.Form("hpiIDTot");
  TH1F *hPionIDTot=new TH1F(hisname.Data(),hisname.Data(),100,0.,10.0);
  fOutputMC->Add(hPionIDTot);
  hisname.Form("hbpiIDTot");
  TH1F *hBPionIDTot=new TH1F(hisname.Data(),hisname.Data(),100,0.,10.0);
  fOutputMC->Add(hBPionIDTot);
  hisname.Form("hbpiIDGood");
  TH1F *hBPionIDGood=new TH1F(hisname.Data(),hisname.Data(),100,0.,10.0);
  fOutputMC->Add(hBPionIDGood);
  hisname.Form("hpinonIDTot");
  TH1F *hPionnonIDTot=new TH1F(hisname.Data(),hisname.Data(),100,0.,10.0);
  fOutputMC->Add(hPionnonIDTot);
  hisname.Form("hbpinonIDTot");
  TH1F *hBPionnonIDTot=new TH1F(hisname.Data(),hisname.Data(),100,0.,10.0);
  fOutputMC->Add(hBPionnonIDTot);


  hisname.Form("hpid0Prong");
  TH1F *hPiond0Prong=new TH1F(hisname.Data(),hisname.Data(),100,-0.1,-0.1);
  fOutputMC->Add(hPiond0Prong);
  hisname.Form("hbpid0Prong");
  TH1F *hBPiond0Prong=new TH1F(hisname.Data(),hisname.Data(),100,-0.1,0.1);
  fOutputMC->Add(hBPiond0Prong);
  hisname.Form("hbpiSignalVspTOF");
  TH2F *hbpiSignalVspTOF=new TH2F(hisname.Data(),hisname.Data(),100,0.,5.0,100,12000.,50000.0);
  fOutputMC->Add(hbpiSignalVspTOF);
  hisname.Form("hbpiSignalVspTPC");
  TH2F *hbpiSignalVspTPC=new TH2F(hisname.Data(),hisname.Data(),100,0.,5.0,150,0.,150.0);
  fOutputMC->Add(hbpiSignalVspTPC);

  hisname.Form("hpiSignalVspTOF");
  TH2F *hpiSignalVspTOF=new TH2F(hisname.Data(),hisname.Data(),100,0.,5.0,100,12000.,50000.0);
  fOutputMC->Add(hpiSignalVspTOF);
  hisname.Form("hpiSignalVspTPC");
  TH2F *hpiSignalVspTPC=new TH2F(hisname.Data(),hisname.Data(),100,0.,5.0,150,0.,150.0);
  fOutputMC->Add(hpiSignalVspTPC);

  hisname.Form("hbpiSigmaVspTOF");
  TH2F *hBpiSigmaVspTOF=new TH2F(hisname.Data(),hisname.Data(),100,0.,5.0,100,-10.0,10.0);
  fOutputMC->Add(hBpiSigmaVspTOF);

  hisname.Form("hbpiSigmaVspTPC");
  TH2F *hBpiSigmaVspTPC=new TH2F(hisname.Data(),hisname.Data(),100,0.,5.0,100,-10.0,10.0);
  fOutputMC->Add(hBpiSigmaVspTPC);
  hisname.Form("hpiSigmaVspTOF");
  TH2F *hpiSigmaVspTOF=new TH2F(hisname.Data(),hisname.Data(),100,0.,5.0,100,-10.0,10.0);
  fOutputMC->Add(hpiSigmaVspTOF);

  hisname.Form("hpiSigmaVspTPC");
  TH2F *hpiSigmaVspTPC=new TH2F(hisname.Data(),hisname.Data(),100,0.,5.0,100,-10.0,10.0);
  fOutputMC->Add(hpiSigmaVspTPC);

  // other generic 
  hisname.Form("hLcpt");
  TH1F *hLcPt=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.0);
  fOutputMC->Add(hLcPt);
  hisname.Form("hbLcpt");
  TH1F *hBLcPt=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.0);
  fOutputMC->Add(hBLcPt);
  hisname.Form("hDist12toPrim");
  TH1F *hDist12Prim=new TH1F(hisname.Data(),hisname.Data(),100,0.,1.0);
  fOutputMC->Add(hDist12Prim);
  hisname.Form("hbDist12toPrim");
  TH1F *hBDist12Prim=new TH1F(hisname.Data(),hisname.Data(),100,0.,1.0);
  fOutputMC->Add(hBDist12Prim);

  hisname.Form("hSigmaVert");
  TH1F *hSigmaVert=new TH1F(hisname.Data(),hisname.Data(),60,0.,0.06);
  fOutputMC->Add(hSigmaVert);
  hisname.Form("hbSigmaVert");
  TH1F *hBSigmaVert=new TH1F(hisname.Data(),hisname.Data(),60,0.,0.06);
  fOutputMC->Add(hBSigmaVert);

  hisname.Form("hDCAs");
  TH1F *hdcas=new TH1F(hisname.Data(),hisname.Data(),200,0.,0.1);
  fOutputMC->Add(hdcas);
  hisname.Form("hbDCAs");
  TH1F *hBdcas=new TH1F(hisname.Data(),hisname.Data(),200,0.,0.1);
  fOutputMC->Add(hBdcas);

  hisname.Form("hCosPointingAngle");
  TH1F *hCosPointingAngle=new TH1F(hisname.Data(),hisname.Data(),40,0.,1.);
  fOutputMC->Add(hCosPointingAngle);
  hisname.Form("hbCosPointingAngle");
  TH1F *hBCosPointingAngle=new TH1F(hisname.Data(),hisname.Data(),40,0.,1.);
  fOutputMC->Add(hBCosPointingAngle);

  hisname.Form("hDecayLength");
  TH1F *hDecayLength=new TH1F(hisname.Data(),hisname.Data(),100,0.,0.1);
  fOutputMC->Add(hDecayLength);
  hisname.Form("hbDecayLength");
  TH1F *hBDecayLength=new TH1F(hisname.Data(),hisname.Data(),100,0.,0.1);
  fOutputMC->Add(hBDecayLength);

  hisname.Form("hSum2");
  TH1F *hSum2=new TH1F(hisname.Data(),hisname.Data(),100,0.,0.1);
  fOutputMC->Add(hSum2);
  hisname.Form("hbSum2");
  TH1F *hBSum2=new TH1F(hisname.Data(),hisname.Data(),100,0.,0.1);
  fOutputMC->Add(hBSum2);

  hisname.Form("hptmax");
  TH1F *hPtMax=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.);
  fOutputMC->Add(hPtMax);
  hisname.Form("hbptmax");
  TH1F *hBPtMax=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.);
  fOutputMC->Add(hBPtMax);



  if(fRDCutsProduction->GetIsUsePID()){
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler=(AliInputEventHandler*)mgr->GetInputEventHandler();
   AliPIDResponse *pidResp=inputHandler->GetPIDResponse();
   fRDCutsProduction->GetPidHF()->SetPidResponse(pidResp);
   fRDCutsProduction->GetPidpion()->SetPidResponse(pidResp);
   fRDCutsProduction->GetPidprot()->SetPidResponse(pidResp);
   fUtilPid=new AliAODpidUtil(pidResp);
   fUtilPid=new AliAODpidUtil();
   }
  if(fRDCutsAnalysis->GetIsUsePID()){
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   AliInputEventHandler *inputHandler=(AliInputEventHandler*)mgr->GetInputEventHandler();
   AliPIDResponse *pidResp=inputHandler->GetPIDResponse();
  fRDCutsAnalysis->GetPidHF()->SetPidResponse(pidResp);
   fRDCutsAnalysis->GetPidpion()->SetPidResponse(pidResp);
   fRDCutsAnalysis->GetPidprot()->SetPidResponse(pidResp);
  }


  PostData(1,fOutput);
  if(fFillVarHists) PostData(3,fOutputMC);
  PostData(4,fNentries);
  if(fFillNtuple){
    //OpenFile(3); // 2 is the slot number of the ntuple
   
    fNtupleLambdac = new TNtuple("fNtupleLambdac","D +","pdg:Px:Py:Pz:PtTrue:VxTrue:VyTrue:VzTrue:Ptpi:PtK:Ptpi2:PtRec:PointingAngle:DecLeng:VxRec:VyRec:VzRec:InvMass:sigvert:d0Pi:d0K:d0Pi2:dca:d0square");  
  PostData(5,fNtupleLambdac);
    
  }
  
  return;
}

//________________________________________________________________________
void AliAnalysisTaskSELambdac::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event:
  // heavy flavor candidates association to MC truth

   AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());
   //tmp
 fHistNEvents->Fill(0); // count event
  // Post the data already here

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
  } else if(aod) {
    array3Prong=(TClonesArray*)aod->GetList()->FindObject("Charm3Prong");
    arrayLikeSign=(TClonesArray*)aod->GetList()->FindObject("LikeSign3Prong");
  }

  if(!array3Prong || !aod) {
    printf("AliAnalysisTaskSELambdac::UserExec: Charm3Prong branch not found!\n");
    return;
  }
  if(!arrayLikeSign) {
    printf("AliAnalysisTaskSELambdac::UserExec: LikeSign3Prong branch not found!\n");
  //  return;
  }

 
  TClonesArray *arrayMC=0;
  AliAODMCHeader *mcHeader=0;

  fNentries->Fill(0);
  TString trigclass=aod->GetFiredTriggerClasses();
    if(trigclass.Contains("C0SMH-B-NOPF-ALLNOTRD") || trigclass.Contains("C0SMH-B-NOPF-ALL")) fNentries->Fill(14);
    if(!fRDCutsProduction->IsEventSelected(aod)) {
     if(fRDCutsProduction->GetWhyRejection()==1) // rejected for pileup
      fNentries->Fill(13);
      return;
    }

  // AOD primary vertex
  AliAODVertex *vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();
  if(!vtx1) return;

  
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
  
  
  for (Int_t i3Prong = 0; i3Prong < n3Prong; i3Prong++) {
    AliAODRecoDecayHF3Prong *d = (AliAODRecoDecayHF3Prong*)array3Prong->UncheckedAt(i3Prong);
    
    
    Bool_t unsetvtx=kFALSE;
    if(!d->GetOwnPrimaryVtx()){
      d->SetOwnPrimaryVtx(vtx1);
      unsetvtx=kTRUE;
    }

//    if(d->GetSelectionMap()) {if(!d->HasSelectionBit(AliRDHFCuts::kLcCuts)) continue;}
//   if(d->GetSelectionMap()) if(!d->HasSelectionBit(AliRDHFCuts::kLcPID)) continue;

    Int_t isSelectedTracks = fRDCutsProduction->IsSelected(d,AliRDHFCuts::kTracks,aod);
    if(!isSelectedTracks) continue;

    if (fRDCutsProduction->IsInFiducialAcceptance(d->Pt(),d->Y(4122))) fNentries->Fill(6);

    Int_t ptbin=fRDCutsProduction->PtBin(d->Pt());
    if(ptbin==-1) {fNentries->Fill(4); continue;} //out of bounds

    FillMassHists(aod,d,arrayMC,fRDCutsProduction);
    if(fFillVarHists) FillVarHists(d,arrayMC,fRDCutsProduction,fOutputMC,aod);
    
      /*
      //start OS analysis
      if(labDp<0)fHistOSbkg->Fill(d->InvMassDplus());
      fHistOS->Fill(d->InvMassDplus());
      */

    if(unsetvtx) d->UnsetOwnPrimaryVtx();
  }

  PostData(1,fOutput); 
  if(fFillVarHists) PostData(3,fOutputMC);
  PostData(4,fNentries);
      
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
 //fHistNEvents = dynamic_cast<TH1F*>(fOutput->FindObject("fHistNEvents"));

 TString hisname;
  if(fFillNtuple){
    fNtupleLambdac = dynamic_cast<TNtuple*>(GetOutputData(3));
  }

 
 return;
}

//________________________________________________________________________
Int_t AliAnalysisTaskSELambdac::MatchToMCLambdac(AliAODRecoDecayHF3Prong *d,TClonesArray *arrayMC) const{
  // check if the candidate is a Lambdac decaying in pKpi or in the resonant channels
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
Bool_t AliAnalysisTaskSELambdac::GetLambdacDaugh(AliAODMCParticle *part,TClonesArray *arrayMC) const{
 // check if the particle is a lambdac and if its decay mode is the correct one 
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
Bool_t AliAnalysisTaskSELambdac::IspKpiMC(AliAODRecoDecayHF3Prong *d,TClonesArray *arrayMC) const{
  // Apply MC PID
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
Bool_t AliAnalysisTaskSELambdac::IspiKpMC(AliAODRecoDecayHF3Prong *d,TClonesArray *arrayMC) const{

  // Apply MC PID
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
Bool_t AliAnalysisTaskSELambdac::VertexingKF(AliAODRecoDecayHF3Prong *d,Int_t *pdgs,Double_t field) const{
 // apply vertexing KF 
   Int_t iprongs[3]={0,1,2};
   Double_t mass[2]={0.,0.};
 //topological constr
   AliKFParticle *lambdac=d->ApplyVertexingKF(iprongs,3,pdgs,kTRUE,field,mass);
   if(!lambdac) return kFALSE;
//  Double_t probTot=TMath::Prob(Lambdac->GetChi2(),Lambdac->GetNDF());
//  if(probTot<fCutsKF[0]) return kFALSE;
  if(lambdac->GetChi2()>fCutsKF[0]) return kFALSE;
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
   AliKFParticle *kappaStar=d->ApplyVertexingKF(ipRes,2,pdgres,kFALSE,field,mass);

  Double_t probKstar=TMath::Prob(kappaStar->GetChi2(),kappaStar->GetNDF());
  if(probKstar>fCutsKF[1]) {
    AliAODTrack *esdProng1=(AliAODTrack*)d->GetDaughter(ipRes[0]);
    AliAODTrack *esdProng2=(AliAODTrack*)d->GetDaughter(ipRes[1]);
    AliKFParticle prong1(*esdProng1,pdgres[0]);
    AliKFParticle prong2(*esdProng2,pdgres[1]);
    if(kappaStar->GetPt()<fCutsKF[2] && prong1.GetAngle(prong2)>fCutsKF[3]) return kFALSE;
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
   AliKFParticle *lambda1520=d->ApplyVertexingKF(ipRes,2,pdgres,kFALSE,field,mass);
  Double_t probLa=TMath::Prob(lambda1520->GetChi2(),lambda1520->GetNDF());
  if(probLa>fCutsKF[4]) {
    AliAODTrack *esdProng1=(AliAODTrack*)d->GetDaughter(ipRes[0]);
    AliAODTrack *esdProng2=(AliAODTrack*)d->GetDaughter(ipRes[1]);
    AliKFParticle prong1(*esdProng1,pdgres[0]);
    AliKFParticle prong2(*esdProng2,pdgres[1]);
    if(lambda1520->GetPt()<fCutsKF[5] && prong1.GetAngle(prong2)>fCutsKF[6]) return kFALSE;
  } 
 //mass constr for Delta
   mass[0]=1.2;mass[1]=0.15;
   ipRes[0]=0;ipRes[1]=2;
   pdgres[0]=pdgs[0];pdgres[1]=pdgs[2];
   AliKFParticle *delta=d->ApplyVertexingKF(ipRes,2,pdgres,kFALSE,field,mass);
  Double_t probDelta=TMath::Prob(delta->GetChi2(),delta->GetNDF());
  if(probDelta>fCutsKF[7]) {
    AliAODTrack *esdProng1=(AliAODTrack*)d->GetDaughter(ipRes[0]);
    AliAODTrack *esdProng2=(AliAODTrack*)d->GetDaughter(ipRes[1]);
    AliKFParticle prong1(*esdProng1,pdgres[0]);
    AliKFParticle prong2(*esdProng2,pdgres[1]);
    if(delta->GetPt()<fCutsKF[8] && prong1.GetAngle(prong2)>fCutsKF[9]) return kFALSE;
  } 
 return kTRUE;
}
//-------------------------------------
Int_t AliAnalysisTaskSELambdac::IspiKpResonant(AliAODRecoDecayHF3Prong *d,Double_t field) const{
  
 // apply PID using the resonant channels 
 //if lambda* -> pk
        Double_t mass[2]={1.520,0.005};
        Int_t ipRes[2]={1,2};
        Int_t pdgres[2]={321,2212};
        AliKFParticle *lambda1520=d->ApplyVertexingKF(ipRes,2,pdgres,kFALSE,field,mass);
	Double_t probLa=TMath::Prob(lambda1520->GetChi2(),lambda1520->GetNDF());
	if(probLa>0.1) return 1;
 //K* -> kpi
        mass[0]=0.8961;mass[1]=0.03;
        ipRes[0]=0;ipRes[1]=1;
        pdgres[0]=211;pdgres[1]=321;
        AliKFParticle *kstar=d->ApplyVertexingKF(ipRes,2,pdgres,kFALSE,field,mass);
	Double_t probKa=TMath::Prob(kstar->GetChi2(),kstar->GetNDF());
       if(probKa>0.1) return 2;

 return 0;

}
//-------------------------------------
Int_t AliAnalysisTaskSELambdac::IspKpiResonant(AliAODRecoDecayHF3Prong *d,Double_t field) const{
   
 // apply PID using the resonant channels 
 //if lambda* -> pk
        Double_t mass[2]={1.520,0.005};
        Int_t ipRes[2]={0,1};
        Int_t pdgres[2]={2212,321};
        AliKFParticle *lambda1520=d->ApplyVertexingKF(ipRes,2,pdgres,kFALSE,field,mass);
	Double_t probLa=TMath::Prob(lambda1520->GetChi2(),lambda1520->GetNDF());
	if(probLa>0.1) return 1;
 //K* -> kpi
        mass[0]=0.8961;mass[1]=0.03;
        ipRes[0]=1;ipRes[1]=2;
        pdgres[1]=211;pdgres[0]=321;
        AliKFParticle *kstar=d->ApplyVertexingKF(ipRes,2,pdgres,kFALSE,field,mass);
	Double_t probKa=TMath::Prob(kstar->GetChi2(),kstar->GetNDF());
	if(probKa>0.1) return 2;

 return 0;

}
//---------------------------
void AliAnalysisTaskSELambdac::FillMassHists(AliAODEvent *aod,AliAODRecoDecayHF3Prong *part, TClonesArray *arrayMC, AliRDHFCutsLctopKpi *cuts){

 //if MC PID or no PID, unset PID
 if(!fRealPid) cuts->SetUsePID(kFALSE);
 Int_t selection=cuts->IsSelected(part,AliRDHFCuts::kCandidate,aod);
 if(selection>0){
 Int_t iPtBin = -1;
 Double_t ptCand = part->Pt();
 Int_t index=0;
 
 for(Int_t ibin=0;ibin<fNPtBins&&iPtBin<0&&ptCand>fArrayBinLimits[0]&&ptCand<fArrayBinLimits[fNPtBins];ibin++){
  if(ptCand<fArrayBinLimits[ibin+1])iPtBin=ibin;
 }
 
 
 Float_t pdgCode=-2;
 Float_t pdgCode1=-1;
 Float_t pdgCode2=-1;
 Int_t labDp=-1;
 Float_t deltaPx=0.;
 Float_t deltaPy=0.;
 Float_t deltaPz=0.;
 Float_t truePt=0.;
 Float_t xDecay=0.;
 Float_t yDecay=0.;
 Float_t zDecay=0.;

 if(fReadMC){
  labDp = MatchToMCLambdac(part,arrayMC);
  if(labDp>0){
   AliAODMCParticle *partDp = (AliAODMCParticle*)arrayMC->At(labDp);
          AliAODMCParticle *dg0 = (AliAODMCParticle*)arrayMC->At(partDp->GetDaughter(0));
          AliAODMCParticle *dg1 = (AliAODMCParticle*)arrayMC->At(partDp->GetDaughter(1));
          deltaPx=partDp->Px()-part->Px();
          deltaPy=partDp->Py()-part->Py();
          deltaPz=partDp->Pz()-part->Pz();
          truePt=partDp->Pt();
          xDecay=dg0->Xv();
          yDecay=dg0->Yv();
          zDecay=dg0->Zv();
          pdgCode=TMath::Abs(partDp->GetPdgCode());
          pdgCode1=TMath::Abs(dg0->GetPdgCode());
          pdgCode2=TMath::Abs(dg1->GetPdgCode());

  }else{
   pdgCode=-1;
  }
 }

  Double_t invMasspKpi=-1.;
  Double_t invMasspiKp=-1.;
  Double_t invMasspKpiLpi=-1.;
  Double_t invMasspiKpLpi=-1.;
  Double_t invMasspKpiKp=-1.;
  Double_t invMasspiKpKp=-1.;
  Int_t pdgs[3]={0,0,0};
  Double_t field=aod->GetMagneticField();
//apply MC PID
  if(fReadMC && fMCPid){

  if(IspKpiMC(part,arrayMC)) {
   invMasspKpi=part->InvMassLcpKpi();
   if(fUseKF){
    pdgs[0]=2212;pdgs[1]=321;pdgs[2]=211;
    if(!VertexingKF(part,pdgs,field)) invMasspKpi=-1.;
   }
  }
  if(IspiKpMC(part,arrayMC)) {
   invMasspiKp=part->InvMassLcpiKp();
   if(fUseKF){
    pdgs[0]=211;pdgs[1]=321;pdgs[2]=2212;
    if(!VertexingKF(part,pdgs,field)) invMasspiKp=-1.;
   }
  }
 }
  
 // apply realistic PID
  if(fRealPid){
   if(selection==1 || selection==3) {
    invMasspKpi=part->InvMassLcpKpi();
    pdgs[0]=2212;pdgs[1]=321;pdgs[2]=211;
    if(fUseKF){
     pdgs[0]=2212;pdgs[1]=321;pdgs[2]=211;
     if(!VertexingKF(part,pdgs,field)) invMasspKpi=-1.;
    }
   }
   if(selection>=2) {
    invMasspiKp=part->InvMassLcpiKp();
    pdgs[0]=211;pdgs[1]=321;pdgs[2]=2212;
    if(fUseKF){
     pdgs[0]=211;pdgs[1]=321;pdgs[2]=2212;
     if(!VertexingKF(part,pdgs,field)) invMasspiKp=-1.;
    }
   }
  } 
  //apply PID using resonances
  if(fResPid && fRealPid){
   if((selection==3 || selection==1)) {
    if(IspKpiResonant(part,field)==1){
     invMasspKpiLpi=part->InvMassLcpKpi();
    }
    if(IspKpiResonant(part,field)==2){
     invMasspKpiKp=part->InvMassLcpKpi();
    }
   }
   if(selection>=2) {
    if(IspiKpResonant(part,field)==1){
     invMasspiKpLpi=part->InvMassLcpiKp();
    }
    if(IspiKpResonant(part,field)==2){
     invMasspiKpKp=part->InvMassLcpiKp();
    }
   }
  }
  
  // no PID
  if(!fResPid && !fRealPid && !fMCPid){
   if(selection==2 || selection==3) {
    invMasspiKp=part->InvMassLcpiKp();
    if(fUseKF){
     pdgs[0]=211;pdgs[1]=321;pdgs[2]=2212;
     if(!VertexingKF(part,pdgs,field)) invMasspiKp=-1.;
    }
   }
   if(selection==1 || selection==3){
    invMasspKpi=part->InvMassLcpKpi();
    if(fUseKF){
     pdgs[0]=2212;pdgs[1]=321;pdgs[2]=211;
     if(!VertexingKF(part,pdgs,field)) invMasspKpi=-1.;
     }
    }
   }

   if(invMasspiKp<0. && invMasspKpi<0.) return;
   
   Int_t passTightCuts=0;
   if(fAnalysis) passTightCuts=fRDCutsAnalysis->IsSelected(part,AliRDHFCuts::kCandidate,aod);

   

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
        if(pdgCode1==2212) {tmp[8]=part->PtProng(0);}else{tmp[8]=0.;}
        if(pdgCode1==211) {tmp[10]=part->PtProng(0);}else{tmp[10]=0.;}
        tmp[9]=part->PtProng(1);
        if(pdgCode2==211) {tmp[10]=part->PtProng(2);}else{tmp[10]=0.;}
        tmp[11]=part->Pt();
        tmp[12]=part->CosPointingAngle();
        tmp[13]=part->DecayLength();
        tmp[14]=part->Xv();
        tmp[15]=part->Yv();
        tmp[16]=part->Zv();
        if(invMasspiKp>0.) tmp[17]=invMasspiKp;
        if(invMasspKpi>0.) tmp[17]=invMasspKpi;
        tmp[18]=part->GetSigmaVert();
        tmp[19]=part->Getd0Prong(0);
        tmp[20]=part->Getd0Prong(1);
        tmp[21]=part->Getd0Prong(2);
        tmp[22]=part->GetDCA();
        tmp[23]=part->Prodd0d0();
        fNtupleLambdac->Fill(tmp);
        PostData(5,fNtupleLambdac);
      }
    
    if(part->Pt()>3.){
     if(invMasspiKp>0. && invMasspKpi>0.){
      if(invMasspiKp>0.) fhMassPtGreater3->Fill(invMasspiKp,0.5);
      if(invMasspKpi>0.) fhMassPtGreater3->Fill(invMasspKpi,0.5);
     }else{
      if(invMasspiKp>0.) fhMassPtGreater3->Fill(invMasspiKp);
      if(invMasspKpi>0.) fhMassPtGreater3->Fill(invMasspKpi);
     }
     if(invMasspiKpLpi>0. && invMasspKpiLpi>0.){
      if(invMasspiKpLpi>0.) fhMassPtGreater3Lpi->Fill(invMasspiKpLpi,0.5);
      if(invMasspKpiLpi>0.) fhMassPtGreater3Lpi->Fill(invMasspKpiLpi,0.5);
     }else{
      if(invMasspiKpLpi>0.) fhMassPtGreater3Lpi->Fill(invMasspiKpLpi);
      if(invMasspKpiLpi>0.) fhMassPtGreater3Lpi->Fill(invMasspKpiLpi);
     }
     if(invMasspiKpKp>0. && invMasspKpiKp>0.){
      if(invMasspiKpKp>0.) fhMassPtGreater3Kp->Fill(invMasspiKpKp,0.5);
      if(invMasspKpiKp>0.) fhMassPtGreater3Kp->Fill(invMasspKpiKp,0.5);
     }else{
      if(invMasspiKpKp>0.) fhMassPtGreater3Kp->Fill(invMasspiKpKp);
      if(invMasspKpiKp>0.) fhMassPtGreater3Kp->Fill(invMasspKpiKp);
     }
     if(passTightCuts>0){
      if(invMasspiKp>0. && invMasspKpi>0.){
       if(invMasspiKp>0.) fhMassPtGreater3TC->Fill(invMasspiKp,0.5);
       if(invMasspKpi>0.) fhMassPtGreater3TC->Fill(invMasspKpi,0.5);
      }else{
       if(invMasspiKp>0.) fhMassPtGreater3TC->Fill(invMasspiKp);
       if(invMasspKpi>0.) fhMassPtGreater3TC->Fill(invMasspKpi);
      }
     if(invMasspiKpLpi>0. && invMasspKpiLpi>0.){
      if(invMasspiKpLpi>0.) fhMassPtGreater3LpiTC->Fill(invMasspiKpLpi,0.5);
      if(invMasspKpiLpi>0.) fhMassPtGreater3LpiTC->Fill(invMasspKpiLpi,0.5);
     }else{
      if(invMasspiKpLpi>0.) fhMassPtGreater3LpiTC->Fill(invMasspiKpLpi);
      if(invMasspKpiLpi>0.) fhMassPtGreater3LpiTC->Fill(invMasspKpiLpi);
     }
     if(invMasspiKpKp>0. && invMasspKpiKp>0.){
      if(invMasspiKpKp>0.) fhMassPtGreater3KpTC->Fill(invMasspiKpKp,0.5);
      if(invMasspKpiKp>0.) fhMassPtGreater3KpTC->Fill(invMasspKpiKp,0.5);
     }else{
      if(invMasspiKpKp>0.) fhMassPtGreater3KpTC->Fill(invMasspiKpKp);
      if(invMasspKpiKp>0.) fhMassPtGreater3KpTC->Fill(invMasspKpiKp);
     }
     }
    }
   if(part->Pt()>2.){
     if(invMasspiKp>0. && invMasspKpi>0.){
      if(invMasspiKp>0.) fhMassPtGreater2->Fill(invMasspiKp,0.5);
      if(invMasspKpi>0.) fhMassPtGreater2->Fill(invMasspKpi,0.5);
     }else{
      if(invMasspiKp>0.) fhMassPtGreater2->Fill(invMasspiKp);
      if(invMasspKpi>0.) fhMassPtGreater2->Fill(invMasspKpi);
     }
     if(invMasspiKpLpi>0. && invMasspKpiLpi>0.){
      if(invMasspiKpLpi>0.) fhMassPtGreater2Lpi->Fill(invMasspiKpLpi,0.5);
      if(invMasspKpiLpi>0.) fhMassPtGreater2Lpi->Fill(invMasspKpiLpi,0.5);
     }else{
      if(invMasspiKpLpi>0.) fhMassPtGreater2Lpi->Fill(invMasspiKpLpi);
      if(invMasspKpiLpi>0.) fhMassPtGreater2Lpi->Fill(invMasspKpiLpi);
     }
     if(invMasspiKpKp>0. && invMasspKpiKp>0.){
      if(invMasspiKpKp>0.) fhMassPtGreater2Kp->Fill(invMasspiKpKp,0.5);
      if(invMasspKpiKp>0.) fhMassPtGreater2Kp->Fill(invMasspKpiKp,0.5);
     }else{
      if(invMasspiKpKp>0.) fhMassPtGreater2Kp->Fill(invMasspiKpKp);
      if(invMasspKpiKp>0.) fhMassPtGreater2Kp->Fill(invMasspKpiKp);
     }
     if(passTightCuts>0){
      if(invMasspiKp>0. && invMasspKpi>0.){
       if(invMasspiKp>0.) fhMassPtGreater2TC->Fill(invMasspiKp,0.5);
       if(invMasspKpi>0.) fhMassPtGreater2TC->Fill(invMasspKpi,0.5);
      }else{
       if(invMasspiKp>0.) fhMassPtGreater2TC->Fill(invMasspiKp);
       if(invMasspKpi>0.) fhMassPtGreater2TC->Fill(invMasspKpi);
      }
     if(invMasspiKpLpi>0. && invMasspKpiLpi>0.){
      if(invMasspiKpLpi>0.) fhMassPtGreater2LpiTC->Fill(invMasspiKpLpi,0.5);
      if(invMasspKpiLpi>0.) fhMassPtGreater2LpiTC->Fill(invMasspKpiLpi,0.5);
     }else{
      if(invMasspiKpLpi>0.) fhMassPtGreater2LpiTC->Fill(invMasspiKpLpi);
      if(invMasspKpiLpi>0.) fhMassPtGreater2LpiTC->Fill(invMasspKpiLpi);
     }
     if(invMasspiKpKp>0. && invMasspKpiKp>0.){
      if(invMasspiKpKp>0.) fhMassPtGreater2KpTC->Fill(invMasspiKpKp,0.5);
      if(invMasspKpiKp>0.) fhMassPtGreater2KpTC->Fill(invMasspKpiKp,0.5);
     }else{
      if(invMasspiKpKp>0.) fhMassPtGreater2KpTC->Fill(invMasspiKpKp);
      if(invMasspKpiKp>0.) fhMassPtGreater2KpTC->Fill(invMasspKpiKp);
     }
     }
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
   if(invMasspiKpLpi>0. && invMasspKpiLpi>0.){
    if(invMasspiKpLpi>0.) fMassHistLpi[index]->Fill(invMasspiKpLpi,0.5);
    if(invMasspKpiLpi>0.) fMassHistLpi[index]->Fill(invMasspKpiLpi,0.5);
   }else{
    if(invMasspiKpLpi>0.) fMassHistLpi[index]->Fill(invMasspiKpLpi);
    if(invMasspKpiLpi>0.) fMassHistLpi[index]->Fill(invMasspKpiLpi);
   }
   if(invMasspiKpKp>0. && invMasspKpiKp>0.){
    if(invMasspiKpKp>0.) fMassHistKp[index]->Fill(invMasspiKpKp,0.5);
    if(invMasspKpiKp>0.) fMassHistKp[index]->Fill(invMasspKpiKp,0.5);
   }else{
    if(invMasspiKpKp>0.) fMassHistKp[index]->Fill(invMasspiKpKp);
    if(invMasspKpiKp>0.) fMassHistKp[index]->Fill(invMasspKpiKp);
   }
   if(passTightCuts>0){
    if(invMasspiKp>0. && invMasspKpi>0. && passTightCuts==3){
     if(invMasspiKp>0.) fMassHistTC[index]->Fill(invMasspiKp,0.5);
     if(invMasspKpi>0.) fMassHistTC[index]->Fill(invMasspKpi,0.5);
    }else{
     if(invMasspiKp>0. && passTightCuts==2) fMassHistTC[index]->Fill(invMasspiKp);
     if(invMasspKpi>0. && passTightCuts==1) fMassHistTC[index]->Fill(invMasspKpi);
    }
   if(invMasspiKpLpi>0. && invMasspKpiLpi>0. && passTightCuts==3){
    if(invMasspiKpLpi>0.) fMassHistLpiTC[index]->Fill(invMasspiKpLpi,0.5);
    if(invMasspKpiLpi>0.) fMassHistLpiTC[index]->Fill(invMasspKpiLpi,0.5);
   }else{
    if(invMasspiKpLpi>0. && passTightCuts==2) fMassHistLpiTC[index]->Fill(invMasspiKpLpi);
    if(invMasspKpiLpi>0.&& passTightCuts==1) fMassHistLpiTC[index]->Fill(invMasspKpiLpi);
   }
   if(invMasspiKpKp>0. && invMasspKpiKp>0. && passTightCuts==3){
    if(invMasspiKpKp>0.) fMassHistKpTC[index]->Fill(invMasspiKpKp,0.5);
    if(invMasspKpiKp>0.) fMassHistKpTC[index]->Fill(invMasspKpiKp,0.5);
   }else{
    if(invMasspiKpKp>0. && passTightCuts==2) fMassHistKpTC[index]->Fill(invMasspiKpKp);
    if(invMasspKpiKp>0.&& passTightCuts==1) fMassHistKpTC[index]->Fill(invMasspKpiKp);
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
    if(invMasspiKpLpi>0. && invMasspKpiLpi>0.){
     if(invMasspiKpLpi>0.) fMassHistLpi[index]->Fill(invMasspiKpLpi,0.5);
     if(invMasspKpiLpi>0.) fMassHistLpi[index]->Fill(invMasspKpiLpi,0.5);
    }else{
     if(invMasspiKpLpi>0.) fMassHistLpi[index]->Fill(invMasspiKpLpi);
     if(invMasspKpiLpi>0.) fMassHistLpi[index]->Fill(invMasspKpiLpi);
    }
    if(invMasspiKpKp>0. && invMasspKpiKp>0.){
     if(invMasspiKpKp>0.) fMassHistKp[index]->Fill(invMasspiKpKp,0.5);
     if(invMasspKpiKp>0.) fMassHistKp[index]->Fill(invMasspKpiKp,0.5);
    }else{
     if(invMasspiKpKp>0.) fMassHistKp[index]->Fill(invMasspiKpKp);
     if(invMasspKpiKp>0.) fMassHistKp[index]->Fill(invMasspKpiKp);
    }
     if(passTightCuts>0){
      if(invMasspiKp>0. && invMasspKpi>0. && passTightCuts==3){
       if(invMasspiKp>0.) fMassHistTC[index]->Fill(invMasspiKp,0.5);
       if(invMasspKpi>0.) fMassHistTC[index]->Fill(invMasspKpi,0.5);
      }else{
       if(invMasspiKp>0. && passTightCuts==2) fMassHistTC[index]->Fill(invMasspiKp);
       if(invMasspKpi>0.&& passTightCuts==1) fMassHistTC[index]->Fill(invMasspKpi);
      }
    if(invMasspiKpLpi>0. && invMasspKpiLpi>0. && passTightCuts==3){
     if(invMasspiKpLpi>0.) fMassHistLpiTC[index]->Fill(invMasspiKpLpi,0.5);
     if(invMasspKpiLpi>0.) fMassHistLpiTC[index]->Fill(invMasspKpiLpi,0.5);
    }else{
     if(invMasspiKpLpi>0. && passTightCuts==2) fMassHistLpiTC[index]->Fill(invMasspiKpLpi);
     if(invMasspKpiLpi>0.&& passTightCuts==1) fMassHistLpiTC[index]->Fill(invMasspKpiLpi);
    } 
    if(invMasspiKpKp>0. && invMasspKpiKp>0. && passTightCuts==3){
     if(invMasspiKpKp>0.) fMassHistKpTC[index]->Fill(invMasspiKpKp,0.5);
     if(invMasspKpiKp>0.) fMassHistKpTC[index]->Fill(invMasspKpiKp,0.5);
    }else{
     if(invMasspiKpKp>0. && passTightCuts==2) fMassHistKpTC[index]->Fill(invMasspiKpKp);
     if(invMasspKpiKp>0.&& passTightCuts==1) fMassHistKpTC[index]->Fill(invMasspKpiKp);
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
    if(invMasspiKpLpi>0. && invMasspKpiLpi>0.){
     if(invMasspiKpLpi>0.) fMassHistLpi[index]->Fill(invMasspiKpLpi,0.5);
     if(invMasspKpiLpi>0.) fMassHistLpi[index]->Fill(invMasspKpiLpi,0.5);
    }else{
     if(invMasspiKpLpi>0.) fMassHistLpi[index]->Fill(invMasspiKpLpi);
     if(invMasspKpiLpi>0.) fMassHistLpi[index]->Fill(invMasspKpiLpi);
    }
    if(invMasspiKpKp>0. && invMasspKpiKp>0.){
     if(invMasspiKpKp>0.) fMassHistKp[index]->Fill(invMasspiKpKp,0.5);
     if(invMasspKpiKp>0.) fMassHistKp[index]->Fill(invMasspKpiKp,0.5);
    }else{
     if(invMasspiKpKp>0.) fMassHistKp[index]->Fill(invMasspiKpKp);
     if(invMasspKpiKp>0.) fMassHistKp[index]->Fill(invMasspKpiKp);
    }
    if(invMasspiKp>0. && invMasspKpi>0. && passTightCuts==3){
     fMassHistTC[index]->Fill(invMasspiKp,0.5);
     fMassHistTC[index]->Fill(invMasspKpi,0.5);
    }else{
     if(invMasspiKp>0. && passTightCuts==2) fMassHistTC[index]->Fill(invMasspiKp);
     if(invMasspKpi>0. && passTightCuts==1) fMassHistTC[index]->Fill(invMasspKpi);
    }
    if(invMasspiKpLpi>0. && invMasspKpiLpi>0. && passTightCuts==3){
     if(invMasspiKpLpi>0.) fMassHistLpiTC[index]->Fill(invMasspiKpLpi,0.5);
     if(invMasspKpiLpi>0.) fMassHistLpiTC[index]->Fill(invMasspKpiLpi,0.5);
    }else{
     if(invMasspiKpLpi>0. && passTightCuts==2) fMassHistLpiTC[index]->Fill(invMasspiKpLpi);
     if(invMasspKpiLpi>0.&& passTightCuts==1) fMassHistLpiTC[index]->Fill(invMasspKpiLpi);
    } 
    if(invMasspiKpKp>0. && invMasspKpiKp>0. && passTightCuts==3){
     if(invMasspiKpKp>0.) fMassHistKpTC[index]->Fill(invMasspiKpKp,0.5);
     if(invMasspKpiKp>0.) fMassHistKpTC[index]->Fill(invMasspKpiKp,0.5);
    }else{
     if(invMasspiKpKp>0. && passTightCuts==2) fMassHistKpTC[index]->Fill(invMasspiKpKp);
     if(invMasspKpiKp>0.&& passTightCuts==1) fMassHistKpTC[index]->Fill(invMasspKpiKp);
    }
   }

  }
 }
}
 return;
}
//-----------------------
void AliAnalysisTaskSELambdac::FillVarHists(AliAODRecoDecayHF3Prong *part, TClonesArray *arrMC, AliRDHFCutsLctopKpi *cuts, TList *listout,AliAODEvent* aod){
  //
  // function used in UserExec to fill variable histograms analysing MC
  //


  Int_t pdgDgLctopKpi[3]={2212,321,211};
  Int_t lab=-9999;
  if(fReadMC) lab=part->MatchToMC(4122,arrMC,3,pdgDgLctopKpi); //return MC particle label if the array corresponds to a Lc, -1 if not (cf. AliAODRecoDecay.cxx)
  Int_t isSelectedPID=cuts->IsSelectedPID(part); //0 rejected,1 Lc -> p K- pi+,2 Lc -> pi+ K- p, (K at center because different sign is there)
                                                 //3 both (it should never happen...)

  if (isSelectedPID==0)fNentries->Fill(7);
  if (isSelectedPID==1)fNentries->Fill(8);
  if (isSelectedPID==2)fNentries->Fill(9);
  if (isSelectedPID==3)fNentries->Fill(10);


  Double_t mPDG=TDatabasePDG::Instance()->GetParticle(4122)->Mass();
  Double_t minvLcpKpi = part->InvMassLcpKpi();
  Double_t minvLcpiKp = part->InvMassLcpiKp();

  Double_t invmasscut=0.05;

  TString fillthis="";


  AliAODTrack *prong0=(AliAODTrack*)part->GetDaughter(0);
  AliAODTrack *prong1=(AliAODTrack*)part->GetDaughter(1);
  AliAODTrack *prong2=(AliAODTrack*)part->GetDaughter(2);
  if (!prong0 || !prong1 || !prong2) {
    fNentries->Fill(5);
    return;
  }

  //check pdg of the prongs
  Int_t labprong[3];
  if(fReadMC){
    labprong[0]=prong0->GetLabel();
    labprong[1]=prong1->GetLabel();
    labprong[2]=prong2->GetLabel();
  }
  AliAODMCParticle *mcprong=0;

  Int_t pdgProng[3]={0,0,0};
  Int_t pdgProngID[3]={0,0,0};
    if(fReadMC) {
     for (Int_t iprong=0;iprong<3;iprong++){
      if(labprong[iprong]<0) continue;
      mcprong= (AliAODMCParticle*)arrMC->At(labprong[iprong]);
      pdgProng[iprong]=TMath::Abs(mcprong->GetPdgCode());
     }
     if(isSelectedPID>0 && fReadMC){
      pdgProngID[1]=321;
      if(isSelectedPID==1) {pdgProngID[0]=2212;pdgProngID[2]=211;}
      if(isSelectedPID==2) {pdgProngID[0]=211;pdgProngID[2]=2212;}
     }
    } else {
      if (isSelectedPID>0) { 
	pdgProng[1]=321;
	if(isSelectedPID==1) {pdgProng[0]=2212;pdgProng[2]=211;}
	if(isSelectedPID==2) {pdgProng[0]=211;pdgProng[2]=2212;}
      }
    }

 Int_t selection=cuts->IsSelected(part,AliRDHFCuts::kCandidate,aod);

  if((lab>=0 && fReadMC) || (!fReadMC && (isSelectedPID>0)) ){ //signal (from MC or PID)

    fillthis="hMass";

     if ((fReadMC && ((AliAODMCParticle*)arrMC->At(lab))->GetPdgCode() == 4122)
	|| (!fReadMC && (isSelectedPID>0 && part->Charge()>0))){    //Lc
	 if ( (pdgProng[0]==2212) && (pdgProng[1]==321) && (pdgProng[2]==211) ) ((TH1F*)listout->FindObject(fillthis))->Fill(minvLcpKpi);
	 else if ( (pdgProng[0]==211) && (pdgProng[1]==321) && (pdgProng[2]==2212) ) ((TH1F*)listout->FindObject(fillthis))->Fill(minvLcpiKp);
    }
    else if ((fReadMC && ((AliAODMCParticle*)arrMC->At(lab))->GetPdgCode() == -4122)
	  || (!fReadMC && (isSelectedPID>0 && part->Charge()<0))){ //anti-Lc
      if ( (pdgProng[0]==2212) && (pdgProng[1]==321) && (pdgProng[2]==211) ) ((TH1F*)listout->FindObject(fillthis))->Fill(minvLcpKpi);
      else if ( (pdgProng[0]==211) && (pdgProng[1]==321) && (pdgProng[2]==2212) ) ((TH1F*)listout->FindObject(fillthis))->Fill(minvLcpiKp);
    }


    //apply cut on invariant mass on the pair
    if(selection>0){

      Double_t ptmax=0.;
      for (Int_t iprong=0; iprong<3; iprong++) {
	AliAODTrack *prong=(AliAODTrack*)part->GetDaughter(iprong);
	if (fReadMC) {
	 labprong[iprong]=prong->GetLabel();
	 if(pdgProngID[iprong]==2212) {
	  fillthis="hpIDTot";
	  ((TH1F*)listout->FindObject(fillthis))->Fill(part->PtProng(iprong));
	 }
	 if(pdgProngID[iprong]==321) {
	  fillthis="hKIDTot";
	  ((TH1F*)listout->FindObject(fillthis))->Fill(part->PtProng(iprong));
	 }
	 if(pdgProngID[iprong]==211) {
	  fillthis="hpiIDTot";
	  ((TH1F*)listout->FindObject(fillthis))->Fill(part->PtProng(iprong));
	}
       }
	AliAODPid *pidObjtrk=prong->GetDetPid();
        AliAODPidHF *pidObj=new AliAODPidHF();
	Bool_t hasTOF=pidObj->CheckStatus(prong,"TOF");
	Bool_t hasTPC=pidObj->CheckStatus(prong,"TPC");
	delete pidObj;
	Double_t tofSignal=0.;
	Double_t dedxTPC=0.;
	Double_t momTOF=0.;
	Double_t momTPC=0.;
	if(hasTOF) {
	 momTOF = prong->P();
	 tofSignal=pidObjtrk->GetTOFsignal();
	}
	if(hasTPC) {
	 momTPC = pidObjtrk->GetTPCmomentum();
	 dedxTPC=pidObjtrk->GetTPCsignal();
	}
	switch (pdgProng[iprong]) {
	case 2212:
	    fillthis="hpTOFSignal";
	    ((TH1F*)listout->FindObject(fillthis))->Fill(tofSignal);
	    fillthis="hpTPCSignal";
	    ((TH1F*)listout->FindObject(fillthis))->Fill(dedxTPC);
	    fillthis="hpptProng";
	    ((TH1F*)listout->FindObject(fillthis))->Fill(part->PtProng(iprong));
	    fillthis="hpd0Prong";
	    ((TH1F*)listout->FindObject(fillthis))->Fill(part->Getd0Prong(iprong));
	    fillthis="hpSignalVspTPC";
	    ((TH2F*)listout->FindObject(fillthis))->Fill(momTPC,dedxTPC);
	    fillthis="hpSignalVspTOF";
	    ((TH2F*)listout->FindObject(fillthis))->Fill(momTOF,tofSignal);
	    ((TH1F*)listout->FindObject(fillthis))->Fill(part->Getd0Prong(iprong));
	    AliPID::EParticleType typep;
	    typep=AliPID::EParticleType(4);
	    if(hasTPC) {
	     Double_t nsigmap = fUtilPid->NumberOfSigmasTPC(prong,typep);
	     fillthis="hpSigmaVspTPC";
	     ((TH2F*)listout->FindObject(fillthis))->Fill(momTPC,nsigmap);
	    }
	    if(hasTOF){
	     Double_t nsigma=fUtilPid->NumberOfSigmasTOF(prong,typep);
	     fillthis="hpSigmaVspTOF";
	     ((TH2F*)listout->FindObject(fillthis))->Fill(momTOF,nsigma);
	    }
	    if(fReadMC){
              // real protons
	     fillthis="hpRealTot";
	     ((TH1F*)listout->FindObject(fillthis))->Fill(part->PtProng(iprong));
	      // id protons
	      if(pdgProngID[iprong]==2212) {
	       fillthis="hpIDGood";
	       ((TH1F*)listout->FindObject(fillthis))->Fill(part->PtProng(iprong));
	      }
	      // misidentified kaons, pions
	      if(pdgProngID[iprong]==321) {
	       fillthis="hKnonIDTot";
	       ((TH1F*)listout->FindObject(fillthis))->Fill(part->PtProng(iprong));
	      }
	      if(pdgProngID[iprong]==211) {
	       fillthis="hpinonIDTot";
	       ((TH1F*)listout->FindObject(fillthis))->Fill(part->PtProng(iprong));
	      }
	    }

	   break;
	case 321:
	    fillthis="hKTOFSignal";
	    ((TH1F*)listout->FindObject(fillthis))->Fill(tofSignal);
	    fillthis="hKTPCSignal";
	    ((TH1F*)listout->FindObject(fillthis))->Fill(dedxTPC);
	    fillthis="hKptProng";
	    ((TH1F*)listout->FindObject(fillthis))->Fill(part->PtProng(iprong));
	    fillthis="hKd0Prong";
	    ((TH1F*)listout->FindObject(fillthis))->Fill(part->Getd0Prong(iprong));
	    fillthis="hKSignalVspTPC";
	    ((TH2F*)listout->FindObject(fillthis))->Fill(momTPC,dedxTPC);
	    fillthis="hKSignalVspTOF";
	    ((TH2F*)listout->FindObject(fillthis))->Fill(momTOF,tofSignal);
	    AliPID::EParticleType typek;
	    typek=AliPID::EParticleType(3);
	    if(hasTPC) {
	     Double_t nsigmap = fUtilPid->NumberOfSigmasTPC(prong,typek);
	     fillthis="hKSigmaVspTPC";
	     ((TH2F*)listout->FindObject(fillthis))->Fill(momTPC,nsigmap);
	    }
	    if(hasTOF){
	     Double_t nsigma = fUtilPid->NumberOfSigmasTOF(prong,typek);
	     fillthis="hKSigmaVspTOF";
	     ((TH2F*)listout->FindObject(fillthis))->Fill(momTOF,nsigma);
	    }

	    if(fReadMC){
              // real kaons
	     fillthis="hKRealTot";
	     ((TH1F*)listout->FindObject(fillthis))->Fill(part->PtProng(iprong));
	      // id kaons
	      if(pdgProngID[iprong]==321) {
	       fillthis="hKIDGood";
	       ((TH1F*)listout->FindObject(fillthis))->Fill(part->PtProng(iprong));
	      }
	      // misidentified protons, pions
	      if(pdgProngID[iprong]==2212) {
	       fillthis="hpnonIDTot";
	       ((TH1F*)listout->FindObject(fillthis))->Fill(part->PtProng(iprong));
	      }
	      if(pdgProngID[iprong]==211) {
	       fillthis="hpinonIDTot";
	       ((TH1F*)listout->FindObject(fillthis))->Fill(part->PtProng(iprong));
	      }
	    }
	    break;
	case 211:
	    fillthis="hpiTOFSignal";
	    ((TH1F*)listout->FindObject(fillthis))->Fill(tofSignal);
	    fillthis="hpiTPCSignal";
	    ((TH1F*)listout->FindObject(fillthis))->Fill(dedxTPC);
	    fillthis="hpiptProng";
	    ((TH1F*)listout->FindObject(fillthis))->Fill(part->PtProng(iprong));
	    fillthis="hpid0Prong";
	    ((TH1F*)listout->FindObject(fillthis))->Fill(part->Getd0Prong(iprong));
	    fillthis="hpiSignalVspTPC";
	    ((TH2F*)listout->FindObject(fillthis))->Fill(momTPC,dedxTPC);
	    fillthis="hpiSignalVspTOF";
	    ((TH2F*)listout->FindObject(fillthis))->Fill(momTOF,tofSignal);
	    AliPID::EParticleType typepi;
	    typepi=AliPID::EParticleType(2);
	    if(hasTPC) {
	     Double_t nsigmap = fUtilPid->NumberOfSigmasTPC(prong,typepi);
	    fillthis="hpiSigmaVspTPC";
	    ((TH2F*)listout->FindObject(fillthis))->Fill(momTPC,nsigmap);
	    }
	    if(hasTOF){
	     Double_t nsigma = fUtilPid->NumberOfSigmasTOF(prong,typepi);
	     fillthis="hpiSigmaVspTOF";
	     ((TH2F*)listout->FindObject(fillthis))->Fill(momTOF,nsigma);
	    }

	    if(fReadMC){
              // real pions
	     fillthis="hpiRealTot";
	     ((TH1F*)listout->FindObject(fillthis))->Fill(part->PtProng(iprong));
	      // id pions
	      if(pdgProngID[iprong]==211) {
	       fillthis="hpiIDGood";
	       ((TH1F*)listout->FindObject(fillthis))->Fill(part->PtProng(iprong));
	      }
	      // misidentified protons, kaons
	      if(pdgProngID[iprong]==2212) {
	       fillthis="hpnonIDTot";
	       ((TH1F*)listout->FindObject(fillthis))->Fill(part->PtProng(iprong));
	      }
	      if(pdgProngID[iprong]==321) {
	       fillthis="hKnonIDTot";
	       ((TH1F*)listout->FindObject(fillthis))->Fill(part->PtProng(iprong));
	      }
	    }
	    break;

	default:
	  break;
	}
	if(part->PtProng(iprong)>ptmax)ptmax=part->PtProng(iprong);

      } //end loop on prongs
	//now histograms where we don't need to check identity
      fillthis = "hDist12toPrim";
      ((TH1F*)listout->FindObject(fillthis))->Fill(part->GetDist12toPrim());
      ((TH1F*)listout->FindObject(fillthis))->Fill(part->GetDist23toPrim());
      fillthis = "hSigmaVert";
      ((TH1F*)listout->FindObject(fillthis))->Fill(part->GetSigmaVert());
      fillthis = "hDCAs";
      Double_t dcas[3];
      part->GetDCAs(dcas);
      for (Int_t idca=0;idca<3;idca++) ((TH1F*)listout->FindObject(fillthis))->Fill(dcas[idca]);
      fillthis = "hCosPointingAngle";
      ((TH1F*)listout->FindObject(fillthis))->Fill(part->CosPointingAngle());
      fillthis = "hDecayLength";
      ((TH1F*)listout->FindObject(fillthis))->Fill(part->DecayLength());
      Double_t sum2=part->Getd0Prong(0)*part->Getd0Prong(0)+part->Getd0Prong(1)*part->Getd0Prong(1)+part->Getd0Prong(2)*part->Getd0Prong(2);
      fillthis = "hSum2";
      ((TH1F*)listout->FindObject(fillthis))->Fill(sum2);
      fillthis = "hptmax";
      ((TH1F*)listout->FindObject(fillthis))->Fill(ptmax);
      fillthis="hLcpt";
      ((TH1F*)listout->FindObject(fillthis))->Fill(part->Pt());

    } //end mass cut


  } else if( (lab<0) && fReadMC) { // background **ONLY MC**

    fillthis="hbMass";

     if (part->Charge()>0){    //Lc background
	 if ( (pdgProng[0]==2212) && (pdgProng[1]==321) && (pdgProng[2]==211) ) ((TH1F*)listout->FindObject(fillthis))->Fill(minvLcpKpi);
	 else if ( (pdgProng[0]==211) && (pdgProng[1]==321) && (pdgProng[2]==2212) ) ((TH1F*)listout->FindObject(fillthis))->Fill(minvLcpiKp);
     }
     else if (part->Charge()<0){ //anti-Lc background
      if ( (pdgProng[0]==2212) && (pdgProng[1]==321) && (pdgProng[2]==211) ) ((TH1F*)listout->FindObject(fillthis))->Fill(minvLcpKpi);
      else if ( (pdgProng[0]==211) && (pdgProng[1]==321) && (pdgProng[2]==2212) ) ((TH1F*)listout->FindObject(fillthis))->Fill(minvLcpiKp);
    }


    //apply cut on invariant mass on the pair
    if(TMath::Abs(minvLcpKpi-mPDG)<invmasscut || TMath::Abs(minvLcpiKp-mPDG)<invmasscut){
    if(selection>0){


      Double_t ptmax=0.;

      for (Int_t iprong=0; iprong<3; iprong++) {
	AliAODTrack *prong=(AliAODTrack*)part->GetDaughter(iprong);
	if (fReadMC) {
	labprong[iprong]=prong->GetLabel();
	 if(pdgProngID[iprong]==2212) {
	  fillthis="hbpIDTot";
	  ((TH1F*)listout->FindObject(fillthis))->Fill(part->PtProng(iprong));
	 }
	 if(pdgProngID[iprong]==321) {
	  fillthis="hbKIDTot";
	  ((TH1F*)listout->FindObject(fillthis))->Fill(part->PtProng(iprong));
	 }
	 if(pdgProngID[iprong]==211) {
	  fillthis="hbpiIDTot";
	  ((TH1F*)listout->FindObject(fillthis))->Fill(part->PtProng(iprong));
	 }
	}
	AliAODPid *pidObjtrk=prong->GetDetPid();
        AliAODPidHF *pidObj=new AliAODPidHF();
	Bool_t hasTOF=pidObj->CheckStatus(prong,"TOF");
	Bool_t hasTPC=pidObj->CheckStatus(prong,"TPC");
	delete pidObj;
	Double_t tofSignal=0.;
	Double_t dedxTPC=0.;
	Double_t momTOF=0.;
	Double_t momTPC=0.;
	if(hasTOF) {
	 momTOF = prong->P();
	 tofSignal=pidObjtrk->GetTOFsignal();
	}
	if(hasTPC) {
	 momTPC = pidObjtrk->GetTPCmomentum();
	 dedxTPC=pidObjtrk->GetTPCsignal();
	}

	switch (pdgProng[iprong]) {
	case 2212:
	    fillthis="hbpTOFSignal";
	    ((TH1F*)listout->FindObject(fillthis))->Fill(tofSignal);
	    fillthis="hbpTPCSignal";
	    ((TH1F*)listout->FindObject(fillthis))->Fill(dedxTPC);
	    fillthis="hbpptProng";
	    ((TH1F*)listout->FindObject(fillthis))->Fill(part->PtProng(iprong));
	    fillthis="hbpd0Prong";
	    ((TH1F*)listout->FindObject(fillthis))->Fill(part->Getd0Prong(iprong));
	    fillthis="hbpSignalVspTPC";
	    ((TH2F*)listout->FindObject(fillthis))->Fill(momTPC,dedxTPC);
	    fillthis="hbpSignalVspTOF";
	    ((TH2F*)listout->FindObject(fillthis))->Fill(momTOF,tofSignal);
	    AliPID::EParticleType typep;
	    typep=AliPID::EParticleType(4);
	    if(hasTPC) {
	     Double_t nsigmap = fUtilPid->NumberOfSigmasTPC(prong,typep);
	     fillthis="hbpSigmaVspTPC";
	     ((TH2F*)listout->FindObject(fillthis))->Fill(momTPC,nsigmap);
	    }
	    if(hasTOF){
	     Double_t nsigma = fUtilPid->NumberOfSigmasTOF(prong,typep);
	     fillthis="hbpSigmaVspTOF";
	     ((TH2F*)listout->FindObject(fillthis))->Fill(momTOF,nsigma);
	    }
	    if(fReadMC){
              // real protons
	     fillthis="hbpRealTot";
	     ((TH1F*)listout->FindObject(fillthis))->Fill(part->PtProng(iprong));
	      // id protons
	      if(pdgProngID[iprong]==2212) {
	       fillthis="hbpIDGood";
	       ((TH1F*)listout->FindObject(fillthis))->Fill(part->PtProng(iprong));
	      }
	      // misidentified kaons, pions
	      if(pdgProngID[iprong]==321) {
	       fillthis="hbKnonIDTot";
	       ((TH1F*)listout->FindObject(fillthis))->Fill(part->PtProng(iprong));
	      }
	      if(pdgProngID[iprong]==211) {
	       fillthis="hbpinonIDTot";
	       ((TH1F*)listout->FindObject(fillthis))->Fill(part->PtProng(iprong));
	      }
	    }
	    break;
	case 321:
	    fillthis="hbKTOFSignal";
	    ((TH1F*)listout->FindObject(fillthis))->Fill(tofSignal);
	    fillthis="hbKTPCSignal";
	    ((TH1F*)listout->FindObject(fillthis))->Fill(dedxTPC);
	    fillthis="hbKptProng";
	    ((TH1F*)listout->FindObject(fillthis))->Fill(part->PtProng(iprong));
	    fillthis="hbKd0Prong";
	    ((TH1F*)listout->FindObject(fillthis))->Fill(part->Getd0Prong(iprong));
	    fillthis="hbKSignalVspTPC";
	    ((TH2F*)listout->FindObject(fillthis))->Fill(momTPC,dedxTPC);
	    fillthis="hbKSignalVspTOF";
	    ((TH2F*)listout->FindObject(fillthis))->Fill(momTOF,tofSignal);
	    AliPID::EParticleType typek;
	    typek=AliPID::EParticleType(3);
	    if(hasTPC) {
	     Double_t nsigmap = fUtilPid->NumberOfSigmasTPC(prong,typek);
	     fillthis="hbKSigmaVspTPC";
	     ((TH2F*)listout->FindObject(fillthis))->Fill(momTPC,nsigmap);
	    }
	    if(hasTOF){
	     Double_t nsigma = fUtilPid->NumberOfSigmasTOF(prong,typek);
	     fillthis="hbKSigmaVspTOF";
	     ((TH2F*)listout->FindObject(fillthis))->Fill(momTOF,nsigma);
	    }
	    if(fReadMC){
              // real kaons
	     fillthis="hbKRealTot";
	     ((TH1F*)listout->FindObject(fillthis))->Fill(part->PtProng(iprong));
	      // id kaons
	      if(pdgProngID[iprong]==321) {
	       fillthis="hbKIDGood";
	       ((TH1F*)listout->FindObject(fillthis))->Fill(part->PtProng(iprong));
	      }
	      // misidentified protons, pions
	      if(pdgProngID[iprong]==2212) {
	       fillthis="hbpnonIDTot";
	       ((TH1F*)listout->FindObject(fillthis))->Fill(part->PtProng(iprong));
	      }
	      if(pdgProngID[iprong]==211) {
	       fillthis="hbpinonIDTot";
	       ((TH1F*)listout->FindObject(fillthis))->Fill(part->PtProng(iprong));
	      }
	    }
	    break;
	case 211:
	    fillthis="hbpiTOFSignal";
	    ((TH1F*)listout->FindObject(fillthis))->Fill(tofSignal);
	    fillthis="hbpiTPCSignal";
	    ((TH1F*)listout->FindObject(fillthis))->Fill(dedxTPC);
	    fillthis="hbpiptProng";
	    ((TH1F*)listout->FindObject(fillthis))->Fill(part->PtProng(iprong));
	    fillthis="hbpid0Prong";
	    ((TH1F*)listout->FindObject(fillthis))->Fill(part->Getd0Prong(iprong));
	    fillthis="hbpiSignalVspTPC";
	    ((TH2F*)listout->FindObject(fillthis))->Fill(momTPC,dedxTPC);
	    fillthis="hbpiSignalVspTOF";
	    ((TH2F*)listout->FindObject(fillthis))->Fill(momTOF,tofSignal);
	    AliPID::EParticleType typepi;
	    typepi=AliPID::EParticleType(2);
	    if(hasTPC) {
	     Double_t nsigmap = fUtilPid->NumberOfSigmasTPC(prong,typepi);
	     fillthis="hbpiSigmaVspTPC";
	     ((TH2F*)listout->FindObject(fillthis))->Fill(momTPC,nsigmap);
	    }
	    if(hasTOF){
	     Double_t nsigma = fUtilPid->NumberOfSigmasTOF(prong,typepi);
	     fillthis="hbpiSigmaVspTOF";
	     ((TH2F*)listout->FindObject(fillthis))->Fill(momTOF,nsigma);
	    }
	    if(fReadMC){
              // real pions
	     fillthis="hbpiRealTot";
	     ((TH1F*)listout->FindObject(fillthis))->Fill(part->PtProng(iprong));
	      // id pions
	      if(pdgProngID[iprong]==211) {
	       fillthis="hbpiIDGood";
	       ((TH1F*)listout->FindObject(fillthis))->Fill(part->PtProng(iprong));
	      }
	      // misidentified protons, kaons
	      if(pdgProngID[iprong]==2212) {
	       fillthis="hbpnonIDTot";
	       ((TH1F*)listout->FindObject(fillthis))->Fill(part->PtProng(iprong));
	      }
	      if(pdgProngID[iprong]==321) {
	       fillthis="hbKnonIDTot";
	       ((TH1F*)listout->FindObject(fillthis))->Fill(part->PtProng(iprong));
	      }
	    }
	    break;
	default:
	  break;
	}
	if(part->PtProng(iprong)>ptmax)ptmax=part->PtProng(iprong);

      } //end loop on prongs
	//now histograms where we don't need to check identity
	fillthis="hbLcpt";
	((TH1F*)listout->FindObject(fillthis))->Fill(part->Pt());
	fillthis = "hbDist12toPrim";
	((TH1F*)listout->FindObject(fillthis))->Fill(part->GetDist12toPrim());
	((TH1F*)listout->FindObject(fillthis))->Fill(part->GetDist23toPrim());
	fillthis = "hbSigmaVert";
	((TH1F*)listout->FindObject(fillthis))->Fill(part->GetSigmaVert());
	fillthis = "hbDCAs";
	Double_t dcas[3];
	part->GetDCAs(dcas);
	for (Int_t idca=0;idca<3;idca++) ((TH1F*)listout->FindObject(fillthis))->Fill(dcas[idca]);
	fillthis = "hbCosPointingAngle";
	((TH1F*)listout->FindObject(fillthis))->Fill(part->CosPointingAngle());
	fillthis = "hbDecayLength";
        ((TH1F*)listout->FindObject(fillthis))->Fill(part->DecayLength());
      Double_t sum2=part->Getd0Prong(0)*part->Getd0Prong(0)+part->Getd0Prong(1)*part->Getd0Prong(1)+part->Getd0Prong(2)*part->Getd0Prong(2);
      fillthis = "hbSum2";
      ((TH1F*)listout->FindObject(fillthis))->Fill(sum2);
      fillthis = "hbptmax";
      ((TH1F*)listout->FindObject(fillthis))->Fill(ptmax);

     }

    } //end mass cut

  } // end background case
  return;
}

