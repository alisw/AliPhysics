/**************************************************************************
 * Copyright(c) 1998-2013, ALICE Experiment at CERN, All rights reserved. *
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
// AliAnalysisTaskSE for the extraction of signal(e.g Lambdac) of heavy flavor
// decay candidates with the MC truth. Specific for ITS upgrade studies
// Authors: r.romita@liv.ac.uk
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
#include "AliAnalysisTaskSELambdacUp.h"
#include "AliKFParticle.h"
#include "AliAODPidHF.h"
#include "AliRDHFCutsLctopKpi.h"
#include "AliRDHFCuts.h"
#include "AliKFVertex.h"
#include "AliESDVertex.h"
//#include "AliAODpidUtil.h"
#include "AliAODPid.h"
#include "AliInputEventHandler.h"
#include "AliPID.h"
#include "AliNormalizationCounter.h"

ClassImp(AliAnalysisTaskSELambdacUp)


//________________________________________________________________________
  AliAnalysisTaskSELambdacUp::AliAnalysisTaskSELambdacUp():
    AliAnalysisTaskSE(),
    fOutput(0), 
    fHistNEvents(0),
    fhChi2(0),
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
    fFillVarHists(kFALSE),
    fNentries(0),
    fOutputMC(0),
    fPIDResponse(0),
    fCounter(0)
{
  // Default constructor
  Float_t ptlims[7]={0.,2.,4.,6.,8.,12.,24.};
  SetPtBinLimit(7,ptlims);
   for(Int_t icut=0; icut<2; icut++) fCutsKF[icut]=0.;
   for(Int_t j=0; j<3*kMaxPtBins; j++){
    fMassHist[j]=0x0;
    fMassHistTC[j]=0x0;
    fMassHistLpi[j]=0x0;
    fMassHistLpiTC[j]=0x0;
    fMassHistKp[j]=0x0;
    fMassHistKpTC[j]=0x0;
    fMassHistDk[j]=0x0;
    fMassHistDkTC[j]=0x0;
    fMassHist3Pr[j]=0x0;
    fMassHist3PrTC[j]=0x0;
   }
}

//________________________________________________________________________
AliAnalysisTaskSELambdacUp::AliAnalysisTaskSELambdacUp(const char *name,Bool_t fillNtuple,AliRDHFCutsLctopKpi *lccutsana,AliRDHFCutsLctopKpi *lccutsprod):
  AliAnalysisTaskSE(name),
  fOutput(0),
  fHistNEvents(0),
  fhChi2(0),
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
  fFillVarHists(kFALSE),
  fNentries(0),
  fOutputMC(0),
  fPIDResponse(0),
  fCounter(0)
{
  SetPtBinLimit(fRDCutsAnalysis->GetNPtBins()+1,
		fRDCutsAnalysis->GetPtBinLimits());
   for(Int_t icut=0; icut<2; icut++) fCutsKF[icut]=0.;
   for(Int_t j=0; j<3*kMaxPtBins; j++){
    fMassHist[j]=0x0;
    fMassHistTC[j]=0x0;
    fMassHistLpi[j]=0x0;
    fMassHistLpiTC[j]=0x0;
    fMassHistKp[j]=0x0;
    fMassHistKpTC[j]=0x0;
    fMassHistDk[j]=0x0;
    fMassHistDkTC[j]=0x0;
    fMassHist3Pr[j]=0x0;
    fMassHist3PrTC[j]=0x0;
   }
  // Default constructor
  // Output slot #1 writes into a TList container
  DefineOutput(1,TList::Class());  //My private output
  DefineOutput(2,TList::Class());
  DefineOutput(3,TList::Class());
  DefineOutput(4,TH1F::Class());
  DefineOutput(5,AliNormalizationCounter::Class());
  if (fFillNtuple) {
    // Output slot #2 writes into a TNtuple container
    DefineOutput(6,TNtuple::Class());  //My private output
  }
}

//________________________________________________________________________
AliAnalysisTaskSELambdacUp::~AliAnalysisTaskSELambdacUp()
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
  /*
    if (fUtilPid){
    delete fUtilPid;
    fUtilPid = 0;
    }
  */
  if (fPIDResponse) {
    delete  fPIDResponse;
  }
  if(fCounter){
   delete fCounter;
   fCounter = 0;
  }

}  
//_________________________________________________________________
void  AliAnalysisTaskSELambdacUp::SetMassLimits(Float_t range){
  fUpmasslimit = 2.286+range;
  fLowmasslimit = 2.286-range;
}
//_________________________________________________________________
void  AliAnalysisTaskSELambdacUp::SetMassLimits(Float_t lowlimit, Float_t uplimit){
  if(uplimit>lowlimit)
    {
      fUpmasslimit = lowlimit;
      fLowmasslimit = uplimit;
    }
}


//________________________________________________________________________
void AliAnalysisTaskSELambdacUp::SetPtBinLimit(Int_t n, Float_t* lim){
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
Double_t  AliAnalysisTaskSELambdacUp::GetPtBinLimit(Int_t ibin) const{
  if(ibin>fNPtBins)return -1;
  return fArrayBinLimits[ibin];
} 

//_________________________________________________________________
void AliAnalysisTaskSELambdacUp::Init()
{
  // Initialization

  if (fDebug > 1) printf("AnalysisTaskSELambdac::Init() \n");

  fListCuts=new TList();
  fListCuts->SetOwner();

  fListCuts->Add(new AliRDHFCutsLctopKpi(*fRDCutsAnalysis));
  fListCuts->Add(new AliRDHFCutsLctopKpi(*fRDCutsProduction));
  PostData(2,fListCuts);
  return;
}

//________________________________________________________________________
void AliAnalysisTaskSELambdacUp::UserCreateOutputObjects()
{
  // Create the output container
  //
  if (fDebug > 1) printf("AnalysisTaskSELambdac::UserCreateOutputObjects() \n");

  // Several histograms are more conveniently managed in a TList
  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("OutputHistos");

  TString hisname;
  Int_t index=0;
  for(Int_t i=0;i<fNPtBins;i++){

    index=GetHistoIndex(i);

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
    hisname.Form("hMassPtDk%d",i);
    fMassHistDk[index]=new TH1F(hisname.Data(),hisname.Data(),100,fLowmasslimit,fUpmasslimit);
    fMassHistDk[index]->Sumw2();
    hisname.Form("hMassPtDk%dTC",i);
    fMassHistDkTC[index]=new TH1F(hisname.Data(),hisname.Data(),100,fLowmasslimit,fUpmasslimit);
    fMassHistDkTC[index]->Sumw2();

    hisname.Form("hMassPt3Pr%d",i);
    fMassHist3Pr[index]=new TH1F(hisname.Data(),hisname.Data(),100,fLowmasslimit,fUpmasslimit);
    fMassHist3Pr[index]->Sumw2();
    hisname.Form("hMassPt3Pr%dTC",i);
    fMassHist3PrTC[index]=new TH1F(hisname.Data(),hisname.Data(),100,fLowmasslimit,fUpmasslimit);
    fMassHist3PrTC[index]->Sumw2();
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

    hisname.Form("hSigPtDk%d",i);
    fMassHistDk[index]=new TH1F(hisname.Data(),hisname.Data(),100,fLowmasslimit,fUpmasslimit);
    fMassHistDk[index]->Sumw2();
    hisname.Form("hSigPtDk%dTC",i);
    fMassHistDkTC[index]=new TH1F(hisname.Data(),hisname.Data(),100,fLowmasslimit,fUpmasslimit);
    fMassHistDkTC[index]->Sumw2();

    hisname.Form("hSigPt3Pr%d",i);
    fMassHist3Pr[index]=new TH1F(hisname.Data(),hisname.Data(),100,fLowmasslimit,fUpmasslimit);
    fMassHist3Pr[index]->Sumw2();
    hisname.Form("hSigPt3Pr%dTC",i);
    fMassHist3PrTC[index]=new TH1F(hisname.Data(),hisname.Data(),100,fLowmasslimit,fUpmasslimit);
    fMassHist3PrTC[index]->Sumw2();

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

    hisname.Form("hBkgPtDk%d",i);
    fMassHistDk[index]=new TH1F(hisname.Data(),hisname.Data(),100,fLowmasslimit,fUpmasslimit);
    fMassHistDk[index]->Sumw2();
    hisname.Form("hBkgPtDk%dTC",i);
    fMassHistDkTC[index]=new TH1F(hisname.Data(),hisname.Data(),100,fLowmasslimit,fUpmasslimit);
    fMassHistDkTC[index]->Sumw2();

    hisname.Form("hBkgPt3Pr%d",i);
    fMassHist3Pr[index]=new TH1F(hisname.Data(),hisname.Data(),100,fLowmasslimit,fUpmasslimit);
    fMassHist3Pr[index]->Sumw2();
    hisname.Form("hBkgPt3Pr%dTC",i);
    fMassHist3PrTC[index]=new TH1F(hisname.Data(),hisname.Data(),100,fLowmasslimit,fUpmasslimit);
    fMassHist3PrTC[index]->Sumw2();


    index=GetPYTHIAHistoIndex(i);
    printf("index=%d , i=%d \n",index,i);
    hisname.Form("hPYTHIAPt%d",i);
    fMassHist[index]=new TH1F(hisname.Data(),hisname.Data(),100,fLowmasslimit,fUpmasslimit);
    fMassHist[index]->Sumw2();
    hisname.Form("hPYTHIAPt%dTC",i);
    fMassHistTC[index]=new TH1F(hisname.Data(),hisname.Data(),100,fLowmasslimit,fUpmasslimit);
    fMassHistTC[index]->Sumw2();
    hisname.Form("hPYTHIAPtLpi%d",i);
    fMassHistLpi[index]=new TH1F(hisname.Data(),hisname.Data(),100,fLowmasslimit,fUpmasslimit);
    fMassHistLpi[index]->Sumw2();
    hisname.Form("hPYTHIAPtLpi%dTC",i);
    fMassHistLpiTC[index]=new TH1F(hisname.Data(),hisname.Data(),100,fLowmasslimit,fUpmasslimit);
    fMassHistLpiTC[index]->Sumw2();

    hisname.Form("hPYTHIAPtKp%d",i);
    fMassHistKp[index]=new TH1F(hisname.Data(),hisname.Data(),100,fLowmasslimit,fUpmasslimit);
    fMassHistKp[index]->Sumw2();
    hisname.Form("hPYTHIAPtKp%dTC",i);
    fMassHistKpTC[index]=new TH1F(hisname.Data(),hisname.Data(),100,fLowmasslimit,fUpmasslimit);
    fMassHistKpTC[index]->Sumw2();

  }
  
  for(Int_t i=0; i<4*fNPtBins; i++){
    fOutput->Add(fMassHist[i]);
    fOutput->Add(fMassHistTC[i]);
    fOutput->Add(fMassHistLpi[i]);
    fOutput->Add(fMassHistLpiTC[i]);
    fOutput->Add(fMassHistKp[i]);
    fOutput->Add(fMassHistKpTC[i]);
    fOutput->Add(fMassHistDk[i]);
    fOutput->Add(fMassHistDkTC[i]);
    fOutput->Add(fMassHist3Pr[i]);
    fOutput->Add(fMassHist3PrTC[i]);
  }

  fHistNEvents = new TH1F("fHistNEvents", "Number of processed events; ; Events",3,-1.5,1.5);
  fHistNEvents->Sumw2();
  fHistNEvents->SetMinimum(0);
  fOutput->Add(fHistNEvents);

  fhChi2 = new TH1F("fhChi2", "Chi2",100,0.,10.);
  fhChi2->Sumw2();
  fOutput->Add(fhChi2);

  
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
  hisname.Form("hnopIDp");
  TH1F *hnoProtonIDpTot=new TH1F(hisname.Data(),hisname.Data(),100,0.,10.0);
  fOutputMC->Add(hnoProtonIDpTot);
  hisname.Form("hbnopIDp");
  TH1F *hBnoProtonIDpTot=new TH1F(hisname.Data(),hisname.Data(),100,0.,10.0);
  fOutputMC->Add(hBnoProtonIDpTot);

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
  hisname.Form("hnokIDk");
  TH1F *hnoKaonIDkTot=new TH1F(hisname.Data(),hisname.Data(),100,0.,10.0);
  fOutputMC->Add(hnoKaonIDkTot);
  hisname.Form("hbnokIDk");
  TH1F *hBnoKaonIDkTot=new TH1F(hisname.Data(),hisname.Data(),100,0.,10.0);
  fOutputMC->Add(hBnoKaonIDkTot);


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
  hisname.Form("hnopiIDpi");
  TH1F *hnoPionIDpiTot=new TH1F(hisname.Data(),hisname.Data(),100,0.,10.0);
  fOutputMC->Add(hnoPionIDpiTot);
  hisname.Form("hbnopiIDpi");
  TH1F *hBnoPionIDpiTot=new TH1F(hisname.Data(),hisname.Data(),100,0.,10.0);
  fOutputMC->Add(hBnoPionIDpiTot);

  hisname.Form("hpid0Prong");
  TH1F *hPiond0Prong=new TH1F(hisname.Data(),hisname.Data(),100,-0.1,-0.1);
  fOutputMC->Add(hPiond0Prong);
  hisname.Form("hbpid0Prong");
  TH1F *hBPiond0Prong=new TH1F(hisname.Data(),hisname.Data(),100,-0.1,0.1);
  fOutputMC->Add(hBPiond0Prong);

  hisname.Form("hpiSignalVspTOF");
  TH2F *hpiSignalVspTOF=new TH2F(hisname.Data(),hisname.Data(),100,0.,5.0,100,12000.,50000.0);
  fOutputMC->Add(hpiSignalVspTOF);
  hisname.Form("hpiSignalVspTPC");
  TH2F *hpiSignalVspTPC=new TH2F(hisname.Data(),hisname.Data(),100,0.,5.0,150,0.,150.0);
  fOutputMC->Add(hpiSignalVspTPC);
  hisname.Form("hbpiSignalVspTOF");
  TH2F *hbpiSignalVspTOF=new TH2F(hisname.Data(),hisname.Data(),100,0.,5.0,100,12000.,50000.0);
  fOutputMC->Add(hbpiSignalVspTOF);
  hisname.Form("hbpiSignalVspTPC");
  TH2F *hbpiSignalVspTPC=new TH2F(hisname.Data(),hisname.Data(),100,0.,5.0,150,0.,150.0);
  fOutputMC->Add(hbpiSignalVspTPC);
  hisname.Form("hpiSigmaVspTOF");
  TH2F *hpiSigmaVspTOF=new TH2F(hisname.Data(),hisname.Data(),100,0.,5.0,100,-10.0,10.0);
  fOutputMC->Add(hpiSigmaVspTOF);
  hisname.Form("hpiSigmaVspTPC");
  TH2F *hpiSigmaVspTPC=new TH2F(hisname.Data(),hisname.Data(),100,0.,5.0,100,-10.0,10.0);
  fOutputMC->Add(hpiSigmaVspTPC);
  hisname.Form("hbpiSigmaVspTOF");
  TH2F *hBpiSigmaVspTOF=new TH2F(hisname.Data(),hisname.Data(),100,0.,5.0,100,-10.0,10.0);
  fOutputMC->Add(hBpiSigmaVspTOF);
  hisname.Form("hbpiSigmaVspTPC");
  TH2F *hBpiSigmaVspTPC=new TH2F(hisname.Data(),hisname.Data(),100,0.,5.0,100,-10.0,10.0);
  fOutputMC->Add(hBpiSigmaVspTPC);

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

  
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();

  if(fRDCutsProduction->GetIsUsePID()){
    fRDCutsProduction->GetPidHF()->SetPidResponse(fPIDResponse);
    fRDCutsProduction->GetPidpion()->SetPidResponse(fPIDResponse);
    fRDCutsProduction->GetPidprot()->SetPidResponse(fPIDResponse);
    //fUtilPid=new AliAODpidUtil(fPIDResponse);
    fRDCutsProduction->GetPidHF()->SetOldPid(kFALSE);
    fRDCutsProduction->GetPidpion()->SetOldPid(kFALSE);
    fRDCutsProduction->GetPidprot()->SetOldPid(kFALSE);
  }
  if(fRDCutsAnalysis->GetIsUsePID()){
    fRDCutsAnalysis->GetPidHF()->SetPidResponse(fPIDResponse);
    fRDCutsAnalysis->GetPidpion()->SetPidResponse(fPIDResponse);
    fRDCutsAnalysis->GetPidprot()->SetPidResponse(fPIDResponse);
    fRDCutsAnalysis->GetPidHF()->SetOldPid(kFALSE);
    fRDCutsAnalysis->GetPidpion()->SetOldPid(kFALSE);
    fRDCutsAnalysis->GetPidprot()->SetOldPid(kFALSE);
  }

  PostData(1,fOutput);
  if (fFillVarHists) PostData(3,fOutputMC);
  PostData(4,fNentries); 
  fCounter = new AliNormalizationCounter("NormalizationCounter");
   fCounter->Init();
    PostData(5,fCounter);
  if (fFillNtuple) {
    //OpenFile(3); // 2 is the slot number of the ntuple
    fNtupleLambdac = new TNtuple("fNtupleLambdac","D +",
				 "pdg:Px:Py:Pz:PtTrue:VxTrue:VyTrue:VzTrue:Ptpi:PtK:Ptpi2:PtRec:PointingAngle:DecLeng:VxRec:VyRec:VzRec:InvMass:sigvert:d0Pi:d0K:d0Pi2:dca:d0square");  
    PostData(6,fNtupleLambdac);
  }

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSELambdacUp::UserExec(Option_t */*option*/)
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

  if(!aod) return;

  TClonesArray *arrayMC=0;
  AliAODMCHeader *mcHeader=0;

  // load MC particles
  if(fReadMC){
    
    arrayMC =  (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
    if(!arrayMC) {
      printf("AliAnalysisTaskSELambdacUp::UserExec: MC particles branch not found!\n");
      return;
    }
 

    // load MC header
    mcHeader =  (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if(!mcHeader) {
      printf("AliAnalysisTaskSELambdacUp::UserExec: MC header branch not found!\n");
      return;
    }
  }

  TString fillthis="";
 


  if(!array3Prong || !aod) {
    printf("AliAnalysisTaskSELambdacUp::UserExec: Charm3Prong branch not found!\n");
    return;
  }
  if(!arrayLikeSign) {
    printf("AliAnalysisTaskSELambdacUp::UserExec: LikeSign3Prong branch not found!\n");
    //  return;
  }

  fNentries->Fill(0);
  fCounter->StoreEvent(aod,fRDCutsProduction,fReadMC);
  TString trigclass=aod->GetFiredTriggerClasses();
  if(trigclass.Contains("C0SMH-B-NOPF-ALLNOTRD") || trigclass.Contains("C0SMH-B-NOPF-ALL")) fNentries->Fill(14);
  Bool_t isEvSelAnCuts,isEvSelProdCuts;
  isEvSelProdCuts=fRDCutsProduction->IsEventSelected(aod);
  isEvSelAnCuts=fRDCutsAnalysis->IsEventSelected(aod);
  if(!isEvSelProdCuts){
    if(fRDCutsProduction->GetWhyRejection()==1) // rejected for pileup
      fNentries->Fill(13);
    return;
  }

  // AOD primary vertex
  AliAODVertex *vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();
  if(!vtx1) return;

  Int_t n3Prong = array3Prong->GetEntriesFast();
  Int_t nSelectedloose[1]={0};
  Int_t nSelectedtight[1]={0};

  for (Int_t i3Prong = 0; i3Prong < n3Prong; i3Prong++) {
    AliAODRecoDecayHF3Prong *d = (AliAODRecoDecayHF3Prong*)array3Prong->UncheckedAt(i3Prong);

    Bool_t unsetvtx=kFALSE;
    if(!d->GetOwnPrimaryVtx()){
      d->SetOwnPrimaryVtx(vtx1);
      unsetvtx=kTRUE;
    }

    //if(d->GetSelectionMap()) {if(!d->HasSelectionBit(AliRDHFCuts::kLcCuts)) continue;}
    //if(d->GetSelectionMap()) {if(!d->HasSelectionBit(AliRDHFCuts::kLcPID)) continue;}

    Int_t isSelectedTracks = fRDCutsProduction->IsSelected(d,AliRDHFCuts::kTracks,aod);
    if(!isSelectedTracks) continue;

    if (fRDCutsProduction->IsInFiducialAcceptance(d->Pt(),d->Y(4122))) 
      {fNentries->Fill(6);}else{continue;}

    Int_t ptbin=fRDCutsProduction->PtBin(d->Pt());
    if(ptbin==-1) {fNentries->Fill(4); continue;} //out of bounds

    FillMassHists(aod,d,arrayMC,fRDCutsProduction,nSelectedloose,nSelectedtight);
    if(fFillVarHists) FillVarHists(d,arrayMC,fRDCutsProduction,/*fOutputMC,*/aod);

    

    /*
    //start OS analysis
    if(labDp<0)fHistOSbkg->Fill(d->InvMassDplus());
    fHistOS->Fill(d->InvMassDplus());
    */

    if(unsetvtx) d->UnsetOwnPrimaryVtx();
  }
  fCounter->StoreCandidates(aod,nSelectedloose[0],kTRUE);
  fCounter->StoreCandidates(aod,nSelectedtight[0],kFALSE);

  
  PostData(1,fOutput); 
  if (fFillVarHists) PostData(3,fOutputMC);
  

  PostData(4,fNentries);
  PostData(5,fCounter);
      
  return;
}



//________________________________________________________________________
void AliAnalysisTaskSELambdacUp::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
  if (fDebug > 1) printf("AnalysisTaskSELambdac: Terminate() \n");

  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {     
    printf("ERROR: fOutput not available\n");
    return;
  }
  //fHistNEvents = dynamic_cast<TH1F*>(fOutput->FindObject("fHistNEvents"));

  if(fFillNtuple){
    fNtupleLambdac = dynamic_cast<TNtuple*>(GetOutputData(3));
  }

 
  return;
}

//________________________________________________________________________
Int_t AliAnalysisTaskSELambdacUp::MatchToMCLambdac(AliAODRecoDecayHF3Prong *d,TClonesArray *arrayMC) const{
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
	Int_t granMotherLabel=motherPart->GetMother();
	if(granMotherLabel<0) return 0;
	AliAODMCParticle *granMotherPart = (AliAODMCParticle*)arrayMC->At(granMotherLabel);
	if(!granMotherPart) continue;
	Int_t granMotherPdg  = TMath::Abs(granMotherPart->GetPdgCode());
	if(granMotherPdg ==4122) {
	  if(GetLambdacDaugh(granMotherPart,arrayMC)) {lambdacLab[i]=granMotherLabel;continue;}
	}
      }
    }
  }

  if(lambdacLab[0]==lambdacLab[1] && lambdacLab[1]==lambdacLab[2]) {return lambdacLab[0];}
  return 0;

}
//------------------------
Bool_t AliAnalysisTaskSELambdacUp::GetLambdacDaugh(AliAODMCParticle *part,TClonesArray *arrayMC) const{
  // check if the particle is a lambdac and if its decay mode is the correct one 
  Int_t numberOfLambdac=0;
  if(TMath::Abs(part->GetPdgCode())!=4122) return kFALSE;
  Int_t daughTmp[2];
  daughTmp[0]=part->GetDaughterLabel(0);
  daughTmp[1]=part->GetDaughterLabel(1);
  Int_t nDaugh = (Int_t)part->GetNDaughters();
  if(nDaugh<2) return kFALSE;
  if(nDaugh>3) return kFALSE;
  AliAODMCParticle* pdaugh1 = (AliAODMCParticle*)arrayMC->At(part->GetDaughterLabel(0));
  if(!pdaugh1) {return kFALSE;}
  Int_t number1 = TMath::Abs(pdaugh1->GetPdgCode());
  AliAODMCParticle* pdaugh2 = (AliAODMCParticle*)arrayMC->At(part->GetDaughterLabel(1));
  if(!pdaugh2) {return kFALSE;}
  Int_t number2 = TMath::Abs(pdaugh2->GetPdgCode());

  if(nDaugh==3){
    Int_t thirdDaugh=part->GetDaughterLabel(1)-1;
    AliAODMCParticle* pdaugh3 = (AliAODMCParticle*)arrayMC->At(thirdDaugh);
    Int_t number3 = TMath::Abs(pdaugh3->GetPdgCode());
    if((number1==321 && number2==211 && number3==2212) ||
       (number1==211 && number2==321 && number3==2212) ||
       (number1==211 && number2==2212 && number3==321) ||
       (number1==321 && number2==2212 && number3==211) ||
       (number1==2212 && number2==321 && number3==211) ||
       (number1==2212 && number2==211 && number3==321)) {
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
      AliAODMCParticle* pdaughK1 = (AliAODMCParticle*)arrayMC->At(pdaugh2->GetDaughterLabel(0));
      AliAODMCParticle* pdaughK2 = (AliAODMCParticle*)arrayMC->At(pdaugh2->GetDaughterLabel(1));
      if(!pdaughK1) return kFALSE;
      if(!pdaughK2) return kFALSE;
      if((TMath::Abs(pdaughK1->GetPdgCode())==211 && TMath::Abs(pdaughK2->GetPdgCode())==321) || (TMath::Abs(pdaughK1->GetPdgCode())==321 && TMath::Abs(pdaughK2->GetPdgCode())==211)) numberOfLambdac++;
    }

    if((number1==313 && number2==2212)){
      nfiglieK=pdaugh1->GetNDaughters();
      if(nfiglieK!=2) return kFALSE;
      AliAODMCParticle* pdaughK1 = (AliAODMCParticle*)arrayMC->At(pdaugh1->GetDaughterLabel(0));
      AliAODMCParticle* pdaughK2 = (AliAODMCParticle*)arrayMC->At(pdaugh1->GetDaughterLabel(1));
      if(!pdaughK1) return kFALSE;
      if(!pdaughK2) return kFALSE;
      if((TMath::Abs(pdaughK1->GetPdgCode())==211 && TMath::Abs(pdaughK2->GetPdgCode())==321) || (TMath::Abs(pdaughK1->GetPdgCode())==321 && TMath::Abs(pdaughK2->GetPdgCode())==211)) numberOfLambdac++;
    }

    //Lambda -> Delta++ k
    Int_t nfiglieDelta=0;
    if(number1==321 && number2==2224){
      nfiglieDelta=pdaugh2->GetNDaughters();
      if(nfiglieDelta!=2) return kFALSE;
      AliAODMCParticle *pdaughD1=(AliAODMCParticle*)arrayMC->At(pdaugh2->GetDaughterLabel(0));
      AliAODMCParticle *pdaughD2=(AliAODMCParticle*)arrayMC->At(pdaugh2->GetDaughterLabel(1));
      if(!pdaughD1) return kFALSE;
      if(!pdaughD2) return kFALSE;
      if((TMath::Abs(pdaughD1->GetPdgCode())==211 && TMath::Abs(pdaughD2->GetPdgCode())==2212) || (TMath::Abs(pdaughD1->GetPdgCode())==2212 && TMath::Abs(pdaughD2->GetPdgCode())==211)) numberOfLambdac++;
    }
    if(number1==2224 && number2==321){
      nfiglieDelta=pdaugh1->GetNDaughters();
      if(nfiglieDelta!=2) return kFALSE;
      AliAODMCParticle* pdaughD1 = (AliAODMCParticle*)arrayMC->At(pdaugh1->GetDaughterLabel(0));
      AliAODMCParticle* pdaughD2 = (AliAODMCParticle*)arrayMC->At(pdaugh1->GetDaughterLabel(1));
      if(!pdaughD1) return kFALSE;
      if(!pdaughD2) return kFALSE;
      if((TMath::Abs(pdaughD1->GetPdgCode())==211 && TMath::Abs(pdaughD2->GetPdgCode())==2212) || (TMath::Abs(pdaughD1->GetPdgCode())==2212 && TMath::Abs(pdaughD2->GetPdgCode())==211)) numberOfLambdac++;
    }
    

    //Lambdac -> Lambda(1520) pi
    Int_t nfiglieLa=0;
    if(number1==3124 && number2==211){
      nfiglieLa=pdaugh1->GetNDaughters();
      if(nfiglieLa!=2) return kFALSE;
      AliAODMCParticle *pdaughL1=(AliAODMCParticle*)arrayMC->At(pdaugh1->GetDaughterLabel(0));
      AliAODMCParticle *pdaughL2=(AliAODMCParticle*)arrayMC->At(pdaugh1->GetDaughterLabel(1));
      if(!pdaughL1) return kFALSE;
      if(!pdaughL2) return kFALSE;
      if((TMath::Abs(pdaughL1->GetPdgCode())==321 && TMath::Abs(pdaughL2->GetPdgCode())==2212) || (TMath::Abs(pdaughL1->GetPdgCode())==2212 && TMath::Abs(pdaughL2->GetPdgCode())==321)) numberOfLambdac++;
    }
    if(number1==211 && number2==3124){
      nfiglieLa=pdaugh2->GetNDaughters();
      if(nfiglieLa!=2) return kFALSE;
      AliAODMCParticle *pdaughL1=(AliAODMCParticle*)arrayMC->At(pdaugh2->GetDaughterLabel(0));
      AliAODMCParticle *pdaughL2=(AliAODMCParticle*)arrayMC->At(pdaugh2->GetDaughterLabel(1));
      if(!pdaughL1) return kFALSE;
      if(!pdaughL2) return kFALSE;
      if((TMath::Abs(pdaughL1->GetPdgCode())==321 && TMath::Abs(pdaughL2->GetPdgCode())==2212) || (TMath::Abs(pdaughL1->GetPdgCode())==2212 && TMath::Abs(pdaughL2->GetPdgCode())==321)) numberOfLambdac++;
   
    }
  }

  if(numberOfLambdac>0) {return kTRUE;}
  return kFALSE;
}
//-----------------------------
Bool_t AliAnalysisTaskSELambdacUp::IspKpiMC(AliAODRecoDecayHF3Prong *d,TClonesArray *arrayMC) const{
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
Bool_t AliAnalysisTaskSELambdacUp::IspiKpMC(AliAODRecoDecayHF3Prong *d,TClonesArray *arrayMC) const{

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
Bool_t AliAnalysisTaskSELambdacUp::VertexingKF(AliAODRecoDecayHF3Prong *d,Int_t *pdgs,Double_t field) const{
  // apply vertexing KF 
  Int_t iprongs[3]={0,1,2};
  Double_t mass[2]={0.,0.};
  Bool_t constraint=kFALSE;
  if(fCutsKF[0]>0.)constraint=kTRUE;
  //topological constr
  AliKFParticle *lambdac=d->ApplyVertexingKF(iprongs,3,pdgs,constraint,field,mass);
  if(!lambdac) return kFALSE;
  if(lambdac->GetChi2()/lambdac->GetNDF()>fCutsKF[1]) return kFALSE;
  return kTRUE;
}
//-------------------------------------
void AliAnalysisTaskSELambdacUp::IspiKpResonant(AliAODRecoDecayHF3Prong *d,Double_t field,Int_t *resNumber) const{
  
  // apply PID using the resonant channels 
  //if lambda* -> pk
  Double_t mass[2]={1.520,0.005};
  Int_t ipRes[2]={1,2};
  Int_t pdgres[2]={321,2212};
  AliKFParticle *lambda1520=d->ApplyVertexingKF(ipRes,2,pdgres,kFALSE,field,mass);
  Double_t probLa=TMath::Prob(lambda1520->GetChi2(),lambda1520->GetNDF());
  if(probLa>0.1) resNumber[0]=1;
  //K* -> kpi
  mass[0]=0.8961;mass[1]=0.03;
  ipRes[0]=0;ipRes[1]=1;
  pdgres[0]=211;pdgres[1]=321;
  AliKFParticle *kstar=d->ApplyVertexingKF(ipRes,2,pdgres,kFALSE,field,mass);
  Double_t probKa=TMath::Prob(kstar->GetChi2(),kstar->GetNDF());
  if(probKa>0.1) resNumber[1]=1;

 // Delta++
  mass[0]=1.232;mass[1]=0.15;
  ipRes[0]=0;ipRes[1]=2;
  pdgres[0]=211;pdgres[1]=2122;
  AliKFParticle *delta=d->ApplyVertexingKF(ipRes,2,pdgres,kFALSE,field,mass);
  Double_t probDe=TMath::Prob(delta->GetChi2(),delta->GetNDF());
  if(probDe>0.1) resNumber[2]=1;

  return ;

}
//-------------------------------------
void AliAnalysisTaskSELambdacUp::IspKpiResonant(AliAODRecoDecayHF3Prong *d,Double_t field,Int_t *resNumber) const{
   
  // apply PID using the resonant channels 
  //if lambda* -> pk
  Double_t mass[2]={1.520,0.005};
  Int_t ipRes[2]={0,1};
  Int_t pdgres[2]={2212,321};
  AliKFParticle *lambda1520=d->ApplyVertexingKF(ipRes,2,pdgres,kFALSE,field,mass);
  Double_t probLa=TMath::Prob(lambda1520->GetChi2(),lambda1520->GetNDF());
  if(probLa>0.1) resNumber[0]=1;
  //K* -> kpi
  mass[0]=0.8961;mass[1]=0.03;
  ipRes[0]=1;ipRes[1]=2;
  pdgres[1]=211;pdgres[0]=321;
  AliKFParticle *kstar=d->ApplyVertexingKF(ipRes,2,pdgres,kFALSE,field,mass);
  Double_t probKa=TMath::Prob(kstar->GetChi2(),kstar->GetNDF());
  if(probKa>0.1) resNumber[1]=1;

  //	Delta++
  mass[0]=1.232;mass[1]=0.15;
  ipRes[0]=0;ipRes[1]=2;
  pdgres[0]=2122;pdgres[1]=211;
  AliKFParticle *delta=d->ApplyVertexingKF(ipRes,2,pdgres,kFALSE,field,mass);
  Double_t probDe=TMath::Prob(delta->GetChi2(),delta->GetNDF());
  if(probDe>0.1) resNumber[2]=1;

  return ;

}
//---------------------------
void AliAnalysisTaskSELambdacUp::FillMassHists(AliAODEvent *aod,AliAODRecoDecayHF3Prong *part,
						TClonesArray *arrayMC, AliRDHFCutsLctopKpi *cuts,Int_t *nSelectedloose,Int_t *nSelectedtight)
{

  //if MC PID or no PID, unset PID
  if(!fRealPid) cuts->SetUsePID(kFALSE);
  Int_t selection=cuts->IsSelected(part,AliRDHFCuts::kCandidate,aod);
  if(selection>0){
    nSelectedloose[0]=nSelectedloose[0]+1;
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
      //labDp = MatchToMCLambdac(part,arrayMC);
      Int_t pdgDgLctopKpi[3]={321,2212,211};
      labDp = part->MatchToMC(4122,arrayMC,3,pdgDgLctopKpi);
      if(labDp>=0){
	AliAODMCParticle *partDp = (AliAODMCParticle*)arrayMC->At(labDp);
	AliAODMCParticle *dg0 = (AliAODMCParticle*)arrayMC->At(partDp->GetDaughterLabel(0));
	AliAODMCParticle *dg1 = (AliAODMCParticle*)arrayMC->At(partDp->GetDaughterLabel(1));
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
    Double_t invMasspiKpDk=-1.;
    Double_t invMasspKpiDk=-1.;
    Double_t invMasspiKp3Pr=-1.;
    Double_t invMasspKpi3Pr=-1.;
    Double_t field=aod->GetMagneticField();
    //apply MC PID
    if(fReadMC && fMCPid){
      if(IspKpiMC(part,arrayMC)) invMasspKpi=part->InvMassLcpKpi();
      if(IspiKpMC(part,arrayMC)) invMasspiKp=part->InvMassLcpiKp();
    }

    // apply realistic PID
    if(fRealPid){
      if(selection==1 || selection==3) invMasspKpi=part->InvMassLcpKpi();
      if(selection>=2) invMasspiKp=part->InvMassLcpiKp();
    }
    //apply PID using resonances
    if (fResPid && fRealPid) {
      Int_t ispKpi[3]={0,0,0};
      IspKpiResonant(part,field,ispKpi);
      Int_t ispiKp[3]={0,0,0};
      IspiKpResonant(part,field,ispiKp);
      if (selection==3 || selection==1) {
       if(ispKpi[0]==1){
        invMasspKpiLpi=part->InvMassLcpKpi();
       }
       if(ispKpi[1]==1){
        invMasspKpiKp=part->InvMassLcpKpi();
       }
       if(ispKpi[2]==1){
        invMasspKpiDk=part->InvMassLcpKpi();
       }
       if(ispKpi[2]==0 && ispKpi[1]==0 && ispKpi[0]==0){
        invMasspKpi3Pr=part->InvMassLcpKpi();
       }
      }
      if(selection>=2) {
       if(ispiKp[0]==1){
        invMasspiKpLpi=part->InvMassLcpiKp();
       }
       if(ispiKp[1]==1){
        invMasspiKpKp=part->InvMassLcpiKp();
       }
       if(ispiKp[2]==1){
        invMasspiKpDk=part->InvMassLcpiKp();
       }
       if(ispiKp[2]==0 && ispiKp[1]==0 && ispiKp[0]==0){
        invMasspiKp3Pr=part->InvMassLcpiKp();
       }
      }

     }
  

    if(invMasspiKp<0. && invMasspKpi<0.) return;
   
    Int_t passTightCuts=0;
    if(fAnalysis) {
     passTightCuts=fRDCutsAnalysis->IsSelected(part,AliRDHFCuts::kCandidate,aod);
     if(fUseKF){
      Int_t pdgs[3]={0,0,0};
      if(invMasspKpi>0.){
       pdgs[0]=2212;pdgs[1]=321;pdgs[2]=211;
       if(!VertexingKF(part,pdgs,field)) {
        invMasspKpi=-1.;
	invMasspKpi3Pr=-1.;
	invMasspKpiDk=-1.;
	invMasspKpiLpi=-1.;
	invMasspKpiKp=-1.;
       }
      }
      if(invMasspiKp>0.){
	pdgs[0]=211;pdgs[1]=321;pdgs[2]=2212;
	if(!VertexingKF(part,pdgs,field)) {
	  invMasspiKp=-1.;
	  invMasspiKp3Pr=-1.;
	  invMasspiKpDk=-1.;
	  invMasspiKpLpi=-1.;
	  invMasspiKpKp=-1.;
       }
     }
   }
    if(passTightCuts>0) nSelectedtight[0]=nSelectedtight[0]+1;
  }

    Float_t tmp[24];
    if (fFillNtuple) {
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
      PostData(6,fNtupleLambdac);
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
      if(invMasspiKpDk>0. && invMasspKpiDk>0.){
       if(invMasspiKpDk>0.) fMassHistDk[index]->Fill(invMasspiKpDk,0.5);
       if(invMasspKpiDk>0.) fMassHistDk[index]->Fill(invMasspKpiDk,0.5);
      }else{
       if(invMasspiKpDk>0.) fMassHistDk[index]->Fill(invMasspiKpDk);
       if(invMasspKpiDk>0.) fMassHistDk[index]->Fill(invMasspKpiDk);
      }
      if(invMasspiKp3Pr>0. && invMasspKpi3Pr>0.){
       if(invMasspiKp3Pr>0.) fMassHist3Pr[index]->Fill(invMasspiKp3Pr,0.5);
       if(invMasspKpi3Pr>0.) fMassHist3Pr[index]->Fill(invMasspKpi3Pr,0.5);
      }else{
       if(invMasspiKp3Pr>0.) fMassHist3Pr[index]->Fill(invMasspiKp3Pr);
       if(invMasspKpi3Pr>0.) fMassHist3Pr[index]->Fill(invMasspKpi3Pr);
     }

      if(passTightCuts>0){
	if(passTightCuts==3){
	  if(invMasspiKp>0.) fMassHistTC[index]->Fill(invMasspiKp);
	  if(invMasspKpi>0.) fMassHistTC[index]->Fill(invMasspKpi);
	}else{
	  if(passTightCuts==2) fMassHistTC[index]->Fill(invMasspiKp);
	  if(passTightCuts==1) fMassHistTC[index]->Fill(invMasspKpi);
	}
	if(passTightCuts==3){
	  if(invMasspiKpLpi>0.) fMassHistLpiTC[index]->Fill(invMasspiKpLpi);
	  if(invMasspKpiLpi>0.) fMassHistLpiTC[index]->Fill(invMasspKpiLpi);
	}else{
	  if(passTightCuts==2) fMassHistLpiTC[index]->Fill(invMasspiKpLpi);
	  if(passTightCuts==1) fMassHistLpiTC[index]->Fill(invMasspKpiLpi);
	}
	if(passTightCuts==3){
	  if(invMasspiKpKp>0.) fMassHistKpTC[index]->Fill(invMasspiKpKp);
	  if(invMasspKpiKp>0.) fMassHistKpTC[index]->Fill(invMasspKpiKp);
	}else{
	  if(passTightCuts==2) fMassHistKpTC[index]->Fill(invMasspiKpKp);
	  if(passTightCuts==1) fMassHistKpTC[index]->Fill(invMasspKpiKp);
	}
       if(passTightCuts==3){
        if(invMasspiKpDk>0.) fMassHistDkTC[index]->Fill(invMasspiKpDk);
        if(invMasspKpiDk>0.) fMassHistDkTC[index]->Fill(invMasspKpiDk);
       }else{
        if(passTightCuts==2) fMassHistDkTC[index]->Fill(invMasspiKpDk);
        if(passTightCuts==1) fMassHistDkTC[index]->Fill(invMasspKpiDk);
       }
       if(passTightCuts==3){
        if(invMasspiKp3Pr>0.) fMassHist3PrTC[index]->Fill(invMasspiKp3Pr,0.5);
        if(invMasspKpi3Pr>0.) fMassHist3PrTC[index]->Fill(invMasspKpi3Pr,0.5);
       }else{
        if(passTightCuts==2) fMassHist3PrTC[index]->Fill(invMasspiKp3Pr);
        if(passTightCuts==1) fMassHist3PrTC[index]->Fill(invMasspKpi3Pr);
      }

      }

      if(fReadMC){
	if(labDp>=0) {
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
	  if(invMasspiKpDk>0. && invMasspKpiDk>0.){
	    if(invMasspiKpDk>0.) fMassHistDk[index]->Fill(invMasspiKpDk,0.5);
	    if(invMasspKpiDk>0.) fMassHistDk[index]->Fill(invMasspKpiDk,0.5);
	  }else{
	    if(invMasspiKpDk>0.) fMassHistDk[index]->Fill(invMasspiKpDk);
	    if(invMasspKpiDk>0.) fMassHistDk[index]->Fill(invMasspKpiDk);
	  }
	  if(invMasspiKp3Pr>0. && invMasspKpi3Pr>0.){
	    if(invMasspiKp3Pr>0.) fMassHist3Pr[index]->Fill(invMasspiKp3Pr,0.5);
	    if(invMasspKpi3Pr>0.) fMassHist3Pr[index]->Fill(invMasspKpi3Pr,0.5);
	  }else{
	    if(invMasspiKp3Pr>0.) fMassHist3Pr[index]->Fill(invMasspiKp3Pr);
	    if(invMasspKpi3Pr>0.) fMassHistDk[index]->Fill(invMasspKpi3Pr);
	  }

	  if(passTightCuts>0){
	    if(passTightCuts==3){
	      if(invMasspiKp>0.) fMassHistTC[index]->Fill(invMasspiKp);
	      if(invMasspKpi>0.) fMassHistTC[index]->Fill(invMasspKpi);
	    }else{
	      if(passTightCuts==2) fMassHistTC[index]->Fill(invMasspiKp);
	      if(passTightCuts==1) fMassHistTC[index]->Fill(invMasspKpi);
	    }
	    if(passTightCuts==3){
	      if(invMasspiKpLpi>0.) fMassHistLpiTC[index]->Fill(invMasspiKpLpi);
	      if(invMasspKpiLpi>0.) fMassHistLpiTC[index]->Fill(invMasspKpiLpi);
	    }else{
	      if(passTightCuts==2) fMassHistLpiTC[index]->Fill(invMasspiKpLpi);
	      if(passTightCuts==1) fMassHistLpiTC[index]->Fill(invMasspKpiLpi);
	    } 
	    if(passTightCuts==3){
	      if(invMasspiKpKp>0.) fMassHistKpTC[index]->Fill(invMasspiKpKp);
	      if(invMasspKpiKp>0.) fMassHistKpTC[index]->Fill(invMasspKpiKp);
	    }else{
	      if(passTightCuts==2) fMassHistKpTC[index]->Fill(invMasspiKpKp);
	      if(passTightCuts==1) fMassHistKpTC[index]->Fill(invMasspKpiKp);
	    }
	    if(passTightCuts==3){
	      if(invMasspiKpDk>0.) fMassHistDkTC[index]->Fill(invMasspiKpDk);
	      if(invMasspKpiDk>0.) fMassHistDkTC[index]->Fill(invMasspKpiDk);
	    }else{
	      if(passTightCuts==2) fMassHistDkTC[index]->Fill(invMasspiKpDk);
	      if(passTightCuts==1) fMassHistDkTC[index]->Fill(invMasspKpiDk);
	    }
	    if(passTightCuts==3){
	      if(invMasspiKp3Pr>0.) fMassHist3PrTC[index]->Fill(invMasspiKp3Pr);
	      if(invMasspKpi3Pr>0.) fMassHist3PrTC[index]->Fill(invMasspKpi3Pr);
	    }else{
	      if(passTightCuts==2) fMassHist3PrTC[index]->Fill(invMasspiKp3Pr);
	      if(passTightCuts==1) fMassHist3PrTC[index]->Fill(invMasspKpi3Pr);
	    }
	  }      
     
	}else{
          Bool_t isPYTHIA=CheckOrigin(part,arrayMC);
          if(!isPYTHIA){
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
	  if(invMasspiKpDk>0. && invMasspKpiDk>0.){
	    if(invMasspiKpDk>0.) fMassHistDk[index]->Fill(invMasspiKpDk,0.5);
	    if(invMasspKpiDk>0.) fMassHistDk[index]->Fill(invMasspKpiDk,0.5);
	  }else{
	    if(invMasspiKpDk>0.) fMassHistDk[index]->Fill(invMasspiKpDk);
	    if(invMasspKpiDk>0.) fMassHistDk[index]->Fill(invMasspKpiDk);
	  }
	  if(invMasspiKp3Pr>0. && invMasspKpi3Pr>0.){
	    if(invMasspiKp3Pr>0.) fMassHist3Pr[index]->Fill(invMasspiKp3Pr,0.5);
	    if(invMasspKpi3Pr>0.) fMassHist3Pr[index]->Fill(invMasspKpi3Pr,0.5);
	  }else{
	    if(invMasspiKp3Pr>0.) fMassHist3Pr[index]->Fill(invMasspiKp3Pr);
	    if(invMasspKpi3Pr>0.) fMassHistDk[index]->Fill(invMasspKpi3Pr);
	  }
	  if(passTightCuts==3){
	    fMassHistTC[index]->Fill(invMasspiKp);
	    fMassHistTC[index]->Fill(invMasspKpi);
	  }else{
	    if(passTightCuts==2) fMassHistTC[index]->Fill(invMasspiKp);
	    if(passTightCuts==1) fMassHistTC[index]->Fill(invMasspKpi);
	  }
	  if(passTightCuts==3){
	    if(invMasspiKpLpi>0.) fMassHistLpiTC[index]->Fill(invMasspiKpLpi);
	    if(invMasspKpiLpi>0.) fMassHistLpiTC[index]->Fill(invMasspKpiLpi);
	  }else{
	    if(passTightCuts==2) fMassHistLpiTC[index]->Fill(invMasspiKpLpi);
	    if(passTightCuts==1) fMassHistLpiTC[index]->Fill(invMasspKpiLpi);
	  } 
	  if(passTightCuts==3){
	    if(invMasspiKpKp>0.) fMassHistKpTC[index]->Fill(invMasspiKpKp);
	    if(invMasspKpiKp>0.) fMassHistKpTC[index]->Fill(invMasspKpiKp);
	  }else{
	    if(passTightCuts==2) fMassHistKpTC[index]->Fill(invMasspiKpKp);
	    if(passTightCuts==1) fMassHistKpTC[index]->Fill(invMasspKpiKp);
	  }
	  if(passTightCuts==3){
	    if(invMasspiKpDk>0.) fMassHistDkTC[index]->Fill(invMasspiKpDk);
	    if(invMasspKpiDk>0.) fMassHistDkTC[index]->Fill(invMasspKpiDk);
	  }else{
	    if(passTightCuts==2) fMassHistDkTC[index]->Fill(invMasspiKpDk);
	    if(passTightCuts==1) fMassHistDkTC[index]->Fill(invMasspKpiDk);
	   }
	   if(passTightCuts==3){
	    if(invMasspiKp3Pr>0.) fMassHist3PrTC[index]->Fill(invMasspiKp3Pr);
	    if(invMasspKpi3Pr>0.) fMassHist3PrTC[index]->Fill(invMasspKpi3Pr);
	   }else{
	    if(passTightCuts==2) fMassHist3PrTC[index]->Fill(invMasspiKp3Pr);
	    if(passTightCuts==1) fMassHist3PrTC[index]->Fill(invMasspKpi3Pr);
	  }
         }else{
          index=GetPYTHIAHistoIndex(iPtBin);
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
   if(passTightCuts==3){
     fMassHistTC[index]->Fill(invMasspiKp);
     fMassHistTC[index]->Fill(invMasspKpi);
    }else{
     if(passTightCuts==2) fMassHistTC[index]->Fill(invMasspiKp);
     if(passTightCuts==1) fMassHistTC[index]->Fill(invMasspKpi);
    }
    if(passTightCuts==3){
     if(invMasspiKpLpi>0.) fMassHistLpiTC[index]->Fill(invMasspiKpLpi);
     if(invMasspKpiLpi>0.) fMassHistLpiTC[index]->Fill(invMasspKpiLpi);
    }else{
     if(passTightCuts==2) fMassHistLpiTC[index]->Fill(invMasspiKpLpi);
     if(passTightCuts==1) fMassHistLpiTC[index]->Fill(invMasspKpiLpi);
    }
    if(passTightCuts==3){
     if(invMasspiKpKp>0.) fMassHistKpTC[index]->Fill(invMasspiKpKp);
     if(invMasspKpiKp>0.) fMassHistKpTC[index]->Fill(invMasspKpiKp);
    }else{
     if(passTightCuts==2) fMassHistKpTC[index]->Fill(invMasspiKpKp);
     if(passTightCuts==1) fMassHistKpTC[index]->Fill(invMasspKpiKp);
    }
   }
	}

      }
    }
  }
  return;
}
//-----------------------
void AliAnalysisTaskSELambdacUp::FillVarHists(AliAODRecoDecayHF3Prong *part,
					       TClonesArray *arrMC,
					       AliRDHFCutsLctopKpi *cuts,
					       /* TList *listout,*/
					       AliAODEvent* aod) {
  //
  // function used in UserExec to fill variable histograms analysing MC
  //

  TString fillthis="";

  Double_t mPDG=TDatabasePDG::Instance()->GetParticle(4122)->Mass();
  Double_t invmasscut=0.05;

  Int_t pdgDgLctopKpi[3]={2212,321,211};
  Int_t lab=-9999;
  if(fReadMC)
    lab=part->MatchToMC(4122,arrMC,3,pdgDgLctopKpi); //return MC particle label if the array corresponds to a Lc, -1 if not (cf. AliAODRecoDecay.cxx)

  Int_t isSelectedPID=cuts->IsSelectedPID(part); // 0 rejected, 1 Lc -> p K- pi+ (K at center because different sign is there),
                                                 // 2 Lc -> pi+ K- p (K at center because different sign is there),
                                                 // 3 both (it should never happen...)

  if (isSelectedPID==0)fNentries->Fill(7);
  if (isSelectedPID==1)fNentries->Fill(8);
  if (isSelectedPID==2)fNentries->Fill(9);
  if (isSelectedPID==3)fNentries->Fill(10);

  AliAODTrack *prong0=(AliAODTrack*)part->GetDaughter(0);
  AliAODTrack *prong1=(AliAODTrack*)part->GetDaughter(1);
  AliAODTrack *prong2=(AliAODTrack*)part->GetDaughter(2);
  if (!prong0 || !prong1 || !prong2) {
    fNentries->Fill(5);
    return;
  }

  Double_t minvLcpKpi = part->InvMassLcpKpi();
  Double_t minvLcpiKp = part->InvMassLcpiKp();

  //check pdg of the prongs
  Int_t labprong[3]={-1,-1,-1};
  if(fReadMC){
    labprong[0]=prong0->GetLabel();
    labprong[1]=prong1->GetLabel();
    labprong[2]=prong2->GetLabel();
  }

  AliAODMCParticle *mcprong=0;
  Int_t pdgProngMC[3]={-1,-1,-1};
  if(fReadMC) {
    for (Int_t iprong=0;iprong<3;iprong++){
      if(labprong[iprong]<0) continue;
      mcprong = (AliAODMCParticle*)arrMC->At(labprong[iprong]);
      pdgProngMC[iprong]=TMath::Abs(mcprong->GetPdgCode());
    }
  }

  Int_t pdgProngPID[3]={-1,-1,-1};
  if(isSelectedPID>0){
    pdgProngPID[1]=321;
    if(isSelectedPID==1) {pdgProngPID[0]=2212;pdgProngPID[2]=211;}
    if(isSelectedPID==2) {pdgProngPID[0]=211;pdgProngPID[2]=2212;}
  }

  //fill hRealTot ---> PID efficiency denominators
  cuts->SetUsePID(kFALSE); //PAOLA
  Int_t selectionCandidate=cuts->IsSelected(part,AliRDHFCuts::kCandidate,aod);//PAOLA
  if(fReadMC && selectionCandidate>0) { // 3prongs candidate x Lc (no PID)
    if(lab>0){ // Signal
      for (Int_t iprong=0; iprong<3; iprong++) {
	switch (pdgProngMC[iprong]) {
	case 2212:
	  fillthis="hpRealTot";
	  ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
	  break;
	case 321:
	  fillthis="hKRealTot";
	  ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
	  break;
	case 211:
	  fillthis="hpiRealTot";
	  ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
	  break;
	default:
	  break;
	}
      }
    } else { // bkg
      for (Int_t iprong=0; iprong<3; iprong++) {
	switch (pdgProngMC[iprong]) {
	case 2212:
	  fillthis="hbpRealTot";
	  ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
	  break;
	case 321:
	  fillthis="hbKRealTot";
	  ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
	  break;
	case 211:
	  fillthis="hbpiRealTot";
	  ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
	  break;
	default:
	  break;
	}
      }
    }
  }


  cuts->SetUsePID(kTRUE); //PAOLA
  Int_t selection=cuts->IsSelected(part,AliRDHFCuts::kCandidate,aod);

  if ( (lab>0 && fReadMC) ||             // signal X MC
       (isSelectedPID>0 && !fReadMC) ) { // signal+bkg X real data

    fillthis="hMass";
    if ( (fReadMC && ((AliAODMCParticle*)arrMC->At(lab))->GetPdgCode() == 4122) ||
	 (!fReadMC && (isSelectedPID>0 && part->Charge()>0)) ) { // Lc
      if ( (pdgProngPID[0]==2212) && (pdgProngPID[1]==321) && (pdgProngPID[2]==211) )
	((TH1F*)fOutputMC->FindObject(fillthis))->Fill(minvLcpKpi); // signal+bkg X MC and real data
      else if ( (pdgProngPID[0]==211) && (pdgProngPID[1]==321) && (pdgProngPID[2]==2212) )
	((TH1F*)fOutputMC->FindObject(fillthis))->Fill(minvLcpiKp); // signal+bkg X MC and real data
    }
    else if ( (fReadMC && ((AliAODMCParticle*)arrMC->At(lab))->GetPdgCode() == -4122) ||
	      (!fReadMC && (isSelectedPID>0 && part->Charge()<0)) ) { // anti-Lc
      if ( (pdgProngPID[0]==2212) && (pdgProngPID[1]==321) && (pdgProngPID[2]==211) )
	((TH1F*)fOutputMC->FindObject(fillthis))->Fill(minvLcpKpi); // signal+bkg X MC and real data
      else if ( (pdgProngPID[0]==211) && (pdgProngPID[1]==321) && (pdgProngPID[2]==2212) )
	((TH1F*)fOutputMC->FindObject(fillthis))->Fill(minvLcpiKp); // signal+bkg X MC and real data
    }

    if (selection>0) { // 3prongs candidate x Lc (yes PID)

      Double_t ptmax=0.;

      for (Int_t iprong=0; iprong<3; iprong++) {
	if (part->PtProng(iprong)>ptmax) ptmax=part->PtProng(iprong);

	AliAODTrack *prong = (AliAODTrack*)part->GetDaughter(iprong);
	AliAODPid *pidObjtrk = (AliAODPid*)prong->GetDetPid();
	AliAODPidHF *pidObj = (AliAODPidHF*)cuts->GetPidHF();
	Bool_t hasTOF=pidObj->CheckStatus(prong,"TOF");
	Bool_t hasTPC=pidObj->CheckStatus(prong,"TPC");
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
	switch (pdgProngPID[iprong]) {
	case 2212:
	  fillthis="hpTOFSignal";
	  ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(tofSignal);
	  fillthis="hpTPCSignal";
	  ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(dedxTPC);
	  fillthis="hpptProng";
	  ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
	  fillthis="hpd0Prong";
	  ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->Getd0Prong(iprong));
	  fillthis="hpSignalVspTPC";
	  ((TH2F*)fOutputMC->FindObject(fillthis))->Fill(momTPC,dedxTPC);
	  fillthis="hpSignalVspTOF";
	  ((TH2F*)fOutputMC->FindObject(fillthis))->Fill(momTOF,tofSignal);
	  AliPID::EParticleType typep;
	  typep=AliPID::EParticleType(4);
	  if(hasTPC) {
	    //Double_t nsigmap = fUtilPid->NumberOfSigmasTPC(prong,typep);
	    Double_t nsigmap = fPIDResponse->NumberOfSigmasTPC(prong,typep);
	    fillthis="hpSigmaVspTPC";
	    ((TH2F*)fOutputMC->FindObject(fillthis))->Fill(momTPC,nsigmap);
	  }
	  if(hasTOF){
	    //Double_t nsigma=fUtilPid->NumberOfSigmasTOF(prong,typep);
	    Double_t nsigma=fPIDResponse->NumberOfSigmasTOF(prong,typep);
	    fillthis="hpSigmaVspTOF";
	    ((TH2F*)fOutputMC->FindObject(fillthis))->Fill(momTOF,nsigma);
	  }

	  if (fReadMC) { // fill hpIDTot ---> PID contamination denominator
	                 //      hIDGood, hnoID ---> PID numerators
	    fillthis="hpIDTot";
	    ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
	    if(pdgProngMC[iprong]==2212) { // id protons
	      fillthis="hpIDGood";
	      ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
	    }
	    else { // misidentified electrons, muons, pions and kaons
	      fillthis="hnopIDp";
	      ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
	    }
	  }
	  break;

	case 321:
	  fillthis="hKTOFSignal";
	  ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(tofSignal);
	  fillthis="hKTPCSignal";
	  ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(dedxTPC);
	  fillthis="hKptProng";
	  ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
	  fillthis="hKd0Prong";
	  ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->Getd0Prong(iprong));
	  fillthis="hKSignalVspTPC";
	  ((TH2F*)fOutputMC->FindObject(fillthis))->Fill(momTPC,dedxTPC);
	  fillthis="hKSignalVspTOF";
	  ((TH2F*)fOutputMC->FindObject(fillthis))->Fill(momTOF,tofSignal);
	  AliPID::EParticleType typek;
	  typek=AliPID::EParticleType(3);
	  if(hasTPC) {
	    //Double_t nsigmap = fUtilPid->NumberOfSigmasTPC(prong,typek);
	    Double_t nsigmap = fPIDResponse->NumberOfSigmasTPC(prong,typek);
	    fillthis="hKSigmaVspTPC";
	    ((TH2F*)fOutputMC->FindObject(fillthis))->Fill(momTPC,nsigmap);
	  }
	  if(hasTOF){
	    //Double_t nsigma = fUtilPid->NumberOfSigmasTOF(prong,typek);
	    Double_t nsigma=fPIDResponse->NumberOfSigmasTOF(prong,typek); // protection against 'old' AODs...
	    fillthis="hKSigmaVspTOF";
	    ((TH2F*)fOutputMC->FindObject(fillthis))->Fill(momTOF,nsigma);
	  }

	  if (fReadMC) { // fill hKIDTot ---> PID contamination denominator
	                 //      hIDGood, hnoID ---> PID numerators
	    fillthis="hKIDTot";
	    ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
	    if(pdgProngMC[iprong]==321) { // id kaons
	      fillthis="hKIDGood";
	      ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
	    }
	    else { // misidentified electrons, muons, pions and protons
	      fillthis="hnokIDk";
	      ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
	    }
	  }
	  break;

	case 211:
	  fillthis="hpiTOFSignal";
	  ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(tofSignal);
	  fillthis="hpiTPCSignal";
	  ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(dedxTPC);
	  fillthis="hpiptProng";
	  ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
	  fillthis="hpid0Prong";
	  ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->Getd0Prong(iprong));
	  fillthis="hpiSignalVspTPC";
	  ((TH2F*)fOutputMC->FindObject(fillthis))->Fill(momTPC,dedxTPC);
	  fillthis="hpiSignalVspTOF";
	  ((TH2F*)fOutputMC->FindObject(fillthis))->Fill(momTOF,tofSignal);
	  AliPID::EParticleType typepi;
	  typepi=AliPID::EParticleType(2);
	  if(hasTPC) {
	    //Double_t nsigmap = fUtilPid->NumberOfSigmasTPC(prong,typepi);
	    Double_t nsigmap = fPIDResponse->NumberOfSigmasTPC(prong,typepi);
	    fillthis="hpiSigmaVspTPC";
	    ((TH2F*)fOutputMC->FindObject(fillthis))->Fill(momTPC,nsigmap);
	  }
	  if(hasTOF){
	    //Double_t nsigma = fUtilPid->NumberOfSigmasTOF(prong,typepi);
	    Double_t nsigma=fPIDResponse->NumberOfSigmasTOF(prong,typepi); // protection against 'old' AODs...
	    fillthis="hpiSigmaVspTOF";
	    ((TH2F*)fOutputMC->FindObject(fillthis))->Fill(momTOF,nsigma);
	  }

	  if (fReadMC) { // fill hpiIDTot ---> PID contamination denominator
	                 //      hIDGood, hnoID ---> PID numerators
	    fillthis="hpiIDTot";
	    ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
	    if(pdgProngMC[iprong]==211) { // id pions
	      fillthis="hpiIDGood";
	      ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
	    }
	    else { // misidentified electrons, muons, kaons and protons
	      fillthis="hnopiIDpi";
	      ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
	    }
	  }
	  break;

	default:
	  break;
	}

      } //end loop on prongs

      //now histograms where we don't need to check identity
      fillthis = "hDist12toPrim";
      ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->GetDist12toPrim());
      ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->GetDist23toPrim());
      fillthis = "hSigmaVert";
      ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->GetSigmaVert());
      fillthis = "hDCAs";
      Double_t dcas[3];
      part->GetDCAs(dcas);
      for (Int_t idca=0;idca<3;idca++) ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(dcas[idca]);
      fillthis = "hCosPointingAngle";
      ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->CosPointingAngle());
      fillthis = "hDecayLength";
      ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->DecayLength());
      Double_t sum2=part->Getd0Prong(0)*part->Getd0Prong(0)+
	part->Getd0Prong(1)*part->Getd0Prong(1)+
	part->Getd0Prong(2)*part->Getd0Prong(2);
      fillthis = "hSum2";
      ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(sum2);
      fillthis = "hptmax";
      ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(ptmax);
      fillthis="hLcpt";
      ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->Pt());

    } // end if (selection)

  } else if( lab<=0 && fReadMC) { // bkg x MC

    fillthis="hbMass";

    if (part->Charge()>0) {      //Lc background
      if ( (pdgProngPID[0]==2212) && (pdgProngPID[1]==321) && (pdgProngPID[2]==211) )
	((TH1F*)fOutputMC->FindObject(fillthis))->Fill(minvLcpKpi); // bkg X MC
      else if ( (pdgProngPID[0]==211) && (pdgProngPID[1]==321) && (pdgProngPID[2]==2212) )
	((TH1F*)fOutputMC->FindObject(fillthis))->Fill(minvLcpiKp); // bkg X MC
    }
    else if (part->Charge()<0){ //anti-Lc background
      if ( (pdgProngPID[0]==2212) && (pdgProngPID[1]==321) && (pdgProngPID[2]==211) )
	((TH1F*)fOutputMC->FindObject(fillthis))->Fill(minvLcpKpi); // bkg X MC
      else if ( (pdgProngPID[0]==211) && (pdgProngPID[1]==321) && (pdgProngPID[2]==2212) )
	((TH1F*)fOutputMC->FindObject(fillthis))->Fill(minvLcpiKp); // bkg X MC
    }


    //apply cut on invariant mass on the pair
    if (TMath::Abs(minvLcpKpi-mPDG)<invmasscut || TMath::Abs(minvLcpiKp-mPDG)<invmasscut) {

      if (selection>0) { // 3prongs candidate x Lc (yes PID)

	Double_t ptmax=0.;
	for (Int_t iprong=0; iprong<3; iprong++) {
	  if(part->PtProng(iprong)>ptmax)ptmax=part->PtProng(iprong);

	  AliAODTrack *prong = (AliAODTrack*)part->GetDaughter(iprong);
	  AliAODPid *pidObjtrk = (AliAODPid*)prong->GetDetPid();
	  AliAODPidHF *pidObj = (AliAODPidHF*)cuts->GetPidHF();
	  Bool_t hasTOF=pidObj->CheckStatus(prong,"TOF");
	  Bool_t hasTPC=pidObj->CheckStatus(prong,"TPC");
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

	  switch (pdgProngPID[iprong]) {
	  case 2212:
	    fillthis="hbpTOFSignal";
	    ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(tofSignal);
	    fillthis="hbpTPCSignal";
	    ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(dedxTPC);
	    fillthis="hbpptProng";
	    ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
	    fillthis="hbpd0Prong";
	    ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->Getd0Prong(iprong));
	    fillthis="hbpSignalVspTPC";
	    ((TH2F*)fOutputMC->FindObject(fillthis))->Fill(momTPC,dedxTPC);
	    fillthis="hbpSignalVspTOF";
	    ((TH2F*)fOutputMC->FindObject(fillthis))->Fill(momTOF,tofSignal);
	    AliPID::EParticleType typep;
	    typep=AliPID::EParticleType(4);
	    if(hasTPC) {
	      //Double_t nsigmap = fUtilPid->NumberOfSigmasTPC(prong,typep);
	      Double_t nsigmap = fPIDResponse->NumberOfSigmasTPC(prong,typep);
	      fillthis="hbpSigmaVspTPC";
	      ((TH2F*)fOutputMC->FindObject(fillthis))->Fill(momTPC,nsigmap);
	    }
	    if(hasTOF){
	      //Double_t nsigma = fUtilPid->NumberOfSigmasTOF(prong,typep);
	      Double_t nsigma=fPIDResponse->NumberOfSigmasTOF(prong,typep); // protection against 'old' AODs...
	      fillthis="hbpSigmaVspTOF";
	      ((TH2F*)fOutputMC->FindObject(fillthis))->Fill(momTOF,nsigma);
	    }
	    if(fReadMC){ // fill hbpIDTot ---> PID contamination denominator
	                 //      hbIDGood, hbnoID ---> PID numerators
	      fillthis="hbpIDTot";
	      ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
	      if(pdgProngMC[iprong]==2212) { // id protons
		fillthis="hbpIDGood";
		((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
	      }
	      else { // misidentified electrons, muons, pions and kaons
		fillthis="hbnopIDp";
		((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
	      }
	    }
	    break;

	  case 321:
	    fillthis="hbKTOFSignal";
	    ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(tofSignal);
	    fillthis="hbKTPCSignal";
	    ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(dedxTPC);
	    fillthis="hbKptProng";
	    ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
	    fillthis="hbKd0Prong";
	    ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->Getd0Prong(iprong));
	    fillthis="hbKSignalVspTPC";
	    ((TH2F*)fOutputMC->FindObject(fillthis))->Fill(momTPC,dedxTPC);
	    fillthis="hbKSignalVspTOF";
	    ((TH2F*)fOutputMC->FindObject(fillthis))->Fill(momTOF,tofSignal);
	    AliPID::EParticleType typek;
	    typek=AliPID::EParticleType(3);
	    if(hasTPC) {
	      //Double_t nsigmap = fUtilPid->NumberOfSigmasTPC(prong,typek);
	      Double_t nsigmap = fPIDResponse->NumberOfSigmasTPC(prong,typek);
	      fillthis="hbKSigmaVspTPC";
	      ((TH2F*)fOutputMC->FindObject(fillthis))->Fill(momTPC,nsigmap);
	    }
	    if(hasTOF){
	      //Double_t nsigma = fUtilPid->NumberOfSigmasTOF(prong,typek);
	      Double_t nsigma=fPIDResponse->NumberOfSigmasTOF(prong,typek); // protection against 'old' AODs...
	      fillthis="hbKSigmaVspTOF";
	      ((TH2F*)fOutputMC->FindObject(fillthis))->Fill(momTOF,nsigma);
	    }
	    if (fReadMC) { // fill hbKIDTot ---> PID contamination denominator
	                   //      hIDGood, hnoID ---> PID numerators
	      fillthis="hbKIDTot";
	      ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
	      if(pdgProngMC[iprong]==321) { // id kaons
		fillthis="hbKIDGood";
		((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
	      }
	      else { // misidentified electrons, muons, pions and protons
		fillthis="hbnokIDk";
		((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
	      }
	    }
	    break;

	  case 211:
	    fillthis="hbpiTOFSignal";
	    ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(tofSignal);
	    fillthis="hbpiTPCSignal";
	    ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(dedxTPC);
	    fillthis="hbpiptProng";
	    ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
	    fillthis="hbpid0Prong";
	    ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->Getd0Prong(iprong));
	    fillthis="hbpiSignalVspTPC";
	    ((TH2F*)fOutputMC->FindObject(fillthis))->Fill(momTPC,dedxTPC);
	    fillthis="hbpiSignalVspTOF";
	    ((TH2F*)fOutputMC->FindObject(fillthis))->Fill(momTOF,tofSignal);
	    AliPID::EParticleType typepi;
	    typepi=AliPID::EParticleType(2);
	    if(hasTPC) {
	      //Double_t nsigmap = fUtilPid->NumberOfSigmasTPC(prong,typepi);
	      Double_t nsigmap = fPIDResponse->NumberOfSigmasTPC(prong,typepi);
	      fillthis="hbpiSigmaVspTPC";
	      ((TH2F*)fOutputMC->FindObject(fillthis))->Fill(momTPC,nsigmap);
	    }
	    if(hasTOF){
	      //Double_t nsigma = fUtilPid->NumberOfSigmasTOF(prong,typepi);
	      Double_t nsigma=fPIDResponse->NumberOfSigmasTOF(prong,typepi); // protection against 'old' AODs...
	      fillthis="hbpiSigmaVspTOF";
	      ((TH2F*)fOutputMC->FindObject(fillthis))->Fill(momTOF,nsigma);
	    }
	    if (fReadMC) { // fill hbpiIDTot ---> PID contamination denominator
	                   //      hIDGood, hnonIDTot ---> PID numerators
	      fillthis="hbpiIDTot";
	      ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
	      if(pdgProngMC[iprong]==211) { // id pions
		fillthis="hbpiIDGood";
		((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
	      }
	      else { // misidentified electrons, muons, kaons and protons
		fillthis="hbnopiIDpi";
		((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->PtProng(iprong));
	      }
	    }
	    break;

	  default:
	    break;
	  }

	} //end loop on prongs

	//now histograms where we don't need to check identity
	fillthis="hbLcpt";
	((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->Pt());

	fillthis = "hbDist12toPrim";
	((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->GetDist12toPrim());
	((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->GetDist23toPrim());
	fillthis = "hbSigmaVert";
	((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->GetSigmaVert());
	fillthis = "hbDCAs";
	Double_t dcas[3];
	part->GetDCAs(dcas);
	for (Int_t idca=0;idca<3;idca++)
	  ((TH1F*)fOutputMC->FindObject(fillthis))->Fill(dcas[idca]);
	fillthis = "hbCosPointingAngle";
	((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->CosPointingAngle());
	fillthis = "hbDecayLength";
	((TH1F*)fOutputMC->FindObject(fillthis))->Fill(part->DecayLength());
	Double_t sum2=part->Getd0Prong(0)*part->Getd0Prong(0)+
	  part->Getd0Prong(1)*part->Getd0Prong(1)+
	  part->Getd0Prong(2)*part->Getd0Prong(2);
	fillthis = "hbSum2";
	((TH1F*)fOutputMC->FindObject(fillthis))->Fill(sum2);
	fillthis = "hbptmax";
	((TH1F*)fOutputMC->FindObject(fillthis))->Fill(ptmax);

      } // end if (selection)

    } // end mass cut

  } // end background case

  return;

}


//--------------
Bool_t AliAnalysisTaskSELambdacUp::CheckOrigin(AliAODRecoDecayHF3Prong *part, TClonesArray *arrayMC){

 for(Int_t i=0;i<3;i++){
  AliAODTrack *daugh=(AliAODTrack*)part->GetDaughter(i);
  Int_t lab=daugh->GetLabel();
  if(lab<0) return 0;
  AliAODMCParticle *mcpart= (AliAODMCParticle*)arrayMC->At(lab);
  Int_t pdgGranma = 0;
  Int_t mother = 0;
  mother = mcpart->GetMother();
  Int_t istep = 0;
  Int_t abspdgGranma =0;
  Bool_t isQuarkFound=kFALSE;
  while (mother >0 ){
                istep++;
                AliDebug(2,Form("mother at step %d = %d", istep, mother));
                AliAODMCParticle* mcGranma = dynamic_cast<AliAODMCParticle*>(arrayMC->At(mother));
                if (mcGranma){
                        pdgGranma = mcGranma->GetPdgCode();
                        AliDebug(2,Form("Pdg mother at step %d = %d", istep, pdgGranma));
                        abspdgGranma = TMath::Abs(pdgGranma);
                        if(abspdgGranma==4 || abspdgGranma==5) {isQuarkFound=kTRUE;return isQuarkFound;}
                        mother = mcGranma->GetMother();
                }else{
                        AliError("Failed casting the mother particle!");
                        break;
                }
        }
  }

 return kFALSE;
}
