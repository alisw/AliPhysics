#include <stdio.h>
#include <stdlib.h>

#include "TH1F.h"
#include "TH2F.h"
#include "TChain.h"
#include "TFile.h"
#include "TParameter.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include "AliESDInputHandler.h"
#include "AliAnaVZEROTrigger.h"
#include "AliCentrality.h"

// VZERO includes
#include "AliESDVZERO.h"

ClassImp(AliAnaVZEROTrigger)

AliAnaVZEROTrigger::AliAnaVZEROTrigger() 
  : AliAnalysisTaskSE("AliAnaVZEROTrigger"), fESD(0), fOutputList(0),
  fMinThr(400.),
  fMaxThr(14000.),
  fRatio(2.0),
  fNThr(30),
  fMBTrigName("CPBI"),
  fUsePhysSel(kFALSE),
  fV0Percent(0),
  fV0PercentAll(0),
  fZvtx(0),
  fXYvtx(0),
  fV0Mult1d(0),
  fV0Charge2d(0),
  fV0Charge2dPercent(0),
  fV0Charge2dAll(0),
  fV0PercentBins(0),
  fV0PercentBinsAll(0),
  fV0Cent(0),
  fV0CentAll(0),
  fV0SemiCent(0),
  fV0SemiCentAll(0),
  fV0CentHw(0),
  fV0CentHwAll(0),
  fV0SemiCentHw(0),
  fV0SemiCentHwAll(0),
  fV0CentTr(0),
  fV0CentTrAll(0),
  fV0SemiCentTr(0),
  fV0SemiCentTrAll(0),
  fV0Percent63(0),
  fV0Percent63All(0),
  fV0MultAll(0),
  fV0Mult63(0),
  fV0CentVsMult(0),
  fV0CentVsCharge(0),
  fV0CentVsTrCharge(0),
  fV0TrChargeVsChargeA(0),
  fV0TrChargeVsChargeC(0),
  fAdcPmt(0),
  fAdcPmt0(0),
  fAdcPmt1(0),
  fAdcPmt2(0),
  fAdcPmt3(0)
{
  // Constructor
  // Init default thr
  fCentCuts[0] = 4152;
  fCentCuts[1] = 8637;
  fSemiCentCuts[0] = 617;
  fSemiCentCuts[1] = 1283;

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliAnaVZEROTrigger::AliAnaVZEROTrigger(const char *name) 
  : AliAnalysisTaskSE(name), fESD(0), fOutputList(0),
  fMinThr(400.),
  fMaxThr(14000.),
  fRatio(2.0),
  fNThr(30),
  fMBTrigName("CPBI"),
  fUsePhysSel(kFALSE),
  fV0Percent(0),
  fV0PercentAll(0),
  fZvtx(0),
  fXYvtx(0),
  fV0Mult1d(0),
  fV0Charge2d(0),
  fV0Charge2dPercent(0),
  fV0Charge2dAll(0),
  fV0PercentBins(0),
  fV0PercentBinsAll(0),
  fV0Cent(0),
  fV0CentAll(0),
  fV0SemiCent(0),
  fV0SemiCentAll(0),
  fV0CentHw(0),
  fV0CentHwAll(0),
  fV0SemiCentHw(0),
  fV0SemiCentHwAll(0),
  fV0CentTr(0),
  fV0CentTrAll(0),
  fV0SemiCentTr(0),
  fV0SemiCentTrAll(0),
  fV0Percent63(0),
  fV0Percent63All(0),
  fV0MultAll(0),
  fV0Mult63(0),
  fV0CentVsMult(0),
  fV0CentVsCharge(0),
  fV0CentVsTrCharge(0),
  fV0TrChargeVsChargeA(0),
  fV0TrChargeVsChargeC(0),
  fAdcPmt(0),
  fAdcPmt0(0),
  fAdcPmt1(0),
  fAdcPmt2(0),
  fAdcPmt3(0)
{
  // Constructor
  // Init default thr
  fCentCuts[0] = 4152;
  fCentCuts[1] = 8637;
  fSemiCentCuts[0] = 617;
  fSemiCentCuts[1] = 1283;

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
void AliAnaVZEROTrigger::Setup(const char *filename)
{
  // Open the config file and
  // read the parameters
  FILE *fin = fopen(filename,"r");
  if (!fin) {
    AliInfo(Form("File %s is not found running with the default parameters.",filename));
  }
  else {
    Int_t res = fscanf(fin,"%f %f %f %d %f %f %f %f",
		       &fMinThr,&fMaxThr,&fRatio,&fNThr,
		       &fSemiCentCuts[0],&fSemiCentCuts[1],&fCentCuts[0],&fCentCuts[1]);
    if(res!=8) {
      AliFatal("Failed to get values from the config file.\n");
    }
    fclose(fin);
  }

  AliInfo(Form("MinThr=%.1f MaxThr=%.1f Ratio=%.4f NThr=%d",fMinThr,fMaxThr,fRatio,fNThr));
  AliInfo(Form("CTA1=%.1f CTC1=%.1f CTA2=%.1f CTC2=%.1f",fSemiCentCuts[0],fSemiCentCuts[1],fCentCuts[0],fCentCuts[1]));
}

//________________________________________________________________________
void AliAnaVZEROTrigger::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);

  fV0Percent    = new TH1F("fV0Percent","Centrality percentile based on v0 mult",400,0,100);
  fOutputList->Add(fV0Percent);
  fV0PercentAll    = new TH1F("fV0PercentAll","Centrality percentile based on v0 mult (no event selection)",400,0,100);
  fOutputList->Add(fV0PercentAll);
  fZvtx    = new TH1F("fZvtx","",500,-20,20);
  fOutputList->Add(fZvtx);
  fXYvtx    = new TH2F("fXYvtx","",250,-0.5,0.5,250,-0.5,0.5);
  fOutputList->Add(fXYvtx);
  fV0Mult1d     = new TH1F("fV0Mult1d","Total v0 mult",500,0,25000);
  fOutputList->Add(fV0Mult1d);

  fV0Charge2d     = new TH2F("fV0Charge2d","V0C vs V0A charge",125,0,25000,125,0,25000);
  fOutputList->Add(fV0Charge2d);
  fV0Charge2dPercent     = new TH2F("fV0Charge2dPercent","V0C vs V0A charge vs centrality percentile",125,0,25000,125,0,25000);
  fOutputList->Add(fV0Charge2dPercent);
  fV0Charge2dAll  = new TH2F("fV0Charge2dAll","V0C vs V0A charge (all events)",125,0,25000,125,0,25000);
  fOutputList->Add(fV0Charge2dAll);

  fV0PercentBins = new TH1F*[fNThr];
  fV0PercentBinsAll = new TH1F*[fNThr];
  for(Int_t j = 0; j < fNThr; ++j) {
    fV0PercentBins[j] = new TH1F(Form("fV0PercentBins_%d",j),"Centrality percentile with V0 charge thresholds",400,0,100);
    fOutputList->Add(fV0PercentBins[j]);
    fV0PercentBinsAll[j] = new TH1F(Form("fV0PercentBinsAll_%d",j),"Centrality percentile with V0 charge thresholds (no event selection)",100,0,100);
    fOutputList->Add(fV0PercentBinsAll[j]);
  }

  fV0Cent = new TH1F("fV0Cent","Centrality percentile with custom V0 charge thresholds",400,0,100);
  fOutputList->Add(fV0Cent);
  fV0CentAll = new TH1F("fV0CentAll","Centrality percentile with custom V0 charge thresholds (no event selection)",100,0,100);
  fOutputList->Add(fV0CentAll);
  fV0SemiCent = new TH1F("fV0SemiCent","Centrality percentile with custom V0 charge thresholds",400,0,100);
  fOutputList->Add(fV0SemiCent);
  fV0SemiCentAll = new TH1F("fV0SemiCentAll","Centrality percentile with custom V0 charge thresholds (no event selection)",100,0,100);
  fOutputList->Add(fV0SemiCentAll);

  fV0CentHw = new TH1F("fV0CentHw","Centrality percentile with hardware V0 charge thresholds",400,0,100);
  fOutputList->Add(fV0CentHw);
  fV0CentHwAll = new TH1F("fV0CentHwAll","Centrality percentile with hardware V0 charge thresholds (no event selection)",100,0,100);
  fOutputList->Add(fV0CentHwAll);
  fV0SemiCentHw = new TH1F("fV0SemiCentHw","Centrality percentile with hardware V0 charge thresholds",400,0,100);
  fOutputList->Add(fV0SemiCentHw);
  fV0SemiCentHwAll = new TH1F("fV0SemiCentHwAll","Centrality percentile with hardware V0 charge thresholds (no event selection)",100,0,100);
  fOutputList->Add(fV0SemiCentHwAll);

  fV0CentTr = new TH1F("fV0CentTr","Centrality percentile with hardware V0 charge thresholds",400,0,100);
  fOutputList->Add(fV0CentTr);
  fV0CentTrAll = new TH1F("fV0CentTrAll","Centrality percentile with hardware V0 charge thresholds (no event selection)",100,0,100);
  fOutputList->Add(fV0CentTrAll);
  fV0SemiCentTr = new TH1F("fV0SemiCentTr","Centrality percentile with hardware V0 charge thresholds",400,0,100);
  fOutputList->Add(fV0SemiCentTr);
  fV0SemiCentTrAll = new TH1F("fV0SemiCentTrAll","Centrality percentile with hardware V0 charge thresholds (no event selection)",100,0,100);
  fOutputList->Add(fV0SemiCentTrAll);

  fV0Percent63    = new TH1F("fV0Percent63","Centrality percentile based on v0 mult",400,0,100);
  fOutputList->Add(fV0Percent63);
  fV0Percent63All    = new TH1F("fV0Percent63All","Centrality percentile based on v0 mult (no event selection)",400,0,100);
  fOutputList->Add(fV0Percent63All);

  fV0MultAll = new TH2F("fV0MultAll","",250,0,25000,250,0,25000);
  fOutputList->Add(fV0MultAll);
  fV0Mult63 = new TH2F("fV0Mult63","",250,0,25000,250,0,25000);
  fOutputList->Add(fV0Mult63);

  fV0CentVsMult = new TH2F("fV0CentVsMult","",250,0,25000,200,0,100);
  fOutputList->Add(fV0CentVsMult);
  fV0CentVsCharge = new TH2F("fV0CentVsCharge","",250,0,100000,200,0,100);
  fOutputList->Add(fV0CentVsCharge); 
  fV0CentVsTrCharge = new TH2F("fV0CentVsTrCharge","",250,0,50000,200,0,100);
  fOutputList->Add(fV0CentVsTrCharge);
  fV0TrChargeVsChargeA = new TH2F("fV0TrChargeVsChargeA","",250,0,75000,250,0,37500);
  fOutputList->Add(fV0TrChargeVsChargeA);
  fV0TrChargeVsChargeC = new TH2F("fV0TrChargeVsChargeC","",250,0,75000,250,0,37500);
  fOutputList->Add(fV0TrChargeVsChargeC);

  fAdcPmt = new TH2F("fAdcPmt","",64,-0.5,63.5,400,0,4000);
  fOutputList->Add(fAdcPmt);
  fAdcPmt0 = new TH2F("fAdcPmt0","",64,-0.5,63.5,400,0,4000);
  fOutputList->Add(fAdcPmt0);
  fAdcPmt1 = new TH2F("fAdcPmt1","",64,-0.5,63.5,400,0,4000);
  fOutputList->Add(fAdcPmt1);
  fAdcPmt2 = new TH2F("fAdcPmt2","",64,-0.5,63.5,400,0,4000);
  fOutputList->Add(fAdcPmt2);
  fAdcPmt3 = new TH2F("fAdcPmt3","",64,-0.5,63.5,400,0,4000);
  fOutputList->Add(fAdcPmt3);

  PostData(1, fOutputList);
 }
//________________________________________________________________________
void AliAnaVZEROTrigger::Init()
{
  // Nothing here
  // ....
}

//________________________________________________________________________
void AliAnaVZEROTrigger::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event


  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!fESD) {
    printf("ERROR: fESD not available\n");
    return;
  }

  AliESDVZERO* esdV0 = fESD->GetVZEROData();
  if (!esdV0) {
    Printf("ERROR: esd V0  not available");
    return;
  }

  // Trigger
  TString trigStr(fESD->GetFiredTriggerClasses());
  if (!trigStr.Contains("-B-")) return;
  if (!trigStr.Contains(fMBTrigName.Data())) return;

  // Count V0 flags
  Int_t nV0A = 0;
  Int_t nV0C = 0;
  for(Int_t i = 0; i < 32; ++i) {
    if (esdV0->GetBBFlag(i)) nV0C++;
    if (esdV0->GetBBFlag(i+32)) nV0A++;
  }

  // Phys sel
  Bool_t goodEvent = kTRUE;
  Bool_t isSelected;
  if (fUsePhysSel)
    isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kAny);
  else
    isSelected = ((esdV0->GetV0ADecision()==1) && (esdV0->GetV0CDecision()==1));

  if (!isSelected) goodEvent = kFALSE;

  const AliESDVertex *primaryVtx = fESD->GetPrimaryVertex();
  if (!primaryVtx) goodEvent = kFALSE;
  if (!primaryVtx->GetStatus()) goodEvent = kFALSE;
  Double_t tPrimaryVtxPosition[3];
  primaryVtx->GetXYZ(tPrimaryVtxPosition);
  if (goodEvent) {
    fZvtx->Fill(tPrimaryVtxPosition[2]);
    fXYvtx->Fill(tPrimaryVtxPosition[0],tPrimaryVtxPosition[1]);
  }
  if (TMath::Abs(tPrimaryVtxPosition[2]) > 10.0) goodEvent = kFALSE;

  if (goodEvent) {
    for(Int_t i = 0; i < 64; ++i) {
      fAdcPmt->Fill(i,esdV0->GetAdc(i));
      if (tPrimaryVtxPosition[2]>5.0) fAdcPmt0->Fill(i,esdV0->GetAdc(i));
      if (tPrimaryVtxPosition[2]>0.0 && tPrimaryVtxPosition[2]<5.0) fAdcPmt1->Fill(i,esdV0->GetAdc(i));
      if (tPrimaryVtxPosition[2]>-5.0 && tPrimaryVtxPosition[2]<0.0) fAdcPmt2->Fill(i,esdV0->GetAdc(i));
      if (tPrimaryVtxPosition[2]<-5.0) fAdcPmt3->Fill(i,esdV0->GetAdc(i));
    }
  }

  UShort_t chargeA = esdV0->GetTriggerChargeA();
  UShort_t chargeC = esdV0->GetTriggerChargeC();

  Float_t offChargeA = 0;
  Float_t offChargeC = 0;
  for(Int_t i = 0; i < 64; ++i) {
    if (i < 32) offChargeC += esdV0->GetAdc(i);
    else offChargeA += esdV0->GetAdc(i);
  }

  if (goodEvent) {
    fV0TrChargeVsChargeA->Fill(offChargeA,chargeA);
    fV0TrChargeVsChargeC->Fill(offChargeC,chargeC);
  }

  Float_t percentile = 0;
  AliCentrality *centrality = fESD->GetCentrality();
  //  if (centrality->GetQuality()) goodEvent = kFALSE;
  percentile = centrality->GetCentralityPercentile("V0M");
  //  if (percentile < 0) percentile = 0;

  if (goodEvent) {
    fV0CentVsMult->Fill(esdV0->GetMTotV0A()+esdV0->GetMTotV0C(),percentile);
    fV0CentVsCharge->Fill(offChargeA+offChargeC,percentile);
    fV0CentVsTrCharge->Fill(chargeA+chargeC,percentile);
  }

  fV0MultAll->Fill(esdV0->GetMTotV0A(),esdV0->GetMTotV0C());
  if ((nV0A+nV0C)>=63)
    fV0Mult63->Fill(esdV0->GetMTotV0A(),esdV0->GetMTotV0C());

  fV0PercentAll->Fill(percentile);
  if (goodEvent) fV0Percent->Fill(percentile);

  if ((nV0A+nV0C)>=63) {
    fV0Percent63All->Fill(percentile);
    if (goodEvent) fV0Percent63->Fill(percentile);
  }

  Float_t multV0 = esdV0->GetMTotV0A()+esdV0->GetMTotV0C();
  if (goodEvent) fV0Mult1d->Fill(multV0);

  fV0Charge2dAll->Fill(chargeA,chargeC);
  if (goodEvent) {
    fV0Charge2d->Fill(chargeA,chargeC);
    fV0Charge2dPercent->Fill(chargeA,chargeC,percentile);
  }

  for(Int_t j = 0; j < fNThr; ++j) {
    Float_t thrA = GetThrA(j);
    Float_t thrC = GetThrC(j);
    if (chargeA >= thrA && chargeC >= thrC) {
      fV0PercentBinsAll[j]->Fill(percentile);
      if (goodEvent) fV0PercentBins[j]->Fill(percentile);
    }
  }

  if (chargeA >= fCentCuts[0] && chargeC >= fCentCuts[1]) {
    fV0CentAll->Fill(percentile);
    if (goodEvent) fV0Cent->Fill(percentile);
  }
  if (chargeA >= fSemiCentCuts[0] && chargeC >= fSemiCentCuts[1]) {
    fV0SemiCentAll->Fill(percentile);
    if (goodEvent) fV0SemiCent->Fill(percentile);
  }

  if (esdV0->GetTriggerBits() & (1<<8)) {
    fV0CentHwAll->Fill(percentile);
    if (goodEvent) fV0CentHw->Fill(percentile);
  }
  if (esdV0->GetTriggerBits() & (1<<6)) {
    fV0SemiCentHwAll->Fill(percentile);
    if (goodEvent) fV0SemiCentHw->Fill(percentile);
  }

  if (trigStr.Contains("CVHN")) {
    fV0CentTrAll->Fill(percentile);
    if (goodEvent) fV0CentTr->Fill(percentile);
  }
  if (trigStr.Contains("CVLN")) {
    fV0SemiCentTrAll->Fill(percentile);
    if (goodEvent) fV0SemiCentTr->Fill(percentile);
  }

  PostData(1, fOutputList);
}      

//________________________________________________________________________
void AliAnaVZEROTrigger::Terminate(Option_t *) 
{
  // Check the output list and store output config file
  // Called once at the end of the query

  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputList) {
    printf("ERROR: Output list not available\n");
    return;
  }

  printf("\n\n");
  FILE *fout=fopen("./trigger.txt","w");
  if (!fout) {
    printf("Failed to open local result file\n");
    return;
  }
  printf("%.0f %.0f %f %d %.0f %.0f %.0f %.0f\n\n",
	 GetMinThr(),
	 GetMaxThr(),
	 GetRatio(),
	 GetNThr(),
	 GetSemiCentCuts(0),GetSemiCentCuts(1),
	 GetCentCuts(0),GetCentCuts(1));

  fprintf(fout,"%.0f %.0f %f %d %.0f %.0f %.0f %.0f\n\n",
	  GetMinThr(),
	  GetMaxThr(),
	  GetRatio(),
	  GetNThr(),
	  GetSemiCentCuts(0),GetSemiCentCuts(1),
	  GetCentCuts(0),GetCentCuts(1));
  fclose(fout);
}

//________________________________________________________________________
Float_t AliAnaVZEROTrigger::GetThrA(Int_t j) const
{
  // Get the threshold on A side with
  // index i
  if (j<0 || j>= fNThr) return 0;

  return (fMinThr + ((Float_t)j)*(fMaxThr-fMinThr)/((Float_t)fNThr-1.));
}

//________________________________________________________________________
Float_t AliAnaVZEROTrigger::GetThrC(Int_t j) const
{
  // Get the threshold on C side with
  // index i
  if (j<0 || j>= fNThr) return 0;

  Float_t thrA = (fMinThr + ((Float_t)j)*(fMaxThr-fMinThr)/((Float_t)fNThr-1.));
  return (thrA*fRatio);
}

