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

//*************************************************************************
// Class AliAnalysisTaskSEDplus
// AliAnalysisTaskSE for the D+ candidates Invariant Mass Histogram and 
//comparison of heavy-flavour decay candidates
// to MC truth (kinematics stored in the AOD)
// Authors: Renu Bala, bala@to.infn.it
// F. Prino, prino@to.infn.it
// G. Ortona, ortona@to.infn.it
/////////////////////////////////////////////////////////////

#include <TClonesArray.h>
#include <TNtuple.h>
#include <TCanvas.h>
#include <TList.h>
#include <TString.h>
#include <TDatabasePDG.h>

#include "AliAnalysisManager.h"
#include "AliRDHFCutsDplustoKpipi.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliAnalysisVertexingHF.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskSEDplus.h"
#include "AliNormalizationCounter.h"

ClassImp(AliAnalysisTaskSEDplus)


//________________________________________________________________________
AliAnalysisTaskSEDplus::AliAnalysisTaskSEDplus():
AliAnalysisTaskSE(),
  fOutput(0), 
  fHistNEvents(0),
  fPtVsMass(0),
  fPtVsMassTC(0),
  fYVsPt(0),
  fYVsPtTC(0),
  fYVsPtSig(0),
  fYVsPtSigTC(0),
  fNtupleDplus(0),
  fUpmasslimit(1.965),
  fLowmasslimit(1.765),
  fNPtBins(0),
  fBinWidth(0.002),
  fListCuts(0),
  fRDCutsProduction(0),
  fRDCutsAnalysis(0),
  fCounter(0),
  fFillNtuple(kFALSE),
  fReadMC(kFALSE),
  fUseStrangeness(kFALSE),
  fUseBit(kTRUE),
  fCutsDistr(kFALSE),
  fDoImpPar(kFALSE),
  fNImpParBins(400),
  fLowerImpPar(-2000.),
  fHigherImpPar(2000.),
  fDoLS(0)
{
   // Default constructor
}

//________________________________________________________________________
AliAnalysisTaskSEDplus::AliAnalysisTaskSEDplus(const char *name,AliRDHFCutsDplustoKpipi *dpluscutsana,AliRDHFCutsDplustoKpipi *dpluscutsprod,Bool_t fillNtuple):
  AliAnalysisTaskSE(name),
  fOutput(0),
  fHistNEvents(0),
  fPtVsMass(0),
  fPtVsMassTC(0),
  fYVsPt(0),
  fYVsPtTC(0),
  fYVsPtSig(0),
  fYVsPtSigTC(0),
  fNtupleDplus(0),
  fUpmasslimit(1.965),
  fLowmasslimit(1.765),
  fNPtBins(0),
  fBinWidth(0.002),
  fListCuts(0),
  fRDCutsProduction(dpluscutsprod),
  fRDCutsAnalysis(dpluscutsana),
  fCounter(0),
  fFillNtuple(fillNtuple),
  fReadMC(kFALSE),
  fUseStrangeness(kFALSE),
  fUseBit(kTRUE),
  fCutsDistr(kFALSE),
  fDoImpPar(kFALSE),
  fNImpParBins(400),
  fLowerImpPar(-2000.),
  fHigherImpPar(2000.),
  fDoLS(0)
{
  // 
  // Standrd constructor
  //
  fNPtBins=fRDCutsAnalysis->GetNPtBins();
  // Default constructor
   // Output slot #1 writes into a TList container
  DefineOutput(1,TList::Class());  //My private output
 // Output slot #2 writes cut to private output
  //  DefineOutput(2,AliRDHFCutsDplustoKpipi::Class());
  DefineOutput(2,TList::Class());
// Output slot #3 writes cut to private output
  DefineOutput(3,AliNormalizationCounter::Class());

  if(fFillNtuple){
    // Output slot #4 writes into a TNtuple container
    DefineOutput(4,TNtuple::Class());  //My private output
  }
}

//________________________________________________________________________
AliAnalysisTaskSEDplus::~AliAnalysisTaskSEDplus()
{
  //
  // Destructor
  //
  if (fOutput) {
    delete fOutput;
    fOutput = 0;
  }
  if(fHistNEvents){
    delete fHistNEvents;
    fHistNEvents=0;
  }  
  for(Int_t i=0;i<3;i++){
    if(fHistCentrality[i]){delete fHistCentrality[i]; fHistCentrality[i]=0;}
  }

  for(Int_t i=0;i<3*fNPtBins;i++){
    if(fMassHist[i]){ delete fMassHist[i]; fMassHist[i]=0;}
    if(fCosPHist[i]){ delete fCosPHist[i]; fCosPHist[i]=0;}
    if(fDLenHist[i]){ delete fDLenHist[i]; fDLenHist[i]=0;}
    if(fSumd02Hist[i]){ delete fSumd02Hist[i]; fSumd02Hist[i]=0;}
    if(fSigVertHist[i]){ delete fSigVertHist[i]; fSigVertHist[i]=0;}
    if(fPtMaxHist[i]){ delete fPtMaxHist[i]; fPtMaxHist[i]=0;}
    if(fPtKHist[i]){ delete fPtKHist[i]; fPtKHist[i]=0;}
    if(fPtpi1Hist[i]){ delete fPtpi1Hist[i]; fPtpi1Hist[i]=0;}
    if(fPtpi2Hist[i]){ delete fPtpi2Hist[i]; fPtpi2Hist[i]=0;}
    if(fDCAHist[i]){ delete fDCAHist[i]; fDCAHist[i]=0;}
    if(fMassHistTC[i]){ delete fMassHistTC[i]; fMassHistTC[i]=0;}
    if(fMassHistTCPlus[i]){ delete fMassHistTCPlus[i]; fMassHistTCPlus[i]=0;}
    if(fMassHistTCMinus[i]){ delete fMassHistTCMinus[i]; fMassHistTCMinus[i]=0;}

    if(fDLxy[i]){delete fDLxy[i]; fDLxy[i]=0;}
    if(fDLxyTC[i]){delete fDLxyTC[i]; fDLxyTC[i]=0;}
    if(fCosxy[i]){delete fCosxy[i]; fCosxy[i]=0;}
    if(fCosxyTC[i]){delete fCosxyTC[i]; fCosxyTC[i]=0;}
    if(fMassHistLS[i]){ delete fMassHistLS[i]; fMassHistLS[i]=0;}
    if(fCosPHistLS[i]){ delete fCosPHistLS[i]; fCosPHistLS[i]=0;}
    if(fDLenHistLS[i]){ delete fDLenHistLS[i]; fDLenHistLS[i]=0;}
    if(fSumd02HistLS[i]){ delete fSumd02HistLS[i]; fSumd02HistLS[i]=0;}
    if(fSigVertHistLS[i]){ delete fSigVertHistLS[i]; fSigVertHistLS[i]=0;}
    if(fPtMaxHistLS[i]){ delete fPtMaxHistLS[i]; fPtMaxHistLS[i]=0;}
    if(fDCAHistLS[i]){ delete fDCAHistLS[i]; fDCAHistLS[i]=0;}
    if(fMassHistLSTC[i]){ delete fMassHistLSTC[i]; fMassHistLSTC[i]=0;}
  }

  for(Int_t i=0;i<3;i++){
    if(fCorreld0Kd0pi[i]){ delete fCorreld0Kd0pi[i]; fCorreld0Kd0pi[i]=0;}
  }

  if(fPtVsMass){
    delete fPtVsMass;
    fPtVsMass=0;
  }
  if(fPtVsMassTC){
    delete fPtVsMassTC;
    fPtVsMassTC=0;
  }
  if(fYVsPt){
    delete fYVsPt;
    fYVsPt=0;
  }
  if(fYVsPtTC){
    delete fYVsPtTC;
    fYVsPtTC=0;
  }
  if(fYVsPtSig){
    delete fYVsPtSig;
    fYVsPtSig=0;
  }
  if(fYVsPtSigTC){
    delete fYVsPtSigTC;
    fYVsPtSigTC=0;
  }
  if(fNtupleDplus){
    delete fNtupleDplus;
    fNtupleDplus=0;
  }
  if (fListCuts) {
    delete fListCuts;
    fListCuts = 0;
  }
  if(fRDCutsProduction){
    delete fRDCutsProduction;
    fRDCutsProduction = 0;
  }
  if(fRDCutsAnalysis){
    delete fRDCutsAnalysis;
    fRDCutsAnalysis = 0;
  }
  for(Int_t i=0; i<5; i++){
    delete fHistMassPtImpParTC[i];
  }
  if(fCounter){
    delete fCounter;
    fCounter = 0;
  }


}  
//_________________________________________________________________
void  AliAnalysisTaskSEDplus::SetMassLimits(Float_t range){
  // set invariant mass limits
  Float_t bw=GetBinWidth();
  fUpmasslimit = 1.865+range;
  fLowmasslimit = 1.865-range;
  SetBinWidth(bw);
}
//_________________________________________________________________
void  AliAnalysisTaskSEDplus::SetMassLimits(Float_t lowlimit, Float_t uplimit){
  // set invariant mass limits
  if(uplimit>lowlimit)
    {
      Float_t bw=GetBinWidth();
      fUpmasslimit = lowlimit;
      fLowmasslimit = uplimit;
      SetBinWidth(bw);
    }
}
//________________________________________________________________
void AliAnalysisTaskSEDplus::SetBinWidth(Float_t w){
  Float_t width=w;
  Int_t nbins=(Int_t)((fUpmasslimit-fLowmasslimit)/width+0.5);
  Int_t missingbins=4-nbins%4;
  nbins=nbins+missingbins;
  width=(fUpmasslimit-fLowmasslimit)/nbins;
  if(missingbins!=0){
    printf("AliAnalysisTaskSEDplus::SetBinWidth: W-bin width of %f will produce histograms not rebinnable by 4. New width set to %f\n",w,width);
  }
  else{
    if(fDebug>1) printf("AliAnalysisTaskSEDplus::SetBinWidth: width set to %f\n",width);
  }
  fBinWidth=width;
}
//_________________________________________________________________
Int_t AliAnalysisTaskSEDplus::GetNBinsHistos(){
  return (Int_t)((fUpmasslimit-fLowmasslimit)/fBinWidth+0.5);
}
//_________________________________________________________________
void AliAnalysisTaskSEDplus::LSAnalysis(TClonesArray *arrayOppositeSign,TClonesArray *arrayLikeSign,AliAODEvent *aod,AliAODVertex *vtx1, Int_t nDplusOS){
  //
  //
  // Fill the Like Sign histograms
  //
  if(fDebug>1)printf("started LS\n");
  
  //histograms for like sign
  Int_t nbins=GetNBinsHistos();;
  TH1F *histLSPlus = new TH1F("LSPlus","LSPlus",nbins,fLowmasslimit,fUpmasslimit);
  TH1F *histLSMinus = new TH1F("LSMinus","LSMinus",nbins,fLowmasslimit,fUpmasslimit);
  
  Int_t nPosTrks=0,nNegTrks=0;
  Int_t nOStriplets = arrayOppositeSign->GetEntriesFast();
  Int_t nDplusLS=0;
  Int_t nDminusLS=0;
  Int_t nLikeSign = arrayLikeSign->GetEntriesFast();
  Int_t index=0; 
  
  // loop over like sign candidates
  for(Int_t iLikeSign = 0; iLikeSign < nLikeSign; iLikeSign++) {
    AliAODRecoDecayHF3Prong *d = (AliAODRecoDecayHF3Prong*)arrayLikeSign->UncheckedAt(iLikeSign);
    if(fUseBit && !d->HasSelectionBit(AliRDHFCuts::kDplusCuts))continue;
    Bool_t unsetvtx=kFALSE;
    if(!d->GetOwnPrimaryVtx()) {
      d->SetOwnPrimaryVtx(vtx1); // needed to compute all variables
      unsetvtx=kTRUE;
    }
    Double_t ptCand = d->Pt();
    Int_t iPtBin = fRDCutsAnalysis->PtBin(ptCand);
    if(iPtBin<0)continue;
    Int_t sign= d->GetCharge();
    if(sign>0){
      nPosTrks++;
    }else{
      nNegTrks++;
    }
    //    if(fRDCutsProduction->IsSelected(d,AliRDHFCuts::kCandidate,aod)){
    fRDCutsAnalysis->IsSelected(d,AliRDHFCuts::kCandidate,aod);
    Int_t passTightCuts=fRDCutsAnalysis->GetIsSelectedCuts();
    if(passTightCuts>0){
      Float_t invMass = d->InvMassDplus();
      index=GetLSHistoIndex(iPtBin);
      fMassHistLS[index+1]->Fill(invMass);
      if(sign>0){
	histLSPlus->Fill(invMass);
	nDplusLS++;
      }else{
	histLSMinus->Fill(invMass);
	nDminusLS++;
      }
      if(fCutsDistr){
	Double_t dlen=d->DecayLength();
	Double_t cosp=d->CosPointingAngle();
	Double_t sumD02=d->Getd0Prong(0)*d->Getd0Prong(0)+d->Getd0Prong(1)*d->Getd0Prong(1)+d->Getd0Prong(2)*d->Getd0Prong(2);
	Double_t dca=d->GetDCA();   
	Double_t sigvert=d->GetSigmaVert();   
	Double_t ptmax=0;
	for(Int_t i=0;i<3;i++){
	  if(d->PtProng(i)>ptmax)ptmax=d->PtProng(i);
	}
	fCosPHistLS[iPtBin]->Fill(cosp);
	fDLenHistLS[iPtBin]->Fill(dlen);
	fSumd02HistLS[iPtBin]->Fill(sumD02);
	fSigVertHistLS[iPtBin]->Fill(sigvert);
	fPtMaxHistLS[iPtBin]->Fill(ptmax);
	fDCAHistLS[iPtBin]->Fill(dca);
      }
    }
    if(unsetvtx) d->UnsetOwnPrimaryVtx();
  }
  //wei:normalized to the number of combinations (production)
  //wei2:normalized to the number of  LS/OS (production)
  //wei3:normalized to the number of  LS/OS (analysis)
  //wei4:normalized to the number of combinations (analysis)
  Float_t wei2=0;
  if(nLikeSign!=0)wei2 = (Float_t)nOStriplets/(Float_t)nLikeSign;
  Float_t wei3=0;
  if(nDplusLS!=0)wei3 = (Float_t)nDplusOS/(Float_t)(nDplusLS+nDminusLS);
  Float_t weiplus=1.,weiminus=1.;
  Float_t wei4plus=1.,wei4minus=1.;
  //wei* should be automatically protected, since to get a triplet there must be at least 3 good tracks in the event
  if(nPosTrks>2)weiplus=3.*(Float_t)nNegTrks/((Float_t)nPosTrks-2.);
  if(nDplusLS>2)wei4plus=3.*(Float_t)nDminusLS/((Float_t)nDplusLS-2.);
  if(nNegTrks>2)weiminus=3.*(Float_t)nPosTrks/((Float_t)nNegTrks-2.);
  if(nDminusLS>2)wei4minus=3.*(Float_t)nDplusLS/((Float_t)nDminusLS-2.);

  fMassHistLS[index]->Add(histLSPlus,weiplus);
  fMassHistLS[index]->Add(histLSMinus,weiminus);
  fMassHistLS[index+2]->Add(histLSPlus,wei2);
  fMassHistLS[index+2]->Add(histLSMinus,wei2);
  fMassHistLS[index+3]->Add(histLSPlus,wei3);
  fMassHistLS[index+3]->Add(histLSMinus,wei3);
  fMassHistLS[index+4]->Add(histLSPlus,wei4plus);
  fMassHistLS[index+4]->Add(histLSMinus,wei4minus);

  delete histLSPlus;histLSPlus=0;
  delete histLSMinus;histLSMinus=0;
  
  if(fDebug>1) printf("LS analysis terminated\n");  
}

//__________________________________________
void AliAnalysisTaskSEDplus::Init(){
  //
  // Initialization
  //
  if(fDebug > 1) printf("AnalysisTaskSEDplus::Init() \n");
  
  //PostData(2,fRDCutsloose);//we should then put those cuts in a tlist if we have more than 1
  fListCuts=new TList();
  AliRDHFCutsDplustoKpipi *production = new AliRDHFCutsDplustoKpipi(*fRDCutsProduction);
  production->SetName("ProductionCuts");
  AliRDHFCutsDplustoKpipi *analysis = new AliRDHFCutsDplustoKpipi(*fRDCutsAnalysis);
  analysis->SetName("AnalysisCuts");
  
  fListCuts->Add(production);
  fListCuts->Add(analysis);
  PostData(2,fListCuts);
  
  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEDplus::UserCreateOutputObjects()
{
  // Create the output container
  //
  if(fDebug > 1) printf("AnalysisTaskSEDplus::UserCreateOutputObjects() \n");

  // Several histograms are more conveniently managed in a TList
  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("OutputHistos");

  TString hisname;
  Int_t index=0;
  Int_t nbins=GetNBinsHistos();
  fHistCentrality[0]=new TH1F("centrality","centrality",100,0.5,30000.5);
  fHistCentrality[1]=new TH1F("centrality(selectedCent)","centrality(selectedCent)",100,0.5,30000.5);
  fHistCentrality[2]=new TH1F("centrality(OutofCent)","centrality(OutofCent)",100,0.5,30000.5);
  for(Int_t i=0;i<3;i++){
    fHistCentrality[i]->Sumw2();
    fOutput->Add(fHistCentrality[i]);
  }
  for(Int_t i=0;i<fNPtBins;i++){

    index=GetHistoIndex(i);

    hisname.Form("hMassPt%d",i);
    fMassHist[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHist[index]->Sumw2();
    hisname.Form("hCosPAllPt%d",i);
    fCosPHist[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.5,1.);
    fCosPHist[index]->Sumw2();
    hisname.Form("hDLenAllPt%d",i);
    fDLenHist[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.,0.5);
    fDLenHist[index]->Sumw2();
    hisname.Form("hSumd02AllPt%d",i);
    fSumd02Hist[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.,1.);
    fSumd02Hist[index]->Sumw2();
    hisname.Form("hSigVertAllPt%d",i);
    fSigVertHist[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.,0.1);
    fSigVertHist[index]->Sumw2();
    hisname.Form("hPtMaxAllPt%d",i);
    fPtMaxHist[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.5,5.);
    fPtMaxHist[index]->Sumw2();
    hisname.Form("hPtKPt%d",i);
    fPtKHist[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.5,5.);
    fPtKHist[index]->Sumw2();
    hisname.Form("hPtpi1Pt%d",i);
    fPtpi1Hist[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.5,5.);
    fPtpi1Hist[index]->Sumw2();
    hisname.Form("hPtpi2Pt%d",i);
    fPtpi2Hist[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.5,5.);
    fPtpi2Hist[index]->Sumw2();
    hisname.Form("hDCAAllPt%d",i);
    fDCAHist[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.,0.1);
    fDCAHist[index]->Sumw2();

    hisname.Form("hDLxyPt%d",i);
    fDLxy[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.,10.);
    fDLxy[index]->Sumw2();
    hisname.Form("hCosxyPt%d",i);
    fCosxy[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.85,1.);
    fCosxy[index]->Sumw2();
    hisname.Form("hDLxyPt%dTC",i);
    fDLxyTC[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.,10.);
    fDLxyTC[index]->Sumw2();
    hisname.Form("hCosxyPt%dTC",i);
    fCosxyTC[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.85,1.);
    fCosxyTC[index]->Sumw2();
   
    hisname.Form("hMassPt%dTC",i);
    fMassHistTC[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHistTC[index]->Sumw2();
    hisname.Form("hMassPt%dTCPlus",i);
    fMassHistTCPlus[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHistTCPlus[index]->Sumw2();
    hisname.Form("hMassPt%dTCMinus",i);
    fMassHistTCMinus[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHistTCMinus[index]->Sumw2();


       
    index=GetSignalHistoIndex(i);    
    hisname.Form("hSigPt%d",i);
    fMassHist[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHist[index]->Sumw2();
    hisname.Form("hCosPSigPt%d",i);
    fCosPHist[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.5,1.);
    fCosPHist[index]->Sumw2();
    hisname.Form("hDLenSigPt%d",i);
    fDLenHist[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.,0.5);
    fDLenHist[index]->Sumw2();
    hisname.Form("hSumd02SigPt%d",i);
    fSumd02Hist[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.,1.);
    fSumd02Hist[index]->Sumw2();
    hisname.Form("hSigVertSigPt%d",i);
    fSigVertHist[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.,0.1);
    fSigVertHist[index]->Sumw2();
    hisname.Form("hPtMaxSigPt%d",i);
    fPtMaxHist[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.5,5.);
    fPtMaxHist[index]->Sumw2();  
    hisname.Form("hPtKSigPt%d",i);  
    fPtKHist[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.5,5.);
    fPtKHist[index]->Sumw2();
    hisname.Form("hPtpi1SigPt%d",i);
    fPtpi1Hist[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.5,5.);
    fPtpi1Hist[index]->Sumw2();
    hisname.Form("hPtpi2SigPt%d",i);
    fPtpi2Hist[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.5,5.);
    fPtpi2Hist[index]->Sumw2();

    hisname.Form("hDCASigPt%d",i);
    fDCAHist[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.,0.1);
    fDCAHist[index]->Sumw2();    

    hisname.Form("hDLxySigPt%d",i);
    fDLxy[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.,10.);
    fDLxy[index]->Sumw2();
    hisname.Form("hCosxySigPt%d",i);
    fCosxy[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.85,1.);
    fCosxy[index]->Sumw2();
    hisname.Form("hDLxySigPt%dTC",i);
    fDLxyTC[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.,10.);
    fDLxyTC[index]->Sumw2();
    hisname.Form("hCosxySigPt%dTC",i);
    fCosxyTC[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.85,1.);
    fCosxyTC[index]->Sumw2();
    hisname.Form("hSigPt%dTC",i);
    fMassHistTC[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHistTC[index]->Sumw2();
    hisname.Form("hSigPt%dTCPlus",i);
    fMassHistTCPlus[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHistTCPlus[index]->Sumw2();
    hisname.Form("hSigPt%dTCMinus",i);
    fMassHistTCMinus[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHistTCMinus[index]->Sumw2();


    index=GetBackgroundHistoIndex(i); 
    hisname.Form("hBkgPt%d",i);
    fMassHist[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHist[index]->Sumw2();
    hisname.Form("hCosPBkgPt%d",i);
    fCosPHist[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.5,1.);
    fCosPHist[index]->Sumw2();
    hisname.Form("hDLenBkgPt%d",i);
    fDLenHist[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.,0.5);
    fDLenHist[index]->Sumw2();
    hisname.Form("hSumd02BkgPt%d",i);
    fSumd02Hist[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.,1.);
    fSumd02Hist[index]->Sumw2();
    hisname.Form("hSigVertBkgPt%d",i);
    fSigVertHist[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.,0.1);
    fSigVertHist[index]->Sumw2();
    hisname.Form("hPtMaxBkgPt%d",i);
    fPtMaxHist[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.5,5.);
    fPtMaxHist[index]->Sumw2();
    hisname.Form("hPtKBkgPt%d",i);  
    fPtKHist[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.5,5.);
    fPtKHist[index]->Sumw2();
    hisname.Form("hPtpi1BkgPt%d",i);
    fPtpi1Hist[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.5,5.);
    fPtpi1Hist[index]->Sumw2();
    hisname.Form("hPtpi2BkgPt%d",i);
    fPtpi2Hist[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.5,5.);
    fPtpi2Hist[index]->Sumw2();
    hisname.Form("hDCABkgPt%d",i);
    fDCAHist[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.,0.1);
    fDCAHist[index]->Sumw2();

    hisname.Form("hDLxyBkgPt%d",i);
    fDLxy[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.,10.);
    fDLxy[index]->Sumw2();
    hisname.Form("hCosxyBkgPt%d",i);
    fCosxy[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.85,1.);
    fCosxy[index]->Sumw2();
    hisname.Form("hDLxyBkgPt%dTC",i);
    fDLxyTC[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.,10.);
    fDLxyTC[index]->Sumw2();
    hisname.Form("hCosxyBkgPt%dTC",i);
    fCosxyTC[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.85,1.);
    fCosxyTC[index]->Sumw2();
  

    hisname.Form("hBkgPt%dTC",i);
    fMassHistTC[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHistTC[index]->Sumw2();
    hisname.Form("hBkgPt%dTCPlus",i);
    fMassHistTCPlus[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHistTCPlus[index]->Sumw2();
    hisname.Form("hBkgPt%dTCMinus",i);
    fMassHistTCMinus[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHistTCMinus[index]->Sumw2();
  }
    

  for(Int_t i=0; i<3*fNPtBins; i++){
    fOutput->Add(fMassHist[i]);
    if(fCutsDistr){  
      fOutput->Add(fCosPHist[i]);
      fOutput->Add(fDLenHist[i]);
      fOutput->Add(fSumd02Hist[i]);
      fOutput->Add(fSigVertHist[i]);
      fOutput->Add(fPtMaxHist[i]);
      fOutput->Add(fPtKHist[i]);
      fOutput->Add(fPtpi1Hist[i]);
      fOutput->Add(fPtpi2Hist[i]);
      fOutput->Add(fDCAHist[i]);
      fOutput->Add(fDLxy[i]);
      fOutput->Add(fDLxyTC[i]);
      fOutput->Add(fCosxy[i]);
      fOutput->Add(fCosxyTC[i]);
    }  
    fOutput->Add(fMassHistTC[i]);
    fOutput->Add(fMassHistTCPlus[i]);
    fOutput->Add(fMassHistTCMinus[i]);

  }
  
  if(fCutsDistr){  
    fCorreld0Kd0pi[0]=new TH2F("hCorreld0Kd0piAll","",100,-0.02,0.02,100,-0.02,0.02);
    fCorreld0Kd0pi[1]=new TH2F("hCorreld0Kd0piSig","",100,-0.02,0.02,100,-0.02,0.02);
    fCorreld0Kd0pi[2]=new TH2F("hCorreld0Kd0piBkg","",100,-0.02,0.02,100,-0.02,0.02);
    for(Int_t i=0; i<3; i++){
      fCorreld0Kd0pi[i]->Sumw2();
      fOutput->Add(fCorreld0Kd0pi[i]);
    }
  }
  
  fHistNEvents = new TH1F("fHistNEvents", "number of events ",9,-0.5,8.5);
  fHistNEvents->GetXaxis()->SetBinLabel(1,"nEventsAnal");
  fHistNEvents->GetXaxis()->SetBinLabel(2,"nEvents with good vertex");
  fHistNEvents->GetXaxis()->SetBinLabel(3,"nEvents with PbPb HM trigger");
  fHistNEvents->GetXaxis()->SetBinLabel(4,"no. of Rejected pileup events");
  fHistNEvents->GetXaxis()->SetBinLabel(5,"no. of candidate");
  fHistNEvents->GetXaxis()->SetBinLabel(6,"no. of D+ after loose cuts");
  fHistNEvents->GetXaxis()->SetBinLabel(7,"no. of D+ after tight cuts");
  fHistNEvents->GetXaxis()->SetBinLabel(8,"no. of out centrality events");
  fHistNEvents->GetXaxis()->SetBinLabel(9,"no. of cand wo bitmask");
  fHistNEvents->GetXaxis()->SetNdivisions(1,kFALSE);  
  fHistNEvents->Sumw2();
  fHistNEvents->SetMinimum(0);
  fOutput->Add(fHistNEvents);

  fPtVsMass=new TH2F("hPtVsMass","PtVsMass (prod. cuts)",nbins,fLowmasslimit,fUpmasslimit,200,0.,20.);
  fPtVsMassTC=new TH2F("hPtVsMassTC","PtVsMass (analysis cuts)",nbins,fLowmasslimit,fUpmasslimit,200,0.,20.);  
  fYVsPt=new TH2F("hYVsPt","YvsPt (prod. cuts)",40,0.,20.,80,-2.,2.);
  fYVsPtTC=new TH2F("hYVsPtTC","YvsPt (analysis cuts)",40,0.,20.,80,-2.,2.);
  fYVsPtSig=new TH2F("hYVsPtSig","YvsPt (MC, only sig., prod. cuts)",40,0.,20.,80,-2.,2.);
  fYVsPtSigTC=new TH2F("hYVsPtSigTC","YvsPt (MC, only Sig, analysis cuts)",40,0.,20.,80,-2.,2.);

  fOutput->Add(fPtVsMass);
  fOutput->Add(fPtVsMassTC);
  fOutput->Add(fYVsPt);
  fOutput->Add(fYVsPtTC);
  fOutput->Add(fYVsPtSig);
  fOutput->Add(fYVsPtSigTC);


  //Counter for Normalization
  TString normName="NormalizationCounter";
  AliAnalysisDataContainer *cont = GetOutputSlot(3)->GetContainer();
  if(cont)normName=(TString)cont->GetName();
  fCounter = new AliNormalizationCounter(normName.Data());
  fCounter->Init();

  if(fDoLS) CreateLikeSignHistos();
  if(fDoImpPar) CreateImpactParameterHistos();

  if(fFillNtuple){
    OpenFile(4); // 4 is the slot number of the ntuple
   
    fNtupleDplus = new TNtuple("fNtupleDplus","D +","pdg:Px:Py:Pz:PtTrue:VxTrue:VyTrue:VzTrue:Ptpi:PtK:Ptpi2:PtRec:PointingAngle:DecLeng:VxRec:VyRec:VzRec:InvMass:sigvert:d0Pi:d0K:d0Pi2:dca:d0square");  
    
  }
  
  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEDplus::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event:
  // heavy flavor candidates association to MC truth

  AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());
  
 

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

  if(!aod || !array3Prong) {
    printf("AliAnalysisTaskSEDplus::UserExec: Charm3Prong branch not found!\n");
    return;
  }
  if(!arrayLikeSign && fDoLS) {
    printf("AliAnalysisTaskSEDplus::UserExec: LikeSign3Prong branch not found!\n");
    return;
  }


  // fix for temporary bug in ESDfilter 
  // the AODs with null vertex pointer didn't pass the PhysSel
  if(!aod->GetPrimaryVertex()||TMath::Abs(aod->GetMagneticField())<0.001) return;
  fCounter->StoreEvent(aod,fRDCutsAnalysis,fReadMC);
  fHistNEvents->Fill(0); // count event

  Bool_t isEvSel=fRDCutsAnalysis->IsEventSelected(aod);
  Bool_t isEvSelP=kTRUE;
  isEvSelP=fRDCutsProduction->IsEventSelected(aod); // to have proper PID object settings

  Float_t centrality=aod->GetNTracks();//fRDCutsAnalysis->GetCentrality(aod);
  fHistCentrality[0]->Fill(centrality);
  // trigger class for PbPb C0SMH-B-NOPF-ALLNOTRD
  TString trigclass=aod->GetFiredTriggerClasses();
  if(trigclass.Contains("C0SMH-B-NOPF-ALLNOTRD")||trigclass.Contains("C0SMH-B-NOPF-ALL")) fHistNEvents->Fill(2);
  if(fRDCutsAnalysis->GetWhyRejection()==1)fHistNEvents->Fill(3); 
  if(fRDCutsAnalysis->GetWhyRejection()==2){fHistNEvents->Fill(7);fHistCentrality[2]->Fill(centrality);}

  // Post the data already here  
  PostData(1,fOutput);
  if(!isEvSel)return;

  fHistCentrality[1]->Fill(centrality);
  fHistNEvents->Fill(1);

  TClonesArray *arrayMC=0;
  AliAODMCHeader *mcHeader=0;

  // AOD primary vertex
  AliAODVertex *vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();
  //    vtx1->Print();
   TString primTitle = vtx1->GetTitle();
   //if(primTitle.Contains("VertexerTracks") && vtx1->GetNContributors()>0)fHistNEvents->Fill(2);
 
  // load MC particles
  if(fReadMC){
    
    arrayMC =  (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
    if(!arrayMC) {
      printf("AliAnalysisTaskSEDplus::UserExec: MC particles branch not found!\n");
      return;
    }
    
  // load MC header
    mcHeader =  (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if(!mcHeader) {
    printf("AliAnalysisTaskSEDplus::UserExec: MC header branch not found!\n");
    return;
    }
  }
  
  Int_t n3Prong = array3Prong->GetEntriesFast();
  //  printf("Number of D+->Kpipi: %d and of tracks: %d\n",n3Prong,aod->GetNumberOfTracks());
  
  
  Int_t nOS=0;
  Int_t index;
  Int_t pdgDgDplustoKpipi[3]={321,211,211};

  if(fDoLS>1){
    for (Int_t i3Prong = 0; i3Prong < n3Prong; i3Prong++) {
      AliAODRecoDecayHF3Prong *d = (AliAODRecoDecayHF3Prong*)array3Prong->UncheckedAt(i3Prong);
      if(fUseBit && !d->HasSelectionBit(AliRDHFCuts::kDplusCuts)){
	if(fRDCutsAnalysis->IsSelected(d,AliRDHFCuts::kCandidate,aod))nOS++;
      }
    }
  }else{
  // Double_t cutsDplus[12]={0.2,0.4,0.4,0.,0.,0.01,0.06,0.02,0.,0.85,0.,10000000000.};//TO REMOVE
  //Double_t *cutsDplus = new (Double_t*)fRDCuts->GetCuts();
  Int_t nSelectedloose=0,nSelectedtight=0;
  for (Int_t i3Prong = 0; i3Prong < n3Prong; i3Prong++) {
    AliAODRecoDecayHF3Prong *d = (AliAODRecoDecayHF3Prong*)array3Prong->UncheckedAt(i3Prong);
    fHistNEvents->Fill(4);
    if(fUseBit && !d->HasSelectionBit(AliRDHFCuts::kDplusCuts)){
      fHistNEvents->Fill(8);
      continue;
    }
    Bool_t unsetvtx=kFALSE;
    if(!d->GetOwnPrimaryVtx()){
      d->SetOwnPrimaryVtx(vtx1);
      unsetvtx=kTRUE;
    }

    if(fRDCutsProduction->IsSelected(d,AliRDHFCuts::kCandidate,aod)) {

      Double_t ptCand = d->Pt();
      Int_t iPtBin = fRDCutsProduction->PtBin(ptCand);
      
      Int_t passTightCuts=fRDCutsAnalysis->IsSelected(d,AliRDHFCuts::kCandidate,aod);
      Bool_t recVtx=kFALSE;
      AliAODVertex *origownvtx=0x0;
      if(fRDCutsProduction->GetIsPrimaryWithoutDaughters()){
	if(d->GetOwnPrimaryVtx()) origownvtx=new AliAODVertex(*d->GetOwnPrimaryVtx());	
	if(fRDCutsProduction->RecalcOwnPrimaryVtx(d,aod))recVtx=kTRUE;
	else fRDCutsProduction->CleanOwnPrimaryVtx(d,aod,origownvtx);
      }
      
      Int_t labDp=-1;
      Bool_t isPrimary=kTRUE;
      Float_t deltaPx=0.;
      Float_t deltaPy=0.;
      Float_t deltaPz=0.;
      Float_t truePt=0.;
      Float_t xDecay=0.;
      Float_t yDecay=0.;
      Float_t zDecay=0.;
      Float_t pdgCode=-2;
      Float_t trueImpParXY=0.;
      if(fReadMC){
	labDp = d->MatchToMC(411,arrayMC,3,pdgDgDplustoKpipi);
	if(labDp>=0){
	  AliAODMCParticle *partDp = (AliAODMCParticle*)arrayMC->At(labDp);
	  if(CheckOrigin(arrayMC,partDp)==5) isPrimary=kFALSE;
	  AliAODMCParticle *dg0 = (AliAODMCParticle*)arrayMC->At(partDp->GetDaughter(0));
	  deltaPx=partDp->Px()-d->Px();
	  deltaPy=partDp->Py()-d->Py();
	  deltaPz=partDp->Pz()-d->Pz();
	  truePt=partDp->Pt();
	  xDecay=dg0->Xv();	  
	  yDecay=dg0->Yv();	  
	  zDecay=dg0->Zv();
	  pdgCode=TMath::Abs(partDp->GetPdgCode());
	  if(!isPrimary){
	    trueImpParXY=GetTrueImpactParameter(mcHeader,arrayMC,partDp)*10000.;
	  }
	}else{
	  pdgCode=-1;
	}
      }

      Double_t invMass=d->InvMassDplus();
      Double_t rapid=d->YDplus();
      fYVsPt->Fill(ptCand,rapid);
      if(passTightCuts) {fYVsPtTC->Fill(ptCand,rapid);nOS++;}
      Bool_t isFidAcc=fRDCutsAnalysis->IsInFiducialAcceptance(ptCand,rapid);
      if(isFidAcc){
	fPtVsMass->Fill(invMass,ptCand);
	if(passTightCuts) fPtVsMassTC->Fill(invMass,ptCand);
      }
      Float_t tmp[24];
      Double_t  dlen=d->DecayLength();
      Double_t cosp=d->CosPointingAngle();
      Double_t sumD02=d->Getd0Prong(0)*d->Getd0Prong(0)+d->Getd0Prong(1)*d->Getd0Prong(1)+d->Getd0Prong(2)*d->Getd0Prong(2);
      Double_t dca=d->GetDCA();
      Double_t sigvert=d->GetSigmaVert();         
      Double_t ptmax=0;
      for(Int_t i=0;i<3;i++){
	if(d->PtProng(i)>ptmax)ptmax=d->PtProng(i);
      }
      Double_t impparXY=d->ImpParXY()*10000.;
      Double_t arrayForSparse[3]={invMass,ptCand,impparXY};
      Double_t arrayForSparseTrue[3]={invMass,ptCand,trueImpParXY};
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
	tmp[12]=cosp;
	tmp[13]=dlen;
	tmp[14]=d->Xv();
	tmp[15]=d->Yv();
	tmp[16]=d->Zv();
	tmp[17]=d->InvMassDplus();
	tmp[18]=sigvert;
	tmp[19]=d->Getd0Prong(0);
	tmp[20]=d->Getd0Prong(1);
	tmp[21]=d->Getd0Prong(2);
	tmp[22]=dca;
	tmp[23]=d->Prodd0d0(); 
	fNtupleDplus->Fill(tmp);
	PostData(4,fNtupleDplus);
      }
      if(iPtBin>=0){
	Float_t dlxy=d->NormalizedDecayLengthXY();
	Float_t cxy=d->CosPointingAngleXY();
	index=GetHistoIndex(iPtBin);
	if(isFidAcc){
	  fHistNEvents->Fill(5);
	  nSelectedloose++;
	  fMassHist[index]->Fill(invMass);
	  if(fCutsDistr){	  
	    fCosPHist[index]->Fill(cosp);
	    fDLenHist[index]->Fill(dlen);
	    fSumd02Hist[index]->Fill(sumD02);
	    fSigVertHist[index]->Fill(sigvert);
	    fPtMaxHist[index]->Fill(ptmax);
	    fPtKHist[index]->Fill(d->PtProng(1));
	    fPtpi1Hist[index]->Fill(d->PtProng(0));
	    fPtpi2Hist[index]->Fill(d->PtProng(2));
	    fDCAHist[index]->Fill(dca);
	    fDLxy[index]->Fill(dlxy);
	    fCosxy[index]->Fill(cxy);
	    fCorreld0Kd0pi[0]->Fill(d->Getd0Prong(0)*d->Getd0Prong(1),
				    d->Getd0Prong(2)*d->Getd0Prong(1));
	  }
	  if(passTightCuts){ fHistNEvents->Fill(6);
	    nSelectedtight++;
	    fMassHistTC[index]->Fill(invMass);
	    if(fCutsDistr){  
	      fDLxyTC[index]->Fill(dlxy);
	      fCosxyTC[index]->Fill(cxy);
	    }
	    if(d->GetCharge()>0) fMassHistTCPlus[index]->Fill(invMass);
	    else if(d->GetCharge()<0) fMassHistTCMinus[index]->Fill(invMass);
	    if(fDoImpPar){
	      fHistMassPtImpParTC[0]->Fill(arrayForSparse);
	    }
 	  }
	}
      
	if(fReadMC){
	  //  if(fCutsDistr){
	  if(labDp>=0) {
	    index=GetSignalHistoIndex(iPtBin);
	    if(isFidAcc){
	      Float_t factor[3]={1.,1.,1.};
	      if(fUseStrangeness){
		for(Int_t iprong=0;iprong<3;iprong++){
		  AliAODTrack *trad = (AliAODTrack*)d->GetDaughter(iprong);
		  Int_t labd= trad->GetLabel();
		  if(labd>=0){
		    AliAODMCParticle *dau = (AliAODMCParticle*)arrayMC->At(labd);
		    if(dau){
		      Int_t labm = dau->GetMother();
		      if(labm>=0){
			AliAODMCParticle *mot = (AliAODMCParticle*)arrayMC->At(labm);
			if(mot){
			  if(TMath::Abs(mot->GetPdgCode())==310 || TMath::Abs(mot->GetPdgCode())==130 || TMath::Abs(mot->GetPdgCode())==321){ //K0_S, K0_L, K^+-
			    if(d->PtProng(iprong)<=1)factor[iprong]=1./.7;
			    else factor[iprong]=1./.6;
			    //          fNentries->Fill(12);
			  }
			  if(TMath::Abs(mot->GetPdgCode())==3122) { //Lambda
			    factor[iprong]=1./0.25;
			    //		  fNentries->Fill(13);
			  }//if 3122
			}//if(mot)
		      }//if labm>0
		    }//if(dau)
		  }//if labd>=0
		}//prong loop
	      }
	      Float_t fact=1.;for(Int_t k=0;k<3;k++)fact=fact*factor[k];
	      fMassHist[index]->Fill(invMass);
	      if(fCutsDistr){
		fCosPHist[index]->Fill(cosp,fact);
		fDLenHist[index]->Fill(dlen,fact);
		fDLxy[index]->Fill(dlxy);
		fCosxy[index]->Fill(cxy);
	    
		Float_t sumd02s=d->Getd0Prong(0)*d->Getd0Prong(0)*factor[0]*factor[0]+d->Getd0Prong(1)*d->Getd0Prong(1)*factor[1]*factor[1]+d->Getd0Prong(2)*d->Getd0Prong(2)*factor[2]*factor[2];
		fSumd02Hist[index]->Fill(sumd02s);
		fSigVertHist[index]->Fill(sigvert,fact);
		fPtMaxHist[index]->Fill(ptmax,fact);
		fPtKHist[index]->Fill(d->PtProng(1),fact);
		fPtpi1Hist[index]->Fill(d->PtProng(0),fact);
		fPtpi2Hist[index]->Fill(d->PtProng(2),fact);
		fDCAHist[index]->Fill(dca,fact);
		fCorreld0Kd0pi[1]->Fill(d->Getd0Prong(0)*d->Getd0Prong(1),
					d->Getd0Prong(2)*d->Getd0Prong(1));
	      }
	      if(passTightCuts){
		fMassHistTC[index]->Fill(invMass);	      
		if(fCutsDistr){
		  fDLxyTC[index]->Fill(dlxy);
		  fCosxyTC[index]->Fill(cxy);
		}	      
		if(d->GetCharge()>0) fMassHistTCPlus[index]->Fill(invMass);
		else if(d->GetCharge()<0) fMassHistTCMinus[index]->Fill(invMass);
		if(fDoImpPar){
		  if(isPrimary) fHistMassPtImpParTC[1]->Fill(arrayForSparse);
		  else{
		    fHistMassPtImpParTC[2]->Fill(arrayForSparse);
		    fHistMassPtImpParTC[3]->Fill(arrayForSparseTrue);
		  }
		}
	      }
	    }	    
	    fYVsPtSig->Fill(ptCand,rapid);
	    if(passTightCuts) fYVsPtSigTC->Fill(ptCand,rapid);
	  }else{
	    index=GetBackgroundHistoIndex(iPtBin);
	    if(isFidAcc){
	      Float_t factor[3]={1.,1.,1.};
	      if(fUseStrangeness){
		for(Int_t iprong=0;iprong<3;iprong++){
		  AliAODTrack *trad = (AliAODTrack*)d->GetDaughter(iprong);
		  Int_t labd= trad->GetLabel();
		  if(labd>=0){
		    AliAODMCParticle *dau = (AliAODMCParticle*)arrayMC->At(labd);
		    if(dau){
		      Int_t labm = dau->GetMother();
		      if(labm>=0){
			AliAODMCParticle *mot = (AliAODMCParticle*)arrayMC->At(labm);
			if(mot){
			  if(TMath::Abs(mot->GetPdgCode())==310 || TMath::Abs(mot->GetPdgCode())==130 || TMath::Abs(mot->GetPdgCode())==321){ //K0_S, K0_L, K^+-
			    if(d->PtProng(iprong)<=1)factor[iprong]=1./.7;
			    else factor[iprong]=1./.6;
			    //          fNentries->Fill(12);
			  }
			  if(TMath::Abs(mot->GetPdgCode())==3122) { //Lambda
			    factor[iprong]=1./0.25;
			    //		  fNentries->Fill(13);
			  }//if 3122
			}//if(mot)
		      }//if labm>0
		    }//if(dau)
		  }//if labd>=0
		}//prong loop
	      }
	    
	      Float_t fact=1.;for(Int_t k=0;k<3;k++)fact=fact*factor[k];
	      fMassHist[index]->Fill(invMass);
	      if(fCutsDistr){
		fCosPHist[index]->Fill(cosp,fact);
		fDLenHist[index]->Fill(dlen,fact);
		fDLxy[index]->Fill(dlxy);
		fCosxy[index]->Fill(cxy);
	    
		Float_t sumd02s=d->Getd0Prong(0)*d->Getd0Prong(0)*factor[0]*factor[0]+d->Getd0Prong(1)*d->Getd0Prong(1)*factor[1]*factor[1]+d->Getd0Prong(2)*d->Getd0Prong(2)*factor[2]*factor[2];
		fSumd02Hist[index]->Fill(sumd02s);
		fSigVertHist[index]->Fill(sigvert,fact);
		fPtMaxHist[index]->Fill(ptmax,fact);
		fPtKHist[index]->Fill(d->PtProng(1),fact);
		fPtpi1Hist[index]->Fill(d->PtProng(0),fact);
		fPtpi2Hist[index]->Fill(d->PtProng(2),fact);
		fDCAHist[index]->Fill(dca,fact);
		fCorreld0Kd0pi[2]->Fill(d->Getd0Prong(0)*d->Getd0Prong(1),
					d->Getd0Prong(2)*d->Getd0Prong(1));
	      }
	      if(passTightCuts){
		fMassHistTC[index]->Fill(invMass);
		if(fCutsDistr){
		  fDLxyTC[index]->Fill(dlxy);
		  fCosxyTC[index]->Fill(cxy);
		}
		if(d->GetCharge()>0) fMassHistTCPlus[index]->Fill(invMass);
		else if(d->GetCharge()<0) fMassHistTCMinus[index]->Fill(invMass);
		if(fDoImpPar){
		  fHistMassPtImpParTC[4]->Fill(arrayForSparse);
		}
	      }
	    }
	  }
	  
	}
      }
    
      if(recVtx)fRDCutsProduction->CleanOwnPrimaryVtx(d,aod,origownvtx);
    }
    if(unsetvtx) d->UnsetOwnPrimaryVtx();
  }
  fCounter->StoreCandidates(aod,nSelectedloose,kTRUE);
  fCounter->StoreCandidates(aod,nSelectedtight,kFALSE);
  }
  //start LS analysis
  if(fDoLS && arrayLikeSign) LSAnalysis(array3Prong,arrayLikeSign,aod,vtx1,nOS);
  
  PostData(1,fOutput); 
  PostData(2,fListCuts);
  PostData(3,fCounter);    
  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEDplus::CreateLikeSignHistos(){
  // Histos for Like Sign bckground

  TString hisname;
  Int_t indexLS=0;
  Int_t index=0;
  Int_t nbins=GetNBinsHistos();
  for(Int_t i=0;i<fNPtBins;i++){

    index=GetHistoIndex(i);
    indexLS=GetLSHistoIndex(i);

    hisname.Form("hLSPt%dLC",i);
    fMassHistLS[indexLS] = new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHistLS[indexLS]->Sumw2();
    hisname.Form("hLSPt%dTC",i);
    fMassHistLSTC[indexLS] = new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHistLSTC[indexLS]->Sumw2();

    hisname.Form("hCosPAllPt%dLS",i);
    fCosPHistLS[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.5,1.);
    fCosPHistLS[index]->Sumw2();
    hisname.Form("hDLenAllPt%dLS",i);
    fDLenHistLS[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.,0.5);
    fDLenHistLS[index]->Sumw2();
    hisname.Form("hSumd02AllPt%dLS",i);
    fSumd02HistLS[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.,1.);
    fSumd02HistLS[index]->Sumw2();
    hisname.Form("hSigVertAllPt%dLS",i);
    fSigVertHistLS[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.,0.1);
    fSigVertHistLS[index]->Sumw2();
    hisname.Form("hPtMaxAllPt%dLS",i);
    fPtMaxHistLS[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.5,5.);
    fPtMaxHistLS[index]->Sumw2();
    hisname.Form("hDCAAllPt%dLS",i);
    fDCAHistLS[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.,0.1);
    fDCAHistLS[index]->Sumw2();    

    index=GetSignalHistoIndex(i);    
    indexLS++;
 
    hisname.Form("hLSPt%dLCnw",i);
    fMassHistLS[indexLS]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHistLS[indexLS]->Sumw2();
    hisname.Form("hLSPt%dTCnw",i);
    fMassHistLSTC[indexLS]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHistLSTC[indexLS]->Sumw2();

    hisname.Form("hCosPSigPt%dLS",i);
    fCosPHistLS[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.5,1.);
    fCosPHistLS[index]->Sumw2();
    hisname.Form("hDLenSigPt%dLS",i);
    fDLenHistLS[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.,0.5);
    fDLenHistLS[index]->Sumw2();
    hisname.Form("hSumd02SigPt%dLS",i);
    fSumd02HistLS[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.,1.);
    fSumd02HistLS[index]->Sumw2();
    hisname.Form("hSigVertSigPt%dLS",i);
    fSigVertHistLS[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.,0.1);
    fSigVertHistLS[index]->Sumw2();
    hisname.Form("hPtMaxSigPt%dLS",i);
    fPtMaxHistLS[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.5,5.);
    fPtMaxHistLS[index]->Sumw2();
    hisname.Form("hDCASigPt%dLS",i);
    fDCAHistLS[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.,0.1);
    fDCAHistLS[index]->Sumw2();

    index=GetBackgroundHistoIndex(i); 
    indexLS++;

    hisname.Form("hLSPt%dLCntrip",i);
    fMassHistLS[indexLS]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHistLS[indexLS]->Sumw2();
    hisname.Form("hLSPt%dTCntrip",i);
    fMassHistLSTC[indexLS]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHistLSTC[indexLS]->Sumw2();

    hisname.Form("hCosPBkgPt%dLS",i);
    fCosPHistLS[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.5,1.);
    fCosPHistLS[index]->Sumw2();
    hisname.Form("hDLenBkgPt%dLS",i);
    fDLenHistLS[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.,0.5);
    fDLenHistLS[index]->Sumw2();
    hisname.Form("hSumd02BkgPt%dLS",i);
    fSumd02HistLS[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.,1.);
    fSumd02HistLS[index]->Sumw2();
    hisname.Form("hSigVertBkgPt%dLS",i);
    fSigVertHistLS[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.,0.1);
    fSigVertHistLS[index]->Sumw2();
    hisname.Form("hPtMaxBkgPt%dLS",i);
    fPtMaxHistLS[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.5,5.);
    fPtMaxHistLS[index]->Sumw2();
    hisname.Form("hDCABkgPt%dLS",i);
    fDCAHistLS[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,0.,0.1);
    fDCAHistLS[index]->Sumw2();

    indexLS++;
    hisname.Form("hLSPt%dLCntripsinglecut",i);
    fMassHistLS[indexLS]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHistLS[indexLS]->Sumw2();
    hisname.Form("hLSPt%dTCntripsinglecut",i);
    fMassHistLSTC[indexLS]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHistLSTC[indexLS]->Sumw2();

    indexLS++;
    hisname.Form("hLSPt%dLCspc",i);
    fMassHistLS[indexLS]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHistLS[indexLS]->Sumw2();
    hisname.Form("hLSPt%dTCspc",i);
    fMassHistLSTC[indexLS]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHistLSTC[indexLS]->Sumw2();
  }

  for(Int_t i=0; i<3*fNPtBins; i++){
    fOutput->Add(fCosPHistLS[i]);
    fOutput->Add(fDLenHistLS[i]);
    fOutput->Add(fSumd02HistLS[i]);
    fOutput->Add(fSigVertHistLS[i]);
    fOutput->Add(fPtMaxHistLS[i]);  
    fOutput->Add(fDCAHistLS[i]);  
  }
  for(Int_t i=0; i<5*fNPtBins; i++){
    fOutput->Add(fMassHistLS[i]);
    fOutput->Add(fMassHistLSTC[i]);
  }
}

//________________________________________________________________________
void AliAnalysisTaskSEDplus::CreateImpactParameterHistos(){
  // Histos for impact paramter study

  Int_t nmassbins=GetNBinsHistos();
  Int_t nbins[3]={nmassbins,200,fNImpParBins};
  Double_t xmin[3]={fLowmasslimit,0.,fLowerImpPar};
  Double_t xmax[3]={fUpmasslimit,20.,fHigherImpPar};

  fHistMassPtImpParTC[0]=new THnSparseF("hMassPtImpParAll",
					"Mass vs. pt vs.imppar - All",
					3,nbins,xmin,xmax);
  fHistMassPtImpParTC[1]=new THnSparseF("hMassPtImpParPrompt",
					"Mass vs. pt vs.imppar - promptD",
					3,nbins,xmin,xmax);
  fHistMassPtImpParTC[2]=new THnSparseF("hMassPtImpParBfeed",
					"Mass vs. pt vs.imppar - DfromB",
					3,nbins,xmin,xmax);
  fHistMassPtImpParTC[3]=new THnSparseF("hMassPtImpParTrueBfeed",
					"Mass vs. pt vs.true imppar -DfromB",
					3,nbins,xmin,xmax);
  fHistMassPtImpParTC[4]=new THnSparseF("hMassPtImpParBkg",
				        "Mass vs. pt vs.imppar - backgr.",
					3,nbins,xmin,xmax);

  for(Int_t i=0; i<5;i++){
    fOutput->Add(fHistMassPtImpParTC[i]);
  }
}

//________________________________________________________________________
void AliAnalysisTaskSEDplus::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
  if(fDebug > 1) printf("AnalysisTaskSEDplus: Terminate() \n");

  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {     
    printf("ERROR: fOutput not available\n");
    return;
  }
  fHistNEvents = dynamic_cast<TH1F*>(fOutput->FindObject("fHistNEvents"));

  TString hisname;
  Int_t index=0;

  for(Int_t i=0;i<fNPtBins;i++){
    index=GetHistoIndex(i);
    
    hisname.Form("hMassPt%dTC",i);
    fMassHistTC[index]=dynamic_cast<TH1F*>(fOutput->FindObject(hisname.Data()));
  } 
    
  TCanvas *c1=new TCanvas("c1","D+ invariant mass distribution",500,500);
  c1->cd();
  TH1F *hMassPt=(TH1F*)fOutput->FindObject("hMassPt3TC");
  hMassPt->SetLineColor(kBlue);
  hMassPt->SetXTitle("M[GeV/c^{2}]"); 
  hMassPt->Draw();
 
  return;
}
//_________________________________________________________________________________________________
Int_t AliAnalysisTaskSEDplus::CheckOrigin(TClonesArray* arrayMC, AliAODMCParticle *mcPartCandidate) const {		
  //
  // checking whether the mother of the particles come from a charm or a bottom quark
  //
	
  Int_t pdgGranma = 0;
  Int_t mother = 0;
  mother = mcPartCandidate->GetMother();
  Int_t istep = 0;
  Int_t abspdgGranma =0;
  Bool_t isFromB=kFALSE;
  Bool_t isQuarkFound=kFALSE;
  while (mother >0 ){
    istep++;
    AliAODMCParticle* mcGranma = dynamic_cast<AliAODMCParticle*>(arrayMC->At(mother));
    if (mcGranma){
      pdgGranma = mcGranma->GetPdgCode();
      abspdgGranma = TMath::Abs(pdgGranma);
      if ((abspdgGranma > 500 && abspdgGranma < 600) || (abspdgGranma > 5000 && abspdgGranma < 6000)){
	isFromB=kTRUE;
      }
      if(abspdgGranma==4 || abspdgGranma==5) isQuarkFound=kTRUE;
      mother = mcGranma->GetMother();
    }else{
      AliError("Failed casting the mother particle!");
      break;
    }
  }
  
  if(isFromB) return 5;
  else return 4;
}
//_________________________________________________________________________________________________
Float_t AliAnalysisTaskSEDplus::GetTrueImpactParameter(AliAODMCHeader *mcHeader, TClonesArray* arrayMC, AliAODMCParticle *partDp) const {
  // true impact parameter calculation

  Double_t vtxTrue[3];
  mcHeader->GetVertex(vtxTrue);
  Double_t origD[3];
  partDp->XvYvZv(origD);
  Short_t charge=partDp->Charge();
  Double_t pXdauTrue[3],pYdauTrue[3],pZdauTrue[3];
  for(Int_t iDau=0; iDau<3; iDau++){
    pXdauTrue[iDau]=0.;
    pYdauTrue[iDau]=0.;
    pZdauTrue[iDau]=0.;
  }

  Int_t nDau=partDp->GetNDaughters();
  Int_t labelFirstDau = partDp->GetDaughter(0); 
  if(nDau==3){
    for(Int_t iDau=0; iDau<3; iDau++){
      Int_t ind = labelFirstDau+iDau;
      AliAODMCParticle* part = dynamic_cast<AliAODMCParticle*>(arrayMC->At(ind));
      if(!part){
	AliError("Daughter particle not found in MC array");
	return 99999.;
      } 
      pXdauTrue[iDau]=part->Px();
      pYdauTrue[iDau]=part->Py();
      pZdauTrue[iDau]=part->Pz();
    }
  }else if(nDau==2){
    Int_t theDau=0;
    for(Int_t iDau=0; iDau<2; iDau++){
      Int_t ind = labelFirstDau+iDau;
      AliAODMCParticle* part = dynamic_cast<AliAODMCParticle*>(arrayMC->At(ind));
      if(!part){
	AliError("Daughter particle not found in MC array");
	return 99999.;
      } 
      Int_t pdgCode=TMath::Abs(part->GetPdgCode());
      if(pdgCode==211 || pdgCode==321){
	pXdauTrue[theDau]=part->Px();
	pYdauTrue[theDau]=part->Py();
	pZdauTrue[theDau]=part->Pz();
	++theDau;
      }else{
	Int_t nDauRes=part->GetNDaughters();
	if(nDauRes==2){
	  Int_t labelFirstDauRes = part->GetDaughter(0); 	
	  for(Int_t iDauRes=0; iDauRes<2; iDauRes++){
	    Int_t indDR = labelFirstDauRes+iDauRes;
	    AliAODMCParticle* partDR = dynamic_cast<AliAODMCParticle*>(arrayMC->At(indDR));
	    if(!partDR){
	      AliError("Daughter particle not found in MC array");
	      return 99999.;
	    } 
	    
	    Int_t pdgCodeDR=TMath::Abs(partDR->GetPdgCode());
 	    if(pdgCodeDR==211 || pdgCodeDR==321){
	      pXdauTrue[theDau]=partDR->Px();
	      pYdauTrue[theDau]=partDR->Py();
	      pZdauTrue[theDau]=partDR->Pz();
	      ++theDau;
	    }
	  }
	}
      }
    }
    if(theDau!=3){
      AliError("Wrong number of decay prongs");
      return 99999.;
    }
  }

  Double_t d0dummy[3]={0.,0.,0.};
  AliAODRecoDecayHF aodDplusMC(vtxTrue,origD,3,charge,pXdauTrue,pYdauTrue,pZdauTrue,d0dummy);
  return aodDplusMC.ImpParXY();

}
