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
// comparison of heavy-flavour decay candidates
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
#include "AliVertexingHFUtils.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskSEDplus);
/// \endcond

//________________________________________________________________________
AliAnalysisTaskSEDplus::AliAnalysisTaskSEDplus():
  AliAnalysisTaskSE(),
  fOutput(0),
  fHistNEvents(0),
  fHistNCandidates(0),
  fMassHist(0x0),
  fMassHistPlus(0x0),
  fMassHistMinus(0x0),
  fMassHistNoPid(0x0),
  fCosPHist(0x0),
  fDLenHist(0x0),
  fSumd02Hist(0x0),
  fSigVertHist(0x0),
  fPtMaxHist(0x0),
  fPtKHist(0x0),
  fPtpi1Hist(0x0),
  fPtpi2Hist(0x0),
  fDCAHist(0x0),
  fDLxy(0x0),
  fCosxy(0x0),
  fHistTrackVar(0),
  fMCAccPrompt(0),
  fMCAccBFeed(0),
  fPtVsMassNoPid(0),
  fPtVsMass(0),
  fPtVsMassPlus(0),
  fPtVsMassMinus(0),
  fPtVsMassBadDaus(0),
  fPtVsMassGoodDaus(0),
  fYVsPtNoPid(0),
  fYVsPt(0),
  fYVsPtSigNoPid(0),
  fYVsPtSig(0),
  fPhiEtaCand(0),
  fPhiEtaCandSigReg(0),
  fSPDMult(0),
  fDaughterClass(0),
  fDeltaID(0),
  fIDDauVsIDTra(0),
  fMassHistLS(0x0),
  fCosPHistLS(0x0),
  fDLenHistLS(0x0),
  fSumd02HistLS(0x0),
  fSigVertHistLS(0x0),
  fPtMaxHistLS(0x0),
  fDCAHistLS(0x0),
  fNtupleDplus(0),
  fUpmasslimit(1.965),
  fLowmasslimit(1.765),
  fNPtBins(1),
  fBinWidth(0.002),
  fListCuts(0),
  fRDCutsAnalysis(0),
  fCounter(0),
  fFillNtuple(0),
  fAODProtection(1),
  fReadMC(kFALSE),
  fUseStrangeness(kFALSE),
  fUseBit(kTRUE),
  fCutsDistr(kFALSE),
  fDoImpPar(kFALSE),
  fDoSparse(kFALSE),
  fDoTrackVarHist(kFALSE),
  fStepMCAcc(kFALSE),
  fUseQuarkTagInKine(kTRUE),
  fNImpParBins(400),
  fLowerImpPar(-1000.),
  fHigherImpPar(1000.),
  fDoLS(0),
  fEtaSelection(0),
  fSystem(0),
  fNtrcklMin(0),
  fNtrcklMax(10000),
  fCutOnTrckl(kFALSE),
  fFillOnlySignalSparses(kFALSE),
  fUseFinPtBinsForSparse(kFALSE)
{
  /// Default constructor

  for(Int_t i=0;i<3;i++){
    fHistCentrality[i]=0;
    fCorreld0Kd0pi[i]=0;
  }

  for(Int_t i=0; i<5; i++)fHistMassPtImpPar[i]=0;
  for(Int_t i=0; i<3; i++)fSparseCutVars[i]=0;
}

//________________________________________________________________________
AliAnalysisTaskSEDplus::AliAnalysisTaskSEDplus(const char *name,AliRDHFCutsDplustoKpipi *dpluscutsana,Int_t fillNtuple):
  AliAnalysisTaskSE(name),
  fOutput(0),
  fHistNEvents(0),
  fHistNCandidates(0),
  fMassHist(0x0),
  fMassHistPlus(0x0),
  fMassHistMinus(0x0),
  fMassHistNoPid(0x0),
  fCosPHist(0x0),
  fDLenHist(0x0),
  fSumd02Hist(0x0),
  fSigVertHist(0x0),
  fPtMaxHist(0x0),
  fPtKHist(0x0),
  fPtpi1Hist(0x0),
  fPtpi2Hist(0x0),
  fDCAHist(0x0),
  fDLxy(0x0),
  fCosxy(0x0),
  fHistTrackVar(0),
  fMCAccPrompt(0),
  fMCAccBFeed(0),
  fPtVsMassNoPid(0),
  fPtVsMass(0),
  fPtVsMassPlus(0),
  fPtVsMassMinus(0),
  fPtVsMassBadDaus(0),
  fPtVsMassGoodDaus(0),
  fYVsPtNoPid(0),
  fYVsPt(0),
  fYVsPtSigNoPid(0),
  fYVsPtSig(0),
  fPhiEtaCand(0),
  fPhiEtaCandSigReg(0),
  fSPDMult(0),
  fDaughterClass(0),
  fDeltaID(0),
  fIDDauVsIDTra(0),
  fMassHistLS(0x0),
  fCosPHistLS(0x0),
  fDLenHistLS(0x0),
  fSumd02HistLS(0x0),
  fSigVertHistLS(0x0),
  fPtMaxHistLS(0x0),
  fDCAHistLS(0x0),
  fNtupleDplus(0),
  fUpmasslimit(1.965),
  fLowmasslimit(1.765),
  fNPtBins(1),
  fBinWidth(0.002),
  fListCuts(0),
  fRDCutsAnalysis(dpluscutsana),
  fCounter(0),
  fFillNtuple(fillNtuple),
  fAODProtection(1),
  fReadMC(kFALSE),
  fUseStrangeness(kFALSE),
  fUseBit(kTRUE),
  fCutsDistr(kFALSE),
  fDoImpPar(kFALSE),
  fDoSparse(kFALSE),
  fDoTrackVarHist(kFALSE),
  fStepMCAcc(kFALSE),
  fUseQuarkTagInKine(kTRUE),
  fNImpParBins(400),
  fLowerImpPar(-1000.),
  fHigherImpPar(1000.),
  fDoLS(0),
  fEtaSelection(0),
  fSystem(0),
  fNtrcklMin(0),
  fNtrcklMax(10000),
  fCutOnTrckl(kFALSE),
  fFillOnlySignalSparses(kFALSE),
  fUseFinPtBinsForSparse(kFALSE)
{
  //
  /// Standrd constructor
  //
  fNPtBins=fRDCutsAnalysis->GetNPtBins();

  for(Int_t i=0;i<3;i++){
    fHistCentrality[i]=0;
    fCorreld0Kd0pi[i]=0;
  }

  for(Int_t i=0; i<5; i++)fHistMassPtImpPar[i]=0;
  for(Int_t i=0; i<3; i++)fSparseCutVars[i]=0;

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
  /// Destructor
  //
  if(fOutput && !fOutput->IsOwner()){
    delete fHistNEvents;
    delete fHistNCandidates;
    for(Int_t i=0;i<3*fNPtBins;i++){
      delete fMassHist[i];
      delete fMassHistPlus[i];
      delete fMassHistMinus[i];
      delete fMassHistNoPid[i];
      delete fCosPHist[i];
      delete fDLenHist[i];
      delete fSumd02Hist[i];
      delete fSigVertHist[i];
      delete fPtMaxHist[i];
      delete fPtKHist[i];
      delete fPtpi1Hist[i];
      delete fPtpi2Hist[i];
      delete fDCAHist[i];
      delete fDLxy[i];
      delete fCosxy[i];
      delete fCosPHistLS[i];
      delete fDLenHistLS[i];
      delete fSumd02HistLS[i];
      delete fSigVertHistLS[i];
      delete fPtMaxHistLS[i];
      delete fDCAHistLS[i];
    }
    delete [] fMassHist;
    delete [] fMassHistPlus;
    delete [] fMassHistMinus;
    delete [] fMassHistNoPid;
    delete [] fCosPHist;
    delete [] fDLenHist;
    delete [] fSumd02Hist;
    delete [] fSigVertHist;
    delete [] fPtMaxHist;
    delete [] fPtKHist;
    delete [] fPtpi1Hist;
    delete [] fPtpi2Hist;
    delete [] fDCAHist;
    delete [] fDLxy;
    delete [] fCosxy;

    for(Int_t i=0;i<5*fNPtBins;i++){
      delete fMassHistLS[i];
    }

    delete [] fMassHistLS;
    delete [] fCosPHistLS;
    delete [] fDLenHistLS;
    delete [] fSumd02HistLS;
    delete [] fSigVertHistLS;
    delete [] fPtMaxHistLS;
    delete [] fDCAHistLS;

    for(Int_t i=0;i<3;i++){
      delete fCorreld0Kd0pi[i];
      delete fHistCentrality[i];
    }
    for(Int_t i=0;i<5;i++){
      delete fHistMassPtImpPar[i];
    }
    for(Int_t i=0;i<3;i++){
      delete fSparseCutVars[i];
    }
    delete fPtVsMassNoPid;
    delete fPtVsMass;
    delete fPtVsMassPlus;
    delete fPtVsMassMinus;
    delete fPtVsMassBadDaus;
    delete fPtVsMassGoodDaus;
    delete fYVsPtNoPid;
    delete fYVsPt;
    delete fYVsPtSigNoPid;
    delete fYVsPtSig;
    delete fPhiEtaCand;
    delete fPhiEtaCandSigReg;
    delete fSPDMult;
    delete fDaughterClass;
    delete fDeltaID;
    delete fIDDauVsIDTra;
    delete fHistTrackVar;
    delete fMCAccPrompt;
    delete fMCAccBFeed;
  }

  delete fOutput;
  delete fNtupleDplus;
  delete fListCuts;
  delete fRDCutsAnalysis;
  delete fCounter;

}
//_________________________________________________________________
void  AliAnalysisTaskSEDplus::SetMassLimits(Float_t range){
  /// set invariant mass limits
  Float_t bw=GetBinWidth();
  fUpmasslimit = 1.865+range;
  fLowmasslimit = 1.865-range;
  SetBinWidth(bw);
}
//_________________________________________________________________
void  AliAnalysisTaskSEDplus::SetMassLimits(Float_t lowlimit, Float_t uplimit){
  /// set invariant mass limits
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
  /// Fill the Like Sign histograms
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
    fRDCutsAnalysis->IsSelected(d,AliRDHFCuts::kAll,aod);
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
  /// Initialization
  //
  if(fDebug > 1) printf("AnalysisTaskSEDplus::Init() \n");

  //PostData(2,fRDCutsloose);//we should then put those cuts in a tlist if we have more than 1
  fListCuts=new TList();
  fListCuts->SetOwner();
  AliRDHFCutsDplustoKpipi *analysis = new AliRDHFCutsDplustoKpipi(*fRDCutsAnalysis);
  analysis->SetName("AnalysisCuts");

  fListCuts->Add(analysis);
  PostData(2,fListCuts);

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEDplus::UserCreateOutputObjects()
{
  /// Create the output container
  //
  if(fDebug > 1) printf("AnalysisTaskSEDplus::UserCreateOutputObjects() \n");

  // Several histograms are more conveniently managed in a TList
  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("OutputHistos");

  fHistNEvents = new TH1F("fHistNEvents", "number of events ",10,-0.5,9.5);
  fHistNEvents->GetXaxis()->SetBinLabel(1,"nEvents read");
  fHistNEvents->GetXaxis()->SetBinLabel(2,"Rejected due to mismatch in trees");
  fHistNEvents->GetXaxis()->SetBinLabel(3,"nEvents with good AOD");
  fHistNEvents->GetXaxis()->SetBinLabel(4,"Rejected due to trigger");
  fHistNEvents->GetXaxis()->SetBinLabel(5,"Rejected due to vertex reco");
  fHistNEvents->GetXaxis()->SetBinLabel(6,"Rejected due to pileup");
  fHistNEvents->GetXaxis()->SetBinLabel(7,"Rejected due to centrality");
  fHistNEvents->GetXaxis()->SetBinLabel(8,"Rejected due to vtxz");
  fHistNEvents->GetXaxis()->SetBinLabel(9,"Rejected due to Physics Sel");
  fHistNEvents->GetXaxis()->SetBinLabel(10,"nEvents accepted");
  fHistNEvents->GetXaxis()->SetNdivisions(1,kFALSE);
  fHistNEvents->SetMinimum(0);
  fOutput->Add(fHistNEvents);

  fHistNCandidates = new TH1F("hNCandidates","number of candidates",7,-0.5,6.5);
  fHistNCandidates->GetXaxis()->SetBinLabel(1,"no. of 3prong candidates");
  fHistNCandidates->GetXaxis()->SetBinLabel(2,"no. of cand with D+ bitmask");
  fHistNCandidates->GetXaxis()->SetBinLabel(3,"D+ not on-the-fly reco");
  fHistNCandidates->GetXaxis()->SetBinLabel(4,"D+ after topological cuts");
  fHistNCandidates->GetXaxis()->SetBinLabel(5,"D+ after Topological+SingleTrack cuts");
  fHistNCandidates->GetXaxis()->SetBinLabel(6,"D+ after Topological+SingleTrack+PID cuts");
  fHistNCandidates->GetXaxis()->SetBinLabel(7,"D+ rejected by preselect");
  fHistNCandidates->GetXaxis()->SetNdivisions(1,kFALSE);
  fHistNCandidates->SetMinimum(0);
  fOutput->Add(fHistNCandidates);

  fMassHist = new TH1F*[3*fNPtBins];
  fMassHistPlus = new TH1F*[3*fNPtBins];
  fMassHistMinus  = new TH1F*[3*fNPtBins];
  fMassHistNoPid = new TH1F*[3*fNPtBins];
  fCosPHist = new TH1F*[3*fNPtBins];
  fDLenHist = new TH1F*[3*fNPtBins];
  fSumd02Hist = new TH1F*[3*fNPtBins];
  fSigVertHist = new TH1F*[3*fNPtBins];
  fPtMaxHist = new TH1F*[3*fNPtBins];
  fPtKHist = new TH1F*[3*fNPtBins];
  fPtpi1Hist = new TH1F*[3*fNPtBins];
  fPtpi2Hist = new TH1F*[3*fNPtBins];
  fDCAHist = new TH1F*[3*fNPtBins];
  fDLxy = new TH1F*[3*fNPtBins];
  fCosxy = new TH1F*[3*fNPtBins];

  TString hisname;
  Int_t index=0;
  Int_t nbins=GetNBinsHistos();
  fHistCentrality[0]=new TH2F("hCentrMult","centrality",100,0.5,30000.5,40,0.,100.);
  fHistCentrality[1]=new TH2F("hCentrMult(selectedCent)","centrality(selectedCent)",100,0.5,30000.5,40,0.,100.);
  fHistCentrality[2]=new TH2F("hCentrMult(OutofCent)","centrality(OutofCent)",100,0.5,30000.5,40,0.,100.);
  for(Int_t i=0;i<3;i++){
    fHistCentrality[i]->Sumw2();
    fOutput->Add(fHistCentrality[i]);
  }
  for(Int_t i=0;i<fNPtBins;i++){

    index=GetHistoIndex(i);

    hisname.Form("hMassNoPidPt%d",i);
    fMassHistNoPid[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHistNoPid[index]->Sumw2();
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
    hisname.Form("hPtKPt%d",i);
    fPtKHist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.);
    fPtKHist[index]->Sumw2();
    hisname.Form("hPtpi1Pt%d",i);
    fPtpi1Hist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.);
    fPtpi1Hist[index]->Sumw2();
    hisname.Form("hPtpi2Pt%d",i);
    fPtpi2Hist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.);
    fPtpi2Hist[index]->Sumw2();
    hisname.Form("hDCAAllPt%d",i);
    fDCAHist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.,0.1);
    fDCAHist[index]->Sumw2();

    hisname.Form("hDLxyPt%d",i);
    fDLxy[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.,10.);
    fDLxy[index]->Sumw2();
    hisname.Form("hCosxyPt%d",i);
    fCosxy[index]=new TH1F(hisname.Data(),hisname.Data(),100,-1,1.);
    fCosxy[index]->Sumw2();

    hisname.Form("hMassPt%dTC",i);
    fMassHist[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHist[index]->Sumw2();
    hisname.Form("hMassPt%dTCPlus",i);
    fMassHistPlus[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHistPlus[index]->Sumw2();
    hisname.Form("hMassPt%dTCMinus",i);
    fMassHistMinus[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHistMinus[index]->Sumw2();

    index=GetSignalHistoIndex(i);
    hisname.Form("hSigNoPidPt%d",i);
    fMassHistNoPid[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHistNoPid[index]->Sumw2();
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
    hisname.Form("hPtKSigPt%d",i);
    fPtKHist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.);
    fPtKHist[index]->Sumw2();
    hisname.Form("hPtpi1SigPt%d",i);
    fPtpi1Hist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.);
    fPtpi1Hist[index]->Sumw2();
    hisname.Form("hPtpi2SigPt%d",i);
    fPtpi2Hist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.);
    fPtpi2Hist[index]->Sumw2();

    hisname.Form("hDCASigPt%d",i);
    fDCAHist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.,0.1);
    fDCAHist[index]->Sumw2();

    hisname.Form("hDLxySigPt%d",i);
    fDLxy[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.,10.);
    fDLxy[index]->Sumw2();
    hisname.Form("hCosxySigPt%d",i);
    fCosxy[index]=new TH1F(hisname.Data(),hisname.Data(),100,-1,1.);
    fCosxy[index]->Sumw2();
    hisname.Form("hSigPt%dTC",i);
    fMassHist[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHist[index]->Sumw2();
    hisname.Form("hSigPt%dTCPlus",i);
    fMassHistPlus[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHistPlus[index]->Sumw2();
    hisname.Form("hSigPt%dTCMinus",i);
    fMassHistMinus[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHistMinus[index]->Sumw2();


    index=GetBackgroundHistoIndex(i);
    hisname.Form("hBkgNoPidPt%d",i);
    fMassHistNoPid[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHistNoPid[index]->Sumw2();
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
    hisname.Form("hPtKBkgPt%d",i);
    fPtKHist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.);
    fPtKHist[index]->Sumw2();
    hisname.Form("hPtpi1BkgPt%d",i);
    fPtpi1Hist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.);
    fPtpi1Hist[index]->Sumw2();
    hisname.Form("hPtpi2BkgPt%d",i);
    fPtpi2Hist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.,5.);
    fPtpi2Hist[index]->Sumw2();
    hisname.Form("hDCABkgPt%d",i);
    fDCAHist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.,0.1);
    fDCAHist[index]->Sumw2();

    hisname.Form("hDLxyBkgPt%d",i);
    fDLxy[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.,10.);
    fDLxy[index]->Sumw2();
    hisname.Form("hCosxyBkgPt%d",i);
    fCosxy[index]=new TH1F(hisname.Data(),hisname.Data(),100,-1,1.);
    fCosxy[index]->Sumw2();


    hisname.Form("hBkgPt%dTC",i);
    fMassHist[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHist[index]->Sumw2();
    hisname.Form("hBkgPt%dTCPlus",i);
    fMassHistPlus[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHistPlus[index]->Sumw2();
    hisname.Form("hBkgPt%dTCMinus",i);
    fMassHistMinus[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHistMinus[index]->Sumw2();
  }


  for(Int_t i=0; i<3*fNPtBins; i++){
    fOutput->Add(fMassHistNoPid[i]);
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
      fOutput->Add(fCosxy[i]);
    }
    fOutput->Add(fMassHist[i]);
    fOutput->Add(fMassHistPlus[i]);
    fOutput->Add(fMassHistMinus[i]);

  }

  fCorreld0Kd0pi[0]=new TH2F("hCorreld0Kd0piAll","",100,-0.02,0.02,100,-0.02,0.02);
  fCorreld0Kd0pi[1]=new TH2F("hCorreld0Kd0piSig","",100,-0.02,0.02,100,-0.02,0.02);
  fCorreld0Kd0pi[2]=new TH2F("hCorreld0Kd0piBkg","",100,-0.02,0.02,100,-0.02,0.02);
  if(fCutsDistr){
    for(Int_t i=0; i<3; i++){
      fCorreld0Kd0pi[i]->Sumw2();
      fOutput->Add(fCorreld0Kd0pi[i]);
    }
  }


  const Int_t nPtBins=440;
  Double_t ptBinLims[nPtBins+1];
  for(Int_t jb=0; jb<=300; jb++) ptBinLims[jb]=0.1*(Double_t)jb; // 100 MeV bins in 0<pt<30
  for(Int_t jb=301; jb<=440; jb++) ptBinLims[jb]=ptBinLims[300]+0.5*(Double_t)(jb-300); // 500 MeV bins in 30<pt<100

  fPtVsMassNoPid=new TH2F("hPtVsMassNoPid","PtVsMass (no PID)",nbins,fLowmasslimit,fUpmasslimit,nPtBins,ptBinLims);
  fPtVsMass=new TH2F("hPtVsMass","PtVsMass",nbins,fLowmasslimit,fUpmasslimit,nPtBins,ptBinLims);
  fPtVsMassPlus=new TH2F("hPtVsMassPlus","PtVsMass",nbins,fLowmasslimit,fUpmasslimit,nPtBins,ptBinLims);
  fPtVsMassMinus=new TH2F("hPtVsMassMinus","PtVsMass",nbins,fLowmasslimit,fUpmasslimit,nPtBins,ptBinLims);
  fPtVsMassGoodDaus=new TH2F("hPtVsMassGoodDaus","PtVsMassGoodDaus",nbins,fLowmasslimit,fUpmasslimit,nPtBins,ptBinLims);
  fPtVsMassBadDaus=new TH2F("hPtVsMassBadDaus","PtVsMassBadDaus",nbins,fLowmasslimit,fUpmasslimit,nPtBins,ptBinLims);

  fYVsPtNoPid=new TH3F("hYVsPtNoPid","YvsPt (no PID)",40,0.,20.,80,-2.,2.,nbins,fLowmasslimit,fUpmasslimit);
  fYVsPt=new TH3F("hYVsPt","YvsPt",40,0.,20.,80,-2.,2.,nbins,fLowmasslimit,fUpmasslimit);
  fYVsPtSigNoPid=new TH2F("hYVsPtSigNoPid","YvsPt (MC, only sig., no PID)",40,0.,20.,80,-2.,2.);
  fYVsPtSig=new TH2F("hYVsPtSig","YvsPt (MC, only Sig)",40,0.,20.,80,-2.,2.);
  fPhiEtaCand=new TH2F("hPhiEtaCand","phi vs. eta candidates",20,-1.,1.,50,0.,2*TMath::Pi());
  fPhiEtaCandSigReg=new TH2F("hPhiEtaCandSigReg","phi vs. eta candidates",20,-1.,1.,50,0.,2*TMath::Pi());

  Double_t maxmult;
  if(fSystem==1) maxmult=5000;
  else maxmult=200;
  fSPDMult = new TH1F("hSPDMult", "Tracklets multiplicity; Tracklets ; Entries",200,0.,maxmult);
  fOutput->Add(fPtVsMassNoPid);
  fOutput->Add(fPtVsMass);
  fOutput->Add(fPtVsMassPlus);
  fOutput->Add(fPtVsMassMinus);
  fOutput->Add(fPtVsMassGoodDaus);
  fOutput->Add(fPtVsMassBadDaus);
  fOutput->Add(fYVsPtNoPid);
  fOutput->Add(fYVsPt);
  fOutput->Add(fYVsPtSigNoPid);
  fOutput->Add(fYVsPtSig);
  fOutput->Add(fPhiEtaCand);
  fOutput->Add(fPhiEtaCandSigReg);
  fOutput->Add(fSPDMult);

  fDaughterClass = new TH1F("hDaughterClass","",10,-0.5,9.5);
  fDaughterClass->GetXaxis()->SetBinLabel(1,"AliAODTrack - good ID");
  fDaughterClass->GetXaxis()->SetBinLabel(2,"AliAODTrack - charge0");
  fDaughterClass->GetXaxis()->SetBinLabel(3,"AliAODTrack - ID=0");
  fDaughterClass->GetXaxis()->SetBinLabel(4,"AliAODTrack - neg ID");
  fDaughterClass->GetXaxis()->SetBinLabel(5,"AliAODTrack - different ID");
  fDaughterClass->GetXaxis()->SetBinLabel(6,"AliAODRecoDecayHF2Prong");
  fDaughterClass->GetXaxis()->SetBinLabel(7,"AliAODRecoDecayHF3Prong");
  fDaughterClass->GetXaxis()->SetBinLabel(8,"AliAODRecoCascadeHF");
  fDaughterClass->GetXaxis()->SetBinLabel(9,"Other");
  fDeltaID = new TH1F("hDeltaID"," ; GetDaughter->GetID() - GetProngID()",20001,-10000.5,10000.5);
  if(fSystem==0){
    fIDDauVsIDTra = new TH2F("hIDDauVsIDTra"," ; GetProngID() ; GetDaughter->GetID()",1001,-500.5,500.5,1001,-500.5,500.5);
  }else{
    fIDDauVsIDTra = new TH2F("hIDDauVsIDTra"," ; GetProngID() ; GetDaughter->GetID()",1000,-30000,30000,1000,-30000,30000);
  }
  fOutput->Add(fDaughterClass);
  fOutput->Add(fDeltaID);
  fOutput->Add(fIDDauVsIDTra);

  //Counter for Normalization
  TString normName="NormalizationCounter";
  AliAnalysisDataContainer *cont = GetOutputSlot(3)->GetContainer();
  if(cont)normName=(TString)cont->GetName();
  fCounter = new AliNormalizationCounter(normName.Data());
  fCounter->Init();

  if(fDoLS) CreateLikeSignHistos();
  if(fDoImpPar) CreateImpactParameterHistos();
  if(fDoSparse) CreateCutVarsSparses();
  if(fReadMC && fStepMCAcc) CreateMCAcceptanceHistos();
  if(fDoTrackVarHist) CreateTrackVarHistos();

  PostData(1,fOutput);

  if(fFillNtuple==1){
    OpenFile(4); // 4 is the slot number of the ntuple


    fNtupleDplus = new TNtuple("fNtupleDplus","D +","pdg:Px:Py:Pz:Pt:charge:piddau0:piddau1:piddau2:Ptpi:PtK:Ptpi2:mompi:momK:mompi2:cosp:cospxy:DecLen:NormDecLen:DecLenXY:NormDecLenXY:InvMass:sigvert:d0Pi:d0K:d0Pi2:maxdca:ntracks:centr:RunNumber:BadDau");

  }

  if(fFillNtuple==2){
    OpenFile(4); // 4 is the slot number of the ntuple


    fNtupleDplus = new TNtuple("fNtupleDplus","D +","pdg:Pt:InvMass:d0:origin");

  }

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEDplus::UserExec(Option_t */*option*/)
{
  /// Execute analysis for current event:
  /// heavy flavor candidates association to MC truth

  AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());

  fHistNEvents->Fill(0); // count event

  if(fAODProtection>=0){
    //   Protection against different number of events in the AOD and deltaAOD
    //   In case of discrepancy the event is rejected.
    Int_t matchingAODdeltaAODlevel = AliRDHFCuts::CheckMatchingAODdeltaAODevents();
    if (matchingAODdeltaAODlevel<0 || (matchingAODdeltaAODlevel==0 && fAODProtection==1)) {
      // AOD/deltaAOD trees have different number of entries || TProcessID do not match while it was required
      fHistNEvents->Fill(1);
      return;
    }
  }

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

  //Store the event in AliNormalizationCounter->To disable for Pb-Pb? Add a flag to disable it?
  fCounter->StoreEvent(aod,fRDCutsAnalysis,fReadMC);

  fHistNEvents->Fill(2);

  Int_t runNumber=aod->GetRunNumber();

  //Event selection
  Bool_t isEvSel=fRDCutsAnalysis->IsEventSelected(aod);
  Float_t ntracks=aod->GetNumberOfTracks();
  Int_t tracklets=AliVertexingHFUtils::GetNumberOfTrackletsInEtaRange(aod,-1.,1.);
  Float_t evCentr=0;
  if(fRDCutsAnalysis->GetUseCentrality()>0) evCentr=fRDCutsAnalysis->GetCentrality(aod);
  fHistCentrality[0]->Fill(ntracks,evCentr);
  if(fRDCutsAnalysis->GetWhyRejection()==5) fHistNEvents->Fill(3);
  if(!isEvSel && fRDCutsAnalysis->GetWhyRejection()==0) fHistNEvents->Fill(4);
  if(fRDCutsAnalysis->GetWhyRejection()==1) fHistNEvents->Fill(5);
  if(fRDCutsAnalysis->GetWhyRejection()==2){fHistNEvents->Fill(6);fHistCentrality[2]->Fill(ntracks,evCentr);}
  if(fRDCutsAnalysis->GetWhyRejection()==6)fHistNEvents->Fill(7);
  if(fRDCutsAnalysis->GetWhyRejection()==7)fHistNEvents->Fill(8);

  TClonesArray *arrayMC=0;
  AliAODMCHeader *mcHeader=0;
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
  if(fReadMC && fStepMCAcc){
    if(aod->GetTriggerMask()==0 &&
       (runNumber>=195344 && runNumber<=195677)){
      // protection for events with empty trigger mask in p-Pb
      return;
    }
    if(fRDCutsAnalysis->GetUseCentrality()>0 && fRDCutsAnalysis->IsEventSelectedInCentrality(aod)!=0) return;
    // events not passing the centrality selection can be removed immediately.

    FillMCAcceptanceHistos(arrayMC, mcHeader, tracklets);
  }
  // trigger class for PbPb C0SMH-B-NOPF-ALLNOTRD
  //TString trigclass=aod->GetFiredTriggerClasses();
  // Post the data already here
  PostData(1,fOutput);
  if(!isEvSel)return;
  // printf("ntracklet===%d\n",tracklets);
  if(fCutOnTrckl && (tracklets<fNtrcklMin || tracklets>fNtrcklMax)) {return;}
  fSPDMult->Fill(tracklets);

  fHistCentrality[1]->Fill(ntracks,evCentr);
  fHistNEvents->Fill(9);

  // AOD primary vertex
  AliAODVertex *vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();
  //    vtx1->Print();
  //   TString primTitle = vtx1->GetTitle();
  //if(primTitle.Contains("VertexerTracks") && vtx1->GetNContributors()>0)fHistNEvents->Fill(2);


  Int_t n3Prong = array3Prong->GetEntriesFast();
  //  printf("Number of D+->Kpipi: %d and of tracks: %d\n",n3Prong,aod->GetNumberOfTracks());

  Int_t nOS=0;
  Int_t index;
  Int_t pdgDgDplustoKpipi[3]={321,211,211};

  // vHF object is needed to call the method that refills the missing info of the candidates
  // if they have been deleted in dAOD reconstruction phase
  // in order to reduce the size of the file
  AliAnalysisVertexingHF *vHF=new AliAnalysisVertexingHF();
  if(fDoLS>1){//Normalizations for LS
    for (Int_t i3Prong = 0; i3Prong < n3Prong; i3Prong++) {
      AliAODRecoDecayHF3Prong *d = (AliAODRecoDecayHF3Prong*)array3Prong->UncheckedAt(i3Prong);
      if(fUseBit && !d->HasSelectionBit(AliRDHFCuts::kDplusCuts)){
	if(fRDCutsAnalysis->IsSelected(d,AliRDHFCuts::kAll,aod))nOS++;
      }
    }
  }else{//Standard analysis
    // Double_t cutsDplus[12]={0.2,0.4,0.4,0.,0.,0.01,0.06,0.02,0.,0.85,0.,10000000000.};//TO REMOVE
    //Double_t *cutsDplus = new (Double_t*)fRDCuts->GetCuts();
    Int_t nSelectednopid=0,nSelected=0;
    for (Int_t i3Prong = 0; i3Prong < n3Prong; i3Prong++) {
      AliAODRecoDecayHF3Prong *d = (AliAODRecoDecayHF3Prong*)array3Prong->UncheckedAt(i3Prong);
      fHistNCandidates->Fill(0);
      if(fUseBit && !d->HasSelectionBit(AliRDHFCuts::kDplusCuts)){
	continue;
      }
      fHistNCandidates->Fill(1);

      TObjArray arrTracks(3);
      for(Int_t jdau=0; jdau<3; jdau++){
        AliAODTrack *tr=vHF->GetProng(aod,d,jdau);
        arrTracks.AddAt(tr,jdau);
      }
      if(!fRDCutsAnalysis->PreSelect(arrTracks)){
        fHistNCandidates->Fill(6);
        continue;
      }

      if(!(vHF->FillRecoCand(aod,d))) { //Fill the data members of the candidate only if they are empty.
        fHistNCandidates->Fill(2); //monitor how often this fails
        continue;
      }

      Bool_t goodDaus=kTRUE;
      for(Int_t jdau=0; jdau<3; jdau++){
	Int_t idStored=d->GetProngID(jdau);
	TObject* odau=(TObject*)d->GetDaughter(jdau);
	TString cname=odau->ClassName();
	if(cname.Contains("AliAODTrack")){
	  AliAODTrack* tr=(AliAODTrack*)d->GetDaughter(jdau);
	  if(tr->Charge()==0){
	    fDaughterClass->Fill(1);
	    goodDaus=kFALSE;
	  }
	  Int_t idRecov=tr->GetID();
	  Int_t dId=idRecov-idStored;
	  fIDDauVsIDTra->Fill(idStored,idRecov);
	  fDeltaID->Fill(dId);
	  if(dId!=0){
	    goodDaus=kFALSE;
	    fDaughterClass->Fill(4);
	  }else{
	    if(idStored>0) fDaughterClass->Fill(0);
	    else if(idStored==0) fDaughterClass->Fill(2);
	    else{
	      goodDaus=kFALSE;
	      fDaughterClass->Fill(3);
	    }
	  }
	}else if(cname.Contains("AliAODRecoDecayHF2")){
	  fDaughterClass->Fill(5);
	  goodDaus=kFALSE;
	}else if(cname.Contains("AliAODRecoDecayHF3")){
	  fDaughterClass->Fill(6);
	  goodDaus=kFALSE;
	}else if(cname.Contains("AliAODRecoCascadeHF")){
	  fDaughterClass->Fill(7);
	  goodDaus=kFALSE;
	}else{
	  fDaughterClass->Fill(8);
	  goodDaus=kFALSE;
	}
      }
      Double_t ptCand = d->Pt();
      Double_t invMass=d->InvMassDplus();
      if(goodDaus) fPtVsMassGoodDaus->Fill(invMass,ptCand);
      else fPtVsMassBadDaus->Fill(invMass,ptCand);

      Int_t passTopolAndPIDCuts=fRDCutsAnalysis->IsSelected(d,AliRDHFCuts::kAll,aod);
      Bool_t passSingleTrackCuts=kTRUE;
      if(fRDCutsAnalysis->GetIsSelectedCuts() && fRDCutsAnalysis->GetIsSelectedPID() && !passTopolAndPIDCuts) passSingleTrackCuts=kFALSE;

      if(!fRDCutsAnalysis->GetIsSelectedCuts()) continue;
      fHistNCandidates->Fill(3);
      if(!passSingleTrackCuts) continue;
      fHistNCandidates->Fill(4);
      if(passTopolAndPIDCuts) fHistNCandidates->Fill(5);

      Double_t etaD=d->Eta();
      Double_t phiD=d->Phi();
      if(fEtaSelection!=0){
	if(fEtaSelection==1 && etaD<0) continue;
	if(fEtaSelection==-1 && etaD>0) continue;
      }

      Bool_t unsetvtx=kFALSE;
      if(!d->GetOwnPrimaryVtx()){
	d->SetOwnPrimaryVtx(vtx1);
	unsetvtx=kTRUE;
      }

      Int_t iPtBin = fRDCutsAnalysis->PtBin(ptCand);

      Bool_t recVtx=kFALSE;
      AliAODVertex *origownvtx=0x0;
      if(fRDCutsAnalysis->GetIsPrimaryWithoutDaughters()){
	if(d->GetOwnPrimaryVtx()) origownvtx=new AliAODVertex(*d->GetOwnPrimaryVtx());
	if(fRDCutsAnalysis->RecalcOwnPrimaryVtx(d,aod))recVtx=kTRUE;
	else fRDCutsAnalysis->CleanOwnPrimaryVtx(d,aod,origownvtx);
      }

      Int_t labDp=-1;
      Bool_t isPrimary=kFALSE;
      Bool_t isFeeddown=kFALSE;
      Float_t pdgCode=-2;
      Float_t trueImpParXY=0.;
      Double_t ptB=-1.5;
      if(fReadMC){
	labDp = d->MatchToMC(411,arrayMC,3,pdgDgDplustoKpipi);
	if(labDp>=0){
	  AliAODMCParticle *partDp = (AliAODMCParticle*)arrayMC->At(labDp);
	  Int_t orig=AliVertexingHFUtils::CheckOrigin(arrayMC,partDp,fUseQuarkTagInKine);//Prompt = 4, FeedDown = 5
	  pdgCode=TMath::Abs(partDp->GetPdgCode());
	  if(orig==4){
	    isPrimary=kTRUE;
	    isFeeddown=kFALSE;
	  }else if(orig==5){
	    isPrimary=kFALSE;
	    isFeeddown=kTRUE;
	    trueImpParXY=GetTrueImpactParameter(mcHeader,arrayMC,partDp)*10000.;
	    ptB=AliVertexingHFUtils::GetBeautyMotherPt(arrayMC,partDp);
	  }else{
	    pdgCode=-3;
	  }
	}else{
	  pdgCode=-1;
	}
      }

      Double_t rapid=d->YDplus();
      fYVsPtNoPid->Fill(ptCand,rapid,invMass);
      if(passTopolAndPIDCuts) {
	fYVsPt->Fill(ptCand,rapid,invMass);
	nOS++;
      }
      if(fReadMC && labDp>=0){
	fYVsPtSigNoPid->Fill(ptCand,rapid, invMass);
	if(passTopolAndPIDCuts)fYVsPtSig->Fill(ptCand,rapid, invMass);
      }

      Bool_t isFidAcc=fRDCutsAnalysis->IsInFiducialAcceptance(ptCand,rapid);
      if(isFidAcc){

 	Double_t  minPtDau=999.,ptmax=0,maxdca=0,sigvert=0,sumD02=0;
	Double_t dlen=0,cosp=0,dlenxy=0,cospxy=0, ndlenxy=0, dd0max=0;
	if(fCutsDistr||fFillNtuple||fDoImpPar||fDoSparse){
	  dlen=d->DecayLength();
	  cosp=d->CosPointingAngle();
	  sumD02=d->Getd0Prong(0)*d->Getd0Prong(0)+d->Getd0Prong(1)*d->Getd0Prong(1)+d->Getd0Prong(2)*d->Getd0Prong(2);
	  maxdca=-9999.;
	  for(Int_t idau=0;idau<3;idau++) if(d->GetDCA(idau)>maxdca) maxdca=d->GetDCA(idau);
	  sigvert=d->GetSigmaVert();
	  ptmax=0;
	  dlenxy = d->DecayLengthXY();
	  ndlenxy=d->NormalizedDecayLengthXY();
	  cospxy=d->CosPointingAngleXY();
	  for(Int_t i=0; i<3; i++) {
	    if(d->PtProng(i)>ptmax)ptmax=d->PtProng(i);
	    if(d->PtProng(i)<minPtDau) minPtDau=d->PtProng(i);
	    Double_t diffIP, errdiffIP;
	    d->Getd0MeasMinusExpProng(i,aod->GetMagneticField(),diffIP,errdiffIP);
	    Double_t normdd0= diffIP/errdiffIP;
	    if(i==0) dd0max=normdd0;
	    else if(TMath::Abs(normdd0)>TMath::Abs(dd0max)) dd0max=normdd0;
	  }
	}
	Double_t impparXY=d->ImpParXY()*10000.;
	Double_t resSel=0;
	if(fRDCutsAnalysis->GetIsSelectedCuts()) resSel=1;
	if(passTopolAndPIDCuts) resSel=2;
	if(fSystem==1) dd0max = TMath::Abs(dd0max);

	//for all THnSparses except for FD
	Double_t arrayForSparse[kVarForSparse]={invMass,ptCand,TMath::Abs(impparXY),resSel,minPtDau,sigvert,cosp,cospxy,dlen,dlenxy,ndlenxy,dd0max,(Double_t)tracklets};
	//for THnSparses for FD
	Double_t arrayForSparseFD[kVarForSparseFD]={invMass,ptCand,TMath::Abs(impparXY),resSel,minPtDau,sigvert,cosp,cospxy,dlen,dlenxy,ndlenxy,dd0max,(Double_t)tracklets,ptB};

	//for imppar THnSparses
	Double_t arrayForImpPar[kVarForImpPar]={invMass,ptCand,impparXY};
	//for imppar THnSparse with true FD imppar
	Double_t arrayForImpParFDTrue[kVarForImpPar]={invMass,ptCand,trueImpParXY};

	Double_t flagOrigin = 0;

	AliAODTrack *track;
	if(fDoTrackVarHist) {
	  for(int i = 0; i < 3; i++) {
	    Double_t ptTrack = 0, nCrossedRowsTPC = 0, nClustersTPC = 0, ratioCRowsFClu = 0, isSig = 0;
	    track = (AliAODTrack*)d->GetDaughter(i);
	    AliESDtrack esdTrack(track);
	    esdTrack.SetTPCClusterMap(track->GetTPCClusterMap());
	    esdTrack.SetTPCSharedMap(track->GetTPCSharedMap());
	    esdTrack.SetTPCPointsF(track->GetTPCNclsF());
	    ptTrack = track->Pt();
	    nCrossedRowsTPC = esdTrack.GetTPCCrossedRows();
	    nClustersTPC = esdTrack.GetTPCNcls();
	    if(esdTrack.GetTPCNclsF()>0) {
	      ratioCRowsFClu = nCrossedRowsTPC / esdTrack.GetTPCNclsF();
	      if(labDp>=0){
		if(isPrimary) isSig = 1.;
		else if(isFeeddown) isSig = 2.;
	      }
	    }
	    Double_t arrayForTrackSparse[kVarForTrackSparse]={ptCand,invMass,ptTrack,nClustersTPC,nCrossedRowsTPC,ratioCRowsFClu,isSig};
	    if(passTopolAndPIDCuts){
	      if(fDoTrackVarHist) fHistTrackVar->Fill(arrayForTrackSparse);
	    }
	  }
	}

	//Fill histos
	index=GetHistoIndex(iPtBin);
	nSelectednopid++;
	fPtVsMassNoPid->Fill(invMass,ptCand);
	fMassHistNoPid[index]->Fill(invMass);
	if(fDoImpPar && passTopolAndPIDCuts){
	  fHistMassPtImpPar[0]->Fill(arrayForImpPar);
	}
	if(fDoSparse && (!fReadMC || !fFillOnlySignalSparses)){ //fill in case of false fReadMC or false fFillOnlySignalSparses
	  fSparseCutVars[0]->Fill(arrayForSparse);
	}
	if(passTopolAndPIDCuts){
	  nSelected++;
	  fPtVsMass->Fill(invMass,ptCand);
	  fMassHist[index]->Fill(invMass);
	  if(d->GetCharge()>0){
	    fPtVsMassPlus->Fill(invMass,ptCand);
	    fMassHistPlus[index]->Fill(invMass);
	  }
	  else if(d->GetCharge()<0){
	    fPtVsMassMinus->Fill(invMass,ptCand);
	    fMassHistMinus[index]->Fill(invMass);
	  }
	  fPhiEtaCand->Fill(etaD,phiD);
	  if(TMath::Abs(invMass-1.8696)<0.05) fPhiEtaCandSigReg->Fill(etaD,phiD);
	  if(fCutsDistr){
	    fCosPHist[index]->Fill(cosp);
	    fDLenHist[index]->Fill(dlen);
	    fSumd02Hist[index]->Fill(sumD02);
	    fSigVertHist[index]->Fill(sigvert);
	    fPtMaxHist[index]->Fill(ptmax);
	    fPtKHist[index]->Fill(d->PtProng(1));
	    fPtpi1Hist[index]->Fill(d->PtProng(0));
	    fPtpi2Hist[index]->Fill(d->PtProng(2));
	    fDCAHist[index]->Fill(maxdca);
	    fDLxy[index]->Fill(ndlenxy);
	    fCosxy[index]->Fill(cospxy);
	    fCorreld0Kd0pi[0]->Fill(d->Getd0Prong(0)*d->Getd0Prong(1),d->Getd0Prong(2)*d->Getd0Prong(1));
	  }
	}

	// fill ntuple
	if(fFillNtuple==1){

	  Float_t tmp[31];

	  tmp[0]=pdgCode;
	  if(isFeeddown) tmp[0]+=5000.;
	  tmp[1]=d->Px();
	  tmp[2]=d->Py();
	  tmp[3]=d->Pz();
	  tmp[4]=d->Pt();
	  tmp[5]=d->GetCharge();
	  //	tmp[5]=fRDCutsAnalysis->GetPIDBitMask(d);
	  tmp[6]=fRDCutsAnalysis->GetPIDTrackTPCTOFBitMap((AliAODTrack*)d->GetDaughter(0));
	  tmp[7]=fRDCutsAnalysis->GetPIDTrackTPCTOFBitMap((AliAODTrack*)d->GetDaughter(1));
	  tmp[8]=fRDCutsAnalysis->GetPIDTrackTPCTOFBitMap((AliAODTrack*)d->GetDaughter(2));
	  tmp[9]=d->PtProng(0);
	  tmp[10]=d->PtProng(1);
	  tmp[11]=d->PtProng(2);
	  tmp[12]=d->PProng(0);
	  tmp[13]=d->PProng(1);
	  tmp[14]=d->PProng(2);
	  tmp[15]=cosp;
	  tmp[16]=cospxy;
	  tmp[17]=dlen;
	  tmp[18]=d->NormalizedDecayLength();
	  tmp[19]=dlenxy;
	  tmp[20]=ndlenxy;
	  tmp[21]=d->InvMassDplus();
	  tmp[22]=sigvert;
	  tmp[23]=d->Getd0Prong(0);
	  tmp[24]=d->Getd0Prong(1);
	  tmp[25]=d->Getd0Prong(2);
	  tmp[26]=maxdca;
	  tmp[27]=ntracks;
	  tmp[28]=fRDCutsAnalysis->GetCentrality(aod);
	  tmp[29]=runNumber;
	  tmp[30]=d->HasBadDaughters();
	  fNtupleDplus->Fill(tmp);
	  PostData(4,fNtupleDplus);
	}

	if(fFillNtuple==2 && passTopolAndPIDCuts){
	  Float_t tmp[5];
	  tmp[0]=pdgCode;
	  if(isFeeddown) tmp[0]+=5000.;
	  tmp[1]=d->Pt();
	  tmp[2]=d->InvMassDplus();
	  tmp[3]=impparXY;
	  if(!fReadMC){flagOrigin=0;};
	  if(fReadMC){
	    if(isPrimary&&labDp>=0)flagOrigin=1;
	    if(isFeeddown&&labDp>=0)flagOrigin=2;
	    if(!(labDp>=0))flagOrigin=3;

	  }

	  tmp[4]=flagOrigin;
	  fNtupleDplus->Fill(tmp);
	  PostData(4,fNtupleDplus);
	}

	if(fReadMC){
	  if(labDp>=0) {
	    index=GetSignalHistoIndex(iPtBin);
	    if(fDoImpPar && passTopolAndPIDCuts){
	      if(isPrimary) fHistMassPtImpPar[1]->Fill(arrayForImpPar);
	      else if(isFeeddown){
		fHistMassPtImpPar[2]->Fill(arrayForImpPar);
		fHistMassPtImpPar[3]->Fill(arrayForImpParFDTrue);
	      }
	    }
	    if(fDoSparse){
	      if(isPrimary) fSparseCutVars[1]->Fill(arrayForSparse);
	      else if(isFeeddown){
		fSparseCutVars[2]->Fill(arrayForSparseFD);
	      }
	    }
	  }else{
	    index=GetBackgroundHistoIndex(iPtBin);
	    if(fDoImpPar && passTopolAndPIDCuts)fHistMassPtImpPar[4]->Fill(arrayForImpPar);
	  }
	  fMassHistNoPid[index]->Fill(invMass);
	  if(passTopolAndPIDCuts){
	    if(fCutsDistr){
	      Float_t fact=1.;
	      Float_t factor[3]={1.,1.,1.};
	      if(fUseStrangeness) fact=GetStrangenessWeights(d,arrayMC,factor);
	      fCosPHist[index]->Fill(cosp,fact);
	      fDLenHist[index]->Fill(dlen,fact);
	      fDLxy[index]->Fill(ndlenxy);
	      fCosxy[index]->Fill(cospxy);
	      Float_t sumd02s=d->Getd0Prong(0)*d->Getd0Prong(0)*factor[0]*factor[0]+d->Getd0Prong(1)*d->Getd0Prong(1)*factor[1]*factor[1]+d->Getd0Prong(2)*d->Getd0Prong(2)*factor[2]*factor[2];
	      fSumd02Hist[index]->Fill(sumd02s);
	      fSigVertHist[index]->Fill(sigvert,fact);
	      fPtMaxHist[index]->Fill(ptmax,fact);
	      fPtKHist[index]->Fill(d->PtProng(1),fact);
	      fPtpi1Hist[index]->Fill(d->PtProng(0),fact);
	      fPtpi2Hist[index]->Fill(d->PtProng(2),fact);
	      fDCAHist[index]->Fill(maxdca,fact);
	      fCorreld0Kd0pi[1]->Fill(d->Getd0Prong(0)*d->Getd0Prong(1),d->Getd0Prong(2)*d->Getd0Prong(1));
	    }
	    fMassHist[index]->Fill(invMass);
	    if(d->GetCharge()>0) fMassHistPlus[index]->Fill(invMass);
	    else if(d->GetCharge()<0) fMassHistMinus[index]->Fill(invMass);
	  }
	}
      }

      if(recVtx)fRDCutsAnalysis->CleanOwnPrimaryVtx(d,aod,origownvtx);

      if(unsetvtx) d->UnsetOwnPrimaryVtx();
    }
    fCounter->StoreCandidates(aod,nSelectednopid,kTRUE);
    fCounter->StoreCandidates(aod,nSelected,kFALSE);
  }
  delete vHF;
  //start LS analysis
  if(fDoLS && arrayLikeSign) LSAnalysis(array3Prong,arrayLikeSign,aod,vtx1,nOS);

  PostData(1,fOutput);
  PostData(2,fListCuts);
  PostData(3,fCounter);
  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEDplus::CreateLikeSignHistos(){
  /// Histos for Like Sign bckground

  TString hisname;
  Int_t indexLS=0;
  Int_t index=0;
  Int_t nbins=GetNBinsHistos();

  fMassHistLS = new TH1F*[5*fNPtBins];
  fCosPHistLS = new TH1F*[3*fNPtBins];
  fDLenHistLS = new TH1F*[3*fNPtBins];
  fSumd02HistLS = new TH1F*[3*fNPtBins];
  fSigVertHistLS = new TH1F*[3*fNPtBins];
  fPtMaxHistLS = new TH1F*[3*fNPtBins];
  fDCAHistLS = new TH1F*[3*fNPtBins];

  for(Int_t i=0;i<fNPtBins;i++){

    index=GetHistoIndex(i);
    indexLS=GetLSHistoIndex(i);

    hisname.Form("hLSPt%d",i);
    fMassHistLS[indexLS] = new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHistLS[indexLS]->Sumw2();

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

    hisname.Form("hLSPt%dnw",i);
    fMassHistLS[indexLS]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHistLS[indexLS]->Sumw2();

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

    indexLS++;
    hisname.Form("hLSPt%dspc",i);
    fMassHistLS[indexLS]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHistLS[indexLS]->Sumw2();
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
  }
}

//________________________________________________________________________
void AliAnalysisTaskSEDplus::CreateImpactParameterHistos(){
  /// Histos for impact parameter study

  Int_t nmassbins=GetNBinsHistos();

  Int_t nptbins=80;
  Double_t ptmin=0.;
  Double_t ptmax=40.;

  //dimensions for THnSparse
  TString axTit[kVarForImpPar]={"M_{K#pi#pi} (GeV/c^{2})","p_{T} (GeV/c)","Imp Par (#mum)"};

  Int_t nbins[kVarForImpPar]={nmassbins,nptbins,fNImpParBins};
  Double_t xmin[kVarForImpPar]={fLowmasslimit,ptmin,fLowerImpPar};
  Double_t xmax[kVarForImpPar]={fUpmasslimit,ptmax,fHigherImpPar};

  //mass, pt, imppar, PIDsel
  //mass, pt, imppar, PIDsel
  fHistMassPtImpPar[0]=new THnSparseF("hMassPtImpParAll",
					"Mass vs. pt vs.imppar - All",
					kVarForImpPar,nbins,xmin,xmax);
  fHistMassPtImpPar[1]=new THnSparseF("hMassPtImpParPrompt",
					"Mass vs. pt vs.imppar - promptD",
					kVarForImpPar,nbins,xmin,xmax);
  fHistMassPtImpPar[2]=new THnSparseF("hMassPtImpParBfeed",
					"Mass vs. pt vs.imppar - DfromB",
					kVarForImpPar,nbins,xmin,xmax);
  fHistMassPtImpPar[3]=new THnSparseF("hMassPtImpParTrueBfeed",
					"Mass vs. pt vs.true imppar -DfromB",
					kVarForImpPar,nbins,xmin,xmax);
  fHistMassPtImpPar[4]=new THnSparseF("hMassPtImpParBkg",
				        "Mass vs. pt vs.imppar - backgr.",
					kVarForImpPar,nbins,xmin,xmax);

  for(Int_t i=0; i<5; i++){
    for(Int_t iax=0; iax<kVarForImpPar; iax++) fHistMassPtImpPar[i]->GetAxis(iax)->SetTitle(axTit[iax].Data());
    fOutput->Add(fHistMassPtImpPar[i]);
  }
}

//________________________________________________________________________
void AliAnalysisTaskSEDplus::CreateCutVarsSparses(){
  /// Sparses for cut variation study

  Int_t nmassbins=GetNBinsHistos();

  Int_t nptbins=100;
  Double_t ptmin=0.;
  Double_t ptmax=50.;
  if(fUseFinPtBinsForSparse)
    nptbins = 500;

  Int_t nselbins=2;
  Double_t minsel=0.5;
  Double_t maxsel=2.5;

  Int_t nnormdlbins=30;
  Double_t minnormdl=0.;
  Double_t maxnormdl=30.;

  Int_t nmultbins;
  Double_t maxmult;
  Double_t minmult=-0.5;

  Int_t nd0bins;
  Double_t d0min;
  Double_t d0max;

  Int_t nptmindaubins;
  Double_t minptmindau;
  Double_t maxptmindau;

  Int_t nsigvertbins;
  Double_t minsigvert;
  Double_t maxsigvert;

  Int_t ndeclbins;
  Double_t mindecl;
  Double_t maxdecl;

  Int_t ndeclxybins;
  Double_t mindeclxy;
  Double_t maxdeclxy;

  Int_t ncospbins;
  Double_t mincosp;
  Double_t maxcosp;

  Int_t ncospxybins;
  Double_t mincospxy;
  Double_t maxcospxy;

  Int_t nd0d0expbins;
  Double_t mind0d0;
  Double_t maxd0d0;

  if(fSystem==1) {
    nd0bins=18;
    d0min=0.;
    d0max=180;

    nptmindaubins=1;//dummy axis
    minptmindau=0.2;
    maxptmindau=1.2;

    nsigvertbins=10;
    minsigvert=0.012;
    maxsigvert=0.032;

    ndeclbins=35;
    mindecl=0.;
    maxdecl=0.35;

    ndeclxybins=35;
    mindeclxy=0.;
    maxdeclxy=0.35;

    ncospbins=30;
    mincosp=0.97;
    maxcosp=1.;

    ncospxybins=30;
    mincospxy=0.97;
    maxcospxy=1.;

    nd0d0expbins=12;
    mind0d0=0.;
    maxd0d0=6.;

    maxmult=5000.5;//dummy axis
    nmultbins=1;
  }
  else {
    nd0bins=60;
    d0min=0.;
    d0max=300;

    nptmindaubins=10;
    minptmindau=0.2;
    maxptmindau=1.2;

    nsigvertbins=25;
    minsigvert=0.010;
    maxsigvert=0.035;

    ndeclbins=70;
    mindecl=0.;
    maxdecl=0.70;

    ndeclxybins=70;
    mindeclxy=0.;
    maxdeclxy=0.70;

    ncospbins=100;
    mincosp=0.90;
    maxcosp=1.;

    ncospxybins=30;
    mincospxy=0.97;
    maxcospxy=1.;

    nd0d0expbins=40;
    mind0d0=-10.;
    maxd0d0=10.;
    
    maxmult=200.5;
    nmultbins = maxmult-minmult;
  }
  
  //dimensions for THnSparse which are NOT for BFeed
  TString axTit[kVarForSparse]={"M_{K#pi#pi} (GeV/c^{2})","p_{T} (GeV/c)","Imp Par (#mum)","passTopolPID","min. daughter p_{T} (GeV/c)","sigmaVertex","cos(#theta_{P})","cos(#theta_{P}^{xy})","decL (cm)","decL XY (cm)","Norm decL XY","Norm max d0-d0exp","N_{trkls}"};

  Int_t nbins[kVarForSparse]={nmassbins,nptbins,nd0bins,nselbins,nptmindaubins,nsigvertbins,ncospbins,ncospxybins,ndeclbins,ndeclxybins,nnormdlbins,nd0d0expbins,nmultbins};
  Double_t xmin[kVarForSparse]={fLowmasslimit,ptmin,d0min,minsel,minptmindau,minsigvert,mincosp,mincospxy,mindecl,mindeclxy,minnormdl,mind0d0,minmult};
  Double_t xmax[kVarForSparse]={fUpmasslimit,ptmax,d0max,maxsel,maxptmindau,maxsigvert,maxcosp,maxcospxy,maxdecl,maxdeclxy,maxnormdl,maxd0d0,maxmult};

  //dimensions for THnSparse for BFeed
  Int_t nbinsFD[kVarForSparseFD]={nmassbins,nptbins,nd0bins,nselbins,nptmindaubins,nsigvertbins,ncospbins,ncospxybins,ndeclbins,ndeclbins,nnormdlbins,nd0d0expbins,nmultbins,84};
  Double_t xminFD[kVarForSparseFD]={fLowmasslimit,ptmin,d0min,minsel,minptmindau,minsigvert,mincosp,mincospxy,mindecl,mindeclxy,minnormdl,mind0d0,minmult,-2};
  Double_t xmaxFD[kVarForSparseFD]={fUpmasslimit,ptmax,d0max,maxsel,maxptmindau,maxsigvert,maxcosp,maxcospxy,maxdecl,maxdeclxy,maxnormdl,maxd0d0,maxmult,40};

  //mass, pt, imppar, cosPoinXY, decLXY, norm decLXY (for BFeed also ptB)
  fSparseCutVars[0]=new THnSparseF("hMassPtCutVarsAll",
					"Mass vs. pt vs. cut vars - All",
					kVarForSparse,nbins,xmin,xmax);
  fSparseCutVars[1]=new THnSparseF("hMassPtCutVarsPrompt",
					"Mass vs. pt vs. cut vars - promptD",
					kVarForSparse,nbins,xmin,xmax);
  fSparseCutVars[2]=new THnSparseF("hMassPtCutVarsBfeed",
					"Mass vs. pt vs. cut vars - DfromB",
					kVarForSparseFD,nbinsFD,xminFD,xmaxFD);

  for(Int_t i=0; i<3; i++){
    for(Int_t iax=0; iax<kVarForSparse; iax++) fSparseCutVars[i]->GetAxis(iax)->SetTitle(axTit[iax].Data());
    if(i == 2 || i == 3) fSparseCutVars[i]->GetAxis(kVarForSparseFD-1)->SetTitle("p_{T}^{B} (GeV/c)");
    fOutput->Add(fSparseCutVars[i]);
  }
}

//________________________________________________________________________
void AliAnalysisTaskSEDplus::CreateTrackVarHistos(){
    /// Histos for single track variable cuts studies

    Int_t nmassbins = GetNBinsHistos();

    Int_t nbins[kVarForTrackSparse]   = {40,nmassbins,50,170,170,100,3};
    Double_t xmin[kVarForTrackSparse] = {0.,fLowmasslimit,0.,0.5,0.5,0.,-0.5};
    Double_t xmax[kVarForTrackSparse] = {40.,fUpmasslimit,5.,170.5,170.5,1.,2.5};

    //pt, y
    fHistTrackVar = new THnSparseF("hHistTrackVar","hHistTrackVar",kVarForTrackSparse,nbins,xmin,xmax);
    fHistTrackVar->GetAxis(0)->SetTitle("D^{+} p_{T} (GeV/c)");
    fHistTrackVar->GetAxis(1)->SetTitle("k#pi#pi inv. mass (GeV/c^{2})");
    fHistTrackVar->GetAxis(2)->SetTitle("track p_{T} (GeV/c)");
    fHistTrackVar->GetAxis(3)->SetTitle("N TPC clusters");
    fHistTrackVar->GetAxis(4)->SetTitle("N TPC cross. rows");
    fHistTrackVar->GetAxis(5)->SetTitle("TPC clust./cross.rows");
    fHistTrackVar->GetAxis(6)->SetTitle("Bkg = 0, Signal = 1");

    fOutput->Add(fHistTrackVar);
}


//________________________________________________________________________
void AliAnalysisTaskSEDplus::CreateMCAcceptanceHistos(){
  /// Histos for MC Acceptance histos

  const Int_t nVarPrompt = 3;
  const Int_t nVarFD = 4;

  Double_t multmin=-0.5;
  Double_t multmax;
  Int_t nmultbins;
  if(fSystem==0) {
    multmax=200.5;
    nmultbins=multmax-multmin;
  }
  else {
    multmax=5000.5;
    nmultbins=1;
  }
  
  Int_t nptbins = 100;
  Double_t ptmin = 0.;
  Double_t ptmax = 50.;
  if(fUseFinPtBinsForSparse)
    nptbins = 500;

  Int_t nbinsPrompt[nVarPrompt]={nptbins,100,nmultbins};
  Int_t nbinsFD[nVarFD]={nptbins,100,nmultbins,200};

  Double_t xminPrompt[nVarPrompt] = {ptmin,-1.,multmin};
  Double_t xmaxPrompt[nVarPrompt] = {ptmax,1.,multmax};

  Double_t xminFD[nVarFD] = {ptmin,-1.,multmin,0.};
  Double_t xmaxFD[nVarFD] = {ptmax,1.,multmax,40.};

  //pt, y
  fMCAccPrompt = new THnSparseF("hMCAccPrompt","kStepMCAcceptance pt vs. y vs. Ntracklets - promptD",nVarPrompt,nbinsPrompt,xminPrompt,xmaxPrompt);
  fMCAccPrompt->GetAxis(0)->SetTitle("p_{T} (GeV/c)");
  fMCAccPrompt->GetAxis(1)->SetTitle("y");
  fMCAccPrompt->GetAxis(2)->SetTitle("N_{trklts}");

  //pt,y,ptB
  fMCAccBFeed = new THnSparseF("hMCAccBFeed","kStepMCAcceptance pt vs. y vs. Ntracklets vs. ptB - DfromB",nVarFD,nbinsFD,xminFD,xmaxFD);
  fMCAccBFeed->GetAxis(0)->SetTitle("p_{T} (GeV/c)");
  fMCAccBFeed->GetAxis(1)->SetTitle("y");
  fMCAccBFeed->GetAxis(2)->SetTitle("N_{trklts}");
  fMCAccBFeed->GetAxis(3)->SetTitle("p_{T}^{B} (GeV/c)");

  fOutput->Add(fMCAccPrompt);
  fOutput->Add(fMCAccBFeed);
}

//________________________________________________________________________
void AliAnalysisTaskSEDplus::FillMCAcceptanceHistos(TClonesArray *arrayMC, AliAODMCHeader *mcHeader, Int_t tracklets){
  /// Fill MC acceptance histos for cuts study

  const Int_t nProng = 3;

  Double_t zMCVertex = mcHeader->GetVtxZ(); //vertex MC

  for(Int_t iPart=0; iPart<arrayMC->GetEntriesFast(); iPart++){
    AliAODMCParticle* mcPart = dynamic_cast<AliAODMCParticle*>(arrayMC->At(iPart));
    if (TMath::Abs(mcPart->GetPdgCode()) == 411){

      Int_t orig=AliVertexingHFUtils::CheckOrigin(arrayMC,mcPart,fUseQuarkTagInKine);//Prompt = 4, FeedDown = 5

      Int_t deca = 0;
      Bool_t isGoodDecay=kFALSE;
      Int_t labDau[4]={-1,-1,-1,-1};
      Bool_t isInAcc = kFALSE;
      Bool_t isFidAcc = kFALSE;

      deca=AliVertexingHFUtils::CheckDplusDecay(arrayMC,mcPart,labDau);
      if(deca > 0) isGoodDecay=kTRUE;

      if(labDau[0]==-1){
	continue; //protection against unfilled array of labels
      }

      isFidAcc=fRDCutsAnalysis->IsInFiducialAcceptance(mcPart->Pt(),mcPart->Y());
      isInAcc=CheckAcc(arrayMC,nProng,labDau);

      if(isGoodDecay && TMath::Abs(zMCVertex) < fRDCutsAnalysis->GetMaxVtxZ() && isFidAcc && isInAcc) {
	//for prompt
	if(orig == 4){
	  //fill histo for prompt
	  Double_t arrayMCprompt[3] = {mcPart->Pt(),mcPart->Y(),(Double_t)tracklets};
	  fMCAccPrompt->Fill(arrayMCprompt);
	}
	//for FD
	else if(orig == 5){
	  Double_t ptB = AliVertexingHFUtils::GetBeautyMotherPt(arrayMC,mcPart);
	  //fill histo for FD
	  Double_t arrayMCFD[4] = {mcPart->Pt(),mcPart->Y(),(Double_t)tracklets,ptB};
	  fMCAccBFeed->Fill(arrayMCFD);
	}
	else
	  continue;
      }
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskSEDplus::Terminate(Option_t */*option*/)
{
  /// Terminate analysis
  //
  if(fDebug > 1) printf("AnalysisTaskSEDplus: Terminate() \n");

  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {
    printf("ERROR: fOutput not available\n");
    return;
  }

  fHistNEvents = dynamic_cast<TH1F*>(fOutput->FindObject("fHistNEvents"));
  if(fHistNEvents){
    printf("Number of analyzed events = %d\n",(Int_t)fHistNEvents->GetBinContent(10));
  }else{
    printf("ERROR: fHistNEvents not available\n");
    return;
  }

  return;
}
//_________________________________________________________________________________________________
Float_t AliAnalysisTaskSEDplus::GetTrueImpactParameter(const AliAODMCHeader *mcHeader, TClonesArray* arrayMC, const AliAODMCParticle *partDp) const {
  /// true impact parameter calculation

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
  Int_t labelFirstDau = partDp->GetDaughterLabel(0);
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
	  Int_t labelFirstDauRes = part->GetDaughterLabel(0);
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
//_________________________________________________________________________________________________
Float_t AliAnalysisTaskSEDplus::GetStrangenessWeights(const AliAODRecoDecayHF3Prong* d, TClonesArray* arrayMC, Float_t factor[3]) const {
  /// Computes weights to adapt strangeness in MC to data

  for(Int_t iprong=0;iprong<3;iprong++){
    factor[iprong]=1;
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

  Float_t fact=1.;
  for(Int_t k=0;k<3;k++)fact=fact*factor[k];
  return fact;

}

//_________________________________________________________________
Bool_t AliAnalysisTaskSEDplus::CheckAcc(TClonesArray* arrayMC,Int_t nProng, Int_t *labDau){
  /// check if the decay products are in the good eta and pt range
  for (Int_t iProng = 0; iProng<nProng; iProng++){
    AliAODMCParticle* mcPartDaughter=dynamic_cast<AliAODMCParticle*>(arrayMC->At(labDau[iProng]));
    if(!mcPartDaughter) return kFALSE;
    Double_t eta = mcPartDaughter->Eta();
    Double_t pt = mcPartDaughter->Pt();
    if (TMath::Abs(eta) > 0.9 || pt < 0.1) return kFALSE;
  }
  return kTRUE;
}
