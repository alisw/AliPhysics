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
#include <THnSparse.h>
#include <TDatabasePDG.h>
#include <Riostream.h>

#include "AliAnalysisManager.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODRecoDecay.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliAnalysisVertexingHF.h"
#include "AliRDHFCutsDstoKKpi.h"
#include "AliAnalysisTaskSE.h"
#include "AliNormalizationCounter.h"
#include "AliAnalysisTaskSEDs.h"
#include "AliVertexingHFUtils.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskSEDs);
/// \endcond

//________________________________________________________________________
AliAnalysisTaskSEDs::AliAnalysisTaskSEDs():
  AliAnalysisTaskSE(),
  fOutput(0),
  fHistNEvents(0),
  fHistoPtWeight(0x0),
  fPtVsMass(0),
  fPtVsMassPhi(0),
  fPtVsMassK0st(0),
  fYVsPt(0),
  fYVsPtSig(0),
  fHistAllV0multNTPCout(0),
  fHistSelV0multNTPCout(0),
  fCosPHist3D(0x0),
  fCosPxyHist3D(0x0),
  fDLenHist3D(0x0),
  fDLenxyHist3D(0x0),
  fNDLenxyHist3D(0x0),
  fSigVertHist3D(0x0),
  fDCAHist3D(0x0),
  fNormIPHist3D(0x0),
  fCosPiDsHist3D(0x0),
  fCosPiKPhiHist3D(0x0),
  fPtProng0Hist3D(0x0),
  fPtProng1Hist3D(0x0),
  fPtProng2Hist3D(0x0),
  fNtupleDs(0),
  fFillNtuple(0),
  fSystem(0),
  fReadMC(kFALSE),
  fWriteOnlySignal(kFALSE),
  fDoCutVarHistos(kTRUE),
  fUseSelectionBit(kFALSE),
  fFillSparse(kTRUE),
  fFillSparseDplus(kFALSE),
  fFillImpParSparse(kFALSE),
  fDoRotBkg(kFALSE),
  fDoBkgPhiSB(kFALSE),
  fDoCutV0multTPCout(kFALSE),
  fUseWeight(kFALSE),
  fAODProtection(1),
  fNPtBins(0),
  fListCuts(0),
  fMassRange(0.8),
  fMassBinSize(0.002),
  fminMass(1.6),
  fmaxMass(2.5),
  fMaxDeltaPhiMass4Rot(0.010),
  fCounter(0),
  fAnalysisCuts(0),
  fnSparse(0),
  fnSparseIP(0),
  fImpParSparse(0x0),
  fImpParSparseMC(0x0)
{
  /// Default constructor
    
  for(Int_t i=0;i<3;i++){
    fHistCentrality[i]=0;
    fHistCentralityMult[i]=0;
  }
  for(Int_t i=0;i<4;i++) {
    fChanHist[i]=0;
  }
  for(Int_t i=0;i<4*kMaxPtBins;i++){
    fPtCandHist[i]=0;
    fMassHist[i]=0;
    fMassHistPhi[i]=0;
    fMassHistK0st[i]=0;
    fCosPHist[i]=0;
    fDLenHist[i]=0;
    fSumd02Hist[i]=0;
    fSigVertHist[i]=0;
    fPtMaxHist[i]=0;
    fDCAHist[i]=0;
    fPtProng0Hist[i]=0;
    fPtProng1Hist[i]=0;
    fPtProng2Hist[i]=0;
    fDalitz[i]=0;
    fDalitzPhi[i]=0;
    fDalitzK0st[i]=0;
  }
  for(Int_t i=0;i<kMaxPtBins;i++){
    fMassHistKK[i]=0;
    fMassHistKpi[i]=0;
    fMassRotBkgHistPhi[i]=0;
    fMassLSBkgHistPhi[i]=0;
    fMassRSBkgHistPhi[i]=0;
  }
  for(Int_t i=0;i<kMaxPtBins+1;i++){
    fPtLimits[i]=0;
  }
  for (Int_t i=0; i<4; i++) {
    fnSparseMC[i]=0;
    fnSparseMCDplus[i]=0;
  }
}

//________________________________________________________________________
AliAnalysisTaskSEDs::AliAnalysisTaskSEDs(const char *name,AliRDHFCutsDstoKKpi* analysiscuts,Int_t fillNtuple):
  AliAnalysisTaskSE(name),
  fOutput(0),
  fHistNEvents(0),
  fHistoPtWeight(0x0),
  fPtVsMass(0),
  fPtVsMassPhi(0),
  fPtVsMassK0st(0),
  fYVsPt(0),
  fYVsPtSig(0),
  fHistAllV0multNTPCout(0),
  fHistSelV0multNTPCout(0),
  fCosPHist3D(0x0),
  fCosPxyHist3D(0x0),
  fDLenHist3D(0x0),
  fDLenxyHist3D(0x0),
  fNDLenxyHist3D(0x0),
  fSigVertHist3D(0x0),
  fDCAHist3D(0x0),
  fNormIPHist3D(0x0),
  fCosPiDsHist3D(0x0),
  fCosPiKPhiHist3D(0x0),
  fPtProng0Hist3D(0x0),
  fPtProng1Hist3D(0x0),
  fPtProng2Hist3D(0x0),
  fNtupleDs(0),
  fFillNtuple(fillNtuple),
  fSystem(0),
  fReadMC(kFALSE),
  fWriteOnlySignal(kFALSE),
  fDoCutVarHistos(kTRUE),
  fUseSelectionBit(kFALSE),
  fFillSparse(kTRUE),
  fFillSparseDplus(kFALSE),
  fFillImpParSparse(kFALSE),
  fDoRotBkg(kTRUE),
  fDoBkgPhiSB(kTRUE),
  fDoCutV0multTPCout(kFALSE),
  fUseWeight(kFALSE),
  fAODProtection(1),
  fNPtBins(0),
  fListCuts(0),
  fMassRange(0.8),
  fMassBinSize(0.002),
  fminMass(1.6),
  fmaxMass(2.5),
  fMaxDeltaPhiMass4Rot(0.010),
  fCounter(0),
  fAnalysisCuts(analysiscuts),
  fnSparse(0),
  fnSparseIP(0),
  fImpParSparse(0x0),
  fImpParSparseMC(0x0)
{
  /// Default constructor
  /// Output slot #1 writes into a TList container
    
  for(Int_t i=0;i<3;i++){
    fHistCentrality[i]=0;
    fHistCentralityMult[i]=0;
  }
  for(Int_t i=0;i<4;i++) {
    fChanHist[i]=0;
  }
  for(Int_t i=0;i<4*kMaxPtBins;i++){
    fPtCandHist[i]=0;
    fMassHist[i]=0;
    fMassHistPhi[i]=0;
    fMassHistK0st[i]=0;
    fCosPHist[i]=0;
    fDLenHist[i]=0;
    fSumd02Hist[i]=0;
    fSigVertHist[i]=0;
    fPtMaxHist[i]=0;
    fDCAHist[i]=0;
    fPtProng0Hist[i]=0;
    fPtProng1Hist[i]=0;
    fPtProng2Hist[i]=0;
    fDalitz[i]=0;
    fDalitzPhi[i]=0;
    fDalitzK0st[i]=0;
  }
  for(Int_t i=0;i<kMaxPtBins;i++){
    fMassHistKK[i]=0;
    fMassHistKpi[i]=0;
    fMassRotBkgHistPhi[i]=0;
    fMassLSBkgHistPhi[i]=0;
    fMassRSBkgHistPhi[i]=0;
  }
  for(Int_t i=0;i<kMaxPtBins+1;i++){
    fPtLimits[i]=0;
  }
    
  for (Int_t i=0; i<4; i++) {
    fnSparseMC[i]=0;
    fnSparseMCDplus[i]=0;
  }
    
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
  /// define pt bins for analysis
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
  if(fOutput && !fOutput->IsOwner()){
    delete  fHistNEvents;
    delete  fHistoPtWeight;
    delete  fHistAllV0multNTPCout;
    delete  fHistSelV0multNTPCout;
    delete  fCosPHist3D;
    delete  fCosPxyHist3D;
    delete  fDLenHist3D;
    delete  fDLenxyHist3D;
    delete  fNDLenxyHist3D;
    delete  fSigVertHist3D;
    delete  fDCAHist3D;
    delete  fNormIPHist3D;
    delete  fCosPiDsHist3D;
    delete  fCosPiKPhiHist3D;
    delete  fPtProng0Hist3D;
    delete  fPtProng1Hist3D;
    delete  fPtProng2Hist3D;
      
    for(Int_t i=0;i<4;i++){
      delete fChanHist[i];
    }
    for(Int_t i=0;i<4*fNPtBins;i++){
      delete fMassHist[i];
      delete fMassHistPhi[i];
      delete fMassHistK0st[i];
      delete fCosPHist[i];
      delete fDLenHist[i];
      delete fSumd02Hist[i];
      delete fSigVertHist[i];
      delete fPtMaxHist[i];
      delete fPtCandHist[i];
      delete fDCAHist[i];
      delete fPtProng0Hist[i];
      delete fPtProng1Hist[i];
      delete fPtProng2Hist[i];
      delete fDalitz[i];
      delete fDalitzPhi[i];
      delete fDalitzK0st[i];
    }
    for(Int_t i=0;i<fNPtBins;i++){
      delete fMassHistKK[i];
      delete fMassHistKpi[i];
      delete fMassRotBkgHistPhi[i];
      delete fMassLSBkgHistPhi[i];
      delete fMassRSBkgHistPhi[i];
    }
    delete fPtVsMass;
    delete fPtVsMassPhi;
    delete fPtVsMassK0st;
    delete fYVsPt;
    delete fYVsPtSig;
    for(Int_t i=0;i<3;i++){
      delete fHistCentrality[i];
      delete fHistCentralityMult[i];
    }
        
    delete fnSparse;
    delete fnSparseIP;
    if(fFillImpParSparse) {
      delete fImpParSparse;
      delete fImpParSparseMC;
    }
    for (Int_t i=0; i<4; i++) {
      delete fnSparseMC[i];
      if(fFillSparseDplus)delete fnSparseMCDplus[i];
    }
  }
  delete fOutput;
  delete fListCuts;
  delete fNtupleDs;
  delete fCounter;
  delete fAnalysisCuts;
    
}

//________________________________________________________________________
void AliAnalysisTaskSEDs::Init()
{
  /// Initialization
    
  if(fDebug > 1) printf("AnalysisTaskSEDs::Init() \n");
    
  fListCuts=new TList();
  fListCuts->SetOwner();
  fListCuts->SetName("CutObjects");
    
  AliRDHFCutsDstoKKpi *analysis = new AliRDHFCutsDstoKKpi(*fAnalysisCuts);
  analysis->SetName("AnalysisCuts");
    
  fListCuts->Add(analysis);
  PostData(2,fListCuts);
  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEDs::UserCreateOutputObjects()
{
  /// Create the output container
  //
  if(fDebug > 1) printf("AnalysisTaskSEDs::UserCreateOutputObjects() \n");
    
  // Several histograms are more conveniently managed in a TList
  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("OutputHistos");
    
  fHistNEvents = new TH1F("hNEvents", "number of events ",15,-0.5,14.5);
  fHistNEvents->GetXaxis()->SetBinLabel(1,"nEventsRead");
  fHistNEvents->GetXaxis()->SetBinLabel(2,"nEvents Matched dAOD");
  fHistNEvents->GetXaxis()->SetBinLabel(3,"nEvents Mismatched dAOD");
  fHistNEvents->GetXaxis()->SetBinLabel(4,"nEventsAnal");
  fHistNEvents->GetXaxis()->SetBinLabel(5,"n. passing IsEvSelected");
  fHistNEvents->GetXaxis()->SetBinLabel(6,"n. rejected due to trigger");
  fHistNEvents->GetXaxis()->SetBinLabel(7,"n. rejected due to not reco vertex");
  fHistNEvents->GetXaxis()->SetBinLabel(8,"n. rejected for contr vertex");
  fHistNEvents->GetXaxis()->SetBinLabel(9,"n. rejected for vertex out of accept");
  fHistNEvents->GetXaxis()->SetBinLabel(10,"n. rejected for pileup events");
  fHistNEvents->GetXaxis()->SetBinLabel(11,"no. of out centrality events");
  fHistNEvents->GetXaxis()->SetBinLabel(12,"no. of 3 prong candidates");
  fHistNEvents->GetXaxis()->SetBinLabel(13,"no. of Ds after filtering cuts");
  fHistNEvents->GetXaxis()->SetBinLabel(14,"no. of Ds after selection cuts");
  fHistNEvents->GetXaxis()->SetBinLabel(15,"no. of not on-the-fly rec Ds");
    
  fHistNEvents->GetXaxis()->SetNdivisions(1,kFALSE);
    
  fHistNEvents->Sumw2();
  fHistNEvents->SetMinimum(0);
  fOutput->Add(fHistNEvents);
    
  if(fHistoPtWeight) {
    fHistoPtWeight->Sumw2();
    fOutput->Add(fHistoPtWeight);
  }
    
  fHistCentrality[0]=new TH1F("hCentr","centrality",10000,0.,100.);
  fHistCentrality[1]=new TH1F("hCentr(selectedCent)","centrality(selectedCent)",10000,0.,100.);
  fHistCentrality[2]=new TH1F("hCentr(OutofCent)","centrality(OutofCent)",10000,0.,100.);
  fHistCentralityMult[0]=new TH2F("hCentrMult","centrality vs mult",100,0.5,30000.5,40,0.,100.);
  fHistCentralityMult[1]=new TH2F("hCentrMult(selectedCent)","centrality vs mult(selectedCent)",100,0.5,30000.5,40,0.,100.);
  fHistCentralityMult[2]=new TH2F("hCentrMult(OutofCent)","centrality vs mult(OutofCent)",100,0.5,30000.5,40,0.,100.);
  for(Int_t i=0;i<3;i++){
    fHistCentrality[i]->Sumw2();
    fOutput->Add(fHistCentrality[i]);
    fHistCentralityMult[i]->Sumw2();
    fOutput->Add(fHistCentralityMult[i]);
  }
  if(fDoCutV0multTPCout) {
    fHistAllV0multNTPCout = new TH2F("HistAllV0multNTPCout", "V0mult vs # TPCout (all) ;V0mult ;# TPCout", 1000, 0., 40000, 1000, 0, 30000);
    fHistSelV0multNTPCout = new TH2F("HistSelV0multNTPCout", "V0mult vs # TPCout (sel) ;V0mult ;# TPCout", 1000, 0., 40000, 1000, 0, 30000);
    fHistAllV0multNTPCout->Sumw2();
    fHistSelV0multNTPCout->Sumw2();
    fOutput->Add(fHistAllV0multNTPCout);
    fOutput->Add(fHistSelV0multNTPCout);
  }
    
  Double_t massDs=TDatabasePDG::Instance()->GetParticle(431)->Mass();
    
  Int_t nInvMassBins=(Int_t)(fMassRange/fMassBinSize+0.5);
  if(nInvMassBins%2==1) nInvMassBins++;
  // Double_t minMass=massDs-1.0*nInvMassBins*fMassBinSize;
  Double_t minMass=massDs-0.5*nInvMassBins*fMassBinSize;
  //  Double_t maxMass=massDs+1.0*nInvMassBins*fMassBinSize;
  Double_t maxMass=massDs+0.5*nInvMassBins*fMassBinSize;
  fminMass = minMass;
  fmaxMass = maxMass;
    
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
      hisname.Form("hCosPxy%sPt%d",htype.Data(),i);
      fCosPxyHist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.5,1.);
      fCosPxyHist[index]->Sumw2();
      hisname.Form("hDLen%sPt%d",htype.Data(),i);
      fDLenHist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.,0.5);
      fDLenHist[index]->Sumw2();
      hisname.Form("hDLenxy%sPt%d",htype.Data(),i);
      fDLenxyHist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.,0.5);
      fDLenxyHist[index]->Sumw2();
      hisname.Form("hNDLenxy%sPt%d",htype.Data(),i);
      fNDLenxyHist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.,11.);
      fNDLenxyHist[index]->Sumw2();
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
      hisname.Form("hNormIP%sPt%d",htype.Data(),i);
      fNormIPHist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.,6.);
      fNormIPHist[index]->Sumw2();
      hisname.Form("hCosPiDs%sPt%d",htype.Data(),i);
      fCosPiDsHist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.5,1.);
      fCosPiDsHist[index]->Sumw2();
      hisname.Form("hCosPiKPhi%sPt%d",htype.Data(),i);
      fCosPiKPhiHist[index]=new TH1F(hisname.Data(),hisname.Data(),100,0.,0.5);
      fCosPiKPhiHist[index]->Sumw2();
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
    fOutput->Add(fPtCandHist[i]);
    if(fDoCutVarHistos){
      fOutput->Add(fCosPHist[i]);
      fOutput->Add(fCosPxyHist[i]);
      fOutput->Add(fDLenHist[i]);
      fOutput->Add(fDLenxyHist[i]);
      fOutput->Add(fNDLenxyHist[i]);
      fOutput->Add(fSumd02Hist[i]);
      fOutput->Add(fSigVertHist[i]);
      fOutput->Add(fPtMaxHist[i]);
      fOutput->Add(fDCAHist[i]);
      fOutput->Add(fNormIPHist[i]);
      fOutput->Add(fCosPiDsHist[i]);
      fOutput->Add(fCosPiKPhiHist[i]);
      fOutput->Add(fPtProng0Hist[i]);
      fOutput->Add(fPtProng1Hist[i]);
      fOutput->Add(fPtProng2Hist[i]);
      fOutput->Add(fDalitz[i]);
      fOutput->Add(fDalitzPhi[i]);
      fOutput->Add(fDalitzK0st[i]);
    }
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
    
  fCosPHist3D     = new TH3F("fCosPHist3D","CosP vs Ds mass",nInvMassBins,minMass,maxMass,20,0.,20.,100,0.5,1.);
  fCosPxyHist3D   = new TH3F("fCosPxyHist3D","CosPxy vs Ds mass",nInvMassBins,minMass,maxMass,20,0.,20.,100,0.5,1.);
  fDLenHist3D     = new TH3F("fDLenHist3D","DLen vs Ds mass",nInvMassBins,minMass,maxMass,20,0.,20.,100,0.,0.5);
  fDLenxyHist3D   = new TH3F("fDLenxyHist3D","DLenxy vs Ds mass",nInvMassBins,minMass,maxMass,20,0.,20.,100,0.,0.5);
  fNDLenxyHist3D  = new TH3F("fNDLenxyHist3D","NDLenxy vs Ds mass",nInvMassBins,minMass,maxMass,20,0.,20.,100,0.,11.);
  fSigVertHist3D  = new TH3F("fSigVertHist3D","SigVert vs Ds mass",nInvMassBins,minMass,maxMass,20,0.,20.,100,0.,0.1);
  fDCAHist3D      = new TH3F("fDCAHist3D","DCA vs Ds mass",nInvMassBins,minMass,maxMass,20,0.,20.,100,0.,0.1);
  fNormIPHist3D   = new TH3F("fNormIPHist3D","nIP vs Ds mass",nInvMassBins,minMass,maxMass,20,0.,20.,100,0.,6.);
  fCosPiDsHist3D  = new TH3F("fCosPiDsHist3D","CosPiDs vs Ds mass",nInvMassBins,minMass,maxMass,20,0.,20.,100,0.5,1.);
  fCosPiKPhiHist3D = new TH3F("fCosPiKPhiHist3D","CosPiKPhi vs Ds mass",nInvMassBins,minMass,maxMass,20,0.,20.,100,0.,0.5);
  fPtProng0Hist3D  = new TH3F("fPtProng0Hist3D","Pt prong0 vs Ds mass",nInvMassBins,minMass,maxMass,20,0.,20.,100,0.0,20.);
  fPtProng1Hist3D  = new TH3F("fPtProng1Hist3D","Pt prong1 vs Ds mass",nInvMassBins,minMass,maxMass,20,0.,20.,100,0.0,20.);
  fPtProng2Hist3D  = new TH3F("fPtProng2Hist3D","Pt prong2 vs Ds mass",nInvMassBins,minMass,maxMass,20,0.,20.,100,0.0,20.);
    
  if(!fReadMC && fDoCutVarHistos) {
    fOutput->Add(fCosPHist3D);
    fOutput->Add(fCosPxyHist3D);
    fOutput->Add(fDLenHist3D);
    fOutput->Add(fDLenxyHist3D);
    fOutput->Add(fNDLenxyHist3D);
    fOutput->Add(fSigVertHist3D);
    fOutput->Add(fDCAHist3D);
    fOutput->Add(fNormIPHist3D);
    fOutput->Add(fCosPiDsHist3D);
    fOutput->Add(fCosPiKPhiHist3D);
    fOutput->Add(fPtProng0Hist3D);
    fOutput->Add(fPtProng1Hist3D);
    fOutput->Add(fPtProng2Hist3D);
  }

  fPtVsMass=new TH2F("hPtVsMass","PtVsMass (prod. cuts)",nInvMassBins,minMass,maxMass,40,0.,20.);
  fPtVsMassPhi=new TH2F("hPtVsMassPhi","PtVsMass (phi selection)",nInvMassBins,minMass,maxMass,200,0.,20.);
  fPtVsMassK0st=new TH2F("hPtVsMassK0st","PtVsMass (K0* selection)",nInvMassBins,minMass,maxMass,200,0.,20.);
  fYVsPt=new TH2F("hYVsPt","YvsPt (prod. cuts)",40,0.,20.,80,-2.,2.);
  fYVsPtSig=new TH2F("hYVsPtSig","YvsPt (MC, only sig., prod. cuts)",40,0.,20.,80,-2.,2.);
    
  for(Int_t i=0;i<fNPtBins;i++){
    hisname.Form("hMassKKPt%d",i);
    fMassHistKK[i]=new TH1F(hisname.Data(),hisname.Data(),200,0.95,1.35);
    fMassHistKK[i]->Sumw2();
    fOutput->Add(fMassHistKK[i]);
    hisname.Form("hMassKpiPt%d",i);
    fMassHistKpi[i]=new TH1F(hisname.Data(),hisname.Data(),200,0.7,1.1);
    fMassHistKpi[i]->Sumw2();
    fOutput->Add(fMassHistKpi[i]);
    if(fDoRotBkg) {
      hisname.Form("hMassAllPt%dphi_RotBkg",i);
      fMassRotBkgHistPhi[i]=new TH1F(hisname.Data(),hisname.Data(),nInvMassBins,minMass,maxMass);
      fMassRotBkgHistPhi[i]->Sumw2();
      fOutput->Add(fMassRotBkgHistPhi[i]);
    }
    if(fDoBkgPhiSB) {
      hisname.Form("fMassLSBkgHistPhiPt%d",i);
      fMassLSBkgHistPhi[i]=new TH1F(hisname.Data(),hisname.Data(),nInvMassBins,minMass,maxMass);
      fMassLSBkgHistPhi[i]->Sumw2();
      fOutput->Add(fMassLSBkgHistPhi[i]);
      hisname.Form("fMassRSBkgHistPhiPt%d",i);
      fMassRSBkgHistPhi[i]=new TH1F(hisname.Data(),hisname.Data(),nInvMassBins,minMass,maxMass);
      fMassRSBkgHistPhi[i]->Sumw2();
      fOutput->Add(fMassRSBkgHistPhi[i]);
    }
  }
    
  fOutput->Add(fPtVsMass);
  fOutput->Add(fPtVsMassPhi);
  fOutput->Add(fPtVsMassK0st);
  fOutput->Add(fYVsPt);
  fOutput->Add(fYVsPtSig);

  nInvMassBins=(Int_t)(0.7/fMassBinSize+0.5);
  minMass=massDs-0.5*nInvMassBins*fMassBinSize;
  maxMass=massDs+0.5*nInvMassBins*fMassBinSize;

  Int_t nBinsReco[knVarForSparse]   = {nInvMassBins,  20,     30,     14,    14,   20,    10,   10,     14,     6,     6,   12};
  Double_t xminReco[knVarForSparse] = {minMass,       0.,     0.,     0.,    0.,   0.,   90.,   90.,    0.,    7.,    0.,   0.};
  Double_t xmaxReco[knVarForSparse] = {maxMass,      20.,     15,    70.,   70.,  10.,  100.,  100.,   70.,   10.,    3.,   6.};
  TString  axis[knVarForSparse]     = {"invMassDsAllPhi","p_{T}","#Delta Mass(KK)","dlen","dlen_{xy}","normdl_{xy}","cosP","cosP_{xy}","sigVert","cosPiDs","|cosPiKPhi^{3}|","normIP"};
  if(fSystem == 1) { //pPb,PbPb
    nInvMassBins=(Int_t)(0.45/fMassBinSize+0.5);
    minMass=massDs-0.5*nInvMassBins*fMassBinSize;
    maxMass=massDs+0.5*nInvMassBins*fMassBinSize;
    nBinsReco[0] = nInvMassBins; //Ds mass
    xminReco[0]  = minMass;
    xmaxReco[0]  = maxMass;
        
    nBinsReco[1] = 16; //pt
    xminReco[1]  = 0.;
    xmaxReco[1]  = 16.;
        
    nBinsReco[2] =  12; //#Delta Mass(KK)
    xmaxReco[2]  = 12.;
        
    nBinsReco[3] = 7; //dlen
    nBinsReco[4] = 7; //dlenxy
    nBinsReco[5] = 10; //ndlenxy
        
    nBinsReco[6] =    6; //cosP
    xminReco[6]  =  97.;
    xmaxReco[6]  = 100.;
        
    nBinsReco[7] =    6; //cosPxy
    xminReco[7]  =  97.;
    xmaxReco[7]  = 100.;
  }
    
  Int_t nBinsAcc[knVarForSparseAcc]   = {20,   20};
  Double_t xminAcc[knVarForSparseAcc] = {0., -10.};
  Double_t xmaxAcc[knVarForSparseAcc] = {20,  10.};
    
  Int_t nBinsIP[knVarForSparseIP]   = { 20,  400,  400,  400,  400,  3};
  Double_t xminIP[knVarForSparseIP] = { 0., -10., -10., -10., -10., 0.};
  Double_t xmaxIP[knVarForSparseIP] = {20.,  10.,  10.,  10.,  10., 3.};
  TString axisIP[knVarForSparseIP]  = {"motherPt","maxNormImp","IP0","IP1","IP2","candType"};
    
  if(fFillSparse) {
        
    if(fReadMC) {
      TString label[knVarForSparseAcc] = {"fromC","fromB"};
      for (Int_t i=0; i<2; i++) {
	fnSparseMC[i] = new THnSparseF(Form("fnSparseAcc_%s",label[i].Data()),Form("MC nSparse (Acc.Step)- %s",label[i].Data()),
				       knVarForSparseAcc, nBinsAcc, xminAcc, xmaxAcc);
	fnSparseMC[i]->GetAxis(0)->SetTitle("p_{T} (GeV/c)");
	fnSparseMC[i]->GetAxis(1)->SetTitle("y");
	fOutput->Add(fnSparseMC[i]);
                
	//Dplus
	if(fFillSparseDplus) {
	  fnSparseMCDplus[i] = new THnSparseF(Form("fnSparseAccDplus_%s",label[i].Data()),Form("MC nSparse D^{+} (Acc.Step)- %s",label[i].Data()),
					      knVarForSparseAcc, nBinsAcc, xminAcc, xmaxAcc);
	  fnSparseMCDplus[i]->GetAxis(0)->SetTitle("p_{T} (GeV/c)");
	  fnSparseMCDplus[i]->GetAxis(1)->SetTitle("y");
	  fOutput->Add(fnSparseMCDplus[i]);
	}
      }
      for (Int_t i=2; i<4; i++) {
	fnSparseMC[i] = new THnSparseF(Form("fnSparseReco_%s",label[i-2].Data()),Form("MC nSparse (Reco Step)- %s",label[i-2].Data()),
				       knVarForSparse, nBinsReco, xminReco, xmaxReco);
	for (Int_t j=0; j<knVarForSparse; j++) {
	  fnSparseMC[i]->GetAxis(j)->SetTitle(Form("%s",axis[j].Data()));
	}
	fOutput->Add(fnSparseMC[i]);
                
	//Dplus
	if(fFillSparseDplus) {
	  fnSparseMCDplus[i] = new THnSparseF(Form("fnSparseRecoDplus_%s",label[i-2].Data()),Form("MC nSparse D^{+} (Reco Step)- %s",label[i-2].Data()),
					      knVarForSparse, nBinsReco, xminReco, xmaxReco);
	  for (Int_t j=0; j<knVarForSparse; j++) {
	    fnSparseMCDplus[i]->GetAxis(j)->SetTitle(Form("%s",axis[j].Data()));
	  }
	  fOutput->Add(fnSparseMCDplus[i]);
	}
      }
            
      fnSparseIP = new THnSparseF("fnSparseIP","nSparseIP", knVarForSparseIP, nBinsIP, xminIP, xmaxIP);
      for (Int_t j=0; j<knVarForSparseIP; j++) {
	fnSparseIP->GetAxis(j)->SetTitle(Form("%s",axisIP[j].Data()));
      }
      fnSparseIP->GetAxis(5)->SetTitle("candType (0.5=bkg; 1.5=prompt; 2.5=FD)");
      fOutput->Add(fnSparseIP);
    }
    else {
      fnSparse = new THnSparseF("fnSparse","nSparse", knVarForSparse, nBinsReco, xminReco, xmaxReco);
      for (Int_t j=0; j<knVarForSparse; j++) {
	fnSparse->GetAxis(j)->SetTitle(Form("%s",axis[j].Data()));
      }
      fOutput->Add(fnSparse);
    }
  }
    
  if(fFillImpParSparse) {
    Int_t nBinsImpPar[3]   = { 20,   200,  350};
    Double_t xminImpPar[3] = { 0.,    0.,  1.6};
    Double_t xmaxImpPar[3] = {20., 1000.,  2.3};
    TString axisImpPar[3]  = {"Pt","imp.par. (#mum)","invMassDsAllPhi"};
    if(!fReadMC) {
      fImpParSparse = new THnSparseF("fImpParSparse","ImpParSparse", 3, nBinsImpPar, xminImpPar, xmaxImpPar);
      for (Int_t j=0; j<3; j++) {
	fImpParSparse->GetAxis(j)->SetTitle(Form("%s",axisImpPar[j].Data()));
      }
      fOutput->Add(fImpParSparse);
    }
    else {
      nBinsImpPar[2] = 3;
      xminImpPar[2]  = 0.;
      xmaxImpPar[2]  = 3.;
      axisImpPar[2]  = "candType (0.5=bkg; 1.5=prompt; 2.5=FD)";
      fImpParSparseMC  = new THnSparseF("fImpParSparseMC","ImpParSparseMC", 3, nBinsImpPar, xminImpPar, xmaxImpPar);
      for (Int_t j=0; j<3; j++) {
	fImpParSparseMC->GetAxis(j)->SetTitle(Form("%s",axisImpPar[j].Data()));
      }
      fOutput->Add(fImpParSparseMC);
    }
  }
    
  //Counter for Normalization
  fCounter = new AliNormalizationCounter("NormalizationCounter");
  fCounter->Init();
    
  PostData(1,fOutput);
  PostData(3,fCounter);
    
  if(fFillNtuple>0){
    OpenFile(4); // 4 is the slot number of the ntuple
        
    fNtupleDs = new TNtuple("fNtupleDs","Ds","labDs:retcode:pdgcode0:Pt0:Pt1:Pt2:PtRec:P0:P1:P2:PidTrackBit0:PidTrackBit1:PidTrackBit2:PointingAngle:PointingAngleXY:DecLeng:DecLengXY:NorDecLeng:NorDecLengXY:InvMassKKpi:InvMasspiKK:sigvert:d00:d01:d02:dca:d0square:InvMassPhiKKpi:InvMassPhipiKK:InvMassK0starKKpi:InvMassK0starpiKK:cosinePiDsFrameKKpi:cosinePiDsFramepiKK:cosineKPhiFrameKKpi:cosineKPhiFramepiKK:centrality:runNumber");
        
  }
    
    
  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEDs::UserExec(Option_t */*option*/)
{
  /// Ds selection for current event, fill mass histos and selecetion variable histo
  /// separate signal and backgound if fReadMC is activated
    
  AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());
    
  fHistNEvents->Fill(0); // all events
  if(fAODProtection>=0){
    //   Protection against different number of events in the AOD and deltaAOD
    //   In case of discrepancy the event is rejected.
    Int_t matchingAODdeltaAODlevel = AliRDHFCuts::CheckMatchingAODdeltaAODevents();
    if (matchingAODdeltaAODlevel<0 || (matchingAODdeltaAODlevel==0 && fAODProtection==1)) {
      // AOD/deltaAOD trees have different number of entries || TProcessID do not match while it was required
      fHistNEvents->Fill(2);
      return;
    }
    fHistNEvents->Fill(1);
  }
    
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
    
    
  fHistNEvents->Fill(3); // count event
  // Post the data already here
  PostData(1,fOutput);
    
  fCounter->StoreEvent(aod,fAnalysisCuts,fReadMC);
    
    
  Bool_t isEvSel=fAnalysisCuts->IsEventSelected(aod);
  Float_t ntracks=aod->GetNumberOfTracks();
  Float_t evCentr=fAnalysisCuts->GetCentrality(aod);
    
  fHistCentrality[0]->Fill(evCentr);
  fHistCentralityMult[0]->Fill(ntracks,evCentr);
  if(fAnalysisCuts->IsEventRejectedDueToTrigger())fHistNEvents->Fill(5);
  if(fAnalysisCuts->IsEventRejectedDueToNotRecoVertex())fHistNEvents->Fill(6);
  if(fAnalysisCuts->IsEventRejectedDueToVertexContributors())fHistNEvents->Fill(7);
  if(fAnalysisCuts->IsEventRejectedDueToZVertexOutsideFiducialRegion())fHistNEvents->Fill(8);
  if(fAnalysisCuts->IsEventRejectedDueToPileup())fHistNEvents->Fill(9);
  if(fAnalysisCuts->IsEventRejectedDueToCentrality()){
    fHistNEvents->Fill(10);
    fHistCentrality[2]->Fill(evCentr);
    fHistCentralityMult[2]->Fill(ntracks,evCentr);
  }
    
  Float_t centrality=fAnalysisCuts->GetCentrality(aod);
  Int_t runNumber=aod->GetRunNumber();
    
    
    
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
    
    
  if(fReadMC && fFillSparse){
    if(aod->GetTriggerMask()==0 && (runNumber>=195344 && runNumber<=195677))
      // protection for events with empty trigger mask in p-Pb
      return;
    if(fAnalysisCuts->GetUseCentrality()>0 && fAnalysisCuts->IsEventSelectedInCentrality(aod)!=0)
      // events not passing the centrality selection can be removed immediately.
      return;
    Double_t zMCVertex = mcHeader->GetVtxZ();
    if (TMath::Abs(zMCVertex) > fAnalysisCuts->GetMaxVtxZ())
      return;
    FillMCGenAccHistos(arrayMC, mcHeader);
  }
    
  if(!isEvSel)return;
  fHistNEvents->Fill(4);
  fHistCentrality[1]->Fill(evCentr);
  fHistCentralityMult[1]->Fill(ntracks,evCentr);
    
  if(fDoCutV0multTPCout) {
    //cut on V0mult. vs #tracks kTPCout
    Float_t V0mult=0.;
    AliAODVZERO *aodVZERO = (AliAODVZERO*)aod->GetVZEROData();
    if(aodVZERO) {
      for(int ich=0; ich<64;ich++) V0mult += aodVZERO->GetMultiplicity(ich);
    }
    Int_t nTPCout=0;
    for (Int_t i=0;i<aod->GetNumberOfTracks();++i) {
      AliVTrack *track = aod->GetTrack(i);
      if (!track) continue;
      if((track->GetStatus() & AliVTrack::kTPCout)) nTPCout++;
    }
    fHistAllV0multNTPCout->Fill(V0mult,nTPCout);
    if(nTPCout > (0.32*V0mult+750)) return;
    else fHistSelV0multNTPCout->Fill(V0mult,nTPCout);
  }
    
  Int_t n3Prong = array3Prong->GetEntriesFast();
  if(fDebug>1) printf("Number of Ds->KKpi: %d\n",n3Prong);
    
    
  Int_t pdgDstoKKpi[3]={321,321,211};
  Int_t nSelected=0;
  Int_t nFiltered=0;
  Double_t massPhi=TDatabasePDG::Instance()->GetParticle(333)->Mass();
    
  // vHF object is needed to call the method that refills the missing info of the candidates
  // if they have been deleted in dAOD reconstruction phase
  // in order to reduce the size of the file
  AliAnalysisVertexingHF *vHF=new AliAnalysisVertexingHF();
    
  for (Int_t i3Prong = 0; i3Prong < n3Prong; i3Prong++) {
        
        
    AliAODRecoDecayHF3Prong *d = (AliAODRecoDecayHF3Prong*)array3Prong->UncheckedAt(i3Prong);
    fHistNEvents->Fill(11);
        
    if(fUseSelectionBit && !(d->HasSelectionBit(AliRDHFCuts::kDsCuts))){
      continue;
    }
    nFiltered++;
    fHistNEvents->Fill(12);
        
    if(!(vHF->FillRecoCand(aod,d))) {////Fill the data members of the candidate only if they are empty.
      fHistNEvents->Fill(14); //monitor how often this fails
      continue;
    }
        
    Bool_t unsetvtx=kFALSE;
    if(!d->GetOwnPrimaryVtx()){
      d->SetOwnPrimaryVtx(vtx1);
      unsetvtx=kTRUE;
      // NOTE: the ow primary vertex should be unset, otherwise there is a memory leak
      // Pay attention if you use continue inside this loop!!!
    }
        
    Bool_t recVtx=kFALSE;
    AliAODVertex *origownvtx=0x0;
        
    Double_t ptCand = d->Pt();
    Int_t iPtBin=TMath::BinarySearch(fNPtBins,fPtLimits,(Float_t)ptCand);
    Double_t rapid=d->YDs();
    fYVsPt->Fill(ptCand,rapid);
    Bool_t isFidAcc=fAnalysisCuts->IsInFiducialAcceptance(ptCand,rapid);
        
    if(isFidAcc){
            
      Int_t retCodeAnalysisCuts=fAnalysisCuts->IsSelected(d,AliRDHFCuts::kAll,aod);
      Int_t retCodeNoRes=retCodeAnalysisCuts;
      Bool_t origRes=fAnalysisCuts->IsCutOnResonancesApplied();
      if(origRes){
	fAnalysisCuts->ApplyCutOnResonances(kFALSE);
	retCodeNoRes=fAnalysisCuts->IsSelected(d,AliRDHFCuts::kAll,aod);
	fAnalysisCuts->ApplyCutOnResonances(origRes);
      }
      Double_t massKK = 0.;
            
      if(retCodeNoRes&1){ //KKpi
	massKK=d->InvMass2Prongs(0,1,321,321);
	Double_t massKp=d->InvMass2Prongs(1,2,321,211);
	fMassHistKK[iPtBin]->Fill(massKK);
	fMassHistKpi[iPtBin]->Fill(massKp);
      }
      if(retCodeNoRes&2){ //piKK
	massKK=d->InvMass2Prongs(1,2,321,321);
	Double_t massKp=d->InvMass2Prongs(0,1,211,321);
	fMassHistKK[iPtBin]->Fill(massKK);
	fMassHistKpi[iPtBin]->Fill(massKp);
      }
            
      Int_t isKKpi=retCodeAnalysisCuts&1;
      Int_t ispiKK=retCodeAnalysisCuts&2;
      Int_t isPhiKKpi=retCodeAnalysisCuts&4;
      Int_t isPhipiKK=retCodeAnalysisCuts&8;
      Int_t isK0starKKpi=retCodeAnalysisCuts&16;
      Int_t isK0starpiKK=retCodeAnalysisCuts&32;
            
      if(retCodeAnalysisCuts>0){
	if(fAnalysisCuts->GetIsPrimaryWithoutDaughters()){
	  if(d->GetOwnPrimaryVtx()) origownvtx=new AliAODVertex(*d->GetOwnPrimaryVtx());
	  if(fAnalysisCuts->RecalcOwnPrimaryVtx(d,aod))recVtx=kTRUE;
	  else fAnalysisCuts->CleanOwnPrimaryVtx(d,aod,origownvtx);
	}
                
                
	fHistNEvents->Fill(13);
	nSelected++;
                
	Int_t index=GetHistoIndex(iPtBin);
	fPtCandHist[index]->Fill(ptCand);
                
                
	Double_t weightKKpi=1.;
	Double_t weightpiKK=1.;
	if(fAnalysisCuts->GetPidOption()==AliRDHFCutsDstoKKpi::kBayesianWeights){
	  weightKKpi=fAnalysisCuts->GetWeightForKKpi();
	  weightpiKK=fAnalysisCuts->GetWeightForpiKK();
	  if(weightKKpi>1. || weightKKpi<0.) weightKKpi=0.;
	  if(weightpiKK>1. || weightpiKK<0.) weightpiKK=0.;
	}
                
	fChanHist[0]->Fill(retCodeAnalysisCuts);
                
                
	Double_t invMass = 0.;
	Int_t indexMCKKpi=-1;
	Int_t indexMCpiKK=-1;
	Int_t labDs=-1;
	Int_t labDplus=-1;
	Int_t pdgCode0=-999;
	Int_t isMCSignal=-1;
                
                
	if(fReadMC){
                    
	  labDs = d->MatchToMC(431,arrayMC,3,pdgDstoKKpi);
	  labDplus = d->MatchToMC(411,arrayMC,3,pdgDstoKKpi);
	  if(labDs>=0){
	    Int_t labDau0=((AliAODTrack*)d->GetDaughter(0))->GetLabel();
	    AliAODMCParticle* p=(AliAODMCParticle*)arrayMC->UncheckedAt(TMath::Abs(labDau0));
	    pdgCode0=TMath::Abs(p->GetPdgCode());
                        
	    if(isKKpi){
	      if(pdgCode0==321) {
		indexMCKKpi=GetSignalHistoIndex(iPtBin);
		fYVsPtSig->Fill(ptCand,rapid);
		fChanHist[1]->Fill(retCodeAnalysisCuts);
		isMCSignal=1;
	      }else{
		indexMCKKpi=GetReflSignalHistoIndex(iPtBin);
		fChanHist[3]->Fill(retCodeAnalysisCuts);
		isMCSignal=0;
	      }
	    }
	    if(ispiKK){
	      if(pdgCode0==211) {
		indexMCpiKK=GetSignalHistoIndex(iPtBin);
		fYVsPtSig->Fill(ptCand,rapid);
		fChanHist[1]->Fill(retCodeAnalysisCuts);
		isMCSignal=1;
	      }else{
		indexMCpiKK=GetReflSignalHistoIndex(iPtBin);
		fChanHist[3]->Fill(retCodeAnalysisCuts);
		isMCSignal=0;
	      }
	    }
	  }else{
	    indexMCpiKK=GetBackgroundHistoIndex(iPtBin);
	    indexMCKKpi=GetBackgroundHistoIndex(iPtBin);
	    fChanHist[2]->Fill(retCodeAnalysisCuts);
	    if(labDplus>=0) {
              Int_t labDau0=((AliAODTrack*)d->GetDaughter(0))->GetLabel();
              AliAODMCParticle* p=(AliAODMCParticle*)arrayMC->UncheckedAt(TMath::Abs(labDau0));
              pdgCode0=TMath::Abs(p->GetPdgCode());
	    }
	  }
	}
          
	Double_t candType = 0.5; //for bkg
	Double_t massPhi=TDatabasePDG::Instance()->GetParticle(333)->Mass();

	if(isKKpi){
	  if(fDoRotBkg && TMath::Abs(massKK-massPhi)<=fMaxDeltaPhiMass4Rot)GenerateRotBkg(d,1,iPtBin);
        
	  invMass=d->InvMassDsKKpi();
	  fMassHist[index]->Fill(invMass,weightKKpi);
	  fPtVsMass->Fill(invMass,ptCand,weightKKpi);
        
	  if(fDoBkgPhiSB && 0.010<TMath::Abs(massKK-massPhi)<0.030) {
            if(massKK<massPhi)fMassLSBkgHistPhi[iPtBin]->Fill(invMass);
            else fMassRSBkgHistPhi[iPtBin]->Fill(invMass);
	  }
        
	  if(isPhiKKpi){
	    fMassHistPhi[index]->Fill(invMass,weightKKpi);
	    fPtVsMassPhi->Fill(invMass,ptCand,weightKKpi);
	  }
	  if(isK0starKKpi){
	    fMassHistK0st[index]->Fill(invMass,weightKKpi);
	    fPtVsMassK0st->Fill(invMass,ptCand,weightKKpi);
	  }
	  if(fReadMC  && indexMCKKpi!=-1){
	    fMassHist[indexMCKKpi]->Fill(invMass,weightKKpi);
	    if(isPhiKKpi) {
	      fMassHistPhi[indexMCKKpi]->Fill(invMass,weightKKpi);
	      if(fFillSparse) {
		if(indexMCKKpi==GetSignalHistoIndex(iPtBin) || labDplus >= 0) {
		  AliAODMCParticle *partDs;
		  if(indexMCKKpi==GetSignalHistoIndex(iPtBin)) partDs = (AliAODMCParticle*)arrayMC->At(labDs);
		  if(labDplus >= 0) partDs = (AliAODMCParticle*)arrayMC->At(labDplus);
		  Int_t orig = AliVertexingHFUtils::CheckOrigin(arrayMC,partDs,kTRUE);
		  if(orig==4) {
		    candType = 1.5;
		  }
		  if(orig==5) {
		    candType = 2.5;
		  }
		}
	      }
	    }
	    if(isK0starKKpi) fMassHistK0st[indexMCKKpi]->Fill(invMass,weightKKpi);
	  }
	  if(isPhiKKpi && fFillImpParSparse) {
            Double_t impParxy = d->ImpParXY()*10000.;
            if(!fReadMC) {
	      Double_t array4ImpPar[3] = {ptCand,impParxy,invMass};
	      fImpParSparse->Fill(array4ImpPar);
            }
            else {
	      Double_t array4ImpPar[3] = {ptCand,impParxy,candType};
	      fImpParSparseMC->Fill(array4ImpPar);
            }
	  }
	}
	if(ispiKK){
	  if(fDoRotBkg && TMath::Abs(massKK-massPhi)<=fMaxDeltaPhiMass4Rot)GenerateRotBkg(d,2,iPtBin);
              
	  invMass=d->InvMassDspiKK();
	  fMassHist[index]->Fill(invMass,weightpiKK);
	  fPtVsMass->Fill(invMass,ptCand,weightpiKK);
              
	  if(fDoBkgPhiSB && 0.010<TMath::Abs(massKK-massPhi)<0.030) {
	    if(massKK<massPhi)fMassLSBkgHistPhi[iPtBin]->Fill(invMass);
	    else fMassRSBkgHistPhi[iPtBin]->Fill(invMass);
	  }
              
	  if(isPhipiKK){
	    fMassHistPhi[index]->Fill(invMass,weightpiKK);
	    fPtVsMassPhi->Fill(invMass,ptCand,weightpiKK);
	  }
	  if(isK0starpiKK){
	    fMassHistK0st[index]->Fill(invMass,weightpiKK);
	    fPtVsMassK0st->Fill(invMass,ptCand,weightpiKK);
	  }
	  if(fReadMC  && indexMCpiKK!=-1){
	    fMassHist[indexMCpiKK]->Fill(invMass,weightpiKK);
	    if(isPhipiKK) {
	      fMassHistPhi[indexMCpiKK]->Fill(invMass,weightpiKK);
	      if(fFillSparse) {
		if(indexMCpiKK==GetSignalHistoIndex(iPtBin) || labDplus >= 0) {
		  AliAODMCParticle *partDs;
		  if(indexMCpiKK==GetSignalHistoIndex(iPtBin)) partDs = (AliAODMCParticle*)arrayMC->At(labDs);
		  if(labDplus >= 0) partDs = (AliAODMCParticle*)arrayMC->At(labDplus);
		  Int_t orig = AliVertexingHFUtils::CheckOrigin(arrayMC,partDs,kTRUE);
		  if(orig==4) {
		    candType = 1.5;
		  }
		  if(orig==5) {
		    candType = 2.5;
		  }
		}
	      }
	    }
	    if(isK0starpiKK) fMassHistK0st[indexMCpiKK]->Fill(invMass,weightpiKK);
	  }
	  if(isPhipiKK && fFillImpParSparse) {
	    Double_t impParxy = d->ImpParXY()*10000.;
	    if(!fReadMC) {
	      Double_t array4ImpPar[3] = {ptCand,impParxy,invMass};
	      fImpParSparse->Fill(array4ImpPar);
	    }
	    else {
	      Double_t array4ImpPar[3] = {ptCand,impParxy,candType};
	      fImpParSparseMC->Fill(array4ImpPar);
	    }
	  }
	}
          
          
	///////////////////// CODE FOR NSPARSE /////////////////////////
          
	const Int_t nProng = 3;
	Double_t deltaMassKK=999.;
	Double_t dlen=d->DecayLength();
	Double_t dlenxy=d->DecayLengthXY();
	Double_t normdl=d->NormalizedDecayLength();
	Double_t normdlxy=d->NormalizedDecayLengthXY();
	Double_t cosp=d->CosPointingAngle();
	Double_t cospxy=d->CosPointingAngleXY();
	Double_t pt0=d->PtProng(0);
	Double_t pt1=d->PtProng(1);
	Double_t pt2=d->PtProng(2);
	Double_t sigvert=d->GetSigmaVert();
	Double_t cosPiDs=-99.;
	Double_t cosPiKPhi=-99.;
	Double_t normIP;                     //to store the maximum topomatic var. among the 3 prongs
	Double_t normIPprong[nProng];        //to store IP of k,k,pi

	Double_t ptWeight = 1.;
	if(fFillSparse) {
	  if (fUseWeight && fHistoPtWeight){
	    AliDebug(2,"Using Histogram as Pt weight function");
	    ptWeight = GetPtWeightFromHistogram(ptCand);
	  }
	  if(isPhiKKpi) {
	    Double_t tmpNormIP[nProng];
                        
	    invMass     = d->InvMassDsKKpi();
	    massKK      = d->InvMass2Prongs(0,1,321,321);
	    deltaMassKK = TMath::Abs(massKK-massPhi);
	    cosPiDs     = d->CosPiDsLabFrameKKpi();
	    cosPiKPhi   = d->CosPiKPhiRFrameKKpi();
	    cosPiKPhi   = TMath::Abs(cosPiKPhi*cosPiKPhi*cosPiKPhi);
                        
	    for(Int_t ip=0; ip<nProng; ip++) {
	      Double_t diffIP, errdiffIP;
	      d->Getd0MeasMinusExpProng(ip,aod->GetMagneticField(),diffIP,errdiffIP);
	      tmpNormIP[ip] = diffIP/errdiffIP;
	      if(ip==0) normIP = tmpNormIP[ip];
	      else if(TMath::Abs(tmpNormIP[ip])>TMath::Abs(normIP)) normIP = tmpNormIP[ip];
	    }
	    normIPprong[0] = tmpNormIP[0];
	    normIPprong[1] = tmpNormIP[1];
	    normIPprong[2] = tmpNormIP[2];
                        
	    Double_t var4nSparse[knVarForSparse] = {invMass,ptCand,deltaMassKK*1000,dlen*1000,dlenxy*1000,normdlxy,cosp*100,cospxy*100,
						    sigvert*1000,cosPiDs*10,cosPiKPhi*10,TMath::Abs(normIP)};
          
	    if(!fReadMC) {
	      fnSparse->Fill(var4nSparse);
	    }
	    else {
	      if(indexMCKKpi==GetSignalHistoIndex(iPtBin)) {
		if(candType==1.5) fnSparseMC[2]->Fill(var4nSparse,ptWeight);
		if(candType==2.5) fnSparseMC[3]->Fill(var4nSparse,ptWeight);
	      }
	      else if(fFillSparseDplus && labDplus>=0 && pdgCode0==321) {
		if(candType==1.5) fnSparseMCDplus[2]->Fill(var4nSparse,ptWeight);
		if(candType==2.5) fnSparseMCDplus[3]->Fill(var4nSparse,ptWeight);
	      }
	    }
	  }
	  if(isPhipiKK) {
	    Double_t tmpNormIP[nProng];
                        
	    invMass     = d->InvMassDspiKK();
	    massKK      = d->InvMass2Prongs(1,2,321,321);
	    deltaMassKK = TMath::Abs(massKK-massPhi);
	    cosPiDs     = d->CosPiDsLabFramepiKK();
	    cosPiKPhi   = d->CosPiKPhiRFramepiKK();
	    cosPiKPhi   = TMath::Abs(cosPiKPhi*cosPiKPhi*cosPiKPhi);
                        
	    for(Int_t ip=0; ip<nProng; ip++) {
	      Double_t diffIP, errdiffIP;
	      d->Getd0MeasMinusExpProng(ip,aod->GetMagneticField(),diffIP,errdiffIP);
	      tmpNormIP[ip] = diffIP/errdiffIP;
	      if(ip==0) normIP = tmpNormIP[ip];
	      else if(TMath::Abs(tmpNormIP[ip])>TMath::Abs(normIP)) normIP = tmpNormIP[ip];
	    }
                        
	    normIPprong[0] = tmpNormIP[2];
	    normIPprong[1] = tmpNormIP[1];
	    normIPprong[2] = tmpNormIP[0];
                        
	    Double_t var4nSparse[knVarForSparse] = {invMass,ptCand,deltaMassKK*1000,dlen*1000,dlenxy*1000,normdlxy,cosp*100,cospxy*100,
						    sigvert*1000,cosPiDs*10,cosPiKPhi*10,TMath::Abs(normIP)};
          
	    if(!fReadMC) {
	      fnSparse->Fill(var4nSparse);
	    }
	    else {
	      if(indexMCpiKK==GetSignalHistoIndex(iPtBin)) {
		if(candType==1.5) fnSparseMC[2]->Fill(var4nSparse,ptWeight);
		if(candType==2.5) fnSparseMC[3]->Fill(var4nSparse,ptWeight);
	      }
	      else if(fFillSparseDplus && labDplus>=0 && pdgCode0==211) {
		if(candType==1.5) fnSparseMCDplus[2]->Fill(var4nSparse,ptWeight);
		if(candType==2.5) fnSparseMCDplus[3]->Fill(var4nSparse,ptWeight);
	      }
	    }
	  }
                    
	  if(fReadMC && (isPhiKKpi || isPhiKKpi)) {
	    Double_t var[6] = {ptCand,normIP,normIPprong[0],normIPprong[1],normIPprong[2],candType};
	    fnSparseIP->Fill(var);
	  }
	}
	////////////////////////////////////////////////////////////////
                
	if(fDoCutVarHistos){
	  Double_t sumD02=d->Getd0Prong(0)*d->Getd0Prong(0)+d->Getd0Prong(1)*d->Getd0Prong(1)+d->Getd0Prong(2)*d->Getd0Prong(2);
	  Double_t dca=d->GetDCA();
	  Double_t ptmax=0;
              
	  for(Int_t i=0;i<3;i++){
	    if(d->PtProng(i)>ptmax)ptmax=d->PtProng(i);
	  }
	  fCosPHist[index]->Fill(cosp);
	  fCosPxyHist[index]->Fill(cospxy);
	  fDLenHist[index]->Fill(dlen);
	  fDLenxyHist[index]->Fill(dlenxy);
	  fNDLenxyHist[index]->Fill(normdlxy);
	  fSigVertHist[index]->Fill(sigvert);
	  fSumd02Hist[index]->Fill(sumD02);
	  fPtMaxHist[index]->Fill(ptmax);
	  fDCAHist[index]->Fill(dca);
	  fNormIPHist[index]->Fill(normIP);
	  fCosPiDsHist[index]->Fill(cosPiDs);
	  fCosPiKPhiHist[index]->Fill(cosPiKPhi);
	  fPtProng0Hist[index]->Fill(pt0);
	  fPtProng1Hist[index]->Fill(pt1);
	  fPtProng2Hist[index]->Fill(pt2);
	  if(!fReadMC) {
	    fCosPHist3D->Fill(invMass,ptCand,cosp);
	    fCosPxyHist3D->Fill(invMass,ptCand,cospxy);
	    fDLenHist3D->Fill(invMass,ptCand,dlen);
	    fDLenxyHist3D->Fill(invMass,ptCand,dlenxy);
	    fNDLenxyHist3D->Fill(invMass,ptCand,normdlxy);
	    fSigVertHist3D->Fill(invMass,ptCand,sigvert);
	    fDCAHist3D->Fill(invMass,ptCand,dca);
	    fNormIPHist3D->Fill(invMass,ptCand,normIP);
	    fCosPiDsHist3D->Fill(invMass,ptCand,cosPiDs);
	    fCosPiKPhiHist3D->Fill(invMass,ptCand,cosPiKPhi);
	    fPtProng0Hist3D->Fill(invMass,ptCand,pt0);
	    fPtProng1Hist3D->Fill(invMass,ptCand,pt1);
	    fPtProng2Hist3D->Fill(invMass,ptCand,pt2);
	  }
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
	      fCosPxyHist[indexMCKKpi]->Fill(cospxy);
	      fDLenHist[indexMCKKpi]->Fill(dlen);
	      fDLenxyHist[indexMCKKpi]->Fill(dlenxy);
	      fNDLenxyHist[indexMCKKpi]->Fill(normdlxy);
	      fSigVertHist[indexMCKKpi]->Fill(sigvert);
	      fSumd02Hist[indexMCKKpi]->Fill(sumD02);
	      fPtMaxHist[indexMCKKpi]->Fill(ptmax);
	      fPtCandHist[indexMCKKpi]->Fill(ptCand);
	      fDCAHist[indexMCKKpi]->Fill(dca);
	      fNormIPHist[indexMCKKpi]->Fill(normIP);
	      fCosPiDsHist[indexMCKKpi]->Fill(cosPiDs);
	      fCosPiKPhiHist[indexMCKKpi]->Fill(cosPiKPhi);
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
	      fCosPxyHist[indexMCpiKK]->Fill(cospxy);
	      fDLenHist[indexMCpiKK]->Fill(dlen);
	      fDLenxyHist[indexMCpiKK]->Fill(dlenxy);
	      fNDLenxyHist[indexMCpiKK]->Fill(normdlxy);
	      fSigVertHist[indexMCpiKK]->Fill(sigvert);
	      fSumd02Hist[indexMCpiKK]->Fill(sumD02);
	      fPtMaxHist[indexMCpiKK]->Fill(ptmax);
	      fPtCandHist[indexMCpiKK]->Fill(ptCand);
	      fDCAHist[indexMCpiKK]->Fill(dca);
	      fNormIPHist[indexMCpiKK]->Fill(normIP);
	      fCosPiDsHist[indexMCpiKK]->Fill(cosPiDs);
	      fCosPiKPhiHist[indexMCpiKK]->Fill(cosPiKPhi);
	      fPtProng0Hist[indexMCpiKK]->Fill(pt0);
	      fPtProng1Hist[indexMCpiKK]->Fill(pt1);
	      fPtProng2Hist[indexMCpiKK]->Fill(pt2);
	    }
	  }
	}
                
	Float_t tmp[37];
	if(fFillNtuple>0){
                    
	  if ((fFillNtuple==1 && (isPhiKKpi || isPhipiKK)) || (fFillNtuple==2 && (isK0starKKpi || isK0starpiKK)) || (fFillNtuple==3 && (isKKpi || ispiKK))){
                        
	    AliAODTrack *track0=(AliAODTrack*)d->GetDaughter(0);
	    AliAODTrack *track1=(AliAODTrack*)d->GetDaughter(1);
	    AliAODTrack *track2=(AliAODTrack*)d->GetDaughter(2);
                        
	    UInt_t bitMapPIDTrack0=fAnalysisCuts->GetPIDTrackTPCTOFBitMap(track0);
	    UInt_t bitMapPIDTrack1=fAnalysisCuts->GetPIDTrackTPCTOFBitMap(track1);
	    UInt_t bitMapPIDTrack2=fAnalysisCuts->GetPIDTrackTPCTOFBitMap(track2);
                        
	    tmp[0]=Float_t(labDs);
	    if(fReadMC && fWriteOnlySignal) tmp[0]=Float_t(isMCSignal);
	    tmp[1]=Float_t(retCodeAnalysisCuts);
	    tmp[2]=Float_t(pdgCode0);
	    tmp[3]=d->PtProng(0);
	    tmp[4]=d->PtProng(1);
	    tmp[5]=d->PtProng(2);
	    tmp[6]=d->Pt();
	    tmp[7]=d->PProng(0);
	    tmp[8]=d->PProng(1);
	    tmp[9]=d->PProng(2);
	    tmp[10]=Int_t(bitMapPIDTrack0);
	    tmp[11]=Int_t(bitMapPIDTrack1);
	    tmp[12]=Int_t(bitMapPIDTrack2);
	    tmp[13]=d->CosPointingAngle();
	    tmp[14]=d->CosPointingAngleXY();
	    tmp[15]=d->DecayLength();
	    tmp[16]=d->DecayLengthXY();
	    tmp[17]=d->NormalizedDecayLength();
	    tmp[18]=d->NormalizedDecayLengthXY();
	    tmp[19]=d->InvMassDsKKpi();
	    tmp[20]=d->InvMassDspiKK();
	    tmp[21]=d->GetSigmaVert();
	    tmp[22]=d->Getd0Prong(0);
	    tmp[23]=d->Getd0Prong(1);
	    tmp[24]=d->Getd0Prong(2);
	    tmp[25]=d->GetDCA();
	    tmp[26]=d->Getd0Prong(0)*d->Getd0Prong(0)+d->Getd0Prong(1)*d->Getd0Prong(1)+d->Getd0Prong(2)*d->Getd0Prong(2);
	    tmp[27]=d->InvMass2Prongs(0,1,321,321);
	    tmp[28]=d->InvMass2Prongs(1,2,321,321);
	    tmp[29]=d->InvMass2Prongs(1,2,321,211);
	    tmp[30]=d->InvMass2Prongs(0,1,211,321);
	    tmp[31]=d->CosPiDsLabFrameKKpi();
	    tmp[32]=d->CosPiDsLabFramepiKK();
	    tmp[33]=d->CosPiKPhiRFrameKKpi();
	    tmp[34]=d->CosPiKPhiRFramepiKK();
	    tmp[35]=(Float_t)(centrality);
	    tmp[36]=(Float_t)(runNumber);
                        
	    if(fReadMC && fWriteOnlySignal){
	      if(isMCSignal>=0) fNtupleDs->Fill(tmp);
	    }else{
	      fNtupleDs->Fill(tmp);
	    }
	    PostData(4,fNtupleDs);
	  }
	}
      } //if(retCodeAnalysisCuts
    } // if(isFidAcc)
        
    if(unsetvtx) d->UnsetOwnPrimaryVtx();
    if(recVtx)fAnalysisCuts->CleanOwnPrimaryVtx(d,aod,origownvtx);
  }
    
  fCounter->StoreCandidates(aod,nFiltered,kTRUE);
  fCounter->StoreCandidates(aod,nSelected,kFALSE);
    
  delete vHF;
    
  PostData(1,fOutput); 
  PostData(3,fCounter);    
    
  return;
}

//_________________________________________________________________

void AliAnalysisTaskSEDs::Terminate(Option_t */*option*/)
{
  /// Terminate analysis
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

//_________________________________________________________________
void AliAnalysisTaskSEDs::FillMCGenAccHistos(TClonesArray *arrayMC, AliAODMCHeader *mcHeader){
  /// Fill MC histos for cuts study at GenLimAccStep and AccStep
    
  Int_t nProng = 3;
  Double_t zMCVertex = mcHeader->GetVtxZ(); //vertex MC
  if(TMath::Abs(zMCVertex) <= fAnalysisCuts->GetMaxVtxZ()) {
    for(Int_t iPart=0; iPart<arrayMC->GetEntriesFast(); iPart++){
            
      AliAODMCParticle* mcPart = dynamic_cast<AliAODMCParticle*>(arrayMC->At(iPart));
            
      if(TMath::Abs(mcPart->GetPdgCode()) == 431) {
	Int_t orig = AliVertexingHFUtils::CheckOrigin(arrayMC,mcPart,kTRUE);//Prompt = 4, FeedDown = 5
                
	Int_t  deca        = 0;
	Bool_t isGoodDecay = kFALSE;
	Int_t  labDau[3]   = {-1,-1,-1};
	Bool_t isFidAcc    = kFALSE;
	Bool_t isDaugInAcc = kFALSE;
                
	deca = AliVertexingHFUtils::CheckDsDecay(arrayMC,mcPart,labDau);
	if(deca == 1) isGoodDecay=kTRUE; // == 1 -> Phi pi -> kkpi
                
	if(labDau[0]==-1) continue; //protection against unfilled array of labels
                
	if(isGoodDecay) {
	  Double_t pt = mcPart->Pt();
	  Double_t rapid = mcPart->Y();
	  isFidAcc = fAnalysisCuts->IsInFiducialAcceptance(pt,rapid);
	  isDaugInAcc = CheckDaugAcc(arrayMC,nProng,labDau);
        
	  Double_t ptWeight = 1.;
	  if (fUseWeight && fHistoPtWeight){
            AliDebug(2,"Using Histogram as Pt weight function");
            ptWeight = GetPtWeightFromHistogram(pt);
	  }

	  if(isFidAcc) {
	    Double_t var4nSparseAcc[2] = {pt,rapid*10};
	    if(isDaugInAcc) {
	      if(orig==4) fnSparseMC[0]->Fill(var4nSparseAcc,ptWeight);
	      if(orig==5) fnSparseMC[1]->Fill(var4nSparseAcc,ptWeight);
	    }
	  }
	}
      }
            
      if(fFillSparseDplus && TMath::Abs(mcPart->GetPdgCode()) == 411) {
	Int_t orig = AliVertexingHFUtils::CheckOrigin(arrayMC,mcPart,kTRUE);//Prompt = 4, FeedDown = 5
                
	Int_t  deca        = 0;
	Bool_t isGoodDecay = kFALSE;
	Int_t  labDau[3]   = {-1,-1,-1};
	Bool_t isFidAcc    = kFALSE;
	Bool_t isDaugInAcc = kFALSE;
                
	deca = AliVertexingHFUtils::CheckDplusKKpiDecay(arrayMC,mcPart,labDau);
	if(deca == 1) isGoodDecay=kTRUE; // == 1 -> Phi pi -> kkpi
                
	if(labDau[0]==-1) continue; //protection against unfilled array of labels
                
	if(isGoodDecay) {
	  Double_t pt = mcPart->Pt();
	  Double_t rapid = mcPart->Y();
	  isFidAcc = fAnalysisCuts->IsInFiducialAcceptance(pt,rapid);
	  isDaugInAcc = CheckDaugAcc(arrayMC,nProng,labDau);
     
	  Double_t ptWeight = 1.;
	  if (fUseWeight && fHistoPtWeight){
            AliDebug(2,"Using Histogram as Pt weight function");
            ptWeight = GetPtWeightFromHistogram(pt);
	  }

	  if(isFidAcc) {
	    Double_t var4nSparseAcc[2] = {pt,rapid*10};
	    if(isDaugInAcc) {
	      if(orig==4) fnSparseMCDplus[0]->Fill(var4nSparseAcc,ptWeight);
	      if(orig==5) fnSparseMCDplus[1]->Fill(var4nSparseAcc,ptWeight);
	    }
	  }
	}
      }
    }
  }
}
//_________________________________________________________________
Bool_t AliAnalysisTaskSEDs::CheckDaugAcc(TClonesArray* arrayMC,Int_t nProng, Int_t *labDau){
  /// check if the decay products are in the good eta and pt range
    
  for (Int_t iProng = 0; iProng<nProng; iProng++){
    AliAODMCParticle* mcPartDaughter=dynamic_cast<AliAODMCParticle*>(arrayMC->At(labDau[iProng]));
    if(!mcPartDaughter) {
      return kFALSE;
    }
    Double_t eta = mcPartDaughter->Eta();
    Double_t pt = mcPartDaughter->Pt();
    if (TMath::Abs(eta) > 0.9 || pt < 0.1) {
      return kFALSE;
    }
  }
  return kTRUE;
}

//_________________________________________________________________

void AliAnalysisTaskSEDs::GenerateRotBkg(AliAODRecoDecayHF3Prong *d, Int_t dec, Int_t iPtBin) {
    
  const Int_t nprongs = 3;
  Double_t PxProng[nprongs], PyProng[nprongs], PtProng[nprongs], PzProng[nprongs], P2Prong[nprongs], mProng[nprongs];
  Double_t Px, Py, Pz, P2;
  UInt_t pdg[3]={321,321,211};
  int idPion = 2;
  if(dec==2) {
    pdg[0]=211;
    pdg[2]=321;
    idPion = 0;
  }
    
  for (Int_t ip=0; ip<nprongs; ip++) {
    PxProng[ip] = d->PxProng(ip);
    PyProng[ip] = d->PxProng(ip);
    PtProng[ip] = d->PtProng(ip);
    PzProng[ip] = d->PzProng(ip);
    P2Prong[ip] = d->P2Prong(ip);
    mProng[ip]  = TDatabasePDG::Instance()->GetParticle(pdg[ip])->Mass();
  }
    
  for(Int_t i=0; i<9; i++) { //9 rotations implemented for the pion track around Pi
    Px = 0.;
    Py = 0.;
    Pz = 0.;
      
    Double_t phirot=TMath::Pi()*(5/6.+1/27.*i);

    PxProng[idPion] = PxProng[idPion]*TMath::Cos(phirot)-PyProng[idPion]*TMath::Sin(phirot);
    PyProng[idPion] = PxProng[idPion]*TMath::Sin(phirot)+PyProng[idPion]*TMath::Cos(phirot);

    for (Int_t j=0; j<nprongs; j++) {
      Px += PxProng[j];
      Py += PyProng[j];
      Pz += PzProng[j];
    }
    P2 = Px*Px + Py*Py + Pz*Pz;
        
    Double_t energysum = 0.;
    for(Int_t j=0; j<nprongs; j++) {
      energysum += TMath::Sqrt(mProng[j]*mProng[j]+P2Prong[j]);
    }
    Double_t mass = TMath::Sqrt(energysum*energysum-P2);
    if(fminMass<=mass<fmaxMass) fMassRotBkgHistPhi[iPtBin]->Fill(mass);
  }
}
//_________________________________________________________________________
void AliAnalysisTaskSEDs::SetPtWeightsFromFONLL5anddataoverLHC16i2a(){
  // weight function from the ratio of the LHC16i2a MC
  // 1.5-14 GeV/c using data and 1-1.5, 14-50 GeV/c using FONLL calculations
    
  if(fHistoPtWeight) delete fHistoPtWeight;
  fHistoPtWeight = new TH1F("histoWeight","histoWeight",500,0.,50.);
  Float_t binc[500]={ 1.695705, 1.743693, 1.790289, 1.835410, 1.878981, 1.920938, 1.961223, 1.999787, 2.036589, 2.071597, 2.104784, 2.136132, 2.165629, 2.193270, 2.219057, 2.174545, 2.064698, 1.959489, 1.858770, 1.762396, 1.670224, 1.582115, 1.497931, 1.417541, 1.340814, 1.267622, 1.197842, 1.131352, 1.068033, 1.007770, 0.950450, 0.895963, 0.844202, 0.795062, 0.748441, 0.704241, 0.662363, 0.622715, 0.585204, 0.549742, 0.516242, 0.484620, 0.454795, 0.426686, 0.400217, 0.375314, 0.351903, 0.329915, 0.309281, 0.289936, 0.271816, 0.254860, 0.239007, 0.224201, 0.210386, 0.197508, 0.185516, 0.174360, 0.163992, 0.154366, 0.145438, 0.137166, 0.129508, 0.122426, 0.115882, 0.109840, 0.104266, 0.099128, 0.094395, 0.090036, 0.086023, 0.082331, 0.078933, 0.075805, 0.072925, 0.070271, 0.067823, 0.065562, 0.063471, 0.061532, 0.059730, 0.058051, 0.056481, 0.055007, 0.053619, 0.052306, 0.051059, 0.049867, 0.048725, 0.047624, 0.046558, 0.045522, 0.044511, 0.043521, 0.042548, 0.041590, 0.040643, 0.039706, 0.038778, 0.037857, 0.036944, 0.036039, 0.035141, 0.034251, 0.033370, 0.032500, 0.031641, 0.030796, 0.029966, 0.029153, 0.028359, 0.027587, 0.026837, 0.026113, 0.025416, 0.024748, 0.024111, 0.023507, 0.022937, 0.022402, 0.021904, 0.021443, 0.021020, 0.020634, 0.020286, 0.019974, 0.019698, 0.019455, 0.019244, 0.019062, 0.018905, 0.018770, 0.018652, 0.018545, 0.018444, 0.018342, 0.018231, 0.018102, 0.017947, 0.017755, 0.017536, 0.017327, 0.017120, 0.016915, 0.016713, 0.016514, 0.016317, 0.016122, 0.015929, 0.015739, 0.015551, 0.015366, 0.015182, 0.015001, 0.014822, 0.014645, 0.014470, 0.014297, 0.014127, 0.013958, 0.013791, 0.013627, 0.013464, 0.013303, 0.013145, 0.012988, 0.012833, 0.012679, 0.012528, 0.012378, 0.012231, 0.012085, 0.011940, 0.011798, 0.011657, 0.011518, 0.011380, 0.011244, 0.011110, 0.010978, 0.010846, 0.010717, 0.010589, 0.010463, 0.010338, 0.010214, 0.010092, 0.009972, 0.009853, 0.009735, 0.009619, 0.009504, 0.009391, 0.009279, 0.009168, 0.009058, 0.008950, 0.008843, 0.008738, 0.008633, 0.008530, 0.008429, 0.008328, 0.008229, 0.008130, 0.008033, 0.007937, 0.007843, 0.007749, 0.007656, 0.007565, 0.007475, 0.007385, 0.007297, 0.007210, 0.007124, 0.007039, 0.006955, 0.006872, 0.006790, 0.006709, 0.006629, 0.006550, 0.006471, 0.006394, 0.006318, 0.006242, 0.006168, 0.006094, 0.006022, 0.005950, 0.005879, 0.005808, 0.005739, 0.005671, 0.005603, 0.005536, 0.005470, 0.005405, 0.005340, 0.005276, 0.005213, 0.005151, 0.005090, 0.005029, 0.004969, 0.004909, 0.004851, 0.004793, 0.004736, 0.004679, 0.004623, 0.004568, 0.004514, 0.004460, 0.004406, 0.004354, 0.004302, 0.004251, 0.004200, 0.004150, 0.004100, 0.004051, 0.004003, 0.003955, 0.003908, 0.003861, 0.003815, 0.003770, 0.003725, 0.003680, 0.003636, 0.003593, 0.003550, 0.003507, 0.003466, 0.003424, 0.003383, 0.003343, 0.003303, 0.003264, 0.003225, 0.003186, 0.003148, 0.003110, 0.003073, 0.003037, 0.003000, 0.002965, 0.002929, 0.002894, 0.002860, 0.002826, 0.002792, 0.002758, 0.002726, 0.002693, 0.002661, 0.002629, 0.002598, 0.002567, 0.002536, 0.002506, 0.002476, 0.002446, 0.002417, 0.002388, 0.002360, 0.002332, 0.002304, 0.002276, 0.002249, 0.002222, 0.002196, 0.002169, 0.002144, 0.002118, 0.002093, 0.002068, 0.002043, 0.002019, 0.001995, 0.001971, 0.001947, 0.001924, 0.001901, 0.001878, 0.001856, 0.001834, 0.001812, 0.001790, 0.001769, 0.001748, 0.001727, 0.001706, 0.001686, 0.001666, 0.001646, 0.001626, 0.001607, 0.001588, 0.001569, 0.001550, 0.001531, 0.001513, 0.001495, 0.001477, 0.001460, 0.001442, 0.001425, 0.001408, 0.001391, 0.001374, 0.001358, 0.001342, 0.001326, 0.001310, 0.001294, 0.001279, 0.001264, 0.001249, 0.001234, 0.001219, 0.001204, 0.001190, 0.001176, 0.001162, 0.001148, 0.001134, 0.001121, 0.001107, 0.001094, 0.001081, 0.001068, 0.001055, 0.001043, 0.001030, 0.001018, 0.001006, 0.000994, 0.000982, 0.000970, 0.000959, 0.000947, 0.000936, 0.000925, 0.000914, 0.000903, 0.000892, 0.000881, 0.000871, 0.000860, 0.000850, 0.000840, 0.000830, 0.000820, 0.000810, 0.000801, 0.000791, 0.000782, 0.000772, 0.000763, 0.000754, 0.000745, 0.000736, 0.000727, 0.000719, 0.000710, 0.000702, 0.000693, 0.000685, 0.000677, 0.000669, 0.000661, 0.000653, 0.000645, 0.000637, 0.000630, 0.000622, 0.000615, 0.000607, 0.000600, 0.000593, 0.000586, 0.000579, 0.000572, 0.000565, 0.000558, 0.000552, 0.000545, 0.000539, 0.000532, 0.000526, 0.000520, 0.000513, 0.000507, 0.000501, 0.000495, 0.000489, 0.000483, 0.000478, 0.000472, 0.000466, 0.000461, 0.000455, 0.000450, 0.000444, 0.000439, 0.000434, 0.000429, 0.000424, 0.000419, 0.000414, 0.000409, 0.000404, 0.000399, 0.000394, 0.000389, 0.000385, 0.000380, 0.000376, 0.000371, 0.000367, 0.000362, 0.000358, 0.000354, 0.000350, 0.000345, 0.000341, 0.000337, 0.000333, 0.000329, 0.000325, 0.000321, 0.000318, 0.000314, 0.000310, 0.000306, 0.000303, 0.000299, 0.000295, 0.000292, 0.000288, 0.000285, 0.000282, 0.000278, 0.000275, 0.000272, 0.000268, 0.000265, 0.000262, 0.000259, 0.000256, 0.000253, 0.000250, 0.000247, 0.000244, 0.000241, 0.000238, 0.000235};
  for(Int_t i=0; i<500; i++){
    fHistoPtWeight->SetBinContent(i+1,binc[i]);
  }
  //SetWeightHistogram();
  fUseWeight=kTRUE;
}

//_________________________________________________________________________
void AliAnalysisTaskSEDs::SetPtWeightsFromFONLL5overLHC16i2abc(){
  // weight function from the ratio of the LHC16i2a+b+c MC
  // and FONLL calculations for pp data
    
  if(fHistoPtWeight) delete fHistoPtWeight;
  fHistoPtWeight = new TH1F("histoWeight","histoWeight",400,0.,40.);
  Float_t binc[400]={1.118416, 1.003458, 0.935514, 0.907222, 0.904359, 0.913668, 0.933906, 0.963898, 0.996388, 1.031708, 1.066404, 1.099683, 1.125805, 1.145181, 1.165910, 1.181905, 1.193425, 1.203891, 1.204726, 1.209411, 1.209943, 1.204763, 1.205291, 1.198912, 1.197390, 1.182005, 1.184194, 1.175994, 1.167881, 1.158348, 1.147190, 1.139833, 1.126940, 1.123322, 1.108389, 1.102199, 1.089464, 1.075874, 1.061964, 1.051429, 1.038113, 1.026668, 1.011441, 0.998567, 0.987658, 0.972434, 0.950068, 0.940758, 0.916880, 0.911931, 0.894512, 0.878691, 0.860589, 0.848025, 0.830774, 0.819399, 0.801134, 0.775276, 0.766382, 0.750495, 0.736935, 0.717529, 0.702637, 0.689152, 0.671334, 0.652030, 0.635696, 0.621365, 0.608362, 0.599019, 0.576024, 0.562136, 0.550938, 0.533587, 0.516410, 0.509744, 0.501655, 0.487402, 0.476469, 0.463762, 0.445979, 0.438088, 0.422214, 0.417467, 0.404357, 0.391450, 0.379996, 0.371201, 0.361497, 0.352912, 0.343189, 0.329183, 0.327662, 0.310783, 0.304525, 0.301007, 0.293306, 0.278332, 0.274419, 0.267361, 0.261459, 0.255514, 0.249293, 0.241129, 0.237600, 0.231343, 0.221982, 0.216872, 0.211094, 0.206954, 0.202333, 0.196572, 0.193274, 0.188240, 0.181817, 0.178364, 0.173614, 0.167135, 0.166055, 0.163423, 0.156557, 0.155821, 0.151985, 0.144909, 0.145062, 0.139720, 0.138873, 0.131892, 0.129969, 0.126509, 0.126978, 0.120451, 0.117661, 0.116300, 0.115604, 0.112215, 0.109237, 0.107720, 0.106419, 0.102050, 0.102777, 0.097406, 0.098447, 0.095964, 0.093868, 0.092430, 0.089329, 0.088249, 0.085881, 0.084417, 0.085498, 0.082444, 0.079151, 0.079565, 0.077811, 0.077293, 0.075218, 0.072445, 0.073054, 0.071545, 0.070279, 0.068046, 0.067854, 0.068092, 0.065378, 0.064405, 0.062060, 0.063391, 0.061718, 0.059616, 0.058913, 0.058895, 0.058311, 0.056320, 0.056527, 0.055349, 0.053701, 0.054735, 0.052264, 0.051277, 0.051554, 0.050545, 0.048995, 0.049507, 0.048466, 0.048156, 0.046809, 0.047600, 0.046078, 0.044801, 0.044113, 0.043700, 0.043530, 0.043396, 0.042556, 0.041048, 0.041657, 0.040394, 0.041314, 0.040720, 0.039656, 0.038478, 0.039276, 0.038777, 0.037730, 0.036918, 0.036466, 0.035827, 0.035285, 0.035963, 0.034371, 0.034757, 0.033205, 0.033666, 0.033266, 0.032583, 0.033570, 0.032102, 0.032107, 0.031464, 0.032160, 0.030091, 0.030564, 0.029464, 0.029613, 0.029626, 0.029512, 0.029324, 0.028607, 0.027628, 0.027251, 0.027072, 0.027077, 0.026724, 0.026961, 0.026303, 0.026237, 0.025454, 0.025133, 0.025365, 0.026014, 0.024807, 0.023901, 0.023459, 0.023405, 0.023654, 0.023981, 0.023675, 0.022493, 0.022781, 0.021801, 0.021704, 0.022372, 0.021189, 0.020681, 0.020779, 0.021324, 0.020558, 0.020901, 0.020586, 0.020808, 0.019276, 0.019516, 0.019706, 0.018935, 0.018632, 0.018516, 0.019187, 0.018916, 0.018039, 0.018208, 0.018045, 0.017628, 0.017916, 0.017711, 0.017838, 0.017222, 0.016565, 0.015733, 0.016264, 0.015826, 0.016090, 0.016622, 0.015802, 0.016621, 0.015441, 0.015309, 0.014860, 0.014935, 0.014968, 0.014443, 0.014485, 0.015136, 0.014078, 0.014414, 0.013908, 0.014071, 0.014078, 0.013766, 0.013436, 0.013507, 0.013480, 0.013224, 0.013041, 0.013935, 0.012885, 0.012453, 0.012528, 0.012492, 0.012225, 0.012542, 0.012706, 0.012136, 0.011902, 0.011560, 0.011448, 0.011861, 0.011271, 0.011831, 0.011159, 0.011171, 0.010966, 0.011311, 0.011002, 0.011130, 0.010995, 0.010450, 0.010663, 0.010678, 0.010492, 0.009861, 0.010507, 0.009916, 0.010121, 0.010029, 0.010046, 0.009370, 0.009647, 0.010104, 0.009282, 0.009830, 0.009403, 0.009148, 0.009172, 0.008893, 0.009158, 0.009019, 0.008780, 0.008579, 0.009063, 0.008634, 0.008988, 0.008265, 0.008581, 0.008575, 0.008690, 0.008181, 0.008352, 0.008150, 0.008430, 0.008256, 0.008119, 0.008453, 0.008447, 0.008021, 0.007938, 0.008025, 0.007718, 0.008127, 0.007651, 0.007590, 0.007316, 0.007839, 0.007504, 0.007341, 0.007527, 0.007263, 0.007668, 0.007306, 0.007271, 0.006910, 0.007257, 0.007260, 0.006810, 0.006967, 0.006887, 0.006867, 0.007202, 0.006829, 0.006370, 0.006710, 0.006417, 0.006361, 0.006800, 0.006410, 0.006323, 0.006790, 0.006322, 0.006673, 0.006547};
  for(Int_t i=0; i<400; i++){
    fHistoPtWeight->SetBinContent(i+1,binc[i]);
  }
  //SetWeightHistogram();
  fUseWeight=kTRUE;
}

//_________________________________________________________________________
void AliAnalysisTaskSEDs::SetPtWeightsFromFONLL5andBAMPSoverLHC16i2abc(){
  // weight function from the ratio of the LHC16i2a+b+c MC
  // and FONLL calculations for pp data
  // corrected by the BAMPS Raa calculation for 30-50% CC
  if(fHistoPtWeight) delete fHistoPtWeight;
  fHistoPtWeight = new TH1F("histoWeight","histoWeight",400,0.,40.);
  Float_t binc[400]={2.166180, 1.866117, 1.667595, 1.547176, 1.486661, 1.457891, 1.426949, 1.399055, 1.383278, 1.349383, 1.317009, 1.282321, 1.234257, 1.181136, 1.136655, 1.087523, 1.037912, 0.993256, 0.944746, 0.900948, 0.865869, 0.827193, 0.794424, 0.757723, 0.733020, 0.700164, 0.682189, 0.659872, 0.637918, 0.615749, 0.593020, 0.574402, 0.556158, 0.542663, 0.525494, 0.516038, 0.503629, 0.490980, 0.479143, 0.469005, 0.457749, 0.447668, 0.436803, 0.427073, 0.418282, 0.407867, 0.395093, 0.387861, 0.374742, 0.369462, 0.360146, 0.351991, 0.342990, 0.336259, 0.327730, 0.322382, 0.314602, 0.303874, 0.299820, 0.293049, 0.287539, 0.280329, 0.274866, 0.269939, 0.263299, 0.256057, 0.249215, 0.242170, 0.235704, 0.230709, 0.220529, 0.213921, 0.208394, 0.202424, 0.196700, 0.194943, 0.192620, 0.187894, 0.184411, 0.180204, 0.172915, 0.169077, 0.162201, 0.159636, 0.153904, 0.148296, 0.143282, 0.139306, 0.135561, 0.132342, 0.128696, 0.123444, 0.122873, 0.116544, 0.114197, 0.112878, 0.110018, 0.104547, 0.103222, 0.100707, 0.098622, 0.096513, 0.094295, 0.091334, 0.090122, 0.087870, 0.084894, 0.083729, 0.082265, 0.081404, 0.080323, 0.078750, 0.078132, 0.076781, 0.074823, 0.074050, 0.072614, 0.070093, 0.069828, 0.068907, 0.066189, 0.066054, 0.064600, 0.061757, 0.061986, 0.059862, 0.059656, 0.056807, 0.055956, 0.054386, 0.054507, 0.051629, 0.050358, 0.049702, 0.049331, 0.047814, 0.046476, 0.045762, 0.045142, 0.043224, 0.043484, 0.041282, 0.041794, 0.040809, 0.039985, 0.039439, 0.038181, 0.037782, 0.036831, 0.036264, 0.036790, 0.035535, 0.034173, 0.034409, 0.033659, 0.033308, 0.032290, 0.030981, 0.031121, 0.030361, 0.029708, 0.028653, 0.028461, 0.028449, 0.027208, 0.026697, 0.025623, 0.026069, 0.025279, 0.024332, 0.024341, 0.024629, 0.024677, 0.024117, 0.024490, 0.024257, 0.023804, 0.024537, 0.023692, 0.023502, 0.023888, 0.023673, 0.023193, 0.023684, 0.023429, 0.023521, 0.023014, 0.023346, 0.022544, 0.021866, 0.021477, 0.021224, 0.021089, 0.020972, 0.020515, 0.019739, 0.019982, 0.019328, 0.019719, 0.019387, 0.018833, 0.018227, 0.018558, 0.018276, 0.017738, 0.017460, 0.017365, 0.017178, 0.017033, 0.017478, 0.016817, 0.017119, 0.016463, 0.016802, 0.016711, 0.016475, 0.017083, 0.016441, 0.016548, 0.016320, 0.016786, 0.015804, 0.016153, 0.015668, 0.015843, 0.015810, 0.015651, 0.015454, 0.014981, 0.014376, 0.014089, 0.013906, 0.013818, 0.013549, 0.013580, 0.013160, 0.013040, 0.012566, 0.012324, 0.012353, 0.012582, 0.011915, 0.011401, 0.011112, 0.011008, 0.011046, 0.011119, 0.010954, 0.010439, 0.010604, 0.010179, 0.010163, 0.010507, 0.009981, 0.009771, 0.009846, 0.010134, 0.009798, 0.009991, 0.009869, 0.010005, 0.009295, 0.009438, 0.009557, 0.009210, 0.009088, 0.009057, 0.009412, 0.009306, 0.008899, 0.009009, 0.008952, 0.008764, 0.008926, 0.008842, 0.008924, 0.008634, 0.008322, 0.007920, 0.008205, 0.008000, 0.008151, 0.008438, 0.008037, 0.008472, 0.007886, 0.007835, 0.007621, 0.007675, 0.007707, 0.007452, 0.007489, 0.007841, 0.007308, 0.007497, 0.007248, 0.007348, 0.007367, 0.007227, 0.007097, 0.007179, 0.007209, 0.007115, 0.007059, 0.007588, 0.007058, 0.006862, 0.006945, 0.006965, 0.006856, 0.007075, 0.007209, 0.006925, 0.006830, 0.006672, 0.006645, 0.006923, 0.006615, 0.006982, 0.006622, 0.006666, 0.006579, 0.006823, 0.006673, 0.006786, 0.006740, 0.006440, 0.006606, 0.006650, 0.006568, 0.006206, 0.006646, 0.006305, 0.006468, 0.006442, 0.006486, 0.006080, 0.006291, 0.006622, 0.006113, 0.006506, 0.006254, 0.006114, 0.006161, 0.006002, 0.006211, 0.006146, 0.006012, 0.005902, 0.006264, 0.005996, 0.006271, 0.005793, 0.006043, 0.006067, 0.006177, 0.005842, 0.005991, 0.005872, 0.006102, 0.006003, 0.005930, 0.006201, 0.006224, 0.005937, 0.005901, 0.005992, 0.005788, 0.006121, 0.005787, 0.005766, 0.005582, 0.006006, 0.005774, 0.005672, 0.005841, 0.005660, 0.006000, 0.005741, 0.005737, 0.005475, 0.005773, 0.005799, 0.005462, 0.005610, 0.005569, 0.005574, 0.005871, 0.005589, 0.005234, 0.005535, 0.005314, 0.005288, 0.005676, 0.005371, 0.005319, 0.005734, 0.005360, 0.005679, 0.005593};
  for(Int_t i=0; i<400; i++){
    fHistoPtWeight->SetBinContent(i+1,binc[i]);
  }
  fUseWeight=kTRUE;
}

//_________________________________________________________________________
void AliAnalysisTaskSEDs::SetPtWeightsFromFONLL5andTAMUoverLHC16i2abc(){
  // weight function from the ratio of the LHC16i2a+b+c MC
  // and FONLL calculations for pp data
  // corrected by the TAMU Raa calculation for 0-10% CC (not available in 30-50% CC)
  if(fHistoPtWeight) delete fHistoPtWeight;
  fHistoPtWeight = new TH1F("histoWeight","histoWeight",400,0.,40.);
  Float_t binc[400]={1.179906, 1.091249, 1.047774, 1.045579, 1.071679, 1.112413, 1.167414, 1.236240, 1.310301, 1.390289, 1.471711, 1.553389, 1.626886, 1.692115, 1.760647, 1.813658, 1.850817, 1.886699, 1.907671, 1.934832, 1.955433, 1.966727, 1.987262, 1.996316, 2.013326, 1.973926, 1.931144, 1.871654, 1.812942, 1.752718, 1.690846, 1.635303, 1.572611, 1.523510, 1.459790, 1.402510, 1.331908, 1.261575, 1.192241, 1.127915, 1.061798, 0.998830, 0.933514, 0.871774, 0.812936, 0.762844, 0.719340, 0.686587, 0.644108, 0.615714, 0.579512, 0.545254, 0.510508, 0.479884, 0.447423, 0.426154, 0.408934, 0.388264, 0.376424, 0.361389, 0.347757, 0.331685, 0.318029, 0.305285, 0.290922, 0.278523, 0.269807, 0.262025, 0.254878, 0.249325, 0.238179, 0.230899, 0.224792, 0.216253, 0.207879, 0.204465, 0.201153, 0.195373, 0.190926, 0.185773, 0.178589, 0.175371, 0.168959, 0.167004, 0.161705, 0.156809, 0.152788, 0.149806, 0.146429, 0.143478, 0.140037, 0.134813, 0.134679, 0.128205, 0.126078, 0.125038, 0.122214, 0.116329, 0.115044, 0.112427, 0.110279, 0.108098, 0.105784, 0.102628, 0.101429, 0.099101, 0.095464, 0.093631, 0.091491, 0.090045, 0.088374, 0.086188, 0.085067, 0.083168, 0.080636, 0.079414, 0.077610, 0.075013, 0.074825, 0.073932, 0.071106, 0.071050, 0.069574, 0.066593, 0.066924, 0.064876, 0.065064, 0.062345, 0.061980, 0.060859, 0.061616, 0.058952, 0.058079, 0.057894, 0.058031, 0.056604, 0.055180, 0.054490, 0.053909, 0.051768, 0.052210, 0.049552, 0.050152, 0.048955, 0.047953, 0.047224, 0.045588, 0.044985, 0.043728, 0.042934, 0.043434, 0.041834, 0.040118, 0.040281, 0.039348, 0.038987, 0.037793, 0.036258, 0.036420, 0.035528, 0.034761, 0.033524, 0.033296, 0.033280, 0.031825, 0.031351, 0.030329, 0.031103, 0.030401, 0.029481, 0.029247, 0.029352, 0.029174, 0.028286, 0.028500, 0.028017, 0.027293, 0.027932, 0.026779, 0.026379, 0.026628, 0.026211, 0.025508, 0.025877, 0.025433, 0.025328, 0.024636, 0.025069, 0.024282, 0.023625, 0.023278, 0.023074, 0.023000, 0.022943, 0.022514, 0.021767, 0.022180, 0.021594, 0.022175, 0.021944, 0.021456, 0.020901, 0.021419, 0.021230, 0.020738, 0.020322, 0.020055, 0.019686, 0.019371, 0.019725, 0.018835, 0.019029, 0.018163, 0.018398, 0.018163, 0.017719, 0.018126, 0.017208, 0.017086, 0.016622, 0.016865, 0.015663, 0.015791, 0.015108, 0.015069, 0.015033, 0.015006, 0.014940, 0.014604, 0.014133, 0.013968, 0.013904, 0.013934, 0.013780, 0.013930, 0.013727, 0.013940, 0.013763, 0.013826, 0.014192, 0.014801, 0.014347, 0.014048, 0.014009, 0.014197, 0.014571, 0.014999, 0.015030, 0.014491, 0.014891, 0.014456, 0.014596, 0.015256, 0.014648, 0.014492, 0.014756, 0.015344, 0.014986, 0.015433, 0.015394, 0.015756, 0.014778, 0.015145, 0.015478, 0.015051, 0.014986, 0.015067, 0.015793, 0.015748, 0.015188, 0.015502, 0.015533, 0.015340, 0.015759, 0.015745, 0.016026, 0.015635, 0.015194, 0.014579, 0.015225, 0.014963, 0.015365, 0.016030, 0.015387, 0.016341, 0.015327, 0.015340, 0.015030, 0.015246, 0.015420, 0.015015, 0.015195, 0.016021, 0.015034, 0.015528, 0.015114, 0.015423, 0.015564, 0.015348, 0.015107, 0.015314, 0.015411, 0.015243, 0.015154, 0.016324, 0.015215, 0.014823, 0.015030, 0.015104, 0.014896, 0.015400, 0.015721, 0.015131, 0.014951, 0.014630, 0.014597, 0.015235, 0.014583, 0.015418, 0.014648, 0.014769, 0.014601, 0.015167, 0.014857, 0.015134, 0.015053, 0.014405, 0.014800, 0.014921, 0.014760, 0.013966, 0.014979, 0.014230, 0.014620, 0.014581, 0.014701, 0.013799, 0.014299, 0.015071, 0.013931, 0.014846, 0.014290, 0.013988, 0.014113, 0.013767, 0.014263, 0.014131, 0.013840, 0.013604, 0.014456, 0.013853, 0.014505, 0.013416, 0.014010, 0.014081, 0.014352, 0.013589, 0.013952, 0.013690, 0.014241, 0.014024, 0.013868, 0.014517, 0.014587, 0.013927, 0.013857, 0.014084, 0.013619, 0.014417, 0.013644, 0.013607, 0.013185, 0.014200, 0.013665, 0.013437, 0.013849, 0.013431, 0.014252, 0.013648, 0.013652, 0.013039, 0.013761, 0.013836, 0.013043, 0.013408, 0.013319, 0.013344, 0.014065, 0.013400, 0.012560, 0.013294, 0.012773, 0.012721, 0.013663, 0.012939, 0.012823, 0.013835, 0.012942, 0.013723, 0.013525};
  for(Int_t i=0; i<400; i++){
    fHistoPtWeight->SetBinContent(i+1,binc[i]);
  }
  fUseWeight=kTRUE;
}

//_________________________________________________________________________
Double_t AliAnalysisTaskSEDs::GetPtWeightFromHistogram(Double_t pt)
{
  //
  // Using an histogram as weight function
  //  weight = 0 in the range outside the function
  //
  Double_t weight = 0.0;
  Int_t histoNbins = fHistoPtWeight->GetNbinsX();
  Int_t bin2 = fHistoPtWeight->FindBin(pt);
  if( (bin2>0) && (bin2<=histoNbins) ) {
    Int_t bin1=bin2-1;
    Int_t bin3=bin2+1;
    if(bin2==1) bin1=bin2+2;
    if(bin2==histoNbins) bin3=bin2-2;
    Float_t x_1=fHistoPtWeight->GetXaxis()->GetBinCenter(bin1);
    Float_t x_2=fHistoPtWeight->GetXaxis()->GetBinCenter(bin2);
    Float_t x_3=fHistoPtWeight->GetXaxis()->GetBinCenter(bin3);
    Float_t y_1=fHistoPtWeight->GetBinContent(bin1);
    Float_t y_2=fHistoPtWeight->GetBinContent(bin2);
    Float_t y_3=fHistoPtWeight->GetBinContent(bin3);
    Double_t a=( (y_3-y_2)*(x_1-x_2) - (y_1-y_2)*(x_3-x_2) )/( (x_3*x_3-x_2*x_2)*(x_1-x_2) - (x_1*x_1-x_2*x_2)*(x_3-x_2) );
    Double_t b=((y_1-y_2)-a*(x_1*x_1-x_2*x_2))/(x_1-x_2);
    Double_t c=y_3-a*(x_3*x_3)-b*x_3;
    weight = a*pt*pt+b*pt+c;
  }
  return weight;
}







