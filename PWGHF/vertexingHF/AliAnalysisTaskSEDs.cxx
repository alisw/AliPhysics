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
    delete fHistNEvents;
    delete fHistAllV0multNTPCout;
    delete fHistSelV0multNTPCout;
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

	if(fFillSparse) {
                    
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
		if(candType==1.5) fnSparseMC[2]->Fill(var4nSparse);
		if(candType==2.5) fnSparseMC[3]->Fill(var4nSparse);
	      }
	      else if(fFillSparseDplus && labDplus>=0 && pdgCode0==321) {
		if(candType==1.5) fnSparseMCDplus[2]->Fill(var4nSparse);
		if(candType==2.5) fnSparseMCDplus[3]->Fill(var4nSparse);
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
		if(candType==1.5) fnSparseMC[2]->Fill(var4nSparse);
		if(candType==2.5) fnSparseMC[3]->Fill(var4nSparse);
	      }
	      else if(fFillSparseDplus && labDplus>=0 && pdgCode0==211) {
		if(candType==1.5) fnSparseMCDplus[2]->Fill(var4nSparse);
		if(candType==2.5) fnSparseMCDplus[3]->Fill(var4nSparse);
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
                    
	  if(isFidAcc) {
	    Double_t var4nSparseAcc[2] = {pt,rapid*10};
	    if(isDaugInAcc) {
	      if(orig==4) fnSparseMC[0]->Fill(var4nSparseAcc);
	      if(orig==5) fnSparseMC[1]->Fill(var4nSparseAcc);
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
                    
	  if(isFidAcc) {
	    Double_t var4nSparseAcc[2] = {pt,rapid*10};
	    if(isDaugInAcc) {
	      if(orig==4) fnSparseMCDplus[0]->Fill(var4nSparseAcc);
	      if(orig==5) fnSparseMCDplus[1]->Fill(var4nSparseAcc);
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









