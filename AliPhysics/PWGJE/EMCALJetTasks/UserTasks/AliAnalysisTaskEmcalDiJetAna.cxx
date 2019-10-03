//
// Dijet analysis task.
//
// Author: M.Verweij

#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TProfile.h>
#include <TChain.h>
#include <TSystem.h>
#include <TFile.h>
#include <TKey.h>

#include "AliVCluster.h"
#include "AliVTrack.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliLog.h"
#include "AliEmcalParticle.h"
#include "AliMCEvent.h"
#include "AliGenPythiaEventHeader.h"
#include "AliAODMCHeader.h"
#include "AliMCEvent.h"
#include "AliAnalysisManager.h"
#include "AliJetContainer.h"
#include "AliCentrality.h"

#include "AliAnalysisTaskEmcalDiJetAna.h"

ClassImp(AliAnalysisTaskEmcalDiJetAna)

//________________________________________________________________________
AliAnalysisTaskEmcalDiJetAna::AliAnalysisTaskEmcalDiJetAna() : 
  AliAnalysisTaskEmcalDiJetBase("AliAnalysisTaskEmcalDiJetAna"),
  fDoMatchFullCharged(kTRUE),
  fNKtBins(30),
  fNDiJetEtaBins(1),
  fNAjBins(1),
  fh2CentRhoCh(0),
  fh2CentRhoScaled(0),
  fh3PtEtaPhiJetFull(0),
  fh3PtEtaPhiJetCharged(0),
  fhnDiJetVarsFull(0),
  fhnDiJetVarsCh(0),
  fhnDiJetVarsFullCharged(0),
  fhnMatchingFullCharged(0),
  fh3PtTrigKt1Kt2Ch(0),
  fh3PtTrigKt1Kt2FuCh(0),
  fh3PtTrigDPhi1DPhi2Ch(0),
  fh3PtTrigDPhi1DPhi2FuCh(0)
{
  // Default constructor.

  for(Int_t i=0; i<4; i++) {
    fh3DiJetKtNEFPtAssoc[i]          = 0;
    fCentCorrPtAssocCh[i]            = 0;
    fCentCorrPtAssocFuCh[i]          = 0;
    fAjPtAssocCentCh[i]              = 0;
    fAjPtAssocCentFuCh[i]            = 0;
    fh3PtAssoc1PtAssoc2DPhi23Ch[i]   = 0;
    fh3PtAssoc1PtAssoc2DPhi23FuCh[i] = 0;
  }  

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskEmcalDiJetAna::AliAnalysisTaskEmcalDiJetAna(const char *name) : 
  AliAnalysisTaskEmcalDiJetBase(name),
  fDoMatchFullCharged(kTRUE),
  fNKtBins(30),
  fNDiJetEtaBins(1),
  fNAjBins(1),
  fh2CentRhoCh(0),
  fh2CentRhoScaled(0),
  fh3PtEtaPhiJetFull(0),
  fh3PtEtaPhiJetCharged(0),
  fhnDiJetVarsFull(0),
  fhnDiJetVarsCh(0),
  fhnDiJetVarsFullCharged(0),
  fhnMatchingFullCharged(0),
  fh3PtTrigKt1Kt2Ch(0),
  fh3PtTrigKt1Kt2FuCh(0),
  fh3PtTrigDPhi1DPhi2Ch(0),
  fh3PtTrigDPhi1DPhi2FuCh(0)
{
  // Standard constructor.

  for(Int_t i=0; i<4; i++) {
    fh3DiJetKtNEFPtAssoc[i]      = 0;
    fCentCorrPtAssocCh[i]        = 0;
    fCentCorrPtAssocFuCh[i]      = 0;
    fAjPtAssocCentCh[i]          = 0;
    fAjPtAssocCentFuCh[i]        = 0;
    fh3PtAssoc1PtAssoc2DPhi23Ch[i]   = 0;
    fh3PtAssoc1PtAssoc2DPhi23FuCh[i] = 0;
  }

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskEmcalDiJetAna::~AliAnalysisTaskEmcalDiJetAna()
{
  // Destructor.
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalDiJetAna::RetrieveEventObjects() {
  //
  // retrieve event objects
  //

  if (!AliAnalysisTaskEmcalDiJetBase::RetrieveEventObjects())
    return kFALSE;

  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskEmcalDiJetAna::UserCreateOutputObjects()
{
  // Create user output.

  AliAnalysisTaskEmcalDiJetBase::UserCreateOutputObjects();

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  const Int_t nBinsCent = 100;
  Double_t minCent = 0.;
  Double_t maxCent = 100.;

  const Int_t nBinsRho = 200;
  Double_t minRho = 0.;
  Double_t maxRho = 20.;
  fh2CentRhoCh = new TH2F("fh2CentRhoCh","fh2CentRhoCh;centrality;#rho_{ch}",nBinsCent,minCent,maxCent,nBinsRho,minRho,maxRho);
  fOutput->Add(fh2CentRhoCh);
  fh2CentRhoScaled = new TH2F("fh2CentRhoScaled","fh2CentRhoScaled;centrality;s_{EMC}#rho_{ch}",nBinsCent,minCent,maxCent,nBinsRho,minRho,maxRho);
  fOutput->Add(fh2CentRhoScaled);

  const Int_t nBinsPt = 150;
  Double_t minPt = -20.;
  Double_t maxPt = 130.;
  const Int_t nBinsEta = 40;
  Double_t minEta = -1.;
  Double_t maxEta = 1.;
  const Int_t nBinsPhi = 18*6;
  Double_t minPhi = 0.;
  Double_t maxPhi = TMath::TwoPi();

  fh3PtEtaPhiJetFull = new TH3F("fh3PtEtaPhiJetFull","fh3PtEtaPhiJetFull;#it{p}_{T}^{jet};#eta;#varphi",nBinsPt,minPt,maxPt,nBinsEta,minEta,maxEta,nBinsPhi,minPhi,maxPhi);
  fOutput->Add(fh3PtEtaPhiJetFull);

  fh3PtEtaPhiJetCharged = new TH3F("fh3PtEtaPhiJetCharged","fh3PtEtaPhiJetCharged;#it{p}_{T}^{jet};#eta;#varphi",nBinsPt,minPt,maxPt,nBinsEta,minEta,maxEta,nBinsPhi,minPhi,maxPhi);
  fOutput->Add(fh3PtEtaPhiJetCharged);

  const Int_t nBinsSparse0 = 7;
  const Int_t nBinsPtW      = 30;
  const Int_t nBinsDPhi     = 72;
  const Int_t nBinsKt       = fNKtBins;
  const Int_t nBinsDiJetEta = fNDiJetEtaBins;//40;
  const Int_t nBinsCentr    = fNcentBins;
  const Int_t nBinsAj       = fNAjBins;//20
  const Int_t nBins0[nBinsSparse0] = {nBinsPtW,nBinsPtW,nBinsDPhi,nBinsKt,nBinsDiJetEta,nBinsCentr,nBinsAj};
  //pT1, pT2, deltaPhi, kT
  const Double_t xmin0[nBinsSparse0]  = {  minPt, minPt, -0.5*TMath::Pi(),   0.,-1.,0.  , 0.};
  const Double_t xmax0[nBinsSparse0]  = {  maxPt, maxPt,  1.5*TMath::Pi(), 120., 1.,100., 1.};
  const Double_t centArrayBins[8] = {0.,2.,5.,10.,20.,40.,60.,100.};

  if(fDoChargedCharged) {
    fhnDiJetVarsCh = new THnSparseF("fhnDiJetVarsCh",
				    "fhnDiJetVarsCh;#it{p}_{T,1} (GeV/#it{c});#it{p}_{T,2} (GeV/#it{c});#Delta#varphi;#it{k}_{T} = #it{p}_{T,1}sin(#Delta#varphi) (GeV/#it{c});(#eta_{1}+#eta_{2})/2);centrality;#it{A}_{j}",
				    nBinsSparse0,nBins0,xmin0,xmax0);
    if(fNcentBins==7) fhnDiJetVarsCh->SetBinEdges(5,centArrayBins);
    fOutput->Add(fhnDiJetVarsCh);
  }

  if(fDoFullCharged) {
    fhnDiJetVarsFullCharged = new THnSparseF("fhnDiJetVarsFullCharged",
				"fhnDiJetVarsFullCharged;#it{p}_{T,1} (GeV/#it{c});#it{p}_{T,2} (GeV/#it{c});#Delta#varphi;#it{k}_{T} = #it{p}_{T,1}sin(#Delta#varphi) (GeV/#it{c});(#eta_{1}+#eta_{2})/2);centrality;#it{A}_{j}",
				nBinsSparse0,nBins0,xmin0,xmax0);
    if(fNcentBins==7) fhnDiJetVarsFullCharged->SetBinEdges(5,centArrayBins);
    fOutput->Add(fhnDiJetVarsFullCharged);
  }

  if(fDoFullFull) {
    fhnDiJetVarsFull = new THnSparseF("fhnDiJetVarsFull",
				    "fhnDiJetVarsFull;#it{p}_{T,1} (GeV/#it{c});#it{p}_{T,2} (GeV/#it{c});#Delta#varphi;#it{k}_{T} = #it{p}_{T,1}sin(#Delta#varphi) (GeV/#it{c});(#eta_{1}+#eta_{2})/2);centrality;#it{A}_{j}",
				    nBinsSparse0,nBins0,xmin0,xmax0);
    fOutput->Add(fhnDiJetVarsFull);
  }

  fh3PtTrigKt1Kt2Ch = new TH3F("fh3PtTrigKt1Kt2Ch","fh3PtTrigKt1Kt2Ch;#it{p}_{T,1} (GeV/#it{c});#it{k}_{T,1};#it{k}_{T,2}",nBinsPt,minPt,maxPt,nBinsKt,0.,120.,nBinsKt,0.,120.);
  fOutput->Add(fh3PtTrigKt1Kt2Ch);  

  fh3PtTrigKt1Kt2FuCh = new TH3F("fh3PtTrigKt1Kt2FuCh","fh3PtTrigKt1Kt2FuCh;#it{p}_{T,1} (GeV/#it{c});#it{k}_{T,1};#it{k}_{T,2}",nBinsPt,minPt,maxPt,nBinsKt,0.,120.,nBinsKt,0.,120.);
  fOutput->Add(fh3PtTrigKt1Kt2FuCh); 

  fh3PtTrigDPhi1DPhi2Ch = new TH3F("fh3PtTrigDPhi1DPhi2Ch","fh3PtTrigDPhi1DPhi2Ch;#it{p}_{T,1} (GeV/#it{c});#it{k}_{T,1};#it{k}_{T,2}",nBinsPt,minPt,maxPt,36,-0.5*TMath::Pi(),1.5*TMath::Pi(),36,-0.5*TMath::Pi(),1.5*TMath::Pi());
  fOutput->Add(fh3PtTrigDPhi1DPhi2Ch);  

  fh3PtTrigDPhi1DPhi2FuCh = new TH3F("fh3PtTrigDPhi1DPhi2FuCh","fh3PtTrigDPhi1DPhi2FuCh;#it{p}_{T,1} (GeV/#it{c});#it{k}_{T,1};#it{k}_{T,2}",nBinsPt,minPt,maxPt,36,-0.5*TMath::Pi(),1.5*TMath::Pi(),36,-0.5*TMath::Pi(),1.5*TMath::Pi());
  fOutput->Add(fh3PtTrigDPhi1DPhi2FuCh);  
 

  for(Int_t i=0; i<4; i++) {
    TString histoName = Form("fh3DiJetKtNEFPtAssoc_TrigBin%d",i);
    fh3DiJetKtNEFPtAssoc[i] = new TH3F(histoName.Data(),histoName.Data(),nBinsKt,0.,120.,50,0.,1.,nBinsPt,minPt,maxPt);
    fOutput->Add(fh3DiJetKtNEFPtAssoc[i]);

    histoName = Form("fCentCorrPtAssocCh_TrigBin%d",i);
    fCentCorrPtAssocCh[i]  = new TH3F(histoName.Data(),histoName.Data(),10,0.,100.,10,0.,100.,nBinsPt,minPt,maxPt);
    fOutput->Add(fCentCorrPtAssocCh[i]);

    histoName = Form("fCentCorrPtAssocFuCh_TrigBin%d",i);
    fCentCorrPtAssocFuCh[i]  = new TH3F(histoName.Data(),histoName.Data(),10,0.,100.,10,0.,100.,nBinsPt,minPt,maxPt);
    fOutput->Add(fCentCorrPtAssocFuCh[i]);

    histoName = Form("fAjPtAssocCentCh_TrigBin%d",i);
    fAjPtAssocCentCh[i]  = new TH3F(histoName.Data(),histoName.Data(),50,0.,1.,nBinsPt,minPt,maxPt,20,0.,100.);
    fOutput->Add(fAjPtAssocCentCh[i]);

    histoName = Form("fAjPtAssocCentFuCh_TrigBin%d",i);
    fAjPtAssocCentFuCh[i]  = new TH3F(histoName.Data(),histoName.Data(),50,0.,1.,nBinsPt,minPt,maxPt,20,0.,100.);
    fOutput->Add(fAjPtAssocCentFuCh[i]);

    histoName = Form("fh3PtAssoc1PtAssoc2DPhi23Ch_TrigBin%d",i);
    fh3PtAssoc1PtAssoc2DPhi23Ch[i] = new TH3F(histoName.Data(),histoName.Data(),nBinsPt,minPt,maxPt,nBinsPt,minPt,maxPt,nBinsDPhi,-0.5*TMath::Pi(),1.5*TMath::Pi());
    fOutput->Add(fh3PtAssoc1PtAssoc2DPhi23Ch[i]);

    histoName = Form("fh3PtAssoc1PtAssoc2DPhi23FuCh_TrigBin%d",i);
    fh3PtAssoc1PtAssoc2DPhi23FuCh[i] = new TH3F(histoName.Data(),histoName.Data(),nBinsPt,minPt,maxPt,nBinsPt,minPt,maxPt,nBinsDPhi,-0.5*TMath::Pi(),1.5*TMath::Pi());
    fOutput->Add(fh3PtAssoc1PtAssoc2DPhi23FuCh[i]);
  }

  const Int_t nBinsSparseMatch = 7;
  const Int_t nBinsDPhiMatch = 80;
  const Int_t nBinsDEtaMatch = 80;
  const Int_t nBinsDR = 20;
  const Int_t nBinsFraction = 21;
  const Int_t nBinsType = 3;
  const Int_t nBinsMatch[nBinsSparseMatch] = {nBinsPt,nBinsPt,nBinsDPhiMatch,nBinsDEtaMatch,nBinsDR,nBinsFraction,nBinsType};
  //pTfull, pTch, deltaPhi, deltaEta, deltaR, fraction, jet type (leading,subleading,other)
  const Double_t xminMatch[nBinsSparseMatch]  = { minPt, minPt, -0.5,-0.5, 0. ,0.  ,0};
  const Double_t xmaxMatch[nBinsSparseMatch]  = { maxPt, maxPt,  0.5, 0.5, 0.5,1.05,3};
  if(fDoMatchFullCharged) {
    fhnMatchingFullCharged = new THnSparseF("fhnMatchingFullCharged","fhnMatchingFullCharged;#it{p}_{T,full} (GeV/#it{c});#it{p}_{T,ch} (GeV/#it{c});#Delta#varphi;#Delta#eta;#Delta R;f_{ch};type",
					  nBinsSparseMatch,nBinsMatch,xminMatch,xmaxMatch);
    fOutput->Add(fhnMatchingFullCharged);
  }
  
  // =========== Switch on Sumw2 for all histos ===========
  for (Int_t i=0; i<fOutput->GetEntries(); ++i) {
    TH1 *h1 = dynamic_cast<TH1*>(fOutput->At(i));
    if (h1){
      h1->Sumw2();
      continue;
    }
    THnSparse *hn = dynamic_cast<THnSparse*>(fOutput->At(i));
    if(hn) {
      hn->Sumw2();
      continue;
    }
  }

  TH1::AddDirectory(oldStatus);

  PostData(1, fOutput); // Post data for ALL output slots > 0 here.
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalDiJetAna::FillHistograms()
{
  // Fill histograms.

  fh2CentRhoCh->Fill(fCent,fRhoChVal);
  fh2CentRhoScaled->Fill(fCent,fRhoFullVal);

  Int_t nJetsFull = GetNJets(fContainerFull);
  for(Int_t ij=0; ij<nJetsFull; ij++) {
    AliEmcalJet *jet = static_cast<AliEmcalJet*>(GetAcceptJetFromArray(ij, fContainerFull));
    if(!jet) continue; //jet not selected

    Double_t jetPt = GetJetPt(jet,0);
    fh3PtEtaPhiJetFull->Fill(jetPt,jet->Eta(),jet->Phi());
  }

  Int_t nJetsCh = GetNJets(fContainerCharged);
  for(Int_t ij=0; ij<nJetsCh; ij++) {
    AliEmcalJet *jet = static_cast<AliEmcalJet*>(GetAcceptJetFromArray(ij, fContainerCharged));
    if(!jet) continue; //jet not selected

      Double_t jetPt = GetJetPt(jet,1);
      fh3PtEtaPhiJetCharged->Fill(jetPt,jet->Eta(),jet->Phi());
  }

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalDiJetAna::Run()
{
  // Run analysis code here, if needed. It will be executed before FillHistograms().

  //Check if event is selected (vertex & pile-up)
  if(!SelectEvent())
    return kFALSE;

  if(fRhoType==0) {
    fRhoFullVal = 0.;
    fRhoChVal = 0.;
  }
  if(fRhoType==1) {
    fRhoFullVal = GetRhoVal(fContainerFull);
    fRhoChVal = GetRhoVal(fContainerCharged);
  }
  
  if(fDoFullFull)
    CorrelateJets(0);

  // MatchFullAndChargedJets();
  if(fDoChargedCharged)   CorrelateJets(1);

  if(fDoMatchFullCharged) FillMatchFullChargedHistos(fContainerFull,fContainerCharged);

  if(fDoFullCharged)      {
    SetChargedFractionIndex();
    CorrelateJets(2);
  }

  return kTRUE;  // If return kFALSE FillHistogram() will NOT be executed.
}

//________________________________________________________________________
void AliAnalysisTaskEmcalDiJetAna::CorrelateJets(const Int_t type) {
  //
  // Correlate jets and fill histos
  //

  if( fJetCorrelationType==kCorrelateAll )
    CorrelateAllJets(type);
  else if( fJetCorrelationType==kCorrelateTwo )
    CorrelateTwoJets(type);
  else if( fJetCorrelationType==kCorrelateLS )
    CorrelateLeadingSubleadingJets(type);

  return;

}

//________________________________________________________________________
void AliAnalysisTaskEmcalDiJetAna::CorrelateTwoJets(const Int_t type) {
  //
  // Correlate jets and fill histos
  //

  Int_t typet = 0;
  Int_t typea = 0;
  if(type==0) { //full-full
    typet = fContainerFull;
    typea = fContainerFull;
  }
  else if(type==1) { //charged-charged
    typet = fContainerCharged;
    typea = fContainerCharged;
  }
  else if(type==2) { //full-charged
    typet = fContainerFull;
    typea = fContainerCharged;
  }
  else {
    AliWarning(Form("%s: type %d of dijet correlation not defined!",GetName(),type));
    return;
  }

  Int_t nJetsTrig  = GetNJets(typet);
  for(Int_t ijt=0; ijt<nJetsTrig; ijt++) {

    AliEmcalJet *jetTrig = NULL; 
    if(type==0) {
      jetTrig = static_cast<AliEmcalJet*>(GetJetFromArray(ijt, typet));
      if(TMath::Abs(jetTrig->Eta())>0.5)
	jetTrig = NULL;
    }
    else
      jetTrig = static_cast<AliEmcalJet*>(GetAcceptJetFromArray(ijt, typet));

    if(!jetTrig)
      continue; //jet not selected
    
    Double_t jetTrigPt = GetJetPt(jetTrig,typet);

    if(jetTrigPt<fPtMinTriggerJet)
      continue;

    AliEmcalJet *jetAssoc = GetLeadingAssociatedJet(typea,jetTrig);
    if(!jetAssoc)
      continue;

    if(fDoPtBias) {
      if(type==0 || type==1) {
	if(GetJetPt(jetAssoc,typea)>jetTrigPt)
	  continue;
      }
    }

    FillDiJetHistos(jetTrig,jetAssoc, type);

    //Look for second jet on away side - 3-jet events
    AliEmcalJet *jetAssoc2 = GetSecondLeadingAssociatedJet(typea,jetTrig);
    if(jetAssoc2)
      FillThreeJetHistos(jetTrig,jetAssoc,jetAssoc2,type);
    
  }
}

//________________________________________________________________________
AliEmcalJet* AliAnalysisTaskEmcalDiJetAna::GetLeadingJet(const Int_t type) {

  //Get associated jet which is the leading jet in the opposite hemisphere

  Int_t cont = 0;
  if(type==0)  //full-full
    cont = fContainerFull;
  else if(type==1)  //charged-charged
    cont = fContainerCharged;
  else if(type==2)  //full-charged
    cont = fContainerFull;

  Int_t nJets = GetNJets(cont);
  Double_t ptLead = -999;
  Int_t    iJetLead = -1;
  for(Int_t ij=0; ij<nJets; ij++) {
    AliEmcalJet *jet = NULL;
    if(type==0) {
      jet = static_cast<AliEmcalJet*>(GetJetFromArray(ij, cont));
      if(TMath::Abs(jet->Eta())>0.5)
	jet = NULL;
    }
    else
      jet = static_cast<AliEmcalJet*>(GetAcceptJetFromArray(ij, cont));

    if(!jet)
      continue;

    Double_t jetPt = GetJetPt(jet,cont);

    if(jetPt>ptLead) {
      ptLead = jetPt;
      iJetLead = ij;
    }
  }

  AliEmcalJet *jetLead = static_cast<AliEmcalJet*>(GetJetFromArray(iJetLead, cont));

  return jetLead;
}

//________________________________________________________________________
AliEmcalJet* AliAnalysisTaskEmcalDiJetAna::GetLeadingAssociatedJet(const Int_t type, AliEmcalJet *jetTrig) {

  //Get associated jet which is the leading jet in the opposite hemisphere

  Int_t typea = 0;
  if(type==0)  //full-full
    typea = fContainerFull;
  else if(type==1)  //charged-charged
    typea = fContainerCharged;
  else if(type==2)  //full-charged
    typea = fContainerCharged;

  AliEmcalJet *jetAssocLead = GetLeadingJetOppositeHemisphere(type, typea, jetTrig);
  
  return jetAssocLead;
}

//________________________________________________________________________
AliEmcalJet* AliAnalysisTaskEmcalDiJetAna::GetSecondLeadingAssociatedJet(const Int_t type, AliEmcalJet *jetTrig) {

  //Get associated jet which is the leading jet in the opposite hemisphere

  Int_t typea = 0;
  if(type==0)  //full-full
    typea = fContainerFull;
  else if(type==1)  //charged-charged
    typea = fContainerCharged;
  else if(type==2)  //full-charged
    typea = fContainerCharged;

  AliEmcalJet *jetAssocLead2 = GetSecondLeadingJetOppositeHemisphere(type, typea, jetTrig);
  
  return jetAssocLead2;
}

//________________________________________________________________________
void AliAnalysisTaskEmcalDiJetAna::CorrelateAllJets(const Int_t type) {
  //
  // Correlate jets and fill histos
  //

  Int_t typet = 0;
  Int_t typea = 0;
  if(type==0) { //full-full
    typet = fContainerFull;
    typea = fContainerFull;
  }
  else if(type==1) { //charged-charged
    typet = fContainerCharged;
    typea = fContainerCharged;
  }
  else if(type==2) { //full-charged
    typet = fContainerFull;
    typea = fContainerCharged;
  }
  else {
    AliWarning(Form("%s: type %d of dijet correlation not defined!",GetName(),type));
    return;
  }

  Int_t nJetsTrig  = GetNJets(typet);
  Int_t nJetsAssoc = GetNJets(typea);

  for(Int_t ijt=0; ijt<nJetsTrig; ijt++) {
    AliEmcalJet *jetTrig = NULL; 
    if(type==0) {
      jetTrig = static_cast<AliEmcalJet*>(GetJetFromArray(ijt, typet));
      if(TMath::Abs(jetTrig->Eta())>0.5)
	jetTrig = NULL;
    }
    else
      jetTrig = static_cast<AliEmcalJet*>(GetAcceptJetFromArray(ijt, typet));

    if(!jetTrig)
      continue; //jet not selected
    
    Double_t jetTrigPt = GetJetPt(jetTrig,typet);

    if(jetTrigPt<fPtMinTriggerJet)
      continue;

    for(Int_t ija=0; ija<nJetsAssoc; ija++) {
      if(IsSameJet(ijt,ija,type)) continue;

      AliEmcalJet *jetAssoc = NULL;
      if(type==0) {
	jetAssoc = static_cast<AliEmcalJet*>(GetJetFromArray(ija, typea));
	if(TMath::Abs(jetAssoc->Eta())>0.5)
	  jetAssoc = NULL;
      }
      else
	jetAssoc = static_cast<AliEmcalJet*>(GetAcceptJetFromArray(ija, typea));

      if(!jetAssoc)
	continue;
	
      Double_t jetAssocPt = GetJetPt(jetAssoc,typea);
      
      if(jetTrigPt>jetAssocPt)
	FillDiJetHistos(jetTrig,jetAssoc, type);
      
    } // associate jet loop
  }//trigger jet loop

}

//________________________________________________________________________
void AliAnalysisTaskEmcalDiJetAna::CorrelateLeadingSubleadingJets(const Int_t type) {
  //
  // Correlate leading jet in event with leading jet in opposite hemisphere
  //

  Int_t typet = 0;
  Int_t typea = 0;
  if(type==0) { //full-full
    typet = fContainerFull;
    typea = fContainerFull;
  }
  else if(type==1) { //charged-charged
    typet = fContainerCharged;
    typea = fContainerCharged;
  }
  else if(type==2) { //full-charged
    typet = fContainerFull;
    typea = fContainerCharged;
  }
  else {
    AliWarning(Form("%s: type %d of dijet correlation not defined!",GetName(),type));
    return;
  }

  AliEmcalJet *jetTrig = GetLeadingJet(typet);
  if(!jetTrig)
    return;

  Double_t jetTrigPt = GetJetPt(jetTrig,typet);

  if(jetTrigPt<fPtMinTriggerJet)
    return;
  
  AliEmcalJet *jetAssoc = GetLeadingAssociatedJet(typea,jetTrig);
  if(!jetAssoc)
    return;
  
  FillDiJetHistos(jetTrig,jetAssoc, type);
}

//________________________________________________________________________
void AliAnalysisTaskEmcalDiJetAna::FillDiJetHistos(const AliEmcalJet *jet1, const AliEmcalJet *jet2, const Int_t mode) {
  //
  // Fill histos
  // mode: full vs full        = 0
  //       charged vs charged  = 1
  //       full vs charged     = 2
  //

  Int_t typet = 0;
  Int_t typea = 0;
  if(mode==0) { //full-full
    typet = fContainerFull;
    typea = fContainerFull;
  }
  else if(mode==1) { //charged-charged
    typet = fContainerCharged;
    typea = fContainerCharged;
  }
  else if(mode==2) { //full-charged
    typet = fContainerFull;
    typea = fContainerCharged;
  }
  else {
    AliWarning(Form("%s: mode %d of dijet correlation not defined!",GetName(),mode));
    return;
  }

  Double_t jetTrigPt = GetJetPt(jet1,typet);
  Double_t jetAssocPt = GetJetPt(jet2,typea);

  Double_t deltaPhi = GetDeltaPhi(jet1->Phi(),jet2->Phi());

  Double_t kT = TMath::Abs(jetTrigPt*TMath::Sin(deltaPhi));

  Double_t dijetEta = (jet1->Eta()+jet2->Eta())/2.;

  Double_t aj = 0.;
  if((jetTrigPt+jetAssocPt)>0.) aj = (jetTrigPt-jetAssocPt)/(jetTrigPt+jetAssocPt);

  Double_t diJetVars[7] = {jetTrigPt,jetAssocPt,deltaPhi,kT,dijetEta,fCent,aj};

  if(mode==0)
    fhnDiJetVarsFull->Fill(diJetVars);
  else if(mode==1)
    fhnDiJetVarsCh->Fill(diJetVars);
  else if(mode==2)
    fhnDiJetVarsFullCharged->Fill(diJetVars);

  Double_t dPhiMin = TMath::Pi() - 1./3.*TMath::Pi();
  Double_t dPhiMax = TMath::Pi() + 1./3.*TMath::Pi();
  Int_t trigBin = GetPtTriggerBin(jetTrigPt);
  if(mode==2) {
    if(trigBin>-1 && trigBin<4) {
      if(deltaPhi>dPhiMin && deltaPhi<dPhiMax)
	fh3DiJetKtNEFPtAssoc[trigBin]->Fill(kT, jet1->NEF(), jetAssocPt);
    }
  }

  //Fill centrality correlation histos in case a dijet is present in acceptance
  Double_t centZNA = -1.;
  AliCentrality *aliCent = InputEvent()->GetCentrality();
  if (aliCent) {
    centZNA = aliCent->GetCentralityPercentile("ZNA");
    if(trigBin>-1 && trigBin<4) {
      if(deltaPhi>dPhiMin && deltaPhi<dPhiMax) {
	if(mode==1) {
	  fCentCorrPtAssocCh[trigBin]->Fill(fCent,centZNA,jetAssocPt);
	  fAjPtAssocCentCh[trigBin]->Fill(aj,jetAssocPt,fCent);
	}
	else if(mode==2) {
	  fCentCorrPtAssocFuCh[trigBin]->Fill(fCent,centZNA,jetAssocPt);
	  fAjPtAssocCentFuCh[trigBin]->Fill(aj,jetAssocPt,fCent);
	}
      }
    }
  }

}

//________________________________________________________________________
void AliAnalysisTaskEmcalDiJetAna::FillThreeJetHistos(const AliEmcalJet *jet1, const AliEmcalJet *jet2, const AliEmcalJet *jet3, const Int_t mode) {
  //
  // Fill histos
  // mode: full vs full        = 0
  //       charged vs charged  = 1
  //       full vs charged     = 2
  //

  Int_t typet = 0;
  Int_t typea = 0;
  if(mode==0) { //full-full
    typet = fContainerFull;
    typea = fContainerFull;
  }
  else if(mode==1) { //charged-charged
    typet = fContainerCharged;
    typea = fContainerCharged;
  }
  else if(mode==2) { //full-charged
    typet = fContainerFull;
    typea = fContainerCharged;
  }
  else {
    AliWarning(Form("%s: mode %d of dijet correlation not defined!",GetName(),mode));
    return;
  }

  Double_t jetTrigPt = GetJetPt(jet1,typet);
  Double_t jetAssoc2Pt = GetJetPt(jet2,typea);
  Double_t jetAssoc3Pt = GetJetPt(jet3,typea);

  Double_t deltaPhi12 = GetDeltaPhi(jet1->Phi(),jet2->Phi());
  Double_t deltaPhi13 = GetDeltaPhi(jet1->Phi(),jet3->Phi());
  Double_t deltaPhi23 = GetDeltaPhi(jet2->Phi(),jet3->Phi());

  Double_t kT12 = TMath::Abs(jetTrigPt*TMath::Sin(deltaPhi12));
  Double_t kT13 = TMath::Abs(jetTrigPt*TMath::Sin(deltaPhi13));
  
  Double_t dPhiMin = TMath::Pi() - 1./3.*TMath::Pi();
  Double_t dPhiMax = TMath::Pi() + 1./3.*TMath::Pi();

  Int_t trigBin = GetPtTriggerBin(jetTrigPt);

  if(jetAssoc2Pt>20. && jetAssoc3Pt>20.) {
    if(mode==1) {
      fh3PtTrigDPhi1DPhi2Ch->Fill(jetTrigPt,deltaPhi12,deltaPhi13);
      fh3PtAssoc1PtAssoc2DPhi23Ch[trigBin]->Fill(jetAssoc2Pt,jetAssoc3Pt,deltaPhi23);    
    }
    else if(mode==1) {
      fh3PtTrigDPhi1DPhi2FuCh->Fill(jetTrigPt,deltaPhi12,deltaPhi13);
      fh3PtAssoc1PtAssoc2DPhi23FuCh[trigBin]->Fill(jetAssoc2Pt,jetAssoc3Pt,deltaPhi23);    
    }

    if(deltaPhi12>dPhiMin && deltaPhi12<dPhiMax) {
      if(mode==1)
	fh3PtTrigKt1Kt2Ch->Fill(jetTrigPt,kT12,kT13);
      else if(mode==2)
	fh3PtTrigKt1Kt2FuCh->Fill(jetTrigPt,kT12,kT13);
    }
  }

}

//________________________________________________________________________
Int_t AliAnalysisTaskEmcalDiJetAna::GetPtTriggerBin(Double_t pt) {

  Int_t binTrig = -1;
  if(pt>=20 && pt<40)
    binTrig = 0;
  else if(pt>=40 && pt<60)
    binTrig = 1;
  else if(pt>=60 && pt<80)
    binTrig = 2;
  else if(pt>=80 && pt<100)
    binTrig = 3;

  return binTrig;
}

//________________________________________________________________________
void AliAnalysisTaskEmcalDiJetAna::FillMatchFullChargedHistos(Int_t cFull,Int_t cCharged) {
  //
  // Match full to charged jets and fill histo
  //

  Int_t match = MatchFullAndChargedJets(cFull,cCharged);
  if(match==0) {
    AliDebug(11,Form("%s: matching failed",GetName()));
    return;
  }
  
  for(int ig = 0;ig < GetNJets(cFull);++ig){
    AliEmcalJet *jetFull = static_cast<AliEmcalJet*>(GetAcceptJetFromArray(ig, cFull));
    if(!jetFull) continue;

    AliEmcalJet *jetCh = jetFull->ClosestJet();
    if(!jetCh) continue;

    Double_t shFraction = GetFractionSharedPt(jetFull,jetCh);
    Double_t matchVars[7] = {
      jetFull->Pt(),
      jetCh->Pt(),
      GetDeltaPhi(jetFull->Phi(),jetCh->Phi()),
      jetFull->Eta()-jetCh->Eta(),
      GetDeltaR(jetFull,jetCh),
      shFraction,TMath::Min((Float_t)ig+0.5,2.5)
    };
      fhnMatchingFullCharged->Fill(matchVars);
    
  }//loop over full jets

}


//________________________________________________________________________
Int_t AliAnalysisTaskEmcalDiJetAna::MatchFullAndChargedJets(Int_t cFull, Int_t cCharged) {
  //
  // Match charged jets to full jets
  //

  if(GetNJets(cFull)<1) {
    AliDebug(2,Form("%s: no full jets: %d", GetName(),GetNJets(cFull)));
    return 0;
  }

  if(GetNJets(cCharged)<1) {
    AliDebug(2,Form("%s: no charged jets: %d", GetName(),GetNJets(cCharged)));
    return 0;
  }

  TClonesArray *cJetsFull = GetJetArray(cFull);
  TClonesArray *cJetsCharged = GetJetArray(cCharged);

  if(!cJetsFull) {
    AliDebug(2,Form("%s: no full jet array",GetName()));
    return 0;
  }

  if(!cJetsCharged) {
    AliDebug(2,Form("%s: no charged jet array",GetName()));
    return 0;
  }

  if(!fMatchingDone) {
      MatchJetsGeo(cFull, cCharged, 0);
      return 1;  
  } else {
    AliDebug(11,Form("%s: Matching already done before",GetName()));
    return 1;
  }

}

//_______________________________________________________________________
void AliAnalysisTaskEmcalDiJetAna::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
}



