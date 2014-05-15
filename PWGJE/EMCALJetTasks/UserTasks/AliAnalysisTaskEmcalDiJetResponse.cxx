//
// Dijet response analysis task.
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

#include "AliAnalysisTaskEmcalDiJetResponse.h"

ClassImp(AliAnalysisTaskEmcalDiJetResponse)

//________________________________________________________________________
AliAnalysisTaskEmcalDiJetResponse::AliAnalysisTaskEmcalDiJetResponse() : 
  AliAnalysisTaskEmcalDiJetBase("AliAnalysisTaskEmcalDiJetResponse"),
  fDoMatchFullCharged(kTRUE),
  fhnDiJetResponseCharged(0),
  fhnDiJetResponseFullCharged(0),
  fh1TriggersLostCharged(0),
  fh1TriggersLostFull(0),
  fh3AssocLostPtDeltaPhiCharged(0),
  fh3AssocLostPtDeltaPhiFull(0),
  fhnMatchingCharged(0),
  fhnMatchingFull(0),
  fnUsedResponseVar(0)
{
  // Default constructor.
  
  SetMakeGeneralHistograms(kTRUE);

  for(Int_t i = 0; i<2; i++) {
    fh1TriggersCharged[i] = 0;
    fh1TriggersFull[i] = 0;
  }

}

//________________________________________________________________________
AliAnalysisTaskEmcalDiJetResponse::AliAnalysisTaskEmcalDiJetResponse(const char *name) : 
  AliAnalysisTaskEmcalDiJetBase(name),
  fDoMatchFullCharged(kTRUE),
  fhnDiJetResponseCharged(0),
  fhnDiJetResponseFullCharged(0),
  fh1TriggersLostCharged(0),
  fh1TriggersLostFull(0),
  fh3AssocLostPtDeltaPhiCharged(0),
  fh3AssocLostPtDeltaPhiFull(0),
  fhnMatchingCharged(0),
  fhnMatchingFull(0),
  fnUsedResponseVar(0)
{
  // Standard constructor.

  SetMakeGeneralHistograms(kTRUE);

  for(Int_t i = 0; i<2; i++) {
    fh1TriggersCharged[i] = 0;
    fh1TriggersFull[i] = 0;
  }

}

//________________________________________________________________________
AliAnalysisTaskEmcalDiJetResponse::~AliAnalysisTaskEmcalDiJetResponse()
{
  // Destructor.
}


//________________________________________________________________________
void AliAnalysisTaskEmcalDiJetResponse::UserCreateOutputObjects()
{
  // Create user output.

  AliAnalysisTaskEmcalDiJetBase::UserCreateOutputObjects();

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  //Store dijet vars: pt,trig MC, pt,trig DET, pt,ass MC, pt,ass DET, dPhi MC, dPhi Det, kT MC, kT Det 
  const Int_t nBinsSparse0 = 10;
  const Int_t nBinsPt = 250;
  const Int_t nBinsDPhi     = 36;
  const Int_t nBinsKt       = 25;
  const Int_t nBinsDiJetEta = 40;
  const Int_t nBinsAj       = 50;
  const Int_t nBinsVar[2] = {nBinsKt,nBinsDiJetEta};

  const Int_t nBins0[nBinsSparse0] = {nBinsPt,nBinsPt,nBinsPt,nBinsPt,nBinsDPhi,nBinsDPhi,nBinsVar[fnUsedResponseVar],nBinsVar[fnUsedResponseVar],nBinsAj,nBinsAj};

  const Double_t minPt = 0.;
  const Double_t maxPt = 250.;
  const Double_t minVar[2] = {   0.,-1.};
  const Double_t maxVar[2] = { 100., 1.};

  const Double_t xmin0[nBinsSparse0]  = {  minPt, minPt, minPt, minPt, 0.5*TMath::Pi(), 0.5*TMath::Pi(), minVar[fnUsedResponseVar], minVar[fnUsedResponseVar],0.,0.};
  const Double_t xmax0[nBinsSparse0]  = {  maxPt, maxPt, maxPt, maxPt, 1.5*TMath::Pi(), 1.5*TMath::Pi(), maxVar[fnUsedResponseVar], maxVar[fnUsedResponseVar],1.,1.};

  fhnDiJetResponseCharged = new THnSparseF("fhnDiJetResponseCharged","fhnDiJetResponseCharged;p_{T,trig}^{part};p_{T,trig}^{det};p_{T,ass}^{part};p_{T,ass}^{det};#Delta#varphi_{part};#Delta#varphi_{det};k_{T}^{part},k_{T}^{det};A_{j}^{part}A_{j}^{det}",nBinsSparse0,nBins0,xmin0,xmax0);

  fhnDiJetResponseFullCharged = new THnSparseF("fhnDiJetResponseFullCharged","fhnDiJetResponseFullCharged;p_{T,trig}^{part};p_{T,trig}^{det};p_{T,ass}^{part};p_{T,ass}^{det};#Delta#varphi_{part};#Delta#varphi_{det};k_{T}^{part},k_{T}^{det};A_{j}^{part}A_{j}^{det}",nBinsSparse0,nBins0,xmin0,xmax0);

  if(fnUsedResponseVar==1) {
    fhnDiJetResponseCharged->SetTitle("fhnDiJetResponseCharged DiJetEta"); 
    fhnDiJetResponseCharged->GetAxis(6)->SetTitle("#eta_{dijet}^{part}");
    fhnDiJetResponseCharged->GetAxis(7)->SetTitle("#eta_{dijet}^{det}");

    fhnDiJetResponseFullCharged->SetTitle("fhnDiJetResponseFullCharged DiJetEta"); 
    fhnDiJetResponseFullCharged->GetAxis(6)->SetTitle("#eta_{dijet}^{part}");
    fhnDiJetResponseFullCharged->GetAxis(7)->SetTitle("#eta_{dijet}^{det}");
  }

  fOutput->Add(fhnDiJetResponseCharged);
  fOutput->Add(fhnDiJetResponseFullCharged);

  TString strType = "";
  for(Int_t i = 0; i<2; i++) {
    if(i==0)      strType="Part";
    else if(i==1) strType="Det";
    fh1TriggersCharged[i] = new TH1F(Form("fh1TriggersCharged%s",strType.Data()),Form("fh1TriggersCharged%s;p_{T,trig}^{ch}",strType.Data()),nBinsPt,minPt,maxPt);
    fOutput->Add(fh1TriggersCharged[i]);
  
    fh1TriggersFull[i] = new TH1F(Form("fh1TriggersFull%s",strType.Data()),Form("fh1TriggersFull%s;p_{T,trig}^{ch}",strType.Data()),nBinsPt,minPt,maxPt);
    fOutput->Add(fh1TriggersFull[i]);
  }

  fh1TriggersLostCharged = new TH1F("fh1TriggersLostCharged","fh1TriggersLostCharged;p_{T,trig}^{ch}",nBinsPt,minPt,maxPt);
  fOutput->Add(fh1TriggersLostCharged);

  fh1TriggersLostFull = new TH1F("fh1TriggersLostFull","fh1TriggersLostFull;p_{T,trig}^{ch}",nBinsPt,minPt,maxPt);
  fOutput->Add(fh1TriggersLostFull);

  fh3AssocLostPtDeltaPhiCharged = new TH3F("fh3AssocLostPtDeltaPhiCharged","fh3AssocLostPtDeltaPhiCharged;p_{T,trig}^{ch};p_{T,assoc}^{ch};#Delta#varphi",nBinsPt,minPt,maxPt,nBinsPt,minPt,maxPt,nBinsDPhi,-0.5*TMath::Pi(),1.5*TMath::Pi());
  fOutput->Add(fh3AssocLostPtDeltaPhiCharged);

  fh3AssocLostPtDeltaPhiFull = new TH3F("fh3AssocLostPtDeltaPhiFull","fh3AssocLostPtDeltaPhiFull;p_{T,trig}^{ch};p_{T,assoc}^{ch};#Delta#varphi",nBinsPt,minPt,maxPt,nBinsPt,minPt,maxPt,nBinsDPhi,-0.5*TMath::Pi(),1.5*TMath::Pi());
  fOutput->Add(fh3AssocLostPtDeltaPhiFull);

  const Int_t nBinsSparseMatch = 6;
  const Int_t nBinsDPhiMatch = 80;
  const Int_t nBinsDEtaMatch = 80;
  const Int_t nBinsDR = 20;
  const Int_t nBinsType = 3;
  const Int_t nBinsMatch[nBinsSparseMatch] = {nBinsPt,nBinsPt,nBinsDPhiMatch,nBinsDEtaMatch,nBinsDR,nBinsType};
  //pTpart, pTdet, deltaPhi, deltaEta, deltaR, jet type (leading,subleading,other)
  const Double_t xminMatch[nBinsSparseMatch]  = { minPt, minPt, -0.5,-0.5, 0., 0};
  const Double_t xmaxMatch[nBinsSparseMatch]  = { maxPt, maxPt,  0.5, 0.5, 0.5,3};
  fhnMatchingCharged = new THnSparseF("fhnMatchingCharged","fhnMatchingCharged;#it{p}_{T,part} (GeV/#it{c});#it{p}_{T,det} (GeV/#it{c});#Delta#varphi;#Delta#eta;#Delta R;type",
					  nBinsSparseMatch,nBinsMatch,xminMatch,xmaxMatch);
  fOutput->Add(fhnMatchingCharged);

  fhnMatchingFull = new THnSparseF("fhnMatchingFull","fhnMatchingFull;#it{p}_{T,part} (GeV/#it{c});#it{p}_{T,det} (GeV/#it{c});#Delta#varphi;#Delta#eta;#Delta R;type",
					  nBinsSparseMatch,nBinsMatch,xminMatch,xmaxMatch);
  fOutput->Add(fhnMatchingFull);


  // =========== Switch on Sumw2 for all histos ===========
  for (Int_t i=0; i<fOutput->GetEntries(); ++i) {
    TH1 *h1 = dynamic_cast<TH1*>(fOutput->At(i));
    if (h1){
      h1->Sumw2();
      continue;
    }
    THnSparse *hn = dynamic_cast<THnSparse*>(fOutput->At(i));
    if(hn)hn->Sumw2();
  }

  TH1::AddDirectory(oldStatus);

  PostData(1, fOutput); // Post data for ALL output slots > 0 here.
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalDiJetResponse::Run()
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
  

  //Do matching
  MatchJetsGeo(fContainerCharged,fContainerChargedMC,0,0.3,1);
  MatchJetsGeo(fContainerFull,fContainerFullMC,0,0.3,2);

  //Fill particle-detector level matching histos
  
  if(fDoChargedCharged)   CorrelateJets(1);

  if(fDoFullCharged) {
    SetChargedFractionIndexMC();
    CorrelateJets(2);
  }
  
  return kTRUE;  // If return kFALSE FillHistogram() will NOT be executed.
}

//________________________________________________________________________
void AliAnalysisTaskEmcalDiJetResponse::CorrelateJets(const Int_t type) {
  //
  // Correlate jets and fill histos
  //

  if( fJetCorrelationType==kCorrelateAll )
    CorrelateAllJets(type);
  else if( fJetCorrelationType==kCorrelateTwo )
    CorrelateTwoJets(type);
  else if( fJetCorrelationType==kCorrelateLS )
    AliWarning(Form("%s: leading-subleading correlation not implemented for response!",GetName()));

  return;

}


//________________________________________________________________________
void AliAnalysisTaskEmcalDiJetResponse::CorrelateAllJets(const Int_t type) {
  //
  // Correlate jets and fill histos
  //

  Int_t typet = 0;
  Int_t typetMC = 0;
  Int_t typeaMC = 0;
  if(type==0) { //full-full
    typetMC = fContainerFullMC;
    typeaMC = fContainerFullMC;
    typet = fContainerFull;
  }
  else if(type==1) { //charged-charged
    typetMC = fContainerChargedMC;
    typeaMC = fContainerChargedMC;
    typet = fContainerCharged;
  }
  else if(type==2) { //full-charged
    typetMC = fContainerFullMC;
    typeaMC = fContainerChargedMC;
    typet = fContainerFull;
  }
  else {
    AliWarning(Form("%s: type %d of dijet correlation not defined!",GetName(),type));
    return;
  }

  Int_t nJetsTrig  = 0;
  Int_t nJetsAssoc = 0;
  if(type==0) {
    nJetsTrig  = GetNJets(fContainerFullMC);
    nJetsAssoc = nJetsTrig;
  }
  else if(type==1) {
    nJetsTrig  = GetNJets(fContainerChargedMC);
    nJetsAssoc = nJetsTrig;
  }
  else if(type==2) {
    nJetsTrig  = GetNJets(fContainerFullMC);
    nJetsAssoc = GetNJets(fContainerChargedMC);
  }


  for(Int_t ijt=0; ijt<nJetsTrig; ijt++) {

    AliEmcalJet *jetTrigMC = static_cast<AliEmcalJet*>(GetAcceptJetFromArray(ijt, typetMC));
    if(!jetTrigMC) continue; //jet not selected

    Double_t jetTrigPtMC = GetJetPt(jetTrigMC,typetMC);

    if(jetTrigPtMC<fPtMinTriggerJet)
      continue;

    if(type==1)
      fh1TriggersCharged[0]->Fill(jetTrigPtMC);
    if(type==2)
      fh1TriggersFull[0]->Fill(jetTrigPtMC);

    AliEmcalJet *jetTrigDet = jetTrigMC->ClosestJet();
    if(!jetTrigDet) {
      //trigger is lost
      if(type==1)
	fh1TriggersLostCharged->Fill(jetTrigPtMC);
      if(type==2)
	fh1TriggersLostFull->Fill(jetTrigPtMC);
      
      continue;
    }

    if(type==1)
      fh1TriggersCharged[1]->Fill(GetJetPt(jetTrigDet,typet));
    if(type==2)
      fh1TriggersFull[1]->Fill(GetJetPt(jetTrigDet,typet));

    for(Int_t ija=0; ija<nJetsAssoc; ija++) {
      if(IsSameJet(ijt,ija,type,kTRUE)) continue;

      AliEmcalJet *jetAssocMC = static_cast<AliEmcalJet*>(GetAcceptJetFromArray(ija, typeaMC));
      if(!jetAssocMC) continue;

      Double_t jetAssocPtMC = GetJetPt(jetAssocMC,typeaMC);

      //Now check if jets are also there on detector level
      AliEmcalJet *jetAssocDet = jetAssocMC->ClosestJet();
      if(!jetAssocDet) {
	//dijet is lost
      if(type==1)
	fh3AssocLostPtDeltaPhiCharged->Fill(jetTrigPtMC,jetAssocPtMC,GetDeltaPhi(jetTrigMC,jetAssocMC));
      if(type==2)
	fh3AssocLostPtDeltaPhiFull->Fill(jetTrigPtMC,jetAssocPtMC,GetDeltaPhi(jetTrigMC,jetAssocMC));
	continue;
      }

      FillDiJetResponse(jetTrigMC,jetAssocMC,jetTrigDet,jetAssocDet,type);

    } // associate jet loop
  }//trigger jet loop

}

//________________________________________________________________________
void AliAnalysisTaskEmcalDiJetResponse::CorrelateTwoJets(const Int_t type) {
  //
  // Correlate jets and fill histos
  //

  Int_t typet = 0;
  Int_t typea = 0;
  Int_t typetMC = 0;
  Int_t typeaMC = 0;
  if(type==0) { //full-full
    typetMC = fContainerFullMC;
    typeaMC = fContainerFullMC;
    typet = fContainerFull;
    typea = fContainerFull;
  }
  else if(type==1) { //charged-charged
    typetMC = fContainerChargedMC;
    typeaMC = fContainerChargedMC;
    typet = fContainerCharged;
    typea = fContainerCharged;
  }
  else if(type==2) { //full-charged
    typetMC = fContainerFullMC;
    typeaMC = fContainerChargedMC;
    typet = fContainerFull;
    typea = fContainerCharged;
  }
  else {
    AliWarning(Form("%s: type %d of dijet correlation not defined!",GetName(),type));
    return;
  }

  Int_t nJetsTrig  = 0;
  if(type==0) {
    nJetsTrig  = GetNJets(fContainerFullMC);
  }
  else if(type==1) {
    nJetsTrig  = GetNJets(fContainerChargedMC);
  }
  else if(type==2) {
    nJetsTrig  = GetNJets(fContainerFullMC);
  }

  for(Int_t ijt=0; ijt<nJetsTrig; ijt++) {

    AliEmcalJet *jetTrigMC = static_cast<AliEmcalJet*>(GetAcceptJetFromArray(ijt, typetMC));
    if(!jetTrigMC) continue; //jet not selected

    Double_t jetTrigPtMC = GetJetPt(jetTrigMC,typetMC);

    if(jetTrigPtMC<fPtMinTriggerJet)
      continue;

    if(type==1)
      fh1TriggersCharged[0]->Fill(jetTrigPtMC);
    if(type==2)
      fh1TriggersFull[0]->Fill(jetTrigPtMC);

    AliEmcalJet *jetTrigDet = jetTrigMC->ClosestJet();
    if(!jetTrigDet) {
      //trigger is lost
      if(type==1)
	fh1TriggersLostCharged->Fill(jetTrigPtMC);
      if(type==2)
	fh1TriggersLostFull->Fill(jetTrigPtMC);
      continue;
    }

    if(type==1)
      fh1TriggersCharged[1]->Fill(GetJetPt(jetTrigDet,typet));
    if(type==2)
      fh1TriggersFull[1]->Fill(GetJetPt(jetTrigDet,typet));


    AliEmcalJet *jetAssocMC = GetLeadingJetOppositeHemisphere(type,typeaMC,jetTrigMC);
    if(!jetAssocMC) continue;

    Double_t jetAssocPtMC = GetJetPt(jetAssocMC,typeaMC);
      
    //Now check if jets are also there on detector level
    AliEmcalJet *jetAssocDet = jetAssocMC->ClosestJet();
    if(!jetAssocDet) {
      //dijet is lost
      if(type==1)
	fh3AssocLostPtDeltaPhiCharged->Fill(jetTrigPtMC,jetAssocPtMC,GetDeltaPhi(jetTrigMC,jetAssocMC));
      if(type==2)
	fh3AssocLostPtDeltaPhiFull->Fill(jetTrigPtMC,jetAssocPtMC,GetDeltaPhi(jetTrigMC,jetAssocMC));
      continue;
    }

    if(fDoPtBias) {
      if(type==0 || type==1) {
	if(GetJetPt(jetAssocDet,typea)>GetJetPt(jetTrigDet,typet))
	  continue;
      }
    }

    FillDiJetResponse(jetTrigMC,jetAssocMC,jetTrigDet,jetAssocDet,type);


  }//trigger jet loop

}

//________________________________________________________________________
void AliAnalysisTaskEmcalDiJetResponse::FillDiJetResponse(const AliEmcalJet *jetTrigMC, const AliEmcalJet *jetAssocMC, const AliEmcalJet *jetTrigDet, const AliEmcalJet *jetAssocDet, Int_t type) {

  //Fill dijet response

  Int_t typet = 0;
  Int_t typea = 0;
  Int_t typetMC = 0;
  Int_t typeaMC = 0;
  if(type==0) { //full-full
    typetMC = fContainerFullMC;
    typeaMC = fContainerFullMC;
    typet = fContainerFull;
    typea = fContainerFull;
  }
  else if(type==1) { //charged-charged
    typetMC = fContainerChargedMC;
    typeaMC = fContainerChargedMC;
    typet = fContainerCharged;
    typea = fContainerCharged;
  }
  else if(type==2) { //full-charged
    typetMC = fContainerFullMC;
    typeaMC = fContainerChargedMC;
    typet = fContainerFull;
    typea = fContainerCharged;
  }
  else {
    AliWarning(Form("%s: type %d of dijet correlation not defined!",GetName(),type));
    return;
  }

  Double_t jetTrigPtMC  = GetJetPt(jetTrigMC,typetMC);
  Double_t jetAssocPtMC = GetJetPt(jetAssocMC,typeaMC);

  Double_t varDet[2] = {TMath::Abs(GetJetPt(jetTrigDet,typet)*TMath::Sin(GetDeltaPhi(jetTrigDet,jetAssocDet))),(jetTrigDet->Eta()+jetAssocDet->Eta())/2.};
  Double_t varPart[2] = {TMath::Abs(jetTrigPtMC*TMath::Sin(GetDeltaPhi(jetTrigMC,jetAssocMC))),(jetTrigMC->Eta()+jetAssocMC->Eta())/2.};

  Double_t ajDet  = (GetJetPt(jetTrigDet,typet)-GetJetPt(jetAssocDet,typea))/(GetJetPt(jetTrigDet,typet)+GetJetPt(jetAssocDet,typea));
  Double_t ajPart = (jetTrigPtMC-jetAssocPtMC)/(jetTrigPtMC+jetAssocPtMC);

  //Store dijet vars: pt,trig MC; pt,trig DET; pt,ass MC; pt,ass DET; dPhi MC; dPhi Det; kT MC; kT Det;
  Double_t diJetVars[10] = {
    jetTrigPtMC,
    GetJetPt(jetTrigDet,typet),
    jetAssocPtMC,
    GetJetPt(jetAssocDet,typea),
    GetDeltaPhi(jetTrigMC,jetAssocMC),
    GetDeltaPhi(jetTrigDet,jetAssocDet),
    varPart[fnUsedResponseVar],
    varDet[fnUsedResponseVar],
    ajDet,
    ajPart
  }; 
  
  if(type==1)
    fhnDiJetResponseCharged->Fill(diJetVars);
  else if(type==2)
    fhnDiJetResponseFullCharged->Fill(diJetVars);


}

//________________________________________________________________________
void AliAnalysisTaskEmcalDiJetResponse::FillMatchHistos() {
  //
  // Fill Particle-Detector level matching histos
  //

  for(int i = 0; i < GetNJets(fContainerFull);++i) {
    AliEmcalJet *jetDet = static_cast<AliEmcalJet*>(GetAcceptJetFromArray(i, fContainerFull));
    if(!jetDet) continue;

    AliEmcalJet *jetPart = jetDet->ClosestJet();
    if(!jetPart) continue;

    Double_t matchVars[6] = {
      jetPart->Pt(),
      jetDet->Pt(),
      GetDeltaPhi(jetPart->Phi(),jetDet->Phi()),
      jetPart->Eta()-jetDet->Eta(),
      GetDeltaR(jetPart,jetDet),
      TMath::Min((Float_t)i+0.5,2.5)
    };
    fhnMatchingFull->Fill(matchVars);

  }//loop over full jets

  for(int i = 0; i < GetNJets(fContainerCharged);++i) {
    AliEmcalJet *jetDet = static_cast<AliEmcalJet*>(GetAcceptJetFromArray(i, fContainerCharged));
    if(!jetDet) continue;

    AliEmcalJet *jetPart = jetDet->ClosestJet();
    if(!jetPart) continue;

    Double_t matchVars[6] = {
      jetPart->Pt(),
      jetDet->Pt(),
      GetDeltaPhi(jetPart->Phi(),jetDet->Phi()),
      jetPart->Eta()-jetDet->Eta(),
      GetDeltaR(jetPart,jetDet),
      TMath::Min((Float_t)i+0.5,2.5)
    };
    fhnMatchingCharged->Fill(matchVars);

  }//loop over charged jets

}


//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalDiJetResponse::RetrieveEventObjects() {
  //
  // retrieve event objects
  //

  if (!AliAnalysisTaskEmcalDiJetBase::RetrieveEventObjects())
    return kFALSE;

  return kTRUE;

}

//_______________________________________________________________________
void AliAnalysisTaskEmcalDiJetResponse::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
}
