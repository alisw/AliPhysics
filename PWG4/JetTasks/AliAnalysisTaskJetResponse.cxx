#include "TChain.h"
#include "TTree.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"

#include "AliLog.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliCentrality.h"
#include "AliAnalysisHelperJetTasks.h"
#include "AliInputEventHandler.h"

#include "AliAODEvent.h"
#include "AliAODJet.h"

#include "AliAnalysisTaskJetResponse.h"

ClassImp(AliAnalysisTaskJetResponse)

AliAnalysisTaskJetResponse::AliAnalysisTaskJetResponse() :
  AliAnalysisTaskSE(),
  fESD(0x0),
  fAOD(0x0),
  fOfflineTrgMask(AliVEvent::kAny),
  fMinContribVtx(1),
  fVtxZMin(-8.),
  fVtxZMax(8.),
  fEvtClassMin(0),
  fEvtClassMax(10),
  fCentMin(0.),
  fCentMax(100.),
  fJetEtaMin(-.5),
  fJetEtaMax(.5),
  fJetDeltaEta(.4),
  fJetDeltaPhi(.4),
  fkNbranches(2),
  fkEvtClasses(10),
  fOutputList(0x0),
  fHistEvtSelection(0x0),
  fHistPtLeadingJet(new TH1F*[fkNbranches]),
  fHistEtaPhiLeadingJet(new TH2F*[fkNbranches]),
  fHistDeltaEtaDeltaPhiLeadingJet(0x0),
  fHistPtPtExtra(0x0),
  fHistPtResponse(new TH2F*[fkEvtClasses]),
  fHistPtSmearing(new TH2F*[fkEvtClasses])
{
  // default Constructor

  fJetBranchName[0] = "";
  fJetBranchName[1] = "";

  fListJets[0] = new TList;
  fListJets[1] = new TList;
}

AliAnalysisTaskJetResponse::AliAnalysisTaskJetResponse(const char *name) :
  AliAnalysisTaskSE(name),
  fESD(0x0),
  fAOD(0x0),
  fOfflineTrgMask(AliVEvent::kAny),
  fMinContribVtx(1),
  fVtxZMin(-8.),
  fVtxZMax(8.),
  fEvtClassMin(0),
  fEvtClassMax(10),
  fCentMin(0.),
  fCentMax(100.),
  fJetEtaMin(-.5),
  fJetEtaMax(.5),
  fJetDeltaEta(.4),
  fJetDeltaPhi(.4),
  fkNbranches(2),
  fkEvtClasses(10),
  fOutputList(0x0),
  fHistEvtSelection(0x0),
  fHistPtLeadingJet(new TH1F*[fkNbranches]),
  fHistEtaPhiLeadingJet(new TH2F*[fkNbranches]),
  fHistDeltaEtaDeltaPhiLeadingJet(0x0),
  fHistPtPtExtra(0x0),
  fHistPtResponse(new TH2F*[fkEvtClasses]),
  fHistPtSmearing(new TH2F*[fkEvtClasses])
{
  // Constructor

  fJetBranchName[0] = "";
  fJetBranchName[1] = "";

  fListJets[0] = new TList;
  fListJets[1] = new TList;

  DefineOutput(1, TList::Class());
}

AliAnalysisTaskJetResponse::~AliAnalysisTaskJetResponse()
{
  delete fListJets[0];
  delete fListJets[1];
}

void AliAnalysisTaskJetResponse::SetBranchNames(const TString &branch1, const TString &branch2)
{
  fJetBranchName[0] = branch1;
  fJetBranchName[1] = branch2;
}

void AliAnalysisTaskJetResponse::Init()
{

   // check for jet branches
   if(!strlen(fJetBranchName[0].Data()) || !strlen(fJetBranchName[1].Data())){
      AliError("Jet branch name not set.");
   }

}

void AliAnalysisTaskJetResponse::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

  fHistEvtSelection = new TH1F("fHistEvtSelection", "event selection", 5, -0.5, 5.5);
  fHistEvtSelection->GetXaxis()->SetBinLabel(1,"ACCEPTED");
  fHistEvtSelection->GetXaxis()->SetBinLabel(2,"events IN");
  fHistEvtSelection->GetXaxis()->SetBinLabel(3,"event selection (rejected)");
  fHistEvtSelection->GetXaxis()->SetBinLabel(4,"vertex cut (rejected)");
  fHistEvtSelection->GetXaxis()->SetBinLabel(5,"centrality (rejected)");

  for (Int_t iBranch = 0; iBranch < fkNbranches; iBranch++) {
    fHistPtLeadingJet[iBranch] = new TH1F(Form("fHistPtLeadingJet%i", iBranch),
					  Form("p_{T} of leading jet from branch %i;p_{T} (GeV/c);counts", iBranch),
					  300, 0., 300.);
    fHistPtLeadingJet[iBranch]->SetMarkerStyle(kFullCircle);

    fHistEtaPhiLeadingJet[iBranch] = new TH2F(Form("fHistEtaPhiLeadingJet%i", iBranch),
					      Form("#eta - #phi of leading jet from branch %i;#eta;#phi", iBranch),
					      100, -2., 2., 100, 0., 2*TMath::Pi());

  }
  fHistDeltaEtaDeltaPhiLeadingJet = new TH2F("fHistDeltaEtaDeltaPhiLeadingJet",
					     "#Delta#eta - #Delta#phi of leading jet;#Delta#eta;#Delta#phi",
					     100, -4., 4., 100, -2.*TMath::Pi(), 2*TMath::Pi());

  fHistPtPtExtra = new TH2F("fHistPtPtExtra", "p_{T} response;p_{T} (GeV/c);p_{T} (GeV/c)",
			    300, 0., 300., 300, 0., 300.);

  for (Int_t iEvtClass =0; iEvtClass < fkEvtClasses; iEvtClass++) {
    fHistPtResponse[iEvtClass] = new TH2F(Form("pt_response%i",iEvtClass), Form("pt_response%i;p_{T} (GeV/c);p_{T} (GeV/c)",iEvtClass),
					  300, 0., 300., 300, 0., 300.);
    fHistPtSmearing[iEvtClass] = new TH2F(Form("pt_smearing%i",iEvtClass), Form("pt_smearing%i;p_{T} (GeV/c);p_{T} (GeV/c)",iEvtClass),
					  200, -50., 150., 300, 0., 300.);
  }

  OpenFile(1);
  fOutputList = new TList;
  fOutputList->Add(fHistEvtSelection);
  for (Int_t iBranch = 0; iBranch < fkNbranches; iBranch++) {
    fOutputList->Add(fHistPtLeadingJet[iBranch]);
    fOutputList->Add(fHistEtaPhiLeadingJet[iBranch]);
  }
  fOutputList->Add(fHistDeltaEtaDeltaPhiLeadingJet);
  fOutputList->Add(fHistPtPtExtra);
  for (Int_t iEvtClass =0; iEvtClass < fkEvtClasses; iEvtClass++) {
    fOutputList->Add(fHistPtResponse[iEvtClass]);
    fOutputList->Add(fHistPtSmearing[iEvtClass]);
  }
}

void AliAnalysisTaskJetResponse::UserExec(Option_t *)
{
  // load events, apply event cuts, then compare leading jets

  if(!strlen(fJetBranchName[0].Data()) || !strlen(fJetBranchName[1].Data())){
      AliError("Jet branch name not set.");
      return;
   }

  fESD=dynamic_cast<AliESDEvent*>(InputEvent());
  if (!fESD) {
    AliError("ESD not available");
    return;
  }
  fAOD = dynamic_cast<AliAODEvent*>(AODEvent());
  if (!fAOD) {
    AliError("AOD not available");
    return;
  }

  fHistEvtSelection->Fill(1); // number of events before event selection

  // physics selection
  AliInputEventHandler* inputHandler = (AliInputEventHandler*)
    ((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
  if(!(inputHandler->IsEventSelected() & fOfflineTrgMask)){
    Printf(" Trigger Selection: event REJECTED ... ");
    fHistEvtSelection->Fill(2);
    PostData(1, fOutputList);
    return;
  }

  // vertex selection
  AliAODVertex* primVtx = fAOD->GetPrimaryVertex();
  Int_t nTracksPrim = primVtx->GetNContributors();
  if ((nTracksPrim < fMinContribVtx) ||
      (primVtx->GetZ() < fVtxZMin) ||
      (primVtx->GetZ() > fVtxZMax) ){
    Printf("%s:%d primary vertex z = %f: event REJECTED...",(char*)__FILE__,__LINE__,primVtx->GetZ());
    fHistEvtSelection->Fill(3);
    PostData(1, fOutputList);
    return;
  }

  // event class selection (from jet helper task)
  Int_t eventClass = AliAnalysisHelperJetTasks::EventClass();
  Printf("Event class %d", eventClass);
  if (eventClass < fEvtClassMin || eventClass > fEvtClassMax){
    fHistEvtSelection->Fill(4);
    PostData(1, fOutputList);
    return;
  }

  // centrality selection
  AliCentrality *cent = fESD->GetCentrality();
  Float_t centValue = cent->GetCentralityPercentile("TRK");
  printf("centrality: %f\n", centValue);
  if (centValue < fCentMin || centValue > fCentMax){
    fHistEvtSelection->Fill(5);
    PostData(1, fOutputList);
    return;
  }

  fHistEvtSelection->Fill(0); // accepted events


  TClonesArray *aodJets[2];
  aodJets[0] = dynamic_cast<TClonesArray*>(fAOD->FindListObject(fJetBranchName[0].Data()));
  aodJets[1] = dynamic_cast<TClonesArray*>(fAOD->FindListObject(fJetBranchName[1].Data()));
  AliAODJet *leadingJet[2] = { 0x0, 0x0 };

  for (Int_t iJetType = 0; iJetType < 2; iJetType++) {
    fListJets[iJetType]->Clear();
    if (!aodJets[iJetType])
      continue;
    for (Int_t iJet = 0; iJet < aodJets[iJetType]->GetEntriesFast(); iJet++) {
      AliAODJet *jet = dynamic_cast<AliAODJet*>((*aodJets[iJetType])[iJet]);
      if (jet)
	fListJets[iJetType]->Add(jet);
    }
    fListJets[iJetType]->Sort();
    leadingJet[iJetType] = (AliAODJet*) fListJets[iJetType]->First();

    if (leadingJet[iJetType])
      fHistEtaPhiLeadingJet[iJetType]->Fill(leadingJet[iJetType]->Eta(), leadingJet[iJetType]->Phi());
  }

  // check if leading jets are close in eta-phi and compare their Pt
  if (leadingJet[0] && leadingJet[1]) {

    // leading jets in "jet" acceptance
    if(leadingJet[0]->Eta()>fJetEtaMax || leadingJet[0]->Eta()<fJetEtaMin ||
       leadingJet[1]->Eta()>fJetEtaMax || leadingJet[1]->Eta()<fJetEtaMin){
       Printf("Jet not in eta acceptance.");
    }
    else{
       // check association of jets
       Float_t deltaEta = leadingJet[0]->Eta() - leadingJet[1]->Eta();
       Float_t deltaPhi = leadingJet[0]->Phi() - leadingJet[1]->Phi();
       fHistDeltaEtaDeltaPhiLeadingJet->Fill(deltaEta, deltaPhi);

       if (TMath::Abs(deltaEta) > fJetDeltaEta || (TMath::Pi() - TMath::Abs(TMath::Abs(deltaPhi) - TMath::Pi())) > fJetDeltaPhi)
          printf("leading jets two far apart\n");
       else {
          fHistPtLeadingJet[0]->Fill(leadingJet[0]->Pt());
          fHistPtLeadingJet[1]->Fill(leadingJet[1]->Pt());

          fHistPtPtExtra->Fill(leadingJet[0]->Pt(), leadingJet[1]->Pt());

          if (eventClass > -1 && eventClass < fkEvtClasses){
	      fHistPtResponse[eventClass]->Fill(leadingJet[1]->Pt(), leadingJet[0]->Pt());
              fHistPtSmearing[eventClass]->Fill(leadingJet[1]->Pt()-leadingJet[0]->Pt(), leadingJet[0]->Pt());
          }
       }
    }
  }

  PostData(1, fOutputList);
}

void AliAnalysisTaskJetResponse::Terminate(const Option_t *)
{
  // Draw result to the screen
  // Called once at the end of the query

  if (!GetOutputData(1))
    return;
}
