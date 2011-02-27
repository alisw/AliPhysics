#include "TChain.h"
#include "TTree.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
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
  fEvtClassMin(1),
  fEvtClassMax(4),
  fCentMin(0.),
  fCentMax(100.),
  fJetEtaMin(-.5),
  fJetEtaMax(.5),
  fJetPtMin(20.),
  fJetPtFractionMin(0.5),
  fNMatchJets(4),
  //fJetDeltaEta(.4),
  //fJetDeltaPhi(.4),
  fkNbranches(2),
  fkEvtClasses(5),
  fOutputList(0x0),
  fHistEvtSelection(0x0),
  fHistEvtClass(0x0),
  fHistCentrality(0x0),
  fHistPtJet(new TH1F*[fkNbranches]),
  fHistEtaPhiJet(new TH2F*[fkNbranches]),
  fHistEtaPhiJetCut(new TH2F*[fkNbranches]),
  fHistDeltaEtaDeltaPhiJet(new TH2F*[fkEvtClasses]),
  fHistDeltaEtaDeltaPhiJetCut(new TH2F*[fkEvtClasses]),
  fHistDeltaEtaDeltaPhiJetNOMatching(new TH2F*[fkEvtClasses]),
  fHistDeltaEtaEtaJet(new TH2F*[fkEvtClasses]),
  fHistDeltaPtEtaJet(new TH2F*[fkEvtClasses]),
  fHistPtFraction(new TH2F*[fkEvtClasses]),
  fHistPtPtExtra(0x0),
  fHistPtResponse(new TH2F*[fkEvtClasses]),
  fHistPtSmearing(new TH2F*[fkEvtClasses]),
  fHistDeltaR(new TH2F*[fkEvtClasses]),
  fHistArea(new TH2F*[fkEvtClasses]),
  fHistDeltaArea(new TH2F*[fkEvtClasses])
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
  fEvtClassMin(1),
  fEvtClassMax(4),
  fCentMin(0.),
  fCentMax(100.),
  fJetEtaMin(-.5),
  fJetEtaMax(.5),
  fJetPtMin(20.),
  fJetPtFractionMin(0.5),
  fNMatchJets(4),
  //fJetDeltaEta(.4),
  //fJetDeltaPhi(.4),
  fkNbranches(2),
  fkEvtClasses(5),
  fOutputList(0x0),
  fHistEvtSelection(0x0),
  fHistEvtClass(0x0),
  fHistCentrality(0x0),
  fHistPtJet(new TH1F*[fkNbranches]),
  fHistEtaPhiJet(new TH2F*[fkNbranches]),
  fHistEtaPhiJetCut(new TH2F*[fkNbranches]),
  fHistDeltaEtaDeltaPhiJet(new TH2F*[fkEvtClasses]),
  fHistDeltaEtaDeltaPhiJetCut(new TH2F*[fkEvtClasses]),
  fHistDeltaEtaDeltaPhiJetNOMatching(new TH2F*[fkEvtClasses]),
  fHistDeltaEtaEtaJet(new TH2F*[fkEvtClasses]),
  fHistDeltaPtEtaJet(new TH2F*[fkEvtClasses]),
  fHistPtFraction(new TH2F*[fkEvtClasses]),
  fHistPtPtExtra(0x0),
  fHistPtResponse(new TH2F*[fkEvtClasses]),
  fHistPtSmearing(new TH2F*[fkEvtClasses]),
  fHistDeltaR(new TH2F*[fkEvtClasses]),
  fHistArea(new TH2F*[fkEvtClasses]),
  fHistDeltaArea(new TH2F*[fkEvtClasses])
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
  OpenFile(1);
  if(!fOutputList) fOutputList = new TList;
  fOutputList->SetOwner(kTRUE);

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);


  fHistEvtSelection = new TH1I("fHistEvtSelection", "event selection", 5, -0.5, 4.5);
  fHistEvtSelection->GetXaxis()->SetBinLabel(1,"ACCEPTED");
  fHistEvtSelection->GetXaxis()->SetBinLabel(2,"events IN");
  fHistEvtSelection->GetXaxis()->SetBinLabel(3,"event selection (rejected)");
  fHistEvtSelection->GetXaxis()->SetBinLabel(4,"vertex cut (rejected)");
  fHistEvtSelection->GetXaxis()->SetBinLabel(5,"centrality (rejected)");
  
  fHistEvtClass   = new TH1I("fHistEvtClass", "event classes", fkEvtClasses, -0.5, fkEvtClasses-0.5);
  fHistCentrality = new TH1F("fHistCentrality", "event centrality", 100, 0., 100.);

  for (Int_t iBranch = 0; iBranch < fkNbranches; iBranch++) {
    fHistPtJet[iBranch] = new TH1F(Form("fHistPtJet%i", iBranch),
					  Form("p_{T} of jets from branch %i;p_{T} (GeV/c);counts", iBranch),
					  250, 0., 250.);
    fHistPtJet[iBranch]->SetMarkerStyle(kFullCircle);

    fHistEtaPhiJet[iBranch]    = new TH2F(Form("fHistEtaPhiJet%i", iBranch),
					                 Form("#eta - #phi of jets from branch %i (before cut);#eta;#phi", iBranch),
					                 48, -1.2, 1.2, 100, 0., 2*TMath::Pi());
	fHistEtaPhiJetCut[iBranch] = new TH2F(Form("fHistEtaPhiJetCut%i", iBranch),
					                 Form("#eta - #phi of jets from branch %i (before cut);#eta;#phi", iBranch),
					                 48, -1.2, 1.2, 100, 0., 2*TMath::Pi());

  }
  

  fHistPtPtExtra = new TH2F("fHistPtPtExtra", "p_{T} response;p_{T} (GeV/c);p_{T} (GeV/c)",
			    250, 0., 250., 250, 0., 250.);

  for (Int_t iEvtClass =0; iEvtClass < fkEvtClasses; iEvtClass++) {
    fHistDeltaEtaDeltaPhiJet[iEvtClass] = new TH2F(Form("fHistDeltaEtaDeltaPhiJet%i",iEvtClass),
					                              "#Delta#eta - #Delta#phi of jet;#Delta#eta;#Delta#phi",
					                              101, -1.01, 1.01, 101, -1.01, 1.01);
    fHistDeltaEtaDeltaPhiJetCut[iEvtClass] = new TH2F(Form("fHistDeltaEtaDeltaPhiJetCut%i",iEvtClass),
					                              "#Delta#eta - #Delta#phi of jet;#Delta#eta;#Delta#phi",
					                              101, -1.01, 1.01, 101, -1.01, 1.01);												  
    fHistDeltaEtaDeltaPhiJetNOMatching[iEvtClass] = new TH2F("fHistDeltaEtaDeltaPhiJetNOMatching",
					                                         "#Delta#eta - #Delta#phi of jets which do not match;#Delta#eta;#Delta#phi",
					                                         100, -2., 2., 100, TMath::Pi(), TMath::Pi());
    fHistPtResponse[iEvtClass] = new TH2F(Form("pt_response%i",iEvtClass), Form("pt_response%i;p_{T} (GeV/c);p_{T} (GeV/c)",iEvtClass),
					                      250, 0., 250., 250, 0., 250.);
    fHistPtSmearing[iEvtClass] = new TH2F(Form("pt_smearing%i",iEvtClass), Form("pt_smearing%i;#Deltap_{T} (GeV/c);p_{T} (GeV/c)",iEvtClass),
					                      200, -50., 150., 250, 0., 250.);
    fHistDeltaR[iEvtClass]     = new TH2F(Form("hist_DeltaR%i",iEvtClass), "#DeltaR of matched jets;#Deltap_{T};#DeltaR", 200, -50.,150., 60, 0.,.6);
    fHistArea[iEvtClass]       = new TH2F(Form("hist_Area%i",iEvtClass), "jet area;#Deltap_{T};jet area", 200, -50.,150., 100, 0.,1.0);
    fHistDeltaArea[iEvtClass]  = new TH2F(Form("hist_DeltaArea%i",iEvtClass), "#DeltaArea of matched jets;#Deltap_{T};#Deltaarea", 200, -50., 150., 60, 0., .3); 
	
    fHistDeltaEtaEtaJet[iEvtClass] = new TH2F(Form("fHistDeltaEtaEtaJet%i", iEvtClass),
                                                Form("#eta - #Delta#eta of matched jets from event class %i;#eta;#Delta#eta", iEvtClass),
                                                60, -.6, .6, 100, -.5, .5);
    fHistDeltaPtEtaJet[iEvtClass] = new TH2F(Form("fHistDeltaPtEtaJet%i", iEvtClass),
                                                Form("#eta - #Deltap_{T} of matched jets from event class %i;#eta;#Deltap_{T}", iEvtClass),
                                                60, -.6, .6, 200, -50., 150);
    fHistPtFraction[iEvtClass] = new TH2F(Form("fHistPtFraction%i", iEvtClass),
	                                      Form("fraction from embedded jet in reconstructed jet;fraction (event class %i);p_{T}^{emb}", iEvtClass),
										  100, 0., 1., 250, 0., 250.);
  }

  
  fOutputList->Add(fHistEvtSelection);
  fOutputList->Add(fHistEvtClass);
  fOutputList->Add(fHistCentrality);

  for (Int_t iBranch = 0; iBranch < fkNbranches; iBranch++) {
    fOutputList->Add(fHistPtJet[iBranch]);
    fOutputList->Add(fHistEtaPhiJet[iBranch]);
	fOutputList->Add(fHistEtaPhiJetCut[iBranch]);
  }
  
  fOutputList->Add(fHistPtPtExtra);
  
  for (Int_t iEvtClass =0; iEvtClass < fkEvtClasses; iEvtClass++) {
    fOutputList->Add(fHistDeltaEtaDeltaPhiJet[iEvtClass]);
    fOutputList->Add(fHistDeltaEtaDeltaPhiJetCut[iEvtClass]);
    fOutputList->Add(fHistDeltaEtaDeltaPhiJetNOMatching[iEvtClass]);
	fOutputList->Add(fHistDeltaEtaEtaJet[iEvtClass]);
    fOutputList->Add(fHistDeltaPtEtaJet[iEvtClass]);
	fOutputList->Add(fHistPtFraction[iEvtClass]);
	
    fOutputList->Add(fHistPtResponse[iEvtClass]);
    fOutputList->Add(fHistPtSmearing[iEvtClass]);
    fOutputList->Add(fHistDeltaR[iEvtClass]);
    fOutputList->Add(fHistArea[iEvtClass]);
    fOutputList->Add(fHistDeltaArea[iEvtClass]);
  }

  // =========== Switch on Sumw2 for all histos ===========
  for (Int_t i=0; i<fOutputList->GetEntries(); ++i) {
    TH1 *h1 = dynamic_cast<TH1*>(fOutputList->At(i));
    if (h1){
      h1->Sumw2();
      continue;
    }
  }
  TH1::AddDirectory(oldStatus);

  PostData(1, fOutputList);
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

  // -- event selection --
  fHistEvtSelection->Fill(1); // number of events before event selection

  // physics selection
  AliInputEventHandler* inputHandler = (AliInputEventHandler*)
    ((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
  if(!(inputHandler->IsEventSelected() & fOfflineTrgMask)){
    if(fDebug) Printf(" Trigger Selection: event REJECTED ... ");
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
    if(fDebug) Printf("%s:%d primary vertex z = %f: event REJECTED...",(char*)__FILE__,__LINE__,primVtx->GetZ());
    fHistEvtSelection->Fill(3);
    PostData(1, fOutputList);
    return;
  }

  // event class selection (from jet helper task)
  Int_t eventClass = AliAnalysisHelperJetTasks::EventClass();
  if(fDebug) Printf("Event class %d", eventClass);
  if (eventClass < fEvtClassMin || eventClass > fEvtClassMax){
    fHistEvtSelection->Fill(4);
    PostData(1, fOutputList);
    return;
  }

  // centrality selection
  AliCentrality *cent = fESD->GetCentrality();
  Float_t centValue = cent->GetCentralityPercentile("TRK");
  if(fDebug) printf("centrality: %f\n", centValue);
  if (centValue < fCentMin || centValue > fCentMax){
    fHistEvtSelection->Fill(4);
    PostData(1, fOutputList);
    return;
  }

  fHistEvtSelection->Fill(0); // accepted events
  fHistEvtClass->Fill(eventClass);
  fHistCentrality->Fill(centValue);
  // -- end event selection --

  // fetch jets
  TClonesArray *aodJets[2];
  aodJets[0] = dynamic_cast<TClonesArray*>(fAOD->FindListObject(fJetBranchName[0].Data())); // in general: embedded jet
  aodJets[1] = dynamic_cast<TClonesArray*>(fAOD->FindListObject(fJetBranchName[1].Data())); // in general: embedded jet + UE

  for (Int_t iJetType = 0; iJetType < 2; iJetType++) {
    fListJets[iJetType]->Clear();
    if (!aodJets[iJetType]) continue;
	
    for (Int_t iJet = 0; iJet < aodJets[iJetType]->GetEntriesFast(); iJet++) {
      AliAODJet *jet = dynamic_cast<AliAODJet*>((*aodJets[iJetType])[iJet]);
      if (jet) fListJets[iJetType]->Add(jet);
    }
    fListJets[iJetType]->Sort();
  }
  
  // jet matching 
  static TArrayI aMatchIndex(fListJets[0]->GetEntries());
  static TArrayF aPtFraction(fListJets[0]->GetEntries());
  if(aMatchIndex.GetSize()<fListJets[0]->GetEntries()) aMatchIndex.Set(fListJets[0]->GetEntries());
  if(aPtFraction.GetSize()<fListJets[0]->GetEntries()) aPtFraction.Set(fListJets[0]->GetEntries());
  
  // stores matched jets in 'aMatchIndex' and fraction of pT in 'aPtFraction'
  AliAnalysisHelperJetTasks::GetJetMatching(fListJets[0], TMath::Min((Int_t)fNMatchJets,(Int_t)fListJets[0]->GetEntries()),
                                            fListJets[1], TMath::Min((Int_t)fNMatchJets,(Int_t)fListJets[1]->GetEntries()),
                                            aMatchIndex, aPtFraction, fDebug);
											
  // loop over matched jets
  Int_t ir = -1; // index of matched reconstruced jet
  Float_t fraction = 0.;
  AliAODJet *jet[2]  = { 0x0, 0x0 };
  Float_t jetEta[2]  = { 0., 0. };
  Float_t jetPhi[2]  = { 0., 0. };
  Float_t jetPt[2]   = { 0., 0. };
  Float_t jetArea[2] = { 0., 0. };
   
  for(Int_t ig=0; ig<fListJets[0]->GetEntries(); ++ig){
     ir = aMatchIndex[ig];
	 if(ir<0) continue;
	 fraction = aPtFraction[ig];
	 
	 // fetch jets
	 jet[0] = (AliAODJet*)(fListJets[0]->At(ig));
	 jet[1] = (AliAODJet*)(fListJets[1]->At(ir));
	 if(!jet[0] || !jet[1]) continue;
	 
	 for(Int_t i=0; i<fkNbranches; ++i){
	    jetEta[i]  = jet[i]->Eta();
		jetPhi[i]  = jet[i]->Phi();
		jetPt[i]   = jet[i]->Pt();
		jetArea[i] = jet[i]->EffectiveAreaCharged();
	 }
	 if(eventClass > -1 && eventClass < fkEvtClasses){
        fHistPtFraction[eventClass] -> Fill(fraction, jetPt[0]);
	 }
	 
	 if(fraction<fJetPtFractionMin) continue;
	 
	 // calculate parameters of associated jets
	 Float_t deltaPt    = jetPt[1]-jetPt[0];
	 Float_t deltaEta   = jetEta[1]-jetEta[0];
	 Float_t deltaPhi   = TVector2::Phi_mpi_pi(jetPhi[1]-jetPhi[0]);
	 Float_t deltaR     = TMath::Sqrt(deltaEta*deltaEta+deltaPhi*deltaPhi);
	 Float_t deltaArea  = jetArea[1]-jetArea[0];
	 
	 // fill histograms before acceptance cut
	 for(Int_t i=0; i<fkNbranches; ++i){
	     fHistEtaPhiJet[i] -> Fill(jetEta[i], jetPhi[i]);
         }
	 if(eventClass > -1 && eventClass < fkEvtClasses){
             fHistDeltaEtaDeltaPhiJet[eventClass] -> Fill(deltaEta, deltaPhi);
	 }
	 
	 // jet acceptance + minimum pT check
	 if(jetEta[0]>fJetEtaMax || jetEta[0]<fJetEtaMin ||
	    jetEta[1]>fJetEtaMax || jetEta[1]<fJetEtaMin){
		
		if(fDebug){
     		Printf("Jet not in eta acceptance.");
			Printf("[0]: jet %d eta %.2f", ig, jetEta[0]);
			Printf("[1]: jet %d eta %.2f", ir, jetEta[1]);
		}
		continue;
     }
	 if(jetPt[1] < fJetPtMin){
	    if(fDebug) Printf("Jet %d (pT %.1f GeV/c) has less than required pT.", ir, jetPt[1]);
	    continue;
	 }
	 
	 
	 
	 // fill histograms
	 for(Int_t i=0; i<fkNbranches; ++i){
	     fHistPtJet[i]        -> Fill(jetPt[i]);
             fHistEtaPhiJetCut[i] -> Fill(jetEta[i], jetPhi[i]);
	 }
	 
	 fHistPtPtExtra->Fill(jetPt[0], jetPt[1]);
	 
	 if(eventClass > -1 && eventClass < fkEvtClasses){
		 fHistDeltaEtaDeltaPhiJetCut[eventClass] -> Fill(deltaEta, deltaPhi);
		 
		 fHistPtResponse[eventClass]             -> Fill(jetPt[1], jetPt[0]);
		 fHistPtSmearing[eventClass]             -> Fill(deltaPt,  jetPt[0]);
		 
		 fHistDeltaPtEtaJet[eventClass]          -> Fill(jetEta[0], deltaPt);
                 fHistDeltaEtaEtaJet[eventClass]         -> Fill(jetEta[0], deltaEta);
		 
		 fHistDeltaR[eventClass]                 -> Fill(deltaPt, deltaR);
		 fHistArea[eventClass]                   -> Fill(deltaPt, jetArea[1]);
                 fHistDeltaArea[eventClass]              -> Fill(deltaPt, deltaArea);
		 
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

