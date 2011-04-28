#include "TChain.h"
#include "TTree.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THnSparse.h"
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

#include "AliAnalysisTaskJetResponseV2.h"

ClassImp(AliAnalysisTaskJetResponseV2)

AliAnalysisTaskJetResponseV2::AliAnalysisTaskJetResponseV2() :
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
  fNInputTracksMin(0),
  fNInputTracksMax(-1),
  fReactionPlaneBin(-1),
  fJetEtaMin(-.5),
  fJetEtaMax(.5),
  fJetPtMin(20.),
  fJetPtFractionMin(0.5),
  fNMatchJets(4),
  fMatchMaxDist(0.8),
  fkNbranches(2),
  fkEvtClasses(12),
  fOutputList(0x0),
  fHistEvtSelection(0x0),
  fhnEvent(0x0),
  fhnJetsRp(0x0),
  fhnJetsDeltaPt(0x0),
  fhnJetsEta(0x0),
  fhnJetsArea(0x0),
  fhnJetsBeforeCut1(0x0),
  fhnJetsBeforeCut2(0x0)
{
  // default Constructor

  fJetBranchName[0] = "";
  fJetBranchName[1] = "";

  fListJets[0] = new TList;
  fListJets[1] = new TList;
}

AliAnalysisTaskJetResponseV2::AliAnalysisTaskJetResponseV2(const char *name) :
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
  fNInputTracksMin(0),
  fNInputTracksMax(-1),
  fReactionPlaneBin(-1),
  fJetEtaMin(-.5),
  fJetEtaMax(.5),
  fJetPtMin(20.),
  fJetPtFractionMin(0.5),
  fNMatchJets(4),
  fMatchMaxDist(0.8),
  fkNbranches(2),
  fkEvtClasses(12),
  fOutputList(0x0),
  fHistEvtSelection(0x0),
  fhnEvent(0x0),
  fhnJetsRp(0x0),
  fhnJetsDeltaPt(0x0),
  fhnJetsEta(0x0),
  fhnJetsArea(0x0),
  fhnJetsBeforeCut1(0x0),
  fhnJetsBeforeCut2(0x0)
{
  // Constructor

  fJetBranchName[0] = "";
  fJetBranchName[1] = "";

  fListJets[0] = new TList;
  fListJets[1] = new TList;

  DefineOutput(1, TList::Class());
}

AliAnalysisTaskJetResponseV2::~AliAnalysisTaskJetResponseV2()
{
  delete fListJets[0];
  delete fListJets[1];
}

void AliAnalysisTaskJetResponseV2::SetBranchNames(const TString &branch1, const TString &branch2)
{
  fJetBranchName[0] = branch1;
  fJetBranchName[1] = branch2;
}

void AliAnalysisTaskJetResponseV2::Init()
{

   // check for jet branches
   if(!strlen(fJetBranchName[0].Data()) || !strlen(fJetBranchName[1].Data())){
      AliError("Jet branch name not set.");
   }

}

void AliAnalysisTaskJetResponseV2::UserCreateOutputObjects()
{
  // Create histograms
  // Called once
  OpenFile(1);
  if(!fOutputList) fOutputList = new TList;
  fOutputList->SetOwner(kTRUE);

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);


  fHistEvtSelection = new TH1I("fHistEvtSelection", "event selection", 6, -0.5, 5.5);
  fHistEvtSelection->GetXaxis()->SetBinLabel(1,"ACCEPTED");
  fHistEvtSelection->GetXaxis()->SetBinLabel(2,"events IN");
  fHistEvtSelection->GetXaxis()->SetBinLabel(3,"event selection (rejected)");
  fHistEvtSelection->GetXaxis()->SetBinLabel(4,"vertex cut (rejected)");
  fHistEvtSelection->GetXaxis()->SetBinLabel(5,"centrality (rejected)");
  fHistEvtSelection->GetXaxis()->SetBinLabel(6,"multiplicity (rejected)");
  
  Float_t pi = TMath::Pi();
  
  const Int_t dim1 = 3;
  // cent : nInpTrks : rp
  Int_t    nbins1[dim1] = { 100,   400,  30  };
  Double_t xmin1[dim1]  = {   0.,    0.,  0. };
  Double_t xmax1[dim1]  = { 100., 4000., pi  };
  
  TString hnTitle("variables per event;centrality;nb. of input tracks;reaction plane #phi");
  
  fhnEvent = new THnSparseF("fhnEvent", hnTitle.Data(), dim1, nbins1, xmin1, xmax1);
  
  /* original thnsparse
  const Int_t dim2 = 18;
  // cent : nInpTrks : rp bins : rp wrt jet : fraction : 
  // jetPt(2x) : jetEta(2x) : jetPhi (2x) : jetArea(2x)
  // deltaPt : deltaEta : deltaPhi : deltaR : deltaArea
  
  Int_t    nbins2[dim2] = {  16,   400,    3,  48,    52, 125,  125,    56,   56,   25,   25, 100, 100,    201,   101,    51,  50,   81 };
  Double_t xmin2[dim2]  = {   0.,    0., -.5, -pi,  0.  ,   0.,   0., -0.7, -0.7,   0.,   0., 0. , 0. , -100.5, -1.01, -1.02,  0., -.81 };
  Double_t xmax2[dim2]  = {  80., 4000., 2.5,  pi,  1.04, 250., 250.,  0.7,  0.7, 2*pi, 2*pi, 1.0, 1.0,  100.5,  1.01,  1.02,  1.,  .81 };
  
  fhnJets = new THnSparseF("fhnJets", "variables per jet", dim2, nbins2, xmin2, xmax2);
  */
  
  
  const Int_t dim2 = 6;
  // cent : nInpTrks : rp bins : rp wrt jet : probe pT
  Int_t    nbins2[dim2] = {   8 ,    40,   3,  48,  48,  50 };
  Double_t xmin2[dim2]  = {   0.,    0., -.5, -pi, -pi,  0. };
  Double_t xmax2[dim2]  = {  80., 4000., 2.5,  pi,  pi, 250.};
  hnTitle = "variables per jet;centrality;nb. of input tracks; reaction plane bin;#Delta#phi(RP-jet);probe p_{T}";
  fhnJetsRp = new THnSparseF("fhnJetsRp", hnTitle.Data(), dim2, nbins2, xmin2, xmax2);
  
  
  const Int_t dim3 = 6;
  // cent : nInpTrks : rp bins:  deltaPt : jetPt(2x) (hr delta pt)
  Int_t    nbins3[dim3] = {  16,   400,    3,    241,  250,  250 };
  Double_t xmin3[dim3]  = {   0.,    0., -.5, -120.5,   0.,   0. };
  Double_t xmax3[dim3]  = {  80., 4000., 2.5,  120.5, 250., 250. };
  hnTitle = "variables per jet;centrality;nb. of input tracks; reaction plane bin;#delta p_{T};probe p_{T};rec p_{T}";
  fhnJetsDeltaPt = new THnSparseF("fhnJetsDeltaPt", hnTitle.Data(), dim3, nbins3, xmin3, xmax3);
    
  const Int_t dim4 = 10;	
  // cent : nInpTrks : rp bins : deltaPt : jetPt(2x) : deltaR : deltaEta : jetEta(2x) (hr for eta)
  Int_t    nbins4[dim4] = {   8 ,   40 ,   3,  101 ,  50 ,  50 , 50,    51,   56,   56 };
  Double_t xmin4[dim4]  = {   0.,    0., -.5, -101.,   0.,   0., 0., -1.02, -0.7, -0.7 };
  Double_t xmax4[dim4]  = {  80., 4000., 2.5,  101., 250., 250., 1.,  1.02,  0.7,  0.7 };
  hnTitle = "variables per jet;centrality;nb. of input tracks; reaction plane bin;#delta p_{T};probe p_{T};rec p_{T};#DeltaR;#Delta#eta;#eta(probe);#eta(rec)";
  fhnJetsEta = new THnSparseF("fhnJetsEta", hnTitle.Data(), dim4, nbins4, xmin4, xmax4);
  
  const Int_t dim5 = 13; 
  // cent : nInpTrks : rp bins : deltaArea : jetArea(2x) : deltaR : fraction : distance next rec jet : pT next jet : deltaPt : jetPt(2x) (hr for area) 
  hnTitle = "variables per jet;centrality;nb. of input tracks; reaction plane bin;#Deltaarea;probe area;rec area;#DeltaR;fraction;distance to closest rec jet;p_{T} of closest rec jet;#delta p_{T};probe p_{T};rec p_{T}";
  Int_t    nbins5[dim5] = {   8 ,   40 ,   3,   81, 100, 100, 50,   52,   51 , 100 ,  101 ,  50 ,  50  };
  Double_t xmin5[dim5]  = {   0.,    0., -.5, -.81,  0.,  0., 0., 0.  , -0.02,   0., -101.,   0.,   0. };
  Double_t xmax5[dim5]  = {  80., 4000., 2.5,  .81,  1.,  1., 1., 1.04,  1.  , 200.,  101., 250., 250. };
  fhnJetsArea = new THnSparseF("fhnJetsArea", hnTitle.Data(), dim5, nbins5, xmin5, xmax5);
  
  
  //before cut
  const Int_t dim6 = 10;
  // cent : nInpTrks : rp bins : fraction : jetPt(2x) : jetEta(2x) : jetPhi(2x) (low resolution) (with fraction, eta, phi, pt cuts possible)
  Int_t    nbins6[dim6] = {  8 ,   40 ,   3, 52 ,  50 ,  50 ,   28,   28,   25,   25 };
  Double_t xmin6[dim6]  = {  0.,    0., -.5,  0.,   0.,   0., -0.7, -0.7,   0.,   0. };
  Double_t xmax6[dim6]  = { 80., 4000., 2.5, 1.04, 250., 250.,  0.7,  0.7, 2*pi, 2*pi };
  hnTitle = "variables before cut;centrality;nb. of input tracks; reaction plane bin;fraction;probe p_{T};rec p_{T};probe #eta; rec #eta;probe #phi;rec #phi";
  fhnJetsBeforeCut1 = new THnSparseF("fhnJetsBeforeCut1", hnTitle.Data(), dim6, nbins6, xmin6, xmax6);
  
  const Int_t dim7 = 10;
  // cent : nInpTrks : rp bins : deltaPt : jetPt(2x) : deltaR : deltaEta : jetEta(2x) (low resolution)
  Int_t    nbins7[dim7] = {   8 ,   40 ,   3,  101 ,  50 ,  50 , 50,    51,   28,   28 };
  Double_t xmin7[dim7]  = {   0.,    0., -.5, -101.,   0.,   0., 0., -1.02, -0.7, -0.7 };
  Double_t xmax7[dim7]  = {  80., 4000., 2.5,  101., 250., 250., 1.,  1.02,  0.7,  0.7 };
  hnTitle = "variables before cut;centrality;nb. of input tracks; reaction plane bin;#delta p_{T};probe p_{T};rec p_{T};#Delta R;#Delta #eta;probe #eta;rec #eta";
  fhnJetsBeforeCut2 = new THnSparseF("fhnJetsBeforeCut2", hnTitle.Data(), dim7, nbins7, xmin7, xmax7);
  
  
  
  fOutputList->Add(fHistEvtSelection);
  fOutputList->Add(fhnEvent);
  fOutputList->Add(fhnJetsRp);
  fOutputList->Add(fhnJetsDeltaPt);
  fOutputList->Add(fhnJetsEta);
  fOutputList->Add(fhnJetsArea);
  fOutputList->Add(fhnJetsBeforeCut1);
  fOutputList->Add(fhnJetsBeforeCut2);

  // =========== Switch on Sumw2 for all histos ===========
  for (Int_t i=0; i<fOutputList->GetEntries(); ++i) {
    TH1 *h1 = dynamic_cast<TH1*>(fOutputList->At(i));
    if (h1){
      h1->Sumw2();
      continue;
    }
	THnSparse *hn = dynamic_cast<THnSparse*>(fOutputList->At(i));
    if (hn){
      hn->Sumw2();
    }	  
  }
  TH1::AddDirectory(oldStatus);

  PostData(1, fOutputList);
}

void AliAnalysisTaskJetResponseV2::UserExec(Option_t *)
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
  Float_t centValue = cent->GetCentralityPercentile("V0M");
  if(fDebug) printf("centrality: %f\n", centValue);
  if (centValue < fCentMin || centValue > fCentMax){
    fHistEvtSelection->Fill(4);
    PostData(1, fOutputList);
    return;
  }


  // multiplicity due to input tracks
  Int_t nInputTracks = GetNInputTracks();
  
  if (nInputTracks < fNInputTracksMin || (fNInputTracksMax > -1 && nInputTracks > fNInputTracksMax)){
       fHistEvtSelection->Fill(5);
       PostData(1, fOutputList);
       return;
  }
  
   
  fHistEvtSelection->Fill(0); // accepted events  
  // -- end event selection --

  Double_t rp = AliAnalysisHelperJetTasks::ReactionPlane(kFALSE);
  Double_t eventEntries[3] = { (Double_t)centValue, (Double_t)nInputTracks, rp };
  fhnEvent->Fill(eventEntries);
  
  
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
                                            aMatchIndex, aPtFraction, fDebug, fMatchMaxDist);
											
  // loop over matched jets
  Int_t ir = -1; // index of matched reconstruced jet
  Float_t fraction = -1.;
  AliAODJet *jet[2]  = { 0x0, 0x0 };
  Float_t jetEta[2]  = { -990., -990. };
  Float_t jetPhi[2]  = { -990., -990. };
  Float_t jetPt[2]   = { -990., -990. };
  Float_t jetArea[2] = { -990., -990. };
  Float_t rpJet[2]   = { -990., -990. };
   
  for(Int_t ig=0; ig<fListJets[0]->GetEntries(); ++ig){
     ir = aMatchIndex[ig];
	 if(ir<0) continue;
	 fraction = aPtFraction[ig];
	 
	 // fetch jets
	 jet[0] = (AliAODJet*)(fListJets[0]->At(ig));
	 jet[1] = (AliAODJet*)(fListJets[1]->At(ir));
	 if(!jet[0] || !jet[1]) continue;
	 
	 // look for distance to next rec jet
	 Float_t distNextJet = -0.01; // no neighbor
	 Float_t ptNextJet = -1.;
	 for(Int_t r=0; r<fListJets[1]->GetEntries(); ++r){
	    if(r==ir) continue;
		Float_t tmpDeltaR = jet[1]->DeltaR((AliAODJet*)fListJets[1]->At(r));
		if(distNextJet<0. || distNextJet>tmpDeltaR){
     		distNextJet = tmpDeltaR;
		    ptNextJet   = ((AliAODJet*)fListJets[1]->At(r))->Pt();
		}
	 }
	 
	 // get parameters
	 for(Int_t i=0; i<fkNbranches; ++i){
	    jetEta[i]  = jet[i]->Eta();
		jetPhi[i]  = jet[i]->Phi();
		jetPt[i]   = jet[i]->Pt();
		jetArea[i] = jet[i]->EffectiveAreaCharged();
		rpJet[i]   = TVector2::Phi_mpi_pi(rp-jetPhi[i]);
	 }
	 Int_t rpBin = AliAnalysisHelperJetTasks::GetPhiBin(TVector2::Phi_mpi_pi(rp-jetPhi[1]), 3);

	 
	 // calculate parameters of associated jets
	 Float_t deltaPt    = jetPt[1]-jetPt[0];
	 Float_t deltaEta   = jetEta[1]-jetEta[0];
	 Float_t deltaPhi   = TVector2::Phi_mpi_pi(jetPhi[1]-jetPhi[0]);
	 Float_t deltaR     = TMath::Sqrt(deltaEta*deltaEta+deltaPhi*deltaPhi);
	 Float_t deltaArea  = jetArea[1]-jetArea[0];
	 
	
	 // fill thnsparse before acceptance cut
	 Double_t jetBeforeCutEntries1[10] = { 
	    (Double_t)centValue, (Double_t)nInputTracks, (Double_t)rpBin,
		(Double_t)fraction,	(Double_t)jetPt[0], (Double_t)jetPt[1], 
		(Double_t)jetEta[0], (Double_t)jetEta[1], (Double_t)jetPhi[0], (Double_t)jetPhi[1] };
	 
	 Double_t jetBeforeCutEntries2[10] = {
        (Double_t)centValue, (Double_t)nInputTracks, (Double_t)rpBin, 
		(Double_t)deltaPt, (Double_t)jetPt[0], (Double_t)jetPt[1],
		(Double_t)deltaR, (Double_t)deltaEta,
		(Double_t)jetEta[0], (Double_t)jetEta[1] };
	 
	 fhnJetsBeforeCut1->Fill(jetBeforeCutEntries1);
	 fhnJetsBeforeCut2->Fill(jetBeforeCutEntries2);
	 
	 
	 // minimum fraction required
	 if(fraction<fJetPtFractionMin) continue;
	 
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
	 
	 // fill thnsparse
	 Double_t jetEntriesRp[6] = {
          (Double_t)centValue, (Double_t)nInputTracks, (Double_t)rpBin, 
          (Double_t)rpJet[0], (Double_t)rpJet[1], (Double_t)jetPt[0]
          };
     fhnJetsRp->Fill(jetEntriesRp);
	 
	 Double_t jetEntriesDeltaPt[6] = {
           (Double_t)centValue, (Double_t)nInputTracks, (Double_t)rpBin,
		   (Double_t)deltaPt, (Double_t)jetPt[0], (Double_t)jetPt[1]
		   };		 
     fhnJetsDeltaPt->Fill(jetEntriesDeltaPt);
	 
	 Double_t jetEntriesEta[10] = {
           (Double_t)centValue, (Double_t)nInputTracks, (Double_t)rpBin,
		   (Double_t)deltaPt, (Double_t)jetPt[0], (Double_t)jetPt[1],
		   (Double_t)deltaR, (Double_t)deltaEta, (Double_t)jetEta[0], (Double_t)jetEta[1]
		   };				 
     fhnJetsEta->Fill(jetEntriesEta);
	 
	 Double_t jetEntriesArea[13] = {
           (Double_t)centValue, (Double_t)nInputTracks, (Double_t)rpBin,
           (Double_t)deltaArea, (Double_t)jetArea[0], (Double_t)jetArea[1],
		   (Double_t)deltaR, (Double_t)fraction, (Double_t)distNextJet, (Double_t)ptNextJet,
		   (Double_t)deltaPt, (Double_t)jetPt[0], (Double_t)jetPt[1]
		   };				 
     fhnJetsArea->Fill(jetEntriesArea);
	 
  }

  PostData(1, fOutputList);
}

void AliAnalysisTaskJetResponseV2::Terminate(const Option_t *)
{
  // Draw result to the screen
  // Called once at the end of the query

  if (!GetOutputData(1))
    return;
}

Int_t AliAnalysisTaskJetResponseV2::GetNInputTracks()
{

  Int_t nInputTracks = 0;
  
  TString jbname(fJetBranchName[1]);
  //needs complete event, use jets without background subtraction
  for(Int_t i=1; i<=3; ++i){
    if(jbname.Contains(Form("B%d",i))) jbname.ReplaceAll(Form("B%d",i),"B0");
  }
  // use only HI event
  if(jbname.Contains("AODextraonly")) jbname.ReplaceAll("AODextraonly","AOD");
  if(jbname.Contains("AODextra")) jbname.ReplaceAll("AODextra","AOD");
  
  if(fDebug) Printf("Multiplicity from jet branch %s", jbname.Data());
  TClonesArray *tmpAODjets = dynamic_cast<TClonesArray*>(fAOD->FindListObject(jbname.Data()));
  if(!tmpAODjets){
    Printf("Jet branch %s not found", jbname.Data());
	Printf("AliAnalysisTaskJetResponseV2::GetNInputTracks FAILED");
	return -1;
  }
      
  for (Int_t iJet=0; iJet<tmpAODjets->GetEntriesFast(); iJet++){
      AliAODJet *jet = dynamic_cast<AliAODJet*>((*tmpAODjets)[iJet]);
	  if(!jet) continue;
	  TRefArray *trackList = jet->GetRefTracks();
	  Int_t nTracks = trackList->GetEntriesFast();
	  nInputTracks += nTracks;
	  if(fDebug) Printf("#jet%d: %d tracks", iJet, nTracks);
  }
  if(fDebug) Printf("---> input tracks: %d", nInputTracks);
  
  return nInputTracks;  
}
