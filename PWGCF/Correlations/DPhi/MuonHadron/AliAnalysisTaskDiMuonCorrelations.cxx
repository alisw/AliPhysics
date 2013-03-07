#include "AliLog.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
#include "TAxis.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "TMath.h"
#include "TString.h"

#include "AliAnalysisTaskDiMuonCorrelations.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliEventPoolManager.h"

ClassImp(AliAnalysisTaskDiMuonCorrelations)

//====================================================================================================================================================

AliAnalysisTaskDiMuonCorrelations::AliAnalysisTaskDiMuonCorrelations() : 
  AliAnalysisTaskSE(), 
  fAOD(0x0),
  fPoolMgr(0x0),
  fMaxChi2Muon(5.), 
  fMinEtaMuon(-4.0), 
  fMaxEtaMuon(-2.5),
  fTriggerMatchLevelMuon(0),
  fNbinsCent(1), 
  fNbinsPt(1),
  fCentAxis(0x0), 
  fPtAxis(0x0),
  fEtaAxis(0x0),
  fHistV0Multiplicity(0x0), 
  fHistITSMultiplicity(0x0),
  fHistCentrality(0x0),
  fHistEvStat(0x0),
  fCentMethod(0),
  fOutputList(0x0)
{

  // Default constructor

  fMuonTrack[0] = NULL;
  fMuonTrack[1] = NULL;

  for (Int_t iCent=0; iCent<fNMaxBinsCentrality; iCent++) {
    for (Int_t iPtBinMuon1=0; iPtBinMuon1<fNMaxBinsPt; iPtBinMuon1++) {
      for (Int_t iPtBinMuon2=0; iPtBinMuon2<fNMaxBinsPt; iPtBinMuon2++) {
	fHistDeltaPhi[iCent][iPtBinMuon1][iPtBinMuon2]       = NULL;
	fHistDeltaPhiMix[iCent][iPtBinMuon1][iPtBinMuon2]    = NULL;
	fHistEtaDeltaPhi[iCent][iPtBinMuon1][iPtBinMuon2]    = NULL;
	fHistEtaDeltaPhiMix[iCent][iPtBinMuon1][iPtBinMuon2] = NULL;
      }
    }
    fHistNMuons_vs_NMuons[iCent]         = NULL;
    fHistNMuons_vs_NMuons_Mixed[iCent]   = NULL;
    fHistTracksEtavsEta[iCent]           = NULL;
    fHistTracksEtavsEta_Mixed[iCent]     = NULL;
    fHistSingleMuonsPt[iCent]            = NULL;
    fHistSingleMuonsPt_Mixed[iCent]      = NULL;
    fHistSingleMuonsEtaPt[iCent]         = NULL;
    fHistSingleMuonsEtaPt_Mixed[iCent]   = NULL;
  }  
  
}


//====================================================================================================================================================

AliAnalysisTaskDiMuonCorrelations::AliAnalysisTaskDiMuonCorrelations(const char *name) : 
  AliAnalysisTaskSE(name), 
  fAOD(0x0),
  fPoolMgr(0x0),
  fMaxChi2Muon(5.), 
  fMinEtaMuon(-4.0), 
  fMaxEtaMuon(-2.5),
  fTriggerMatchLevelMuon(0),
  fNbinsCent(1), 
  fNbinsPt(1),
  fCentAxis(0x0), 
  fPtAxis(0x0),
  fEtaAxis(0x0),
  fHistV0Multiplicity(0x0), 
  fHistITSMultiplicity(0x0),
  fHistCentrality(0x0),
  fHistEvStat(0x0),
  fCentMethod(0),
  fOutputList(0x0)
{

  // Constructor

  fMuonTrack[0] = NULL;
  fMuonTrack[1] = NULL;

  for (Int_t iCent=0; iCent<fNMaxBinsCentrality; iCent++) {
    for (Int_t iPtBinMuon1=0; iPtBinMuon1<fNMaxBinsPt; iPtBinMuon1++) {
      for (Int_t iPtBinMuon2=0; iPtBinMuon2<fNMaxBinsPt; iPtBinMuon2++) {
	fHistDeltaPhi[iCent][iPtBinMuon1][iPtBinMuon2]    = NULL;
	fHistDeltaPhiMix[iCent][iPtBinMuon1][iPtBinMuon2] = NULL;
	fHistEtaDeltaPhi[iCent][iPtBinMuon1][iPtBinMuon2]    = NULL;
	fHistEtaDeltaPhiMix[iCent][iPtBinMuon1][iPtBinMuon2] = NULL;
      }
    }
    fHistNMuons_vs_NMuons[iCent]         = NULL;
    fHistNMuons_vs_NMuons_Mixed[iCent]   = NULL;
    fHistTracksEtavsEta[iCent]           = NULL;
    fHistTracksEtavsEta_Mixed[iCent]     = NULL;
    fHistSingleMuonsPt[iCent]            = NULL;
    fHistSingleMuonsPt_Mixed[iCent]      = NULL;
    fHistSingleMuonsEtaPt[iCent]         = NULL;
    fHistSingleMuonsEtaPt_Mixed[iCent]   = NULL;
  }  
  
  // Define input and output slots here
  DefineOutput(1, TList::Class());
  
}

//====================================================================================================================================================

AliAnalysisTaskDiMuonCorrelations::~AliAnalysisTaskDiMuonCorrelations() {
  
  delete fCentAxis;
  delete fPtAxis;
  delete fEtaAxis;

  if (fOutputList  && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) delete fOutputList;

}

//====================================================================================================================================================

void AliAnalysisTaskDiMuonCorrelations::UserCreateOutputObjects() {
  
  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);

  for (Int_t iCent=0; iCent<fNbinsCent; iCent++) {
    for (Int_t iPtBinMuon1=0; iPtBinMuon1<fNbinsPt; iPtBinMuon1++) {
      for (Int_t iPtBinMuon2=0; iPtBinMuon2<fNbinsPt; iPtBinMuon2++) {
	
	fHistDeltaPhi[iCent][iPtBinMuon1][iPtBinMuon2]    = new TH1D(Form("fHistDeltaPhi_Cent%02d_PtBin%02d_%02d",iCent,iPtBinMuon1,iPtBinMuon2), 
								     Form("%d-%d %%, %3.1f<p_{T}^{#mu1}<%3.1f, %3.1f<p_{T}^{#mu2}<%3.1f",
									  Int_t(fCentAxis->GetBinLowEdge(iCent+1)),
									  Int_t(fCentAxis->GetBinUpEdge(iCent+1)),
									  fPtAxis->GetBinLowEdge(iPtBinMuon1+1),
									  fPtAxis->GetBinUpEdge(iPtBinMuon1+1),
									  fPtAxis->GetBinLowEdge(iPtBinMuon2+1),
									  fPtAxis->GetBinUpEdge(iPtBinMuon2+1)),
								     100, -0.5*TMath::RadToDeg()*TMath::Pi(), 1.5*TMath::RadToDeg()*TMath::Pi());
	
	fHistDeltaPhiMix[iCent][iPtBinMuon1][iPtBinMuon2]    = new TH1D(Form("fHistDeltaPhiMix_Cent%02d_PtBin%02d_%02d",iCent,iPtBinMuon1,iPtBinMuon2), 
									Form("%d-%d %%, %3.1f<p_{T}^{#mu1}<%3.1f, %3.1f<p_{T}^{#mu2}<%3.1f MIXED",
									     Int_t(fCentAxis->GetBinLowEdge(iCent+1)),
									     Int_t(fCentAxis->GetBinUpEdge(iCent+1)),
									     fPtAxis->GetBinLowEdge(iPtBinMuon1+1),
									     fPtAxis->GetBinUpEdge(iPtBinMuon1+1),
									     fPtAxis->GetBinLowEdge(iPtBinMuon2+1),
									     fPtAxis->GetBinUpEdge(iPtBinMuon2+1)),
									100, -0.5*TMath::RadToDeg()*TMath::Pi(), 1.5*TMath::RadToDeg()*TMath::Pi());
	
	fHistDeltaPhi[iCent][iPtBinMuon1][iPtBinMuon2]    -> SetXTitle("#Delta#varphi  [degrees]");
	fHistDeltaPhiMix[iCent][iPtBinMuon1][iPtBinMuon2] -> SetXTitle("#Delta#varphi  [degrees]");
	
	fHistDeltaPhi[iCent][iPtBinMuon1][iPtBinMuon2]    -> Sumw2();
	fHistDeltaPhiMix[iCent][iPtBinMuon1][iPtBinMuon2] -> Sumw2();
	
	fOutputList -> Add(fHistDeltaPhi[iCent][iPtBinMuon1][iPtBinMuon2]);
	fOutputList -> Add(fHistDeltaPhiMix[iCent][iPtBinMuon1][iPtBinMuon2]);
	
	fHistEtaDeltaPhi[iCent][iPtBinMuon1][iPtBinMuon2]    = new TH2D(Form("fHistEtaDeltaPhi_Cent%02d_PtBin%02d_%02d",iCent,iPtBinMuon1,iPtBinMuon2), 
									Form("%d-%d %%, %3.1f<p_{T}^{#mu1}<%3.1f, %3.1f<p_{T}^{#mu2}<%3.1f",
									     Int_t(fCentAxis->GetBinLowEdge(iCent+1)),
									     Int_t(fCentAxis->GetBinUpEdge(iCent+1)),
									     fPtAxis->GetBinLowEdge(iPtBinMuon1+1),
									     fPtAxis->GetBinUpEdge(iPtBinMuon1+1),
									     fPtAxis->GetBinLowEdge(iPtBinMuon2+1),
									     fPtAxis->GetBinUpEdge(iPtBinMuon2+1)),
									100, -0.5*TMath::RadToDeg()*TMath::Pi(), 1.5*TMath::RadToDeg()*TMath::Pi(),
									fEtaAxis->GetNbins(),(Double_t*)fEtaAxis->GetXbins()->GetArray());
	
	fHistEtaDeltaPhiMix[iCent][iPtBinMuon1][iPtBinMuon2]    = new TH2D(Form("fHistEtaDeltaPhiMix_Cent%02d_PtBin%02d_%02d",iCent,iPtBinMuon1,iPtBinMuon2), 
									   Form("%d-%d %%, %3.1f<p_{T}^{#mu1}<%3.1f, %3.1f<p_{T}^{#mu2}<%3.1f MIXED",
										Int_t(fCentAxis->GetBinLowEdge(iCent+1)),
										Int_t(fCentAxis->GetBinUpEdge(iCent+1)),
										fPtAxis->GetBinLowEdge(iPtBinMuon1+1),
										fPtAxis->GetBinUpEdge(iPtBinMuon1+1),
										fPtAxis->GetBinLowEdge(iPtBinMuon2+1),
										fPtAxis->GetBinUpEdge(iPtBinMuon2+1)),
									   100, -0.5*TMath::RadToDeg()*TMath::Pi(), 1.5*TMath::RadToDeg()*TMath::Pi(),
									   fEtaAxis->GetNbins(),(Double_t*)fEtaAxis->GetXbins()->GetArray());
	
	fHistEtaDeltaPhi[iCent][iPtBinMuon1][iPtBinMuon2]    -> SetXTitle("#Delta#varphi  [degrees]");
	fHistEtaDeltaPhiMix[iCent][iPtBinMuon1][iPtBinMuon2] -> SetXTitle("#Delta#varphi  [degrees]");
	fHistEtaDeltaPhi[iCent][iPtBinMuon1][iPtBinMuon2]    -> SetYTitle("#eta");
	fHistEtaDeltaPhiMix[iCent][iPtBinMuon1][iPtBinMuon2] -> SetYTitle("#eta");
	
	fHistEtaDeltaPhi[iCent][iPtBinMuon1][iPtBinMuon2]    -> Sumw2();
	fHistEtaDeltaPhiMix[iCent][iPtBinMuon1][iPtBinMuon2] -> Sumw2();
	
	fOutputList -> Add(fHistEtaDeltaPhi[iCent][iPtBinMuon1][iPtBinMuon2]);
	fOutputList -> Add(fHistEtaDeltaPhiMix[iCent][iPtBinMuon1][iPtBinMuon2]);

      }
    }

    fHistNMuons_vs_NMuons[iCent] = new TH2D(Form("fHistNMuons_vs_NMuons_Cent%02d",iCent),
						  Form("%d-%d %%",Int_t(fCentAxis->GetBinLowEdge(iCent+1)),Int_t(fCentAxis->GetBinUpEdge(iCent+1))),
						  100, 0, 500, 20, 0, 20);
    fHistNMuons_vs_NMuons[iCent] -> SetXTitle("N_{tracks} Muon Arm");
    fHistNMuons_vs_NMuons[iCent] -> SetYTitle("N_{tracks} Muon Arm");
    fHistNMuons_vs_NMuons[iCent] -> Sumw2();

    fHistNMuons_vs_NMuons_Mixed[iCent] = new TH2D(Form("fHistNMuons_vs_NMuons_Mixed_Cent%02d",iCent),
						  Form("%d-%d %% MIXED",Int_t(fCentAxis->GetBinLowEdge(iCent+1)),Int_t(fCentAxis->GetBinUpEdge(iCent+1))),
						  100, 0, 500, 20, 0, 20);
    fHistNMuons_vs_NMuons_Mixed[iCent] -> SetXTitle("N_{tracks} Muon Arm");
    fHistNMuons_vs_NMuons_Mixed[iCent] -> SetYTitle("N_{tracks} Muon Arm");
    fHistNMuons_vs_NMuons_Mixed[iCent] -> Sumw2();

    fHistTracksEtavsEta[iCent]       = new TH2D(Form("fHistTracksEtavsEta_%02d",iCent),       "#eta 1st muon vs #eta 2nd muon", 100,-4.5,-2.,100,-1.5,1.5);
    fHistTracksEtavsEta_Mixed[iCent] = new TH2D(Form("fHistTracksEtavsEta_Mixed_%02d",iCent), "#eta 1st muon vs #eta 2nd muon", 100,-4.5,-2.,100,-1.5,1.5);

    fOutputList -> Add(fHistNMuons_vs_NMuons[iCent]);
    fOutputList -> Add(fHistNMuons_vs_NMuons_Mixed[iCent]);
    fOutputList -> Add(fHistTracksEtavsEta[iCent]);
    fOutputList -> Add(fHistTracksEtavsEta_Mixed[iCent]);

    fHistSingleMuonsPt[iCent]       = new TH1D(Form("fHistSingleMuonPt_Cent%02d",iCent),      "p_{T} for single muons", fNbinsPt, (Double_t*)fPtAxis->GetXbins()->GetArray());
    fHistSingleMuonsPt_Mixed[iCent] = new TH1D(Form("fHistSingleMuonPtmixed_Cent%02d",iCent), "p_{T} for single muons", fNbinsPt, (Double_t*)fPtAxis->GetXbins()->GetArray());
    fOutputList -> Add(fHistSingleMuonsPt[iCent]);
    fOutputList -> Add(fHistSingleMuonsPt_Mixed[iCent]);

    fHistSingleMuonsEtaPt[iCent] = new TH2D(Form("fHistSingleMuonEtaPt_Cent%02d",iCent), "#eta vs p_{T} for single muons",
					    fNbinsPt, (Double_t*)fPtAxis->GetXbins()->GetArray(),
					    fEtaAxis->GetNbins(),(Double_t*)fEtaAxis->GetXbins()->GetArray());
    fHistSingleMuonsEtaPt_Mixed[iCent] = new TH2D(Form("fHistSingleMuonEtaPtmixed_Cent%02d",iCent), "#eta vs p_{T} for single muons",
						  fNbinsPt, (Double_t*)fPtAxis->GetXbins()->GetArray(),
						  fEtaAxis->GetNbins(),(Double_t*)fEtaAxis->GetXbins()->GetArray());
    fOutputList -> Add(fHistSingleMuonsEtaPt[iCent]);
    fOutputList -> Add(fHistSingleMuonsEtaPt_Mixed[iCent]);

  }
  
  fHistV0Multiplicity = new TH1D("fHistV0Multiplicity", "V0 Multiplicity", 500, 0, 1000);
  fHistV0Multiplicity -> SetXTitle("Multiplicity");
  fHistV0Multiplicity -> Sumw2();

  fHistITSMultiplicity = new TH1D("fHistITSMultiplicity", "ITS Multiplicity", 500, 0, 500);
  fHistITSMultiplicity -> SetXTitle("N_{Clusters}");
  fHistITSMultiplicity -> Sumw2();

  fHistCentrality = new TH1D("fHistCentrality", Form("%s Centrality",fCentMethod.Data()), 300, -100, 200);
  fHistCentrality -> SetXTitle("Centrality  [%]");
  fHistCentrality -> Sumw2();

  fOutputList -> Add(fHistV0Multiplicity);
  fOutputList -> Add(fHistITSMultiplicity);
  fOutputList -> Add(fHistCentrality);

  fHistEvStat = new TH1D("fHistEvStat","Event cuts statistics",20,-0.5,19.5);
  fHistEvStat->SetXTitle("Cut index");
  fOutputList->Add(fHistEvStat);

  const Int_t kNZvtxBins  = 10;
  // bins for further buffers are shifted by 100 cm
  Double_t vertexBins[kNZvtxBins+1] = { -10,   -8,  -6,  -4,  -2,   0,   2,   4,   6,   8,  10 };
  Int_t nZvtxBins  = kNZvtxBins;
  Double_t* zvtxbin = vertexBins;

  fPoolMgr = new AliEventPoolManager(1000, 20000, fNbinsCent, (Double_t*)fCentAxis->GetXbins()->GetArray(), nZvtxBins, zvtxbin);

  PostData(1, fOutputList); 

}

//====================================================================================================================================================

void AliAnalysisTaskDiMuonCorrelations::UserExec(Option_t *) {

  fAOD = dynamic_cast<AliAODEvent *>(InputEvent());  
  if (!fAOD) return;  

  Int_t cutIndex = 0;
  fHistEvStat->Fill(cutIndex++);
  // Trigger selection
  if (!(IsTriggerFired())) return;
  fHistEvStat->Fill(cutIndex++);
  
  fHistV0Multiplicity  -> Fill(GetV0Multiplicity());
  fHistITSMultiplicity -> Fill(GetITSMultiplicity());

  Int_t centBin = GetCentBin();
  if (centBin<0) return;
  fHistEvStat->Fill(cutIndex++);
  Double_t percentile = fAOD->GetCentrality()->GetCentralityPercentile(fCentMethod.Data());
  fHistCentrality->Fill(percentile);

  // Vertex selection
  const AliAODVertex* trkVtx = fAOD->GetPrimaryVertex();
  if (!trkVtx || trkVtx->GetNContributors()<=0) return;
  TString vtxTtl = trkVtx->GetTitle();
  if (!vtxTtl.Contains("VertexerTracks")) return;
  fHistEvStat->Fill(cutIndex++);
  Double_t zvtx = trkVtx->GetZ();
  const AliAODVertex* spdVtx = fAOD->GetPrimaryVertexSPD();
  if (spdVtx->GetNContributors()<=0) return;
  TString vtxTyp = spdVtx->GetTitle();
  Double_t cov[6]={0};
  spdVtx->GetCovarianceMatrix(cov);
  Double_t zRes = TMath::Sqrt(cov[5]);
  if (vtxTyp.Contains("vertexer:Z") && (zRes>0.25)) return;
  if (TMath::Abs(spdVtx->GetZ() - trkVtx->GetZ())>0.5) return;
  fHistEvStat->Fill(cutIndex++);

  if (TMath::Abs(zvtx) > 10.) return;
  fHistEvStat->Fill(cutIndex++);

  TObjArray *tracksMuonArm = GetAcceptedTracksMuonArm(fAOD);
  if (tracksMuonArm->GetEntriesFast() < 2) {
    delete tracksMuonArm;
    return;
  }
  fHistEvStat->Fill(cutIndex++);
  
  fHistNMuons_vs_NMuons[centBin] -> Fill(tracksMuonArm->GetEntries(), tracksMuonArm->GetEntries());

  AliDebug(1, Form("Single Event analysis : nTracksMuonArm = %4d", tracksMuonArm->GetEntries()));

  // Same event  
  for (Int_t iTrMuon1=0; iTrMuon1<tracksMuonArm->GetEntriesFast(); iTrMuon1++) {
    fMuonTrack[0] = (AliAODTrack*) tracksMuonArm->At(iTrMuon1);
    fHistSingleMuonsPt[centBin]->Fill(fMuonTrack[0]->Pt());
    fHistSingleMuonsEtaPt[centBin]->Fill(fMuonTrack[0]->Pt(),fMuonTrack[0]->Eta());
    for (Int_t iTrMuon2=iTrMuon1+1; iTrMuon2<tracksMuonArm->GetEntriesFast(); iTrMuon2++) {
      fMuonTrack[1] = (AliAODTrack*) tracksMuonArm -> At(iTrMuon2);
      FillHistograms(centBin, kSingleEvent);
    }
  }

  // Mixed event
  AliEventPool* pool = fPoolMgr->GetEventPool(percentile, zvtx);
  //pool->PrintInfo();
  if (pool->IsReady() || pool->NTracksInPool() > 2000 || pool->GetCurrentNEvents() >= 5) {
    for (Int_t jMix=0; jMix<pool->GetCurrentNEvents(); jMix++) {
      TObjArray *mixedTracks = pool->GetEvent(jMix);     // Cvetan, here we would need to retrieve the MUON tracks of the event we mix...
      fHistNMuons_vs_NMuons_Mixed[centBin]->Fill(mixedTracks->GetEntriesFast(), tracksMuonArm->GetEntriesFast());
      for (Int_t iTrMuon1=0; iTrMuon1<tracksMuonArm->GetEntriesFast(); iTrMuon1++) {
	fMuonTrack[0] = (AliAODTrack*) tracksMuonArm->At(iTrMuon1);
	fHistSingleMuonsPt_Mixed[centBin]->Fill(fMuonTrack[0]->Pt());
	fHistSingleMuonsEtaPt_Mixed[centBin]->Fill(fMuonTrack[0]->Pt(),fMuonTrack[0]->Eta());
	for (Int_t iTrMuon2=0; iTrMuon2<mixedTracks->GetEntriesFast(); iTrMuon2++) {
	  fMuonTrack[1] = (AliAODTrack*) mixedTracks -> At(iTrMuon2);
	  FillHistograms(centBin, kMixedEvent);
	}
      }
    }
  }
  //  pool->UpdatePool(tracksCentralBarrel);    // Cvetan, I think I do not fully understand what this line does, and how it should be modified in a Di-Muon analysis

  delete tracksMuonArm;

  PostData(1, fOutputList); 

}

//====================================================================================================================================================

void AliAnalysisTaskDiMuonCorrelations::FillHistograms(Int_t centrality, Int_t option) {

  Int_t ptBinTrackMuon1 = fPtAxis -> FindBin(fMuonTrack[0]->Pt());
  Int_t ptBinTrackMuon2 = fPtAxis -> FindBin(fMuonTrack[1]->Pt());

  if (ptBinTrackMuon1<1 || ptBinTrackMuon1>fNbinsPt || ptBinTrackMuon2<1 || ptBinTrackMuon2>fNbinsPt) return;

  Double_t deltaPhi = fMuonTrack[0]->Phi() - fMuonTrack[1]->Phi();
  if (deltaPhi >  1.5*TMath::Pi()) deltaPhi -= TMath::TwoPi();
  if (deltaPhi < -0.5*TMath::Pi()) deltaPhi += TMath::TwoPi();

  if (option==kSingleEvent)     fHistDeltaPhi[centrality][ptBinTrackMuon1-1][ptBinTrackMuon2-1]    -> Fill(TMath::RadToDeg()*deltaPhi);
  else if (option==kMixedEvent) fHistDeltaPhiMix[centrality][ptBinTrackMuon1-1][ptBinTrackMuon2-1] -> Fill(TMath::RadToDeg()*deltaPhi);

  if (option==kSingleEvent)     fHistEtaDeltaPhi[centrality][ptBinTrackMuon1-1][ptBinTrackMuon2-1]    -> Fill(TMath::RadToDeg()*deltaPhi,fMuonTrack[0]->Eta());
  else if (option==kMixedEvent) fHistEtaDeltaPhiMix[centrality][ptBinTrackMuon1-1][ptBinTrackMuon2-1] -> Fill(TMath::RadToDeg()*deltaPhi,fMuonTrack[0]->Eta());

  if (option==kSingleEvent)     fHistTracksEtavsEta[centrality]->Fill(fMuonTrack[0]->Eta(),fMuonTrack[1]->Eta());
  else if (option==kMixedEvent) fHistTracksEtavsEta_Mixed[centrality]->Fill(fMuonTrack[0]->Eta(),fMuonTrack[1]->Eta());

}

//====================================================================================================================================================

Bool_t AliAnalysisTaskDiMuonCorrelations::IsTriggerFired() {
  
  Bool_t isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kINT7); 

  return isSelected;
}

//====================================================================================================================================================

Float_t AliAnalysisTaskDiMuonCorrelations::GetV0Multiplicity() {
  
  Float_t multiplicity=0;
  for (Int_t iChannel=0; iChannel<64; iChannel++) multiplicity += fAOD->GetVZEROData()->GetMultiplicity(iChannel);
  return multiplicity;

}

//====================================================================================================================================================

TObjArray* AliAnalysisTaskDiMuonCorrelations::GetAcceptedTracksMuonArm(AliAODEvent *aodEvent) {

  // fills the array of muon tracks that pass the cuts

  TObjArray *tracks = new TObjArray;
  tracks->SetOwner(kFALSE);

  Int_t nTracks = aodEvent->GetNTracks();

  AliAODTrack *track = 0;
  
  for (Int_t iTrack=0; iTrack<nTracks; iTrack++) {
    track = aodEvent->GetTrack(iTrack);
    if (track->IsMuonTrack() && track->GetMatchTrigger()>=fTriggerMatchLevelMuon) {
      if (track->Chi2perNDF() < fMaxChi2Muon) {
	if (track->Eta() > fMinEtaMuon && track->Eta() < fMaxEtaMuon) {
	  tracks->Add(new AliAODTrack(*track));
	}
      }
    }
  }

  return tracks;

}

//====================================================================================================================================================

void AliAnalysisTaskDiMuonCorrelations::SetCentBinning(Int_t nBins, Double_t *limits) {

  if (nBins>fNMaxBinsCentrality) {
    AliInfo(Form("WARNING : only %d centrality bins (out of the %d proposed) will be considered",fNMaxBinsCentrality,nBins));
    nBins = fNMaxBinsCentrality;
  }
  if (nBins<=0) {
    AliInfo("WARNING : at least one centrality bin must be considered");
    nBins = 1;
  }
  
  fNbinsCent = nBins;
  fCentAxis  = new TAxis(fNbinsCent, limits);

}

//====================================================================================================================================================

void AliAnalysisTaskDiMuonCorrelations::SetPtBinning(Int_t nBins, Double_t *limits) {

  if (nBins>fNMaxBinsPt) {
    AliInfo(Form("WARNING : only %d pt bins (out of the %d proposed) will be considered",fNMaxBinsPt,nBins));
    nBins = fNMaxBinsPt;
  }
  if (nBins<=0) {
    AliInfo("WARNING : at least one pt bin must be considered");
    nBins = 1;
  }
  
  fNbinsPt = nBins;
  fPtAxis  = new TAxis(fNbinsPt, limits);

}

//====================================================================================================================================================

void AliAnalysisTaskDiMuonCorrelations::SetEtaBinning(Int_t nBins, Double_t *limits) {

  if (nBins>fNMaxBinsEta) {
    AliInfo(Form("WARNING : only %d pt bins (out of the %d proposed) will be considered",fNMaxBinsEta,nBins));
    nBins = fNMaxBinsEta;
  }
  if (nBins<=0) {
    AliInfo("WARNING : at least one pt bin must be considered");
    nBins = 1;
  }
  
  fEtaAxis  = new TAxis(nBins, limits);

}

//====================================================================================================================================================

Int_t AliAnalysisTaskDiMuonCorrelations::GetCentBin() {

  Double_t percentile = fAOD->GetCentrality()->GetCentralityPercentile(fCentMethod.Data());

  Int_t bin = fCentAxis->FindBin(percentile) - 1;
  if (bin >= fNbinsCent) bin = -1;
  return bin;
  
}

//====================================================================================================================================================

Double_t AliAnalysisTaskDiMuonCorrelations::GetITSMultiplicity() {

  Double_t multiplicity = fAOD->GetHeader()->GetNumberOfITSClusters(1);

  return multiplicity;

}

//====================================================================================================================================================

void AliAnalysisTaskDiMuonCorrelations::Terminate(Option_t *) {

  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputList) {
    printf("ERROR: Output list not available\n");
    return;
  }  

}

//====================================================================================================================================================

