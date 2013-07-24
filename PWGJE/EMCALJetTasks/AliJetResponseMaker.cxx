// $Id$
//
// Emcal jet response matrix maker task.
//
// Author: S. Aiola

#include "AliJetResponseMaker.h"

#include <TClonesArray.h>
#include <TH2F.h>
#include <THnSparse.h>
#include <TLorentzVector.h>

#include "AliAnalysisManager.h"
#include "AliVCluster.h"
#include "AliVTrack.h"
#include "AliEmcalJet.h"
#include "AliLog.h"
#include "AliRhoParameter.h"
#include "AliNamedArrayI.h"

ClassImp(AliJetResponseMaker)

//________________________________________________________________________
AliJetResponseMaker::AliJetResponseMaker() : 
  AliAnalysisTaskEmcalJet("AliJetResponseMaker", kTRUE),
  fTracks2Name(""),
  fCalo2Name(""),
  fJets2Name(""),
  fRho2Name(""),
  fJet2Radius(0.4),
  fJet2AreaCut(-1),
  fPtBiasJet2Track(0),
  fPtBiasJet2Clus(0),
  fJet2MinEta(-999),
  fJet2MaxEta(-999),
  fJet2MinPhi(-999),
  fJet2MaxPhi(-999),
  fMaxClusterPt2(1000),
  fMaxTrackPt2(1000),
  fAreCollections1MC(kFALSE),  
  fAreCollections2MC(kTRUE),
  fMatching(kNoMatching),
  fMatchingPar1(0),
  fMatchingPar2(0),
  fUseCellsToMatch(kFALSE),
  fMinJetMCPt(1),
  fHistoType(0),
  fDeltaPtAxis(0),
  fDeltaEtaDeltaPhiAxis(0),
  fNEFAxis(0),
  fZAxis(0),
  fDoJet2Histogram(0),
  fTracks2(0),
  fCaloClusters2(0),
  fJets2(0),
  fRho2(0),
  fRho2Val(0),
  fTracks2Map(0),
  fHistLeadingJets1PtArea(0),
  fHistLeadingJets1CorrPtArea(0),
  fHistLeadingJets2PtArea(0),
  fHistLeadingJets2CorrPtArea(0),
  fHistLeadingJets2PtAreaAcceptance(0),
  fHistLeadingJets2CorrPtAreaAcceptance(0),
  fHistJets1(0),
  fHistJets2(0),
  fHistJets2Acceptance(0),
  fHistMatching(0),
  fHistJets1PhiEta(0),
  fHistJets1PtArea(0),
  fHistJets1CorrPtArea(0),
  fHistJets1NEFvsPt(0),
  fHistJets1CEFvsCEFPt(0),
  fHistJets1ZvsPt(0),
  fHistJets2PhiEta(0),
  fHistJets2PtArea(0),
  fHistJets2CorrPtArea(0),
  fHistJets2PhiEtaAcceptance(0),
  fHistJets2PtAreaAcceptance(0),
  fHistJets2CorrPtAreaAcceptance(0),
  fHistJets2NEFvsPt(0),
  fHistJets2CEFvsCEFPt(0),
  fHistJets2ZvsPt(0),
  fHistCommonEnergy1vsJet1Pt(0),
  fHistCommonEnergy2vsJet2Pt(0),
  fHistDistancevsJet1Pt(0),
  fHistDistancevsJet2Pt(0),
  fHistDistancevsCommonEnergy1(0),
  fHistDistancevsCommonEnergy2(0),
  fHistCommonEnergy1vsCommonEnergy2(0),
  fHistDeltaEtaDeltaPhi(0),
  fHistJet2PtOverJet1PtvsJet2Pt(0),
  fHistJet1PtOverJet2PtvsJet1Pt(0),
  fHistDeltaPtvsJet1Pt(0),
  fHistDeltaPtvsJet2Pt(0),
  fHistDeltaPtOverJet1PtvsJet1Pt(0),
  fHistDeltaPtOverJet2PtvsJet2Pt(0),
  fHistDeltaPtvsDistance(0),
  fHistDeltaPtvsCommonEnergy1(0),
  fHistDeltaPtvsCommonEnergy2(0),
  fHistDeltaPtvsArea1(0),
  fHistDeltaPtvsArea2(0),
  fHistDeltaPtvsDeltaArea(0),
  fHistJet1PtvsJet2Pt(0),
  fHistDeltaCorrPtOverJet1CorrPtvsJet1CorrPt(0),
  fHistDeltaCorrPtOverJet2CorrPtvsJet2CorrPt(0),
  fHistDeltaCorrPtvsJet1CorrPt(0),
  fHistDeltaCorrPtvsJet2CorrPt(0),
  fHistDeltaCorrPtvsDistance(0),
  fHistDeltaCorrPtvsCommonEnergy1(0),
  fHistDeltaCorrPtvsCommonEnergy2(0),
  fHistDeltaCorrPtvsArea1(0),
  fHistDeltaCorrPtvsArea2(0),
  fHistDeltaCorrPtvsDeltaArea(0),
  fHistJet1CorrPtvsJet2CorrPt(0),
  fHistDeltaMCPtOverJet1MCPtvsJet1MCPt(0),
  fHistDeltaMCPtOverJet2PtvsJet2Pt(0),
  fHistDeltaMCPtvsJet1MCPt(0),
  fHistDeltaMCPtvsJet2Pt(0),
  fHistDeltaMCPtvsDistance(0),
  fHistDeltaMCPtvsCommonEnergy1(0),
  fHistDeltaMCPtvsCommonEnergy2(0),
  fHistDeltaMCPtvsArea1(0),
  fHistDeltaMCPtvsArea2(0),
  fHistDeltaMCPtvsDeltaArea(0),
  fHistJet1MCPtvsJet2Pt(0)
{
  // Default constructor.

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliJetResponseMaker::AliJetResponseMaker(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fTracks2Name("MCParticles"),
  fCalo2Name(""),
  fJets2Name("MCJets"),
  fRho2Name(""),
  fJet2Radius(0.4),
  fJet2AreaCut(-1),
  fPtBiasJet2Track(0),
  fPtBiasJet2Clus(0),
  fJet2MinEta(-999),
  fJet2MaxEta(-999),
  fJet2MinPhi(-999),
  fJet2MaxPhi(-999),
  fMaxClusterPt2(1000),
  fMaxTrackPt2(1000),
  fAreCollections1MC(kFALSE),  
  fAreCollections2MC(kTRUE),
  fMatching(kNoMatching),
  fMatchingPar1(0),
  fMatchingPar2(0),
  fUseCellsToMatch(kFALSE),
  fMinJetMCPt(1),
  fHistoType(0),
  fDeltaPtAxis(0),
  fDeltaEtaDeltaPhiAxis(0),
  fNEFAxis(0),
  fZAxis(0),
  fDoJet2Histogram(0),
  fTracks2(0),
  fCaloClusters2(0),
  fJets2(0),
  fRho2(0),
  fRho2Val(0),
  fTracks2Map(0),
  fHistLeadingJets1PtArea(0),
  fHistLeadingJets1CorrPtArea(0),
  fHistLeadingJets2PtArea(0),
  fHistLeadingJets2CorrPtArea(0),
  fHistLeadingJets2PtAreaAcceptance(0),
  fHistLeadingJets2CorrPtAreaAcceptance(0),
  fHistJets1(0),
  fHistJets2(0),
  fHistJets2Acceptance(0),
  fHistMatching(0),
  fHistJets1PhiEta(0),
  fHistJets1PtArea(0),
  fHistJets1CorrPtArea(0),
  fHistJets1NEFvsPt(0),
  fHistJets1CEFvsCEFPt(0),
  fHistJets1ZvsPt(0),
  fHistJets2PhiEta(0),
  fHistJets2PtArea(0),
  fHistJets2CorrPtArea(0),
  fHistJets2PhiEtaAcceptance(0),
  fHistJets2PtAreaAcceptance(0),
  fHistJets2CorrPtAreaAcceptance(0),
  fHistJets2NEFvsPt(0),
  fHistJets2CEFvsCEFPt(0),
  fHistJets2ZvsPt(0),
  fHistCommonEnergy1vsJet1Pt(0),
  fHistCommonEnergy2vsJet2Pt(0),
  fHistDistancevsJet1Pt(0),
  fHistDistancevsJet2Pt(0),
  fHistDistancevsCommonEnergy1(0),
  fHistDistancevsCommonEnergy2(0),
  fHistCommonEnergy1vsCommonEnergy2(0),
  fHistDeltaEtaDeltaPhi(0),
  fHistJet2PtOverJet1PtvsJet2Pt(0),
  fHistJet1PtOverJet2PtvsJet1Pt(0),
  fHistDeltaPtvsJet1Pt(0),
  fHistDeltaPtvsJet2Pt(0),
  fHistDeltaPtOverJet1PtvsJet1Pt(0),
  fHistDeltaPtOverJet2PtvsJet2Pt(0),
  fHistDeltaPtvsDistance(0),
  fHistDeltaPtvsCommonEnergy1(0),
  fHistDeltaPtvsCommonEnergy2(0),
  fHistDeltaPtvsArea1(0),
  fHistDeltaPtvsArea2(0),
  fHistDeltaPtvsDeltaArea(0),
  fHistJet1PtvsJet2Pt(0),
  fHistDeltaCorrPtOverJet1CorrPtvsJet1CorrPt(0),
  fHistDeltaCorrPtOverJet2CorrPtvsJet2CorrPt(0),
  fHistDeltaCorrPtvsJet1CorrPt(0),
  fHistDeltaCorrPtvsJet2CorrPt(0),
  fHistDeltaCorrPtvsDistance(0),
  fHistDeltaCorrPtvsCommonEnergy1(0),
  fHistDeltaCorrPtvsCommonEnergy2(0),
  fHistDeltaCorrPtvsArea1(0),
  fHistDeltaCorrPtvsArea2(0),
  fHistDeltaCorrPtvsDeltaArea(0),
  fHistJet1CorrPtvsJet2CorrPt(0),
  fHistDeltaMCPtOverJet1MCPtvsJet1MCPt(0),
  fHistDeltaMCPtOverJet2PtvsJet2Pt(0),
  fHistDeltaMCPtvsJet1MCPt(0),
  fHistDeltaMCPtvsJet2Pt(0),
  fHistDeltaMCPtvsDistance(0),
  fHistDeltaMCPtvsCommonEnergy1(0),
  fHistDeltaMCPtvsCommonEnergy2(0),
  fHistDeltaMCPtvsArea1(0),
  fHistDeltaMCPtvsArea2(0),
  fHistDeltaMCPtvsDeltaArea(0),
  fHistJet1MCPtvsJet2Pt(0)
{
  // Standard constructor.

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliJetResponseMaker::~AliJetResponseMaker()
{
  // Destructor
}


//________________________________________________________________________
void AliJetResponseMaker::AllocateTH2()
{
  // Allocate TH2 histograms.

  fHistJets1PhiEta = new TH2F("fHistJets1PhiEta", "fHistJets1PhiEta", 40, -1, 1, 40, 0, TMath::Pi()*2);
  fHistJets1PhiEta->GetXaxis()->SetTitle("#eta");
  fHistJets1PhiEta->GetYaxis()->SetTitle("#phi");
  fOutput->Add(fHistJets1PhiEta);
  
  fHistJets1PtArea = new TH2F("fHistJets1PtArea", "fHistJets1PtArea", fNbins/2, 0, 1.5, fNbins, fMinBinPt, fMaxBinPt);
  fHistJets1PtArea->GetXaxis()->SetTitle("area");
  fHistJets1PtArea->GetYaxis()->SetTitle("p_{T,1} (GeV/c)");
  fHistJets1PtArea->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistJets1PtArea);
 
  if (!fRhoName.IsNull()) {
    fHistJets1CorrPtArea = new TH2F("fHistJets1CorrPtArea", "fHistJets1CorrPtArea", 
				    fNbins/2, 0, 1.5, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    fHistJets1CorrPtArea->GetXaxis()->SetTitle("area");
    fHistJets1CorrPtArea->GetYaxis()->SetTitle("p_{T,1}^{corr} (GeV/c)");
    fHistJets1CorrPtArea->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistJets1CorrPtArea);
  }

  fHistJets1ZvsPt = new TH2F("fHistJets1ZvsPt", "fHistJets1ZvsPt", 120, 0, 1.2, fNbins, fMinBinPt, fMaxBinPt);
  fHistJets1ZvsPt->GetXaxis()->SetTitle("Z");
  fHistJets1ZvsPt->GetYaxis()->SetTitle("p_{T,1} (GeV/c)");
  fHistJets1ZvsPt->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistJets1ZvsPt);
  
  fHistJets1NEFvsPt = new TH2F("fHistJets1NEFvsPt", "fHistJets1NEFvsPt", 120, 0, 1.2, fNbins, fMinBinPt, fMaxBinPt);
  fHistJets1NEFvsPt->GetXaxis()->SetTitle("NEF");
  fHistJets1NEFvsPt->GetYaxis()->SetTitle("p_{T,1} (GeV/c)");
  fHistJets1NEFvsPt->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistJets1NEFvsPt);
  
  fHistJets1CEFvsCEFPt = new TH2F("fHistJets1CEFvsCEFPt", "fHistJets1CEFvsCEFPt", 120, 0, 1.2, fNbins, fMinBinPt, fMaxBinPt);
  fHistJets1CEFvsCEFPt->GetXaxis()->SetTitle("1-NEF");
  fHistJets1CEFvsCEFPt->GetYaxis()->SetTitle("(1-NEF)*p_{T,1} (GeV/c)");
  fHistJets1CEFvsCEFPt->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistJets1CEFvsCEFPt);

  // Jets 2 histograms

  fHistJets2PhiEta = new TH2F("fHistJets2PhiEta", "fHistJets2PhiEta", 40, -1, 1, 40, 0, TMath::Pi()*2);
  fHistJets2PhiEta->GetXaxis()->SetTitle("#eta");
  fHistJets2PhiEta->GetYaxis()->SetTitle("#phi");
  fOutput->Add(fHistJets2PhiEta);

  fHistJets2PtArea = new TH2F("fHistJets2PtArea", "fHistJets2PtArea", fNbins/2, 0, 1.5, fNbins, fMinBinPt, fMaxBinPt);
  fHistJets2PtArea->GetXaxis()->SetTitle("area");
  fHistJets2PtArea->GetYaxis()->SetTitle("p_{T,2} (GeV/c)");
  fHistJets2PtArea->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistJets2PtArea);

  if (!fRho2Name.IsNull()) {
    fHistJets2CorrPtArea = new TH2F("fHistJets2CorrPtArea", "fHistJets2CorrPtArea", 
				    fNbins/2, 0, 1.5, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    fHistJets2CorrPtArea->GetXaxis()->SetTitle("area");
    fHistJets2CorrPtArea->GetYaxis()->SetTitle("p_{T,2}^{corr} (GeV/c)");
    fHistJets2CorrPtArea->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistJets2CorrPtArea);
  }

  fHistJets2PhiEtaAcceptance = new TH2F("fHistJets2PhiEtaAcceptance", "fHistJets2PhiEtaAcceptance", 40, -1, 1, 40, 0, TMath::Pi()*2);
  fHistJets2PhiEtaAcceptance->GetXaxis()->SetTitle("#eta");
  fHistJets2PhiEtaAcceptance->GetYaxis()->SetTitle("#phi");
  fOutput->Add(fHistJets2PhiEtaAcceptance);

  fHistJets2PtAreaAcceptance = new TH2F("fHistJets2PtAreaAcceptance", "fHistJets2PtAreaAcceptance", 
					fNbins/2, 0, 1.5, fNbins, fMinBinPt, fMaxBinPt);
  fHistJets2PtAreaAcceptance->GetXaxis()->SetTitle("area");
  fHistJets2PtAreaAcceptance->GetYaxis()->SetTitle("p_{T,2} (GeV/c)");
  fHistJets2PtAreaAcceptance->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistJets2PtAreaAcceptance);

  if (!fRho2Name.IsNull()) {
    fHistJets2CorrPtAreaAcceptance = new TH2F("fHistJets2CorrPtAreaAcceptance", "fHistJets2CorrPtAreaAcceptance", 
					      fNbins/2, 0, 1.5, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    fHistJets2CorrPtAreaAcceptance->GetXaxis()->SetTitle("area");
    fHistJets2CorrPtAreaAcceptance->GetYaxis()->SetTitle("p_{T,2}^{corr} (GeV/c)");
    fHistJets2CorrPtAreaAcceptance->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistJets2CorrPtAreaAcceptance);
  }

  fHistJets2ZvsPt = new TH2F("fHistJets2ZvsPt", "fHistJets2ZvsPt", 120, 0, 1.2, fNbins, fMinBinPt, fMaxBinPt);
  fHistJets2ZvsPt->GetXaxis()->SetTitle("Z");
  fHistJets2ZvsPt->GetYaxis()->SetTitle("p_{T,2} (GeV/c)");
  fHistJets2ZvsPt->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistJets2ZvsPt);
  
  fHistJets2NEFvsPt = new TH2F("fHistJets2NEFvsPt", "fHistJets2NEFvsPt", 120, 0, 1.2, fNbins, fMinBinPt, fMaxBinPt);
  fHistJets2NEFvsPt->GetXaxis()->SetTitle("NEF");
  fHistJets2NEFvsPt->GetYaxis()->SetTitle("p_{T,2} (GeV/c)");
  fHistJets2NEFvsPt->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistJets2NEFvsPt);
  
  fHistJets2CEFvsCEFPt = new TH2F("fHistJets2CEFvsCEFPt", "fHistJets2CEFvsCEFPt", 120, 0, 1.2, fNbins, fMinBinPt, fMaxBinPt);
  fHistJets2CEFvsCEFPt->GetXaxis()->SetTitle("1-NEF");
  fHistJets2CEFvsCEFPt->GetYaxis()->SetTitle("(1-NEF)*p_{T,2} (GeV/c)");
  fHistJets2CEFvsCEFPt->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistJets2CEFvsCEFPt);

  // Matching histograms

  fHistCommonEnergy1vsJet1Pt = new TH2F("fHistCommonEnergy1vsJet1Pt", "fHistCommonEnergy1vsJet1Pt", fNbins/2, 0, 1.2, fNbins/2, fMinBinPt, fMaxBinPt);
  fHistCommonEnergy1vsJet1Pt->GetXaxis()->SetTitle("Common energy 1 (%)");
  fHistCommonEnergy1vsJet1Pt->GetYaxis()->SetTitle("p_{T,1}");  
  fHistCommonEnergy1vsJet1Pt->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistCommonEnergy1vsJet1Pt);

  fHistCommonEnergy2vsJet2Pt = new TH2F("fHistCommonEnergy2vsJet2Pt", "fHistCommonEnergy2vsJet2Pt", fNbins/2, 0, 1.2, fNbins/2, fMinBinPt, fMaxBinPt);
  fHistCommonEnergy2vsJet2Pt->GetXaxis()->SetTitle("Common energy 2 (%)");
  fHistCommonEnergy2vsJet2Pt->GetYaxis()->SetTitle("p_{T,2}");  
  fHistCommonEnergy2vsJet2Pt->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistCommonEnergy2vsJet2Pt);

  fHistDistancevsJet1Pt = new TH2F("fHistDistancevsJet1Pt", "fHistDistancevsJet1Pt", fNbins/2, 0, 1.2, fNbins/2, fMinBinPt, fMaxBinPt);
  fHistDistancevsJet1Pt->GetXaxis()->SetTitle("Distance");
  fHistDistancevsJet1Pt->GetYaxis()->SetTitle("p_{T,1}");  
  fHistDistancevsJet1Pt->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistDistancevsJet1Pt);

  fHistDistancevsJet2Pt = new TH2F("fHistDistancevsJet2Pt", "fHistDistancevsJet2Pt", fNbins/2, 0, 1.2, fNbins/2, fMinBinPt, fMaxBinPt);
  fHistDistancevsJet2Pt->GetXaxis()->SetTitle("Distance");
  fHistDistancevsJet2Pt->GetYaxis()->SetTitle("p_{T,2}");  
  fHistDistancevsJet2Pt->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistDistancevsJet2Pt);

  fHistDistancevsCommonEnergy1 = new TH2F("fHistDistancevsCommonEnergy1", "fHistDistancevsCommonEnergy1", fNbins/2, 0, 1.2, fNbins/2, 0, 1.2);
  fHistDistancevsCommonEnergy1->GetXaxis()->SetTitle("Distance");
  fHistDistancevsCommonEnergy1->GetYaxis()->SetTitle("Common energy 1 (%)");  
  fHistDistancevsCommonEnergy1->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistDistancevsCommonEnergy1);

  fHistDistancevsCommonEnergy2 = new TH2F("fHistDistancevsCommonEnergy2", "fHistDistancevsCommonEnergy2", fNbins/2, 0, 1.2, fNbins/2, 0, 1.2);
  fHistDistancevsCommonEnergy2->GetXaxis()->SetTitle("Distance");
  fHistDistancevsCommonEnergy2->GetYaxis()->SetTitle("Common energy 2 (%)");  
  fHistDistancevsCommonEnergy2->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistDistancevsCommonEnergy2);

  fHistCommonEnergy1vsCommonEnergy2 = new TH2F("fHistCommonEnergy1vsCommonEnergy2", "fHistCommonEnergy1vsCommonEnergy2", fNbins/2, 0, 1.2, fNbins/2, 0, 1.2);
  fHistCommonEnergy1vsCommonEnergy2->GetXaxis()->SetTitle("Common energy 1 (%)");
  fHistCommonEnergy1vsCommonEnergy2->GetYaxis()->SetTitle("Common energy 2 (%)");  
  fHistCommonEnergy1vsCommonEnergy2->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistCommonEnergy1vsCommonEnergy2);

  fHistDeltaEtaDeltaPhi = new TH2F("fHistDeltaEtaDeltaPhi", "fHistDeltaEtaDeltaPhi", fNbins/4, -1, 1, fNbins/4, -TMath::Pi()/2, TMath::Pi()*3/2);
  fHistDeltaEtaDeltaPhi->GetXaxis()->SetTitle("Common energy 1 (%)");
  fHistDeltaEtaDeltaPhi->GetYaxis()->SetTitle("Common energy 2 (%)");  
  fHistDeltaEtaDeltaPhi->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistDeltaEtaDeltaPhi);

  fHistJet2PtOverJet1PtvsJet2Pt = new TH2F("fHistJet2PtOverJet1PtvsJet2Pt", "fHistJet2PtOverJet1PtvsJet2Pt", fNbins, fMinBinPt, fMaxBinPt, 300, 0, 1.5);
  fHistJet2PtOverJet1PtvsJet2Pt->GetXaxis()->SetTitle("p_{T,2}");  
  fHistJet2PtOverJet1PtvsJet2Pt->GetYaxis()->SetTitle("p_{T,2} / p_{T,1}");
  fHistJet2PtOverJet1PtvsJet2Pt->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistJet2PtOverJet1PtvsJet2Pt);

  fHistJet1PtOverJet2PtvsJet1Pt = new TH2F("fHistJet1PtOverJet2PtvsJet1Pt", "fHistJet1PtOverJet2PtvsJet1Pt", fNbins, fMinBinPt, fMaxBinPt, 300, 0, 1.5);
  fHistJet1PtOverJet2PtvsJet1Pt->GetXaxis()->SetTitle("p_{T,1}");  
  fHistJet1PtOverJet2PtvsJet1Pt->GetYaxis()->SetTitle("p_{T,1} / p_{T,2}");
  fHistJet1PtOverJet2PtvsJet1Pt->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistJet1PtOverJet2PtvsJet1Pt);

  fHistDeltaPtvsJet1Pt = new TH2F("fHistDeltaPtvsJet1Pt", "fHistDeltaPtvsJet1Pt", 
				  fNbins, fMinBinPt, fMaxBinPt, 2*fNbins, -fMaxBinPt, fMaxBinPt);
  fHistDeltaPtvsJet1Pt->GetXaxis()->SetTitle("p_{T,1}");  
  fHistDeltaPtvsJet1Pt->GetYaxis()->SetTitle("#deltap_{T} (GeV/c)");
  fHistDeltaPtvsJet1Pt->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistDeltaPtvsJet1Pt);

  fHistDeltaPtvsJet2Pt = new TH2F("fHistDeltaPtvsJet2Pt", "fHistDeltaPtvsJet2Pt", 
				  fNbins, fMinBinPt, fMaxBinPt, 2*fNbins, -fMaxBinPt, fMaxBinPt);
  fHistDeltaPtvsJet2Pt->GetXaxis()->SetTitle("p_{T,2}");  
  fHistDeltaPtvsJet2Pt->GetYaxis()->SetTitle("#deltap_{T} (GeV/c)");
  fHistDeltaPtvsJet2Pt->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistDeltaPtvsJet2Pt);

  fHistDeltaPtOverJet1PtvsJet1Pt = new TH2F("fHistDeltaPtOverJet1PtvsJet1Pt", "fHistDeltaPtOverJet1PtvsJet1Pt", 
					    fNbins, fMinBinPt, fMaxBinPt, fNbins, -5, 5);
  fHistDeltaPtOverJet1PtvsJet1Pt->GetXaxis()->SetTitle("p_{T,1}");  
  fHistDeltaPtOverJet1PtvsJet1Pt->GetYaxis()->SetTitle("#deltap_{T} / p_{T,1}");
  fHistDeltaPtOverJet1PtvsJet1Pt->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistDeltaPtOverJet1PtvsJet1Pt);

  fHistDeltaPtOverJet2PtvsJet2Pt = new TH2F("fHistDeltaPtOverJet2PtvsJet2Pt", "fHistDeltaPtOverJet2PtvsJet2Pt", 
					    fNbins, fMinBinPt, fMaxBinPt, fNbins, -5, 5);
  fHistDeltaPtOverJet2PtvsJet2Pt->GetXaxis()->SetTitle("p_{T,2}");  
  fHistDeltaPtOverJet2PtvsJet2Pt->GetYaxis()->SetTitle("#deltap_{T} / p_{T,2}");
  fHistDeltaPtOverJet2PtvsJet2Pt->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistDeltaPtOverJet2PtvsJet2Pt);

  fHistDeltaPtvsDistance = new TH2F("fHistDeltaPtvsDistance", "fHistDeltaPtvsDistance", 
				    fNbins/2, 0, 1.2, 2*fNbins, -fMaxBinPt, fMaxBinPt);
  fHistDeltaPtvsDistance->GetXaxis()->SetTitle("Distance");  
  fHistDeltaPtvsDistance->GetYaxis()->SetTitle("#deltap_{T} (GeV/c)");
  fHistDeltaPtvsDistance->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistDeltaPtvsDistance);

  fHistDeltaPtvsCommonEnergy1 = new TH2F("fHistDeltaPtvsCommonEnergy1", "fHistDeltaPtvsCommonEnergy1",
					 fNbins/2, 0, 1.2, 2*fNbins, -fMaxBinPt, fMaxBinPt);
  fHistDeltaPtvsCommonEnergy1->GetXaxis()->SetTitle("Common energy 1 (%)");  
  fHistDeltaPtvsCommonEnergy1->GetYaxis()->SetTitle("#deltap_{T} (GeV/c)");
  fHistDeltaPtvsCommonEnergy1->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistDeltaPtvsCommonEnergy1);

  fHistDeltaPtvsCommonEnergy2 = new TH2F("fHistDeltaPtvsCommonEnergy2", "fHistDeltaPtvsCommonEnergy2", 
					 fNbins/2, 0, 1.2, 2*fNbins, -fMaxBinPt, fMaxBinPt);
  fHistDeltaPtvsCommonEnergy2->GetXaxis()->SetTitle("Common energy 2 (%)");  
  fHistDeltaPtvsCommonEnergy2->GetYaxis()->SetTitle("#deltap_{T} (GeV/c)");
  fHistDeltaPtvsCommonEnergy2->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistDeltaPtvsCommonEnergy2);

  fHistDeltaPtvsArea1 = new TH2F("fHistDeltaPtvsArea1", "fHistDeltaPtvsArea1", 
				 fNbins/2, 0, 1.5, 2*fNbins, -fMaxBinPt, fMaxBinPt);
  fHistDeltaPtvsArea1->GetXaxis()->SetTitle("A_{jet,1}");
  fHistDeltaPtvsArea1->GetYaxis()->SetTitle("#deltap_{T} (GeV/c)");
  fHistDeltaPtvsArea1->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistDeltaPtvsArea1);

  fHistDeltaPtvsArea2 = new TH2F("fHistDeltaPtvsArea2", "fHistDeltaPtvsArea2", 
				 fNbins/2, 0, 1.5, 2*fNbins, -fMaxBinPt, fMaxBinPt);
  fHistDeltaPtvsArea2->GetXaxis()->SetTitle("A_{jet,2}");
  fHistDeltaPtvsArea2->GetYaxis()->SetTitle("#deltap_{T} (GeV/c)");
  fHistDeltaPtvsArea2->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistDeltaPtvsArea2);

  fHistDeltaPtvsDeltaArea = new TH2F("fHistDeltaPtvsDeltaArea", "fHistDeltaPtvsDeltaArea", 
				     fNbins, -1.+1./fNbins, 1.+1./fNbins, 2*fNbins, -fMaxBinPt, fMaxBinPt);
  fHistDeltaPtvsDeltaArea->GetXaxis()->SetTitle("#deltaA_{jet}");
  fHistDeltaPtvsDeltaArea->GetYaxis()->SetTitle("#deltap_{T} (GeV/c)");
  fHistDeltaPtvsDeltaArea->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistDeltaPtvsDeltaArea);

  fHistJet1PtvsJet2Pt = new TH2F("fHistJet1PtvsJet2Pt", "fHistJet1PtvsJet2Pt", fNbins, fMinBinPt, fMaxBinPt, fNbins, fMinBinPt, fMaxBinPt);
  fHistJet1PtvsJet2Pt->GetXaxis()->SetTitle("p_{T,1}");
  fHistJet1PtvsJet2Pt->GetYaxis()->SetTitle("p_{T,2}");
  fHistJet1PtvsJet2Pt->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistJet1PtvsJet2Pt);

  if (!fRhoName.IsNull() || !fRho2Name.IsNull()) {
    if (fRhoName.IsNull()) 
      fHistDeltaCorrPtvsJet1CorrPt = new TH2F("fHistDeltaCorrPtvsJet1CorrPt", "fHistDeltaCorrPtvsJet1CorrPt", 
					      fNbins, fMinBinPt, fMaxBinPt, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    else
      fHistDeltaCorrPtvsJet1CorrPt = new TH2F("fHistDeltaCorrPtvsJet1CorrPt", "fHistDeltaCorrPtvsJet1CorrPt", 
					      2*fNbins, -fMaxBinPt, fMaxBinPt, 2*fNbins, -fMaxBinPt, fMaxBinPt);
      
    fHistDeltaCorrPtvsJet1CorrPt->GetXaxis()->SetTitle("p_{T,1}^{corr}");  
    fHistDeltaCorrPtvsJet1CorrPt->GetYaxis()->SetTitle("#deltap_{T}^{corr} (GeV/c)");
    fHistDeltaCorrPtvsJet1CorrPt->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistDeltaCorrPtvsJet1CorrPt);

    if (fRho2Name.IsNull()) 
      fHistDeltaCorrPtvsJet2CorrPt = new TH2F("fHistDeltaCorrPtvsJet2CorrPt", "fHistDeltaCorrPtvsJet2CorrPt", 
					      fNbins, fMinBinPt, fMaxBinPt, 2*fNbins, -fMaxBinPt, fMaxBinPt);					      
    else
      fHistDeltaCorrPtvsJet2CorrPt = new TH2F("fHistDeltaCorrPtvsJet2CorrPt", "fHistDeltaCorrPtvsJet2CorrPt", 
					      2*fNbins, -fMaxBinPt, fMaxBinPt, 2*fNbins, -fMaxBinPt, fMaxBinPt);

    fHistDeltaCorrPtvsJet2CorrPt->GetXaxis()->SetTitle("p_{T,2}^{corr}");  
    fHistDeltaCorrPtvsJet2CorrPt->GetYaxis()->SetTitle("#deltap_{T}^{corr} (GeV/c)");
    fHistDeltaCorrPtvsJet2CorrPt->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistDeltaCorrPtvsJet2CorrPt);

    fHistDeltaCorrPtOverJet1CorrPtvsJet1CorrPt = new TH2F("fHistDeltaCorrPtOverJet1CorrPtvsJet1CorrPt", "fHistDeltaCorrPtOverJet1CorrPtvsJet1CorrPt", 
							  2*fNbins, -fMaxBinPt, fMaxBinPt, fNbins, -5, 5);
    fHistDeltaCorrPtOverJet1CorrPtvsJet1CorrPt->GetXaxis()->SetTitle("p_{T,1}^{corr}");  
    fHistDeltaCorrPtOverJet1CorrPtvsJet1CorrPt->GetYaxis()->SetTitle("#deltap_{T}^{corr} / p_{T,1}^{corr}");
    fHistDeltaCorrPtOverJet1CorrPtvsJet1CorrPt->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistDeltaCorrPtOverJet1CorrPtvsJet1CorrPt);

    fHistDeltaCorrPtOverJet2CorrPtvsJet2CorrPt = new TH2F("fHistDeltaCorrPtOverJet2CorrPtvsJet2CorrPt", "fHistDeltaCorrPtOverJet2CorrPtvsJet2CorrPt", 
							  2*fNbins, -fMaxBinPt, fMaxBinPt, fNbins, -5, 5);
    fHistDeltaCorrPtOverJet2CorrPtvsJet2CorrPt->GetXaxis()->SetTitle("p_{T,2}^{corr}");  
    fHistDeltaCorrPtOverJet2CorrPtvsJet2CorrPt->GetYaxis()->SetTitle("#deltap_{T}^{corr} / p_{T,2}^{corr}");
    fHistDeltaCorrPtOverJet2CorrPtvsJet2CorrPt->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistDeltaCorrPtOverJet2CorrPtvsJet2CorrPt);

    fHistDeltaCorrPtvsDistance = new TH2F("fHistDeltaCorrPtvsDistance", "fHistDeltaCorrPtvsDistance", 
					  fNbins/2, 0, 1.2, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    fHistDeltaCorrPtvsDistance->GetXaxis()->SetTitle("Distance");  
    fHistDeltaCorrPtvsDistance->GetYaxis()->SetTitle("#deltap_{T}^{corr} (GeV/c)");
    fHistDeltaCorrPtvsDistance->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistDeltaCorrPtvsDistance);

    fHistDeltaCorrPtvsCommonEnergy1 = new TH2F("fHistDeltaCorrPtvsCommonEnergy1", "fHistDeltaCorrPtvsCommonEnergy1", 
					       fNbins/2, 0, 1.2, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    fHistDeltaCorrPtvsCommonEnergy1->GetXaxis()->SetTitle("Common energy 1 (%)");  
    fHistDeltaCorrPtvsCommonEnergy1->GetYaxis()->SetTitle("#deltap_{T}^{corr} (GeV/c)");
    fHistDeltaCorrPtvsCommonEnergy1->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistDeltaCorrPtvsCommonEnergy1);

    fHistDeltaCorrPtvsCommonEnergy2 = new TH2F("fHistDeltaCorrPtvsCommonEnergy2", "fHistDeltaCorrPtvsCommonEnergy2", 
					       fNbins/2, 0, 1.2, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    fHistDeltaCorrPtvsCommonEnergy2->GetXaxis()->SetTitle("Common energy 2 (%)");  
    fHistDeltaCorrPtvsCommonEnergy2->GetYaxis()->SetTitle("#deltap_{T}^{corr} (GeV/c)");
    fHistDeltaCorrPtvsCommonEnergy2->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistDeltaCorrPtvsCommonEnergy2);

    fHistDeltaCorrPtvsArea1 = new TH2F("fHistDeltaCorrPtvsArea1", "fHistDeltaCorrPtvsArea1", 
				       fNbins/2, 0, 1.5, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    fHistDeltaCorrPtvsArea1->GetXaxis()->SetTitle("A_{jet,1}");
    fHistDeltaCorrPtvsArea1->GetYaxis()->SetTitle("#deltap_{T}^{corr} (GeV/c)");
    fHistDeltaCorrPtvsArea1->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistDeltaCorrPtvsArea1);
    
    fHistDeltaCorrPtvsArea2 = new TH2F("fHistDeltaCorrPtvsArea2", "fHistDeltaCorrPtvsArea2", 
				       fNbins/2, 0, 1.5, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    fHistDeltaCorrPtvsArea2->GetXaxis()->SetTitle("A_{jet,2}");
    fHistDeltaCorrPtvsArea2->GetYaxis()->SetTitle("#deltap_{T}^{corr} (GeV/c)");
    fHistDeltaCorrPtvsArea2->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistDeltaCorrPtvsArea2);
    
    fHistDeltaCorrPtvsDeltaArea = new TH2F("fHistDeltaCorrPtvsDeltaArea", "fHistDeltaCorrPtvsDeltaArea", 
					   fNbins, -1.+1./fNbins, 1.+1./fNbins, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    fHistDeltaCorrPtvsDeltaArea->GetXaxis()->SetTitle("#deltaA_{jet}");
    fHistDeltaCorrPtvsDeltaArea->GetYaxis()->SetTitle("#deltap_{T}^{corr} (GeV/c)");
    fHistDeltaCorrPtvsDeltaArea->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistDeltaCorrPtvsDeltaArea);
    
    if (fRhoName.IsNull()) 
      fHistJet1CorrPtvsJet2CorrPt = new TH2F("fHistJet1CorrPtvsJet2CorrPt", "fHistJet1CorrPtvsJet2CorrPt", 
					     fNbins, fMinBinPt, fMaxBinPt, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    else if (fRho2Name.IsNull()) 
      fHistJet1CorrPtvsJet2CorrPt = new TH2F("fHistJet1CorrPtvsJet2CorrPt", "fHistJet1CorrPtvsJet2CorrPt", 
					     2*fNbins, -fMaxBinPt, fMaxBinPt, fNbins, fMinBinPt, fMaxBinPt);
    else
      fHistJet1CorrPtvsJet2CorrPt = new TH2F("fHistJet1CorrPtvsJet2CorrPt", "fHistJet1CorrPtvsJet2CorrPt", 
					     2*fNbins, -fMaxBinPt, fMaxBinPt, 
					     2*fNbins, -fMaxBinPt, fMaxBinPt);
    fHistJet1CorrPtvsJet2CorrPt->GetXaxis()->SetTitle("p_{T,1}^{corr}");
    fHistJet1CorrPtvsJet2CorrPt->GetYaxis()->SetTitle("p_{T,2}^{corr}");
    fHistJet1CorrPtvsJet2CorrPt->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistJet1CorrPtvsJet2CorrPt);
  }

  if (fIsEmbedded) {
    fHistDeltaMCPtvsJet1MCPt = new TH2F("fHistDeltaMCPtvsJet1MCPt", "fHistDeltaMCPtvsJet1MCPt", 
					fNbins, fMinBinPt, fMaxBinPt, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    fHistDeltaMCPtvsJet1MCPt->GetXaxis()->SetTitle("p_{T,1}^{MC}");  
    fHistDeltaMCPtvsJet1MCPt->GetYaxis()->SetTitle("#deltap_{T} (GeV/c)");
    fHistDeltaMCPtvsJet1MCPt->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistDeltaMCPtvsJet1MCPt);

    fHistDeltaMCPtvsJet2Pt = new TH2F("fHistDeltaMCPtvsJet2Pt", "fHistDeltaMCPtvsJet2Pt", 
				      fNbins, fMinBinPt, fMaxBinPt, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    fHistDeltaMCPtvsJet2Pt->GetXaxis()->SetTitle("p_{T,2}");  
    fHistDeltaMCPtvsJet2Pt->GetYaxis()->SetTitle("#deltap_{T} (GeV/c)");
    fHistDeltaMCPtvsJet2Pt->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistDeltaMCPtvsJet2Pt);

    fHistDeltaMCPtOverJet1MCPtvsJet1MCPt = new TH2F("fHistDeltaMCPtOverJet1MCPtvsJet1MCPt", "fHistDeltaMCPtOverJet1MCPtvsJet1MCPt", 
						    fNbins, fMinBinPt, fMaxBinPt, fNbins, -5, 5);
    fHistDeltaMCPtOverJet1MCPtvsJet1MCPt->GetXaxis()->SetTitle("p_{T,1}^{MC}");  
    fHistDeltaMCPtOverJet1MCPtvsJet1MCPt->GetYaxis()->SetTitle("#deltap_{T} / p_{T,1}^{MC}");
    fHistDeltaMCPtOverJet1MCPtvsJet1MCPt->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistDeltaMCPtOverJet1MCPtvsJet1MCPt);

    fHistDeltaMCPtOverJet2PtvsJet2Pt = new TH2F("fHistDeltaMCPtOverJet2PtvsJet2Pt", "fHistDeltaMCPtOverJet2PtvsJet2Pt", 
						fNbins, fMinBinPt, fMaxBinPt, fNbins, -5, 5);
    fHistDeltaMCPtOverJet2PtvsJet2Pt->GetXaxis()->SetTitle("p_{T,2}");  
    fHistDeltaMCPtOverJet2PtvsJet2Pt->GetYaxis()->SetTitle("#deltap_{T} / p_{T,2}");
    fHistDeltaMCPtOverJet2PtvsJet2Pt->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistDeltaMCPtOverJet2PtvsJet2Pt);

    fHistDeltaMCPtvsDistance = new TH2F("fHistDeltaMCPtvsDistance", "fHistDeltaMCPtvsDistance", 
					fNbins/2, 0, 1.2, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    fHistDeltaMCPtvsDistance->GetXaxis()->SetTitle("Distance");  
    fHistDeltaMCPtvsDistance->GetYaxis()->SetTitle("#deltap_{T} (GeV/c)");
    fHistDeltaMCPtvsDistance->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistDeltaMCPtvsDistance);

    fHistDeltaMCPtvsCommonEnergy1 = new TH2F("fHistDeltaMCPtvsCommonEnergy1", "fHistDeltaMCPtvsCommonEnergy1", 
					     fNbins/2, 0, 1.2, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    fHistDeltaMCPtvsCommonEnergy1->GetXaxis()->SetTitle("Common energy 1 (%)");  
    fHistDeltaMCPtvsCommonEnergy1->GetYaxis()->SetTitle("#deltap_{T} (GeV/c)");
    fHistDeltaMCPtvsCommonEnergy1->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistDeltaMCPtvsCommonEnergy1);

    fHistDeltaMCPtvsCommonEnergy2 = new TH2F("fHistDeltaMCPtvsCommonEnergy2", "fHistDeltaMCPtvsCommonEnergy2", 
					     fNbins/2, 0, 1.2, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    fHistDeltaMCPtvsCommonEnergy2->GetXaxis()->SetTitle("Common energy 2 (%)");  
    fHistDeltaMCPtvsCommonEnergy2->GetYaxis()->SetTitle("#deltap_{T} (GeV/c)");
    fHistDeltaMCPtvsCommonEnergy2->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistDeltaMCPtvsCommonEnergy2);

    fHistDeltaMCPtvsArea1 = new TH2F("fHistDeltaMCPtvsArea1", "fHistDeltaMCPtvsArea1", 
				     fNbins/2, 0, 1.5, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    fHistDeltaMCPtvsArea1->GetXaxis()->SetTitle("A_{jet,1}");
    fHistDeltaMCPtvsArea1->GetYaxis()->SetTitle("#deltap_{T} (GeV/c)");
    fHistDeltaMCPtvsArea1->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistDeltaMCPtvsArea1);

    fHistDeltaMCPtvsArea2 = new TH2F("fHistDeltaMCPtvsArea2", "fHistDeltaMCPtvsArea2", 
				     fNbins/2, 0, 1.5, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    fHistDeltaMCPtvsArea2->GetXaxis()->SetTitle("A_{jet,2}");
    fHistDeltaMCPtvsArea2->GetYaxis()->SetTitle("#deltap_{T} (GeV/c)");
    fHistDeltaMCPtvsArea2->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistDeltaMCPtvsArea2);

    fHistDeltaMCPtvsDeltaArea = new TH2F("fHistDeltaMCPtvsDeltaArea", "fHistDeltaMCPtvsDeltaArea", 
					 fNbins, -1.+1./fNbins, 1.+1./fNbins, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    fHistDeltaMCPtvsDeltaArea->GetXaxis()->SetTitle("#deltaA_{jet}");
    fHistDeltaMCPtvsDeltaArea->GetYaxis()->SetTitle("#deltap_{T} (GeV/c)");
    fHistDeltaMCPtvsDeltaArea->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistDeltaMCPtvsDeltaArea);

    fHistJet1MCPtvsJet2Pt = new TH2F("fHistJet1MCPtvsJet2Pt", "fHistJet1MCPtvsJet2Pt", 
				     fNbins, fMinBinPt, fMaxBinPt, fNbins, fMinBinPt, fMaxBinPt);
    fHistJet1MCPtvsJet2Pt->GetXaxis()->SetTitle("p_{T,1}^{MC}");
    fHistJet1MCPtvsJet2Pt->GetYaxis()->SetTitle("p_{T,2}");
    fHistJet1MCPtvsJet2Pt->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistJet1MCPtvsJet2Pt);
  }
}

//________________________________________________________________________
void AliJetResponseMaker::AllocateTHnSparse()
{
  // Allocate THnSparse histograms.

  TString title[20]= {""};
  Int_t nbins[20]  = {0};
  Double_t min[20] = {0.};
  Double_t max[20] = {0.};
  Int_t dim = 0;

  title[dim] = "#phi";
  nbins[dim] = fNbins/4;
  min[dim] = 0;
  max[dim] = 2*TMath::Pi()*(1 + 1./(nbins[dim]-1));
  dim++;

  title[dim] = "#eta";
  nbins[dim] = fNbins/4;
  min[dim] = -1;
  max[dim] = 1;
  dim++;

  title[dim] = "p_{T}";
  nbins[dim] = fNbins;
  min[dim] = 0;
  max[dim] = 250;
  dim++;

  title[dim] = "A_{jet}";
  nbins[dim] = fNbins/4;
  min[dim] = 0;
  max[dim] = 1.5;
  dim++;

  title[dim] = "NEF";
  nbins[dim] = fNbins/4;
  min[dim] = 0;
  max[dim] = 1.2;
  dim++;

  title[dim] = "Z";
  nbins[dim] = fNbins/4;
  min[dim] = 0;
  max[dim] = 1.2;
  dim++;

  Int_t dim1 = dim, dim2 = dim;

  if (!fRhoName.IsNull()) {
    title[dim1] = "p_{T}^{corr}";
    nbins[dim1] = fNbins*2;
    min[dim1] = -250;
    max[dim1] = 250;
    dim1++;
  }

  if (fIsEmbedded) {
    title[dim1] = "p_{T}^{MC}";
    nbins[dim1] = fNbins;
    min[dim1] = 0;
    max[dim1] = 250;
    dim1++;
  }

  fHistJets1 = new THnSparseD("fHistJets1","fHistJets1",dim1,nbins,min,max);
  for (Int_t i = 0; i < dim1; i++) 
    fHistJets1->GetAxis(i)->SetTitle(title[i]);
  fOutput->Add(fHistJets1);

  if (!fRho2Name.IsNull()) {
    title[dim2] = "p_{T}^{corr}";
    nbins[dim2] = fNbins*2;
    min[dim2] = -250;
    max[dim2] = 250;
    dim2++;
  }

  if (fDoJet2Histogram) {
    fHistJets2 = new THnSparseD("fHistJets2","fHistJets2",dim2,nbins,min,max);
    for (Int_t i = 0; i < dim2; i++) 
      fHistJets2->GetAxis(i)->SetTitle(title[i]);
    fOutput->Add(fHistJets2);
  }

  fHistJets2Acceptance = new THnSparseD("fHistJets2Acceptance","fHistJets2Acceptance",dim2,nbins,min,max);
  for (Int_t i = 0; i < dim2; i++) 
    fHistJets2Acceptance->GetAxis(i)->SetTitle(title[i]);
  fOutput->Add(fHistJets2Acceptance);
  
  // Matching

  dim = 0;

  title[dim] = "p_{T,1}";
  nbins[dim] = fNbins;
  min[dim] = 0;
  max[dim] = 250;
  dim++;

  title[dim] = "p_{T,2}";
  nbins[dim] = fNbins;
  min[dim] = 0;
  max[dim] = 250;
  dim++;

  title[dim] = "A_{jet,1}";
  nbins[dim] = fNbins/4;
  min[dim] = 0;
  max[dim] = 1.5;
  dim++;

  title[dim] = "A_{jet,2}";
  nbins[dim] = fNbins/4;
  min[dim] = 0;
  max[dim] = 1.5;
  dim++;

  title[dim] = "distance";
  nbins[dim] = fNbins/2;
  min[dim] = 0;
  max[dim] = 1.2;
  dim++;

  title[dim] = "CE1";
  nbins[dim] = fNbins/2;
  min[dim] = 0;
  max[dim] = 1.2;
  dim++;

  title[dim] = "CE2";
  nbins[dim] = fNbins/2;
  min[dim] = 0;
  max[dim] = 1.2;
  dim++;

  if (fDeltaPtAxis) {
    title[dim] = "#deltaA_{jet}";
    nbins[dim] = fNbins/2;
    min[dim] = -1;
    max[dim] = 1;
    dim++;

    title[dim] = "#deltap_{T}";
    nbins[dim] = fNbins*2;
    min[dim] = -250;
    max[dim] = 250;
    dim++;
  }
  if (!fRhoName.IsNull()) {
    title[dim] = "p_{T,1}^{corr}";
    nbins[dim] = fNbins*2;
    min[dim] = -250;
    max[dim] = 250;
    dim++;
  }
  if (!fRho2Name.IsNull()) {
    title[dim] = "p_{T,2}^{corr}";
    nbins[dim] = fNbins*2;
    min[dim] = -250;
    max[dim] = 250;
    dim++;
  }
  if (fDeltaPtAxis && (!fRhoName.IsNull() || !fRho2Name.IsNull())) {
    title[dim] = "#deltap_{T}^{corr}";
    nbins[dim] = fNbins*2;
    min[dim] = -250;
    max[dim] = 250;
    dim++;
  }
  if (fDeltaEtaDeltaPhiAxis) {
    title[dim] = "#delta#eta";
    nbins[dim] = fNbins/2;
    min[dim] = -1;
    max[dim] = 1;
    dim++;

    title[dim] = "#delta#phi";
    nbins[dim] = fNbins/2;
    min[dim] = -TMath::Pi()/2;
    max[dim] = TMath::Pi()*3/2;
    dim++;
  }
  if (fIsEmbedded) {
    title[dim] = "p_{T,1}^{MC}";
    nbins[dim] = fNbins;
    min[dim] = 0;
    max[dim] = 250;
    dim++;

    if (fDeltaPtAxis) {
      title[dim] = "#deltap_{T}^{MC}";
      nbins[dim] = fNbins*2;
      min[dim] = -250;
      max[dim] = 250;
      dim++;
    }
  }

  if (fNEFAxis) {
    title[dim] = "NEF_{1}";
    nbins[dim] = fNbins/4;
    min[dim] = 0;
    max[dim] = 1.2;
    dim++;

    title[dim] = "NEF_{2}";
    nbins[dim] = fNbins/4;
    min[dim] = 0;
    max[dim] = 1.2;
    dim++;
  }

  if (fZAxis) {
    title[dim] = "Z_{1}";
    nbins[dim] = fNbins/4;
    min[dim] = 0;
    max[dim] = 1.2;
    dim++;

    title[dim] = "Z_{2}";
    nbins[dim] = fNbins/4;
    min[dim] = 0;
    max[dim] = 1.2;
    dim++;
  }
      
  fHistMatching = new THnSparseD("fHistMatching","fHistMatching",dim,nbins,min,max);
    
  for (Int_t i = 0; i < dim; i++)
    fHistMatching->GetAxis(i)->SetTitle(title[i]);
   
  fOutput->Add(fHistMatching);
}

//________________________________________________________________________
void AliJetResponseMaker::UserCreateOutputObjects()
{
  // Create user objects.

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  // Jets 1 histograms
  fHistLeadingJets1PtArea = new TH2F("fHistLeadingJets1PtArea", "fHistLeadingJets1PtArea", 
				     fNbins/2, 0, 1, fNbins, fMinBinPt, fMaxBinPt);
  fHistLeadingJets1PtArea->GetXaxis()->SetTitle("area");
  fHistLeadingJets1PtArea->GetYaxis()->SetTitle("p_{T,1} (GeV/c)");
  fHistLeadingJets1PtArea->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistLeadingJets1PtArea);

  if (!fRhoName.IsNull()) {
    fHistLeadingJets1CorrPtArea = new TH2F("fHistLeadingJets1CorrPtArea", "fHistLeadingJets1CorrPtArea", 
					   fNbins/2, 0, 1, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    fHistLeadingJets1CorrPtArea->GetXaxis()->SetTitle("area");
    fHistLeadingJets1CorrPtArea->GetYaxis()->SetTitle("p_{T,1}^{corr} (GeV/c)");
    fHistLeadingJets1CorrPtArea->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistLeadingJets1CorrPtArea);
  }

  fHistLeadingJets2PtAreaAcceptance = new TH2F("fHistLeadingJets2PtAreaAcceptance", "fHistLeadingJets2PtAreaAcceptance", 
					       fNbins/2, 1, 2, fNbins, fMinBinPt, fMaxBinPt);
  fHistLeadingJets2PtAreaAcceptance->GetXaxis()->SetTitle("area");
  fHistLeadingJets2PtAreaAcceptance->GetYaxis()->SetTitle("p_{T,2} (GeV/c)");
  fHistLeadingJets2PtAreaAcceptance->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistLeadingJets2PtAreaAcceptance);

  fHistLeadingJets2PtArea = new TH2F("fHistLeadingJets2PtArea", "fHistLeadingJets2PtArea", 
				     fNbins/2, 0, 1, fNbins, fMinBinPt, fMaxBinPt);
  fHistLeadingJets2PtArea->GetXaxis()->SetTitle("area");
  fHistLeadingJets2PtArea->GetYaxis()->SetTitle("p_{T,2} (GeV/c)");
  fHistLeadingJets2PtArea->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistLeadingJets2PtArea);

  if (!fRho2Name.IsNull()) {
    fHistLeadingJets2CorrPtAreaAcceptance = new TH2F("fHistLeadingJets2CorrPtAreaAcceptance", "fHistLeadingJets2CorrPtAreaAcceptance", 
						     fNbins/2, 0, 1, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    fHistLeadingJets2CorrPtAreaAcceptance->GetXaxis()->SetTitle("area");
    fHistLeadingJets2CorrPtAreaAcceptance->GetYaxis()->SetTitle("p_{T,2}^{corr} (GeV/c)");
    fHistLeadingJets2CorrPtAreaAcceptance->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistLeadingJets2CorrPtAreaAcceptance);

    fHistLeadingJets2CorrPtArea = new TH2F("fHistLeadingJets2CorrPtArea", "fHistLeadingJets2CorrPtArea", 
					   fNbins/2, 0, 1, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    fHistLeadingJets2CorrPtArea->GetXaxis()->SetTitle("area");
    fHistLeadingJets2CorrPtArea->GetYaxis()->SetTitle("p_{T,2}^{corr} (GeV/c)");
    fHistLeadingJets2CorrPtArea->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistLeadingJets2CorrPtArea);
  }

  if (fHistoType==0)
    AllocateTH2();
  else 
    AllocateTHnSparse();

  PostData(1, fOutput); // Post data for ALL output slots > 0 here, to get at least an empty histogram
}

//________________________________________________________________________
void AliJetResponseMaker::FillJetHisto(Double_t Phi, Double_t Eta, Double_t Pt, Double_t A, Double_t NEF, Double_t Z, Double_t CorrPt, Double_t MCPt, Int_t Set)
{
  if (fHistoType==1) {
    THnSparse *histo = 0;
    if (Set==1)
      histo = fHistJets1;
    else if (Set==2)
      histo = fHistJets2;
    else if (Set==3)
      histo = fHistJets2Acceptance;

    if (!histo)
      return;

    Double_t contents[20]={0};

    for (Int_t i = 0; i < histo->GetNdimensions(); i++) {
      TString title(histo->GetAxis(i)->GetTitle());
      if (title=="#phi")
	contents[i] = Phi;
      else if (title=="#eta")
	contents[i] = Eta;
      else if (title=="p_{T}")
	contents[i] = Pt;
      else if (title=="A_{jet}")
	contents[i] = A;
      else if (title=="NEF")
	contents[i] = NEF;
      else if (title=="Z")
	contents[i] = Z;
      else if (title=="p_{T}^{corr}")
	contents[i] = CorrPt;
      else if (title=="p_{T}^{MC}")
	contents[i] = MCPt;
      else 
	AliWarning(Form("Unable to fill dimension %s!",title.Data()));
    }

    histo->Fill(contents);
  }
  else {
    if (Set == 1) {
      fHistJets1PtArea->Fill(A, Pt);
      fHistJets1PhiEta->Fill(Eta, Phi);
      fHistJets1ZvsPt->Fill(Z, Pt);
      fHistJets1NEFvsPt->Fill(NEF, Pt);
      fHistJets1CEFvsCEFPt->Fill(1-NEF, (1-NEF)*Pt);
      if (!fRhoName.IsNull())
	fHistJets1CorrPtArea->Fill(A, CorrPt);
    }
    else if (Set == 2) {
      fHistJets2PtArea->Fill(A, Pt);
      fHistJets2PhiEta->Fill(Eta, Phi);
      if (!fRho2Name.IsNull())
	fHistJets2CorrPtArea->Fill(A, CorrPt);
    }
    else if (Set == 3) {
      fHistJets2PtAreaAcceptance->Fill(A, Pt);
      fHistJets2PhiEtaAcceptance->Fill(Eta, Phi);
      fHistJets2ZvsPt->Fill(Z, Pt);
      fHistJets2NEFvsPt->Fill(NEF, Pt);
      fHistJets2CEFvsCEFPt->Fill(1-NEF, (1-NEF)*Pt);
      if (!fRho2Name.IsNull())
	fHistJets2CorrPtAreaAcceptance->Fill(A, CorrPt);
    }
  }
}

//________________________________________________________________________
void AliJetResponseMaker::FillMatchingHistos(Double_t Pt1, Double_t Pt2, Double_t Eta1, Double_t Eta2, Double_t Phi1, Double_t Phi2, 
					     Double_t A1, Double_t A2, Double_t d, Double_t CE1, Double_t CE2, Double_t CorrPt1, Double_t CorrPt2, 
					     Double_t MCPt1, Double_t NEF1, Double_t NEF2, Double_t Z1, Double_t Z2)
{
  if (fHistoType==1) {
    Double_t contents[20]={0};

    for (Int_t i = 0; i < fHistMatching->GetNdimensions(); i++) {
      TString title(fHistMatching->GetAxis(i)->GetTitle());
      if (title=="p_{T,1}")
	contents[i] = Pt1;
      else if (title=="p_{T,2}")
	contents[i] = Pt2;
      else if (title=="A_{jet,1}")
	contents[i] = A1;
      else if (title=="A_{jet,2}")
	contents[i] = A2;
      else if (title=="distance")
	contents[i] = d;
      else if (title=="CE1")
	contents[i] = CE1;
      else if (title=="CE2")
	contents[i] = CE2;
      else if (title=="#deltaA_{jet}")
	contents[i] = A1-A2;
      else if (title=="#deltap_{T}")
	contents[i] = Pt1-Pt2;
      else if (title=="#delta#eta")
	contents[i] = Eta1-Eta2;
      else if (title=="#delta#phi")
	contents[i] = Phi1-Phi2;
      else if (title=="p_{T,1}^{corr}")
	contents[i] = CorrPt1;
      else if (title=="p_{T,2}^{corr}")
	contents[i] = CorrPt2;
      else if (title=="#deltap_{T}^{corr}")
	contents[i] = CorrPt1-CorrPt2;
      else if (title=="p_{T,1}^{MC}")
	contents[i] = MCPt1;
      else if (title=="#deltap_{T}^{MC}")
	contents[i] = MCPt1-Pt2;
      else if (title=="NEF_{1}")
	contents[i] = NEF1;
      else if (title=="NEF_{2}")
	contents[i] = NEF2;
      else if (title=="Z_{1}")
	contents[i] = Z1;
      else if (title=="Z_{2}")
	contents[i] = Z2;
      else 
	AliWarning(Form("Unable to fill dimension %s!",title.Data()));
    }

    fHistMatching->Fill(contents);
  }
  else {
    fHistCommonEnergy1vsJet1Pt->Fill(CE1, Pt1);
    fHistCommonEnergy2vsJet2Pt->Fill(CE2, Pt2);

    fHistDistancevsJet1Pt->Fill(d, Pt1);
    fHistDistancevsJet2Pt->Fill(d, Pt2);

    fHistDistancevsCommonEnergy1->Fill(d, CE1);
    fHistDistancevsCommonEnergy2->Fill(d, CE2);
    fHistCommonEnergy1vsCommonEnergy2->Fill(CE1, CE2);

    fHistDeltaEtaDeltaPhi->Fill(Eta1-Eta2,Phi1-Phi2);

    fHistJet2PtOverJet1PtvsJet2Pt->Fill(Pt2, Pt2 / Pt1);
    fHistJet1PtOverJet2PtvsJet1Pt->Fill(Pt1, Pt1 / Pt2);

    Double_t dpt = Pt1 - Pt2;
    fHistDeltaPtvsJet1Pt->Fill(Pt1, dpt);
    fHistDeltaPtvsJet2Pt->Fill(Pt2, dpt);
    fHistDeltaPtOverJet1PtvsJet1Pt->Fill(Pt1, dpt/Pt1);
    fHistDeltaPtOverJet2PtvsJet2Pt->Fill(Pt2, dpt/Pt2);

    fHistDeltaPtvsDistance->Fill(d, dpt);
    fHistDeltaPtvsCommonEnergy1->Fill(CE1, dpt);
    fHistDeltaPtvsCommonEnergy2->Fill(CE2, dpt);

    fHistDeltaPtvsArea1->Fill(A1, dpt);
    fHistDeltaPtvsArea2->Fill(A2, dpt);

    Double_t darea = A1 - A2;
    fHistDeltaPtvsDeltaArea->Fill(darea, dpt);

    fHistJet1PtvsJet2Pt->Fill(Pt1, Pt2);

    if (!fRhoName.IsNull() || !fRho2Name.IsNull()) {
      Double_t dcorrpt = CorrPt1 - CorrPt2;
      fHistDeltaCorrPtvsJet1CorrPt->Fill(CorrPt1, dcorrpt);
      fHistDeltaCorrPtvsJet2CorrPt->Fill(CorrPt2, dcorrpt);
      fHistDeltaCorrPtOverJet1CorrPtvsJet1CorrPt->Fill(CorrPt1, dcorrpt/CorrPt1);
      fHistDeltaCorrPtOverJet2CorrPtvsJet2CorrPt->Fill(CorrPt2, dcorrpt/CorrPt2);
      fHistDeltaCorrPtvsDistance->Fill(d, dcorrpt);
      fHistDeltaCorrPtvsCommonEnergy1->Fill(CE1, dcorrpt);
      fHistDeltaCorrPtvsCommonEnergy2->Fill(CE2, dcorrpt);
      fHistDeltaCorrPtvsArea1->Fill(A1, dcorrpt);
      fHistDeltaCorrPtvsArea2->Fill(A2, dcorrpt);
      fHistDeltaCorrPtvsDeltaArea->Fill(darea, dcorrpt);
      fHistJet1CorrPtvsJet2CorrPt->Fill(CorrPt1, CorrPt2);
    }

    if (fIsEmbedded) {
      Double_t dmcpt = MCPt1 - Pt2;
      fHistDeltaMCPtvsJet1MCPt->Fill(MCPt1, dmcpt);
      fHistDeltaMCPtvsJet2Pt->Fill(Pt2, dmcpt);
      fHistDeltaMCPtOverJet1MCPtvsJet1MCPt->Fill(MCPt1, dmcpt/MCPt1);
      fHistDeltaMCPtOverJet2PtvsJet2Pt->Fill(Pt2, dmcpt/Pt2);
      fHistDeltaMCPtvsDistance->Fill(d, dmcpt);
      fHistDeltaMCPtvsCommonEnergy1->Fill(CE1, dmcpt);
      fHistDeltaMCPtvsCommonEnergy2->Fill(CE2, dmcpt);
      fHistDeltaMCPtvsArea1->Fill(A1, dmcpt);
      fHistDeltaMCPtvsArea2->Fill(A2, dmcpt);
      fHistDeltaMCPtvsDeltaArea->Fill(darea, dmcpt);
      fHistJet1MCPtvsJet2Pt->Fill(MCPt1, Pt2);
    }
  }
}

//________________________________________________________________________
Bool_t AliJetResponseMaker::AcceptJet(AliEmcalJet *jet) const
{   
  // Return true if jet is accepted.

  if (jet->Pt() < fJetPtCut)
    return kFALSE;
  if (jet->Area() < fJetAreaCut)
    return kFALSE;
  if (jet->MCPt() < fMinJetMCPt)
    return kFALSE;

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliJetResponseMaker::AcceptBiasJet2(AliEmcalJet *jet) const
{ 
  // Accept jet with a bias.

  if (fLeadingHadronType == 0) {
    if (jet->MaxTrackPt() < fPtBiasJet2Track) return kFALSE;
  }
  else if (fLeadingHadronType == 1) {
    if (jet->MaxClusterPt() < fPtBiasJet2Clus) return kFALSE;
  }
  else {
    if (jet->MaxTrackPt() < fPtBiasJet2Track && jet->MaxClusterPt() < fPtBiasJet2Clus) return kFALSE;
  }

  return kTRUE;
}

//________________________________________________________________________
void AliJetResponseMaker::ExecOnce()
{
  // Execute once.

  if (!fJets2Name.IsNull() && !fJets2) {
    fJets2 = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fJets2Name));
    if (!fJets2) {
      AliError(Form("%s: Could not retrieve jets2 %s!", GetName(), fJets2Name.Data()));
      return;
    }
    else if (!fJets2->GetClass()->GetBaseClass("AliEmcalJet")) {
      AliError(Form("%s: Collection %s does not contain AliEmcalJet objects!", GetName(), fJets2Name.Data())); 
      fJets2 = 0;
      return;
    }
  }

  if (!fTracks2Name.IsNull() && !fTracks2) {
    fTracks2 = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTracks2Name));
    if (!fTracks2) {
      AliError(Form("%s: Could not retrieve tracks2 %s!", GetName(), fTracks2Name.Data())); 
      return;
    }
    else {
      TClass *cl = fTracks2->GetClass();
      if (!cl->GetBaseClass("AliVParticle") && !cl->GetBaseClass("AliEmcalParticle")) {
	AliError(Form("%s: Collection %s does not contain AliVParticle nor AliEmcalParticle objects!", GetName(), fTracks2Name.Data())); 
	fTracks2 = 0;
	return;
      }
    }

    if (fAreCollections2MC) {
      fTracks2Map = dynamic_cast<AliNamedArrayI*>(InputEvent()->FindListObject(fTracks2Name + "_Map"));
      // this is needed to map the MC labels with the indexes of the MC particle collection
      // if teh map is not given, the MC labels are assumed to be consistent with the indexes (which is not the case if AliEmcalMCTrackSelector is used)
      if (!fTracks2Map) {
	AliWarning(Form("%s: Could not retrieve map for tracks2 %s! Will assume MC labels consistent with indexes...", GetName(), fTracks2Name.Data())); 
	fTracks2Map = new AliNamedArrayI("tracksMap",9999);
	for (Int_t i = 0; i < 9999; i++) {
	  fTracks2Map->AddAt(i,i);
	}
      }
    }
  }

  if (!fCalo2Name.IsNull() && !fCaloClusters2) {
    fCaloClusters2 =  dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fCalo2Name));
    if (!fCaloClusters2) {
      AliError(Form("%s: Could not retrieve clusters %s!", GetName(), fCalo2Name.Data())); 
      return;
    } else {
      TClass *cl = fCaloClusters2->GetClass();
      if (!cl->GetBaseClass("AliVCluster") && !cl->GetBaseClass("AliEmcalParticle")) {
	AliError(Form("%s: Collection %s does not contain AliVCluster nor AliEmcalParticle objects!", GetName(), fCalo2Name.Data())); 
	fCaloClusters2 = 0;
	return;
      }
    }
  }

  if (!fRho2Name.IsNull() && !fRho2) {
    fRho2 = dynamic_cast<AliRhoParameter*>(InputEvent()->FindListObject(fRho2Name));
    if (!fRho2) {
      AliError(Form("%s: Could not retrieve rho %s!", GetName(), fRho2Name.Data()));
      fInitialized = kFALSE;
      return;
    }
  }

  if (fPercAreaCut >= 0) {
    if (fJet2AreaCut >= 0)
      AliInfo(Form("%s: jet 2 area cut will be calculated as a percentage of the average area, given value will be overwritten", GetName()));
    fJet2AreaCut = fPercAreaCut * fJet2Radius * fJet2Radius * TMath::Pi();
  }
  if (fJet2AreaCut < 0)
    fJet2AreaCut = 0;

  if (fJet2MinEta == -999)
    fJet2MinEta = fJetMinEta - fJetRadius;
  if (fJet2MaxEta == -999)
    fJet2MaxEta = fJetMaxEta + fJetRadius;
  if (fJet2MinPhi == -999)
    fJet2MinPhi = fJetMinPhi - fJetRadius;
  if (fJet2MaxPhi == -999)
    fJet2MaxPhi = fJetMaxPhi + fJetRadius;

  if (fMatching == kMCLabel && (!fAreCollections2MC || fAreCollections1MC)) {
    if (fAreCollections2MC == fAreCollections1MC) {
      AliWarning("Changing matching type from MC label to same collection...");
      fMatching = kSameCollections;
    }
    else {
      AliWarning("Changing matching type from MC label to geometrical...");
      fMatching = kGeometrical;
    }
  }
  else if (fMatching == kSameCollections && fAreCollections2MC != fAreCollections1MC) {
    if (fAreCollections2MC && !fAreCollections1MC) {
      AliWarning("Changing matching type from same collection to MC label...");
      fMatching = kMCLabel;
    }
    else {
      AliWarning("Changing matching type from same collection to geometrical...");
      fMatching = kGeometrical;
    }
  }

  AliAnalysisTaskEmcalJet::ExecOnce();
}

//________________________________________________________________________
Bool_t AliJetResponseMaker::RetrieveEventObjects()
{
  // Retrieve event objects.

  if (!AliAnalysisTaskEmcalJet::RetrieveEventObjects())
    return kFALSE;

  if (fRho2)
    fRho2Val = fRho2->GetVal();

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliJetResponseMaker::Run()
{
  // Find the closest jets

  if (fMatching == kNoMatching) 
    return kTRUE;
  else
    return DoJetMatching();
}

//________________________________________________________________________
Bool_t AliJetResponseMaker::DoJetMatching()
{
  DoJetLoop(kFALSE);

  const Int_t nJets = fJets->GetEntriesFast();

  for (Int_t i = 0; i < nJets; i++) {

    AliEmcalJet* jet1 = static_cast<AliEmcalJet*>(fJets->At(i));

    if (!jet1) {
      AliError(Form("Could not receive jet %d", i));
      continue;
    }  

    if (!AcceptJet(jet1))
      continue;

    if (jet1->Eta() < fJetMinEta || jet1->Eta() > fJetMaxEta || jet1->Phi() < fJetMinPhi || jet1->Phi() > fJetMaxPhi)
      continue;

    if (jet1->ClosestJet() && jet1->ClosestJet()->ClosestJet() == jet1 && 
	jet1->ClosestJetDistance() < fMatchingPar1 && jet1->ClosestJet()->ClosestJetDistance() < fMatchingPar2) {    // Matched jet found
      jet1->SetMatchedToClosest(fMatching);
      jet1->ClosestJet()->SetMatchedToClosest(fMatching);
      AliDebug(2,Form("Found matching: jet1 pt = %f, eta = %f, phi = %f, jet2 pt = %f, eta = %f, phi = %f", 
		      jet1->Pt(), jet1->Eta(), jet1->Phi(), 
		      jet1->MatchedJet()->Pt(), jet1->MatchedJet()->Eta(), jet1->MatchedJet()->Phi()));
    }
  }

  return kTRUE;
}

//________________________________________________________________________
void AliJetResponseMaker::DoJetLoop(Bool_t order)
{
  // Do the jet loop.

  TClonesArray *jets1 = 0;
  TClonesArray *jets2 = 0;

  if (order) {
    jets1 = fJets2;
    jets2 = fJets;
  }
  else {
    jets1 = fJets;
    jets2 = fJets2;
  }

  Int_t nJets1 = jets1->GetEntriesFast();
  Int_t nJets2 = jets2->GetEntriesFast();

  for (Int_t j = 0; j < nJets2; j++) {
      
    AliEmcalJet* jet2 = static_cast<AliEmcalJet*>(jets2->At(j));
      
    if (!jet2) {
      AliError(Form("Could not receive jet %d", j));
      continue;
    }

    jet2->ResetMatching();


    if (!AcceptJet(jet2))
      continue;

    if (jet2->Eta() < fJet2MinEta || jet2->Eta() > fJet2MaxEta || jet2->Phi() < fJet2MinPhi || jet2->Phi() > fJet2MaxPhi)
      continue;
  }
    
  for (Int_t i = 0; i < nJets1; i++) {

    AliEmcalJet* jet1 = static_cast<AliEmcalJet*>(jets1->At(i));

    if (!jet1) {
      AliError(Form("Could not receive jet %d", i));
      continue;
    }

    jet1->ResetMatching();

    if (!AcceptJet(jet1))
      continue;

    if (order) {
     if (jet1->Eta() < fJet2MinEta || jet1->Eta() > fJet2MaxEta || jet1->Phi() < fJet2MinPhi || jet1->Phi() > fJet2MaxPhi)
	continue;
    }
    else {
      if (jet1->Eta() < fJetMinEta || jet1->Eta() > fJetMaxEta || jet1->Phi() < fJetMinPhi || jet1->Phi() > fJetMaxPhi)
	continue;
    }

    for (Int_t j = 0; j < nJets2; j++) {
      
      AliEmcalJet* jet2 = static_cast<AliEmcalJet*>(jets2->At(j));
      
      if (!jet2) {
	AliError(Form("Could not receive jet %d", j));
	continue;
      }
      if (!AcceptJet(jet2))
	continue;

      if (order) {
	if (jet2->Eta() < fJetMinEta || jet2->Eta() > fJetMaxEta || jet2->Phi() < fJetMinPhi || jet2->Phi() > fJetMaxPhi)
	  continue;
      }
      else {
	if (jet2->Eta() < fJet2MinEta || jet2->Eta() > fJet2MaxEta || jet2->Phi() < fJet2MinPhi || jet2->Phi() > fJet2MaxPhi)
	  continue;
      }

      SetMatchingLevel(jet1, jet2, fMatching);
    } // jet2 loop

  } // jet1 loop
}

//________________________________________________________________________
void AliJetResponseMaker::GetGeometricalMatchingLevel(AliEmcalJet *jet1, AliEmcalJet *jet2, Double_t &d) const
{
  Double_t deta = jet2->Eta() - jet1->Eta();
  Double_t dphi = jet2->Phi() - jet1->Phi();
  d = TMath::Sqrt(deta * deta + dphi * dphi);
}

//________________________________________________________________________
void AliJetResponseMaker::GetMCLabelMatchingLevel(AliEmcalJet *jet1, AliEmcalJet *jet2, Double_t &d1, Double_t &d2) const
{ 
  // d1 and d2 represent the matching level: 0 = maximum level of matching, 1 = the two jets are completely unrelated
  d1 = jet1->Pt();
  d2 = jet2->Pt();
  Double_t totalPt1 = d1; // the total pt of the reconstructed jet will be cleaned from the background

  for (Int_t iTrack2 = 0; iTrack2 < jet2->GetNumberOfTracks(); iTrack2++) {
    Bool_t track2Found = kFALSE;
    Int_t index2 = jet2->TrackAt(iTrack2);
    for (Int_t iTrack = 0; iTrack < jet1->GetNumberOfTracks(); iTrack++) {
      AliVParticle *track = jet1->TrackAt(iTrack,fTracks);
      if (!track) {
	AliWarning(Form("Could not find track %d!", iTrack));
	continue;
      }
      Int_t MClabel = TMath::Abs(track->GetLabel());
      if (MClabel > fMCLabelShift)
	MClabel -= fMCLabelShift;
      Int_t index = -1;
	  
      if (MClabel == 0) {// this is not a MC particle; remove it completely
	AliDebug(3,Form("Track %d (pT = %f) is not a MC particle (MClabel = %d)!",iTrack,track->Pt(),MClabel));
	totalPt1 -= track->Pt();
	d1 -= track->Pt();
	continue;
      }
      else if (MClabel < fTracks2Map->GetSize()) {
	index = fTracks2Map->At(MClabel);
      }
	  
      if (index < 0) {
	AliDebug(2,Form("Track %d (pT = %f) does not have an associated MC particle (MClabel = %d)!",iTrack,track->Pt(),MClabel));
	continue;
      }

      if (index2 == index) { // found common particle
	track2Found = kTRUE;
	d1 -= track->Pt();
	AliVParticle *MCpart = static_cast<AliVParticle*>(fTracks2->At(index2));
	AliDebug(3,Form("Track %d (pT = %f, eta = %f, phi = %f) is associated with the MC particle %d (pT = %f, eta = %f, phi = %f)!",
			iTrack,track->Pt(),track->Eta(),track->Phi(),MClabel,MCpart->Pt(),MCpart->Eta(),MCpart->Phi()));
	d2 -= MCpart->Pt();
	break;
      }
    }
    for (Int_t iClus = 0; iClus < jet1->GetNumberOfClusters(); iClus++) {
      AliVCluster *clus = jet1->ClusterAt(iClus,fCaloClusters);
      if (!clus) {
	AliWarning(Form("Could not find cluster %d!", iClus));
	continue;
      }
      TLorentzVector part;
      clus->GetMomentum(part, const_cast<Double_t*>(fVertex));
	  
      if (fUseCellsToMatch && fCaloCells) { // if the cell colection is available, look for cells with a matched MC particle
	for (Int_t iCell = 0; iCell < clus->GetNCells(); iCell++) {
	  Int_t cellId = clus->GetCellAbsId(iCell);
	  Double_t cellFrac = clus->GetCellAmplitudeFraction(iCell);

	  Int_t MClabel = TMath::Abs(fCaloCells->GetCellMCLabel(cellId));
	  if (MClabel > fMCLabelShift)
	    MClabel -= fMCLabelShift;
	  Int_t index = -1;
	  
	  if (MClabel == 0) {// this is not a MC particle; remove it completely
	    AliDebug(3,Form("Cell %d (frac = %f) is not a MC particle (MClabel = %d)!",iCell,cellFrac,MClabel));
	    totalPt1 -= part.Pt() * cellFrac;
	    d1 -= part.Pt() * cellFrac;
	    continue;
	  }
	  else if (MClabel < fTracks2Map->GetSize()) {
	    index = fTracks2Map->At(MClabel);
	  }

	  if (index < 0) {
	    AliDebug(3,Form("Cell %d (frac = %f) does not have an associated MC particle (MClabel = %d)!",iCell,cellFrac,MClabel));
	    continue;
	  }
	  if (index2 == index) { // found common particle
	    d1 -= part.Pt() * cellFrac;
		
	    if (!track2Found) {// only if it is not already found among charged tracks (charged particles are most likely already found)
	      AliVParticle *MCpart = static_cast<AliVParticle*>(fTracks2->At(index2));
	      AliDebug(3,Form("Cell %d belonging to cluster %d (pT = %f, eta = %f, phi = %f) is associated with the MC particle %d (pT = %f, eta = %f, phi = %f)!",
			      iCell,iClus,part.Pt(),part.Eta(),part.Phi(),MClabel,MCpart->Pt(),MCpart->Eta(),MCpart->Phi()));		  
	      d2 -= MCpart->Pt() * cellFrac;
	    }
	    break;
	  }
	}
      }
      else { //otherwise look for the first contributor to the cluster, and if matched to a MC label remove it
	Int_t MClabel = TMath::Abs(clus->GetLabel());
	if (MClabel > fMCLabelShift)
	  MClabel -= fMCLabelShift;
	Int_t index = -1;
	    
	if (MClabel == 0) {// this is not a MC particle; remove it completely
	  AliDebug(3,Form("Cluster %d (pT = %f) is not a MC particle (MClabel = %d)!",iClus,part.Pt(),MClabel));
	  totalPt1 -= part.Pt();
	  d1 -= part.Pt();
	  continue;
	}
	else if (MClabel < fTracks2Map->GetSize()) {
	  index = fTracks2Map->At(MClabel);
	}
	 
	if (index < 0) {
	  AliDebug(3,Form("Cluster %d (pT = %f) does not have an associated MC particle (MClabel = %d)!",iClus,part.Pt(),MClabel));
	  continue;
	}
	if (index2 == index) { // found common particle
	  d1 -= part.Pt();

	  if (!track2Found) {// only if it is not already found among charged tracks (charged particles are most likely already found)
	    AliVParticle *MCpart = static_cast<AliVParticle*>(fTracks2->At(index2));
	    AliDebug(3,Form("Cluster %d (pT = %f, eta = %f, phi = %f) is associated with the MC particle %d (pT = %f, eta = %f, phi = %f)!",
			    iClus,part.Pt(),part.Eta(),part.Phi(),MClabel,MCpart->Pt(),MCpart->Eta(),MCpart->Phi()));
		
	    d2 -= MCpart->Pt();
	  }
	  break;
	}
      }
    }
  }

  if (d1 < 0)
    d1 = 0;

  if (d2 < 0)
    d2 = 0;

  if (totalPt1 < 1)
    d1 = -1;
  else
    d1 /= totalPt1;

  if (jet2->Pt() < 1)
    d2 = -1;
  else
    d2 /= jet2->Pt();
}

//________________________________________________________________________
void AliJetResponseMaker::GetSameCollectionsMatchingLevel(AliEmcalJet *jet1, AliEmcalJet *jet2, Double_t &d1, Double_t &d2) const
{ 
  // d1 and d2 represent the matching level: 0 = maximum level of matching, 1 = the two jets are completely unrelated
  d1 = jet1->Pt();
  d2 = jet2->Pt();

  if (fTracks && fTracks2) {

    for (Int_t iTrack2 = 0; iTrack2 < jet2->GetNumberOfTracks(); iTrack2++) {
      Int_t index2 = jet2->TrackAt(iTrack2);
      for (Int_t iTrack = 0; iTrack < jet1->GetNumberOfTracks(); iTrack++) {
	Int_t index = jet1->TrackAt(iTrack);
	if (index2 == index) { // found common particle
	  AliVParticle *part = static_cast<AliVParticle*>(fTracks->At(index));
	  if (!part) {
	    AliWarning(Form("Could not find track %d!", index));
	    continue;
	  }
	  AliVParticle *part2 = static_cast<AliVParticle*>(fTracks2->At(index2));
	  if (!part2) {
	    AliWarning(Form("Could not find track %d!", index2));
	    continue;
	  }

	  d1 -= part->Pt();
	  d2 -= part2->Pt();
	  break;
	}
      }
    }

  }

  if (fCaloClusters && fCaloClusters2) {

    if (fUseCellsToMatch) {
      const Int_t nClus1 = jet1->GetNumberOfClusters();

      Int_t ncells1[nClus1];
      UShort_t *cellsId1[nClus1];
      Double_t *cellsFrac1[nClus1];
      Int_t *sortedIndexes1[nClus1];
      Double_t ptClus1[nClus1];
      for (Int_t iClus1 = 0; iClus1 < nClus1; iClus1++) {
	Int_t index1 = jet1->ClusterAt(iClus1);
	AliVCluster *clus1 =  static_cast<AliVCluster*>(fCaloClusters->At(index1));
	if (!clus1) {
	  AliWarning(Form("Could not find cluster %d!", index1));
	  ncells1[iClus1] = 0;
	  cellsId1[iClus1] = 0;
	  cellsFrac1[iClus1] = 0;
	  sortedIndexes1[iClus1] = 0;
	  ptClus1[iClus1] = 0;
	  continue;
	}
	TLorentzVector part1;
	clus1->GetMomentum(part1, const_cast<Double_t*>(fVertex));

	ncells1[iClus1] = clus1->GetNCells();
	cellsId1[iClus1] = clus1->GetCellsAbsId();
	cellsFrac1[iClus1] = clus1->GetCellsAmplitudeFraction();
	sortedIndexes1[iClus1] = new Int_t[ncells1[iClus1]];
	ptClus1[iClus1] = part1.Pt();

	TMath::Sort(ncells1[iClus1], cellsId1[iClus1], sortedIndexes1[iClus1], kFALSE);
      }
      
      const Int_t nClus2 = jet2->GetNumberOfClusters();

      const Int_t maxNcells2 = 11520;
      Int_t sortedIndexes2[maxNcells2];
      for (Int_t iClus2 = 0; iClus2 < nClus2; iClus2++) {
	Int_t index2 = jet2->ClusterAt(iClus2);
	AliVCluster *clus2 =  static_cast<AliVCluster*>(fCaloClusters2->At(index2));
	if (!clus2) {
	  AliWarning(Form("Could not find cluster %d!", index2));
	  continue;
	}
	Int_t ncells2 = clus2->GetNCells();
	if (ncells2 >= maxNcells2) {
	  AliError(Form("Number of cells in the cluster %d >= %d",ncells2,maxNcells2));
	  continue;
	}
	UShort_t *cellsId2 = clus2->GetCellsAbsId();
	Double_t *cellsFrac2 = clus2->GetCellsAmplitudeFraction();

	TLorentzVector part2;
	clus2->GetMomentum(part2, const_cast<Double_t*>(fVertex));
	Double_t ptClus2 = part2.Pt();

	TMath::Sort(ncells2, cellsId2, sortedIndexes2, kFALSE);

	for (Int_t iClus1 = 0; iClus1 < nClus1; iClus1++) {
	  if (sortedIndexes1[iClus1] == 0)
	    continue;
	  Int_t iCell1 = 0, iCell2 = 0;
	  while (iCell1 < ncells1[iClus1] && iCell2 < ncells2) {
	    if (cellsId1[iClus1][sortedIndexes1[iClus1][iCell1]] == cellsId2[sortedIndexes2[iCell2]]) { // found a common cell
	      d1 -= cellsFrac1[iClus1][sortedIndexes1[iClus1][iCell1]] * ptClus1[iClus1];
	      d2 -= cellsFrac2[sortedIndexes2[iCell2]] * ptClus2;
	      iCell1++;
	      iCell2++;
	    }
	    else if (cellsId1[iClus1][sortedIndexes1[iClus1][iCell1]] > cellsId2[sortedIndexes2[iCell2]]) { 
	      iCell2++;
	    }
	    else {
	      iCell1++;
	    }
	  }
	}
      }
    }
    else {
      for (Int_t iClus2 = 0; iClus2 < jet2->GetNumberOfClusters(); iClus2++) {
	Int_t index2 = jet2->ClusterAt(iClus2);
	for (Int_t iClus = 0; iClus < jet1->GetNumberOfClusters(); iClus++) {
	  Int_t index = jet1->ClusterAt(iClus);
	  if (index2 == index) { // found common particle
	    AliVCluster *clus =  static_cast<AliVCluster*>(fCaloClusters->At(index));
	    if (!clus) {
	      AliWarning(Form("Could not find cluster %d!", index));
	      continue;
	    }
	    AliVCluster *clus2 =  static_cast<AliVCluster*>(fCaloClusters2->At(index2));
	    if (!clus2) {
	      AliWarning(Form("Could not find cluster %d!", index2));
	      continue;
	    }
	    TLorentzVector part, part2;
	    clus->GetMomentum(part, const_cast<Double_t*>(fVertex));
	    clus2->GetMomentum(part2, const_cast<Double_t*>(fVertex));
	    
	    d1 -= part.Pt();
	    d2 -= part2.Pt();
	    break;
	  }
	}
      }
    }
  }

  if (d1 < 0)
    d1 = 0;

  if (d2 < 0)
    d2 = 0;

  if (jet1->Pt() > 0)
    d1 /= jet1->Pt();
  else
    d1 = -1;

  if (jet2->Pt() > 0)
    d2 /= jet2->Pt();
  else
    d2 = -1;
}

//________________________________________________________________________
void AliJetResponseMaker::SetMatchingLevel(AliEmcalJet *jet1, AliEmcalJet *jet2, MatchingType matching) 
{
  Double_t d1 = -1;
  Double_t d2 = -1;

  switch (matching) {
  case kGeometrical:
    GetGeometricalMatchingLevel(jet1,jet2,d1);
    d2 = d1;
    break;
  case kMCLabel: // jet1 = detector level and jet2 = particle level!
    GetMCLabelMatchingLevel(jet1,jet2,d1,d2);
    break;
  case kSameCollections:
    GetSameCollectionsMatchingLevel(jet1,jet2,d1,d2);
    break;
  default:
    ;
  }

  if (d1 >= 0) {

    if (d1 < jet1->ClosestJetDistance()) {
      jet1->SetSecondClosestJet(jet1->ClosestJet(), jet1->ClosestJetDistance());
      jet1->SetClosestJet(jet2, d1);
    }
    else if (d1 < jet1->SecondClosestJetDistance()) {
      jet1->SetSecondClosestJet(jet2, d1);
    }
  }
  
  if (d2 >= 0) {
    
    if (d2 < jet2->ClosestJetDistance()) {
      jet2->SetSecondClosestJet(jet2->ClosestJet(), jet2->ClosestJetDistance());
      jet2->SetClosestJet(jet1, d2);
    }
    else if (d2 < jet2->SecondClosestJetDistance()) {
      jet2->SetSecondClosestJet(jet1, d2);
    }
  }
}

//________________________________________________________________________
Bool_t AliJetResponseMaker::FillHistograms()
{
  // Fill histograms.

  static Int_t indexes[9999] = {-1};

  const Int_t nJets2 = GetSortedArray(indexes, fJets2, fRho2Val);

  Int_t naccJets2 = 0;
  Int_t naccJets2Acceptance = 0;

  for (Int_t i = 0; i < nJets2; i++) {
    AliDebug(2,Form("Processing jet (2) %d", indexes[i]));

    AliEmcalJet* jet2 = static_cast<AliEmcalJet*>(fJets2->At(indexes[i]));

    if (!jet2) {
      AliError(Form("Could not receive jet2 %d", i));
      continue;
    }

    // Jet 2 cuts
    if (!AcceptJet(jet2))
      continue;
    if (!AcceptBiasJet2(jet2))
      continue;
    if (jet2->Eta() < fJet2MinEta || jet2->Eta() > fJet2MaxEta || jet2->Phi() < fJet2MinPhi || jet2->Phi() > fJet2MaxPhi)
      continue;
    if (jet2->MaxTrackPt() > fMaxTrackPt2 || jet2->MaxClusterPt() > fMaxClusterPt2)
      continue;

    if (naccJets2 < fNLeadingJets) {
      fHistLeadingJets2PtArea->Fill(jet2->Area(), jet2->Pt());
      if (!fRho2Name.IsNull())
	fHistLeadingJets2CorrPtArea->Fill(jet2->Area(), jet2->Pt() - fRho2Val * jet2->Area());
    }

    if (fDoJet2Histogram)
      FillJetHisto(jet2->Phi(), jet2->Eta(), jet2->Pt(), jet2->Area(), jet2->NEF(), jet2->MaxPartPt()/jet2->Pt(), jet2->Pt() - fRho2Val * jet2->Area(), jet2->MCPt(), 2);

    naccJets2++;

    // Verify also jet cuts 1 on jet 2
    if (AcceptBiasJet(jet2) &&
	(jet2->Eta() > fJetMinEta && jet2->Eta() < fJetMaxEta && jet2->Phi() > fJetMinPhi && jet2->Phi() < fJetMaxPhi)) {
      
      if (naccJets2Acceptance < fNLeadingJets) {
	fHistLeadingJets2PtAreaAcceptance->Fill(jet2->Area(), jet2->Pt());
	if (!fRho2Name.IsNull()) 
	  fHistLeadingJets2CorrPtAreaAcceptance->Fill(jet2->Area(), jet2->Pt() - fRho2Val * jet2->Area());
      }
      
      FillJetHisto(jet2->Phi(), jet2->Eta(), jet2->Pt(), jet2->Area(), jet2->NEF(), jet2->MaxPartPt()/jet2->Pt(), jet2->Pt() - fRho2Val * jet2->Area(), jet2->MCPt(), 3);
      naccJets2Acceptance++;
    }

    if (jet2->MatchedJet()) {

      if (AcceptJet(jet2->MatchedJet()) &&
	  jet2->MatchedJet()->Eta() > fJetMinEta && jet2->MatchedJet()->Eta() < fJetMaxEta &&
	  jet2->MatchedJet()->Phi() > fJetMinPhi && jet2->MatchedJet()->Phi() < fJetMaxPhi &&
	  AcceptBiasJet(jet2->MatchedJet()) &&
	  jet2->MatchedJet()->MaxTrackPt() < fMaxTrackPt && jet2->MatchedJet()->MaxClusterPt() < fMaxClusterPt) {

	Double_t d=-1, ce1=-1, ce2=-1;
	if (jet2->GetMatchingType() == kGeometrical) {
	  if (fAreCollections2MC && !fAreCollections1MC) // the other way around is not supported
	    GetMCLabelMatchingLevel(jet2->MatchedJet(), jet2, ce1, ce2);
	  else if (fAreCollections1MC == fAreCollections2MC)
	    GetSameCollectionsMatchingLevel(jet2->MatchedJet(), jet2, ce1, ce2);

	  d = jet2->ClosestJetDistance();
	}
	else if (jet2->GetMatchingType() == kMCLabel || jet2->GetMatchingType() == kSameCollections) {
	  GetGeometricalMatchingLevel(jet2->MatchedJet(), jet2, d);

	  ce1 = jet2->MatchedJet()->ClosestJetDistance();
	  ce2 = jet2->ClosestJetDistance();
	}

	Double_t corrpt1 = jet2->MatchedJet()->Pt() - fRhoVal * jet2->MatchedJet()->Area();
	Double_t corrpt2 = jet2->Pt() - fRho2Val * jet2->Area();

	FillMatchingHistos(jet2->MatchedJet()->Pt(), jet2->Pt(), jet2->MatchedJet()->Eta(), jet2->Eta(), jet2->MatchedJet()->Phi(), jet2->Phi(), 
			   jet2->MatchedJet()->Area(), jet2->Area(), d, ce1, ce2, corrpt1, corrpt2, jet2->MatchedJet()->MCPt(), 
			   jet2->MatchedJet()->NEF(), jet2->NEF(), jet2->MatchedJet()->MaxPartPt()/jet2->MatchedJet()->Pt(), jet2->MaxPartPt()/jet2->Pt());
      }
    }
  }

  const Int_t nJets1 = GetSortedArray(indexes, fJets, fRhoVal);

  Int_t naccJets1 = 0;

  for (Int_t i = 0; i < nJets1; i++) {
    AliDebug(2,Form("Processing jet (1) %d", indexes[i]));
    AliEmcalJet* jet1 = static_cast<AliEmcalJet*>(fJets->At(indexes[i]));

    if (!jet1) {
      AliError(Form("Could not receive jet %d", i));
      continue;
    }  

    if (!AcceptJet(jet1))
      continue;
    if (!AcceptBiasJet(jet1))
      continue;
    if (jet1->MaxTrackPt() > fMaxTrackPt || jet1->MaxClusterPt() > fMaxClusterPt)
      continue;
    if (jet1->Eta() < fJetMinEta || jet1->Eta() > fJetMaxEta || jet1->Phi() < fJetMinPhi || jet1->Phi() > fJetMaxPhi)
      continue;

    if (naccJets1 < fNLeadingJets){
      fHistLeadingJets1PtArea->Fill(jet1->Area(), jet1->Pt());
      if (!fRhoName.IsNull())
	fHistLeadingJets1CorrPtArea->Fill(jet1->Area(), jet1->Pt() - fRhoVal * jet1->Area());
    }

    FillJetHisto(jet1->Phi(), jet1->Eta(), jet1->Pt(), jet1->Area(), jet1->NEF(), jet1->MaxPartPt()/jet1->Pt(), jet1->Pt() - fRhoVal * jet1->Area(), jet1->MCPt(), 1);
    naccJets1++;
  }

  return kTRUE;
}
