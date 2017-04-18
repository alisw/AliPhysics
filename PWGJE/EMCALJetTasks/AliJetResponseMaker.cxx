/*************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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
//
// Emcal jet response matrix maker task.
//
// Author: Salvatore Aiola, Yale University, salvatore.aiola@cern.ch

#include "AliJetResponseMaker.h"

#include <TClonesArray.h>
#include <TH2F.h>
#include <THnSparse.h>

#include "AliTLorentzVector.h"
#include "AliAnalysisManager.h"
#include "AliVCluster.h"
#include "AliVTrack.h"
#include "AliEmcalJet.h"
#include "AliLog.h"
#include "AliRhoParameter.h"
#include "AliNamedArrayI.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"

ClassImp(AliJetResponseMaker)

//________________________________________________________________________
AliJetResponseMaker::AliJetResponseMaker() : 
  AliAnalysisTaskEmcalJet("AliJetResponseMaker", kTRUE),
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
  fZgAxis(0),
  fdRAxis(0),
  fPtgAxis(0),
  fDBCAxis(0),
  fFlavourZAxis(0),
  fFlavourPtAxis(0),
  fJetRelativeEPAngle(0),
  fIsJet1Rho(kFALSE),
  fIsJet2Rho(kFALSE),
  fHistRejectionReason1(0),
  fHistRejectionReason2(0),
  fHistJets1(0),
  fHistJets2(0),
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
  fZgAxis(0),
  fdRAxis(0),
  fPtgAxis(0),
  fDBCAxis(0),
  fFlavourZAxis(0),
  fFlavourPtAxis(0),
  fJetRelativeEPAngle(0),
  fIsJet1Rho(kFALSE),
  fIsJet2Rho(kFALSE),
  fHistRejectionReason1(0),
  fHistRejectionReason2(0),
  fHistJets1(0),
  fHistJets2(0),
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

  if (fIsJet1Rho) {
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

  if (fIsJet2Rho) {
    fHistJets2CorrPtArea = new TH2F("fHistJets2CorrPtArea", "fHistJets2CorrPtArea", 
        fNbins/2, 0, 1.5, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    fHistJets2CorrPtArea->GetXaxis()->SetTitle("area");
    fHistJets2CorrPtArea->GetYaxis()->SetTitle("p_{T,2}^{corr} (GeV/c)");
    fHistJets2CorrPtArea->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistJets2CorrPtArea);
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

  if (fIsJet1Rho || fIsJet2Rho) {
    if (!fIsJet1Rho) 
      fHistDeltaCorrPtvsJet1CorrPt = new TH2F("fHistDeltaCorrPtvsJet1CorrPt", "fHistDeltaCorrPtvsJet1CorrPt", 
          fNbins, fMinBinPt, fMaxBinPt, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    else
      fHistDeltaCorrPtvsJet1CorrPt = new TH2F("fHistDeltaCorrPtvsJet1CorrPt", "fHistDeltaCorrPtvsJet1CorrPt", 
          2*fNbins, -fMaxBinPt, fMaxBinPt, 2*fNbins, -fMaxBinPt, fMaxBinPt);

    fHistDeltaCorrPtvsJet1CorrPt->GetXaxis()->SetTitle("p_{T,1}^{corr}");  
    fHistDeltaCorrPtvsJet1CorrPt->GetYaxis()->SetTitle("#deltap_{T}^{corr} (GeV/c)");
    fHistDeltaCorrPtvsJet1CorrPt->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistDeltaCorrPtvsJet1CorrPt);

    if (!fIsJet2Rho) 
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

    if (!fIsJet1Rho) 
      fHistJet1CorrPtvsJet2CorrPt = new TH2F("fHistJet1CorrPtvsJet2CorrPt", "fHistJet1CorrPtvsJet2CorrPt", 
          fNbins, fMinBinPt, fMaxBinPt, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    else if (!fIsJet2Rho) 
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

  TString title[25]= {""};
  Int_t nbins[25]  = {0};
  Double_t min[25] = {0.};
  Double_t max[25] = {0.};
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

  if (fNEFAxis) {
    title[dim] = "NEF";
    nbins[dim] = 120;
    min[dim] = 0.0;
    max[dim] = 1.2;
    dim++;
  }

  if (fZAxis) {
    title[dim] = "Z";
    nbins[dim] = 120;
    min[dim] = 0.0;
    max[dim] = 1.2;
    dim++;
  }

  if (fFlavourZAxis) {
    title[dim] = "z_{flavour}";
    nbins[dim] = 100;
    min[dim] = 0.0;
    max[dim] = 2.0;
    dim++;
  }

  if (fFlavourPtAxis) {
    title[dim] = "p_{T}^{D}";
    nbins[dim] = fNbins/2;
    min[dim] = 0;
    max[dim] = 125;
    dim++;
  }

  if (fJetRelativeEPAngle) {
    title[dim] = "#theta_{jet}^{EP}";
    nbins[dim] = 3;
    min[dim] = 0;
    max[dim] = TMath::Pi()/2;
    dim++;
  }

  title[dim] = "p_{T,particle}^{leading} (GeV/c)";
  nbins[dim] = 120;
  min[dim] = 0;
  max[dim] = 120;
  dim++;

  Int_t dim1 = dim, dim2 = dim;

  if (fIsJet1Rho) {
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

  if (fIsJet2Rho) {
    title[dim2] = "p_{T}^{corr}";
    nbins[dim2] = fNbins*2;
    min[dim2] = -250;
    max[dim2] = 250;
    dim2++;
  }

  fHistJets2 = new THnSparseD("fHistJets2","fHistJets2",dim2,nbins,min,max);
  for (Int_t i = 0; i < dim2; i++) 
    fHistJets2->GetAxis(i)->SetTitle(title[i]);
  fOutput->Add(fHistJets2);

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

  title[dim] = "p_{T,particle,1}^{leading} (GeV/c)";
  nbins[dim] = 120;
  min[dim] = 0;
  max[dim] = 120;
  dim++;

  title[dim] = "p_{T,particle,2}^{leading} (GeV/c)";
  nbins[dim] = 120;
  min[dim] = 0;
  max[dim] = 120;
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
  if (fIsJet1Rho) {
    title[dim] = "p_{T,1}^{corr}";
    nbins[dim] = fNbins*2;
    min[dim] = -250;
    max[dim] = 250;
    dim++;
  }
  if (fIsJet2Rho) {
    title[dim] = "p_{T,2}^{corr}";
    nbins[dim] = fNbins*2;
    min[dim] = -250;
    max[dim] = 250;
    dim++;
  }
  if (fDeltaPtAxis && (fIsJet1Rho || fIsJet2Rho)) {
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
    nbins[dim] = 120;
    min[dim] = 0.0;
    max[dim] = 1.2;
    dim++;

    title[dim] = "NEF_{2}";
    nbins[dim] = 120;
    min[dim] = 0.0;
    max[dim] = 1.2;
    dim++;
  }

  if (fZAxis) {
    title[dim] = "Z_{1}";
    nbins[dim] = 120;
    min[dim] = 0.0;
    max[dim] = 1.2;
    dim++;

    title[dim] = "Z_{2}";
    nbins[dim] = 120;
    min[dim] = 0.0;
    max[dim] = 1.2;
    dim++;
  }

  if (fFlavourZAxis) {
    title[dim] = "z_{flavour,1}";
    nbins[dim] = 100;
    min[dim] = 0.0;
    max[dim] = 2.0;
    dim++;

    title[dim] = "z_{flavour,2}";
    nbins[dim] = 100;
    min[dim] = 0.0;
    max[dim] = 2.0;
    dim++;
  }

  if (fFlavourPtAxis) {
    title[dim] = "p_{T,1}^{D}";
    nbins[dim] = fNbins/2;
    min[dim] = 0;
    max[dim] = 125;
    dim++;

    title[dim] = "p_{T,2}^{D}";
    nbins[dim] = fNbins/2;
    min[dim] = 0;
    max[dim] = 125;
    dim++;
  }

  if (fZgAxis) {
    title[dim] = "Z_{g,1}";
    nbins[dim] = 20;
    min[dim] = 0.0;
    max[dim] = 1.0;
    dim++;
    title[dim] = "Z_{g,2}";
    nbins[dim] = 20;
    min[dim] = 0.0;
    max[dim] = 1.0;
    dim++;
  }

  if (fdRAxis) {
    title[dim] = "dR_{1}";
    nbins[dim] = 40;
    min[dim] = 0.0;
    max[dim] = 0.5;
    dim++;
    title[dim] = "dR_{2}";
    nbins[dim] = 40;
    min[dim] = 0.0;
    max[dim] = 0.5;
    dim++;
  }

  if (fPtgAxis) {
    title[dim] = "p_{T,g,1}";
    nbins[dim] = 16;
    min[dim] = 0.0;
    max[dim] = 160.0;
    dim++;
    title[dim] = "p_{T,g,2}";
    nbins[dim] = 16;
    min[dim] = 0.0;
    max[dim] = 160.0;
    dim++;
  }

  if (fDBCAxis) {
    title[dim] = "DBC_{1}";
    nbins[dim] = 20;
    min[dim] = 0.0;
    max[dim] = 20.0;
    dim++;
    title[dim] = "DBC_{2}";
    nbins[dim] = 20;
    min[dim] = 0.0;
    max[dim] = 20.0;
    dim++;
  }

  if (fJetRelativeEPAngle) {
    title[dim] = "#theta_{jet,1}^{EP}";
    nbins[dim] = 3;
    min[dim] = 0;
    max[dim] = TMath::Pi()/2;
    dim++;

    title[dim] = "#theta_{jet,2}^{EP}";
    nbins[dim] = 3;
    min[dim] = 0;
    max[dim] = TMath::Pi()/2;
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

  AliJetContainer *jets1 = static_cast<AliJetContainer*>(fJetCollArray.At(0));
  AliJetContainer *jets2 = static_cast<AliJetContainer*>(fJetCollArray.At(1));

  if (!jets1 || !jets2) return;

  if (jets1->GetRhoName().IsNull()) fIsJet1Rho = kFALSE;
  else fIsJet1Rho = kTRUE;
  if (jets2->GetRhoName().IsNull()) fIsJet2Rho = kFALSE;
  else fIsJet2Rho = kTRUE;

  fHistRejectionReason1 = new TH2F("fHistRejectionReason1", "fHistRejectionReason1", 32, 0, 32, 100, 0, 250);
  fHistRejectionReason1->GetXaxis()->SetTitle("Rejection reason");
  fHistRejectionReason1->GetYaxis()->SetTitle("p_{T,jet} (GeV/c)");
  fHistRejectionReason1->GetZaxis()->SetTitle("counts");
  SetRejectionReasonLabels(fHistRejectionReason1->GetXaxis());
  fOutput->Add(fHistRejectionReason1);

  fHistRejectionReason2 = new TH2F("fHistRejectionReason2", "fHistRejectionReason2", 32, 0, 32, 100, 0, 250);
  fHistRejectionReason2->GetXaxis()->SetTitle("Rejection reason");
  fHistRejectionReason2->GetYaxis()->SetTitle("p_{T,jet} (GeV/c)");
  fHistRejectionReason2->GetZaxis()->SetTitle("counts");
  SetRejectionReasonLabels(fHistRejectionReason2->GetXaxis());
  fOutput->Add(fHistRejectionReason2);

  if (fHistoType==0)
    AllocateTH2();
  else 
    AllocateTHnSparse();

  PostData(1, fOutput); // Post data for ALL output slots > 0 here, to get at least an empty histogram
}

//________________________________________________________________________
void AliJetResponseMaker::FillJetHisto(AliEmcalJet* jet, Int_t Set)
{
  AliJetContainer* jets = GetJetContainer(Set-1);

  AliTLorentzVector leadPart;
  jets->GetLeadingHadronMomentum(leadPart, jet);
  Double_t zleading = GetParallelFraction(leadPart.Vect(), jet);
  if (zleading == 1 || (zleading > 1 && zleading - 1 < 1e-3)) zleading = 0.999; // so that it will contribute to the bin 0.9-1 rather than 1-1.1

  Double_t corrpt = jet->Pt() - jets->GetRhoVal() * jet->Area();
  Double_t zflavour = 0;
  Double_t ptflavour = 0;
  AliVParticle* hftrack = jet->GetFlavourTrack();
  if (hftrack) {
    zflavour = GetParallelFraction(hftrack, jet);
    ptflavour = hftrack->Pt();

    if (zflavour == 1 || (zflavour > 1 && zflavour - 1 < 1e-3)) zflavour = 0.999; // so that it will contribute to the bin 0.9-1 rather than 1-1.1
  }
  Double_t jetRelativeEPAngle = GetRelativeEPAngle(jet->Phi(), fEPV0);

  if (fHistoType==1) {
    THnSparse *histo = 0;
    if (Set==1) {
      histo = fHistJets1;
    }
    else if (Set==2) {
      histo = fHistJets2;
    }

    if (!histo) return;

    Double_t contents[20]={0};

    for (Int_t i = 0; i < histo->GetNdimensions(); i++) {
      TString title(histo->GetAxis(i)->GetTitle());
      if (title=="#phi")
        contents[i] = jet->Phi();
      else if (title=="#eta")
        contents[i] = jet->Eta();
      else if (title=="p_{T}")
        contents[i] = jet->Pt();
      else if (title=="A_{jet}")
        contents[i] = jet->Area();
      else if (title=="NEF")
        contents[i] = jet->NEF();
      else if (title=="Z")
        contents[i] = zleading;
      else if (title=="p_{T}^{corr}")
        contents[i] = corrpt;
      else if (title=="p_{T}^{MC}")
        contents[i] = jet->MCPt();
      else if (title=="p_{T,particle}^{leading} (GeV/c)")
        contents[i] = leadPart.Pt();
      else if (title=="z_{flavour}")
        contents[i] = zflavour;
      else if (title=="p_{T}^{D}")
        contents[i] = ptflavour;
      else if (title=="#theta_{jet}^{EP}")
        contents[i] = jetRelativeEPAngle;
      else 
        AliWarning(Form("Unable to fill dimension %s!",title.Data()));
    }

    histo->Fill(contents);
  }
  else {
    if (Set == 1) {
      fHistJets1PtArea->Fill(jet->Area(), jet->Pt());
      fHistJets1PhiEta->Fill(jet->Eta(), jet->Phi());
      fHistJets1ZvsPt->Fill(zleading, jet->Pt());
      fHistJets1NEFvsPt->Fill(jet->NEF(), jet->Pt());
      fHistJets1CEFvsCEFPt->Fill(1-jet->NEF(), (1-jet->NEF())*jet->Pt());
      if (fIsJet1Rho)
        fHistJets1CorrPtArea->Fill(jet->Area(), corrpt);
    }
    else if (Set == 2) {
      fHistJets2PtArea->Fill(jet->Area(), jet->Pt());
      fHistJets2PhiEta->Fill(jet->Eta(), jet->Phi());
      fHistJets2ZvsPt->Fill(zleading, jet->Pt());
      fHistJets2NEFvsPt->Fill(jet->NEF(), jet->Pt());
      fHistJets2CEFvsCEFPt->Fill(1-jet->NEF(), (1-jet->NEF())*jet->Pt());
      if (fIsJet2Rho)
        fHistJets2CorrPtArea->Fill(jet->Area(), corrpt);
    }
  }
}

//________________________________________________________________________
void AliJetResponseMaker::FillMatchingHistos(AliEmcalJet* jet1, AliEmcalJet* jet2, Double_t d, Double_t CE1, Double_t CE2)
{
  AliJetContainer* jets1 = GetJetContainer(0);
  AliJetContainer* jets2 = GetJetContainer(1);

  /*
  Printf("Detector level jet");
  jet1->Print();

  jet1->PrintConstituents(jets1->GetParticleContainer() != 0 ? jets1->GetParticleContainer()->GetArray() : 0,
                          jets1->GetClusterContainer() != 0 ? jets1->GetClusterContainer()->GetArray() : 0);

  Printf("Matched with particle level jet");
  jet2->Print();
  jet2->PrintConstituents(jets2->GetParticleContainer() != 0 ? jets2->GetParticleContainer()->GetArray() : 0,
                          jets2->GetClusterContainer() != 0 ? jets2->GetClusterContainer()->GetArray() : 0);
   */

  AliTLorentzVector leadPart1;
  jets1->GetLeadingHadronMomentum(leadPart1, jet1);
  Double_t zleading1 = GetParallelFraction(leadPart1.Vect(), jet1);
  if (zleading1 == 1 || (zleading1 > 1 && zleading1 - 1 < 1e-3)) zleading1 = 0.999; // so that it will contribute to the bin 0.9-1 rather than 1-1.1

  Double_t corrpt1 = jet1->Pt() - jets1->GetRhoVal() * jet1->Area();
  Double_t zflavour1 = 0;
  Double_t ptflavour1 = 0;
  AliVParticle* hftrack1 = jet1->GetFlavourTrack();
  if (hftrack1) {
    zflavour1 = GetParallelFraction(hftrack1, jet1);
    ptflavour1 = hftrack1->Pt();

    if (zflavour1 == 1 || (zflavour1 > 1 && zflavour1 - 1 < 1e-3)) zflavour1 = 0.999; // so that it will contribute to the bin 0.9-1 rather than 1-1.1
  }
  Double_t jetRelativeEPAngle1 = GetRelativeEPAngle(jet1->Phi(), fEPV0);


  AliTLorentzVector leadPart2;
  jets2->GetLeadingHadronMomentum(leadPart2, jet2);
  Double_t zleading2 = GetParallelFraction(leadPart2.Vect(), jet2);
  if (zleading2 == 1 || (zleading2 > 1 && zleading2 - 1 < 1e-3)) zleading2 = 0.999; // so that it will contribute to the bin 0.9-1 rather than 1-1.1

  Double_t corrpt2 = jet2->Pt() - jets2->GetRhoVal() * jet2->Area();
  Double_t zflavour2 = 0;
  Double_t ptflavour2 = 0;
  AliVParticle* hftrack2 = jet2->GetFlavourTrack();
  if (hftrack2) {
    zflavour2 = GetParallelFraction(hftrack2, jet2);
    ptflavour2 = hftrack2->Pt();

    if (zflavour2 == 1 || (zflavour2 > 1 && zflavour2 - 1 < 1e-3)) zflavour2 = 0.999; // so that it will contribute to the bin 0.9-1 rather than 1-1.1
  }
  Double_t jetRelativeEPAngle2 = GetRelativeEPAngle(jet2->Phi(), fEPV0);

  if (fHistoType==1) {
    Double_t contents[20]={0};

    for (Int_t i = 0; i < fHistMatching->GetNdimensions(); i++) {
      TString title(fHistMatching->GetAxis(i)->GetTitle());
      if (title=="p_{T,1}")
        contents[i] = jet1->Pt();
      else if (title=="p_{T,2}")
        contents[i] = jet2->Pt();
      else if (title=="A_{jet,1}")
        contents[i] = jet1->Area();
      else if (title=="A_{jet,2}")
        contents[i] = jet2->Area();
      else if (title=="distance")
        contents[i] = d;
      else if (title=="CE1")
        contents[i] = CE1;
      else if (title=="CE2")
        contents[i] = CE2;
      else if (title=="#deltaA_{jet}")
        contents[i] = jet1->Area()-jet2->Area();
      else if (title=="#deltap_{T}")
        contents[i] = jet1->Pt()-jet2->Pt();
      else if (title=="#delta#eta")
        contents[i] = jet1->Eta()-jet2->Eta();
      else if (title=="#delta#phi")
        contents[i] = jet1->Phi()-jet2->Phi();
      else if (title=="p_{T,1}^{corr}")
        contents[i] = corrpt1;
      else if (title=="p_{T,2}^{corr}")
        contents[i] = corrpt2;
      else if (title=="#deltap_{T}^{corr}")
        contents[i] = corrpt1-corrpt2;
      else if (title=="p_{T,1}^{MC}")
        contents[i] = jet1->MCPt();
      else if (title=="#deltap_{T}^{MC}")
        contents[i] = jet1->MCPt()-jet2->Pt();
      else if (title=="NEF_{1}")
        contents[i] = jet1->NEF();
      else if (title=="NEF_{2}")
        contents[i] = jet2->NEF();
      else if (title=="Z_{1}")
        contents[i] = zleading1;
      else if (title=="Z_{2}")
        contents[i] = zleading2;
      else if (title=="p_{T,particle,1}^{leading} (GeV/c)")
        contents[i] = leadPart1.Pt();
      else if (title=="p_{T,particle,2}^{leading} (GeV/c)")
        contents[i] = leadPart2.Pt();
      else if (title=="z_{flavour,1}")
        contents[i] = zflavour1;
      else if (title=="z_{flavour,2}")
        contents[i] = zflavour2;
      else if (title=="p_{T,1}^{D}")
        contents[i] = ptflavour1;
      else if (title=="p_{T,2}^{D}")
        contents[i] = ptflavour2;
      else if (title=="Z_{g,1}")
        contents[i] = jet1->GetShapeProperties()->GetSoftDropZg();
      else if (title=="Z_{g,2}")
        contents[i] = jet2->GetShapeProperties()->GetSoftDropZg();
      else if (title=="dR_{1}")
        contents[i] = jet1->GetShapeProperties()->GetSoftDropdR();
      else if (title=="dR_{2}")
        contents[i] = jet2->GetShapeProperties()->GetSoftDropdR();
      else if (title=="p_{T,g,1}")
        contents[i] = ( jet1->GetShapeProperties()->GetSoftDropPtfrac() )*( jet1->Pt() );
      else if (title=="p_{T,g,2}")
        contents[i] = ( jet2->GetShapeProperties()->GetSoftDropPtfrac() )*( jet2->Pt() );
      else if (title=="DBC_{1}")
        contents[i] = ( jet1->GetShapeProperties()->GetSoftDropDropCount() );
      else if (title=="DBC_{2}")
        contents[i] = ( jet2->GetShapeProperties()->GetSoftDropDropCount() );
      else if (title=="#theta_{jet,1}^{EP}")
        contents[i] = jetRelativeEPAngle1;
      else if (title=="#theta_{jet,2}^{EP}")
        contents[i] = jetRelativeEPAngle2;
      else 
        AliWarning(Form("Unable to fill dimension %s!",title.Data()));
    }

    fHistMatching->Fill(contents);
  }
  else {
    fHistCommonEnergy1vsJet1Pt->Fill(CE1, jet1->Pt());
    fHistCommonEnergy2vsJet2Pt->Fill(CE2, jet2->Pt());

    fHistDistancevsJet1Pt->Fill(d, jet1->Pt());
    fHistDistancevsJet2Pt->Fill(d, jet2->Pt());

    fHistDistancevsCommonEnergy1->Fill(d, CE1);
    fHistDistancevsCommonEnergy2->Fill(d, CE2);
    fHistCommonEnergy1vsCommonEnergy2->Fill(CE1, CE2);

    fHistDeltaEtaDeltaPhi->Fill(jet1->Eta()-jet2->Eta(),jet1->Phi()-jet2->Phi());

    fHistJet2PtOverJet1PtvsJet2Pt->Fill(jet2->Pt(), jet2->Pt() / jet1->Pt());
    fHistJet1PtOverJet2PtvsJet1Pt->Fill(jet1->Pt(), jet1->Pt() / jet2->Pt());

    Double_t dpt = jet1->Pt() - jet2->Pt();
    fHistDeltaPtvsJet1Pt->Fill(jet1->Pt(), dpt);
    fHistDeltaPtvsJet2Pt->Fill(jet2->Pt(), dpt);
    fHistDeltaPtOverJet1PtvsJet1Pt->Fill(jet1->Pt(), dpt/jet1->Pt());
    fHistDeltaPtOverJet2PtvsJet2Pt->Fill(jet2->Pt(), dpt/jet2->Pt());

    fHistDeltaPtvsDistance->Fill(d, dpt);
    fHistDeltaPtvsCommonEnergy1->Fill(CE1, dpt);
    fHistDeltaPtvsCommonEnergy2->Fill(CE2, dpt);

    fHistDeltaPtvsArea1->Fill(jet1->Area(), dpt);
    fHistDeltaPtvsArea2->Fill(jet2->Area(), dpt);

    Double_t darea = jet1->Area() - jet2->Area();
    fHistDeltaPtvsDeltaArea->Fill(darea, dpt);

    fHistJet1PtvsJet2Pt->Fill(jet1->Pt(), jet2->Pt());

    if (fIsJet1Rho || fIsJet2Rho) {
      Double_t dcorrpt = corrpt1 - corrpt2;
      fHistDeltaCorrPtvsJet1CorrPt->Fill(corrpt1, dcorrpt);
      fHistDeltaCorrPtvsJet2CorrPt->Fill(corrpt2, dcorrpt);
      fHistDeltaCorrPtOverJet1CorrPtvsJet1CorrPt->Fill(corrpt1, dcorrpt/corrpt1);
      fHistDeltaCorrPtOverJet2CorrPtvsJet2CorrPt->Fill(corrpt2, dcorrpt/corrpt2);
      fHistDeltaCorrPtvsDistance->Fill(d, dcorrpt);
      fHistDeltaCorrPtvsCommonEnergy1->Fill(CE1, dcorrpt);
      fHistDeltaCorrPtvsCommonEnergy2->Fill(CE2, dcorrpt);
      fHistDeltaCorrPtvsArea1->Fill(jet1->Area(), dcorrpt);
      fHistDeltaCorrPtvsArea2->Fill(jet2->Area(), dcorrpt);
      fHistDeltaCorrPtvsDeltaArea->Fill(darea, dcorrpt);
      fHistJet1CorrPtvsJet2CorrPt->Fill(corrpt1, corrpt2);
    }

    if (fIsEmbedded) {
      Double_t dmcpt = jet1->MCPt() - jet2->Pt();
      fHistDeltaMCPtvsJet1MCPt->Fill(jet1->MCPt(), dmcpt);
      fHistDeltaMCPtvsJet2Pt->Fill(jet2->MCPt(), dmcpt);
      fHistDeltaMCPtOverJet1MCPtvsJet1MCPt->Fill(jet1->MCPt(), dmcpt/jet1->MCPt());
      fHistDeltaMCPtOverJet2PtvsJet2Pt->Fill(jet2->Pt(), dmcpt/jet2->Pt());
      fHistDeltaMCPtvsDistance->Fill(d, dmcpt);
      fHistDeltaMCPtvsCommonEnergy1->Fill(CE1, dmcpt);
      fHistDeltaMCPtvsCommonEnergy2->Fill(CE2, dmcpt);
      fHistDeltaMCPtvsArea1->Fill(jet1->Area(), dmcpt);
      fHistDeltaMCPtvsArea2->Fill(jet2->Area(), dmcpt);
      fHistDeltaMCPtvsDeltaArea->Fill(darea, dmcpt);
      fHistJet1MCPtvsJet2Pt->Fill(jet1->MCPt(), jet2->Pt());
    }
  }
}

//________________________________________________________________________
void AliJetResponseMaker::ExecOnce()
{
  // Execute once.

  AliAnalysisTaskEmcalJet::ExecOnce();

  AliJetContainer *jets1 = static_cast<AliJetContainer*>(fJetCollArray.At(0));
  AliJetContainer *jets2 = static_cast<AliJetContainer*>(fJetCollArray.At(1));

  if (!jets1 || !jets1->GetArray() || !jets2 || !jets2->GetArray()) return;

  if (fMatching == kMCLabel && (!jets2->GetIsParticleLevel() || jets1->GetIsParticleLevel())) {
    if (jets1->GetIsParticleLevel() == jets2->GetIsParticleLevel()) {
      AliWarning("Changing matching type from MC label to same collection...");
      fMatching = kSameCollections;
    }
    else {
      AliWarning("Changing matching type from MC label to geometrical...");
      fMatching = kGeometrical;
    }
  }
  else if (fMatching == kSameCollections && jets1->GetIsParticleLevel() != jets2->GetIsParticleLevel()) {
    if (jets2->GetIsParticleLevel() && !jets1->GetIsParticleLevel()) {
      AliWarning("Changing matching type from same collection to MC label...");
      fMatching = kMCLabel;
    }
    else {
      AliWarning("Changing matching type from same collection to geometrical...");
      fMatching = kGeometrical;
    }
  }
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
  AliJetContainer *jets1 = static_cast<AliJetContainer*>(fJetCollArray.At(0));
  AliJetContainer *jets2 = static_cast<AliJetContainer*>(fJetCollArray.At(1));

  if (!jets1 || !jets1->GetArray() || !jets2 || !jets2->GetArray()) return kFALSE;

  DoJetLoop();

  AliEmcalJet* jet1 = 0;

  jets1->ResetCurrentID();
  while ((jet1 = jets1->GetNextJet())) {

    AliEmcalJet *jet2 = jet1->ClosestJet();

    if (!jet2) continue;
    if (jet2->ClosestJet() != jet1) continue;
    if (jet1->ClosestJetDistance() > fMatchingPar1 || jet2->ClosestJetDistance() > fMatchingPar2) continue;

    // Matched jet found
    jet1->SetMatchedToClosest(fMatching);
    jet2->SetMatchedToClosest(fMatching);
    AliDebug(2,Form("Found matching: jet1 pt = %f, eta = %f, phi = %f, jet2 pt = %f, eta = %f, phi = %f", 
        jet1->Pt(), jet1->Eta(), jet1->Phi(),
        jet2->Pt(), jet2->Eta(), jet2->Phi()));
  }

  return kTRUE;
}

//________________________________________________________________________
void AliJetResponseMaker::DoJetLoop()
{
  // Do the jet loop.

  AliJetContainer *jets1 = static_cast<AliJetContainer*>(fJetCollArray.At(0));
  AliJetContainer *jets2 = static_cast<AliJetContainer*>(fJetCollArray.At(1));

  if (!jets1 || !jets1->GetArray() || !jets2 || !jets2->GetArray()) return;

  AliEmcalJet* jet1 = 0;
  AliEmcalJet* jet2 = 0;

  jets2->ResetCurrentID();
  while ((jet2 = jets2->GetNextJet())) jet2->ResetMatching();

  jets1->ResetCurrentID();
  while ((jet1 = jets1->GetNextJet())) {
    jet1->ResetMatching();

    if (jet1->MCPt() < fMinJetMCPt) continue;

    jets2->ResetCurrentID();
    while ((jet2 = jets2->GetNextJet())) {
      SetMatchingLevel(jet1, jet2, fMatching);
    } // jet2 loop
  } // jet1 loop
}

//________________________________________________________________________
void AliJetResponseMaker::GetGeometricalMatchingLevel(AliEmcalJet *jet1, AliEmcalJet *jet2, Double_t &d) const
{
  d = jet1->DeltaR(jet2);
}

//________________________________________________________________________
void AliJetResponseMaker::GetMCLabelMatchingLevel(AliEmcalJet *jet1, AliEmcalJet *jet2, Double_t &d1, Double_t &d2) const
{ 
  AliJetContainer *jets1 = static_cast<AliJetContainer*>(fJetCollArray.At(0));
  AliJetContainer *jets2 = static_cast<AliJetContainer*>(fJetCollArray.At(1));

  if (!jets1 || !jets1->GetArray() || !jets2 || !jets2->GetArray()) return;

  // tracks1 just serves as a proxy to ensure that tracks are in jets1
  AliParticleContainer *tracks1   = jets1->GetParticleContainer();
  // tracks2 is used to retrieve MC labels associated with tracks in the container
  // NOTE: For multiple containers, this would need to be generalized!
  AliParticleContainer *tracks2   = jets2->GetParticleContainer();

  // d1 and d2 represent the matching level: 0 = maximum level of matching, 1 = the two jets are completely unrelated
  d1 = jet1->Pt();
  d2 = jet2->Pt();
  Double_t totalPt1 = d1; // the total pt of the reconstructed jet will be cleaned from the background

  // remove completely tracks that are not MC particles (label == 0)
  if (tracks1 && tracks1->GetArray()) {
    for (Int_t iTrack = 0; iTrack < jet1->GetNumberOfTracks(); iTrack++) {
      AliVParticle *track = jet1->Track(iTrack);
      if (!track) {
        AliWarning(Form("Could not find track %d!", iTrack));
        continue;
      }

      Int_t MClabel = TMath::Abs(track->GetLabel());
      MClabel -= fMCLabelShift;
      if (MClabel != 0) continue;

      // this is not a MC particle; remove it completely
      AliDebug(3,Form("Track %d (pT = %f) is not a MC particle (MClabel = %d)!",iTrack,track->Pt(),MClabel));
      totalPt1 -= track->Pt();
      d1 -= track->Pt();
    }
  }

  // remove completely clusters that are not MC particles (label == 0)
  if (fUseCellsToMatch && fCaloCells) { 
    for (Int_t iClus = 0; iClus < jet1->GetNumberOfClusters(); iClus++) {
      AliVCluster *clus = jet1->Cluster(iClus);
      if (!clus) {
        AliWarning(Form("Could not find cluster %d!", iClus));
        continue;
      }
      AliTLorentzVector part;
      clus->GetMomentum(part, fVertex);

      for (Int_t iCell = 0; iCell < clus->GetNCells(); iCell++) {
        Int_t cellId = clus->GetCellAbsId(iCell);
        Double_t cellFrac = clus->GetCellAmplitudeFraction(iCell);

        Int_t MClabel = TMath::Abs(fCaloCells->GetCellMCLabel(cellId));
        MClabel -= fMCLabelShift;
        if (MClabel != 0) continue;

        // this is not a MC particle; remove it completely
        AliDebug(3,Form("Cell %d (frac = %f) is not a MC particle (MClabel = %d)!",iCell,cellFrac,MClabel));
        totalPt1 -= part.Pt() * cellFrac;
        d1 -= part.Pt() * cellFrac;
      }
    }
  }
  else {
    for (Int_t iClus = 0; iClus < jet1->GetNumberOfClusters(); iClus++) {
      AliVCluster *clus = jet1->Cluster(iClus);
      if (!clus) {
        AliWarning(Form("Could not find cluster %d!", iClus));
        continue;
      }
      TLorentzVector part;
      clus->GetMomentum(part, fVertex);

      Int_t MClabel = TMath::Abs(clus->GetLabel());
      MClabel -= fMCLabelShift;
      if (MClabel != 0) continue;

      // this is not a MC particle; remove it completely
      AliDebug(3,Form("Cluster %d (pT = %f) is not a MC particle (MClabel = %d)!",iClus,part.Pt(),MClabel));
      totalPt1 -= part.Pt();
      d1 -= part.Pt();
    }
  }

  for (Int_t iTrack2 = 0; iTrack2 < jet2->GetNumberOfTracks(); iTrack2++) {
    Bool_t track2Found = kFALSE;
    Int_t index2 = jet2->TrackAt(iTrack2);

    // now look for common particles in the track array
    for (Int_t iTrack = 0; iTrack < jet1->GetNumberOfTracks(); iTrack++) {
      AliVParticle *track = jet1->Track(iTrack);
      if (!track) {
        AliWarning(Form("Could not find track %d!", iTrack));
        continue;
      }
      Int_t MClabel = TMath::Abs(track->GetLabel());
      MClabel -= fMCLabelShift;	  
      if (MClabel <= 0) continue;

      Int_t index = -1;
      index = tracks2->GetIndexFromLabel(MClabel);
      if (index < 0) {
        AliDebug(2,Form("Track %d (pT = %f) does not have an associated MC particle (MClabel = %d)!",iTrack,track->Pt(),MClabel));
        continue;
      }

      if (index2 != index) continue;

      // found common particle
      d1 -= track->Pt();

      if (!track2Found) {
        AliVParticle *MCpart = jet2->Track(iTrack2);
        AliDebug(3,Form("Track %d (pT = %f, eta = %f, phi = %f) is associated with the MC particle %d (pT = %f, eta = %f, phi = %f)!",
            iTrack,track->Pt(),track->Eta(),track->Phi(),MClabel,MCpart->Pt(),MCpart->Eta(),MCpart->Phi()));
        d2 -= MCpart->Pt();
      }

      track2Found = kTRUE;
    }

    // now look for common particles in the cluster array
    if (fUseCellsToMatch && fCaloCells) { // if the cell colection is available, look for cells with a matched MC particle
      for (Int_t iClus = 0; iClus < jet1->GetNumberOfClusters(); iClus++) {
        AliVCluster *clus = jet1->Cluster(iClus);
        if (!clus) {
          AliWarning(Form("Could not find cluster %d!", iClus));
          continue;
        }
        AliTLorentzVector part;
        clus->GetMomentum(part, fVertex);

        for (Int_t iCell = 0; iCell < clus->GetNCells(); iCell++) {
          Int_t cellId = clus->GetCellAbsId(iCell);
          Double_t cellFrac = clus->GetCellAmplitudeFraction(iCell);

          Int_t MClabel = TMath::Abs(fCaloCells->GetCellMCLabel(cellId));
          MClabel -= fMCLabelShift;
          if (MClabel <= 0) continue;

          Int_t index1 = -1;
          index1 = tracks2->GetIndexFromLabel(MClabel);
          if (index1 < 0) {
            AliDebug(3,Form("Cell %d (frac = %f) does not have an associated MC particle (MClabel = %d)!",iCell,cellFrac,MClabel));
            continue;
          }

          if (index2 != index1) continue;

          // found common particle
          d1 -= part.Pt() * cellFrac;

          if (!track2Found) { // only if it is not already found among charged tracks (charged particles are most likely already found)
            AliVParticle *MCpart = jet2->Track(iTrack2);
            AliDebug(3,Form("Cell %d belonging to cluster %d (pT = %f, eta = %f, phi = %f) is associated with the MC particle %d (pT = %f, eta = %f, phi = %f)!",
                iCell,iClus,part.Pt(),part.Eta(),part.Phi_0_2pi(),MClabel,MCpart->Pt(),MCpart->Eta(),MCpart->Phi()));
            d2 -= MCpart->Pt() * cellFrac;
          }

          track2Found = kTRUE;
        }
      }
    }
    else { //otherwise look for the first contributor to the cluster, and if matched to a MC label remove it
      for (Int_t iClus = 0; iClus < jet1->GetNumberOfClusters(); iClus++) {
        AliVCluster *clus = jet1->Cluster(iClus);
        if (!clus) {
          AliWarning(Form("Could not find cluster %d!", iClus));
          continue;
        }
        AliTLorentzVector part;
        clus->GetMomentum(part, fVertex);

        Int_t MClabel = TMath::Abs(clus->GetLabel());
        MClabel -= fMCLabelShift;
        if (MClabel <= 0) continue;

        Int_t index = -1;
        index = tracks2->GetIndexFromLabel(MClabel);

        if (index < 0) {
          AliDebug(3,Form("Cluster %d (pT = %f) does not have an associated MC particle (MClabel = %d)!",iClus,part.Pt(),MClabel));
          continue;
        }

        if (index2 != index) continue;

        // found common particle
        d1 -= part.Pt();

        if (!track2Found) { // only if it is not already found among charged tracks (charged particles are most likely already found)
          AliVParticle *MCpart = jet2->Track(iTrack2);
          AliDebug(3,Form("Cluster %d (pT = %f, eta = %f, phi = %f) is associated with the MC particle %d (pT = %f, eta = %f, phi = %f)!",
              iClus,part.Pt(),part.Eta(),part.Phi_0_2pi(),MClabel,MCpart->Pt(),MCpart->Eta(),MCpart->Phi()));

          d2 -= MCpart->Pt();
        }

        track2Found = kTRUE;
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
  AliJetContainer *jets1 = static_cast<AliJetContainer*>(fJetCollArray.At(0));
  AliJetContainer *jets2 = static_cast<AliJetContainer*>(fJetCollArray.At(1));

  if (!jets1 || !jets1->GetArray() || !jets2 || !jets2->GetArray()) return;

  // All of the containers are simply used as proxies for whether tracks or clusters are in a jet
  AliParticleContainer *tracks1   = jets1->GetParticleContainer();
  AliClusterContainer  *clusters1 = jets1->GetClusterContainer();
  AliParticleContainer *tracks2   = jets2->GetParticleContainer();
  AliClusterContainer  *clusters2 = jets2->GetClusterContainer();

  // d1 and d2 represent the matching level: 0 = maximum level of matching, 1 = the two jets are completely unrelated
  d1 = jet1->Pt();
  d2 = jet2->Pt();

  if (tracks1 && tracks2) {

    for (Int_t iTrack2 = 0; iTrack2 < jet2->GetNumberOfTracks(); iTrack2++) {
      Int_t index2 = jet2->TrackAt(iTrack2);
      for (Int_t iTrack1 = 0; iTrack1 < jet1->GetNumberOfTracks(); iTrack1++) {
        Int_t index1 = jet1->TrackAt(iTrack1);
        if (index2 == index1) { // found common particle
          AliVParticle *part1 = jet1->Track(iTrack1);
          if (!part1) {
            AliWarning(Form("Could not find track %d!", index1));
            continue;
          }
          AliVParticle *part2 = jet2->Track(iTrack2);
          if (!part2) {
            AliWarning(Form("Could not find track %d!", index2));
            continue;
          }

          d1 -= part1->Pt();
          d2 -= part2->Pt();
          break;
        }
      }
    }

  }

  if (clusters1 && clusters2) {

    if (fUseCellsToMatch && fCaloCells) {
      // Note: this section of the code needs to be revised and tested heavily
      // While fixing some inconsistencies in AliAnalysisTaskEMCALClusterizeFast
      // some issues came up, e.g. memory leaks (fixed) and inconsistent use of
      // fCaloCells. In principle the two cluster collections may use cells
      // from different sources (embedding / non-embedding). This is not handled
      // correctly in the current version of this code.
      AliWarning("ATTENTION ATTENTION ATTENTION: this section of the AliJetResponseMaker code needs to be revised and tested before using it for physics!!!");
      const Int_t nClus1 = jet1->GetNumberOfClusters();

      Int_t ncells1[nClus1];
      UShort_t *cellsId1[nClus1];
      Double_t *cellsFrac1[nClus1];
      Double_t *cellsClusFrac1[nClus1];
      Int_t *sortedIndexes1[nClus1];
      Double_t ptClus1[nClus1];
      for (Int_t iClus1 = 0; iClus1 < nClus1; iClus1++) {
        Int_t index1 = jet1->ClusterAt(iClus1);
        AliVCluster *clus1 = clusters1->GetCluster(index1);
        if (!clus1) {
          AliWarning(Form("Could not find cluster %d!", index1));
          ncells1[iClus1] = 0;
          cellsId1[iClus1] = 0;
          cellsFrac1[iClus1] = 0;
          cellsClusFrac1[iClus1] = 0;
          sortedIndexes1[iClus1] = 0;
          ptClus1[iClus1] = 0;
          continue;
        }
        TLorentzVector part1;
        clus1->GetMomentum(part1, fVertex);

        ncells1[iClus1] = clus1->GetNCells();
        cellsId1[iClus1] = clus1->GetCellsAbsId();
        cellsFrac1[iClus1] = clus1->GetCellsAmplitudeFraction();
        cellsClusFrac1[iClus1] = new Double_t[ncells1[iClus1]];
        sortedIndexes1[iClus1] = new Int_t[ncells1[iClus1]];
        ptClus1[iClus1] = part1.Pt();

        for (Int_t iCell = 0; iCell < ncells1[iClus1]; iCell++) {
          cellsClusFrac1[iClus1][iCell] = fCaloCells->GetCellAmplitude(cellsId1[iClus1][iCell]) / clus1->E();
        }

        TMath::Sort(ncells1[iClus1], cellsId1[iClus1], sortedIndexes1[iClus1], kFALSE);
      }

      const Int_t nClus2 = jet2->GetNumberOfClusters();

      const Int_t maxNcells2 = 11520;
      Int_t sortedIndexes2[maxNcells2];
      for (Int_t iClus2 = 0; iClus2 < nClus2; iClus2++) {
        Int_t index2 = jet2->ClusterAt(iClus2);
        AliVCluster *clus2 =  clusters2->GetCluster(index2);
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
        Double_t *cellsClusFrac2 = new Double_t[ncells2];

        for (Int_t iCell = 0; iCell < ncells2; iCell++) {
          cellsClusFrac2[iCell] = fCaloCells->GetCellAmplitude(cellsId2[iCell]) / clus2->E();
        }

        TLorentzVector part2;
        clus2->GetMomentum(part2, fVertex);
        Double_t ptClus2 = part2.Pt();

        TMath::Sort(ncells2, cellsId2, sortedIndexes2, kFALSE);

        for (Int_t iClus1 = 0; iClus1 < nClus1; iClus1++) {
          if (sortedIndexes1[iClus1] == 0)
            continue;
          Int_t iCell1 = 0, iCell2 = 0;
          while (iCell1 < ncells1[iClus1] && iCell2 < ncells2) {
            if (cellsId1[iClus1][sortedIndexes1[iClus1][iCell1]] == cellsId2[sortedIndexes2[iCell2]]) { // found a common cell
              d1 -= cellsFrac1[iClus1][sortedIndexes1[iClus1][iCell1]] * cellsClusFrac1[iClus1][sortedIndexes1[iClus1][iCell1]] * ptClus1[iClus1];
              d2 -= cellsFrac2[sortedIndexes2[iCell2]] * cellsClusFrac2[sortedIndexes2[iCell2]] * ptClus2;
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
        delete[] cellsClusFrac2;
      }
      for (Int_t iClus1 = 0; iClus1 < nClus1; iClus1++) {
        delete[] cellsClusFrac1[iClus1];
        delete[] sortedIndexes1[iClus1];
      }
    }
    else {
      for (Int_t iClus2 = 0; iClus2 < jet2->GetNumberOfClusters(); iClus2++) {
        Int_t index2 = jet2->ClusterAt(iClus2);
        for (Int_t iClus1 = 0; iClus1 < jet1->GetNumberOfClusters(); iClus1++) {
          Int_t index1 = jet1->ClusterAt(iClus1);
          if (index2 == index1) { // found common particle
            AliVCluster *clus1 = jet1->Cluster(iClus1);
            if (!clus1) {
              AliWarning(Form("Could not find cluster %d!", index1));
              continue;
            }
            AliVCluster *clus2 =  jet2->Cluster(iClus2);
            if (!clus2) {
              AliWarning(Form("Could not find cluster %d!", index2));
              continue;
            }
            TLorentzVector part1, part2;
            clus1->GetMomentum(part1, fVertex);
            clus2->GetMomentum(part2, fVertex);

            d1 -= part1.Pt();
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

  AliJetContainer *jets1 = static_cast<AliJetContainer*>(fJetCollArray.At(0));
  AliJetContainer *jets2 = static_cast<AliJetContainer*>(fJetCollArray.At(1));

  if (!jets1 || !jets1->GetArray() || !jets2 || !jets2->GetArray()) return kFALSE;

  AliEmcalJet* jet1 = 0;  
  AliEmcalJet* jet2 = 0;

  jets2->ResetCurrentID();
  while ((jet2 = jets2->GetNextJet())) {

    AliDebug(2,Form("Processing jet (2) %d", jets2->GetCurrentID()));

    if (jet2->Pt() < jets2->GetJetPtCut()) continue;

    UInt_t rejectionReason = 0;
    if (jets2->AcceptJet(jet2, rejectionReason))
      FillJetHisto(jet2, 2);
    else
      fHistRejectionReason2->Fill(jets2->GetRejectionReasonBitPosition(rejectionReason), jet2->Pt());

    jet1 = jet2->MatchedJet();

    if (!jet1) continue;
    rejectionReason = 0;
    if (!jets1->AcceptJet(jet1, rejectionReason)) continue;
    if (jet1->MCPt() < fMinJetMCPt) continue;

    Double_t d=-1, ce1=-1, ce2=-1;
    if (jet2->GetMatchingType() == kGeometrical) {
      if (jets2->GetIsParticleLevel() && !jets1->GetIsParticleLevel()) // the other way around is not supported
        GetMCLabelMatchingLevel(jet1, jet2, ce1, ce2);
      else if (jets1->GetIsParticleLevel() == jets2->GetIsParticleLevel())
        GetSameCollectionsMatchingLevel(jet1, jet2, ce1, ce2);

      d = jet2->ClosestJetDistance();
    }
    else if (jet2->GetMatchingType() == kMCLabel || jet2->GetMatchingType() == kSameCollections) {
      GetGeometricalMatchingLevel(jet1, jet2, d);

      ce1 = jet1->ClosestJetDistance();
      ce2 = jet2->ClosestJetDistance();
    }

    FillMatchingHistos(jet1, jet2, d, ce1, ce2);
  }

  jets1->ResetCurrentID();
  while ((jet1 = jets1->GetNextJet())) {
    UInt_t rejectionReason = 0;
    if (!jets1->AcceptJet(jet1, rejectionReason)) {
      fHistRejectionReason1->Fill(jets1->GetRejectionReasonBitPosition(rejectionReason), jet1->Pt());
      continue;
    }
    if (jet1->MCPt() < fMinJetMCPt) continue;
    AliDebug(2,Form("Processing jet (1) %d", jets1->GetCurrentID()));

    FillJetHisto(jet1, 1);
  }
  return kTRUE;
}

/**
 * Function to calculate angle between jet and EP in the 1st quadrant (0,Pi/2).
 * Adapted from AliAnalysisTaskEmcalJetHadEPpid.
 *
 * @param jetAngle Phi angle of the jet (could be any particle)
 * @param epAngle Event plane angle
 *
 * @return Angle between jet and EP in the 1st quadrant (0,Pi/2)
 */
Double_t AliJetResponseMaker::GetRelativeEPAngle(Double_t jetAngle, Double_t epAngle) const
{
  Double_t dphi = (epAngle - jetAngle);

  // ran into trouble with a few dEP<-Pi so trying this...
  if( dphi<-1*TMath::Pi()  ){
    dphi = dphi + 1*TMath::Pi();
  } // this assumes we are doing full jets currently

  if( (dphi>0) && (dphi<1*TMath::Pi()/2)  ){
    // Do nothing! we are in quadrant 1
  }else if( (dphi>1*TMath::Pi()/2) && (dphi<1*TMath::Pi())  ){
    dphi = 1*TMath::Pi() - dphi;
  }else if( (dphi<0) && (dphi>-1*TMath::Pi()/2)  ){
    dphi = fabs(dphi);
  }else if( (dphi<-1*TMath::Pi()/2) && (dphi>-1*TMath::Pi())  ){
    dphi = dphi + 1*TMath::Pi();
  }

  // test
  if( dphi < 0 || dphi > TMath::Pi()/2  )
    AliWarning(Form("%s: dPHI not in range [0, 0.5*Pi]!", GetName()));

  return dphi;   // dphi in [0, Pi/2]
}
