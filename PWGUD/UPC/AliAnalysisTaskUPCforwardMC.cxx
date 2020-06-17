/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

// c++ headers
#include <iostream>
#include <fstream>
// #include <vector>
// #include <algorithm>


// root headers
#include <TMath.h>
#include "TH1I.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TH1.h"
#include "THnSparse.h"
#include <TFile.h>
#include <TF2.h>
#include <TF1.h>
#include <TRandom.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TTree.h>
#include <TGraph2D.h>
#include <TStopwatch.h>
#include <TMatrixDSym.h>
#include <TFitResult.h>
#include <TLatex.h>
#include "TClonesArray.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TParticle.h"
#include "TObjString.h"
#include "TList.h"
#include "TChain.h"


// aliroot headers
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliMuonTrackCuts.h"
#include "AliAODVertex.h"


// my headers
#include "AliAnalysisTaskUPCforwardMC.h"


// headers for MC
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliAODMCParticle.h"



class AliAnalysisTaskUPCforwardMC;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskUPCforwardMC) // classimp: necessary for root

//_____________________________________________________________________________
AliAnalysisTaskUPCforwardMC::AliAnalysisTaskUPCforwardMC()
    : AliAnalysisTaskSE(),
      fAOD(0),
      fOutputList(0),
      fMCEvent(0),
      fNumberMuonsH(0),
      fCounterH(0),
      fEtaMuonH(0),
      fRAbsMuonH(0),
      fInvariantMassDistributionH(0),
      fInvariantMassDistributionRapidityBinsH{ 0, 0, 0, 0, 0, 0},
      fEntriesAgainstRunNumberH(0),
      fEntriesAgainstRunNumberProperlyH(0),
      fInvariantMassDistributionCoherentH(0),
      fInvariantMassDistributionCoherentRapidityBinsH{ 0, 0, 0, 0, 0, 0},
      fInvariantMassDistributionIncoherentH(0),
      fInvariantMassDistributionIncoherentRapidityBinsH{ 0, 0, 0, 0, 0, 0},
      fDimuonPtDistributionH(0),
      fTemplatePtDistributionH(0),
      fTemplatePtDistributionRapidityH{ 0, 0, 0 },
      //_______________________________
      //
      // SIDEBANDS
      //
      fTemplatePtDistributionHLowerSide(0),
      fTemplatePtDistributionRapidityHLowerSide{ 0, 0, 0 },
      fTemplatePtDistributionHHigherSide(0),
      fTemplatePtDistributionRapidityHHigherSide{ 0, 0, 0 },
      //_______________________________
      fDcaAgainstInvariantMassH(0),
      fInvariantMassDistributionExtendedH(0),
      fInvariantMassDistributionCoherentExtendedH(0),
      fInvariantMassDistributionIncoherentExtendedH(0),
      fMuonTrackCuts(0x0),
      fRunNum(0),
      fTracklets(0),
      fLumiPerRun(0),
      fL0inputs(0),
      fL1inputs(0),
      fZem1Energy(0),
      fZem2Energy(0),
      fZNCEnergy(0),
      fZNAEnergy(0),
      fZPCEnergy(0),
      fZPAEnergy(0),
      fZNATime(0),
      fZNCTime(0),
      fV0ADecision(-10),
      fV0CDecision(-10),
      fADADecision(-10),
      fADCDecision(-10),
      fIR1Map(0),
      fIR2Map(0),
      fZNATDC{0, 0, 0, 0},
      fZNCTDC{0, 0, 0, 0},
      fZPATDC{0, 0, 0, 0},
      fZPCTDC{0, 0, 0, 0},
      fV0Hits{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fV0TotalNCells(0),
      // fVectorGoodRunNumbers(0),
      fAngularDistribOfPositiveMuonRestFrameJPsiH(0),
      fAngularDistribOfNegativeMuonRestFrameJPsiH(0),
      fMCpdgCodesH(0),
      fMCpdgCodesOnlyPrimaryH(0),
      fMCphiGeneratedTruthH(0),
      fMCetaGeneratedTruthH(0),
      fMCpseudorapidityGeneratedTruthH(0),
      fMCptGeneratedTruthH(0),
      fMCphiDimuonGeneratedTruthH(0),
      fMCetaDimuonGeneratedTruthH(0),
      fMCpseudorapidityDimuonGeneratedTruthH(0),
      fMCptDimuonGeneratedTruthSingleMuonsH(0),
      fMCptDimuonGeneratedTruthH(0),
      fMCinvariantMassDistrJPsiGeneratedTruthH(0),
      fMCinvariantMassDistrJPsiAfterEvtAndTrkSelectionTruthH(0),
      fMCthetaDistribOfPositiveMuonRestFrameJPsiGeneratedTruthH(0),
      fMCthetaDistribOfNegativeMuonRestFrameJPsiGeneratedTruthH(0),
      fBBFlag{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fBGFlag{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fBBAFlags(0),
      fBBCFlags(0),
      fBGAFlags(0),
      fBGCFlags(0),
      fBBFlagAD{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fBGFlagAD{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fBBAFlagsAD(0),
      fBBCFlagsAD(0),
      fBGAFlagsAD(0),
      fBGCFlagsAD(0),
      // fVectorCosThetaGenerated(0),
      // fVectorCosThetaReconstructed(0),
      fCosThetaGeneratedHelicityFrame(0),
      fCosThetaReconHelicityFrame(0),
      fPhiGeneratedHelicityFrame(0),
      fPhiReconHelicityFrame(0),
      fCounterUPCevent(0),
      fBinMigrationHelicityH(0),
      fBinMigrationForPhiHelicityH(0),
      fCheckHelicityRestFrameJPsiH(0),
      fThetaDistribOfPositiveMuonRestFrameJPsiRapidityBinH{ 0, 0, 0, 0, 0, 0, 0, 0},
      fMCthetaDistribOfPositiveMuonRestFrameJPsiGeneratedTruthRapidityBinH{ 0, 0, 0, 0, 0, 0, 0, 0},
      fCosThetaHelicityFrameJPsiH(0),
      fCosThetaCollinsSoperFrameJPsiH(0),
      fPhiHelicityFrameJPsiH(0),
      fPhiCollinsSoperFrameJPsiH(0),
      fMCCosThetaHelicityFrameJPsiH(0),
      fMCCosThetaCollinsSoperFrameJPsiH(0),
      fMCPhiHelicityFrameJPsiH(0),
      fMCPhiCollinsSoperFrameJPsiH(0),
      fCosThetaHelicityFrameJPsiRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0},
      fCosThetaCollinsSoperFrameJPsiRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0},
      fPhiHelicityFrameJPsiRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0},
      fPhiCollinsSoperFrameJPsiRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0},
      fMCCosThetaHelicityFrameJPsiRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0},
      fMCCosThetaCollinsSoperFrameJPsiRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0},
      fMCPhiHelicityFrameJPsiRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0},
      fMCPhiCollinsSoperFrameJPsiRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0},
      fCosThetaHelicityFrameJPsiTenRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      fCosThetaCollinsSoperFrameJPsiTenRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      fPhiHelicityFrameJPsiTenRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      fPhiCollinsSoperFrameJPsiTenRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      fMCCosThetaHelicityFrameJPsiTenRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      fMCCosThetaCollinsSoperFrameJPsiTenRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      fMCPhiHelicityFrameJPsiTenRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      fMCPhiCollinsSoperFrameJPsiTenRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      fInvariantMassDistributionBinsOfCosThetaAndPhiHelicityFrameH(0),
      fMCInvariantMassDistributionBinsOfCosThetaAndPhiHelicityFrameH(0),
      fCosThetaAndPhiHelicityFrameInclusivePeopleBinningH(0),
      fMCCosThetaAndPhiHelicityFrameInclusivePeopleBinningH(0),
      fInvariantMassDistributionForSignalExtractionHelicityFrameH(0),
      fMCInvariantMassDistributionForSignalExtractionHelicityFrameH(0),
      fCosThetaHeFrameForSignalExH(0),
      fPhiHeFrameForSignalExH(0),
      fMCCosThetaHeFrameForSignalExH(0),
      fMCPhiHeFrameForSignalExH(0),
      fEfficiencyPerRunH(0),
      fMCEfficiencyPerRunH(0),
      fEfficiencyPerRunRapidityH{0,0,0,0,0,0},
      fMCEfficiencyPerRunRapidityH{0,0,0,0,0,0},
      fEtaAndPhi(0),
      fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fCosThetaAndPhiHelicityFrameMyBinningH(0),
      fMCCosThetaAndPhiHelicityFrameMyBinningH(0),
      fCosThetaAndPhiCsFrameMyBinningH(0),
      fMCCosThetaAndPhiCsFrameMyBinningH(0),
      fCosThetaAndPhiHelicityFrameMyBinningTriggerH{0,0,0,0,0,0,0},
      fMCCosThetaAndPhiHelicityFrameMyBinningTriggerH{0,0,0,0,0,0,0},
      fCosThetaAndPhiCsFrameMyBinningTriggerH{0,0,0,0,0,0,0},
      fMCCosThetaAndPhiCsFrameMyBinningTriggerH{0,0,0,0,0,0,0},
      fCosThetaAndPhiHelicityFrameMyBinningReweightingH(0),
      fMCCosThetaAndPhiHelicityFrameMyBinningReweightingH(0),
      fCosThetaAndPhiCsFrameMyBinningReweightingH(0),
      fMCCosThetaAndPhiCsFrameMyBinningReweightingH(0),
      fMCCosThetaAndPhiHelicityFrameReweightingH(0),
      fMCCosThetaAndPhiCsFrameReweightingH(0),
      fCosThetaHelicityFrameMyBinningH(0),
      fMCCosThetaHelicityFrameMyBinningH(0),
      fCosThetaHelicityFrameMyBinningSmallH(0),
      fMCCosThetaHelicityFrameMyBinningSmallH(0),
      fCosThetaHelicityFrameMySeventeenBinningH(0),
      fMCCosThetaHelicityFrameMySeventeenBinningH(0),
      fCosThetaHelicityFrameTwentyfiveBinsH(0),
      fMCCosThetaHelicityFrameTwentyfiveBinsH(0),
      fPhiHelicityFrameTwentyfiveBinsH(0),
      fMCPhiHelicityFrameTwentyfiveBinsH(0),
      fTildePhiHelicityFrameTwentyfiveBinsH(0),
      fMCTildePhiHelicityFrameTwentyfiveBinsH(0),
      fCosThetaHelicityFrameTwentyfiveBinsTriggerH{0,0,0,0,0,0,0},
      fMCCosThetaHelicityFrameTwentyfiveBinsTriggerH{0,0,0,0,0,0,0},
      fPhiHelicityFrameTwentyfiveBinsTriggerH{0,0,0,0,0,0,0},
      fMCPhiHelicityFrameTwentyfiveBinsTriggerH{0,0,0,0,0,0,0},
      fTildePhiHelicityFrameTwentyfiveBinsTriggerH{0,0,0,0,0,0,0},
      fMCTildePhiHelicityFrameTwentyfiveBinsTriggerH{0,0,0,0,0,0,0},
      fMCCosThetaHeVsCsH(0),
      fMCCosThetaHeVsCsFlatH(0),
      fCosThetaCsFrameTwentyfiveBinsH(0),
      fMCCosThetaCsFrameTwentyfiveBinsH(0),
      fPhiCsFrameTwentyfiveBinsH(0),
      fMCPhiCsFrameTwentyfiveBinsH(0),
      fTildePhiCsFrameTwentyfiveBinsH(0),
      fMCTildePhiCsFrameTwentyfiveBinsH(0),
      fCosThetaCsFrameTwentyfiveBinsTriggerH{0,0,0,0,0,0,0},
      fMCCosThetaCsFrameTwentyfiveBinsTriggerH{0,0,0,0,0,0,0},
      fPhiCsFrameTwentyfiveBinsTriggerH{0,0,0,0,0,0,0},
      fMCPhiCsFrameTwentyfiveBinsTriggerH{0,0,0,0,0,0,0},
      fTildePhiCsFrameTwentyfiveBinsTriggerH{0,0,0,0,0,0,0},
      fMCTildePhiCsFrameTwentyfiveBinsTriggerH{0,0,0,0,0,0,0},
      fPhiHelicityFrameMyBinningH(0),
      fMCPhiHelicityFrameMyBinningH(0),
      fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameMyVariableBinningH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                                                0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                                                0, 0, 0, 0, 0, 0 },
      fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameMyVariableBinningH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fInvariantMassDistributionOnlyPhiHeFrameTwentyfiveBinsH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                               0, 0, 0, 0, 0 },
      fInvariantMassDistributionOnlyCosThetaHeFrameTwentyfiveBinsH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                    0, 0, 0, 0, 0 },
      fInvariantMassDistributionOnlyTildePhiHeFrameTwentyfiveBinsH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                    0, 0, 0, 0, 0 },
      fInvariantMassDistributionOnlyPhiCsFrameTwentyfiveBinsH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                               0, 0, 0, 0, 0 },
      fInvariantMassDistributionOnlyCosThetaCsFrameTwentyfiveBinsH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                    0, 0, 0, 0, 0 },
      fInvariantMassDistributionOnlyTildePhiCsFrameTwentyfiveBinsH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                    0, 0, 0, 0, 0 },
      fInvariantMassDistributionForSignalExtractionHelicityFrameMyBinningH(0),
      fInvariantMassDistributionForSignalExtractionCsFrameMyBinningH(0)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}

//_____________________________________________________________________________
AliAnalysisTaskUPCforwardMC::AliAnalysisTaskUPCforwardMC( const char* name )
    : AliAnalysisTaskSE(name),
      fAOD(0),
      fOutputList(0),
      fMCEvent(0),
      fNumberMuonsH(0),
      fCounterH(0),
      fEtaMuonH(0),
      fRAbsMuonH(0),
      fInvariantMassDistributionH(0),
      fInvariantMassDistributionRapidityBinsH{ 0, 0, 0, 0, 0, 0},
      fEntriesAgainstRunNumberH(0),
      fEntriesAgainstRunNumberProperlyH(0),
      fInvariantMassDistributionCoherentH(0),
      fInvariantMassDistributionCoherentRapidityBinsH{ 0, 0, 0, 0, 0, 0},
      fInvariantMassDistributionIncoherentH(0),
      fInvariantMassDistributionIncoherentRapidityBinsH{ 0, 0, 0, 0, 0, 0},
      fDimuonPtDistributionH(0),
      fTemplatePtDistributionH(0),
      fTemplatePtDistributionRapidityH{ 0, 0, 0 },
      //_______________________________
      //
      // SIDEBANDS
      //
      fTemplatePtDistributionHLowerSide(0),
      fTemplatePtDistributionRapidityHLowerSide{ 0, 0, 0 },
      fTemplatePtDistributionHHigherSide(0),
      fTemplatePtDistributionRapidityHHigherSide{ 0, 0, 0 },
      //_______________________________
      fDcaAgainstInvariantMassH(0),
      fInvariantMassDistributionExtendedH(0),
      fInvariantMassDistributionCoherentExtendedH(0),
      fInvariantMassDistributionIncoherentExtendedH(0),
      fMuonTrackCuts(0x0),
      fRunNum(0),
      fTracklets(0),
      fLumiPerRun(0),
      fL0inputs(0),
      fL1inputs(0),
      fZem1Energy(0),
      fZem2Energy(0),
      fZNCEnergy(0),
      fZNAEnergy(0),
      fZPCEnergy(0),
      fZPAEnergy(0),
      fZNATime(0),
      fZNCTime(0),
      fV0ADecision(-10),
      fV0CDecision(-10),
      fADADecision(-10),
      fADCDecision(-10),
      fIR1Map(0),
      fIR2Map(0),
      fZNATDC{0, 0, 0, 0},
      fZNCTDC{0, 0, 0, 0},
      fZPATDC{0, 0, 0, 0},
      fZPCTDC{0, 0, 0, 0},
      fV0Hits{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fV0TotalNCells(0),
      // fVectorGoodRunNumbers(0),
      fAngularDistribOfPositiveMuonRestFrameJPsiH(0),
      fAngularDistribOfNegativeMuonRestFrameJPsiH(0),
      fMCpdgCodesH(0),
      fMCpdgCodesOnlyPrimaryH(0),
      fMCphiGeneratedTruthH(0),
      fMCetaGeneratedTruthH(0),
      fMCpseudorapidityGeneratedTruthH(0),
      fMCptGeneratedTruthH(0),
      fMCphiDimuonGeneratedTruthH(0),
      fMCetaDimuonGeneratedTruthH(0),
      fMCpseudorapidityDimuonGeneratedTruthH(0),
      fMCptDimuonGeneratedTruthSingleMuonsH(0),
      fMCptDimuonGeneratedTruthH(0),
      fMCinvariantMassDistrJPsiGeneratedTruthH(0),
      fMCinvariantMassDistrJPsiAfterEvtAndTrkSelectionTruthH(0),
      fMCthetaDistribOfPositiveMuonRestFrameJPsiGeneratedTruthH(0),
      fMCthetaDistribOfNegativeMuonRestFrameJPsiGeneratedTruthH(0),
      fBBFlag{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fBGFlag{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fBBAFlags(0),
      fBBCFlags(0),
      fBGAFlags(0),
      fBGCFlags(0),
      fBBFlagAD{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fBGFlagAD{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fBBAFlagsAD(0),
      fBBCFlagsAD(0),
      fBGAFlagsAD(0),
      fBGCFlagsAD(0),
      // fVectorCosThetaGenerated(0),
      // fVectorCosThetaReconstructed(0),
      fCosThetaGeneratedHelicityFrame(0),
      fCosThetaReconHelicityFrame(0),
      fPhiGeneratedHelicityFrame(0),
      fPhiReconHelicityFrame(0),
      fCounterUPCevent(0),
      fBinMigrationHelicityH(0),
      fBinMigrationForPhiHelicityH(0),
      fCheckHelicityRestFrameJPsiH(0),
      fThetaDistribOfPositiveMuonRestFrameJPsiRapidityBinH{ 0, 0, 0, 0, 0, 0, 0, 0},
      fMCthetaDistribOfPositiveMuonRestFrameJPsiGeneratedTruthRapidityBinH{ 0, 0, 0, 0, 0, 0, 0, 0},
      fCosThetaHelicityFrameJPsiH(0),
      fCosThetaCollinsSoperFrameJPsiH(0),
      fPhiHelicityFrameJPsiH(0),
      fPhiCollinsSoperFrameJPsiH(0),
      fMCCosThetaHelicityFrameJPsiH(0),
      fMCCosThetaCollinsSoperFrameJPsiH(0),
      fMCPhiHelicityFrameJPsiH(0),
      fMCPhiCollinsSoperFrameJPsiH(0),
      fCosThetaHelicityFrameJPsiRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0},
      fCosThetaCollinsSoperFrameJPsiRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0},
      fPhiHelicityFrameJPsiRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0},
      fPhiCollinsSoperFrameJPsiRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0},
      fMCCosThetaHelicityFrameJPsiRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0},
      fMCCosThetaCollinsSoperFrameJPsiRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0},
      fMCPhiHelicityFrameJPsiRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0},
      fMCPhiCollinsSoperFrameJPsiRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0},
      fCosThetaHelicityFrameJPsiTenRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      fCosThetaCollinsSoperFrameJPsiTenRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      fPhiHelicityFrameJPsiTenRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      fPhiCollinsSoperFrameJPsiTenRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      fMCCosThetaHelicityFrameJPsiTenRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      fMCCosThetaCollinsSoperFrameJPsiTenRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      fMCPhiHelicityFrameJPsiTenRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      fMCPhiCollinsSoperFrameJPsiTenRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      fInvariantMassDistributionBinsOfCosThetaAndPhiHelicityFrameH(0),
      fMCInvariantMassDistributionBinsOfCosThetaAndPhiHelicityFrameH(0),
      fCosThetaAndPhiHelicityFrameInclusivePeopleBinningH(0),
      fMCCosThetaAndPhiHelicityFrameInclusivePeopleBinningH(0),
      fInvariantMassDistributionForSignalExtractionHelicityFrameH(0),
      fMCInvariantMassDistributionForSignalExtractionHelicityFrameH(0),
      fCosThetaHeFrameForSignalExH(0),
      fPhiHeFrameForSignalExH(0),
      fMCCosThetaHeFrameForSignalExH(0),
      fMCPhiHeFrameForSignalExH(0),
      fEfficiencyPerRunH(0),
      fMCEfficiencyPerRunH(0),
      fEfficiencyPerRunRapidityH{0,0,0,0,0,0},
      fMCEfficiencyPerRunRapidityH{0,0,0,0,0,0},
      fEtaAndPhi(0),
      fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fCosThetaAndPhiHelicityFrameMyBinningH(0),
      fMCCosThetaAndPhiHelicityFrameMyBinningH(0),
      fCosThetaAndPhiCsFrameMyBinningH(0),
      fMCCosThetaAndPhiCsFrameMyBinningH(0),
      fCosThetaAndPhiHelicityFrameMyBinningTriggerH{0,0,0,0,0,0,0},
      fMCCosThetaAndPhiHelicityFrameMyBinningTriggerH{0,0,0,0,0,0,0},
      fCosThetaAndPhiCsFrameMyBinningTriggerH{0,0,0,0,0,0,0},
      fMCCosThetaAndPhiCsFrameMyBinningTriggerH{0,0,0,0,0,0,0},
      fCosThetaAndPhiHelicityFrameMyBinningReweightingH(0),
      fMCCosThetaAndPhiHelicityFrameMyBinningReweightingH(0),
      fCosThetaAndPhiCsFrameMyBinningReweightingH(0),
      fMCCosThetaAndPhiCsFrameMyBinningReweightingH(0),
      fMCCosThetaAndPhiHelicityFrameReweightingH(0),
      fMCCosThetaAndPhiCsFrameReweightingH(0),
      fCosThetaHelicityFrameMyBinningH(0),
      fMCCosThetaHelicityFrameMyBinningH(0),
      fCosThetaHelicityFrameMyBinningSmallH(0),
      fMCCosThetaHelicityFrameMyBinningSmallH(0),
      fCosThetaHelicityFrameMySeventeenBinningH(0),
      fMCCosThetaHelicityFrameMySeventeenBinningH(0),
      fCosThetaHelicityFrameTwentyfiveBinsH(0),
      fMCCosThetaHelicityFrameTwentyfiveBinsH(0),
      fPhiHelicityFrameTwentyfiveBinsH(0),
      fMCPhiHelicityFrameTwentyfiveBinsH(0),
      fTildePhiHelicityFrameTwentyfiveBinsH(0),
      fMCTildePhiHelicityFrameTwentyfiveBinsH(0),
      fCosThetaHelicityFrameTwentyfiveBinsTriggerH{0,0,0,0,0,0,0},
      fMCCosThetaHelicityFrameTwentyfiveBinsTriggerH{0,0,0,0,0,0,0},
      fPhiHelicityFrameTwentyfiveBinsTriggerH{0,0,0,0,0,0,0},
      fMCPhiHelicityFrameTwentyfiveBinsTriggerH{0,0,0,0,0,0,0},
      fTildePhiHelicityFrameTwentyfiveBinsTriggerH{0,0,0,0,0,0,0},
      fMCTildePhiHelicityFrameTwentyfiveBinsTriggerH{0,0,0,0,0,0,0},
      fMCCosThetaHeVsCsH(0),
      fMCCosThetaHeVsCsFlatH(0),
      fCosThetaCsFrameTwentyfiveBinsH(0),
      fMCCosThetaCsFrameTwentyfiveBinsH(0),
      fPhiCsFrameTwentyfiveBinsH(0),
      fMCPhiCsFrameTwentyfiveBinsH(0),
      fTildePhiCsFrameTwentyfiveBinsH(0),
      fMCTildePhiCsFrameTwentyfiveBinsH(0),
      fCosThetaCsFrameTwentyfiveBinsTriggerH{0,0,0,0,0,0,0},
      fMCCosThetaCsFrameTwentyfiveBinsTriggerH{0,0,0,0,0,0,0},
      fPhiCsFrameTwentyfiveBinsTriggerH{0,0,0,0,0,0,0},
      fMCPhiCsFrameTwentyfiveBinsTriggerH{0,0,0,0,0,0,0},
      fTildePhiCsFrameTwentyfiveBinsTriggerH{0,0,0,0,0,0,0},
      fMCTildePhiCsFrameTwentyfiveBinsTriggerH{0,0,0,0,0,0,0},
      fPhiHelicityFrameMyBinningH(0),
      fMCPhiHelicityFrameMyBinningH(0),
      fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameMyVariableBinningH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                                                0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                                                0, 0, 0, 0, 0, 0 },
      fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameMyVariableBinningH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fInvariantMassDistributionOnlyPhiHeFrameTwentyfiveBinsH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                               0, 0, 0, 0, 0 },
      fInvariantMassDistributionOnlyCosThetaHeFrameTwentyfiveBinsH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                    0, 0, 0, 0, 0 },
      fInvariantMassDistributionOnlyTildePhiHeFrameTwentyfiveBinsH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                    0, 0, 0, 0, 0 },
      fInvariantMassDistributionOnlyPhiCsFrameTwentyfiveBinsH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                               0, 0, 0, 0, 0 },
      fInvariantMassDistributionOnlyCosThetaCsFrameTwentyfiveBinsH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                    0, 0, 0, 0, 0 },
      fInvariantMassDistributionOnlyTildePhiCsFrameTwentyfiveBinsH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                    0, 0, 0, 0, 0 },
      fInvariantMassDistributionForSignalExtractionHelicityFrameMyBinningH(0),
      fInvariantMassDistributionForSignalExtractionCsFrameMyBinningH(0)
{
    // FillGoodRunVector(fVectorGoodRunNumbers);
    for( Int_t iRun = 0; iRun < 60000; iRun++) {
      fCounterGeneratedLevel[iRun]   = 0;
      // fDeadZoneEtaVsPhiPerRunH[iRun] = 0x0;
    }
    // fDeadZoneEtaVsPhiPerRunH[60000] = 0x0;
    for( Int_t iRun = 0; iRun < 364; iRun++) {
      fDeadZoneEtaVsPhiPerRunH[iRun] = 0x0;
    }
    fEtaAndPhi = new Double_t[2];

    // constructor
    DefineInput(0, TChain::Class());    // define the input of the analysis: in this case we take a 'chain' of events
                                        // this chain is created by the analysis manager, so no need to worry about it,
                                        // it does its work automatically
    DefineOutput(1, TList::Class());    // define the ouptut of the analysis: in this case it's a list of histograms
                                        // you can add more output objects by calling DefineOutput(2, classname::Class())
                                        // if you add more output objects, make sure to call PostData for all of them, and to
                                        // make changes to your AddTask macro!
}
//_____________________________________________________________________________
AliAnalysisTaskUPCforwardMC::~AliAnalysisTaskUPCforwardMC()
{
    // destructor
    if(fOutputList)    {delete fOutputList;}     	// at the end of your task, it is deleted
    if(fMuonTrackCuts) {delete fMuonTrackCuts;}   // from memory by calling this function
}
//_____________________________________________________________________________
// void AliAnalysisTaskUPCforwardMC::FillGoodRunVector(std::vector<Int_t> &fVectorGoodRunNumbers)
// {
//   fVectorGoodRunNumbers.clear();
//   Int_t listOfGoodRunNumbersLHC18q[] = { 295585, 295586, 295587, 295588, 295589, 295612,
//                                          295615, 295665, 295666, 295667, 295668, 295671,
//                                          295673, 295675, 295676, 295677, 295714, 295716,
//                                          295717, 295718, 295719, 295723, 295725, 295753,
//                                          295754, 295755, 295758, 295759, 295762, 295763,
//                                          295786, 295788, 295791, 295816, 295818, 295819,
//                                          295822, 295825, 295826, 295829, 295831, 295854,
//                                          295855, 295856, 295859, 295860, 295861, 295863,
//                                          295881, 295908, 295909, 295910, 295913, 295936,
//                                          295937, 295941, 295942, 295943, 295945, 295947,
//                                          296061, 296062, 296063, 296065, 296066, 296068,
//                                          296123, 296128, 296132, 296133, 296134, 296135,
//                                          296142, 296143, 296191, 296192, 296194, 296195,
//                                          296196, 296197, 296198, 296241, 296242, 296243,
//                                          296244, 296246, 296247, 296269, 296270, 296273,
//                                          296279, 296280, 296303, 296304, 296307, 296309,
//                                          296312, 296376, 296377, 296378, 296379, 296380,
//                                          296381, 296383, 296414, 296419, 296420, 296423,
//                                          296424, 296433, 296472, 296509, 296510, 296511,
//                                          296514, 296516, 296547, 296548, 296549, 296550,
//                                          296551, 296552, 296553, 296615, 296616, 296618,
//                                          296619, 296622, 296623 };
//   Int_t listOfGoodRunNumbersLHC18r[] = { 296690, 296691, 296694, 296749, 296750, 296781,
//                                          296784, 296785, 296786, 296787, 296791, 296793,
//                                          296794, 296799, 296836, 296838, 296839, 296848,
//                                          296849, 296850, 296851, 296852, 296890, 296894,
//                                          296899, 296900, 296903, 296930, 296931, 296932,
//                                          296934, 296935, 296938, 296941, 296966, 296967,
//                                          296968, 296969, 296971, 296975, 296976, 296977,
//                                          296979, 297029, 297031, 297035, 297085, 297117,
//                                          297118, 297119, 297123, 297124, 297128, 297129,
//                                          297132, 297133, 297193, 297194, 297196, 297218,
//                                          297219, 297221, 297222, 297278, 297310, 297312,
//                                          297315, 297317, 297363, 297366, 297367, 297372,
//                                          297379, 297380, 297405, 297408, 297413, 297414,
//                                          297415, 297441, 297442, 297446, 297450, 297451,
//                                          297452, 297479, 297481, 297483, 297512, 297537,
//                                          297540, 297541, 297542, 297544, 297558, 297588,
//                                          297590, 297595, 297623, 297624 };
//   Int_t sizeOfLHC18q = 0;
//   Int_t sizeOfLHC18r = 0;
//   for ( Int_t GoodRunNumberLHC18q : listOfGoodRunNumbersLHC18q ) {
//         fVectorGoodRunNumbers.push_back(GoodRunNumberLHC18q);
//         sizeOfLHC18q++;
//   }
//   for ( Int_t GoodRunNumberLHC18r : listOfGoodRunNumbersLHC18r ) {
//         fVectorGoodRunNumbers.push_back(GoodRunNumberLHC18r);
//         sizeOfLHC18r++;
//   }
// }
//_____________________________________________________________________________
void AliAnalysisTaskUPCforwardMC::UserCreateOutputObjects()
{
  // create output objects
  //
  // this function is called ONCE at the start of your analysis (RUNTIME)
  // here you ceate the histograms that you want to use
  //
  // the histograms are in this case added to a tlist, this list is in the end saved
  // to an output file
  //

  //muon track cuts
  fMuonTrackCuts = new AliMuonTrackCuts("StdMuonCuts", "StdMuonCuts");
  fMuonTrackCuts->SetFilterMask(    AliMuonTrackCuts::kMuEta     |
                                    AliMuonTrackCuts::kMuThetaAbs|
                                    AliMuonTrackCuts::kMuPdca    |
                                    AliMuonTrackCuts::kMuMatchLpt   );
  fMuonTrackCuts->SetAllowDefaultParams(kTRUE);
  fMuonTrackCuts->Print("mask");


  fOutputList = new TList();          // this is a list which will contain all
                                      // of your histograms at the end of the
                                      // analysis, the contents of this list
                                      // are written to the output file

  fOutputList->SetOwner(kTRUE);       // memory management: the list is owner
                                      // of all objects it contains and will
                                      // delete them if requested

  //_______________________________
  // - Adding histograms
  fNumberMuonsH = new TH1F("fNumberMuonsH", "fNumberMuonsH", 12, -0.5, 11.5);
  fOutputList->Add(fNumberMuonsH);    // don't forget to add it to the list!

  fCounterH = new TH1F("fCounterH", "fCounterH", 24, -0.5, 23.5);
  fOutputList->Add(fCounterH);

  fEtaMuonH = new TH1F("fEtaMuonH", "fEtaMuonH", 90, -2, -5);
  fOutputList->Add(fEtaMuonH);

  fRAbsMuonH = new TH1F("fRAbsMuonH", "fRAbsMuonH", 100, 0, 100);
  fOutputList->Add(fRAbsMuonH);

  fInvariantMassDistributionH = new TH1F("fInvariantMassDistributionH", "fInvariantMassDistributionH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionH);

  for( Int_t iRapidityBin = 0; iRapidityBin < 6; iRapidityBin++ ){
    fInvariantMassDistributionRapidityBinsH[iRapidityBin]
          = new TH1F( Form("fInvariantMassDistributionRapidityBinsH_%d", iRapidityBin),
                      Form("fInvariantMassDistributionRapidityBinsH_%d", iRapidityBin),
                      2000, 0, 20
                      );
    fOutputList->Add(fInvariantMassDistributionRapidityBinsH[iRapidityBin]);
  }

  fEntriesAgainstRunNumberH = new TH1F("fEntriesAgainstRunNumberH", "fEntriesAgainstRunNumberH", 10000, 290000, 300000);
  fOutputList->Add(fEntriesAgainstRunNumberH);

  /* - Trying to reproduce the histogram for the RunNumbers as they always
     - show it, properly labelled. Inspiration has come from the website:
     - https://root.cern.ch/doc/master/hlabels1_8C.html
     - Let us see if it works properly.
     -
   */
  fEntriesAgainstRunNumberProperlyH = new TH1F("fEntriesAgainstRunNumberProperlyH", "fEntriesAgainstRunNumberProperlyH", 10000, 290000, 300000);
  fEntriesAgainstRunNumberProperlyH->SetStats(0);
  fEntriesAgainstRunNumberProperlyH->SetFillColor(38);
  fEntriesAgainstRunNumberProperlyH->LabelsDeflate();
  fOutputList->Add(fEntriesAgainstRunNumberProperlyH);

  fInvariantMassDistributionCoherentH = new TH1F("fInvariantMassDistributionCoherentH", "fInvariantMassDistributionCoherentH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionCoherentH);

  for( Int_t iRapidityBin = 0; iRapidityBin < 6; iRapidityBin++ ){
    fInvariantMassDistributionCoherentRapidityBinsH[iRapidityBin]
          = new TH1F( Form("fInvariantMassDistributionCoherentRapidityBinsH_%d", iRapidityBin),
                      Form("fInvariantMassDistributionCoherentRapidityBinsH_%d", iRapidityBin),
                      2000, 0, 20
                      );
    fOutputList->Add(fInvariantMassDistributionCoherentRapidityBinsH[iRapidityBin]);
  }

  fInvariantMassDistributionIncoherentH = new TH1F("fInvariantMassDistributionIncoherentH", "fInvariantMassDistributionIncoherentH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionIncoherentH);

  for( Int_t iRapidityBin = 0; iRapidityBin < 6; iRapidityBin++ ){
    fInvariantMassDistributionIncoherentRapidityBinsH[iRapidityBin]
          = new TH1F( Form("fInvariantMassDistributionIncoherentRapidityBinsH_%d", iRapidityBin),
                      Form("fInvariantMassDistributionIncoherentRapidityBinsH_%d", iRapidityBin),
                      2000, 0, 20
                      );
    fOutputList->Add(fInvariantMassDistributionIncoherentRapidityBinsH[iRapidityBin]);
  }

  fDimuonPtDistributionH = new TH1F("fDimuonPtDistributionH", "fDimuonPtDistributionH", 4000, 0, 20);
  fOutputList->Add(fDimuonPtDistributionH);

  fTemplatePtDistributionH = new TH1F("fTemplatePtDistributionH", "fTemplatePtDistributionH", 4000, 0, 20);
  fOutputList->Add(fTemplatePtDistributionH);

  for( Int_t iRapidityBin = 0; iRapidityBin < 3; iRapidityBin++ ){
    fTemplatePtDistributionRapidityH[iRapidityBin] =
          new TH1F( Form( "fTemplatePtDistributionRapidityH_%d", iRapidityBin),
                    Form( "fTemplatePtDistributionRapidityH_%d", iRapidityBin),
                    4000, 0, 20
                    );
    fOutputList->Add(fTemplatePtDistributionRapidityH[iRapidityBin]);
  }

  //_______________________________
  /* -
   * - SIDEBANDS
   */
  fTemplatePtDistributionHLowerSide = new TH1F("fTemplatePtDistributionHLowerSide", "fTemplatePtDistributionHLowerSide", 4000, 0, 20);
  fOutputList->Add(fTemplatePtDistributionHLowerSide);

  for( Int_t iRapidityBin = 0; iRapidityBin < 3; iRapidityBin++ ){
    fTemplatePtDistributionRapidityHLowerSide[iRapidityBin] =
          new TH1F( Form( "fTemplatePtDistributionRapidityHLowerSide_%d", iRapidityBin),
                    Form( "fTemplatePtDistributionRapidityHLowerSide_%d", iRapidityBin),
                    4000, 0, 20
                    );
    fOutputList->Add(fTemplatePtDistributionRapidityHLowerSide[iRapidityBin]);
  }

  fTemplatePtDistributionHHigherSide = new TH1F("fTemplatePtDistributionHHigherSide", "fTemplatePtDistributionHHigherSide", 4000, 0, 20);
  fOutputList->Add(fTemplatePtDistributionHHigherSide);

  for( Int_t iRapidityBin = 0; iRapidityBin < 3; iRapidityBin++ ){
    fTemplatePtDistributionRapidityHHigherSide[iRapidityBin] =
          new TH1F( Form( "fTemplatePtDistributionRapidityHHigherSide_%d", iRapidityBin),
                    Form( "fTemplatePtDistributionRapidityHHigherSide_%d", iRapidityBin),
                    4000, 0, 20
                    );
    fOutputList->Add(fTemplatePtDistributionRapidityHHigherSide[iRapidityBin]);
  }
  //_______________________________

  fDcaAgainstInvariantMassH = new TH2F("fDcaAgainstInvariantMassH", "fDcaAgainstInvariantMassH", 4000, 0, 40, 2000, -100, 100);
  fOutputList->Add(fDcaAgainstInvariantMassH);

  /* - These histograms have an EXTENDED range (0,20)->(0,40)
     -
   */
  fInvariantMassDistributionExtendedH = new TH1F("fInvariantMassDistributionExtendedH", "fInvariantMassDistributionExtendedH", 4000, 0, 40);
  fOutputList->Add(fInvariantMassDistributionExtendedH);

  fInvariantMassDistributionCoherentExtendedH = new TH1F("fInvariantMassDistributionCoherentExtendedH", "fInvariantMassDistributionCoherentExtendedH", 4000, 0, 40);
  fOutputList->Add(fInvariantMassDistributionCoherentExtendedH);

  fInvariantMassDistributionIncoherentExtendedH = new TH1F("fInvariantMassDistributionIncoherentExtendedH", "fInvariantMassDistributionIncoherentExtendedH", 4000, 0, 40);
  fOutputList->Add(fInvariantMassDistributionIncoherentExtendedH);

  /* - Here starts the list of histograms needed for the analysis of the J/Psi's
     - polarization.
     -
   */
  fAngularDistribOfPositiveMuonRestFrameJPsiH = new TH1F("fAngularDistribOfPositiveMuonRestFrameJPsiH", "fAngularDistribOfPositiveMuonRestFrameJPsiH", 1000, -1., 1.);
  fOutputList->Add(fAngularDistribOfPositiveMuonRestFrameJPsiH);

  fAngularDistribOfNegativeMuonRestFrameJPsiH = new TH1F("fAngularDistribOfNegativeMuonRestFrameJPsiH", "fAngularDistribOfNegativeMuonRestFrameJPsiH", 1000, -1., 1.);
  fOutputList->Add(fAngularDistribOfNegativeMuonRestFrameJPsiH);

  fBinMigrationHelicityH = new TH2F("fBinMigrationHelicityH", "fBinMigrationHelicityH", 1000, -1., 1., 1000, -1., 1.);
  fOutputList->Add(fBinMigrationHelicityH);

  fBinMigrationForPhiHelicityH = new TH2F("fBinMigrationForPhiHelicityH", "fBinMigrationForPhiHelicityH", 1000, -3.14, 3.14, 1000, -3.14, 3.14);
  fOutputList->Add(fBinMigrationForPhiHelicityH);

  //_______________________________
  // - MC-only plots
  fMCpdgCodesH = new TH1F("fMCpdgCodesH", "fMCpdgCodesH", 3, 0, 3);
  fMCpdgCodesH ->SetStats(0);
  fMCpdgCodesH ->SetFillColor(38);
  fMCpdgCodesH->LabelsDeflate();
  fOutputList->Add(fMCpdgCodesH);

  fMCpdgCodesOnlyPrimaryH = new TH1F("fMCpdgCodesOnlyPrimaryH", "fMCpdgCodesOnlyPrimaryH", 3, 0, 3);
  fMCpdgCodesOnlyPrimaryH ->SetStats(0);
  fMCpdgCodesOnlyPrimaryH ->SetFillColor(38);
  fMCpdgCodesOnlyPrimaryH->LabelsDeflate();
  fOutputList->Add(fMCpdgCodesOnlyPrimaryH);

  fMCinvariantMassDistrJPsiGeneratedTruthH = new TH1F("fMCinvariantMassDistrJPsiGeneratedTruthH", "fMCinvariantMassDistrJPsiGeneratedTruthH", 2000000, 0, 20);
  fOutputList->Add(fMCinvariantMassDistrJPsiGeneratedTruthH);

  fMCinvariantMassDistrJPsiAfterEvtAndTrkSelectionTruthH = new TH1F("fMCinvariantMassDistrJPsiAfterEvtAndTrkSelectionTruthH", "fMCinvariantMassDistrJPsiAfterEvtAndTrkSelectionTruthH", 20000, 0, 20);
  fOutputList->Add(fMCinvariantMassDistrJPsiAfterEvtAndTrkSelectionTruthH);

  fMCphiGeneratedTruthH = new TH1F("fMCphiGeneratedTruthH", "fMCphiGeneratedTruthH", 2000, 0, 20);
  fOutputList->Add(fMCphiGeneratedTruthH);

  fMCetaGeneratedTruthH = new TH1F("fMCetaGeneratedTruthH", "fMCetaGeneratedTruthH", 2000, 0, 20);
  fOutputList->Add(fMCetaGeneratedTruthH);

  fMCpseudorapidityGeneratedTruthH = new TH1F("fMCpseudorapidityGeneratedTruthH", "fMCpseudorapidityGeneratedTruthH", 2000, 0, 20);
  fOutputList->Add(fMCpseudorapidityGeneratedTruthH);

  fMCptGeneratedTruthH = new TH1F("fMCptGeneratedTruthH", "fMCptGeneratedTruthH", 2000, 0, 20);
  fOutputList->Add(fMCptGeneratedTruthH);

  fMCphiDimuonGeneratedTruthH = new TH1F("fMCphiDimuonGeneratedTruthH", "fMCphiDimuonGeneratedTruthH", 2000, 0, 20);
  fOutputList->Add(fMCphiDimuonGeneratedTruthH);

  fMCetaDimuonGeneratedTruthH = new TH1F("fMCetaDimuonGeneratedTruthH", "fMCetaDimuonGeneratedTruthH", 2000, 0, 20);
  fOutputList->Add(fMCetaDimuonGeneratedTruthH);

  fMCpseudorapidityDimuonGeneratedTruthH = new TH1F("fMCpseudorapidityDimuonGeneratedTruthH", "fMCpseudorapidityDimuonGeneratedTruthH", 2000, 0, 20);
  fOutputList->Add(fMCpseudorapidityDimuonGeneratedTruthH);

  fMCptDimuonGeneratedTruthSingleMuonsH = new TH1F("fMCptDimuonGeneratedTruthSingleMuonsH", "fMCptDimuonGeneratedTruthSingleMuonsH", 4000, 0, 20);
  fOutputList->Add(fMCptDimuonGeneratedTruthSingleMuonsH);

  fMCptDimuonGeneratedTruthH = new TH1F("fMCptDimuonGeneratedTruthH", "fMCptDimuonGeneratedTruthH", 4000, 0, 20);
  fOutputList->Add(fMCptDimuonGeneratedTruthH);


  /* - Here starts the list of histograms needed for the analysis of the J/Psi's
     - polarization. GENERATED MC TRUTH.
     -
   */
  fMCthetaDistribOfPositiveMuonRestFrameJPsiGeneratedTruthH = new TH1F("fMCthetaDistribOfPositiveMuonRestFrameJPsiGeneratedTruthH", "fMCthetaDistribOfPositiveMuonRestFrameJPsiGeneratedTruthH", 1000, -1., 1.);
  fOutputList->Add(fMCthetaDistribOfPositiveMuonRestFrameJPsiGeneratedTruthH);

  fMCthetaDistribOfNegativeMuonRestFrameJPsiGeneratedTruthH = new TH1F("fMCthetaDistribOfNegativeMuonRestFrameJPsiGeneratedTruthH", "fMCthetaDistribOfNegativeMuonRestFrameJPsiGeneratedTruthH", 1000, -1., 1.);
  fOutputList->Add(fMCthetaDistribOfNegativeMuonRestFrameJPsiGeneratedTruthH);

  fCheckHelicityRestFrameJPsiH = new TH1F("fCheckHelicityRestFrameJPsiH", "fCheckHelicityRestFrameJPsiH", 100000, -50., 50.);
  fOutputList->Add(fCheckHelicityRestFrameJPsiH);

  for(Int_t iRapidityBin = 0; iRapidityBin < 8; iRapidityBin++ ){
    fThetaDistribOfPositiveMuonRestFrameJPsiRapidityBinH[iRapidityBin] = new TH1F(
                Form("fThetaDistribOfPositiveMuonRestFrameJPsiRapidityBinH_%d", iRapidityBin),
                Form("fThetaDistribOfPositiveMuonRestFrameJPsiRapidityBinH_%d", iRapidityBin),
                1000, -1., 1.
              );
    fOutputList->Add(fThetaDistribOfPositiveMuonRestFrameJPsiRapidityBinH[iRapidityBin]);

    fMCthetaDistribOfPositiveMuonRestFrameJPsiGeneratedTruthRapidityBinH[iRapidityBin] = new TH1F(
                Form("fMCthetaDistribOfPositiveMuonRestFrameJPsiGeneratedTruthRapidityBinH_%d", iRapidityBin),
                Form("fMCthetaDistribOfPositiveMuonRestFrameJPsiGeneratedTruthRapidityBinH_%d", iRapidityBin),
                1000, -1., 1.
              );
    fOutputList->Add(fMCthetaDistribOfPositiveMuonRestFrameJPsiGeneratedTruthRapidityBinH[iRapidityBin]);

  }

  fCosThetaHelicityFrameJPsiH = new TH1F("fCosThetaHelicityFrameJPsiH", "fCosThetaHelicityFrameJPsiH", 1000, -1., 1.);
  fOutputList->Add(fCosThetaHelicityFrameJPsiH);

  fCosThetaCollinsSoperFrameJPsiH = new TH1F("fCosThetaCollinsSoperFrameJPsiH", "fCosThetaCollinsSoperFrameJPsiH", 1000, -1., 1.);
  fOutputList->Add(fCosThetaCollinsSoperFrameJPsiH);

  fPhiHelicityFrameJPsiH = new TH1F("fPhiHelicityFrameJPsiH", "fPhiHelicityFrameJPsiH", 4000, -4., 4.);
  fOutputList->Add(fPhiHelicityFrameJPsiH);

  fPhiCollinsSoperFrameJPsiH = new TH1F("fPhiCollinsSoperFrameJPsiH", "fPhiCollinsSoperFrameJPsiH", 4000, -4., 4.);
  fOutputList->Add(fPhiCollinsSoperFrameJPsiH);

  fMCCosThetaHelicityFrameJPsiH = new TH1F("fMCCosThetaHelicityFrameJPsiH", "fMCCosThetaHelicityFrameJPsiH", 1000, -1., 1.);
  fOutputList->Add(fMCCosThetaHelicityFrameJPsiH);

  fMCCosThetaCollinsSoperFrameJPsiH = new TH1F("fMCCosThetaCollinsSoperFrameJPsiH", "fMCCosThetaCollinsSoperFrameJPsiH", 1000, -1., 1.);
  fOutputList->Add(fMCCosThetaCollinsSoperFrameJPsiH);

  fMCPhiHelicityFrameJPsiH = new TH1F("fMCPhiHelicityFrameJPsiH", "fMCPhiHelicityFrameJPsiH", 4000, -4., 4.);
  fOutputList->Add(fMCPhiHelicityFrameJPsiH);

  fMCPhiCollinsSoperFrameJPsiH = new TH1F("fMCPhiCollinsSoperFrameJPsiH", "fMCPhiCollinsSoperFrameJPsiH", 4000, -4., 4.);
  fOutputList->Add(fMCPhiCollinsSoperFrameJPsiH);

  /* - Rapidity-dependent analysis for the helicity of the J/Psi.
     - It is divided into RECONSTRUCTED and GENERATED for both
     - helicity frame and Collins-Soper.
     -
   */
  for(Int_t iRapidityBin = 0; iRapidityBin < 8; iRapidityBin++ ){
    fCosThetaHelicityFrameJPsiRapidityBinsH[iRapidityBin] = new TH1F(
                Form("fCosThetaHelicityFrameJPsiRapidityBinsH_%d", iRapidityBin),
                Form("fCosThetaHelicityFrameJPsiRapidityBinsH_%d", iRapidityBin),
                1000, -1., 1.
              );
    fOutputList->Add(fCosThetaHelicityFrameJPsiRapidityBinsH[iRapidityBin]);
  }

  for(Int_t iRapidityBin = 0; iRapidityBin < 8; iRapidityBin++ ){
    fCosThetaCollinsSoperFrameJPsiRapidityBinsH[iRapidityBin] = new TH1F(
                Form("fCosThetaCollinsSoperFrameJPsiRapidityBinsH_%d", iRapidityBin),
                Form("fCosThetaCollinsSoperFrameJPsiRapidityBinsH_%d", iRapidityBin),
                1000, -1., 1.
              );
    fOutputList->Add(fCosThetaCollinsSoperFrameJPsiRapidityBinsH[iRapidityBin]);
  }

  for(Int_t iRapidityBin = 0; iRapidityBin < 8; iRapidityBin++ ){
    fPhiHelicityFrameJPsiRapidityBinsH[iRapidityBin] = new TH1F(
                Form("fPhiHelicityFrameJPsiRapidityBinsH_%d", iRapidityBin),
                Form("fPhiHelicityFrameJPsiRapidityBinsH_%d", iRapidityBin),
                4000, -4., 4.
              );
    fOutputList->Add(fPhiHelicityFrameJPsiRapidityBinsH[iRapidityBin]);
  }

  for(Int_t iRapidityBin = 0; iRapidityBin < 8; iRapidityBin++ ){
    fPhiCollinsSoperFrameJPsiRapidityBinsH[iRapidityBin] = new TH1F(
                Form("fPhiCollinsSoperFrameJPsiRapidityBinsH_%d", iRapidityBin),
                Form("fPhiCollinsSoperFrameJPsiRapidityBinsH_%d", iRapidityBin),
                4000, -4., 4.
              );
    fOutputList->Add(fPhiCollinsSoperFrameJPsiRapidityBinsH[iRapidityBin]);
  }

  for(Int_t iRapidityBin = 0; iRapidityBin < 8; iRapidityBin++ ){
    fMCCosThetaHelicityFrameJPsiRapidityBinsH[iRapidityBin] = new TH1F(
                Form("fMCCosThetaHelicityFrameJPsiRapidityBinsH_%d", iRapidityBin),
                Form("fMCCosThetaHelicityFrameJPsiRapidityBinsH_%d", iRapidityBin),
                1000, -1., 1.
              );
    fOutputList->Add(fMCCosThetaHelicityFrameJPsiRapidityBinsH[iRapidityBin]);
  }

  for(Int_t iRapidityBin = 0; iRapidityBin < 8; iRapidityBin++ ){
    fMCCosThetaCollinsSoperFrameJPsiRapidityBinsH[iRapidityBin] = new TH1F(
                Form("fMCCosThetaCollinsSoperFrameJPsiRapidityBinsH_%d", iRapidityBin),
                Form("fMCCosThetaCollinsSoperFrameJPsiRapidityBinsH_%d", iRapidityBin),
                1000, -1., 1.
              );
    fOutputList->Add(fMCCosThetaCollinsSoperFrameJPsiRapidityBinsH[iRapidityBin]);
  }

  for(Int_t iRapidityBin = 0; iRapidityBin < 8; iRapidityBin++ ){
    fMCPhiHelicityFrameJPsiRapidityBinsH[iRapidityBin] = new TH1F(
                Form("fMCPhiHelicityFrameJPsiRapidityBinsH_%d", iRapidityBin),
                Form("fMCPhiHelicityFrameJPsiRapidityBinsH_%d", iRapidityBin),
                4000, -4., 4.
              );
    fOutputList->Add(fMCPhiHelicityFrameJPsiRapidityBinsH[iRapidityBin]);
  }

  for(Int_t iRapidityBin = 0; iRapidityBin < 8; iRapidityBin++ ){
    fMCPhiCollinsSoperFrameJPsiRapidityBinsH[iRapidityBin] = new TH1F(
                Form("fMCPhiCollinsSoperFrameJPsiRapidityBinsH_%d", iRapidityBin),
                Form("fMCPhiCollinsSoperFrameJPsiRapidityBinsH_%d", iRapidityBin),
                4000, -4., 4.
              );
    fOutputList->Add(fMCPhiCollinsSoperFrameJPsiRapidityBinsH[iRapidityBin]);
  }

  /* - Rapidity-dependent analysis for the helicity of the J/Psi.
     - It is divided into RECONSTRUCTED and GENERATED for both
     - helicity frame and Collins-Soper.
     -
     - The following are the same as above divided in 10 rapidity bins.
     - The reasons for this addition is explained in the file:
     - PWGUD/UPC/AliAnalysisTaskUPCforward.cxx
   */
  for(Int_t iRapidityBin = 0; iRapidityBin < 10; iRapidityBin++ ){
    fCosThetaHelicityFrameJPsiTenRapidityBinsH[iRapidityBin] = new TH1F(
                Form("fCosThetaHelicityFrameJPsiTenRapidityBinsH_%d", iRapidityBin),
                Form("fCosThetaHelicityFrameJPsiTenRapidityBinsH_%d", iRapidityBin),
                // 250, -1., 1.
                1000, -1., 1.
              );
    fOutputList->Add(fCosThetaHelicityFrameJPsiTenRapidityBinsH[iRapidityBin]);
  }

  for(Int_t iRapidityBin = 0; iRapidityBin < 10; iRapidityBin++ ){
    fCosThetaCollinsSoperFrameJPsiTenRapidityBinsH[iRapidityBin] = new TH1F(
                Form("fCosThetaCollinsSoperFrameJPsiTenRapidityBinsH_%d", iRapidityBin),
                Form("fCosThetaCollinsSoperFrameJPsiTenRapidityBinsH_%d", iRapidityBin),
                // 250, -1., 1.
                1000, -1., 1.
              );
    fOutputList->Add(fCosThetaCollinsSoperFrameJPsiTenRapidityBinsH[iRapidityBin]);
  }

  for(Int_t iRapidityBin = 0; iRapidityBin < 10; iRapidityBin++ ){
    fPhiHelicityFrameJPsiTenRapidityBinsH[iRapidityBin] = new TH1F(
                Form("fPhiHelicityFrameJPsiTenRapidityBinsH_%d", iRapidityBin),
                Form("fPhiHelicityFrameJPsiTenRapidityBinsH_%d", iRapidityBin),
                // 1000, -4., 4.
                4000, -4., 4.
              );
    fOutputList->Add(fPhiHelicityFrameJPsiTenRapidityBinsH[iRapidityBin]);
  }

  for(Int_t iRapidityBin = 0; iRapidityBin < 10; iRapidityBin++ ){
    fPhiCollinsSoperFrameJPsiTenRapidityBinsH[iRapidityBin] = new TH1F(
                Form("fPhiCollinsSoperFrameJPsiTenRapidityBinsH_%d", iRapidityBin),
                Form("fPhiCollinsSoperFrameJPsiTenRapidityBinsH_%d", iRapidityBin),
                // 1000, -4., 4.
                4000, -4., 4.
              );
    fOutputList->Add(fPhiCollinsSoperFrameJPsiTenRapidityBinsH[iRapidityBin]);
  }

  for(Int_t iRapidityBin = 0; iRapidityBin < 10; iRapidityBin++ ){
    fMCCosThetaHelicityFrameJPsiTenRapidityBinsH[iRapidityBin] = new TH1F(
                Form("fMCCosThetaHelicityFrameJPsiTenRapidityBinsH_%d", iRapidityBin),
                Form("fMCCosThetaHelicityFrameJPsiTenRapidityBinsH_%d", iRapidityBin),
                // 250, -1., 1.
                1000, -1., 1.
              );
    fOutputList->Add(fMCCosThetaHelicityFrameJPsiTenRapidityBinsH[iRapidityBin]);
  }

  for(Int_t iRapidityBin = 0; iRapidityBin < 10; iRapidityBin++ ){
    fMCCosThetaCollinsSoperFrameJPsiTenRapidityBinsH[iRapidityBin] = new TH1F(
                Form("fMCCosThetaCollinsSoperFrameJPsiTenRapidityBinsH_%d", iRapidityBin),
                Form("fMCCosThetaCollinsSoperFrameJPsiTenRapidityBinsH_%d", iRapidityBin),
                // 250, -1., 1.
                1000, -1., 1.
              );
    fOutputList->Add(fMCCosThetaCollinsSoperFrameJPsiTenRapidityBinsH[iRapidityBin]);
  }

  for(Int_t iRapidityBin = 0; iRapidityBin < 10; iRapidityBin++ ){
    fMCPhiHelicityFrameJPsiTenRapidityBinsH[iRapidityBin] = new TH1F(
                Form("fMCPhiHelicityFrameJPsiTenRapidityBinsH_%d", iRapidityBin),
                Form("fMCPhiHelicityFrameJPsiTenRapidityBinsH_%d", iRapidityBin),
                // 1000, -4., 4.
                4000, -4., 4.
              );
    fOutputList->Add(fMCPhiHelicityFrameJPsiTenRapidityBinsH[iRapidityBin]);
  }

  for(Int_t iRapidityBin = 0; iRapidityBin < 10; iRapidityBin++ ){
    fMCPhiCollinsSoperFrameJPsiTenRapidityBinsH[iRapidityBin] = new TH1F(
                Form("fMCPhiCollinsSoperFrameJPsiTenRapidityBinsH_%d", iRapidityBin),
                Form("fMCPhiCollinsSoperFrameJPsiTenRapidityBinsH_%d", iRapidityBin),
                // 1000, -4., 4.
                4000, -4., 4.
              );
    fOutputList->Add(fMCPhiCollinsSoperFrameJPsiTenRapidityBinsH[iRapidityBin]);
  }

  fInvariantMassDistributionBinsOfCosThetaAndPhiHelicityFrameH =
        new TH2F( "fInvariantMassDistributionBinsOfCosThetaAndPhiHelicityFrameH",
                  "fInvariantMassDistributionBinsOfCosThetaAndPhiHelicityFrameH",
                  80, -1, 1,
                  80, -4, 4
                  );
  fOutputList->Add(fInvariantMassDistributionBinsOfCosThetaAndPhiHelicityFrameH);

  fMCInvariantMassDistributionBinsOfCosThetaAndPhiHelicityFrameH =
        new TH2F( "fMCInvariantMassDistributionBinsOfCosThetaAndPhiHelicityFrameH",
                  "fMCInvariantMassDistributionBinsOfCosThetaAndPhiHelicityFrameH",
                  80, -1, 1,
                  80, -4, 4
                  );
  fOutputList->Add(fMCInvariantMassDistributionBinsOfCosThetaAndPhiHelicityFrameH);

  /* - Variable binning for CosTheta and Phi.
     - Adopting same binning as inclusive people's.
     -
   */
  const Int_t XBINS = 19;
  const Int_t YBINS = 20;
  Double_t CosThetaBinning[ XBINS + 1 ] = { -1. , -0.8, -0.7 , -0.6 , -0.5, -0.4,
                                            -0.3, -0.2, -0.12, -0.04,  0.04, 0.12,
                                             0.2,  0.3,  0.4,   0.5,   0.6,  0.7,
                                             0.8,  1
                                             };
  Double_t PhiBinning[ YBINS + 1 ] = { -3.142, -2.639, -2.136, -1.885, -1.696,
                                       -1.571, -1.445, -1.257, -1.005, -0.502,
                                        0.,     0.502,  1.005,  1.257,  1.445,
                                        1.571,  1.696,  1.885,  2.136,  2.639,
                                        3.142
                                      };
  fCosThetaAndPhiHelicityFrameInclusivePeopleBinningH =
        new TH2F( "fCosThetaAndPhiHelicityFrameInclusivePeopleBinningH",
                  "fCosThetaAndPhiHelicityFrameInclusivePeopleBinningH",
                  XBINS, CosThetaBinning,
                  YBINS, PhiBinning
                  );
  fOutputList->Add(fCosThetaAndPhiHelicityFrameInclusivePeopleBinningH);

  fMCCosThetaAndPhiHelicityFrameInclusivePeopleBinningH =
        new TH2F( "fMCCosThetaAndPhiHelicityFrameInclusivePeopleBinningH",
                  "fMCCosThetaAndPhiHelicityFrameInclusivePeopleBinningH",
                  XBINS, CosThetaBinning,
                  YBINS, PhiBinning
                  );
  fOutputList->Add(fMCCosThetaAndPhiHelicityFrameInclusivePeopleBinningH);


  fInvariantMassDistributionForSignalExtractionHelicityFrameH =
        new TH2F( "fInvariantMassDistributionForSignalExtractionHelicityFrameH",
                  "fInvariantMassDistributionForSignalExtractionHelicityFrameH",
                  10, -1.00, 1.00,
                  10, -3.14, 3.14
                  );
  fOutputList->Add(fInvariantMassDistributionForSignalExtractionHelicityFrameH);

  fMCInvariantMassDistributionForSignalExtractionHelicityFrameH =
        new TH2F( "fMCInvariantMassDistributionForSignalExtractionHelicityFrameH",
                  "fMCInvariantMassDistributionForSignalExtractionHelicityFrameH",
                  10, -1.00, 1.00,
                  10, -3.14, 3.14
                  );
  fOutputList->Add(fMCInvariantMassDistributionForSignalExtractionHelicityFrameH);

  fCosThetaHeFrameForSignalExH = new TH1F("fCosThetaHeFrameForSignalExH", "fCosThetaHeFrameForSignalExH", 40, -1., 1.);
  fOutputList->Add(fCosThetaHeFrameForSignalExH);

  fMCCosThetaHeFrameForSignalExH = new TH1F("fMCCosThetaHeFrameForSignalExH", "fMCCosThetaHeFrameForSignalExH", 40, -1., 1.);
  fOutputList->Add(fMCCosThetaHeFrameForSignalExH);

  fPhiHeFrameForSignalExH = new TH1F("fPhiHeFrameForSignalExH", "fPhiHeFrameForSignalExH", 50, -3.14, 3.14);
  fOutputList->Add(fPhiHeFrameForSignalExH);

  fMCPhiHeFrameForSignalExH = new TH1F("fMCPhiHeFrameForSignalExH", "fMCPhiHeFrameForSignalExH", 50, -3.14, 3.14);
  fOutputList->Add(fMCPhiHeFrameForSignalExH);

  fEfficiencyPerRunH = new TH1F("fEfficiencyPerRunH", "fEfficiencyPerRunH", 3, 0, 3);
  fEfficiencyPerRunH->SetStats(0);
  fEfficiencyPerRunH->SetFillColor(38);
  fEfficiencyPerRunH->LabelsDeflate();
  fOutputList->Add(fEfficiencyPerRunH);

  fMCEfficiencyPerRunH = new TH1F("fMCEfficiencyPerRunH", "fMCEfficiencyPerRunH", 3, 0, 3);
  fMCEfficiencyPerRunH->SetStats(0);
  fMCEfficiencyPerRunH->SetFillColor(38);
  fMCEfficiencyPerRunH->LabelsDeflate();
  fOutputList->Add(fMCEfficiencyPerRunH);

  for(Int_t iRapidityBin = 0; iRapidityBin < 6; iRapidityBin++ ){

        fEfficiencyPerRunRapidityH[iRapidityBin] = new TH1F( Form("fEfficiencyPerRunRapidityH_%d", iRapidityBin),
                                                             Form("fEfficiencyPerRunRapidityH_%d", iRapidityBin),
                                                             3, 0, 3
                                                             );
        fEfficiencyPerRunRapidityH[iRapidityBin]->SetStats(0);
        fEfficiencyPerRunRapidityH[iRapidityBin]->SetFillColor(38);
        fEfficiencyPerRunRapidityH[iRapidityBin]->LabelsDeflate();
        fOutputList->Add(fEfficiencyPerRunRapidityH[iRapidityBin]);



        fMCEfficiencyPerRunRapidityH[iRapidityBin] = new TH1F( Form("fMCEfficiencyPerRunRapidityH_%d", iRapidityBin),
                                                               Form("fMCEfficiencyPerRunRapidityH_%d", iRapidityBin),
                                                               3, 0, 3
                                                               );
        fMCEfficiencyPerRunRapidityH[iRapidityBin]->SetStats(0);
        fMCEfficiencyPerRunRapidityH[iRapidityBin]->SetFillColor(38);
        fMCEfficiencyPerRunRapidityH[iRapidityBin]->LabelsDeflate();
        fOutputList->Add(fMCEfficiencyPerRunRapidityH[iRapidityBin]);
  }

  /* - Eta vs Phi dead zones per Run.
   * -
   */
  Int_t*    bins = new Int_t[2];
  Double_t* xmin = new Double_t[2];
  Double_t* xmax = new Double_t[2];
  /* - [0] refers to Eta.
   * - I am plotting from -5.0 to -2.0,
   * - hence 150 bins are reasonable...  (REBIN 10x)
   * - [1] refers to Phi.
   * - To avoid problems related to TMath::Pi(),
   * - I am plotting from -4 to 4.
   * - Hence I had thought of 200 bins... (REBIN 10x)
   */
  bins[0] = 150;
  bins[1] = 200;
  xmin[0] = -5.0;
  xmin[1] = -4.0;
  xmax[0] = -2.0;
  xmax[1] = +4.0;
  Int_t listOfGoodRunNumbers[]       = { 295585, 295586, 295587, 295588, 295589, 295612,
                                         295615, 295665, 295666, 295667, 295668, 295671,
                                         295673, 295675, 295676, 295677, 295714, 295716,
                                         295717, 295718, 295719, 295723, 295725, 295753,
                                         295754, 295755, 295758, 295759, 295762, 295763,
                                         295786, 295788, 295791, 295816, 295818, 295819,
                                         295822, 295825, 295826, 295829, 295831, 295854,
                                         295855, 295856, 295859, 295860, 295861, 295863,
                                         295881, 295908, 295909, 295910, 295913, 295936,
                                         295937, 295941, 295942, 295943, 295945, 295947,
                                         296061, 296062, 296063, 296065, 296066, 296068,
                                         296123, 296128, 296132, 296133, 296134, 296135,
                                         296142, 296143, 296191, 296192, 296194, 296195,
                                         296196, 296197, 296198, 296241, 296242, 296243,
                                         296244, 296246, 296247, 296269, 296270, 296273,
                                         296279, 296280, 296303, 296304, 296307, 296309,
                                         296312, 296376, 296377, 296378, 296379, 296380,
                                         296381, 296383, 296414, 296419, 296420, 296423,
                                         296424, 296433, 296472, 296509, 296510, 296511,
                                         296514, 296516, 296547, 296548, 296549, 296550,
                                         296551, 296552, 296553, 296615, 296616, 296618,
                                         296619, 296622, 296623,
                                         296690, 296691, 296694, 296749, 296750, 296781,
                                         296784, 296785, 296786, 296787, 296791, 296793,
                                         296794, 296799, 296836, 296838, 296839, 296848,
                                         296849, 296850, 296851, 296852, 296890, 296894,
                                         296899, 296900, 296903, 296930, 296931, 296932,
                                         296934, 296935, 296938, 296941, 296966, 296967,
                                         296968, 296969, 296971, 296975, 296976, 296977,
                                         296979, 297029, 297031, 297035, 297085, 297117,
                                         297118, 297119, 297123, 297124, 297128, 297129,
                                         297132, 297133, 297193, 297194, 297196, 297218,
                                         297219, 297221, 297222, 297278, 297310, 297312,
                                         297315, 297317, 297363, 297366, 297367, 297372,
                                         297379, 297380, 297405, 297408, 297413, 297414,
                                         297415, 297441, 297442, 297446, 297450, 297451,
                                         297452, 297479, 297481, 297483, 297512, 297537,
                                         297540, 297541, 297542, 297544, 297558, 297588,
                                         297590, 297595,/*, 297623, 297624*/
                                         244918, 244980, 244982, 244983, 245064, 245066, 245068, 245145, 245146, 245151,
                                         245152, 245231, 245232, 245233, 245253, 245259, 245343, 245345, 245346, 245347,
                                         245353, 245401, 245407, 245409, 245410, 245446, 245450, 245496, 245501, 245504,
                                         245505, 245507, 245535, 245540, 245542, 245543, 245554, 245683, 245692, 245700,
                                         245705, 245729, 245731, 245738, 245752, 245759, 245766, 245775, 245785, 245793,
                                         245829, 245831, 245833, 245949, 245952, 245954, 245963, 245996, 246001, 246003,
                                         246012, 246036, 246037, 246042, 246048, 246049, 246053, 246087, 246089, 246113,
                                         246115, 246148, 246151, 246152, 246153, 246178, 246181, 246182, 246217, 246220,
                                         246222, 246225, 246272, 246275, 246276, 246390, 246391, 246392, 246424, 246428,
                                         246431, 246433, 246434, 246487, 246488, 246493, 246495, 246675, 246676, 246750,
                                         246751, 246755, 246757, 246758, 246759, 246760, 246763, 246765, 246804, 246805,
                                         246806, 246807, 246808, 246809, 246844, 246845, 246846, 246847, 246851, 246855,
                                         246859, 246864, 246865, 246867, 246871, 246930, 246937, 246942, 246945, 246948,
                                         246949, 246980, 246982, 246984, 246989, 246991, 246994
                                       };
  // for( Int_t iRuns = 0; iRuns < 60001; iRuns++ ) {
  for( Int_t iRuns = 0; iRuns < 364; iRuns++ ) {
    // fDeadZoneEtaVsPhiPerRunH[iRuns] = new THnSparseF( Form( "fDeadZoneEtaVsPhiPerRunH_%d", listOfGoodRunNumbers[iRuns] ),
    //                                                   Form( "fDeadZoneEtaVsPhiPerRunH_%d", listOfGoodRunNumbers[iRuns] ),
    //                                                   2, // number of dimensions
    //                                                   bins,
    //                                                   xmin,
    //                                                   xmax
    //                                                   );
    fDeadZoneEtaVsPhiPerRunH[iRuns] = new TH2F( Form( "fDeadZoneEtaVsPhiPerRunH_%d", listOfGoodRunNumbers[iRuns] ),
                                                Form( "fDeadZoneEtaVsPhiPerRunH_%d", listOfGoodRunNumbers[iRuns] ),
                                                150, -5.0, -2.0,
                                                // 200, -4.0,  4.0
                                                200,  0.0,  8.0
                                                );
    fOutputList->Add(fDeadZoneEtaVsPhiPerRunH[iRuns]);
  }

  //________________________________________
  /* - Templates for polarisation analysis.
   * -
   */
  for(Int_t iCosThetaBins = 0; iCosThetaBins < 40; iCosThetaBins++ ){
    fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameH[iCosThetaBins] = new TH1F(
                Form("fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameH_%d", iCosThetaBins),
                Form("fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameH_%d", iCosThetaBins),
                2000, 0, 20
                );
    fOutputList->Add(fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameH[iCosThetaBins]);
  }

  for(Int_t iPhiBins = 0; iPhiBins < 50; iPhiBins++ ){
    fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameH[iPhiBins] = new TH1F(
                Form("fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameH_%d", iPhiBins),
                Form("fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameH_%d", iPhiBins),
                2000, 0, 20
                );
    fOutputList->Add(fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameH[iPhiBins]);
  }
  //________________________________________

  /* - My Variable binning for CosTheta and Phi.
   * - 2D analysis.
   */
  const Int_t XBINS2 = 7;
  const Int_t YBINS2 = 20;
  Double_t MyVariableCosThetaBinning2[] = { -0.65, -0.35, -0.15, -0.05,
                                             0.05,  0.15,  0.35,  0.65 };
  Double_t MyVariablePhiBinning2[] = { -3.14*1,       -3.14*19*0.05, -3.14*18*0.05, -3.14*17*0.05,
                                       -3.14*13*0.05, -3.14*9*0.05,  -3.14*6*0.05,  -3.14*4*0.05,
                                       -3.14*2*0.05,  -3.14*1*0.05,   0,            +3.14*1*0.05,
                                       +3.14*2*0.05,  +3.14*4*0.05,  +3.14*6*0.05,  +3.14*9*0.05,
                                       +3.14*13*0.05, +3.14*17*0.05, +3.14*18*0.05, +3.14*19*0.05,
                                       +3.14*1 };
  /* -
   * - For flat polarisation.
   * -
   */
  const Int_t XBINS7 = 28;
  const Int_t YBINS7 = 80;
  Double_t MyVariableCosThetaBinning3[] = { -0.65, -0.65+0.25*0.3, -0.65+0.5*0.3, -0.65+0.75*0.3,
                                            -0.35, -0.35+0.25*0.2, -0.35+0.5*0.2, -0.35+0.75*0.2,
                                            -0.15, -0.15+0.25*0.1, -0.15+0.5*0.1, -0.15+0.75*0.1,
                                            -0.05, -0.05+0.25*0.1, -0.05+0.5*0.1, -0.05+0.75*0.1,
                                             0.05,  0.05+0.25*0.1,  0.05+0.5*0.1,  0.05+0.75*0.1,
                                             0.15,  0.15+0.25*0.2,  0.15+0.5*0.2,  0.15+0.75*0.2,
                                             0.35,  0.35+0.25*0.3,  0.35+0.5*0.3,  0.35+0.75*0.3,
                                             0.65 };
  Double_t MyVariablePhiBinning3[] = { -3.14*1.,      -3.14*1.     +3.14*0.05*0.25, -3.14*1.     +3.14*0.05*0.5, -3.14*1.     +3.14*0.05*0.75,
                                       -3.14*19*0.05, -3.14*19*0.05+3.14*0.05*0.25, -3.14*19*0.05+3.14*0.05*0.5, -3.14*19*0.05+3.14*0.05*0.75,
                                       -3.14*18*0.05, -3.14*18*0.05+3.14*0.05*0.25, -3.14*18*0.05+3.14*0.05*0.5, -3.14*18*0.05+3.14*0.05*0.75,
                                       -3.14*17*0.05, -3.14*17*0.05+3.14*0.20*0.25, -3.14*17*0.05+3.14*0.20*0.5, -3.14*17*0.05+3.14*0.20*0.75,
                                       -3.14*13*0.05, -3.14*13*0.05+3.14*0.20*0.25, -3.14*13*0.05+3.14*0.20*0.5, -3.14*13*0.05+3.14*0.20*0.75,
                                       -3.14*9*0.05,  -3.14*9*0.05 +3.14*0.15*0.25, -3.14*9*0.05 +3.14*0.15*0.5, -3.14*9*0.05 +3.14*0.15*0.75,
                                       -3.14*6*0.05,  -3.14*6*0.05 +3.14*0.10*0.25, -3.14*6*0.05 +3.14*0.10*0.5, -3.14*6*0.05 +3.14*0.10*0.75,
                                       -3.14*4*0.05,  -3.14*4*0.05 +3.14*0.10*0.25, -3.14*4*0.05 +3.14*0.10*0.5, -3.14*4*0.05 +3.14*0.10*0.75,
                                       -3.14*2*0.05,  -3.14*2*0.05 +3.14*0.05*0.25, -3.14*2*0.05 +3.14*0.05*0.5, -3.14*2*0.05 +3.14*0.05*0.75,
                                       -3.14*1*0.05,  -3.14*1*0.05 +3.14*0.05*0.25, -3.14*1*0.05 +3.14*0.05*0.5, -3.14*1*0.05 +3.14*0.05*0.75,
                                        0,                          3.14*0.05*0.25,               3.14*0.05*0.5,               3.14*0.05*0.75,
                                       +3.14*1*0.05,  +3.14*1*0.05 +3.14*0.05*0.25, +3.14*1*0.05 +3.14*0.05*0.5, +3.14*1*0.05 +3.14*0.05*0.75,
                                       +3.14*2*0.05,  +3.14*2*0.05 +3.14*0.10*0.25, +3.14*2*0.05 +3.14*0.10*0.5, +3.14*2*0.05 +3.14*0.10*0.75,
                                       +3.14*4*0.05,  +3.14*4*0.05 +3.14*0.10*0.25, +3.14*4*0.05 +3.14*0.10*0.5, +3.14*4*0.05 +3.14*0.10*0.75,
                                       +3.14*6*0.05,  +3.14*6*0.05 +3.14*0.15*0.25, +3.14*6*0.05 +3.14*0.15*0.5, +3.14*6*0.05 +3.14*0.15*0.75,
                                       +3.14*9*0.05,  +3.14*9*0.05 +3.14*0.20*0.25, +3.14*9*0.05 +3.14*0.20*0.5, +3.14*9*0.05 +3.14*0.20*0.75,
                                       +3.14*13*0.05, +3.14*13*0.05+3.14*0.20*0.25, +3.14*13*0.05+3.14*0.20*0.5, +3.14*13*0.05+3.14*0.20*0.75,
                                       +3.14*17*0.05, +3.14*17*0.05+3.14*0.05*0.25, +3.14*17*0.05+3.14*0.05*0.5, +3.14*17*0.05+3.14*0.05*0.75,
                                       +3.14*18*0.05, +3.14*18*0.05+3.14*0.05*0.25, +3.14*18*0.05+3.14*0.05*0.5, +3.14*18*0.05+3.14*0.05*0.75,
                                       +3.14*19*0.05, +3.14*19*0.05+3.14*0.05*0.25, +3.14*19*0.05+3.14*0.05*0.5, +3.14*19*0.05+3.14*0.05*0.75,
                                       +3.14*1 };

  fCosThetaAndPhiHelicityFrameMyBinningH =
        new TH2F( "fCosThetaAndPhiHelicityFrameMyBinningH",
                  "fCosThetaAndPhiHelicityFrameMyBinningH",
                  XBINS2, MyVariableCosThetaBinning2,
                  YBINS2, MyVariablePhiBinning2
                  );
  fOutputList->Add(fCosThetaAndPhiHelicityFrameMyBinningH);

  fMCCosThetaAndPhiHelicityFrameMyBinningH =
        new TH2F( "fMCCosThetaAndPhiHelicityFrameMyBinningH",
                  "fMCCosThetaAndPhiHelicityFrameMyBinningH",
                  XBINS2, MyVariableCosThetaBinning2,
                  YBINS2, MyVariablePhiBinning2
                  );
  fOutputList->Add(fMCCosThetaAndPhiHelicityFrameMyBinningH);

  fCosThetaAndPhiCsFrameMyBinningH =
        new TH2F( "fCosThetaAndPhiCsFrameMyBinningH",
                  "fCosThetaAndPhiCsFrameMyBinningH",
                  XBINS2, MyVariableCosThetaBinning2,
                  YBINS2, MyVariablePhiBinning2
                  );
  fOutputList->Add(fCosThetaAndPhiCsFrameMyBinningH);

  fMCCosThetaAndPhiCsFrameMyBinningH =
        new TH2F( "fMCCosThetaAndPhiCsFrameMyBinningH",
                  "fMCCosThetaAndPhiCsFrameMyBinningH",
                  XBINS2, MyVariableCosThetaBinning2,
                  YBINS2, MyVariablePhiBinning2
                  );
  fOutputList->Add(fMCCosThetaAndPhiCsFrameMyBinningH);


  fCosThetaAndPhiHelicityFrameMyBinningReweightingH =
        new TH2F( "fCosThetaAndPhiHelicityFrameMyBinningReweightingH",
                  "fCosThetaAndPhiHelicityFrameMyBinningReweightingH",
                  XBINS7, MyVariableCosThetaBinning3,
                  YBINS7, MyVariablePhiBinning3
                  );
  fOutputList->Add(fCosThetaAndPhiHelicityFrameMyBinningReweightingH);

  fMCCosThetaAndPhiHelicityFrameMyBinningReweightingH =
        new TH2F( "fMCCosThetaAndPhiHelicityFrameMyBinningReweightingH",
                  "fMCCosThetaAndPhiHelicityFrameMyBinningReweightingH",
                  XBINS7, MyVariableCosThetaBinning3,
                  YBINS7, MyVariablePhiBinning3
                  );
  fOutputList->Add(fMCCosThetaAndPhiHelicityFrameMyBinningReweightingH);

  fCosThetaAndPhiCsFrameMyBinningReweightingH =
        new TH2F( "fCosThetaAndPhiCsFrameMyBinningReweightingH",
                  "fCosThetaAndPhiCsFrameMyBinningReweightingH",
                  XBINS7, MyVariableCosThetaBinning3,
                  YBINS7, MyVariablePhiBinning3
                  );
  fOutputList->Add(fCosThetaAndPhiCsFrameMyBinningReweightingH);

  fMCCosThetaAndPhiCsFrameMyBinningReweightingH =
        new TH2F( "fMCCosThetaAndPhiCsFrameMyBinningReweightingH",
                  "fMCCosThetaAndPhiCsFrameMyBinningReweightingH",
                  XBINS7, MyVariableCosThetaBinning3,
                  YBINS7, MyVariablePhiBinning3
                  );
  fOutputList->Add(fMCCosThetaAndPhiCsFrameMyBinningReweightingH);





  fMCCosThetaAndPhiHelicityFrameReweightingH =
        new TH2F( "fMCCosThetaAndPhiHelicityFrameReweightingH",
                  "fMCCosThetaAndPhiHelicityFrameReweightingH",
                  100, -1., 1.,
                  100, -3.14, 3.14
                  );
  fOutputList->Add(fMCCosThetaAndPhiHelicityFrameReweightingH);

  fMCCosThetaAndPhiCsFrameReweightingH =
        new TH2F( "fMCCosThetaAndPhiCsFrameReweightingH",
                  "fMCCosThetaAndPhiCsFrameReweightingH",
                  100, -1., 1.,
                  100, -3.14, 3.14
                  );
  fOutputList->Add(fMCCosThetaAndPhiCsFrameReweightingH);


  for ( Int_t iTrigger = 0; iTrigger < 7; iTrigger++ ) {
    fCosThetaAndPhiHelicityFrameMyBinningTriggerH[iTrigger] =
          new TH2F( Form("fCosThetaAndPhiHelicityFrameMyBinningTriggerH_%d", iTrigger),
                    Form("fCosThetaAndPhiHelicityFrameMyBinningTriggerH_%d", iTrigger),
                    XBINS2, MyVariableCosThetaBinning2,
                    YBINS2, MyVariablePhiBinning2
                    );
    fOutputList->Add(fCosThetaAndPhiHelicityFrameMyBinningTriggerH[iTrigger]);
  }

  for ( Int_t iTrigger = 0; iTrigger < 7; iTrigger++ ) {
    fCosThetaAndPhiCsFrameMyBinningTriggerH[iTrigger] =
          new TH2F( Form("fCosThetaAndPhiCsFrameMyBinningTriggerH_%d", iTrigger),
                    Form("fCosThetaAndPhiCsFrameMyBinningTriggerH_%d", iTrigger),
                    XBINS2, MyVariableCosThetaBinning2,
                    YBINS2, MyVariablePhiBinning2
                    );
    fOutputList->Add(fCosThetaAndPhiCsFrameMyBinningTriggerH[iTrigger]);
  }


  //_____________________________________________
  /* - My Variable binning for CosTheta and Phi.
   * - 1D analysis.
   */
  const Int_t XBINS3 = 26;
  const Int_t YBINS3 = 30;
  Double_t MyVariableCosThetaBinning1D[] = { -0.65,  -0.5,  -0.4,  -0.35,  -0.3,
                                             -0.25,  -0.2,  -0.15, -0.125, -0.1,
                                             -0.075, -0.05, -0.025, 0,      0.025,
                                              0.05,   0.075, 0.1,   0.125,  0.15,
                                              0.2,    0.25,  0.3,   0.35,   0.4,
                                              0.5,    0.65 };
  Double_t MyVariablePhiBinning1D[] = { -3.14*1,      -3.14*14/15,  -3.14*13/15,  -3.14*12/15,
                                        -3.14*11/15,  -3.14*10/15,  -3.14*9/15,   -3.14*8/15,
                                        -3.14*7/15,   -3.14*6/15,   -3.14*5/15,   -3.14*4/15,
                                        -3.14*3/15,   -3.14*2/15,   -3.14*1/15,    0,
                                        +3.14*1/15,   +3.14*2/15,   +3.14*3/15,   +3.14*4/15,
                                        +3.14*5/15,   +3.14*6/15,   +3.14*7/15,   +3.14*8/15,
                                        +3.14*9/15,   +3.14*10/15,  +3.14*11/15,  +3.14*12/15,
                                        +3.14*13/15,  +3.14*14/15,  +3.14*1 };
  const Int_t XBINS4 = 9;
  Double_t MyVariableCosThetaBinning1Dv2[] = { -0.65, -0.45, -0.3, -0.15,  -0.05,
                                                0.05,  0.15,  0.3,  0.45,   0.65 };
  const Int_t XBINS5 = 17;
  Double_t MyVariableCosThetaBinning1Dv3[] = { -0.65,  -0.5,  -0.425, -0.35,
                                               -0.275, -0.2,  -0.125, -0.075,
                                               -0.025,  0.025, 0.075,  0.125,
                                                0.2,    0.275, 0.35,   0.425,
                                                0.5,    0.65 };
  fCosThetaHelicityFrameMyBinningH =
        new TH1F( "fCosThetaHelicityFrameMyBinningH",
                  "fCosThetaHelicityFrameMyBinningH",
                  XBINS3, MyVariableCosThetaBinning1D
                  );
  fOutputList->Add(fCosThetaHelicityFrameMyBinningH);

  fMCCosThetaHelicityFrameMyBinningH =
        new TH1F( "fMCCosThetaHelicityFrameMyBinningH",
                  "fMCCosThetaHelicityFrameMyBinningH",
                  XBINS3, MyVariableCosThetaBinning1D
                  );
  fOutputList->Add(fMCCosThetaHelicityFrameMyBinningH);

  fCosThetaHelicityFrameMyBinningSmallH =
        new TH1F( "fCosThetaHelicityFrameMyBinningSmallH",
                  "fCosThetaHelicityFrameMyBinningSmallH",
                  XBINS4, MyVariableCosThetaBinning1Dv2
                  );
  fOutputList->Add(fCosThetaHelicityFrameMyBinningSmallH);

  fMCCosThetaHelicityFrameMyBinningSmallH =
        new TH1F( "fMCCosThetaHelicityFrameMyBinningSmallH",
                  "fMCCosThetaHelicityFrameMyBinningSmallH",
                  XBINS4, MyVariableCosThetaBinning1Dv2
                  );
  fOutputList->Add(fMCCosThetaHelicityFrameMyBinningSmallH);

  fCosThetaHelicityFrameMySeventeenBinningH =
        new TH1F( "fCosThetaHelicityFrameMySeventeenBinningH",
                  "fCosThetaHelicityFrameMySeventeenBinningH",
                  XBINS5, MyVariableCosThetaBinning1Dv3
                  );
  fOutputList->Add(fCosThetaHelicityFrameMySeventeenBinningH);

  fMCCosThetaHelicityFrameMySeventeenBinningH =
        new TH1F( "fMCCosThetaHelicityFrameMySeventeenBinningH",
                  "fMCCosThetaHelicityFrameMySeventeenBinningH",
                  XBINS5, MyVariableCosThetaBinning1Dv3
                  );
  fOutputList->Add(fMCCosThetaHelicityFrameMySeventeenBinningH);

  fPhiHelicityFrameMyBinningH =
        new TH1F( "fPhiHelicityFrameMyBinningH",
                  "fPhiHelicityFrameMyBinningH",
                  YBINS3, MyVariablePhiBinning1D
                  );
  fOutputList->Add(fPhiHelicityFrameMyBinningH);

  fMCPhiHelicityFrameMyBinningH =
        new TH1F( "fMCPhiHelicityFrameMyBinningH",
                  "fMCPhiHelicityFrameMyBinningH",
                  YBINS3, MyVariablePhiBinning1D
                  );
  fOutputList->Add(fMCPhiHelicityFrameMyBinningH);

  for(Int_t iCosThetaBins = 0; iCosThetaBins < 26; iCosThetaBins++ ){
    fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameMyVariableBinningH[iCosThetaBins] = new TH1F(
                Form("fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameMyVariableBinningH_%d", iCosThetaBins),
                Form("fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameMyVariableBinningH_%d", iCosThetaBins),
                2000, 0, 20
                );
    fOutputList->Add(fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameMyVariableBinningH[iCosThetaBins]);
  }

  for(Int_t iPhiBins = 0; iPhiBins < 30; iPhiBins++ ){
    fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameMyVariableBinningH[iPhiBins] = new TH1F(
                Form("fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameMyVariableBinningH_%d", iPhiBins),
                Form("fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameMyVariableBinningH_%d", iPhiBins),
                2000, 0, 20
                );
    fOutputList->Add(fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameMyVariableBinningH[iPhiBins]);
  }

  /* - HELICITY FRAME ANALYSIS.
   * -
   */
  fCosThetaHelicityFrameTwentyfiveBinsH =
        new TH1F( "fCosThetaHelicityFrameTwentyfiveBinsH",
                  "fCosThetaHelicityFrameTwentyfiveBinsH",
                  25, -1, 1
                  // 100, -1, 1 // Reweighting
                  );
  fOutputList->Add(fCosThetaHelicityFrameTwentyfiveBinsH);

  fMCCosThetaHelicityFrameTwentyfiveBinsH =
        new TH1F( "fMCCosThetaHelicityFrameTwentyfiveBinsH",
                  "fMCCosThetaHelicityFrameTwentyfiveBinsH",
                  25, -1, 1
                  // 100, -1, 1 // Reweighting
                  );
  fOutputList->Add(fMCCosThetaHelicityFrameTwentyfiveBinsH);

  fPhiHelicityFrameTwentyfiveBinsH =
        new TH1F( "fPhiHelicityFrameTwentyfiveBinsH",
                  "fPhiHelicityFrameTwentyfiveBinsH",
                  25, -3.14, 3.14
                  );
  fOutputList->Add(fPhiHelicityFrameTwentyfiveBinsH);

  fMCPhiHelicityFrameTwentyfiveBinsH =
        new TH1F( "fMCPhiHelicityFrameTwentyfiveBinsH",
                  "fMCPhiHelicityFrameTwentyfiveBinsH",
                  25, -3.14, 3.14
                  );
  fOutputList->Add(fMCPhiHelicityFrameTwentyfiveBinsH);

  fTildePhiHelicityFrameTwentyfiveBinsH =
        new TH1F( "fTildePhiHelicityFrameTwentyfiveBinsH",
                  "fTildePhiHelicityFrameTwentyfiveBinsH",
                  // 25, -3.14*7.0*0.25, 3.14*3.0*0.25
                  25, 0, 2. * 3.14
                  );
  fOutputList->Add(fTildePhiHelicityFrameTwentyfiveBinsH);

  fMCTildePhiHelicityFrameTwentyfiveBinsH =
        new TH1F( "fMCTildePhiHelicityFrameTwentyfiveBinsH",
                  "fMCTildePhiHelicityFrameTwentyfiveBinsH",
                  // 25, -3.14*7.0*0.25, 3.14*3.0*0.25
                  25, 0, 2. * 3.14
                  );
  fOutputList->Add(fMCTildePhiHelicityFrameTwentyfiveBinsH);

  for ( Int_t iTrigger = 0; iTrigger < 7; iTrigger++ ) {
    fCosThetaHelicityFrameTwentyfiveBinsTriggerH[iTrigger] =
          new TH1F( Form("fCosThetaHelicityFrameTwentyfiveBinsTriggerH_%d", iTrigger),
                    Form("fCosThetaHelicityFrameTwentyfiveBinsTriggerH_%d", iTrigger),
                    25, -1, 1
                    );
    fOutputList->Add(fCosThetaHelicityFrameTwentyfiveBinsTriggerH[iTrigger]);
  }

  for ( Int_t iTrigger = 0; iTrigger < 7; iTrigger++ ) {
    fPhiHelicityFrameTwentyfiveBinsTriggerH[iTrigger] =
          new TH1F( Form("fPhiHelicityFrameTwentyfiveBinsTriggerH_%d", iTrigger),
                    Form("fPhiHelicityFrameTwentyfiveBinsTriggerH_%d", iTrigger),
                    25, -3.14, 3.14
                    );
    fOutputList->Add(fPhiHelicityFrameTwentyfiveBinsTriggerH[iTrigger]);
  }

  for ( Int_t iTrigger = 0; iTrigger < 7; iTrigger++ ) {
    fTildePhiHelicityFrameTwentyfiveBinsTriggerH[iTrigger] =
          new TH1F( Form("fTildePhiHelicityFrameTwentyfiveBinsTriggerH_%d", iTrigger),
                    Form("fTildePhiHelicityFrameTwentyfiveBinsTriggerH_%d", iTrigger),
                    25, 0, 2. * 3.14
                    );
    fOutputList->Add(fTildePhiHelicityFrameTwentyfiveBinsTriggerH[iTrigger]);
  }




  fMCCosThetaHeVsCsH = new TH2F( "fMCCosThetaHeVsCsH", "fMCCosThetaHeVsCsH", 100, -1., -1., 100, -1., -1. );
  fOutputList->Add(fMCCosThetaHeVsCsH);

  fMCCosThetaHeVsCsFlatH = new TH2F( "fMCCosThetaHeVsCsFlatH", "fMCCosThetaHeVsCsFlatH", 100, -1., -1., 100, -1., -1. );
  fOutputList->Add(fMCCosThetaHeVsCsFlatH);

  /* - COLLINS-SOPER ANALYSIS
   * -
   */
  fCosThetaCsFrameTwentyfiveBinsH =
        new TH1F( "fCosThetaCsFrameTwentyfiveBinsH",
                  "fCosThetaCsFrameTwentyfiveBinsH",
                  25, -1, 1
                  // 100, -1, 1 // Reweighting
                  );
  fOutputList->Add(fCosThetaCsFrameTwentyfiveBinsH);

  fMCCosThetaCsFrameTwentyfiveBinsH =
        new TH1F( "fMCCosThetaCsFrameTwentyfiveBinsH",
                  "fMCCosThetaCsFrameTwentyfiveBinsH",
                  25, -1, 1
                  // 100, -1, 1 // Reweighting
                  );
  fOutputList->Add(fMCCosThetaCsFrameTwentyfiveBinsH);

  fPhiCsFrameTwentyfiveBinsH =
        new TH1F( "fPhiCsFrameTwentyfiveBinsH",
                  "fPhiCsFrameTwentyfiveBinsH",
                  25, -3.14, 3.14
                  );
  fOutputList->Add(fPhiCsFrameTwentyfiveBinsH);

  fMCPhiCsFrameTwentyfiveBinsH =
        new TH1F( "fMCPhiCsFrameTwentyfiveBinsH",
                  "fMCPhiCsFrameTwentyfiveBinsH",
                  25, -3.14, 3.14
                  );
  fOutputList->Add(fMCPhiCsFrameTwentyfiveBinsH);

  fTildePhiCsFrameTwentyfiveBinsH =
        new TH1F( "fTildePhiCsFrameTwentyfiveBinsH",
                  "fTildePhiCsFrameTwentyfiveBinsH",
                  // 25, -3.14*7.0*0.25, 3.14*3.0*0.25
                  25, 0, 2. * 3.14
                  );
  fOutputList->Add(fTildePhiCsFrameTwentyfiveBinsH);

  fMCTildePhiCsFrameTwentyfiveBinsH =
        new TH1F( "fMCTildePhiCsFrameTwentyfiveBinsH",
                  "fMCTildePhiCsFrameTwentyfiveBinsH",
                  // 25, -3.14*7.0*0.25, 3.14*3.0*0.25
                  25, 0, 2. * 3.14
                  );
  fOutputList->Add(fMCTildePhiCsFrameTwentyfiveBinsH);

  for ( Int_t iTrigger = 0; iTrigger < 7; iTrigger++ ) {
    fCosThetaCsFrameTwentyfiveBinsTriggerH[iTrigger] =
          new TH1F( Form("fCosThetaCsFrameTwentyfiveBinsTriggerH_%d", iTrigger),
                    Form("fCosThetaCsFrameTwentyfiveBinsTriggerH_%d", iTrigger),
                    25, -1, 1
                    );
    fOutputList->Add(fCosThetaCsFrameTwentyfiveBinsTriggerH[iTrigger]);
  }

  for ( Int_t iTrigger = 0; iTrigger < 7; iTrigger++ ) {
    fPhiCsFrameTwentyfiveBinsTriggerH[iTrigger] =
          new TH1F( Form("fPhiCsFrameTwentyfiveBinsTriggerH_%d", iTrigger),
                    Form("fPhiCsFrameTwentyfiveBinsTriggerH_%d", iTrigger),
                    25, -3.14, 3.14
                    );
    fOutputList->Add(fPhiCsFrameTwentyfiveBinsTriggerH[iTrigger]);
  }

  for ( Int_t iTrigger = 0; iTrigger < 7; iTrigger++ ) {
    fTildePhiCsFrameTwentyfiveBinsTriggerH[iTrigger] =
          new TH1F( Form("fTildePhiCsFrameTwentyfiveBinsTriggerH_%d", iTrigger),
                    Form("fTildePhiCsFrameTwentyfiveBinsTriggerH_%d", iTrigger),
                    25, 0, 2. * 3.14
                    );
    fOutputList->Add(fTildePhiCsFrameTwentyfiveBinsTriggerH[iTrigger]);
  }


  /* - Invariant mass distributions for signal extraction for POLARISATION.
   * - The usage will be:    histo[CosTheta][Phi];
   * - My variable binning.
   * -
   * - NOTE: this is the official binning I am using
   * -       for the analysis.
   * -       I am doing both a HELICITY frame analysis
   * -       and a COLLINS-SOPER.
   * -
   */
  fInvariantMassDistributionForSignalExtractionHelicityFrameMyBinningH = new TH1F**[7];
  for( Int_t iCosTheta = 0; iCosTheta < 7; iCosTheta++ ){
    fInvariantMassDistributionForSignalExtractionHelicityFrameMyBinningH[iCosTheta] = new TH1F*[20];
    for( Int_t iPhi = 0; iPhi < 20; iPhi++ ){
      fInvariantMassDistributionForSignalExtractionHelicityFrameMyBinningH[iCosTheta][iPhi] =
          new TH1F( Form("fInvariantMassDistributionForSignalExtractionHelicityFrameMyBinningH_%d_%d", iCosTheta, iPhi),
                    Form("fInvariantMassDistributionForSignalExtractionHelicityFrameMyBinningH_%d_%d", iCosTheta, iPhi),
                    2000, 0, 20
                    );
      fOutputList->Add(fInvariantMassDistributionForSignalExtractionHelicityFrameMyBinningH[iCosTheta][iPhi]);
    }
  }

  fInvariantMassDistributionForSignalExtractionCsFrameMyBinningH = new TH1F**[7];
  for( Int_t iCosTheta = 0; iCosTheta < 7; iCosTheta++ ){
    fInvariantMassDistributionForSignalExtractionCsFrameMyBinningH[iCosTheta] = new TH1F*[20];
    for( Int_t iPhi = 0; iPhi < 20; iPhi++ ){
      fInvariantMassDistributionForSignalExtractionCsFrameMyBinningH[iCosTheta][iPhi] =
          new TH1F( Form("fInvariantMassDistributionForSignalExtractionCsFrameMyBinningH_%d_%d", iCosTheta, iPhi),
                    Form("fInvariantMassDistributionForSignalExtractionCsFrameMyBinningH_%d_%d", iCosTheta, iPhi),
                    2000, 0, 20
                    );
      fOutputList->Add(fInvariantMassDistributionForSignalExtractionCsFrameMyBinningH[iCosTheta][iPhi]);
    }
  }

  /* -
   * - FINAL POLARISATION ANALYSIS.
   * - HELICITY
   */
  for(Int_t iCosThetaBins = 0; iCosThetaBins < 25; iCosThetaBins++ ){
    fInvariantMassDistributionOnlyCosThetaHeFrameTwentyfiveBinsH[iCosThetaBins] = new TH1F(
                Form("fInvariantMassDistributionOnlyCosThetaHeFrameTwentyfiveBinsH_%d", iCosThetaBins),
                Form("fInvariantMassDistributionOnlyCosThetaHeFrameTwentyfiveBinsH_%d", iCosThetaBins),
                2000, 0, 20
                );
    fOutputList->Add(fInvariantMassDistributionOnlyCosThetaHeFrameTwentyfiveBinsH[iCosThetaBins]);
  }

  for(Int_t iPhiBins = 0; iPhiBins < 25; iPhiBins++ ){
    fInvariantMassDistributionOnlyPhiHeFrameTwentyfiveBinsH[iPhiBins] = new TH1F(
                Form("fInvariantMassDistributionOnlyPhiHeFrameTwentyfiveBinsH_%d", iPhiBins),
                Form("fInvariantMassDistributionOnlyPhiHeFrameTwentyfiveBinsH_%d", iPhiBins),
                2000, 0, 20
                );
    fOutputList->Add(fInvariantMassDistributionOnlyPhiHeFrameTwentyfiveBinsH[iPhiBins]);
  }

  for(Int_t iPhiBins = 0; iPhiBins < 25; iPhiBins++ ){
    fInvariantMassDistributionOnlyTildePhiHeFrameTwentyfiveBinsH[iPhiBins] = new TH1F(
                Form("fInvariantMassDistributionOnlyTildePhiHeFrameTwentyfiveBinsH_%d", iPhiBins),
                Form("fInvariantMassDistributionOnlyTildePhiHeFrameTwentyfiveBinsH_%d", iPhiBins),
                2000, 0, 20
                );
    fOutputList->Add(fInvariantMassDistributionOnlyTildePhiHeFrameTwentyfiveBinsH[iPhiBins]);
  }

  /* -
   * - FINAL POLARISATION ANALYSIS.
   * - COLLINS-SOPER
   */
  for(Int_t iCosThetaBins = 0; iCosThetaBins < 25; iCosThetaBins++ ){
    fInvariantMassDistributionOnlyCosThetaCsFrameTwentyfiveBinsH[iCosThetaBins] = new TH1F(
                Form("fInvariantMassDistributionOnlyCosThetaCsFrameTwentyfiveBinsH_%d", iCosThetaBins),
                Form("fInvariantMassDistributionOnlyCosThetaCsFrameTwentyfiveBinsH_%d", iCosThetaBins),
                2000, 0, 20
                );
    fOutputList->Add(fInvariantMassDistributionOnlyCosThetaCsFrameTwentyfiveBinsH[iCosThetaBins]);
  }

  for(Int_t iPhiBins = 0; iPhiBins < 25; iPhiBins++ ){
    fInvariantMassDistributionOnlyPhiCsFrameTwentyfiveBinsH[iPhiBins] = new TH1F(
                Form("fInvariantMassDistributionOnlyPhiCsFrameTwentyfiveBinsH_%d", iPhiBins),
                Form("fInvariantMassDistributionOnlyPhiCsFrameTwentyfiveBinsH_%d", iPhiBins),
                2000, 0, 20
                );
    fOutputList->Add(fInvariantMassDistributionOnlyPhiCsFrameTwentyfiveBinsH[iPhiBins]);
  }

  for(Int_t iPhiBins = 0; iPhiBins < 25; iPhiBins++ ){
    fInvariantMassDistributionOnlyTildePhiCsFrameTwentyfiveBinsH[iPhiBins] = new TH1F(
                Form("fInvariantMassDistributionOnlyTildePhiCsFrameTwentyfiveBinsH_%d", iPhiBins),
                Form("fInvariantMassDistributionOnlyTildePhiCsFrameTwentyfiveBinsH_%d", iPhiBins),
                2000, 0, 20
                );
    fOutputList->Add(fInvariantMassDistributionOnlyTildePhiCsFrameTwentyfiveBinsH[iPhiBins]);
  }


  //_______________________________
  // - End of the function
  PostData(1, fOutputList);
}
//_____________________________________________________________________________
void AliAnalysisTaskUPCforwardMC::NotifyRun()
{
  /// Set run number for cuts
  fMuonTrackCuts->SetRun(fInputHandler);
}
//_____________________________________________________________________________
void AliAnalysisTaskUPCforwardMC::UserExec(Option_t *)
{
  /* - This iSelectionCounter is used as a token. So at every passing step it is
     - increased by one. All the events are supposed to pass the first step
     - obviously, but only a few get to the end. This is effect is clearly
     - noticeable in fCounterH event with the small trial local sample.
     - Almost 160k possible events at the 0-th step, while only 2k at the 4th step.
   */
  Int_t iSelectionCounter = 0; // no selection applied yet
  fCounterH->Fill(iSelectionCounter); // entering UserExec
  iSelectionCounter++;

  // get AOD event
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if(!fAOD) {
      PostData(1, fOutputList);
      return;
  }
  fMCEvent = MCEvent();
  if(!fMCEvent) {
      PostData(1, fOutputList);
      return;
  }
  if(fMCEvent) {
    fRunNum    = fAOD->GetRunNumber();
    SetLuminosityCap();
    fCounterGeneratedLevel[ fRunNum - 240000 ] += 1;
    // cout << "fCounterGeneratedLevel[ " << (fRunNum - 240000) << " ] = " << fCounterGeneratedLevel[ fRunNum - 240000 ] << endl;
    // if( fCounterGeneratedLevel[ fRunNum - 240000 ] > ( (Int_t)fLumiPerRun * (Int_t)40000 ) ) {
    if( fCounterGeneratedLevel[ fRunNum - 240000 ] > ( fLumiPerRun * 40000 ) ) {
          PostData(1, fOutputList);
          return;
    }
    ProcessMCParticles(fMCEvent);
    fCounterUPCevent += 1;
    fMCEfficiencyPerRunH->Fill( Form("%d", fRunNum) , 1 );
  }
  /* - We are now checking if there were any tracks. If there were at least one,
     - then the histogram gets filled again. If not we are returning. There
     - would be no point in going further.
   */
  Int_t nTracks(fAOD->GetNumberOfTracks());
  if(nTracks<1) {
        PostData(1, fOutputList);
        return;
  }
  fCounterH->Fill(iSelectionCounter); // At least one track
  iSelectionCounter++;


  //_______________________________
  // EVENT DATA EXTRACTION
  /* - Eugeny Krishen's event data extraction. I am trying to implement it.
     - The only thing I am a bit worried about is whether it should go before or
     - after the "nTracks<1" check... I will try and switch it if it sounds
     - better. These data are used for the event selection and maybe later on
     - for track selection, but I did not get to that part yet. If after all of
     - this I remember to do so, I will come back to this point and correct this
     - statement. If you find this part, please, keep in mind to check the
     - following.
   */

  /* - Event information:
     - Run Number, maybe to select the GOOD Runs and discard the others;
     - Number of Tracklets, these are in this case the SPD tracklets, so the
     - almost unit vector roughly 2 cm between two pixels of the SPD in different
     - layers.
   */
  fRunNum    = fAOD->GetRunNumber();
  fTracklets = fAOD->GetTracklets()->GetNumberOfTracklets();

  /* - Trigger Inputs:
     - L0: ..... ;
     - L1: ..... .
   */
  // fL0inputs = fAOD->GetHeader()->GetL0TriggerInputs();
  // fL1inputs = fAOD->GetHeader()->GetL1TriggerInputs();

  /* - Past-future protection maps:
     - IR1: .... ;
     - IR2: .... .
   */
  // fIR1Map = fAOD->GetHeader()->GetIRInt1InteractionMap();
  // fIR2Map = fAOD->GetHeader()->GetIRInt2InteractionMap();

  /* - ZDC: we try to find the ZDC object data in the nano-AOD. If we cannot,
     - we return, because there would be no way to actually select the events
     - otherwise! We are here, so we could even check if there is a discrepancy
     - between good events with and without ZDC's information. Or at least, this
     - is my impression of it (filling fCounterH). ZDC information:
     - fZem1Energy:
     - fZem2Energy:
     - fZNAEnergy:
     - fZNCEnergy:
     - fZPAEnergy:
     - fZPCEnergy:
     - fZNATime:
     - fZNCTime:
     - fZNATDC[i]:
     - fZNCTDC[i]:
     - fZPATDC[i]:
     - fZPCTDC[i]:
   */
  AliAODZDC *dataZDC = dynamic_cast<AliAODZDC*>(fAOD->GetZDCData());
  if(!dataZDC) {
        PostData(1, fOutputList);
        return;
  }
  fCounterH->Fill(iSelectionCounter);
  iSelectionCounter++;

  fZem1Energy = dataZDC->GetZEM1Energy();
  fZem2Energy = dataZDC->GetZEM2Energy();
  fZNAEnergy  = dataZDC->GetZNATowerEnergy()[0];
  fZNCEnergy  = dataZDC->GetZNCTowerEnergy()[0];
  fZPAEnergy  = dataZDC->GetZPATowerEnergy()[0];
  fZPCEnergy  = dataZDC->GetZPCTowerEnergy()[0];

  fZNATime    = dataZDC->GetZNATime();
  fZNCTime    = dataZDC->GetZNCTime();

  for (Int_t i=0;i<4;i++) fZNATDC[i] = dataZDC->GetZNATDCm(i);
  for (Int_t i=0;i<4;i++) fZNCTDC[i] = dataZDC->GetZNCTDCm(i);
  for (Int_t i=0;i<4;i++) fZPATDC[i] = dataZDC->GetZPATDCm(i);
  for (Int_t i=0;i<4;i++) fZPCTDC[i] = dataZDC->GetZPCTDCm(i);

  /* - These lines are the calibration for the ZDC as provided by Evgeny Kryshen.
     -
   */
  // Bool_t calibrated = 0;
  // if ( fRunNum <  295726 ) calibrated = 1;
  // if ( fRunNum == 296509 ) calibrated = 1;
  // if ( fRunNum >  296689 ) calibrated = 1;
  // if ( fRunNum >  296695 ) calibrated = 0;
  // if ( fRunNum == 297219 ) calibrated = 1;
  // if ( fRunNum == 297221 ) calibrated = 1;
  // if ( fRunNum == 297415 ) calibrated = 1;
  //
  // if ( !calibrated ) {
  //   fZNAEnergy *= (2500./190.);
  //   fZNCEnergy *= (2500./190.);
  // }

  // V0
  AliVVZERO *dataVZERO = dynamic_cast<AliVVZERO*>(fAOD->GetVZEROData());
  if(!dataVZERO) {
        PostData(1, fOutputList);
        return;
  }
  fCounterH->Fill(iSelectionCounter);
  iSelectionCounter++;

  fV0ADecision = dataVZERO->GetV0ADecision();
  fV0CDecision = dataVZERO->GetV0CDecision();

  // Reset event info
  fBBCFlags = 0;
  fBGCFlags = 0;
  fBBAFlags = 0;
  fBGAFlags = 0;
  for (Int_t i=0; i<64; i++) {
    // get array of fired cells
    fBBFlag[i] = dataVZERO->GetBBFlag(i);
    fBGFlag[i] = dataVZERO->GetBGFlag(i);
  }

  for (Int_t i=0; i<32; i++){ // loop over cells
    fBBCFlags+=fBBFlag[i];
    fBGCFlags+=fBGFlag[i];
    fBBAFlags+=fBBFlag[i+32];
    fBGAFlags+=fBGFlag[i+32];
  }


  //_____________________________________
  // RUN SELECTION
  /* - This part is the run selection. We call the std::find() method of the
     - <algorithm> library STL. This returns a kTRUE if you find the value you
     - are looking for inside the vector. So, what happens is that I look for
     - the Run Numbers inside the vector containing them. If I cannot find
     - them I move on to the next event.
     -
   */
  // auto findRunNumber = std::find(  std::begin(fVectorGoodRunNumbers),
  //                                  std::end(fVectorGoodRunNumbers),
  //                                  fRunNum
  //                                  );
  // if (findRunNumber != std::end(fVectorGoodRunNumbers)) {
  //     // std::cout << "fVectorGoodRunNumbers DOES     contain: " << fRunNum << std::endl;
  // } else {
  //     PostData(1, fOutputList);
  //     return;
  // }

  fCounterH->Fill(15);
  Int_t listOfGoodRunNumbersLHC18q[] = { 295585, 295586, 295587, 295588, 295589, 295612,
                                         295615, 295665, 295666, 295667, 295668, 295671,
                                         295673, 295675, 295676, 295677, 295714, 295716,
                                         295717, 295718, 295719, 295723, 295725, 295753,
                                         295754, 295755, 295758, 295759, 295762, 295763,
                                         295786, 295788, 295791, 295816, 295818, 295819,
                                         295822, 295825, 295826, 295829, 295831, 295854,
                                         295855, 295856, 295859, 295860, 295861, 295863,
                                         295881, 295908, 295909, 295910, 295913, 295936,
                                         295937, 295941, 295942, 295943, 295945, 295947,
                                         296061, 296062, 296063, 296065, 296066, 296068,
                                         296123, 296128, 296132, 296133, 296134, 296135,
                                         296142, 296143, 296191, 296192, 296194, 296195,
                                         296196, 296197, 296198, 296241, 296242, 296243,
                                         296244, 296246, 296247, 296269, 296270, 296273,
                                         296279, 296280, 296303, 296304, 296307, 296309,
                                         296312, /*296376,*/ 296377, 296378, 296379, 296380,
                                         296381, 296383, 296414, 296419, 296420, 296423,
                                         296424, 296433, 296472, 296509, 296510, 296511,
                                         296514, 296516, 296547, 296548, 296549, 296550,
                                         296551, 296552, 296553, 296615, 296616, 296618,
                                         296619, 296622, 296623 };
  Int_t listOfGoodRunNumbersLHC18r[] = { 296690, 296691, 296694, 296749, 296750, 296781,
                                         296784, 296785, 296786, 296787, 296791, 296793,
                                         296794, 296799, 296836, 296838, 296839, 296848,
                                         296849, 296850, 296851, 296852, 296890, 296894,
                                         296899, 296900, 296903, 296930, 296931, 296932,
                                         296934, 296935, 296938, 296941, 296966, 296967,
                                         296968, 296969, 296971, 296975, 296976, /*296977,*/
                                         296979, 297029, 297031, 297035, 297085, 297117,
                                         297118, 297119, 297123, 297124, 297128, 297129,
                                         297132, 297133, 297193, 297194, 297196, 297218,
                                         297219, 297221, 297222, 297278, 297310, 297312,
                                         297315, 297317, 297363, 297366, 297367, 297372,
                                         297379, 297380, 297405, 297408, 297413, 297414,
                                         297415, 297441, 297442, 297446, 297450, 297451,
                                         297452, 297479, 297481, 297483, 297512, 297537,
                                         297540, 297541, 297542, 297544, 297558, 297588,
                                         297590, 297595/*, 297623, 297624*/ };
  /* - This good run number list has been taken from the analysis
     - note of Kay's talk for DIS 2017, see:
     - https://alice-notes.web.cern.ch/system/files/notes/analysis/596/2017-Feb-08-analysis_note-2017-Feb-08-analysis-note.pdf
     -
   */
  Int_t listOfGoodRunNumbersLHC15o[] = { /*244918,*/ 244980, 244982, 244983, 245064, 245066, 245068, 245145, 245146, 245151,
                                         245152, 245231, 245232, 245233, 245253, 245259, 245343, 245345, 245346, 245347,
                                         245353, 245401, 245407, 245409, 245410, 245446, 245450, 245496, 245501, 245504,
                                         245505, 245507, 245535, 245540, 245542, 245543, 245554, 245683, 245692, 245700,
                                         245705, 245729, 245731, 245738, 245752, 245759, 245766, 245775, 245785, 245793,
                                         245829, 245831, 245833, 245949, 245952, 245954, 245963, 245996, 246001, 246003,
                                         246012, 246036, 246037, 246042, 246048, 246049, 246053, 246087, 246089, 246113,
                                         246115, 246148, 246151, 246152, 246153, 246178, 246181, 246182, 246217, 246220,
                                         246222, 246225, 246272, 246275, 246276, 246390, 246391, 246392, 246424, 246428,
                                         246431, 246433, 246434, 246487, 246488, 246493, 246495, 246675, 246676, 246750,
                                         246751, 246755, 246757, 246758, 246759, 246760, 246763, 246765, 246804, 246805,
                                         246806, 246807, 246808, 246809, 246844, 246845, 246846, 246847, 246851, 246855,
                                         246859, 246864, 246865, 246867, 246871, 246930, 246937, 246942, 246945, 246948,
                                         246949, 246980, 246982, 246984, 246989, 246991, 246994
                                       };
  Int_t listOfRunNumbersZDC[] = { 296244, 296750, 296849, 297219, 297481 };
  Bool_t checkIfGoodRun = kFALSE;
  // for( Int_t iRunLHC18q = 0; iRunLHC18q < 129; iRunLHC18q++){
  for( Int_t iRunLHC18q = 0; iRunLHC18q < 128; iRunLHC18q++){
    if( fRunNum == listOfGoodRunNumbersLHC18q[iRunLHC18q] ) checkIfGoodRun = kTRUE;
  }
  // for( Int_t iRunLHC18r = 0; iRunLHC18r <  98; iRunLHC18r++){
  for( Int_t iRunLHC18r = 0; iRunLHC18r <  97; iRunLHC18r++){
    if( fRunNum == listOfGoodRunNumbersLHC18r[iRunLHC18r] ) checkIfGoodRun = kTRUE;
  }
  // for( Int_t iRunLHC15o = 0; iRunLHC15o < 137; iRunLHC15o++){
  for( Int_t iRunLHC15o = 0; iRunLHC15o < 136; iRunLHC15o++){
    if( fRunNum == listOfGoodRunNumbersLHC15o[iRunLHC15o] ) checkIfGoodRun = kTRUE;
  }
  // for( Int_t iRunZDC = 0; iRunZDC < 5; iRunZDC++){
  //   if( fRunNum == listOfRunNumbersZDC[iRunZDC] )           checkIfGoodRun = kTRUE;
  // }
  // cout << "fRunNum = " << fRunNum << "   and   checkIfGoodRun = " << checkIfGoodRun << endl;
  if(checkIfGoodRun != 1) {
       PostData(1, fOutputList);
       // cout << "OPS!" << endl;
       return;
  }
  fCounterH->Fill(17);

  // END RUN SELECTION
  //_____________________________________



  /* - We have to get the number of fired V0C cells. So firstly, we get the
     - boolean information about the hit cells for all V0. This is done through
     - the GetBBFlag(i) method, where 0<i<32 stands for the V0C cells and
     - 32<i<64 for the V0A cells. Then I thought the easiest way to check
     - whether the number of fired V0C cells is above 2 is just to add up the
     - boolean numbers for 0<i<32. Let's see.
   */
  fV0TotalNCells = 0;
  for(Int_t iV0Hits = 0; iV0Hits < 64; iV0Hits++) {
        fV0Hits[iV0Hits] = dataVZERO->GetBBFlag(iV0Hits);
        if(fV0Hits[iV0Hits] == kTRUE) {
              if(iV0Hits < 32) fV0TotalNCells += 1;
        }
  }

  // AD
  AliVAD *dataAD = dynamic_cast<AliVAD*>(fAOD->GetADData());
  if(dataAD) {
        fCounterH->Fill(iSelectionCounter);
        iSelectionCounter++;

        fADADecision = dataAD->GetADADecision();
        fADCDecision = dataAD->GetADCDecision();

        // Reset event info
        fBBCFlagsAD = 0;
        fBGCFlagsAD = 0;
        fBBAFlagsAD = 0;
        fBGAFlagsAD = 0;
        for(Int_t i=0; i<16; i++) {
          // get array of fired pads
          fBBFlagAD[i] = dataAD->GetBBFlag(i);
          fBGFlagAD[i] = dataAD->GetBGFlag(i);
        }

        for(Int_t i=0; i<4; i++) { // loop over pairs of pads
          if ( fBBFlagAD[i]   && fBBFlagAD[i+4]  ) fBBCFlagsAD++;
          if ( fBGFlagAD[i]   && fBGFlagAD[i+4]  ) fBGCFlagsAD++;
          if ( fBBFlagAD[i+8] && fBBFlagAD[i+12] ) fBBAFlagsAD++;
          if ( fBGFlagAD[i+8] && fBGFlagAD[i+12] ) fBGAFlagsAD++;
        }
  }
  // END EVENT DATA EXTRACTION
  //_______________________________
  // APPLY TRIGGER MC!
  if(!IsTriggered()) {
    cout << "Ehm" ;
    PostData(1, fOutputList);
    return;
  }
  fCounterH->Fill(iSelectionCounter); // right trigger found
  iSelectionCounter++;
  //_______________________________
  // EVENT SELECTION
  /* - This is Eugeny Krishen's event selection from the talk in 14/1/2019 for
     - the PWG-UD (UPC oriented) meeting. The event selection requires:
     - CMUP11-B triggers;
     - Maximum 2 V0C cells fired;
     - Empty V0A decision;
     - Empty ADA decision;
     - Empty ADC decision;
     - 0 tracklets in SPD;
     - Exactly 2 unlike-sign muons;
   */
  /* - Empty V0A decision
     - Empty ADA decision
     - Empty ADC decision
   */
  if(fV0ADecision != 0) {
       PostData(1, fOutputList);
       return;
  }
  if(fADADecision != 0) {
       PostData(1, fOutputList);
       return;
  }
  if(fADCDecision != 0) {
       PostData(1, fOutputList);
       return;
  }
  /* - 0 tracklets in SPD
   */
  if(fTracklets != 0) {
       PostData(1, fOutputList);
       return;
  }
  /* - Maximum 2 V0C cells fired.
   */
  if( fV0TotalNCells > 2 ) {
       PostData(1, fOutputList);
       return;
  }

  /* - We are finally at the starting point. We loop over the tracks and select
     - the good muons. Later on everything should happen in this loop. Let us
     - see what the future has in hold.
     -
     - Saturday: I moved the creation of the AliAODTrack* track outside of the
     - loop as it would have been otherwise created for each single iteration.
     - This could have caused massive memory issues especially to grid. I have
     - added a second AliAODTrack* track[2] to hold the second supposed muon.
     - Now this is ready to send the information to two TLorentzVectors to
     - obtain the invariant mass of the J/Psi through the Mag() method of the
     - class. Hope for the best.
   */
  Int_t nGoodMuons = 0;
  AliAODTrack* track[2];
  track[0]         = 0x0;
  track[1]         = 0x0;
  for(Int_t iTrack(0); iTrack < nTracks; iTrack++) {
    /* - This should be another form of event selection.
       - I am basically requesting the presence of TWO good muons only.
       - Later I will be checking whether of they are likesign or unlikesign.
     */
    if(nGoodMuons > 2) {
         PostData(1, fOutputList);
         return;
    }
    track[nGoodMuons] = static_cast<AliAODTrack*>(fAOD->GetTrack(iTrack));
    if(!track[nGoodMuons]) return;

    // is it a good muon track?
    if(!track[nGoodMuons]->IsMuonTrack()) {
        // track[nGoodMuons] = 0x0;
        continue;
    }
    if(!fMuonTrackCuts->IsSelected(track[nGoodMuons])) {
        // track[nGoodMuons] = 0x0;
        continue;
    }

    // MUON SELECTION
    /* - This is Eugeny Krishen's MUON selection from the talk in 14/1/2019 for
       - the PWG-UD (UPC oriented) meeting. The event selection requires:
       - Muon trigger matching >=2 (1 GeV/c threshold);
       - (-4) < eta < (-2.5);
       - (17.5 cm) < R_{abs} < (89.5 cm);
       - p \times DCA cut;
    */

    // increase counter
    nGoodMuons++;
  }
  /* - We need EXACTLY 2 good muons !!!!!
     -
   */
  if( nGoodMuons != 2 ) {
        PostData(1, fOutputList);
        return;
  }
  /* - Implementing the track cut on the unlike muons
   * -
   */
  if( (track[0]->Charge()) == (track[1]->Charge()) ) {
        PostData(1, fOutputList);
        return;
  }
  for(Int_t iFilling = 0; iFilling < nGoodMuons; iFilling++) {
        fEtaMuonH ->Fill(track[iFilling]->Eta());
        fRAbsMuonH->Fill(track[iFilling]->GetRAtAbsorberEnd());
  }
  // store muons
  fNumberMuonsH->Fill(nGoodMuons);
  fEntriesAgainstRunNumberH->Fill(fRunNum);
  /* - This is the last part of my try to obtain a proper RunNumbers histogram...
     -
   */
  fEntriesAgainstRunNumberProperlyH->Fill( Form("%d", fRunNum) , 1 );
  fEfficiencyPerRunH               ->Fill( Form("%d", fRunNum) , 1 );
  if (nGoodMuons>0) fCounterH->Fill(iSelectionCounter); // At least one good muon
  iSelectionCounter++;

  /* - Filling the fDeadZoneEtaVsPhiPerRunH.
   * -
   */
  // fEtaAndPhi[0] = 0;
  // fEtaAndPhi[1] = 0;
  // fEtaAndPhi[0] = track[0]->Eta();
  // fEtaAndPhi[1] = track[0]->Phi();
  // ((THnSparseF*) fOutputList->FindObject(Form( "fDeadZoneEtaVsPhiPerRunH_%d", fRunNum )) )->Fill( fEtaAndPhi );
  ((TH2F*) fOutputList->FindObject(Form( "fDeadZoneEtaVsPhiPerRunH_%d", fRunNum )) )->Fill( track[0]->Eta(), track[0]->Phi() );
  // fEtaAndPhi[0] = 0;
  // fEtaAndPhi[1] = 0;
  // fEtaAndPhi[0] = track[1]->Eta();
  // fEtaAndPhi[1] = track[1]->Phi();
  // ((THnSparseF*) fOutputList->FindObject(Form( "fDeadZoneEtaVsPhiPerRunH_%d", fRunNum )) )->Fill( fEtaAndPhi );
  ((TH2F*) fOutputList->FindObject(Form( "fDeadZoneEtaVsPhiPerRunH_%d", fRunNum )) )->Fill( track[1]->Eta(), track[1]->Phi() );

  /* - Finally the core!!!
   * - What will be happening is that we will instantiate TLorentzVectors to
   * - obtain the invariant mass of the dimuon system. If everything goes fine
   * - after this we should be able to obtain the peak of the J/Psi. But
   * - things never go as expected, so who knows!
   */
  TLorentzVector muons[2];
  TLorentzVector possibleJPsi;
  Double_t       chargeOfMuons[2];
  for(int indexMuon = 0; indexMuon < 2; indexMuon++) {
        muons[indexMuon].SetPtEtaPhiM(   track[indexMuon]->Pt(),
                                         track[indexMuon]->Eta(),
                                         track[indexMuon]->Phi(),
                                         TDatabasePDG::Instance()->GetParticle(13)->Mass()
                                       );
        possibleJPsi += muons[indexMuon];
        chargeOfMuons[indexMuon] = track[indexMuon]->Charge();
  }
  fInvariantMassDistributionH->Fill(possibleJPsi.Mag());
  if (        possibleJPsi.Rapidity() > -4.0  && possibleJPsi.Rapidity() <= -3.75 ) {
    fInvariantMassDistributionRapidityBinsH[0]->Fill(possibleJPsi.Mag());
    fEfficiencyPerRunRapidityH[0]             ->Fill( Form("%d", fRunNum) , 1 );
  } else if ( possibleJPsi.Rapidity() > -3.75 && possibleJPsi.Rapidity() <= -3.50 ) {
    fInvariantMassDistributionRapidityBinsH[1]->Fill(possibleJPsi.Mag());
    fEfficiencyPerRunRapidityH[1]             ->Fill( Form("%d", fRunNum) , 1 );
  } else if ( possibleJPsi.Rapidity() > -3.50 && possibleJPsi.Rapidity() <= -3.25 ) {
    fInvariantMassDistributionRapidityBinsH[2]->Fill(possibleJPsi.Mag());
    fEfficiencyPerRunRapidityH[2]             ->Fill( Form("%d", fRunNum) , 1 );
  } else if ( possibleJPsi.Rapidity() > -3.25 && possibleJPsi.Rapidity() <= -3.00 ) {
    fInvariantMassDistributionRapidityBinsH[3]->Fill(possibleJPsi.Mag());
    fEfficiencyPerRunRapidityH[3]             ->Fill( Form("%d", fRunNum) , 1 );
  } else if ( possibleJPsi.Rapidity() > -3.00 && possibleJPsi.Rapidity() <= -2.75 ) {
    fInvariantMassDistributionRapidityBinsH[4]->Fill(possibleJPsi.Mag());
    fEfficiencyPerRunRapidityH[4]             ->Fill( Form("%d", fRunNum) , 1 );
  } else if ( possibleJPsi.Rapidity() > -2.75 && possibleJPsi.Rapidity() <= -2.50 ) {
    fInvariantMassDistributionRapidityBinsH[5]->Fill(possibleJPsi.Mag());
    fEfficiencyPerRunRapidityH[5]             ->Fill( Form("%d", fRunNum) , 1 );
  }
  fInvariantMassDistributionExtendedH->Fill(possibleJPsi.Mag());

  /* - This is a TH2F histogram filled with DCA against the invariant mass of
     - the dimuon pair. This should be plotted as LEGO or as a heat map...
     - Is this the correct way of doing it?? I have to ask to my supervisor!!
     -
   */
  for(int indexMuonForDCA = 0; indexMuonForDCA < 2; indexMuonForDCA++) {
        fDcaAgainstInvariantMassH->Fill(possibleJPsi.Mag(), track[indexMuonForDCA]->DCA());
  }

  /* - Now we are evaluating the pt of the dimuon pair. Generally speaking,
     - if such a pt is less than 0.25 GeV/c then it fills the coherent
     - component, otherwise the incoherent component. At this point we may fill
     - even the dimuon pt distribution histogram and see if it looks like Kay's.
     -
   */
  Double_t ptOfTheDimuonPair = possibleJPsi.Pt();
  if( ptOfTheDimuonPair < 0.25) {
        fInvariantMassDistributionCoherentH->Fill(possibleJPsi.Mag());
        fInvariantMassDistributionCoherentExtendedH->Fill(possibleJPsi.Mag());
        if (        possibleJPsi.Rapidity() > -4.0  && possibleJPsi.Rapidity() <= -3.75 ) {
          fInvariantMassDistributionCoherentRapidityBinsH[0]->Fill(possibleJPsi.Mag());
        } else if ( possibleJPsi.Rapidity() > -3.75 && possibleJPsi.Rapidity() <= -3.50 ) {
          fInvariantMassDistributionCoherentRapidityBinsH[1]->Fill(possibleJPsi.Mag());
        } else if ( possibleJPsi.Rapidity() > -3.50 && possibleJPsi.Rapidity() <= -3.25 ) {
          fInvariantMassDistributionCoherentRapidityBinsH[2]->Fill(possibleJPsi.Mag());
        } else if ( possibleJPsi.Rapidity() > -3.25 && possibleJPsi.Rapidity() <= -3.00 ) {
          fInvariantMassDistributionCoherentRapidityBinsH[3]->Fill(possibleJPsi.Mag());
        } else if ( possibleJPsi.Rapidity() > -3.00 && possibleJPsi.Rapidity() <= -2.75 ) {
          fInvariantMassDistributionCoherentRapidityBinsH[4]->Fill(possibleJPsi.Mag());
        } else if ( possibleJPsi.Rapidity() > -2.75 && possibleJPsi.Rapidity() <= -2.50 ) {
          fInvariantMassDistributionCoherentRapidityBinsH[5]->Fill(possibleJPsi.Mag());
        }
  } else {
        fInvariantMassDistributionIncoherentH->Fill(possibleJPsi.Mag());
        fInvariantMassDistributionIncoherentExtendedH->Fill(possibleJPsi.Mag());
        if (        possibleJPsi.Rapidity() > -4.0  && possibleJPsi.Rapidity() <= -3.75 ) {
          fInvariantMassDistributionIncoherentRapidityBinsH[0]->Fill(possibleJPsi.Mag());
        } else if ( possibleJPsi.Rapidity() > -3.75 && possibleJPsi.Rapidity() <= -3.50 ) {
          fInvariantMassDistributionIncoherentRapidityBinsH[1]->Fill(possibleJPsi.Mag());
        } else if ( possibleJPsi.Rapidity() > -3.50 && possibleJPsi.Rapidity() <= -3.25 ) {
          fInvariantMassDistributionIncoherentRapidityBinsH[2]->Fill(possibleJPsi.Mag());
        } else if ( possibleJPsi.Rapidity() > -3.25 && possibleJPsi.Rapidity() <= -3.00 ) {
          fInvariantMassDistributionIncoherentRapidityBinsH[3]->Fill(possibleJPsi.Mag());
        } else if ( possibleJPsi.Rapidity() > -3.00 && possibleJPsi.Rapidity() <= -2.75 ) {
          fInvariantMassDistributionIncoherentRapidityBinsH[4]->Fill(possibleJPsi.Mag());
        } else if ( possibleJPsi.Rapidity() > -2.75 && possibleJPsi.Rapidity() <= -2.50 ) {
          fInvariantMassDistributionIncoherentRapidityBinsH[5]->Fill(possibleJPsi.Mag());
        }
  }
  fDimuonPtDistributionH->Fill(ptOfTheDimuonPair);
  if ( (possibleJPsi.Mag() > 2.85) && (possibleJPsi.Mag() < 3.35) ) {
    fTemplatePtDistributionH->Fill(ptOfTheDimuonPair);
    if (        possibleJPsi.Rapidity() > -4.0  && possibleJPsi.Rapidity() <= -3.50 ) {
      fTemplatePtDistributionRapidityH[0]->Fill(possibleJPsi.Mag());
    } else if ( possibleJPsi.Rapidity() > -3.50 && possibleJPsi.Rapidity() <= -3.00 ) {
      fTemplatePtDistributionRapidityH[1]->Fill(possibleJPsi.Mag());
    } else if ( possibleJPsi.Rapidity() > -3.00 && possibleJPsi.Rapidity() <= -2.50 ) {
      fTemplatePtDistributionRapidityH[2]->Fill(possibleJPsi.Mag());
    }
  }

  //_______________________________
  //    SIDEBANDS
  /* -
   * - (LOWER)
   */
  if ( (possibleJPsi.Mag() > 2.4) && (possibleJPsi.Mag() < 2.8) ) {
    fTemplatePtDistributionHLowerSide->Fill(ptOfTheDimuonPair);
    if (        possibleJPsi.Rapidity() > -4.0  && possibleJPsi.Rapidity() <= -3.50 ) {
      fTemplatePtDistributionRapidityHLowerSide[0]->Fill(possibleJPsi.Mag());
    } else if ( possibleJPsi.Rapidity() > -3.50 && possibleJPsi.Rapidity() <= -3.00 ) {
      fTemplatePtDistributionRapidityHLowerSide[1]->Fill(possibleJPsi.Mag());
    } else if ( possibleJPsi.Rapidity() > -3.00 && possibleJPsi.Rapidity() <= -2.50 ) {
      fTemplatePtDistributionRapidityHLowerSide[2]->Fill(possibleJPsi.Mag());
    }
  }
  /* -
   * - (HIGHER)
   */
  if ( (possibleJPsi.Mag() > 4.) && (possibleJPsi.Mag() < 5.5) ) {
    fTemplatePtDistributionHHigherSide->Fill(ptOfTheDimuonPair);
    if (        possibleJPsi.Rapidity() > -4.0  && possibleJPsi.Rapidity() <= -3.50 ) {
      fTemplatePtDistributionRapidityHHigherSide[0]->Fill(possibleJPsi.Mag());
    } else if ( possibleJPsi.Rapidity() > -3.50 && possibleJPsi.Rapidity() <= -3.00 ) {
      fTemplatePtDistributionRapidityHHigherSide[1]->Fill(possibleJPsi.Mag());
    } else if ( possibleJPsi.Rapidity() > -3.00 && possibleJPsi.Rapidity() <= -2.50 ) {
      fTemplatePtDistributionRapidityHHigherSide[2]->Fill(possibleJPsi.Mag());
    }
  }



  /* - Filling the J/Psi's polarization plots.
     -
     - Now we are ordering the muons. The first muon will always be positive.
     - This is useful for the histograms...
   */
  TLorentzVector muonsCopy[2];
  TLorentzVector muonsCopy2[2];
  TLorentzVector possibleJPsiCopy;
  if( chargeOfMuons[0] > 0 ){
    muonsCopy[0]     = muons[0];
    muonsCopy[1]     = muons[1];
  } else if( chargeOfMuons[0] < 0 ){
    muonsCopy[0]     = muons[1];
    muonsCopy[1]     = muons[0];
  }
  muonsCopy2[0]      = muonsCopy[0];
  muonsCopy2[1]      = muonsCopy[1];
  possibleJPsiCopy = possibleJPsi;
  for( Int_t iBoosting = 0; iBoosting < 2; iBoosting++ ) {
    // TLorentzVector boostBack = -(possibleJPsiCopy).BoostVector();
    /* - This snippet has beem taken from the website:
       - http://personalpages.to.infn.it/~puccio/htmldoc/src/AliAODDimuon.cxx.html
       -
     */
    TVector3 beta=(-1./possibleJPsiCopy.E())*possibleJPsiCopy.Vect();
    muonsCopy[iBoosting].Boost( beta );
  }
  Double_t cosThetaMuonsRestFrame[2];
  for( Int_t iAngle = 0; iAngle < 2; iAngle++ ) {
    TVector3 muonsCopyVector        = muonsCopy[iAngle].Vect();
    TVector3 possibleJPsiCopyVector = possibleJPsiCopy.Vect();
    Double_t dotProductMuonJPsi     = muonsCopyVector.Dot(possibleJPsiCopyVector);
    cosThetaMuonsRestFrame[iAngle]  = dotProductMuonJPsi/( muonsCopyVector.Mag() * possibleJPsiCopyVector.Mag() );
  }
  /* - If we are in the J/Psi peak, hence 2.8 < M < 3.3 GeV/c, AND if we are
     - in the coherent regime, so if the Pt < 0.25 GeV/c, we fill the plots.
     -
     - In the following note that the rapidity is well computed, so we are
     - dealing with negative values... -4.0 < Y < -2.5 !!!
     -
   */
  if ( (possibleJPsiCopy.Mag() > 2.85) && (possibleJPsiCopy.Mag() < 3.35) && (possibleJPsiCopy.Pt() < 0.25) ) {
    fAngularDistribOfPositiveMuonRestFrameJPsiH->Fill(cosThetaMuonsRestFrame[0]);
    fAngularDistribOfNegativeMuonRestFrameJPsiH->Fill(cosThetaMuonsRestFrame[1]);
    fCheckHelicityRestFrameJPsiH->Fill( muonsCopy[0].Dot(muonsCopy[1]) );

    /* - New part: filling all possible histograms!
       -
     */
    fCosThetaHelicityFrameJPsiH->Fill( CosThetaHelicityFrame( muonsCopy2[0],
                                                              muonsCopy2[1],
                                                              possibleJPsiCopy
                                                              )
                                                            );
    fCosThetaCollinsSoperFrameJPsiH->Fill( CosThetaCollinsSoper( muonsCopy2[0],
                                                                 muonsCopy2[1],
                                                                 possibleJPsiCopy
                                                                 )
                                                               );
    fPhiHelicityFrameJPsiH->Fill( CosPhiHelicityFrame( muonsCopy2[0],
                                                       muonsCopy2[1],
                                                       possibleJPsiCopy
                                                       )
                                                     );
    fPhiCollinsSoperFrameJPsiH->Fill( CosPhiCollinsSoper( muonsCopy2[0],
                                                          muonsCopy2[1],
                                                          possibleJPsiCopy
                                                          )
                                                        );

    fCosThetaAndPhiHelicityFrameMyBinningH->Fill( CosThetaHelicityFrame( muonsCopy2[0],
                                                                         muonsCopy2[1],
                                                                         possibleJPsiCopy
                                                                         ),
                                                  CosPhiHelicityFrame( muonsCopy2[0],
                                                                       muonsCopy2[1],
                                                                       possibleJPsiCopy
                                                                       )
                                                  );

    fCosThetaAndPhiCsFrameMyBinningH->Fill( CosThetaCollinsSoper( muonsCopy2[0],
                                                                  muonsCopy2[1],
                                                                  possibleJPsiCopy
                                                                  ),
                                            CosPhiCollinsSoper(   muonsCopy2[0],
                                                                  muonsCopy2[1],
                                                                  possibleJPsiCopy
                                                                  )
                                            );

    /* - Now we are filling in terms of rapidity...
       - The easiest way to do so I have envisioned is to simply
       - check everytime if we are below the following threshold
       - in a consecutive sense. This means that if we have not passed
       - the previous check we are at least above it.
       - This readily defines the rapidity bin.
       -
       */
    for(Int_t iRapidityBin = 0; iRapidityBin < 8; iRapidityBin++){
        if( (possibleJPsiCopy.Rapidity() + 4.) < 1.5*((Double_t)iRapidityBin + 1.)/8. ){
          fThetaDistribOfPositiveMuonRestFrameJPsiRapidityBinH[iRapidityBin]->Fill(cosThetaMuonsRestFrame[0]);
          /* - New part: filling all possible histograms!
             -
           */
          fCosThetaHelicityFrameJPsiRapidityBinsH[iRapidityBin]->Fill( CosThetaHelicityFrame( muonsCopy2[0],
                                                                                              muonsCopy2[1],
                                                                                              possibleJPsiCopy
                                                                                              )
                                                                                           );
          fCosThetaCollinsSoperFrameJPsiRapidityBinsH[iRapidityBin]->Fill( CosThetaCollinsSoper( muonsCopy2[0],
                                                                                                 muonsCopy2[1],
                                                                                                 possibleJPsiCopy
                                                                                                 )
                                                                                               );
          fPhiHelicityFrameJPsiRapidityBinsH[iRapidityBin]->Fill( CosPhiHelicityFrame( muonsCopy2[0],
                                                                                       muonsCopy2[1],
                                                                                       possibleJPsiCopy
                                                                                       )
                                                                                     );
          fPhiCollinsSoperFrameJPsiRapidityBinsH[iRapidityBin]->Fill( CosPhiCollinsSoper( muonsCopy2[0],
                                                                                          muonsCopy2[1],
                                                                                          possibleJPsiCopy
                                                                                          )
                                                                                        );
          break;
        }
    }


    /* - NEW:
       -
     */
    Bool_t controlFlag2 = 0;
    Bool_t controlFlag3 = 0;
    if ( possibleJPsiCopy.Pt() < 0.25 ) {
          Double_t CosThetaHelicityFrameValue = CosThetaHelicityFrame( muonsCopy2[0], muonsCopy2[1], possibleJPsiCopy );
          Double_t PhiHelicityFrameValue      =   CosPhiHelicityFrame( muonsCopy2[0], muonsCopy2[1], possibleJPsiCopy );
          for(Int_t iCosThetaBins = 0; iCosThetaBins < 40; iCosThetaBins++) {
            if( controlFlag2 == 1) break;
            if( (CosThetaHelicityFrameValue + 1.) < 2.*((Double_t)iCosThetaBins + 1.)/40. ){
                fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameH[iCosThetaBins]->Fill(possibleJPsiCopy.Mag());
                controlFlag2 = 1;
            }
          }
          for(Int_t iPhiBins = 0; iPhiBins < 50; iPhiBins++) {
            if( controlFlag3 == 1) break;
            if( (PhiHelicityFrameValue + 3.14) < 6.28*((Double_t)iPhiBins + 1.)/50. ){
                fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameH[iPhiBins]->Fill(possibleJPsiCopy.Mag());
                controlFlag3 = 1;
            }
          }

    }

    /* - NEW: analysis with purity of the binning
     * - above 80% in CosTheta.
     * -
     |* - NB: FINAL INCARNATION OF THE ANALYSIS
     */
    if ( (possibleJPsiCopy.Pt() < 0.25) && (possibleJPsiCopy.Mag() < 3.35) && (possibleJPsiCopy.Mag() > 2.85) ) {
          Double_t CosThetaHelicityFrameValue10 = CosThetaHelicityFrame( muonsCopy2[0], muonsCopy2[1], possibleJPsiCopy );
          Double_t PhiHelicityFrameValue10      =   CosPhiHelicityFrame( muonsCopy2[0], muonsCopy2[1], possibleJPsiCopy );
          Double_t TildePhiPositiveCosTheta     = PhiHelicityFrameValue10 - 0.25 * 3.14;
          Double_t TildePhiNegativeCosTheta     = PhiHelicityFrameValue10 - 0.75 * 3.14;

          Double_t CosThetaCollinsSoperValue   = CosThetaCollinsSoper(  muonsCopy2[0], muonsCopy2[1], possibleJPsiCopy );
          Double_t PhiCollinsSoperValue        =   CosPhiCollinsSoper(  muonsCopy2[0], muonsCopy2[1], possibleJPsiCopy );
          Double_t TildePhiPositiveCosThetaCS  = PhiCollinsSoperValue - 0.25 * 3.14;
          Double_t TildePhiNegativeCosThetaCS  = PhiCollinsSoperValue - 0.75 * 3.14;

          /* -
           * - IMPORTANT:
           * -
           * - Comment this when you are not doing the
           * - polarisation check.
           */
          Double_t ReweightingCosThetaHE  = 1.0 / ( 1.0 + CosThetaHelicityFrameValue10 * CosThetaHelicityFrameValue10 );
          Double_t ReweightingCosThetaCS  = 1.0 / ( 1.0 + 0.66 * CosThetaCollinsSoperValue    * CosThetaCollinsSoperValue    );
          Double_t ReweightingCosThetaCS2 = 1.0 / ( 1.0 + CosThetaCollinsSoperValue    * CosThetaCollinsSoperValue    );


          if( TildePhiPositiveCosTheta < 0. ) {
            TildePhiPositiveCosTheta += 2 * TMath::Pi();
          }
          if( TildePhiNegativeCosTheta < 0. ) {
            TildePhiNegativeCosTheta += 2 * TMath::Pi();
          }

          if( TildePhiPositiveCosThetaCS < 0. ) {
            TildePhiPositiveCosThetaCS += 2. * TMath::Pi();
          }
          if( TildePhiNegativeCosThetaCS < 0. ) {
            TildePhiNegativeCosThetaCS += 2. * TMath::Pi();
          }

          /* - HELICITY FRAME ANALYSIS
           * -
           */
          fCosThetaHelicityFrameTwentyfiveBinsH->Fill( CosThetaHelicityFrameValue10 );
          // fCosThetaHelicityFrameTwentyfiveBinsH->Fill( CosThetaHelicityFrameValue10, ReweightingCosThetaHE );
          if ( (track[0]->Pt() > 0.85) && (track[1]->Pt() > 0.85) ) fCosThetaHelicityFrameTwentyfiveBinsTriggerH[0]->Fill( CosThetaHelicityFrameValue10 );
          if ( (track[0]->Pt() > 0.90) && (track[1]->Pt() > 0.90) ) fCosThetaHelicityFrameTwentyfiveBinsTriggerH[1]->Fill( CosThetaHelicityFrameValue10 );
          if ( (track[0]->Pt() > 0.95) && (track[1]->Pt() > 0.95) ) fCosThetaHelicityFrameTwentyfiveBinsTriggerH[2]->Fill( CosThetaHelicityFrameValue10 );
          if ( (track[0]->Pt() > 1.00) && (track[1]->Pt() > 1.00) ) fCosThetaHelicityFrameTwentyfiveBinsTriggerH[3]->Fill( CosThetaHelicityFrameValue10 );
          if ( (track[0]->Pt() > 1.05) && (track[1]->Pt() > 1.05) ) fCosThetaHelicityFrameTwentyfiveBinsTriggerH[4]->Fill( CosThetaHelicityFrameValue10 );
          if ( (track[0]->Pt() > 1.10) && (track[1]->Pt() > 1.10) ) fCosThetaHelicityFrameTwentyfiveBinsTriggerH[5]->Fill( CosThetaHelicityFrameValue10 );
          if ( (track[0]->Pt() > 1.15) && (track[1]->Pt() > 1.15) ) fCosThetaHelicityFrameTwentyfiveBinsTriggerH[6]->Fill( CosThetaHelicityFrameValue10 );

          fPhiHelicityFrameTwentyfiveBinsH     ->Fill( PhiHelicityFrameValue10 );
          if ( (track[0]->Pt() > 0.85) && (track[1]->Pt() > 0.85) ) fPhiHelicityFrameTwentyfiveBinsTriggerH[0]->Fill( PhiHelicityFrameValue10 );
          if ( (track[0]->Pt() > 0.90) && (track[1]->Pt() > 0.90) ) fPhiHelicityFrameTwentyfiveBinsTriggerH[1]->Fill( PhiHelicityFrameValue10 );
          if ( (track[0]->Pt() > 0.95) && (track[1]->Pt() > 0.95) ) fPhiHelicityFrameTwentyfiveBinsTriggerH[2]->Fill( PhiHelicityFrameValue10 );
          if ( (track[0]->Pt() > 1.00) && (track[1]->Pt() > 1.00) ) fPhiHelicityFrameTwentyfiveBinsTriggerH[3]->Fill( PhiHelicityFrameValue10 );
          if ( (track[0]->Pt() > 1.05) && (track[1]->Pt() > 1.05) ) fPhiHelicityFrameTwentyfiveBinsTriggerH[4]->Fill( PhiHelicityFrameValue10 );
          if ( (track[0]->Pt() > 1.10) && (track[1]->Pt() > 1.10) ) fPhiHelicityFrameTwentyfiveBinsTriggerH[5]->Fill( PhiHelicityFrameValue10 );
          if ( (track[0]->Pt() > 1.15) && (track[1]->Pt() > 1.15) ) fPhiHelicityFrameTwentyfiveBinsTriggerH[6]->Fill( PhiHelicityFrameValue10 );

          if( CosThetaHelicityFrameValue10 > 0 ){
            fTildePhiHelicityFrameTwentyfiveBinsH->Fill( TildePhiPositiveCosTheta );
            if ( (track[0]->Pt() > 0.85) && (track[1]->Pt() > 0.85) ) fTildePhiHelicityFrameTwentyfiveBinsTriggerH[0]->Fill( TildePhiPositiveCosTheta );
            if ( (track[0]->Pt() > 0.90) && (track[1]->Pt() > 0.90) ) fTildePhiHelicityFrameTwentyfiveBinsTriggerH[1]->Fill( TildePhiPositiveCosTheta );
            if ( (track[0]->Pt() > 0.95) && (track[1]->Pt() > 0.95) ) fTildePhiHelicityFrameTwentyfiveBinsTriggerH[2]->Fill( TildePhiPositiveCosTheta );
            if ( (track[0]->Pt() > 1.00) && (track[1]->Pt() > 1.00) ) fTildePhiHelicityFrameTwentyfiveBinsTriggerH[3]->Fill( TildePhiPositiveCosTheta );
            if ( (track[0]->Pt() > 1.05) && (track[1]->Pt() > 1.05) ) fTildePhiHelicityFrameTwentyfiveBinsTriggerH[4]->Fill( TildePhiPositiveCosTheta );
            if ( (track[0]->Pt() > 1.10) && (track[1]->Pt() > 1.10) ) fTildePhiHelicityFrameTwentyfiveBinsTriggerH[5]->Fill( TildePhiPositiveCosTheta );
            if ( (track[0]->Pt() > 1.15) && (track[1]->Pt() > 1.15) ) fTildePhiHelicityFrameTwentyfiveBinsTriggerH[6]->Fill( TildePhiPositiveCosTheta );
          } else {
            fTildePhiHelicityFrameTwentyfiveBinsH->Fill( TildePhiNegativeCosTheta );
            if ( (track[0]->Pt() > 0.85) && (track[1]->Pt() > 0.85) ) fTildePhiHelicityFrameTwentyfiveBinsTriggerH[0]->Fill( TildePhiNegativeCosTheta );
            if ( (track[0]->Pt() > 0.90) && (track[1]->Pt() > 0.90) ) fTildePhiHelicityFrameTwentyfiveBinsTriggerH[1]->Fill( TildePhiNegativeCosTheta );
            if ( (track[0]->Pt() > 0.95) && (track[1]->Pt() > 0.95) ) fTildePhiHelicityFrameTwentyfiveBinsTriggerH[2]->Fill( TildePhiNegativeCosTheta );
            if ( (track[0]->Pt() > 1.00) && (track[1]->Pt() > 1.00) ) fTildePhiHelicityFrameTwentyfiveBinsTriggerH[3]->Fill( TildePhiNegativeCosTheta );
            if ( (track[0]->Pt() > 1.05) && (track[1]->Pt() > 1.05) ) fTildePhiHelicityFrameTwentyfiveBinsTriggerH[4]->Fill( TildePhiNegativeCosTheta );
            if ( (track[0]->Pt() > 1.10) && (track[1]->Pt() > 1.10) ) fTildePhiHelicityFrameTwentyfiveBinsTriggerH[5]->Fill( TildePhiNegativeCosTheta );
            if ( (track[0]->Pt() > 1.15) && (track[1]->Pt() > 1.15) ) fTildePhiHelicityFrameTwentyfiveBinsTriggerH[6]->Fill( TildePhiNegativeCosTheta );
          }

          /* - COLLINS-SOPER ANALYSIS
           * -
           */
          fCosThetaCsFrameTwentyfiveBinsH->Fill( CosThetaCollinsSoperValue );
          // fCosThetaCsFrameTwentyfiveBinsH->Fill( CosThetaCollinsSoperValue, ReweightingCosThetaCS  );
          if ( (track[0]->Pt() > 0.85) && (track[1]->Pt() > 0.85) ) fCosThetaCsFrameTwentyfiveBinsTriggerH[0]->Fill( CosThetaCollinsSoperValue );
          if ( (track[0]->Pt() > 0.90) && (track[1]->Pt() > 0.90) ) fCosThetaCsFrameTwentyfiveBinsTriggerH[1]->Fill( CosThetaCollinsSoperValue );
          if ( (track[0]->Pt() > 0.95) && (track[1]->Pt() > 0.95) ) fCosThetaCsFrameTwentyfiveBinsTriggerH[2]->Fill( CosThetaCollinsSoperValue );
          if ( (track[0]->Pt() > 1.00) && (track[1]->Pt() > 1.00) ) fCosThetaCsFrameTwentyfiveBinsTriggerH[3]->Fill( CosThetaCollinsSoperValue );
          if ( (track[0]->Pt() > 1.05) && (track[1]->Pt() > 1.05) ) fCosThetaCsFrameTwentyfiveBinsTriggerH[4]->Fill( CosThetaCollinsSoperValue );
          if ( (track[0]->Pt() > 1.10) && (track[1]->Pt() > 1.10) ) fCosThetaCsFrameTwentyfiveBinsTriggerH[5]->Fill( CosThetaCollinsSoperValue );
          if ( (track[0]->Pt() > 1.15) && (track[1]->Pt() > 1.15) ) fCosThetaCsFrameTwentyfiveBinsTriggerH[6]->Fill( CosThetaCollinsSoperValue );

          fPhiCsFrameTwentyfiveBinsH     ->Fill( PhiCollinsSoperValue );
          if ( (track[0]->Pt() > 0.85) && (track[1]->Pt() > 0.85) ) fPhiCsFrameTwentyfiveBinsTriggerH[0]->Fill( PhiCollinsSoperValue );
          if ( (track[0]->Pt() > 0.90) && (track[1]->Pt() > 0.90) ) fPhiCsFrameTwentyfiveBinsTriggerH[1]->Fill( PhiCollinsSoperValue );
          if ( (track[0]->Pt() > 0.95) && (track[1]->Pt() > 0.95) ) fPhiCsFrameTwentyfiveBinsTriggerH[2]->Fill( PhiCollinsSoperValue );
          if ( (track[0]->Pt() > 1.00) && (track[1]->Pt() > 1.00) ) fPhiCsFrameTwentyfiveBinsTriggerH[3]->Fill( PhiCollinsSoperValue );
          if ( (track[0]->Pt() > 1.05) && (track[1]->Pt() > 1.05) ) fPhiCsFrameTwentyfiveBinsTriggerH[4]->Fill( PhiCollinsSoperValue );
          if ( (track[0]->Pt() > 1.10) && (track[1]->Pt() > 1.10) ) fPhiCsFrameTwentyfiveBinsTriggerH[5]->Fill( PhiCollinsSoperValue );
          if ( (track[0]->Pt() > 1.15) && (track[1]->Pt() > 1.15) ) fPhiCsFrameTwentyfiveBinsTriggerH[6]->Fill( PhiCollinsSoperValue );

          if( CosThetaCollinsSoperValue > 0 ){
            fTildePhiCsFrameTwentyfiveBinsH->Fill( TildePhiPositiveCosThetaCS );
            if ( (track[0]->Pt() > 0.85) && (track[1]->Pt() > 0.85) ) fTildePhiCsFrameTwentyfiveBinsTriggerH[0]->Fill( TildePhiPositiveCosThetaCS );
            if ( (track[0]->Pt() > 0.90) && (track[1]->Pt() > 0.90) ) fTildePhiCsFrameTwentyfiveBinsTriggerH[1]->Fill( TildePhiPositiveCosThetaCS );
            if ( (track[0]->Pt() > 0.95) && (track[1]->Pt() > 0.95) ) fTildePhiCsFrameTwentyfiveBinsTriggerH[2]->Fill( TildePhiPositiveCosThetaCS );
            if ( (track[0]->Pt() > 1.00) && (track[1]->Pt() > 1.00) ) fTildePhiCsFrameTwentyfiveBinsTriggerH[3]->Fill( TildePhiPositiveCosThetaCS );
            if ( (track[0]->Pt() > 1.05) && (track[1]->Pt() > 1.05) ) fTildePhiCsFrameTwentyfiveBinsTriggerH[4]->Fill( TildePhiPositiveCosThetaCS );
            if ( (track[0]->Pt() > 1.10) && (track[1]->Pt() > 1.10) ) fTildePhiCsFrameTwentyfiveBinsTriggerH[5]->Fill( TildePhiPositiveCosThetaCS );
            if ( (track[0]->Pt() > 1.15) && (track[1]->Pt() > 1.15) ) fTildePhiCsFrameTwentyfiveBinsTriggerH[6]->Fill( TildePhiPositiveCosThetaCS );
          } else {
            fTildePhiCsFrameTwentyfiveBinsH->Fill( TildePhiNegativeCosThetaCS );
            if ( (track[0]->Pt() > 0.85) && (track[1]->Pt() > 0.85) ) fTildePhiCsFrameTwentyfiveBinsTriggerH[0]->Fill( TildePhiNegativeCosThetaCS );
            if ( (track[0]->Pt() > 0.90) && (track[1]->Pt() > 0.90) ) fTildePhiCsFrameTwentyfiveBinsTriggerH[1]->Fill( TildePhiNegativeCosThetaCS );
            if ( (track[0]->Pt() > 0.95) && (track[1]->Pt() > 0.95) ) fTildePhiCsFrameTwentyfiveBinsTriggerH[2]->Fill( TildePhiNegativeCosThetaCS );
            if ( (track[0]->Pt() > 1.00) && (track[1]->Pt() > 1.00) ) fTildePhiCsFrameTwentyfiveBinsTriggerH[3]->Fill( TildePhiNegativeCosThetaCS );
            if ( (track[0]->Pt() > 1.05) && (track[1]->Pt() > 1.05) ) fTildePhiCsFrameTwentyfiveBinsTriggerH[4]->Fill( TildePhiNegativeCosThetaCS );
            if ( (track[0]->Pt() > 1.10) && (track[1]->Pt() > 1.10) ) fTildePhiCsFrameTwentyfiveBinsTriggerH[5]->Fill( TildePhiNegativeCosThetaCS );
            if ( (track[0]->Pt() > 1.15) && (track[1]->Pt() > 1.15) ) fTildePhiCsFrameTwentyfiveBinsTriggerH[6]->Fill( TildePhiNegativeCosThetaCS );
          }





          /* =
           * = 2D flat
           */
          fCosThetaAndPhiHelicityFrameMyBinningReweightingH->Fill( CosThetaHelicityFrameValue10, PhiHelicityFrameValue10, ReweightingCosThetaHE  );
          fCosThetaAndPhiCsFrameMyBinningReweightingH->Fill(       CosThetaCollinsSoperValue,    PhiCollinsSoperValue,    ReweightingCosThetaCS2 );


          /* -
           * - 2D trigger
           */
          if ( (track[0]->Pt() > 0.85) && (track[1]->Pt() > 0.85) ) fCosThetaAndPhiHelicityFrameMyBinningTriggerH[0]->Fill( CosThetaHelicityFrameValue10, PhiHelicityFrameValue10 );
          if ( (track[0]->Pt() > 0.90) && (track[1]->Pt() > 0.90) ) fCosThetaAndPhiHelicityFrameMyBinningTriggerH[1]->Fill( CosThetaHelicityFrameValue10, PhiHelicityFrameValue10 );
          if ( (track[0]->Pt() > 0.95) && (track[1]->Pt() > 0.95) ) fCosThetaAndPhiHelicityFrameMyBinningTriggerH[2]->Fill( CosThetaHelicityFrameValue10, PhiHelicityFrameValue10 );
          if ( (track[0]->Pt() > 1.00) && (track[1]->Pt() > 1.00) ) fCosThetaAndPhiHelicityFrameMyBinningTriggerH[3]->Fill( CosThetaHelicityFrameValue10, PhiHelicityFrameValue10 );
          if ( (track[0]->Pt() > 1.05) && (track[1]->Pt() > 1.05) ) fCosThetaAndPhiHelicityFrameMyBinningTriggerH[4]->Fill( CosThetaHelicityFrameValue10, PhiHelicityFrameValue10 );
          if ( (track[0]->Pt() > 1.10) && (track[1]->Pt() > 1.10) ) fCosThetaAndPhiHelicityFrameMyBinningTriggerH[5]->Fill( CosThetaHelicityFrameValue10, PhiHelicityFrameValue10 );
          if ( (track[0]->Pt() > 1.15) && (track[1]->Pt() > 1.15) ) fCosThetaAndPhiHelicityFrameMyBinningTriggerH[6]->Fill( CosThetaHelicityFrameValue10, PhiHelicityFrameValue10 );

          if ( (track[0]->Pt() > 0.85) && (track[1]->Pt() > 0.85) ) fCosThetaAndPhiCsFrameMyBinningTriggerH[0]->Fill( CosThetaCollinsSoperValue,    PhiCollinsSoperValue );
          if ( (track[0]->Pt() > 0.90) && (track[1]->Pt() > 0.90) ) fCosThetaAndPhiCsFrameMyBinningTriggerH[1]->Fill( CosThetaCollinsSoperValue,    PhiCollinsSoperValue );
          if ( (track[0]->Pt() > 0.95) && (track[1]->Pt() > 0.95) ) fCosThetaAndPhiCsFrameMyBinningTriggerH[2]->Fill( CosThetaCollinsSoperValue,    PhiCollinsSoperValue );
          if ( (track[0]->Pt() > 1.00) && (track[1]->Pt() > 1.00) ) fCosThetaAndPhiCsFrameMyBinningTriggerH[3]->Fill( CosThetaCollinsSoperValue,    PhiCollinsSoperValue );
          if ( (track[0]->Pt() > 1.05) && (track[1]->Pt() > 1.05) ) fCosThetaAndPhiCsFrameMyBinningTriggerH[4]->Fill( CosThetaCollinsSoperValue,    PhiCollinsSoperValue );
          if ( (track[0]->Pt() > 1.10) && (track[1]->Pt() > 1.10) ) fCosThetaAndPhiCsFrameMyBinningTriggerH[5]->Fill( CosThetaCollinsSoperValue,    PhiCollinsSoperValue );
          if ( (track[0]->Pt() > 1.15) && (track[1]->Pt() > 1.15) ) fCosThetaAndPhiCsFrameMyBinningTriggerH[6]->Fill( CosThetaCollinsSoperValue,    PhiCollinsSoperValue );



    }




    /* - NEW:
     * - 1D analysis with
     * - my variable binning.
     */
    Bool_t controlFlag5 = 0;
    Bool_t controlFlag6 = 0;
    Double_t MyVariableCosThetaBinning1D[] = { -0.65,  -0.5,  -0.4,  -0.35,  -0.3,
                                               -0.25,  -0.2,  -0.15, -0.125, -0.1,
                                               -0.075, -0.05, -0.025, 0,      0.025,
                                                0.05,   0.075, 0.1,   0.125,  0.15,
                                                0.2,    0.25,  0.3,   0.35,   0.4,
                                                0.5,    0.65 };
    Double_t MyVariablePhiBinning1D[] = { -3.14*1,      -3.14*14/15,  -3.14*13/15,  -3.14*12/15,
                                          -3.14*11/15,  -3.14*10/15,  -3.14*9/15,   -3.14*8/15,
                                          -3.14*7/15,   -3.14*6/15,   -3.14*5/15,   -3.14*4/15,
                                          -3.14*3/15,   -3.14*2/15,   -3.14*1/15,    0,
                                          +3.14*1/15,   +3.14*2/15,   +3.14*3/15,   +3.14*4/15,
                                          +3.14*5/15,   +3.14*6/15,   +3.14*7/15,   +3.14*8/15,
                                          +3.14*9/15,   +3.14*10/15,  +3.14*11/15,  +3.14*12/15,
                                          +3.14*13/15,  +3.14*14/15,  +3.14*1 };
    if ( possibleJPsiCopy.Pt() < 0.25 ) {
          Double_t CosThetaHelicityFrameValue4 = CosThetaHelicityFrame( muonsCopy2[0], muonsCopy2[1], possibleJPsiCopy );
          Double_t PhiHelicityFrameValue4      =   CosPhiHelicityFrame( muonsCopy2[0], muonsCopy2[1], possibleJPsiCopy );
          fCosThetaHelicityFrameMyBinningH         ->Fill( CosThetaHelicityFrameValue4 );
          fCosThetaHelicityFrameMyBinningSmallH    ->Fill( CosThetaHelicityFrameValue4 );
          fCosThetaHelicityFrameMySeventeenBinningH->Fill( CosThetaHelicityFrameValue4 );
          fPhiHelicityFrameMyBinningH              ->Fill( PhiHelicityFrameValue4      );
          for(Int_t iCosThetaBins = 0; iCosThetaBins < 26; iCosThetaBins++) {
            if( controlFlag5 == 1) break;
            if( CosThetaHelicityFrameValue4 < MyVariableCosThetaBinning1D[iCosThetaBins + 1] ){
              fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameMyVariableBinningH[iCosThetaBins]->Fill( possibleJPsiCopy.Mag() );
              controlFlag5 = 1;
            }
          }
          for(Int_t iPhiBins = 0; iPhiBins < 30; iPhiBins++) {
            if( controlFlag6 == 1) break;
            if( PhiHelicityFrameValue4  < MyVariablePhiBinning1D[iPhiBins + 1] ){
              fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameMyVariableBinningH[iPhiBins]->Fill( possibleJPsiCopy.Mag() );
              controlFlag6 = 1;
            }
          }
    }


    /* - Now we are filling in terms of rapidity...
       - The easiest way to do so I have envisioned is to simply
       - check everytime if we are below the following threshold
       - in a consecutive sense. This means that if we have not passed
       - the previous check we are at least above it.
       - This readily defines the rapidity bin.
       -
       -
       - NEW: the following code is for 10 rapidity bins...
       */
    for(Int_t iRapidityBin = 0; iRapidityBin < 10; iRapidityBin++){
        if( (possibleJPsiCopy.Rapidity() + 4.) < 1.5*((Double_t)iRapidityBin + 1.)/10. ){
          fCosThetaHelicityFrameJPsiTenRapidityBinsH[iRapidityBin]->Fill( CosThetaHelicityFrame( muonsCopy2[0],
                                                                                                 muonsCopy2[1],
                                                                                                 possibleJPsiCopy
                                                                                                 )
                                                                                              );
          fCosThetaCollinsSoperFrameJPsiTenRapidityBinsH[iRapidityBin]->Fill( CosThetaCollinsSoper( muonsCopy2[0],
                                                                                                    muonsCopy2[1],
                                                                                                    possibleJPsiCopy
                                                                                                    )
                                                                                               );
          fPhiHelicityFrameJPsiTenRapidityBinsH[iRapidityBin]->Fill( CosPhiHelicityFrame( muonsCopy2[0],
                                                                                          muonsCopy2[1],
                                                                                          possibleJPsiCopy
                                                                                          )
                                                                                     );
          fPhiCollinsSoperFrameJPsiTenRapidityBinsH[iRapidityBin]->Fill( CosPhiCollinsSoper( muonsCopy2[0],
                                                                                             muonsCopy2[1],
                                                                                             possibleJPsiCopy
                                                                                             )
                                                                                        );
          break;
        }
    }


    /* - What we do here is very similar.
       - This time we divide firstly in bins of CosTheta.
       - As many as needed.
       - And then we divide again in terms of Phi.
       - Then we fill.
       - This way we should be able to obtain some kind of map...
       -
     */
    fInvariantMassDistributionBinsOfCosThetaAndPhiHelicityFrameH->Fill( CosThetaHelicityFrame( muonsCopy2[0],
                                                                                               muonsCopy2[1],
                                                                                               possibleJPsiCopy
                                                                                               ),
                                                                        CosPhiHelicityFrame( muonsCopy2[0],
                                                                                             muonsCopy2[1],
                                                                                             possibleJPsiCopy
                                                                                             )
                                                                        );
    fCosThetaAndPhiHelicityFrameInclusivePeopleBinningH->Fill( CosThetaHelicityFrame( muonsCopy2[0],
                                                                                      muonsCopy2[1],
                                                                                      possibleJPsiCopy
                                                                                      ),
                                                               CosPhiHelicityFrame( muonsCopy2[0],
                                                                                    muonsCopy2[1],
                                                                                    possibleJPsiCopy
                                                                                    )
                                                               );
    fInvariantMassDistributionForSignalExtractionHelicityFrameH->Fill( CosThetaHelicityFrame( muonsCopy2[0],
                                                                                              muonsCopy2[1],
                                                                                              possibleJPsiCopy
                                                                                              ),
                                                                       CosPhiHelicityFrame( muonsCopy2[0],
                                                                                            muonsCopy2[1],
                                                                                            possibleJPsiCopy
                                                                                            )
                                                                       );
    fCosThetaHeFrameForSignalExH->Fill( CosThetaHelicityFrame( muonsCopy2[0], muonsCopy2[1], possibleJPsiCopy ) );
    fPhiHeFrameForSignalExH     ->Fill( CosPhiHelicityFrame(   muonsCopy2[0], muonsCopy2[1], possibleJPsiCopy ) );

  }






  /* - NEW:
   * - TEMPLATES needed for signal extraction.
   * -
   */
  Bool_t controlFlag13 = 0;
  Bool_t controlFlag14 = 0;
  Bool_t controlFlag15 = 0;
  Bool_t controlFlag16 = 0;
  Bool_t controlFlag17 = 0;
  Bool_t controlFlag18 = 0;
  if ( possibleJPsiCopy.Pt() < 0.25 ) {
        Double_t CosThetaHelicityFrameValue7 = CosThetaHelicityFrame( muonsCopy2[0], muonsCopy2[1], possibleJPsiCopy );
        Double_t PhiHelicityFrameValue7      =   CosPhiHelicityFrame( muonsCopy2[0], muonsCopy2[1], possibleJPsiCopy );
        // Double_t TildePhiPositiveCosTheta    = CosThetaHelicityFrameValue7 - 0.25 * TMath::Pi();
        // Double_t TildePhiNegativeCosTheta    = CosThetaHelicityFrameValue7 - 0.75 * TMath::Pi();
        Double_t TildePhiPositiveCosTheta    = PhiHelicityFrameValue7 - 0.25 * 3.14;
        Double_t TildePhiNegativeCosTheta    = PhiHelicityFrameValue7 - 0.75 * 3.14;

        Double_t CosThetaCollinsSoperValue   = CosThetaCollinsSoper(  muonsCopy2[0], muonsCopy2[1], possibleJPsiCopy );
        Double_t PhiCollinsSoperValue        =   CosPhiCollinsSoper(  muonsCopy2[0], muonsCopy2[1], possibleJPsiCopy );
        Double_t TildePhiPositiveCosThetaCS  = PhiCollinsSoperValue - 0.25 * 3.14;
        Double_t TildePhiNegativeCosThetaCS  = PhiCollinsSoperValue - 0.75 * 3.14;

        if( TildePhiPositiveCosTheta < 0. ) {
          TildePhiPositiveCosTheta += 2. * TMath::Pi();
        }
        if( TildePhiNegativeCosTheta < 0. ) {
          TildePhiNegativeCosTheta += 2. * TMath::Pi();
        }

        if( TildePhiPositiveCosThetaCS < 0. ) {
          TildePhiPositiveCosThetaCS += 2. * TMath::Pi();
        }
        if( TildePhiNegativeCosThetaCS < 0. ) {
          TildePhiNegativeCosThetaCS += 2. * TMath::Pi();
        }

        /* - HELICITY FRAME ANALYSIS
         * -
         */
        for(Int_t iCosThetaBins = 0; iCosThetaBins < 25; iCosThetaBins++) {
          if( controlFlag13 == 1) break;
          if( (CosThetaHelicityFrameValue7 + 1.) < 2.*((Double_t)iCosThetaBins + 1.)/25. ) {
            fInvariantMassDistributionOnlyCosThetaHeFrameTwentyfiveBinsH[iCosThetaBins]->Fill(possibleJPsiCopy.Mag());
            controlFlag13 = 1;
          }
        }
        for(Int_t iPhiBins = 0; iPhiBins < 25; iPhiBins++) {
          if( controlFlag14 == 1) break;
          if( (PhiHelicityFrameValue7 + 3.14) < 6.28*((Double_t)iPhiBins + 1.)/25. ){
            fInvariantMassDistributionOnlyPhiHeFrameTwentyfiveBinsH[iPhiBins]->Fill(possibleJPsiCopy.Mag());
            controlFlag14 = 1;
          }
        }
        if( (CosThetaHelicityFrameValue7 > 0) || (CosThetaHelicityFrameValue7 == 0) ){
          for(Int_t iTildePhiBins = 0; iTildePhiBins < 25; iTildePhiBins++) {
            if( controlFlag15 == 1) break;
            if( (TildePhiPositiveCosTheta) < 6.28*((Double_t)iTildePhiBins + 1.)/25. ){
              fInvariantMassDistributionOnlyTildePhiHeFrameTwentyfiveBinsH[iTildePhiBins]->Fill(possibleJPsiCopy.Mag());
              controlFlag15 = 1;
            }
          }
        } else if ( CosThetaHelicityFrameValue7 < 0 ){
          for(Int_t iTildePhiBins = 0; iTildePhiBins < 25; iTildePhiBins++) {
            if( controlFlag15 == 1) break;
            if( (TildePhiNegativeCosTheta) < 6.28*((Double_t)iTildePhiBins + 1.)/25. ){
              fInvariantMassDistributionOnlyTildePhiHeFrameTwentyfiveBinsH[iTildePhiBins]->Fill(possibleJPsiCopy.Mag());
              controlFlag15 = 1;
            }
          }
        }

        /* - COLLINS SOPER ANALYSIS
         * -
         */
        for(Int_t iCosThetaBins = 0; iCosThetaBins < 25; iCosThetaBins++) {
          if( controlFlag16 == 1) break;
          if( (CosThetaCollinsSoperValue + 1.) < 2.*((Double_t)iCosThetaBins + 1.)/25. ) {
            fInvariantMassDistributionOnlyCosThetaCsFrameTwentyfiveBinsH[iCosThetaBins]->Fill(possibleJPsiCopy.Mag());
            controlFlag16 = 1;
          }
        }
        for(Int_t iPhiBins = 0; iPhiBins < 25; iPhiBins++) {
          if( controlFlag17 == 1) break;
          if( (PhiCollinsSoperValue + 3.14) < 6.28*((Double_t)iPhiBins + 1.)/25. ){
            fInvariantMassDistributionOnlyPhiCsFrameTwentyfiveBinsH[iPhiBins]->Fill(possibleJPsiCopy.Mag());
            controlFlag17 = 1;
          }
        }
        if( (CosThetaCollinsSoperValue > 0) || (CosThetaCollinsSoperValue == 0) ){
          for(Int_t iTildePhiBins = 0; iTildePhiBins < 25; iTildePhiBins++) {
            if( controlFlag18 == 1) break;
            if( (TildePhiPositiveCosThetaCS) < 6.28*((Double_t)iTildePhiBins + 1.)/25. ){
              fInvariantMassDistributionOnlyTildePhiCsFrameTwentyfiveBinsH[iTildePhiBins]->Fill(possibleJPsiCopy.Mag());
              controlFlag18 = 1;
            }
          }
        } else if ( CosThetaCollinsSoperValue < 0 ){
          for(Int_t iTildePhiBins = 0; iTildePhiBins < 25; iTildePhiBins++) {
            if( controlFlag18 == 1) break;
            if( (TildePhiNegativeCosThetaCS) < 6.28*((Double_t)iTildePhiBins + 1.)/25. ){
              fInvariantMassDistributionOnlyTildePhiCsFrameTwentyfiveBinsH[iTildePhiBins]->Fill(possibleJPsiCopy.Mag());
              controlFlag18 = 1;
            }
          }
        }
  }

  /* - 2D ANALYSIS for POLARISATION.
   * - Comparing HELICITY frame with
   * - COLLINS-SOPER.
   * - Using a really rough binning
   * - to compensate the lack of statistics.
   * -
   */
  Bool_t controlFlag4    = 0;
  Bool_t controlFlag4_CS = 0;
  Double_t MyVariableCosThetaBinning[] = { -0.65, -0.35, -0.15, -0.05,
                                            0.05,  0.15,  0.35,  0.65 };
  Double_t MyVariablePhiBinning[] = { -3.14*1,       -3.14*19*0.05, -3.14*18*0.05, -3.14*17*0.05,
                                      -3.14*13*0.05, -3.14*9*0.05,  -3.14*6*0.05,  -3.14*4*0.05,
                                      -3.14*2*0.05,  -3.14*1*0.05,   0,            +3.14*1*0.05,
                                      +3.14*2*0.05,  +3.14*4*0.05,  +3.14*6*0.05,  +3.14*9*0.05,
                                      +3.14*13*0.05, +3.14*17*0.05, +3.14*18*0.05, +3.14*19*0.05,
                                      +3.14*1 };
  if ( possibleJPsiCopy.Pt() < 0.25 ) {
        Double_t CosThetaHelicityFrameValue3 = CosThetaHelicityFrame( muonsCopy2[0], muonsCopy2[1], possibleJPsiCopy );
        Double_t PhiHelicityFrameValue3      =   CosPhiHelicityFrame( muonsCopy2[0], muonsCopy2[1], possibleJPsiCopy );
        Double_t CosThetaCollinsSoperValue   =  CosThetaCollinsSoper( muonsCopy2[0], muonsCopy2[1], possibleJPsiCopy );
        Double_t PhiCollinsSoperValue        =    CosPhiCollinsSoper( muonsCopy2[0], muonsCopy2[1], possibleJPsiCopy );
        for(Int_t iCosThetaBins = 0; iCosThetaBins < 7; iCosThetaBins++) {
          if( controlFlag4 == 1) break;
          if( CosThetaHelicityFrameValue3 < MyVariableCosThetaBinning[iCosThetaBins + 1] ){
            for(Int_t iPhiBins = 0; iPhiBins < 20; iPhiBins++) {
              if( controlFlag4 == 1) break;
              if( PhiHelicityFrameValue3  < MyVariablePhiBinning[iPhiBins + 1] ){
                  fInvariantMassDistributionForSignalExtractionHelicityFrameMyBinningH[iCosThetaBins][iPhiBins]->Fill(possibleJPsiCopy.Mag());
                  controlFlag4 = 1;
              }
            }
          }
        }
        for(Int_t iCosThetaBins = 0; iCosThetaBins < 7; iCosThetaBins++) {
          if( controlFlag4_CS == 1) break;
          if( CosThetaCollinsSoperValue < MyVariableCosThetaBinning[iCosThetaBins + 1] ){
            for(Int_t iPhiBins = 0; iPhiBins < 20; iPhiBins++) {
              if( controlFlag4_CS == 1) break;
              if( PhiCollinsSoperValue  < MyVariablePhiBinning[iPhiBins + 1] ){
                  fInvariantMassDistributionForSignalExtractionCsFrameMyBinningH[iCosThetaBins][iPhiBins]->Fill(possibleJPsiCopy.Mag());
                  controlFlag4_CS = 1;
              }
            }
          }
        }
  }











  // fVectorCosThetaReconstructed.push_back(cosThetaMuonsRestFrame[0]);
  fCosThetaReconHelicityFrame = 0;
  fCosThetaReconHelicityFrame = cosThetaMuonsRestFrame[0];
  fPhiReconHelicityFrame      = 0;
  fPhiReconHelicityFrame      = CosPhiHelicityFrame( muonsCopy2[0],
                                                     muonsCopy2[1],
                                                     possibleJPsiCopy
                                                     );
  /* - Mind that it could generate segmentation fault without
     - fCounterUPCevent-1, because we are incrementing the counter right after
     - it processes the MC events at Generated level...
     -
     - Comparing old version with vector (unusable on MC due to the huge
     - vector size) to the new easier on resources...
     -
   */
  // if ( fVectorCosThetaGenerated.at(fCounterUPCevent-1) && cosThetaMuonsRestFrame[0] ) {
  //       fBinMigrationHelicityH->Fill( fVectorCosThetaGenerated.at(fCounterUPCevent-1),
  //                                     cosThetaMuonsRestFrame[0]
  //                                   );
  // }
  if ( fCosThetaReconHelicityFrame ) {
        fBinMigrationHelicityH->Fill( fCosThetaGeneratedHelicityFrame,
                                      fCosThetaReconHelicityFrame
                                    );
  }
  if ( fPhiReconHelicityFrame ) {
        fBinMigrationForPhiHelicityH->Fill( fPhiGeneratedHelicityFrame,
                                            fPhiReconHelicityFrame
                                            );
  }




  // post the data
  PostData(1, fOutputList);
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskUPCforwardMC::IsTriggered()
{
  /* - This function IS roughly speaking the trigger for the MC.
     - This code has been greatly inspired by David Horak's work.
     - It returns kTRUE if the trigger has been fired.
     -
   */
  Bool_t is0VBAfired = fBBAFlags > 0;
  Bool_t is0VBCfired = fBBCFlags > 0;
  Bool_t is0UBAfired = fBBAFlagsAD > 0;
  Bool_t is0UBCfired = fBBCFlagsAD > 0;
  // cout << "is0VBAfired = ( fBBAFlags   = " << fBBAFlags   << " ) > 0 => " << is0VBAfired << endl;
  // cout << "is0VBCfired = ( fBBCFlags   = " << fBBCFlags   << " ) > 0 => " << is0VBCfired << endl;
  // cout << "is0UBAfired = ( fBBAFlagsAD = " << fBBAFlagsAD << " ) > 0 => " << is0UBAfired << endl;
  // cout << "is0UBCfired = ( fBBCFlagsAD = " << fBBCFlagsAD << " ) > 0 => " << is0UBCfired << endl;
  if (!is0VBAfired && !is0UBAfired && !is0UBCfired ) return kTRUE;
  // if (!is0VBAfired) return kTRUE;
  else return kFALSE;

}
//_____________________________________________________________________________
void AliAnalysisTaskUPCforwardMC::ProcessMCParticles(AliMCEvent* fMCEventArg)
{
  // Loop over the MC tracks
  for(int iPart = 0; iPart < (fMCEventArg->GetNumberOfTracks()); iPart++) {
    AliAODMCParticle *mcParticle  = (AliAODMCParticle*)fMCEventArg->GetTrack(iPart);
    if (!mcParticle) {
      AliError(Form("Could not receive track %d", iPart));
      continue;
    }
    if (  mcParticle->Charge()    == 0) continue;
    Double_t pT     = mcParticle->Pt();
    Double_t eta    = mcParticle->Eta();
    Double_t pseudo = mcParticle->Y();
    Double_t phi    = mcParticle->Phi();
    Short_t  charge = mcParticle->Charge();
    Int_t    pdg    = mcParticle->PdgCode();
    fMCpdgCodesH->Fill( Form("%d", pdg) , 1 );
    if ( mcParticle->IsPrimary() ) fMCpdgCodesOnlyPrimaryH->Fill( Form("%d", pdg) , 1 );
    fMCphiGeneratedTruthH           ->Fill(phi);
    fMCetaGeneratedTruthH           ->Fill(eta);
    fMCpseudorapidityGeneratedTruthH->Fill(pseudo);
    fMCptGeneratedTruthH            ->Fill(pT);
  }

  /* - This is loop where we record the MC TRUTH stuff BEFORE both
     - the event selection and the track selection.
     - It is nonetheless strange that we can't directly use the
     - J/Psi and that we have to o through such a roundabout
     - road like reconstructing the invariant mass...
     -
   */
  Int_t nGoodMuonsMC = 0;
  TLorentzVector muonsMC[2];
  TLorentzVector possibleJPsiMC;
  Double_t pT[2];
  Double_t eta[2];
  Double_t pseudo[2];
  Double_t phi[2];
  Short_t  charge[2];
  Int_t    pdg[2];
  // Loop over the MC tracks
  for(int iPart = 0; iPart < (fMCEventArg->GetNumberOfTracks()); iPart++) {
    if(nGoodMuonsMC > 2) continue;
    AliAODMCParticle *mcParticle  = (AliAODMCParticle*)fMCEventArg->GetTrack(iPart);
    if (!mcParticle) {
      AliError(Form("Could not receive track %d", iPart));
      continue;
    }
    if (  !mcParticle->IsPrimary()                 ) continue;
    if (   mcParticle->Charge()              == 0  ) continue;
    if (   TMath::Abs(mcParticle->PdgCode()) == 13 ) {
      if ( nGoodMuonsMC < 2 ) {
        // cout << "Ok" << nGoodMuonsMC << endl;
        pT[nGoodMuonsMC]     = mcParticle->Pt();
        eta[nGoodMuonsMC]    = mcParticle->Eta();
        pseudo[nGoodMuonsMC] = mcParticle->Y();
        phi[nGoodMuonsMC]    = mcParticle->Phi();
        charge[nGoodMuonsMC] = mcParticle->Charge();
        pdg[nGoodMuonsMC]    = mcParticle->PdgCode();
        muonsMC[nGoodMuonsMC].SetPtEtaPhiM( pT[nGoodMuonsMC],
                                            eta[nGoodMuonsMC],
                                            phi[nGoodMuonsMC],
                                            TDatabasePDG::Instance()->GetParticle(13)->Mass()
                                           );
        nGoodMuonsMC++;
      }
    }
  }
  // cout << "Tracks in this event:  " << fMCEventArg->GetNumberOfTracks() << endl;
  /* - We need EXACTLY 2 good muons !!!!!
     -
   */
  if( nGoodMuonsMC == 2 ) {
    /* - Unlike-sign muons
       -
     */
    if( charge[0] != charge[1] ) {
      TLorentzVector muonsMCcopy[2];
      muonsMCcopy[0] = muonsMC[0];
      muonsMCcopy[1] = muonsMC[1];
      for( Int_t iMuonsMC = 0; iMuonsMC < 2; iMuonsMC++ ) {
        possibleJPsiMC+=muonsMC[iMuonsMC];
        fMCphiDimuonGeneratedTruthH           ->Fill(phi[iMuonsMC]);
        fMCetaDimuonGeneratedTruthH           ->Fill(eta[iMuonsMC]);
        fMCpseudorapidityDimuonGeneratedTruthH->Fill(pseudo[iMuonsMC]);
        fMCptDimuonGeneratedTruthSingleMuonsH ->Fill(pT[iMuonsMC]);
      }
      fMCinvariantMassDistrJPsiGeneratedTruthH->Fill(possibleJPsiMC.Mag());
      fMCptDimuonGeneratedTruthH->Fill(possibleJPsiMC.Pt());
      if (        possibleJPsiMC.Rapidity() > -4.0  && possibleJPsiMC.Rapidity() <= -3.75 ) {
        fMCEfficiencyPerRunRapidityH[0]->Fill( Form("%d", fRunNum) , 1 );
      } else if ( possibleJPsiMC.Rapidity() > -3.75 && possibleJPsiMC.Rapidity() <= -3.50 ) {
        fMCEfficiencyPerRunRapidityH[1]->Fill( Form("%d", fRunNum) , 1 );
      } else if ( possibleJPsiMC.Rapidity() > -3.50 && possibleJPsiMC.Rapidity() <= -3.25 ) {
        fMCEfficiencyPerRunRapidityH[2]->Fill( Form("%d", fRunNum) , 1 );
      } else if ( possibleJPsiMC.Rapidity() > -3.25 && possibleJPsiMC.Rapidity() <= -3.00 ) {
        fMCEfficiencyPerRunRapidityH[3]->Fill( Form("%d", fRunNum) , 1 );
      } else if ( possibleJPsiMC.Rapidity() > -3.00 && possibleJPsiMC.Rapidity() <= -2.75 ) {
        fMCEfficiencyPerRunRapidityH[4]->Fill( Form("%d", fRunNum) , 1 );
      } else if ( possibleJPsiMC.Rapidity() > -2.75 && possibleJPsiMC.Rapidity() <= -2.50 ) {
        fMCEfficiencyPerRunRapidityH[5]->Fill( Form("%d", fRunNum) , 1 );
      }
      for( Int_t iBoosting = 0; iBoosting < 2; iBoosting++ ) {
        // TLorentzVector boostBack = -(possibleJPsiCopy).BoostVector();
        /* - This snippet has beem taken from the website:
           - http://personalpages.to.infn.it/~puccio/htmldoc/src/AliAODDimuon.cxx.html
           -
         */
        TVector3 betaMC = (-1./possibleJPsiMC.E())*possibleJPsiMC.Vect();
        muonsMC[iBoosting].Boost( betaMC );
      }
      Double_t cosThetaMuonsRestFrameMC[2];
      for( Int_t iAngle = 0; iAngle < 2; iAngle++ ) {
        TVector3 muonsMCVector            = muonsMC[iAngle].Vect();
        TVector3 possibleJPsiMCVector     = possibleJPsiMC.Vect();
        Double_t dotProductMuonJPsiMC     = muonsMCVector.Dot(possibleJPsiMCVector);
        cosThetaMuonsRestFrameMC[iAngle]  = dotProductMuonJPsiMC/( muonsMCVector.Mag() * possibleJPsiMCVector.Mag() );
      }
      if ( (possibleJPsiMC.Mag() > 2.8) && (possibleJPsiMC.Mag() < 3.3) && (possibleJPsiMC.Pt() < 0.25) ) {
          if( charge[0] > 0 ) {
                  /* - This means that [0] is the positive muon while [1]
                     - is the negative muon!
                     -
                   */
                  fMCthetaDistribOfPositiveMuonRestFrameJPsiGeneratedTruthH->Fill(cosThetaMuonsRestFrameMC[0]);
                  fMCthetaDistribOfNegativeMuonRestFrameJPsiGeneratedTruthH->Fill(cosThetaMuonsRestFrameMC[1]);
                  // fVectorCosThetaGenerated.push_back(cosThetaMuonsRestFrameMC[0]);
                  fCosThetaGeneratedHelicityFrame = cosThetaMuonsRestFrameMC[0];
                  fPhiGeneratedHelicityFrame      = CosPhiHelicityFrame( muonsMCcopy[0],
                                                                         muonsMCcopy[1],
                                                                         possibleJPsiMC
                                                                       );
                  /* - New part: filling all possible histograms!
                     -
                   */
                  fMCCosThetaHelicityFrameJPsiH->Fill( CosThetaHelicityFrame( muonsMCcopy[0],
                                                                              muonsMCcopy[1],
                                                                              possibleJPsiMC
                                                                              )
                                                                            );
                  fMCCosThetaCollinsSoperFrameJPsiH->Fill( CosThetaCollinsSoper( muonsMCcopy[0],
                                                                                 muonsMCcopy[1],
                                                                                 possibleJPsiMC
                                                                                 )
                                                                               );
                  fMCPhiHelicityFrameJPsiH->Fill( CosPhiHelicityFrame( muonsMCcopy[0],
                                                                       muonsMCcopy[1],
                                                                       possibleJPsiMC
                                                                       )
                                                                     );
                  fMCPhiCollinsSoperFrameJPsiH->Fill( CosPhiCollinsSoper( muonsMCcopy[0],
                                                                          muonsMCcopy[1],
                                                                          possibleJPsiMC
                                                                          )
                                                                        );
                  fMCCosThetaHeFrameForSignalExH->Fill( CosThetaHelicityFrame( muonsMCcopy[0],
                                                                               muonsMCcopy[1],
                                                                               possibleJPsiMC
                                                                               )
                                                                              );
                  fMCPhiHeFrameForSignalExH->Fill( CosPhiHelicityFrame( muonsMCcopy[0],
                                                                          muonsMCcopy[1],
                                                                          possibleJPsiMC
                                                                          )
                                                                         );
                  fMCCosThetaHelicityFrameMyBinningH->Fill( CosThetaHelicityFrame( muonsMCcopy[0],
                                                                                   muonsMCcopy[1],
                                                                                   possibleJPsiMC
                                                                                   )
                                                                                  );
                  /* - HELICITY FRAME ANALYSIS
                   * -
                   */
                  // fMCCosThetaHelicityFrameTwentyfiveBinsH->Fill( CosThetaHelicityFrame( muonsMCcopy[0],
                  //                                                                       muonsMCcopy[1],
                  //                                                                       possibleJPsiMC
                  //                                                                       )
                  //                                                                      );
                  /* -
                   * - IMPORTANT:
                   * -
                   * - Comment this when you are not doing the
                   * - polarisation check.
                   */
                  Double_t CosThetaHeForTrial   = CosThetaHelicityFrame(muonsMCcopy[0],muonsMCcopy[1],possibleJPsiMC);
                  Double_t CosThetaCsForTrial   = CosThetaCollinsSoper( muonsMCcopy[0],muonsMCcopy[1],possibleJPsiMC);
                  // Double_t ReweightedCosThetaHE = CosThetaHeForTrial / ( 1 + CosThetaHeForTrial * CosThetaHeForTrial );
                  // Double_t ReweightedCosThetaCS = CosThetaCsForTrial / ( 1 + CosThetaCsForTrial * CosThetaCsForTrial );
                  Double_t ReweightedCosThetaHE  = 1.0 / ( 1.0 + CosThetaHeForTrial * CosThetaHeForTrial );
                  Double_t ReweightedCosThetaCS  = 1.0 / ( 1.0 + 0.66 * CosThetaCsForTrial * CosThetaCsForTrial );
                  Double_t ReweightedCosThetaCS2 = 1.0 / ( 1.0 + CosThetaCsForTrial * CosThetaCsForTrial );
                  // fMCCosThetaHelicityFrameTwentyfiveBinsH->Fill( CosThetaHeForTrial*ReweightedCosThetaHE );
                  fMCCosThetaHelicityFrameTwentyfiveBinsH->Fill( CosThetaHeForTrial );
                  // fMCCosThetaHelicityFrameTwentyfiveBinsH->Fill( CosThetaHeForTrial, ReweightedCosThetaHE );
                  fMCPhiHelicityFrameTwentyfiveBinsH->Fill( CosPhiHelicityFrame( muonsMCcopy[0],
                                                                                 muonsMCcopy[1],
                                                                                 possibleJPsiMC
                                                                                 )
                                                                                );
                  if( CosThetaHelicityFrame(muonsMCcopy[0],muonsMCcopy[1],possibleJPsiMC) > 0 ){
                    Double_t PhiHelicityFrameValueTruth  = CosPhiHelicityFrame( muonsMCcopy[0],muonsMCcopy[1],possibleJPsiMC );
                    Double_t TildePhiPositiveCosTheta    = PhiHelicityFrameValueTruth - 0.25 * 3.14;
                    Double_t TildePhiNegativeCosTheta    = PhiHelicityFrameValueTruth - 0.75 * 3.14;
                    if( TildePhiPositiveCosTheta < 0. ) {
                      TildePhiPositiveCosTheta += 2. * TMath::Pi();
                    }
                    if( TildePhiNegativeCosTheta < 0. ) {
                      TildePhiNegativeCosTheta += 2. * TMath::Pi();
                    }
                    fMCTildePhiHelicityFrameTwentyfiveBinsH->Fill( TildePhiPositiveCosTheta );
                  } else {
                    Double_t PhiHelicityFrameValueTruth  = CosPhiHelicityFrame( muonsMCcopy[0],muonsMCcopy[1],possibleJPsiMC );
                    Double_t TildePhiPositiveCosTheta    = PhiHelicityFrameValueTruth - 0.25 * 3.14;
                    Double_t TildePhiNegativeCosTheta    = PhiHelicityFrameValueTruth - 0.75 * 3.14;
                    if( TildePhiPositiveCosTheta < 0. ) {
                      TildePhiPositiveCosTheta += 2. * TMath::Pi();
                    }
                    if( TildePhiNegativeCosTheta < 0. ) {
                      TildePhiNegativeCosTheta += 2. * TMath::Pi();
                    }
                    fMCTildePhiHelicityFrameTwentyfiveBinsH->Fill( TildePhiNegativeCosTheta );
                  }


                  /* =
                   * = 2D flat
                   */
                  fMCCosThetaAndPhiHelicityFrameMyBinningReweightingH->Fill( CosThetaHeForTrial, CosPhiCollinsSoper( muonsMCcopy[0],muonsMCcopy[1],possibleJPsiMC ), ReweightedCosThetaHE  );
                  fMCCosThetaAndPhiCsFrameMyBinningReweightingH->Fill(       CosThetaCsForTrial, CosPhiCollinsSoper( muonsMCcopy[0],muonsMCcopy[1],possibleJPsiMC ), ReweightedCosThetaCS2 );
                  fMCCosThetaAndPhiHelicityFrameReweightingH->Fill(          CosThetaHeForTrial, CosPhiCollinsSoper( muonsMCcopy[0],muonsMCcopy[1],possibleJPsiMC ), ReweightedCosThetaHE  );
                  fMCCosThetaAndPhiCsFrameReweightingH->Fill(                CosThetaCsForTrial, CosPhiCollinsSoper( muonsMCcopy[0],muonsMCcopy[1],possibleJPsiMC ), ReweightedCosThetaCS2 );


                  /* - COLLINS-SOPER ANALYSIS
                   * -
                   */
                  // fMCCosThetaCsFrameTwentyfiveBinsH->Fill( CosThetaCollinsSoper( muonsMCcopy[0],
                  //                                                                muonsMCcopy[1],
                  //                                                                possibleJPsiMC
                  //                                                                )
                  //                                                               );
                  fMCCosThetaCsFrameTwentyfiveBinsH->Fill( CosThetaCsForTrial );
                  // fMCCosThetaCsFrameTwentyfiveBinsH->Fill( CosThetaCsForTrial, ReweightedCosThetaCS );
                  fMCCosThetaHeVsCsH               ->Fill( CosThetaHeForTrial, CosThetaCsForTrial   );
                  fMCCosThetaHeVsCsFlatH           ->Fill( CosThetaHeForTrial, CosThetaCsForTrial, ReweightedCosThetaHE*ReweightedCosThetaCS );
                  // fMCCosThetaCsFrameTwentyfiveBinsH->Fill( CosThetaCsForTrial*ReweightedCosThetaCS );
                  fMCPhiCsFrameTwentyfiveBinsH->Fill( CosPhiCollinsSoper( muonsMCcopy[0],
                                                                          muonsMCcopy[1],
                                                                          possibleJPsiMC
                                                                          )
                                                                         );
                  if( CosThetaCollinsSoper(muonsMCcopy[0],muonsMCcopy[1],possibleJPsiMC) > 0 ){
                    Double_t PhiCsFrameValueTruth          = CosPhiCollinsSoper( muonsMCcopy[0],muonsMCcopy[1],possibleJPsiMC );
                    Double_t TildePhiPositiveCosThetaCS    = PhiCsFrameValueTruth - 0.25 * 3.14;
                    Double_t TildePhiNegativeCosThetaCS    = PhiCsFrameValueTruth - 0.75 * 3.14;
                    if( TildePhiPositiveCosThetaCS < 0. ) {
                      TildePhiPositiveCosThetaCS += 2. * TMath::Pi();
                    }
                    if( TildePhiNegativeCosThetaCS < 0. ) {
                      TildePhiNegativeCosThetaCS += 2. * TMath::Pi();
                    }
                    fMCTildePhiCsFrameTwentyfiveBinsH->Fill( TildePhiPositiveCosThetaCS );
                  } else {
                    Double_t PhiCsFrameValueTruth  = CosPhiCollinsSoper( muonsMCcopy[0],muonsMCcopy[1],possibleJPsiMC );
                    Double_t TildePhiPositiveCosThetaCS    = PhiCsFrameValueTruth - 0.25 * 3.14;
                    Double_t TildePhiNegativeCosThetaCS    = PhiCsFrameValueTruth - 0.75 * 3.14;
                    if( TildePhiPositiveCosThetaCS < 0. ) {
                      TildePhiPositiveCosThetaCS += 2. * TMath::Pi();
                    }
                    if( TildePhiNegativeCosThetaCS < 0. ) {
                      TildePhiNegativeCosThetaCS += 2. * TMath::Pi();
                    }
                    fMCTildePhiCsFrameTwentyfiveBinsH->Fill( TildePhiNegativeCosThetaCS );
                  }
                  //_______________________________
                  fMCCosThetaHelicityFrameMyBinningSmallH->Fill( CosThetaHelicityFrame(  muonsMCcopy[0],
                                                                                         muonsMCcopy[1],
                                                                                         possibleJPsiMC
                                                                                         )
                                                                                        );
                  fMCCosThetaHelicityFrameMySeventeenBinningH->Fill( CosThetaHelicityFrame(  muonsMCcopy[0],
                                                                                             muonsMCcopy[1],
                                                                                             possibleJPsiMC
                                                                                             )
                                                                                            );
                  fMCPhiHelicityFrameMyBinningH->Fill( CosPhiHelicityFrame( muonsMCcopy[0],
                                                                            muonsMCcopy[1],
                                                                            possibleJPsiMC
                                                                            )
                                                                           );
                  fMCInvariantMassDistributionForSignalExtractionHelicityFrameH->Fill( CosThetaHelicityFrame( muonsMCcopy[0],
                                                                                                              muonsMCcopy[1],
                                                                                                              possibleJPsiMC
                                                                                                              ),
                                                                                       CosPhiHelicityFrame( muonsMCcopy[0],
                                                                                                            muonsMCcopy[1],
                                                                                                            possibleJPsiMC
                                                                                                            )
                                                                                       );

                  /* - Now we are filling in terms of rapidity...
                     - The easiest way to do so I have envisioned is to simply
                     - check everytime if we are below the following threshold
                     - in a consecutive sense. This means that if we have not passed
                     - the previous check we are at least above it.
                     - This readily defines the rapidity bin.
                     -
                     */
                  for(Int_t iRapidityBin = 0; iRapidityBin < 8; iRapidityBin++){
                      if( (possibleJPsiMC.Rapidity() + 4.) < 1.5*((Double_t)iRapidityBin + 1.)/8. ){
                        fMCthetaDistribOfPositiveMuonRestFrameJPsiGeneratedTruthRapidityBinH[iRapidityBin]->Fill(cosThetaMuonsRestFrameMC[0]);
                        /* - New part: filling all possible histograms!
                           -
                         */
                        fMCCosThetaHelicityFrameJPsiRapidityBinsH[iRapidityBin]->Fill( CosThetaHelicityFrame( muonsMCcopy[0],
                                                                                                              muonsMCcopy[1],
                                                                                                              possibleJPsiMC
                                                                                                              )
                                                                                                            );
                        fMCCosThetaCollinsSoperFrameJPsiRapidityBinsH[iRapidityBin]->Fill( CosThetaCollinsSoper( muonsMCcopy[0],
                                                                                                                 muonsMCcopy[1],
                                                                                                                 possibleJPsiMC
                                                                                                                 )
                                                                                                               );
                        fMCPhiHelicityFrameJPsiRapidityBinsH[iRapidityBin]->Fill( CosPhiHelicityFrame( muonsMCcopy[0],
                                                                                                       muonsMCcopy[1],
                                                                                                       possibleJPsiMC
                                                                                                       )
                                                                                                     );
                        fMCPhiCollinsSoperFrameJPsiRapidityBinsH[iRapidityBin]->Fill( CosPhiCollinsSoper( muonsMCcopy[0],
                                                                                                          muonsMCcopy[1],
                                                                                                          possibleJPsiMC
                                                                                                          )
                                                                                                        );
                        break;
                      }
                  }

                  /* - Filling for 10 rapidity bins.
                     -
                   */
                  for(Int_t iRapidityBin = 0; iRapidityBin < 10; iRapidityBin++){
                      if( (possibleJPsiMC.Rapidity() + 4.) < 1.5*((Double_t)iRapidityBin + 1.)/10. ){
                        fMCCosThetaHelicityFrameJPsiTenRapidityBinsH[iRapidityBin]->Fill( CosThetaHelicityFrame( muonsMCcopy[0],
                                                                                                                 muonsMCcopy[1],
                                                                                                                 possibleJPsiMC
                                                                                                                 )
                                                                                                            );
                        fMCCosThetaCollinsSoperFrameJPsiTenRapidityBinsH[iRapidityBin]->Fill( CosThetaCollinsSoper( muonsMCcopy[0],
                                                                                                                    muonsMCcopy[1],
                                                                                                                    possibleJPsiMC
                                                                                                                    )
                                                                                                               );
                        fMCPhiHelicityFrameJPsiTenRapidityBinsH[iRapidityBin]->Fill( CosPhiHelicityFrame( muonsMCcopy[0],
                                                                                                          muonsMCcopy[1],
                                                                                                          possibleJPsiMC
                                                                                                          )
                                                                                                     );
                        fMCPhiCollinsSoperFrameJPsiTenRapidityBinsH[iRapidityBin]->Fill( CosPhiCollinsSoper( muonsMCcopy[0],
                                                                                                             muonsMCcopy[1],
                                                                                                             possibleJPsiMC
                                                                                                             )
                                                                                                        );
                        break;
                      }
                  }
                  /* - What we do here is very similar.
                     - This time we divide firstly in bins of CosTheta.
                     - As many as needed.
                     - And then we divide again in terms of Phi.
                     - Then we fill.
                     - This way we should be able to obtain some kind of map...
                     -
                   */
                  fMCInvariantMassDistributionBinsOfCosThetaAndPhiHelicityFrameH->Fill( CosThetaHelicityFrame( muonsMCcopy[0],
                                                                                                               muonsMCcopy[1],
                                                                                                               possibleJPsiMC
                                                                                                               ),
                                                                                        CosPhiHelicityFrame( muonsMCcopy[0],
                                                                                                             muonsMCcopy[1],
                                                                                                             possibleJPsiMC
                                                                                                             )
                                                                                        );
                  fMCCosThetaAndPhiHelicityFrameInclusivePeopleBinningH->Fill( CosThetaHelicityFrame( muonsMCcopy[0],
                                                                                                      muonsMCcopy[1],
                                                                                                      possibleJPsiMC
                                                                                                      ),
                                                                               CosPhiHelicityFrame( muonsMCcopy[0],
                                                                                                    muonsMCcopy[1],
                                                                                                    possibleJPsiMC
                                                                                                    )
                                                                               );
                  fMCCosThetaAndPhiHelicityFrameMyBinningH->Fill( CosThetaHelicityFrame( muonsMCcopy[0],
                                                                                       muonsMCcopy[1],
                                                                                       possibleJPsiMC
                                                                                       ),
                                                                CosPhiHelicityFrame( muonsMCcopy[0],
                                                                                     muonsMCcopy[1],
                                                                                     possibleJPsiMC
                                                                                     )
                                                                );
                  fMCCosThetaAndPhiCsFrameMyBinningH->Fill( CosThetaCollinsSoper( muonsMCcopy[0],
                                                                                  muonsMCcopy[1],
                                                                                  possibleJPsiMC
                                                                                  ),
                                                            CosPhiCollinsSoper(   muonsMCcopy[0],
                                                                                  muonsMCcopy[1],
                                                                                  possibleJPsiMC
                                                                                  )
                                                            );
          } else  {
                  fMCthetaDistribOfNegativeMuonRestFrameJPsiGeneratedTruthH->Fill(cosThetaMuonsRestFrameMC[0]);
                  fMCthetaDistribOfPositiveMuonRestFrameJPsiGeneratedTruthH->Fill(cosThetaMuonsRestFrameMC[1]);
                  // fVectorCosThetaGenerated.push_back(cosThetaMuonsRestFrameMC[1]);
                  fCosThetaGeneratedHelicityFrame = cosThetaMuonsRestFrameMC[1];
                  fPhiGeneratedHelicityFrame      = CosPhiHelicityFrame( muonsMCcopy[1],
                                                                         muonsMCcopy[0],
                                                                         possibleJPsiMC
                                                                         );
                  /* - New part: filling all possible histograms!
                     -
                   */
                  fMCCosThetaHelicityFrameJPsiH->Fill( CosThetaHelicityFrame( muonsMCcopy[1],
                                                                              muonsMCcopy[0],
                                                                              possibleJPsiMC
                                                                              )
                                                                            );
                  fMCCosThetaCollinsSoperFrameJPsiH->Fill( CosThetaCollinsSoper( muonsMCcopy[1],
                                                                                 muonsMCcopy[0],
                                                                                 possibleJPsiMC
                                                                                 )
                                                                               );
                  fMCPhiHelicityFrameJPsiH->Fill( CosPhiHelicityFrame( muonsMCcopy[1],
                                                                       muonsMCcopy[0],
                                                                       possibleJPsiMC
                                                                       )
                                                                     );
                  fMCPhiCollinsSoperFrameJPsiH->Fill( CosPhiCollinsSoper( muonsMCcopy[1],
                                                                          muonsMCcopy[0],
                                                                          possibleJPsiMC
                                                                          )
                                                                        );
                  fMCCosThetaHeFrameForSignalExH->Fill( CosThetaHelicityFrame( muonsMCcopy[1],
                                                                               muonsMCcopy[0],
                                                                               possibleJPsiMC
                                                                               )
                                                                              );
                  fMCPhiHeFrameForSignalExH->Fill( CosPhiHelicityFrame(   muonsMCcopy[1],
                                                                          muonsMCcopy[0],
                                                                          possibleJPsiMC
                                                                          )
                                                                         );
                  fMCCosThetaHelicityFrameMyBinningH->Fill( CosThetaHelicityFrame( muonsMCcopy[1],
                                                                                   muonsMCcopy[0],
                                                                                   possibleJPsiMC
                                                                                   )
                                                                                  );
                  /* - HELICITY FRAME ANALYSIS
                   * -
                   */
                  /* -
                   * - IMPORTANT:
                   * -
                   * - Comment this when you are not doing the
                   * - polarisation check.
                   */
                  Double_t CosThetaHeForTrial   = CosThetaHelicityFrame(muonsMCcopy[1],muonsMCcopy[0],possibleJPsiMC);
                  Double_t CosThetaCsForTrial   = CosThetaCollinsSoper( muonsMCcopy[1],muonsMCcopy[0],possibleJPsiMC);
                  // Double_t ReweightedCosThetaHE = CosThetaHeForTrial / ( 1 + CosThetaHeForTrial * CosThetaHeForTrial );
                  // Double_t ReweightedCosThetaCS = CosThetaCsForTrial / ( 1 + CosThetaCsForTrial * CosThetaCsForTrial );
                  Double_t ReweightedCosThetaHE  = 1.0 / ( 1.0 + CosThetaHeForTrial * CosThetaHeForTrial );
                  Double_t ReweightedCosThetaCS  = 1.0 / ( 1.0 + 0.66 * CosThetaCsForTrial * CosThetaCsForTrial );
                  Double_t ReweightedCosThetaCS2 = 1.0 / ( 1.0 + CosThetaCsForTrial * CosThetaCsForTrial );
                  // fMCCosThetaHelicityFrameTwentyfiveBinsH->Fill( CosThetaHelicityFrame( muonsMCcopy[1],
                  //                                                                       muonsMCcopy[0],
                  //                                                                       possibleJPsiMC
                  //                                                                       )
                  //                                                                      );
                  fMCCosThetaHelicityFrameTwentyfiveBinsH->Fill( CosThetaHeForTrial );
                  // fMCCosThetaHelicityFrameTwentyfiveBinsH->Fill( CosThetaHeForTrial, ReweightedCosThetaHE );
                  // fMCCosThetaHelicityFrameTwentyfiveBinsH->Fill( CosThetaHeForTrial*ReweightedCosThetaHE );
                  fMCPhiHelicityFrameTwentyfiveBinsH->Fill( CosPhiHelicityFrame( muonsMCcopy[1],
                                                                                 muonsMCcopy[0],
                                                                                 possibleJPsiMC
                                                                                 )
                                                                                );
                  if( CosThetaHelicityFrame(muonsMCcopy[1],muonsMCcopy[0],possibleJPsiMC) > 0 ){
                    Double_t PhiHelicityFrameValueTruth  = CosPhiHelicityFrame( muonsMCcopy[1],muonsMCcopy[0],possibleJPsiMC );
                    Double_t TildePhiPositiveCosTheta    = PhiHelicityFrameValueTruth - 0.25 * 3.14;
                    Double_t TildePhiNegativeCosTheta    = PhiHelicityFrameValueTruth - 0.75 * 3.14;
                    if( TildePhiPositiveCosTheta < 0. ) {
                      TildePhiPositiveCosTheta += 2. * TMath::Pi();
                    }
                    if( TildePhiNegativeCosTheta < 0. ) {
                      TildePhiNegativeCosTheta += 2. * TMath::Pi();
                    }
                    fMCTildePhiHelicityFrameTwentyfiveBinsH->Fill( TildePhiPositiveCosTheta );
                  } else {
                    Double_t PhiHelicityFrameValueTruth  = CosPhiHelicityFrame( muonsMCcopy[1],muonsMCcopy[0],possibleJPsiMC );
                    Double_t TildePhiPositiveCosTheta    = PhiHelicityFrameValueTruth - 0.25 * 3.14;
                    Double_t TildePhiNegativeCosTheta    = PhiHelicityFrameValueTruth - 0.75 * 3.14;
                    if( TildePhiPositiveCosTheta < 0. ) {
                      TildePhiPositiveCosTheta += 2. * TMath::Pi();
                    }
                    if( TildePhiNegativeCosTheta < 0. ) {
                      TildePhiNegativeCosTheta += 2. * TMath::Pi();
                    }
                    fMCTildePhiHelicityFrameTwentyfiveBinsH->Fill( TildePhiNegativeCosTheta );
                  }



                  /* =
                   * = 2D flat
                   */
                  fMCCosThetaAndPhiHelicityFrameMyBinningReweightingH->Fill( CosThetaHeForTrial, CosPhiCollinsSoper( muonsMCcopy[1],muonsMCcopy[0],possibleJPsiMC ), ReweightedCosThetaHE  );
                  fMCCosThetaAndPhiCsFrameMyBinningReweightingH->Fill(       CosThetaCsForTrial, CosPhiCollinsSoper( muonsMCcopy[1],muonsMCcopy[0],possibleJPsiMC ), ReweightedCosThetaCS2 );
                  fMCCosThetaAndPhiHelicityFrameReweightingH->Fill(          CosThetaHeForTrial, CosPhiCollinsSoper( muonsMCcopy[1],muonsMCcopy[0],possibleJPsiMC ), ReweightedCosThetaHE  );
                  fMCCosThetaAndPhiCsFrameReweightingH->Fill(                CosThetaCsForTrial, CosPhiCollinsSoper( muonsMCcopy[1],muonsMCcopy[0],possibleJPsiMC ), ReweightedCosThetaCS2 );




                  /* - COLLINS-SOPER ANALYSIS
                   * -
                   */
                  // fMCCosThetaCsFrameTwentyfiveBinsH->Fill( CosThetaCollinsSoper( muonsMCcopy[1],
                  //                                                                muonsMCcopy[0],
                  //                                                                possibleJPsiMC
                  //                                                                )
                  //                                                               );
                  fMCCosThetaCsFrameTwentyfiveBinsH->Fill( CosThetaCsForTrial );
                  // fMCCosThetaCsFrameTwentyfiveBinsH->Fill( CosThetaCsForTrial, ReweightedCosThetaCS );
                  fMCCosThetaHeVsCsH               ->Fill( CosThetaHeForTrial, CosThetaCsForTrial   );
                  fMCCosThetaHeVsCsFlatH           ->Fill( CosThetaHeForTrial, CosThetaCsForTrial, ReweightedCosThetaHE*ReweightedCosThetaCS );
                  // fMCCosThetaCsFrameTwentyfiveBinsH->Fill( CosThetaCsForTrial*ReweightedCosThetaCS );
                  fMCPhiCsFrameTwentyfiveBinsH->Fill( CosPhiCollinsSoper( muonsMCcopy[1],
                                                                          muonsMCcopy[0],
                                                                          possibleJPsiMC
                                                                          )
                                                                         );
                  if( CosThetaCollinsSoper(muonsMCcopy[1],muonsMCcopy[0],possibleJPsiMC) > 0 ){
                      Double_t PhiCsFrameValueTruth          = CosPhiCollinsSoper( muonsMCcopy[1],muonsMCcopy[0],possibleJPsiMC );
                      Double_t TildePhiPositiveCosThetaCS    = PhiCsFrameValueTruth - 0.25 * 3.14;
                      Double_t TildePhiNegativeCosThetaCS    = PhiCsFrameValueTruth - 0.75 * 3.14;
                      if( TildePhiPositiveCosThetaCS < 0. ) {
                        TildePhiPositiveCosThetaCS += 2. * TMath::Pi();
                      }
                      if( TildePhiNegativeCosThetaCS < 0. ) {
                        TildePhiNegativeCosThetaCS += 2. * TMath::Pi();
                      }
                      fMCTildePhiCsFrameTwentyfiveBinsH->Fill( TildePhiPositiveCosThetaCS );
                  } else {
                      Double_t PhiCsFrameValueTruth  = CosPhiCollinsSoper( muonsMCcopy[1],muonsMCcopy[0],possibleJPsiMC );
                      Double_t TildePhiPositiveCosThetaCS    = PhiCsFrameValueTruth - 0.25 * 3.14;
                      Double_t TildePhiNegativeCosThetaCS    = PhiCsFrameValueTruth - 0.75 * 3.14;
                      if( TildePhiPositiveCosThetaCS < 0. ) {
                        TildePhiPositiveCosThetaCS += 2. * TMath::Pi();
                      }
                      if( TildePhiNegativeCosThetaCS < 0. ) {
                        TildePhiNegativeCosThetaCS += 2. * TMath::Pi();
                      }
                      fMCTildePhiCsFrameTwentyfiveBinsH->Fill( TildePhiNegativeCosThetaCS );
                  }
                  //_______________________________
                  fMCCosThetaHelicityFrameMyBinningSmallH->Fill( CosThetaHelicityFrame(  muonsMCcopy[1],
                                                                                         muonsMCcopy[0],
                                                                                         possibleJPsiMC
                                                                                         )
                                                                                        );
                  fMCCosThetaHelicityFrameMySeventeenBinningH->Fill( CosThetaHelicityFrame(  muonsMCcopy[1],
                                                                                             muonsMCcopy[0],
                                                                                             possibleJPsiMC
                                                                                             )
                                                                                            );
                  fMCPhiHelicityFrameMyBinningH->Fill( CosPhiHelicityFrame( muonsMCcopy[1],
                                                                            muonsMCcopy[0],
                                                                            possibleJPsiMC
                                                                            )
                                                                           );
                  fMCInvariantMassDistributionForSignalExtractionHelicityFrameH->Fill( CosThetaHelicityFrame( muonsMCcopy[1],
                                                                                                              muonsMCcopy[0],
                                                                                                              possibleJPsiMC
                                                                                                              ),
                                                                                       CosPhiHelicityFrame( muonsMCcopy[1],
                                                                                                            muonsMCcopy[0],
                                                                                                            possibleJPsiMC
                                                                                                            )
                                                                                       );

                  /* - Now we are filling in terms of rapidity...
                     - The easiest way to do so I have envisioned is to simply
                     - check everytime if we are below the following threshold
                     - in a consecutive sense. This means that if we have not passed
                     - the previous check we are at least above it.
                     - This readily defines the rapidity bin.
                     -
                     */
                  for(Int_t iRapidityBin = 0; iRapidityBin < 8; iRapidityBin++){
                      if( (possibleJPsiMC.Rapidity() + 4.) < 1.5*((Double_t)iRapidityBin + 1.)/8. ){
                        fMCthetaDistribOfPositiveMuonRestFrameJPsiGeneratedTruthRapidityBinH[iRapidityBin]->Fill(cosThetaMuonsRestFrameMC[1]);
                        /* - New part: filling all possible histograms!
                           -
                         */
                        fMCCosThetaHelicityFrameJPsiRapidityBinsH[iRapidityBin]->Fill( CosThetaHelicityFrame( muonsMCcopy[1],
                                                                                                              muonsMCcopy[0],
                                                                                                              possibleJPsiMC
                                                                                                              )
                                                                                                            );
                        fMCCosThetaCollinsSoperFrameJPsiRapidityBinsH[iRapidityBin]->Fill( CosThetaCollinsSoper( muonsMCcopy[1],
                                                                                                                 muonsMCcopy[0],
                                                                                                                 possibleJPsiMC
                                                                                                                 )
                                                                                                               );
                        fMCPhiHelicityFrameJPsiRapidityBinsH[iRapidityBin]->Fill( CosPhiHelicityFrame( muonsMCcopy[1],
                                                                                                       muonsMCcopy[0],
                                                                                                       possibleJPsiMC
                                                                                                       )
                                                                                                     );
                        fMCPhiCollinsSoperFrameJPsiRapidityBinsH[iRapidityBin]->Fill( CosPhiCollinsSoper( muonsMCcopy[1],
                                                                                                          muonsMCcopy[0],
                                                                                                          possibleJPsiMC
                                                                                                          )
                                                                                                        );
                        break;
                      }
                  }

                  /* - Filling for 10 rapidity bins.
                     -
                   */
                  for(Int_t iRapidityBin = 0; iRapidityBin < 10; iRapidityBin++){
                      if( (possibleJPsiMC.Rapidity() + 4.) < 1.5*((Double_t)iRapidityBin + 1.)/10. ){
                        fMCCosThetaHelicityFrameJPsiTenRapidityBinsH[iRapidityBin]->Fill( CosThetaHelicityFrame( muonsMCcopy[1],
                                                                                                                 muonsMCcopy[0],
                                                                                                                 possibleJPsiMC
                                                                                                                 )
                                                                                                            );
                        fMCCosThetaCollinsSoperFrameJPsiTenRapidityBinsH[iRapidityBin]->Fill( CosThetaCollinsSoper( muonsMCcopy[1],
                                                                                                                    muonsMCcopy[0],
                                                                                                                    possibleJPsiMC
                                                                                                                    )
                                                                                                               );
                        fMCPhiHelicityFrameJPsiTenRapidityBinsH[iRapidityBin]->Fill( CosPhiHelicityFrame( muonsMCcopy[1],
                                                                                                          muonsMCcopy[0],
                                                                                                          possibleJPsiMC
                                                                                                          )
                                                                                                     );
                        fMCPhiCollinsSoperFrameJPsiTenRapidityBinsH[iRapidityBin]->Fill( CosPhiCollinsSoper( muonsMCcopy[1],
                                                                                                             muonsMCcopy[0],
                                                                                                             possibleJPsiMC
                                                                                                             )
                                                                                                        );
                        break;
                      }
                  }
                  fMCCosThetaAndPhiHelicityFrameMyBinningH->Fill( CosThetaHelicityFrame( muonsMCcopy[1],
                                                                                         muonsMCcopy[0],
                                                                                         possibleJPsiMC
                                                                                         ),
                                                                  CosPhiHelicityFrame( muonsMCcopy[1],
                                                                                       muonsMCcopy[0],
                                                                                       possibleJPsiMC
                                                                                       )
                                                                  );
                  fMCCosThetaAndPhiCsFrameMyBinningH->Fill( CosThetaCollinsSoper( muonsMCcopy[1],
                                                                                  muonsMCcopy[0],
                                                                                  possibleJPsiMC
                                                                                  ),
                                                            CosPhiCollinsSoper(   muonsMCcopy[1],
                                                                                  muonsMCcopy[0],
                                                                                  possibleJPsiMC
                                                                                  )
                                                            );

                  /* - What we do here is very similar.
                     - This time we divide firstly in bins of CosTheta.
                     - As many as needed.
                     - And then we divide again in terms of Phi.
                     - Then we fill.
                     - This way we should be able to obtain some kind of map...
                     -
                   */
                  fMCInvariantMassDistributionBinsOfCosThetaAndPhiHelicityFrameH->Fill( CosThetaHelicityFrame( muonsMCcopy[1],
                                                                                                               muonsMCcopy[0],
                                                                                                               possibleJPsiMC
                                                                                                               ),
                                                                                        CosPhiHelicityFrame( muonsMCcopy[1],
                                                                                                             muonsMCcopy[0],
                                                                                                             possibleJPsiMC
                                                                                                             )
                                                                                        );
                  fMCCosThetaAndPhiHelicityFrameInclusivePeopleBinningH->Fill( CosThetaHelicityFrame( muonsMCcopy[1],
                                                                                                      muonsMCcopy[0],
                                                                                                      possibleJPsiMC
                                                                                                      ),
                                                                               CosPhiHelicityFrame( muonsMCcopy[1],
                                                                                                    muonsMCcopy[0],
                                                                                                    possibleJPsiMC
                                                                                                    )
                                                                               );
          }
      }
    }
  }

}
//_____________________________________________________________________________
/* - The following are code snippets adapted from the AliAODDimuon class.
   - The problem is that that class was adapted specifically for the
   - inclusive people's analysis, hence it is not fit for the UPC...
   -
 */
Double_t AliAnalysisTaskUPCforwardMC::CosThetaCollinsSoper( TLorentzVector muonPositive,
                                                            TLorentzVector muonNegative,
                                                            TLorentzVector possibleJPsi )
{
  /* - This function computes the Collins-Soper cos(theta) for the
     - helicity of the J/Psi.
     - The idea should be to get back to a reference frame where it
     - is easier to compute and to define the proper z-axis.
     -
   */

  /* - Half of the energy per pair of the colliding nucleons.
     -
   */
  Double_t HalfSqrtSnn   = 2510.;
  Double_t MassOfLead208 = 193.6823;
  Double_t MomentumBeam  = TMath::Sqrt( HalfSqrtSnn*HalfSqrtSnn*208*208 - MassOfLead208*MassOfLead208 );
  /* - Fill the Lorentz vector for projectile and target.
     - For the moment we do not consider the crossing angle.
     - Projectile runs towards the MUON arm.
     -
   */
  TLorentzVector pProjCM(0.,0., -MomentumBeam, HalfSqrtSnn*208); // projectile
  TLorentzVector pTargCM(0.,0.,  MomentumBeam, HalfSqrtSnn*208); // target
  /* - Translate the dimuon parameters in the dimuon rest frame
     -
   */
  TVector3       beta      = ( -1./possibleJPsi.E() ) * possibleJPsi.Vect();
  TLorentzVector pMu1Dimu  = muonPositive;
  TLorentzVector pMu2Dimu  = muonNegative;
  TLorentzVector pProjDimu = pProjCM;
  TLorentzVector pTargDimu = pTargCM;
  pMu1Dimu.Boost(beta);
  pMu2Dimu.Boost(beta);
  pProjDimu.Boost(beta);
  pTargDimu.Boost(beta);
  /* - Determine the z axis for the CS angle.
     -
   */
  TVector3 zaxisCS=(((pProjDimu.Vect()).Unit())-((pTargDimu.Vect()).Unit())).Unit();
  /* - Determine the CS angle (angle between mu+ and the z axis defined above)
     -
   */
  Double_t CosThetaCS = zaxisCS.Dot((pMu1Dimu.Vect()).Unit());
  return   CosThetaCS;
}
//_____________________________________________________________________________
Double_t AliAnalysisTaskUPCforwardMC::CosThetaHelicityFrame( TLorentzVector muonPositive,
                                                             TLorentzVector muonNegative,
                                                             TLorentzVector possibleJPsi )
{
  /* - This function computes the Helicity cos(theta) for the
     - helicity of the J/Psi.
     - The idea should be to get back to a reference frame where it
     - is easier to compute and to define the proper z-axis.
     -
   */

  /* - Half of the energy per pair of the colliding nucleons.
     -
   */
  Double_t HalfSqrtSnn   = 2510.;
  Double_t MassOfLead208 = 193.6823;
  Double_t MomentumBeam  = TMath::Sqrt( HalfSqrtSnn*HalfSqrtSnn*208*208 - MassOfLead208*MassOfLead208 );
  /* - Fill the Lorentz vector for projectile and target.
     - For the moment we do not consider the crossing angle.
     - Projectile runs towards the MUON arm.
     -
   */
  TLorentzVector pProjCM(0.,0., -MomentumBeam, HalfSqrtSnn*208); // projectile
  TLorentzVector pTargCM(0.,0.,  MomentumBeam, HalfSqrtSnn*208); // target
  /* - Translate the dimuon parameters in the dimuon rest frame
     -
   */
  TVector3       beta      = ( -1./possibleJPsi.E() ) * possibleJPsi.Vect();
  TLorentzVector pMu1Dimu  = muonPositive;
  TLorentzVector pMu2Dimu  = muonNegative;
  TLorentzVector pProjDimu = pProjCM;
  TLorentzVector pTargDimu = pTargCM;
  pMu1Dimu.Boost(beta);
  pMu2Dimu.Boost(beta);
  pProjDimu.Boost(beta);
  pTargDimu.Boost(beta);
  //
  // --- Determine the z axis for the calculation of the polarization angle
  // (i.e. the direction of the dimuon in the CM system)
  //
  TVector3 zaxis = (possibleJPsi.Vect()).Unit();
  /* - Determine the He angle (angle between mu+ and the z axis defined above)
     -
   */
  Double_t CosThetaHE = zaxis.Dot((pMu1Dimu.Vect()).Unit());
  return   CosThetaHE;

}
//_____________________________________________________________________________
Double_t AliAnalysisTaskUPCforwardMC::CosPhiCollinsSoper( TLorentzVector muonPositive,
                                                          TLorentzVector muonNegative,
                                                          TLorentzVector possibleJPsi )
{
  /* - This function computes the Collins-Soper PHI for the
     - helicity of the J/Psi.
     - The idea should be to get back to a reference frame where it
     - is easier to compute and to define the proper z-axis.
     -
   */

  /* - Half of the energy per pair of the colliding nucleons.
     -
   */
  Double_t HalfSqrtSnn   = 2510.;
  Double_t MassOfLead208 = 193.6823;
  Double_t MomentumBeam  = TMath::Sqrt( HalfSqrtSnn*HalfSqrtSnn*208*208 - MassOfLead208*MassOfLead208 );
  /* - Fill the Lorentz vector for projectile and target.
     - For the moment we do not consider the crossing angle.
     - Projectile runs towards the MUON arm.
     -
   */
  TLorentzVector pProjCM(0.,0., -MomentumBeam, HalfSqrtSnn*208); // projectile
  TLorentzVector pTargCM(0.,0.,  MomentumBeam, HalfSqrtSnn*208); // target
  /* - Translate the dimuon parameters in the dimuon rest frame
     -
   */
  TVector3       beta      = ( -1./possibleJPsi.E() ) * possibleJPsi.Vect();
  TLorentzVector pMu1Dimu  = muonPositive;
  TLorentzVector pMu2Dimu  = muonNegative;
  TLorentzVector pProjDimu = pProjCM;
  TLorentzVector pTargDimu = pTargCM;
  pMu1Dimu.Boost(beta);
  pMu2Dimu.Boost(beta);
  pProjDimu.Boost(beta);
  pTargDimu.Boost(beta);
  /* - Determine the z axis for the CS angle.
     -
   */
  TVector3 zaxisCS=(((pProjDimu.Vect()).Unit())-((pTargDimu.Vect()).Unit())).Unit();
  //
  // --- Determine the CS angle (angle between mu+ and the z axis defined above)
  //
  TVector3 yaxisCS=(((pProjDimu.Vect()).Unit()).Cross((pTargDimu.Vect()).Unit())).Unit();
  TVector3 xaxisCS=(yaxisCS.Cross(zaxisCS)).Unit();

  Double_t phi = TMath::ATan2((pMu1Dimu.Vect()).Dot(yaxisCS),((pMu1Dimu.Vect()).Dot(xaxisCS)));
  return   phi;
}
//_____________________________________________________________________________
Double_t AliAnalysisTaskUPCforwardMC::CosPhiHelicityFrame(  TLorentzVector muonPositive,
                                                            TLorentzVector muonNegative,
                                                            TLorentzVector possibleJPsi )
{
  /* - This function computes the helicity phi for the
     - helicity of the J/Psi.
     - The idea should be to get back to a reference frame where it
     - is easier to compute and to define the proper z-axis.
     -
   */

  /* - Half of the energy per pair of the colliding nucleons.
     -
  */
  Double_t HalfSqrtSnn   = 2510.;
  Double_t MassOfLead208 = 193.6823;
  Double_t MomentumBeam  = TMath::Sqrt( HalfSqrtSnn*HalfSqrtSnn*208*208 - MassOfLead208*MassOfLead208 );
  /* - Fill the Lorentz vector for projectile and target.
     - For the moment we do not consider the crossing angle.
     - Projectile runs towards the MUON arm.
     -
   */
  TLorentzVector pProjCM(0.,0., -MomentumBeam, HalfSqrtSnn*208); // projectile
  TLorentzVector pTargCM(0.,0.,  MomentumBeam, HalfSqrtSnn*208); // target
  /* - Translate the dimuon parameters in the dimuon rest frame
     -
   */
  TVector3       beta      = ( -1./possibleJPsi.E() ) * possibleJPsi.Vect();
  TLorentzVector pMu1Dimu  = muonPositive;
  TLorentzVector pMu2Dimu  = muonNegative;
  TLorentzVector pProjDimu = pProjCM;
  TLorentzVector pTargDimu = pTargCM;
  pMu1Dimu.Boost(beta);
  pMu2Dimu.Boost(beta);
  pProjDimu.Boost(beta);
  pTargDimu.Boost(beta);
  //
  // --- Determine the z axis for the calculation of the polarization angle
  // (i.e. the direction of the dimuon in the CM system)
  //
  TVector3 zaxis = (possibleJPsi.Vect()).Unit();
  TVector3 yaxis = ((pProjDimu.Vect()).Cross(pTargDimu.Vect())).Unit();
  TVector3 xaxis = (yaxis.Cross(zaxis)).Unit();
  //
  // --- Calculation of the azimuthal angle (Helicity)
  //
  Double_t phi = TMath::ATan2((pMu1Dimu.Vect()).Dot(yaxis),(pMu1Dimu.Vect()).Dot(xaxis));
  return   phi;
}
//_____________________________________________________________________________
void AliAnalysisTaskUPCforwardMC::Terminate(Option_t *)
{
    cout << endl;
    // terminate
    // called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________
/* - There were plenty of ways to do this...
 * - However, recently the STL libraries were
 * - creating confusion on the LEGO framework
 * - (they didn't fire at all).
 * - This problem was not found on local, where
 * - things were working properly...
 * - So I am using the most barbaric C-style
 * - arrays/for...
 */
void AliAnalysisTaskUPCforwardMC::SetLuminosityCap()
{
  fLumiPerRun = 0;
  /* - Here I am rounding up the number for 10k,
   * - so that I don't have to do tedious conversions...
   * - I am adding 1 entry to the number obtained by 40k,
   * - so that I am not missing anything...
   * -
   */
  if      ( fRunNum == 244980 ) { fLumiPerRun = 0.0505; }
  else if ( fRunNum == 244982 ) { fLumiPerRun = 0.0761; }
  else if ( fRunNum == 244983 ) { fLumiPerRun = 0.0291; }
  else if ( fRunNum == 245064 ) { fLumiPerRun = 0.1643; }
  else if ( fRunNum == 245066 ) { fLumiPerRun = 0.0236; }
  else if ( fRunNum == 245068 ) { fLumiPerRun = 0.0202; }
  else if ( fRunNum == 245145 ) { fLumiPerRun = 1.2115; }
  else if ( fRunNum == 245146 ) { fLumiPerRun = 1.3773; }
  else if ( fRunNum == 245151 ) { fLumiPerRun = 0.1469; }
  else if ( fRunNum == 245152 ) { fLumiPerRun = 0.1655; }
  else if ( fRunNum == 245231 ) { fLumiPerRun = 0.3084; }
  else if ( fRunNum == 245232 ) { fLumiPerRun = 1.0146; }
  else if ( fRunNum == 245233 ) { fLumiPerRun = 0.2373; }
  else if ( fRunNum == 245253 ) { fLumiPerRun = 0.3067; }
  else if ( fRunNum == 245259 ) { fLumiPerRun = 0.4893; }
  else if ( fRunNum == 245343 ) { fLumiPerRun = 0.7006; }
  else if ( fRunNum == 245345 ) { fLumiPerRun = 2.2153; }
  else if ( fRunNum == 245346 ) { fLumiPerRun = 0.2785; }
  else if ( fRunNum == 245347 ) { fLumiPerRun = 1.1752; }
  else if ( fRunNum == 245353 ) { fLumiPerRun = 1.6505; }
  else if ( fRunNum == 245401 ) { fLumiPerRun = 0.7485; }
  else if ( fRunNum == 245407 ) { fLumiPerRun = 2.0625; }
  else if ( fRunNum == 245409 ) { fLumiPerRun = 0.8705; }
  else if ( fRunNum == 245410 ) { fLumiPerRun = 0.1819; }
  else if ( fRunNum == 245446 ) { fLumiPerRun = 0.1261; }
  else if ( fRunNum == 245450 ) { fLumiPerRun = 0.2621; }
  else if ( fRunNum == 245496 ) { fLumiPerRun = 1.06;   }
  else if ( fRunNum == 245501 ) { fLumiPerRun = 1.334;  }
  else if ( fRunNum == 245504 ) { fLumiPerRun = 0.6492; }
  else if ( fRunNum == 245505 ) { fLumiPerRun = 0.3624; }
  else if ( fRunNum == 245507 ) { fLumiPerRun = 1.6192; }
  else if ( fRunNum == 245535 ) { fLumiPerRun = 1.3612; }
  else if ( fRunNum == 245540 ) { fLumiPerRun = 0.7121; }
  else if ( fRunNum == 245542 ) { fLumiPerRun = 1.1181; }
  else if ( fRunNum == 245543 ) { fLumiPerRun = 2.0169; }
  else if ( fRunNum == 245554 ) { fLumiPerRun = 1.7248; }
  else if ( fRunNum == 245683 ) { fLumiPerRun = 4.0406; }
  else if ( fRunNum == 245692 ) { fLumiPerRun = 1.9090; }
  else if ( fRunNum == 245700 ) { fLumiPerRun = 1.1167; }
  else if ( fRunNum == 245705 ) { fLumiPerRun = 0.3239; }
  else if ( fRunNum == 245729 ) { fLumiPerRun = 1.1548; }
  else if ( fRunNum == 245731 ) { fLumiPerRun = 3.3932; }
  else if ( fRunNum == 245738 ) { fLumiPerRun = 1.9485; }
  else if ( fRunNum == 245752 ) { fLumiPerRun = 1.2497; }
  else if ( fRunNum == 245759 ) { fLumiPerRun = 1.3785; }
  else if ( fRunNum == 245766 ) { fLumiPerRun = 1.1429; }
  else if ( fRunNum == 245775 ) { fLumiPerRun = 1.7326; }
  else if ( fRunNum == 245785 ) { fLumiPerRun = 0.5102; }
  else if ( fRunNum == 245793 ) { fLumiPerRun = 0.7093; }
  else if ( fRunNum == 245829 ) { fLumiPerRun = 1.958;  }
  else if ( fRunNum == 245831 ) { fLumiPerRun = 1.9939; }
  else if ( fRunNum == 245833 ) { fLumiPerRun = 0.3559; }
  else if ( fRunNum == 245949 ) { fLumiPerRun = 0.5652; }
  else if ( fRunNum == 245952 ) { fLumiPerRun = 3.0759; }
  else if ( fRunNum == 245954 ) { fLumiPerRun = 1.9965; }
  else if ( fRunNum == 245963 ) { fLumiPerRun = 2.2815; }
  else if ( fRunNum == 245996 ) { fLumiPerRun = 0.4644; }
  else if ( fRunNum == 246001 ) { fLumiPerRun = 3.5684; }
  else if ( fRunNum == 246003 ) { fLumiPerRun = 0.5803; }
  else if ( fRunNum == 246012 ) { fLumiPerRun = 0.7302; }
  else if ( fRunNum == 246036 ) { fLumiPerRun = 0.2143; }
  else if ( fRunNum == 246037 ) { fLumiPerRun = 1.7465; }
  else if ( fRunNum == 246042 ) { fLumiPerRun = 4.8713; }
  else if ( fRunNum == 246048 ) { fLumiPerRun = 0.3835; }
  else if ( fRunNum == 246049 ) { fLumiPerRun = 3.2666; }
  else if ( fRunNum == 246053 ) { fLumiPerRun = 1.7691; }
  else if ( fRunNum == 246087 ) { fLumiPerRun = 14.184; }
  else if ( fRunNum == 246089 ) { fLumiPerRun = 0.3296; }
  else if ( fRunNum == 246113 ) { fLumiPerRun = 1.4761; }
  else if ( fRunNum == 246115 ) { fLumiPerRun = 0.4514; }
  else if ( fRunNum == 246148 ) { fLumiPerRun = 5.3175; }
  else if ( fRunNum == 246151 ) { fLumiPerRun = 3.0605; }
  else if ( fRunNum == 246152 ) { fLumiPerRun = 0.4734; }
  else if ( fRunNum == 246153 ) { fLumiPerRun = 4.6676; }
  else if ( fRunNum == 246178 ) { fLumiPerRun = 0.8156; }
  else if ( fRunNum == 246181 ) { fLumiPerRun = 2.7526; }
  else if ( fRunNum == 246182 ) { fLumiPerRun = 2.2047; }
  else if ( fRunNum == 246217 ) { fLumiPerRun = 3.4663; }
  else if ( fRunNum == 246220 ) { fLumiPerRun = 0.682;  }
  else if ( fRunNum == 246222 ) { fLumiPerRun = 3.6826; }
  else if ( fRunNum == 246225 ) { fLumiPerRun = 1.2534; }
  else if ( fRunNum == 246272 ) { fLumiPerRun = 5.5294; }
  else if ( fRunNum == 246275 ) { fLumiPerRun = 1.242;  }
  else if ( fRunNum == 246276 ) { fLumiPerRun = 0.5871; }
  else if ( fRunNum == 246390 ) { fLumiPerRun = 0.0448; }
  else if ( fRunNum == 246391 ) { fLumiPerRun = 0.1446; }
  else if ( fRunNum == 246392 ) { fLumiPerRun = 0.1765; }
  else if ( fRunNum == 246424 ) { fLumiPerRun = 2.866;  }
  else if ( fRunNum == 246428 ) { fLumiPerRun = 0.4417; }
  else if ( fRunNum == 246431 ) { fLumiPerRun = 1.7836; }
  else if ( fRunNum == 246433 ) { fLumiPerRun = 0.4164; }
  else if ( fRunNum == 246434 ) { fLumiPerRun = 4.103;  }
  else if ( fRunNum == 246487 ) { fLumiPerRun = 0.7286; }
  else if ( fRunNum == 246488 ) { fLumiPerRun = 7.5895; }
  else if ( fRunNum == 246493 ) { fLumiPerRun = 1.3534; }
  else if ( fRunNum == 246495 ) { fLumiPerRun = 0.4100; }
  else if ( fRunNum == 246675 ) { fLumiPerRun = 2.3469; }
  else if ( fRunNum == 246676 ) { fLumiPerRun = 0.4794; }
  else if ( fRunNum == 246750 ) { fLumiPerRun = 2.0756; }
  else if ( fRunNum == 246751 ) { fLumiPerRun = 2.0419; }
  else if ( fRunNum == 246755 ) { fLumiPerRun = 1.4197; }
  else if ( fRunNum == 246757 ) { fLumiPerRun = 0.59;   }
  else if ( fRunNum == 246758 ) { fLumiPerRun = 1.626;  }
  else if ( fRunNum == 246759 ) { fLumiPerRun = 0.3335; }
  else if ( fRunNum == 246760 ) { fLumiPerRun = 1.1753; }
  else if ( fRunNum == 246763 ) { fLumiPerRun = 0.549;  }
  else if ( fRunNum == 246765 ) { fLumiPerRun = 0.3274; }
  else if ( fRunNum == 246804 ) { fLumiPerRun = 1.0208; }
  else if ( fRunNum == 246805 ) { fLumiPerRun = 3.1925; }
  else if ( fRunNum == 246806 ) { fLumiPerRun = 2.5555; }
  else if ( fRunNum == 246807 ) { fLumiPerRun = 2.5962; }
  else if ( fRunNum == 246808 ) { fLumiPerRun = 0.3101; }
  else if ( fRunNum == 246809 ) { fLumiPerRun = 2.4707; }
  else if ( fRunNum == 246844 ) { fLumiPerRun = 0.7657; }
  else if ( fRunNum == 246845 ) { fLumiPerRun = 1.4355; }
  else if ( fRunNum == 246846 ) { fLumiPerRun = 0.8986; }
  else if ( fRunNum == 246847 ) { fLumiPerRun = 1.7064; }
  else if ( fRunNum == 246851 ) { fLumiPerRun = 1.2170; }
  else if ( fRunNum == 246855 ) { fLumiPerRun = 1.3014; }
  else if ( fRunNum == 246859 ) { fLumiPerRun = 1.2397; }
  else if ( fRunNum == 246864 ) { fLumiPerRun = 2.4832; }
  else if ( fRunNum == 246865 ) { fLumiPerRun = 0.8111; }
  else if ( fRunNum == 246867 ) { fLumiPerRun = 1.5019; }
  else if ( fRunNum == 246871 ) { fLumiPerRun = 0.8713; }
  else if ( fRunNum == 246930 ) { fLumiPerRun = 0.5641; }
  else if ( fRunNum == 246937 ) { fLumiPerRun = 0.699;  }
  else if ( fRunNum == 246942 ) { fLumiPerRun = 1.0555; }
  else if ( fRunNum == 246945 ) { fLumiPerRun = 2.1676; }
  else if ( fRunNum == 246948 ) { fLumiPerRun = 0.8855; }
  else if ( fRunNum == 246949 ) { fLumiPerRun = 2.8978; }
  else if ( fRunNum == 246980 ) { fLumiPerRun = 7.1999; }
  else if ( fRunNum == 246982 ) { fLumiPerRun = 0.5146; }
  else if ( fRunNum == 246984 ) { fLumiPerRun = 4.143;  }
  else if ( fRunNum == 246989 ) { fLumiPerRun = 3.8342; }
  else if ( fRunNum == 246991 ) { fLumiPerRun = 0.4368; }
  else if ( fRunNum == 246994 ) { fLumiPerRun = 1.2329; }
  //_______________________________
  /* -
   * -
   * - POLARISATION ANALYSIS
   * -
   */
  else if ( fRunNum == 295585 ) { fLumiPerRun = 0.0793; }
  else if ( fRunNum == 295586 ) { fLumiPerRun = 0.2386; }
  else if ( fRunNum == 295587 ) { fLumiPerRun = 0.1095; }
  else if ( fRunNum == 295588 ) { fLumiPerRun = 0.1358; }
  else if ( fRunNum == 295589 ) { fLumiPerRun = 0.2819; }
  else if ( fRunNum == 295612 ) { fLumiPerRun = 0.449;  }
  else if ( fRunNum == 295615 ) { fLumiPerRun = 0.0566; }
  else if ( fRunNum == 295665 ) { fLumiPerRun = 0.3349; }
  else if ( fRunNum == 295666 ) { fLumiPerRun = 0.3239; }
  else if ( fRunNum == 295667 ) { fLumiPerRun = 0.0970; }
  else if ( fRunNum == 295668 ) { fLumiPerRun = 0.1303; }
  else if ( fRunNum == 295671 ) { fLumiPerRun = 0.3259; }
  else if ( fRunNum == 295673 ) { fLumiPerRun = 0.3128; }
  else if ( fRunNum == 295675 ) { fLumiPerRun = 0.132;  }
  else if ( fRunNum == 295676 ) { fLumiPerRun = 0.3213; }
  else if ( fRunNum == 295677 ) { fLumiPerRun = 0.2652; }
  else if ( fRunNum == 295714 ) { fLumiPerRun = 0.3456; }
  else if ( fRunNum == 295716 ) { fLumiPerRun = 0.3389; }
  else if ( fRunNum == 295717 ) { fLumiPerRun = 0.2880; }
  else if ( fRunNum == 295718 ) { fLumiPerRun = 0.2567; }
  else if ( fRunNum == 295719 ) { fLumiPerRun = 0.2947; }
  else if ( fRunNum == 295723 ) { fLumiPerRun = 0.5064; }
  else if ( fRunNum == 295725 ) { fLumiPerRun = 0.8890; }
  else if ( fRunNum == 295753 ) { fLumiPerRun = 0.3846; }
  else if ( fRunNum == 295754 ) { fLumiPerRun = 0.7055; }
  else if ( fRunNum == 295755 ) { fLumiPerRun = 0.7585; }
  else if ( fRunNum == 295758 ) { fLumiPerRun = 1.8934; }
  else if ( fRunNum == 295759 ) { fLumiPerRun = 0.5331; }
  else if ( fRunNum == 295762 ) { fLumiPerRun = 0.2749; }
  else if ( fRunNum == 295763 ) { fLumiPerRun = 1.0282; }
  else if ( fRunNum == 295786 ) { fLumiPerRun = 0.7490; }
  else if ( fRunNum == 295788 ) { fLumiPerRun = 3.0237; }
  else if ( fRunNum == 295791 ) { fLumiPerRun = 0.8580; }
  else if ( fRunNum == 295816 ) { fLumiPerRun = 1.2056; }
  else if ( fRunNum == 295818 ) { fLumiPerRun = 0.1453; }
  else if ( fRunNum == 295819 ) { fLumiPerRun = 2.7474; }
  else if ( fRunNum == 295822 ) { fLumiPerRun = 2.2529; }
  else if ( fRunNum == 295825 ) { fLumiPerRun = 0.2558; }
  else if ( fRunNum == 295826 ) { fLumiPerRun = 1.5814; }
  else if ( fRunNum == 295829 ) { fLumiPerRun = 0.9351; }
  else if ( fRunNum == 295831 ) { fLumiPerRun = 0.7762; }
  else if ( fRunNum == 295854 ) { fLumiPerRun = 1.3119; }
  else if ( fRunNum == 295855 ) { fLumiPerRun = 1.7466; }
  else if ( fRunNum == 295856 ) { fLumiPerRun = 1.4700; }
  else if ( fRunNum == 295859 ) { fLumiPerRun = 1.0510; }
  else if ( fRunNum == 295860 ) { fLumiPerRun = 0.8341; }
  else if ( fRunNum == 295861 ) { fLumiPerRun = 1.0670; }
  else if ( fRunNum == 295863 ) { fLumiPerRun = 0.7279; }
  else if ( fRunNum == 295881 ) { fLumiPerRun = 0.7115; }
  else if ( fRunNum == 295908 ) { fLumiPerRun = 2.9261; }
  else if ( fRunNum == 295909 ) { fLumiPerRun = 0.7875; }
  else if ( fRunNum == 295910 ) { fLumiPerRun = 3.1843; }
  else if ( fRunNum == 295913 ) { fLumiPerRun = 3.1294; }
  else if ( fRunNum == 295936 ) { fLumiPerRun = 1.4736; }
  else if ( fRunNum == 295937 ) { fLumiPerRun = 0.4057; }
  else if ( fRunNum == 295941 ) { fLumiPerRun = 1.6767; }
  else if ( fRunNum == 295942 ) { fLumiPerRun = 1.9237; }
  else if ( fRunNum == 295943 ) { fLumiPerRun = 1.6747; }
  else if ( fRunNum == 295945 ) { fLumiPerRun = 2.0370; }
  else if ( fRunNum == 295947 ) { fLumiPerRun = 2.6337; }
  else if ( fRunNum == 296061 ) { fLumiPerRun = 1.2968; }
  else if ( fRunNum == 296062 ) { fLumiPerRun = 1.8083; }
  else if ( fRunNum == 296063 ) { fLumiPerRun = 2.6876; }
  else if ( fRunNum == 296065 ) { fLumiPerRun = 2.4473; }
  else if ( fRunNum == 296066 ) { fLumiPerRun = 0.7337; }
  else if ( fRunNum == 296068 ) { fLumiPerRun = 1.9812; }
  else if ( fRunNum == 296123 ) { fLumiPerRun = 0.5065; }
  else if ( fRunNum == 296128 ) { fLumiPerRun = 0.4455; }
  else if ( fRunNum == 296132 ) { fLumiPerRun = 1.312;  }
  else if ( fRunNum == 296133 ) { fLumiPerRun = 1.7321; }
  else if ( fRunNum == 296134 ) { fLumiPerRun = 3.9104; }
  else if ( fRunNum == 296135 ) { fLumiPerRun = 2.3412; }
  else if ( fRunNum == 296142 ) { fLumiPerRun = 1.7893; }
  else if ( fRunNum == 296143 ) { fLumiPerRun = 0.5340; }
  else if ( fRunNum == 296191 ) { fLumiPerRun = 5.0507; }
  else if ( fRunNum == 296192 ) { fLumiPerRun = 0.4974; }
  else if ( fRunNum == 296194 ) { fLumiPerRun = 2.8725; }
  else if ( fRunNum == 296195 ) { fLumiPerRun = 0.7376; }
  else if ( fRunNum == 296196 ) { fLumiPerRun = 2.352;  }
  else if ( fRunNum == 296197 ) { fLumiPerRun = 2.0691; }
  else if ( fRunNum == 296198 ) { fLumiPerRun = 0.8140; }
  else if ( fRunNum == 296241 ) { fLumiPerRun = 0.8459; }
  else if ( fRunNum == 296242 ) { fLumiPerRun = 0.9517; }
  else if ( fRunNum == 296243 ) { fLumiPerRun = 1.5674; }
  else if ( fRunNum == 296244 ) { fLumiPerRun = 8.3722; }
  else if ( fRunNum == 296246 ) { fLumiPerRun = 1.8351; }
  else if ( fRunNum == 296247 ) { fLumiPerRun = 1.1765; }
  else if ( fRunNum == 296269 ) { fLumiPerRun = 3.8392; }
  else if ( fRunNum == 296270 ) { fLumiPerRun = 1.5116; }
  else if ( fRunNum == 296273 ) { fLumiPerRun = 7.2237; }
  else if ( fRunNum == 296279 ) { fLumiPerRun = 0.4057; }
  else if ( fRunNum == 296280 ) { fLumiPerRun = 1.5066; }
  else if ( fRunNum == 296303 ) { fLumiPerRun = 2.006;  }
  else if ( fRunNum == 296304 ) { fLumiPerRun = 6.0965; }
  else if ( fRunNum == 296307 ) { fLumiPerRun = 2.9023; }
  else if ( fRunNum == 296309 ) { fLumiPerRun = 2.1026; }
  else if ( fRunNum == 296312 ) { fLumiPerRun = 2.1228; }
  else if ( fRunNum == 296377 ) { fLumiPerRun = 6.0666; }
  else if ( fRunNum == 296378 ) { fLumiPerRun = 5.3897; }
  else if ( fRunNum == 296379 ) { fLumiPerRun = 2.0969; }
  else if ( fRunNum == 296380 ) { fLumiPerRun = 2.8820; }
  else if ( fRunNum == 296381 ) { fLumiPerRun = 1.4418; }
  else if ( fRunNum == 296383 ) { fLumiPerRun = 1.5136; }
  else if ( fRunNum == 296414 ) { fLumiPerRun = 4.8766; }
  else if ( fRunNum == 296419 ) { fLumiPerRun = 2.7523; }
  else if ( fRunNum == 296420 ) { fLumiPerRun = 1.4132; }
  else if ( fRunNum == 296423 ) { fLumiPerRun = 1.5981; }
  else if ( fRunNum == 296424 ) { fLumiPerRun = 0.3864; }
  else if ( fRunNum == 296433 ) { fLumiPerRun = 4.0456; }
  else if ( fRunNum == 296472 ) { fLumiPerRun = 0.8632; }
  else if ( fRunNum == 296509 ) { fLumiPerRun = 2.9592; }
  else if ( fRunNum == 296510 ) { fLumiPerRun = 9.0673; }
  else if ( fRunNum == 296511 ) { fLumiPerRun = 2.5666; }
  else if ( fRunNum == 296514 ) { fLumiPerRun = 0.4898; }
  else if ( fRunNum == 296516 ) { fLumiPerRun = 0.6134; }
  else if ( fRunNum == 296547 ) { fLumiPerRun = 1.0834; }
  else if ( fRunNum == 296548 ) { fLumiPerRun = 1.3771; }
  else if ( fRunNum == 296549 ) { fLumiPerRun = 4.8645; }
  else if ( fRunNum == 296550 ) { fLumiPerRun = 3.9901; }
  else if ( fRunNum == 296551 ) { fLumiPerRun = 2.0214; }
  else if ( fRunNum == 296552 ) { fLumiPerRun = 0.4842; }
  else if ( fRunNum == 296553 ) { fLumiPerRun = 0.7091; }
  else if ( fRunNum == 296615 ) { fLumiPerRun = 1.5676; }
  else if ( fRunNum == 296616 ) { fLumiPerRun = 0.5399; }
  else if ( fRunNum == 296618 ) { fLumiPerRun = 1.7014; }
  else if ( fRunNum == 296619 ) { fLumiPerRun = 1.5613; }
  else if ( fRunNum == 296622 ) { fLumiPerRun = 0.7064; }
  else if ( fRunNum == 296623 ) { fLumiPerRun = 2.1442; }
  else if ( fRunNum == 296690 ) { fLumiPerRun = 6.8615; }
  else if ( fRunNum == 296691 ) { fLumiPerRun = 0.6511; }
  else if ( fRunNum == 296694 ) { fLumiPerRun = 5.1826; }
  else if ( fRunNum == 296749 ) { fLumiPerRun = 9.2413; }
  else if ( fRunNum == 296750 ) { fLumiPerRun = 8.2161; }
  else if ( fRunNum == 296781 ) { fLumiPerRun = 0.8179; }
  else if ( fRunNum == 296784 ) { fLumiPerRun = 2.98;   }
  else if ( fRunNum == 296785 ) { fLumiPerRun = 1.9085; }
  else if ( fRunNum == 296786 ) { fLumiPerRun = 0.7537; }
  else if ( fRunNum == 296787 ) { fLumiPerRun = 3.2190; }
  else if ( fRunNum == 296791 ) { fLumiPerRun = 0.7573; }
  else if ( fRunNum == 296793 ) { fLumiPerRun = 1.3317; }
  else if ( fRunNum == 296794 ) { fLumiPerRun = 3.1335; }
  else if ( fRunNum == 296799 ) { fLumiPerRun = 2.7149; }
  else if ( fRunNum == 296836 ) { fLumiPerRun = 1.5116; }
  else if ( fRunNum == 296838 ) { fLumiPerRun = 0.5432; }
  else if ( fRunNum == 296839 ) { fLumiPerRun = 2.9424; }
  else if ( fRunNum == 296848 ) { fLumiPerRun = 2.1628; }
  else if ( fRunNum == 296849 ) { fLumiPerRun = 11.469; }
  else if ( fRunNum == 296850 ) { fLumiPerRun = 2.7979; }
  else if ( fRunNum == 296851 ) { fLumiPerRun = 0.1392; }
  else if ( fRunNum == 296852 ) { fLumiPerRun = 0.9565; }
  else if ( fRunNum == 296890 ) { fLumiPerRun = 8.0545; }
  else if ( fRunNum == 296894 ) { fLumiPerRun = 4.6472; }
  else if ( fRunNum == 296899 ) { fLumiPerRun = 2.1355; }
  else if ( fRunNum == 296900 ) { fLumiPerRun = 2.7833; }
  else if ( fRunNum == 296903 ) { fLumiPerRun = 1.0391; }
  else if ( fRunNum == 296930 ) { fLumiPerRun = 1.4575; }
  else if ( fRunNum == 296931 ) { fLumiPerRun = 0.5292; }
  else if ( fRunNum == 296932 ) { fLumiPerRun = 1.1863; }
  else if ( fRunNum == 296934 ) { fLumiPerRun = 2.5917; }
  else if ( fRunNum == 296935 ) { fLumiPerRun = 4.4039; }
  else if ( fRunNum == 296938 ) { fLumiPerRun = 1.6678; }
  else if ( fRunNum == 296941 ) { fLumiPerRun = 2.9181; }
  else if ( fRunNum == 296966 ) { fLumiPerRun = 3.3611; }
  else if ( fRunNum == 296967 ) { fLumiPerRun = 0.8051; }
  else if ( fRunNum == 296968 ) { fLumiPerRun = 3.1905; }
  else if ( fRunNum == 296969 ) { fLumiPerRun = 1.8878; }
  else if ( fRunNum == 296971 ) { fLumiPerRun = 0.6907; }
  else if ( fRunNum == 296975 ) { fLumiPerRun = 7.3683; }
  else if ( fRunNum == 296976 ) { fLumiPerRun = 1.1175; }
  else if ( fRunNum == 296979 ) { fLumiPerRun = 1.0995; }
  else if ( fRunNum == 297029 ) { fLumiPerRun = 7.2370; }
  else if ( fRunNum == 297031 ) { fLumiPerRun = 6.0499; }
  else if ( fRunNum == 297035 ) { fLumiPerRun = 0.5705; }
  else if ( fRunNum == 297085 ) { fLumiPerRun = 0.9774; }
  else if ( fRunNum == 297117 ) { fLumiPerRun = 2.3096; }
  else if ( fRunNum == 297118 ) { fLumiPerRun = 2.43;   }
  else if ( fRunNum == 297119 ) { fLumiPerRun = 2.6870; }
  else if ( fRunNum == 297123 ) { fLumiPerRun = 3.2804; }
  else if ( fRunNum == 297124 ) { fLumiPerRun = 0.6395; }
  else if ( fRunNum == 297128 ) { fLumiPerRun = 2.411;  }
  else if ( fRunNum == 297129 ) { fLumiPerRun = 2.8300; }
  else if ( fRunNum == 297132 ) { fLumiPerRun = 2.8179; }
  else if ( fRunNum == 297133 ) { fLumiPerRun = 1.1454; }
  else if ( fRunNum == 297193 ) { fLumiPerRun = 7.5602; }
  else if ( fRunNum == 297194 ) { fLumiPerRun = 8.8428; }
  else if ( fRunNum == 297196 ) { fLumiPerRun = 2.1255; }
  else if ( fRunNum == 297218 ) { fLumiPerRun = 6.42;   }
  else if ( fRunNum == 297219 ) { fLumiPerRun = 10.531; }
  else if ( fRunNum == 297221 ) { fLumiPerRun = 2.8309; }
  else if ( fRunNum == 297222 ) { fLumiPerRun = 1.7175; }
  else if ( fRunNum == 297278 ) { fLumiPerRun = 0.6019; }
  else if ( fRunNum == 297310 ) { fLumiPerRun = 0.6701; }
  else if ( fRunNum == 297312 ) { fLumiPerRun = 2.4002; }
  else if ( fRunNum == 297315 ) { fLumiPerRun = 7.8271; }
  else if ( fRunNum == 297317 ) { fLumiPerRun = 4.3148; }
  else if ( fRunNum == 297363 ) { fLumiPerRun = 1.9122; }
  else if ( fRunNum == 297366 ) { fLumiPerRun = 2.1293; }
  else if ( fRunNum == 297367 ) { fLumiPerRun = 3.1548; }
  else if ( fRunNum == 297372 ) { fLumiPerRun = 3.2003; }
  else if ( fRunNum == 297379 ) { fLumiPerRun = 6.8050; }
  else if ( fRunNum == 297380 ) { fLumiPerRun = 1.5488; }
  else if ( fRunNum == 297405 ) { fLumiPerRun = 0.6007; }
  else if ( fRunNum == 297408 ) { fLumiPerRun = 4.1021; }
  else if ( fRunNum == 297413 ) { fLumiPerRun = 2.9907; }
  else if ( fRunNum == 297414 ) { fLumiPerRun = 2.2140; }
  else if ( fRunNum == 297415 ) { fLumiPerRun = 6.8227; }
  else if ( fRunNum == 297441 ) { fLumiPerRun = 5.0556; }
  else if ( fRunNum == 297442 ) { fLumiPerRun = 1.9878; }
  else if ( fRunNum == 297446 ) { fLumiPerRun = 8.1326; }
  else if ( fRunNum == 297450 ) { fLumiPerRun = 1.9518; }
  else if ( fRunNum == 297451 ) { fLumiPerRun = 1.3327; }
  else if ( fRunNum == 297452 ) { fLumiPerRun = 1.1512; }
  else if ( fRunNum == 297479 ) { fLumiPerRun = 7.7463; }
  else if ( fRunNum == 297481 ) { fLumiPerRun = 10.645; }
  else if ( fRunNum == 297483 ) { fLumiPerRun = 1.9505; }
  else if ( fRunNum == 297512 ) { fLumiPerRun = 1.5848; }
  else if ( fRunNum == 297537 ) { fLumiPerRun = 1.8096; }
  else if ( fRunNum == 297540 ) { fLumiPerRun = 0.6286; }
  else if ( fRunNum == 297541 ) { fLumiPerRun = 4.0120; }
  else if ( fRunNum == 297542 ) { fLumiPerRun = 1.5362; }
  else if ( fRunNum == 297544 ) { fLumiPerRun = 7.2900; }
  else if ( fRunNum == 297558 ) { fLumiPerRun = 0.4783; }
  else if ( fRunNum == 297588 ) { fLumiPerRun = 5.2912; }
  else if ( fRunNum == 297590 ) { fLumiPerRun = 3.06;   }
  //_______________________________
  //_______________________________
  /* -
   * -
   * - CHECKAD ANALYSIS
   * -
   */
  // else if ( fRunNum == 295585 ) { fLumiPerRun = 0.0793352; }
  // else if ( fRunNum == 295586 ) { fLumiPerRun = 0.238599; }
  // else if ( fRunNum == 295587 ) { fLumiPerRun = 0.109518; }
  // else if ( fRunNum == 295588 ) { fLumiPerRun = 0.135709; }
  // else if ( fRunNum == 295589 ) { fLumiPerRun = 0.281897; }
  // else if ( fRunNum == 295612 ) { fLumiPerRun = 0.448985; }
  // else if ( fRunNum == 295615 ) { fLumiPerRun = 0.0565828; }
  // else if ( fRunNum == 295665 ) { fLumiPerRun = 0.334733; }
  // else if ( fRunNum == 295666 ) { fLumiPerRun = 0.323941; }
  // else if ( fRunNum == 295667 ) { fLumiPerRun = 0.0970128; }
  // else if ( fRunNum == 295668 ) { fLumiPerRun = 0.130088; }
  // else if ( fRunNum == 295671 ) { fLumiPerRun = 0.325985; }
  // else if ( fRunNum == 295673 ) { fLumiPerRun = 0.312556; }
  // else if ( fRunNum == 295675 ) { fLumiPerRun = 0.13204; }
  // else if ( fRunNum == 295676 ) { fLumiPerRun = 0.321284; }
  // else if ( fRunNum == 295677 ) { fLumiPerRun = 0.265538; }
  // else if ( fRunNum == 295714 ) { fLumiPerRun = 0.345554; }
  // else if ( fRunNum == 295716 ) { fLumiPerRun = 0.33778; }
  // else if ( fRunNum == 295717 ) { fLumiPerRun = 0.287941; }
  // else if ( fRunNum == 295718 ) { fLumiPerRun = 0.257395; }
  // else if ( fRunNum == 295719 ) { fLumiPerRun = 0.294713; }
  // else if ( fRunNum == 295723 ) { fLumiPerRun = 0.506379; }
  // else if ( fRunNum == 295725 ) { fLumiPerRun = 0.8453; }
  // else if ( fRunNum == 295753 ) { fLumiPerRun = 0.384579; }
  // else if ( fRunNum == 295754 ) { fLumiPerRun = 1.12405; }
  // else if ( fRunNum == 295755 ) { fLumiPerRun = 1.01645; }
  // else if ( fRunNum == 295758 ) { fLumiPerRun = 1.62041; }
  // else if ( fRunNum == 295759 ) { fLumiPerRun = 0.532728; }
  // else if ( fRunNum == 295762 ) { fLumiPerRun = 0.274725; }
  // else if ( fRunNum == 295763 ) { fLumiPerRun = 1.0236; }
  // else if ( fRunNum == 295786 ) { fLumiPerRun = 0.746701; }
  // else if ( fRunNum == 295788 ) { fLumiPerRun = 3.06545; }
  // else if ( fRunNum == 295791 ) { fLumiPerRun = 0.85713; }
  // else if ( fRunNum == 295816 ) { fLumiPerRun = 1.20498; }
  // else if ( fRunNum == 295818 ) { fLumiPerRun = 0.145167; }
  // else if ( fRunNum == 295819 ) { fLumiPerRun = 2.74121; }
  // else if ( fRunNum == 295822 ) { fLumiPerRun = 2.2132; }
  // else if ( fRunNum == 295825 ) { fLumiPerRun = 0.255836; }
  // else if ( fRunNum == 295826 ) { fLumiPerRun = 1.48001; }
  // else if ( fRunNum == 295829 ) { fLumiPerRun = 0.934546; }
  // else if ( fRunNum == 295831 ) { fLumiPerRun = 0.779334; }
  // else if ( fRunNum == 295854 ) { fLumiPerRun = 1.73457; }
  // else if ( fRunNum == 295855 ) { fLumiPerRun = 1.74464; }
  // else if ( fRunNum == 295856 ) { fLumiPerRun = 1.73755; }
  // else if ( fRunNum == 295859 ) { fLumiPerRun = 1.05037; }
  // else if ( fRunNum == 295860 ) { fLumiPerRun = 0.834142; }
  // else if ( fRunNum == 295861 ) { fLumiPerRun = 1.06703; }
  // else if ( fRunNum == 295863 ) { fLumiPerRun = 0.728227; }
  // else if ( fRunNum == 295881 ) { fLumiPerRun = 0.711494; }
  // else if ( fRunNum == 295908 ) { fLumiPerRun = 2.92589; }
  // else if ( fRunNum == 295909 ) { fLumiPerRun = 0.787373; }
  // else if ( fRunNum == 295910 ) { fLumiPerRun = 3.19503; }
  // else if ( fRunNum == 295913 ) { fLumiPerRun = 3.12937; }
  // else if ( fRunNum == 295936 ) { fLumiPerRun = 1.47357; }
  // else if ( fRunNum == 295937 ) { fLumiPerRun = 0.405657; }
  // else if ( fRunNum == 295941 ) { fLumiPerRun = 1.67616; }
  // else if ( fRunNum == 295942 ) { fLumiPerRun = 1.92337; }
  // else if ( fRunNum == 295943 ) { fLumiPerRun = 1.67424; }
  // else if ( fRunNum == 295945 ) { fLumiPerRun = 2.03661; }
  // else if ( fRunNum == 295947 ) { fLumiPerRun = 2.69035; }
  // else if ( fRunNum == 296061 ) { fLumiPerRun = 1.29676; }
  // else if ( fRunNum == 296062 ) { fLumiPerRun = 1.80833; }
  // else if ( fRunNum == 296063 ) { fLumiPerRun = 2.69974; }
  // else if ( fRunNum == 296065 ) { fLumiPerRun = 2.45097; }
  // else if ( fRunNum == 296066 ) { fLumiPerRun = 0.733886; }
  // else if ( fRunNum == 296068 ) { fLumiPerRun = 1.98122; }
  // else if ( fRunNum == 296123 ) { fLumiPerRun = 0.506486; }
  // else if ( fRunNum == 296128 ) { fLumiPerRun = 0.445452; }
  // else if ( fRunNum == 296132 ) { fLumiPerRun = 2.10822; }
  // else if ( fRunNum == 296133 ) { fLumiPerRun = 1.73213; }
  // else if ( fRunNum == 296134 ) { fLumiPerRun = 3.91041; }
  // else if ( fRunNum == 296135 ) { fLumiPerRun = 2.3412; }
  // else if ( fRunNum == 296142 ) { fLumiPerRun = 1.78935; }
  // else if ( fRunNum == 296143 ) { fLumiPerRun = 0.534028; }
  // else if ( fRunNum == 296191 ) { fLumiPerRun = 5.051; }
  // else if ( fRunNum == 296192 ) { fLumiPerRun = 0.497364; }
  // else if ( fRunNum == 296194 ) { fLumiPerRun = 2.84426; }
  // else if ( fRunNum == 296195 ) { fLumiPerRun = 0.737647; }
  // else if ( fRunNum == 296196 ) { fLumiPerRun = 2.42215; }
  // else if ( fRunNum == 296197 ) { fLumiPerRun = 2.07279; }
  // else if ( fRunNum == 296198 ) { fLumiPerRun = 0.813921; }
  // else if ( fRunNum == 296241 ) { fLumiPerRun = 0.845868; }
  // else if ( fRunNum == 296242 ) { fLumiPerRun = 0.95166; }
  // else if ( fRunNum == 296243 ) { fLumiPerRun = 1.56742; }
  // else if ( fRunNum == 296244 ) { fLumiPerRun = 8.37179; }
  // else if ( fRunNum == 296246 ) { fLumiPerRun = 1.82994; }
  // else if ( fRunNum == 296247 ) { fLumiPerRun = 1.1763; }
  // else if ( fRunNum == 296269 ) { fLumiPerRun = 3.8392; }
  // else if ( fRunNum == 296270 ) { fLumiPerRun = 1.51158; }
  // else if ( fRunNum == 296273 ) { fLumiPerRun = 7.23985; }
  // else if ( fRunNum == 296279 ) { fLumiPerRun = 0.405692; }
  // else if ( fRunNum == 296280 ) { fLumiPerRun = 1.61122; }
  // else if ( fRunNum == 296303 ) { fLumiPerRun = 2.01567; }
  // else if ( fRunNum == 296304 ) { fLumiPerRun = 6.11939; }
  // else if ( fRunNum == 296307 ) { fLumiPerRun = 2.93104; }
  // else if ( fRunNum == 296309 ) { fLumiPerRun = 2.10255; }
  // else if ( fRunNum == 296312 ) { fLumiPerRun = 2.12718; }
  // else if ( fRunNum == 296377 ) { fLumiPerRun = 6.08141; }
  // else if ( fRunNum == 296378 ) { fLumiPerRun = 5.34678; }
  // else if ( fRunNum == 296379 ) { fLumiPerRun = 2.09698; }
  // else if ( fRunNum == 296380 ) { fLumiPerRun = 2.88934; }
  // else if ( fRunNum == 296381 ) { fLumiPerRun = 1.44081; }
  // else if ( fRunNum == 296383 ) { fLumiPerRun = 1.51276; }
  // else if ( fRunNum == 296414 ) { fLumiPerRun = 4.8759; }
  // else if ( fRunNum == 296419 ) { fLumiPerRun = 2.74845; }
  // else if ( fRunNum == 296420 ) { fLumiPerRun = 1.40973; }
  // else if ( fRunNum == 296423 ) { fLumiPerRun = 1.57873; }
  // else if ( fRunNum == 296424 ) { fLumiPerRun = 0.386229; }
  // else if ( fRunNum == 296433 ) { fLumiPerRun = 4.46674; }
  // else if ( fRunNum == 296472 ) { fLumiPerRun = 0.862237; }
  // else if ( fRunNum == 296509 ) { fLumiPerRun = 2.99727; }
  // else if ( fRunNum == 296510 ) { fLumiPerRun = 9.07207; }
  // else if ( fRunNum == 296511 ) { fLumiPerRun = 2.5694; }
  // else if ( fRunNum == 296514 ) { fLumiPerRun = 0.48925; }
  // else if ( fRunNum == 296516 ) { fLumiPerRun = 0.613001; }
  // else if ( fRunNum == 296547 ) { fLumiPerRun = 1.08314; }
  // else if ( fRunNum == 296548 ) { fLumiPerRun = 1.37636; }
  // else if ( fRunNum == 296549 ) { fLumiPerRun = 4.86114; }
  // else if ( fRunNum == 296550 ) { fLumiPerRun = 3.96203; }
  // else if ( fRunNum == 296551 ) { fLumiPerRun = 2.02168; }
  // else if ( fRunNum == 296552 ) { fLumiPerRun = 0.483545; }
  // else if ( fRunNum == 296553 ) { fLumiPerRun = 0.707355; }
  // else if ( fRunNum == 296615 ) { fLumiPerRun = 1.56721; }
  // else if ( fRunNum == 296616 ) { fLumiPerRun = 0.53916; }
  // else if ( fRunNum == 296618 ) { fLumiPerRun = 2.28102; }
  // else if ( fRunNum == 296619 ) { fLumiPerRun = 1.5574; }
  // else if ( fRunNum == 296622 ) { fLumiPerRun = 0.703628; }
  // else if ( fRunNum == 296623 ) { fLumiPerRun = 2.16132; }
  // else if ( fRunNum == 296690 ) { fLumiPerRun = 6.55512; }
  // else if ( fRunNum == 296691 ) { fLumiPerRun = 0.651063; }
  // else if ( fRunNum == 296694 ) { fLumiPerRun = 5.02114; }
  // else if ( fRunNum == 296749 ) { fLumiPerRun = 9.27577; }
  // else if ( fRunNum == 296750 ) { fLumiPerRun = 8.21552; }
  // else if ( fRunNum == 296781 ) { fLumiPerRun = 0.817883; }
  // else if ( fRunNum == 296784 ) { fLumiPerRun = 3.0019; }
  // else if ( fRunNum == 296785 ) { fLumiPerRun = 1.9085; }
  // else if ( fRunNum == 296786 ) { fLumiPerRun = 0.753734; }
  // else if ( fRunNum == 296787 ) { fLumiPerRun = 3.24915; }
  // else if ( fRunNum == 296791 ) { fLumiPerRun = 0.757278; }
  // else if ( fRunNum == 296793 ) { fLumiPerRun = 1.39522; }
  // else if ( fRunNum == 296794 ) { fLumiPerRun = 3.13783; }
  // else if ( fRunNum == 296799 ) { fLumiPerRun = 2.77608; }
  // else if ( fRunNum == 296836 ) { fLumiPerRun = 1.51136; }
  // else if ( fRunNum == 296838 ) { fLumiPerRun = 0.542885; }
  // else if ( fRunNum == 296839 ) { fLumiPerRun = 2.95419; }
  // else if ( fRunNum == 296848 ) { fLumiPerRun = 2.19145; }
  // else if ( fRunNum == 296849 ) { fLumiPerRun = 11.8282; }
  // else if ( fRunNum == 296850 ) { fLumiPerRun = 2.8369; }
  // else if ( fRunNum == 296851 ) { fLumiPerRun = 0.890021; }
  // else if ( fRunNum == 296852 ) { fLumiPerRun = 0.948014; }
  // else if ( fRunNum == 296890 ) { fLumiPerRun = 6.96187; }
  // else if ( fRunNum == 296894 ) { fLumiPerRun = 4.47163; }
  // else if ( fRunNum == 296899 ) { fLumiPerRun = 2.11682; }
  // else if ( fRunNum == 296900 ) { fLumiPerRun = 2.78312; }
  // else if ( fRunNum == 296903 ) { fLumiPerRun = 1.03903; }
  // else if ( fRunNum == 296930 ) { fLumiPerRun = 1.42317; }
  // else if ( fRunNum == 296931 ) { fLumiPerRun = 0.528865; }
  // else if ( fRunNum == 296932 ) { fLumiPerRun = 1.19931; }
  // else if ( fRunNum == 296934 ) { fLumiPerRun = 2.62967; }
  // else if ( fRunNum == 296935 ) { fLumiPerRun = 4.49271; }
  // else if ( fRunNum == 296938 ) { fLumiPerRun = 1.65053; }
  // else if ( fRunNum == 296941 ) { fLumiPerRun = 2.93132; }
  // else if ( fRunNum == 296966 ) { fLumiPerRun = 3.37409; }
  // else if ( fRunNum == 296967 ) { fLumiPerRun = 0.804596; }
  // else if ( fRunNum == 296968 ) { fLumiPerRun = 3.18308; }
  // else if ( fRunNum == 296969 ) { fLumiPerRun = 1.88784; }
  // else if ( fRunNum == 296971 ) { fLumiPerRun = 0.6895; }
  // else if ( fRunNum == 296975 ) { fLumiPerRun = 7.47159; }
  // else if ( fRunNum == 296976 ) { fLumiPerRun = 1.11722; }
  // else if ( fRunNum == 296979 ) { fLumiPerRun = 1.09943; }
  // else if ( fRunNum == 297029 ) { fLumiPerRun = 7.20294; }
  // else if ( fRunNum == 297031 ) { fLumiPerRun = 5.94515; }
  // else if ( fRunNum == 297035 ) { fLumiPerRun = 1.30617; }
  // else if ( fRunNum == 297085 ) { fLumiPerRun = 0.977753; }
  // else if ( fRunNum == 297117 ) { fLumiPerRun = 2.34892; }
  // else if ( fRunNum == 297118 ) { fLumiPerRun = 2.44282; }
  // else if ( fRunNum == 297119 ) { fLumiPerRun = 2.68704; }
  // else if ( fRunNum == 297123 ) { fLumiPerRun = 3.3714; }
  // else if ( fRunNum == 297124 ) { fLumiPerRun = 0.639463; }
  // else if ( fRunNum == 297128 ) { fLumiPerRun = 2.41074; }
  // else if ( fRunNum == 297129 ) { fLumiPerRun = 2.82995; }
  // else if ( fRunNum == 297132 ) { fLumiPerRun = 2.81789; }
  // else if ( fRunNum == 297133 ) { fLumiPerRun = 1.16976; }
  // else if ( fRunNum == 297193 ) { fLumiPerRun = 7.64123; }
  // else if ( fRunNum == 297194 ) { fLumiPerRun = 9.86729; }
  // else if ( fRunNum == 297196 ) { fLumiPerRun = 2.1255; }
  // else if ( fRunNum == 297218 ) { fLumiPerRun = 6.39259; }
  // else if ( fRunNum == 297219 ) { fLumiPerRun = 9.29989; }
  // else if ( fRunNum == 297221 ) { fLumiPerRun = 2.83193; }
  // else if ( fRunNum == 297222 ) { fLumiPerRun = 1.69325; }
  // else if ( fRunNum == 297278 ) { fLumiPerRun = 0.601609; }
  // else if ( fRunNum == 297310 ) { fLumiPerRun = 0.670071; }
  // else if ( fRunNum == 297312 ) { fLumiPerRun = 2.40205; }
  // else if ( fRunNum == 297315 ) { fLumiPerRun = 7.93229; }
  // else if ( fRunNum == 297317 ) { fLumiPerRun = 4.31559; }
  // else if ( fRunNum == 297363 ) { fLumiPerRun = 1.89669; }
  // else if ( fRunNum == 297366 ) { fLumiPerRun = 2.05394; }
  // else if ( fRunNum == 297367 ) { fLumiPerRun = 3.11285; }
  // else if ( fRunNum == 297372 ) { fLumiPerRun = 3.22421; }
  // else if ( fRunNum == 297379 ) { fLumiPerRun = 6.92989; }
  // else if ( fRunNum == 297380 ) { fLumiPerRun = 1.46125; }
  // else if ( fRunNum == 297405 ) { fLumiPerRun = 0.51899; }
  // else if ( fRunNum == 297408 ) { fLumiPerRun = 4.05969; }
  // else if ( fRunNum == 297413 ) { fLumiPerRun = 2.99189; }
  // else if ( fRunNum == 297414 ) { fLumiPerRun = 2.21433; }
  // else if ( fRunNum == 297415 ) { fLumiPerRun = 6.76725; }
  // else if ( fRunNum == 297441 ) { fLumiPerRun = 4.92805; }
  // else if ( fRunNum == 297442 ) { fLumiPerRun = 1.98748; }
  // else if ( fRunNum == 297446 ) { fLumiPerRun = 8.1332; }
  // else if ( fRunNum == 297450 ) { fLumiPerRun = 1.95205; }
  // else if ( fRunNum == 297451 ) { fLumiPerRun = 1.33275; }
  // else if ( fRunNum == 297452 ) { fLumiPerRun = 1.15124; }
  // else if ( fRunNum == 297479 ) { fLumiPerRun = 7.71966; }
  // else if ( fRunNum == 297481 ) { fLumiPerRun = 9.55402; }
  // else if ( fRunNum == 297483 ) { fLumiPerRun = 1.99555; }
  // else if ( fRunNum == 297512 ) { fLumiPerRun = 1.5838; }
  // else if ( fRunNum == 297537 ) { fLumiPerRun = 1.80636; }
  // else if ( fRunNum == 297540 ) { fLumiPerRun = 0.628094; }
  // else if ( fRunNum == 297541 ) { fLumiPerRun = 3.96702; }
  // else if ( fRunNum == 297542 ) { fLumiPerRun = 1.5516; }
  // else if ( fRunNum == 297544 ) { fLumiPerRun = 7.29158; }
  // else if ( fRunNum == 297558 ) { fLumiPerRun = 0.478978; }
  // else if ( fRunNum == 297588 ) { fLumiPerRun = 5.26723; }
  // else if ( fRunNum == 297590 ) { fLumiPerRun = 3.14808; }
  else if ( fRunNum == 297595 ) { fLumiPerRun = 0.0; }


}
