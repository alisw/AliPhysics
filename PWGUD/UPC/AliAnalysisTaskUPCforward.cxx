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
#include "AliAODVertex.h"         // My addition, to use Eugeny Krishen's format


// my headers
#include "AliAnalysisTaskUPCforward.h"



class AliAnalysisTaskUPCforward;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskUPCforward) // classimp: necessary for root

//_____________________________________________________________________________
AliAnalysisTaskUPCforward::AliAnalysisTaskUPCforward()
    : AliAnalysisTaskSE(),
      fAOD(0),
      fOutputList(0),
      fSetSingleMuonPt(0),
      fPtSingleMuonH(0),
      fNumberMuonsH(0),
      fCounterH(0),
      fEtaMuonH(0),
      fRAbsMuonH(0),
      fInvariantMassDistributionH(0),
      fInvariantMassDistributionAtDcaH(0),
      fEntriesAgainstRunNumberH(0),
      fEntriesAgainstRunNumberProperlyH(0),
      fRunNumberTriggerCMUP11ClassH(0),
      fRunNumberTriggerCMUP11ClassProperlyH(0),
      fRunNumberTriggerCMUP26ClassH(0),
      fRunNumberTriggerCMUP26ClassProperlyH(0),
      fRunNumberTriggerCMUP6ClassH(0),
      fRunNumberTriggerCMUP6ClassProperlyH(0),
      fRunNumberTriggerCMUP10ClassH(0),
      fRunNumberTriggerCMUP10ClassProperlyH(0),
      fRunNumberTriggerCMUP13ClassH(0),
      fRunNumberTriggerCMUP13ClassProperlyH(0),
      fTriggersVsRunH(0),
      fInvariantMassDistributionCoherentH(0),
      fInvariantMassDistributionCoherentRapidityBinsH{ 0, 0, 0, 0, 0, 0 },
      fInvariantMassDistributionCoherentShiftMinusTwoH(0),
      fInvariantMassDistributionCoherentShiftMinusOneH(0),
      fInvariantMassDistributionCoherentShiftPlusOneH(0),
      fInvariantMassDistributionCoherentShiftPlusTwoH(0),
      fInvariantMassDistributionIncoherentH(0),
      fInvariantMassDistributionIncoherentShiftMinusTwoH(0),
      fInvariantMassDistributionIncoherentShiftMinusOneH(0),
      fInvariantMassDistributionIncoherentShiftPlusOneH(0),
      fInvariantMassDistributionIncoherentShiftPlusTwoH(0),
      fDimuonPtDistributionH(0),
      fDimuonPtDistributionShiftPlusOneH(0),
      fZNCEnergyAgainstEntriesH(0),
      fZNAEnergyAgainstEntriesH(0),
      fZNCEnergyBeforeTimingSelectionH(0),
      fZNAEnergyBeforeTimingSelectionH(0),
      fZNCEnergyAgainstEntriesExtendedH(0),
      fZNAEnergyAgainstEntriesExtendedH(0),
      fZNCEnergyAgainstEntriesExtendedHv2(0),
      fZNAEnergyAgainstEntriesExtendedHv2(0),
      fZNCEnergyBeforeTimingSelectionExtendedH(0),
      fZNAEnergyBeforeTimingSelectionExtendedH(0),
      fZNCEnergyCalibratedH(0),
      fZNAEnergyCalibratedH(0),
      fZNCEnergyUncalibratedH(0),
      fZNAEnergyUncalibratedH(0),
      fZNCEnergyCalibratedHigherGainH(0),
      fZNAEnergyCalibratedHigherGainH(0),
      fZNCTimeAgainstEntriesH(0),
      fZNATimeAgainstEntriesH(0),
      fZNCTimeStrictTimeWindowH(0),
      fZNATimeStrictTimeWindowH(0),
      fZNCTimeWithoutTimingH{0, 0, 0, 0},
      fZNATimeWithoutTimingH{0, 0, 0, 0},
      fZNCTime4FillingH(0),
      fZNATime4FillingH(0),
      fZNCminusZNAtimeVsZNCplusZNAtimeH{0, 0, 0, 0},
      fZNCminusZNAtimeVsZNCplusZNAtime4FillingH(0),
      fCounterZNCH(0),
      fCounterZNAH(0),
      fInvariantMassDistributionNoNeutronsH(0),
      fInvariantMassDistributionOneNeutronH(0),
      fInvariantMassDistributionAtLeastOneNeutronH(0),
      fInvariantMassDistributionCoherentNoNeutronsH(0),
      fInvariantMassDistributionCoherentOneNeutronH(0),
      fInvariantMassDistributionCoherentAtLeastOneNeutronH(0),
      fInvariantMassDistributionIncoherentNoNeutronsH(0),
      fInvariantMassDistributionIncoherentOneNeutronH(0),
      fInvariantMassDistributionIncoherentAtLeastOneNeutronH(0),
      fDcaAgainstInvariantMassH(0),
      fInvariantMassDistributionExtendedH(0),
      fInvariantMassDistributionCoherentExtendedH(0),
      fInvariantMassDistributionIncoherentExtendedH(0),
      fMuonTrackCuts(0x0),
      fRunNum(0),
      fTracklets(0),
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
      fInvariantMassDistributionCoherentZNCzeroZNAzeroH(0),
      fInvariantMassDistributionCoherentZNCzeroZNAzeroShiftMinusTwoH(0),
      fInvariantMassDistributionCoherentZNCzeroZNAzeroShiftMinusOneH(0),
      fInvariantMassDistributionCoherentZNCzeroZNAzeroShiftPlusOneH(0),
      fInvariantMassDistributionCoherentZNCzeroZNAzeroShiftPlusTwoH(0),
      fInvariantMassDistributionCoherentZNCzeroZNAzeroHv2(0),
      fInvariantMassDistributionCoherentZNCzeroZNAzeroRapidityHv2{0, 0, 0},
      fInvariantMassDistributionCoherentZNCzeroZNAanyH(0),
      fInvariantMassDistributionCoherentZNCzeroZNAanyShiftMinusTwoH(0),
      fInvariantMassDistributionCoherentZNCzeroZNAanyShiftMinusOneH(0),
      fInvariantMassDistributionCoherentZNCzeroZNAanyShiftPlusOneH(0),
      fInvariantMassDistributionCoherentZNCzeroZNAanyShiftPlusTwoH(0),
      fInvariantMassDistributionCoherentZNCzeroZNAanyHv2(0),
      fInvariantMassDistributionCoherentZNCzeroZNAanyRapidityHv2{0, 0, 0},
      fInvariantMassDistributionCoherentZNCanyZNAzeroH(0),
      fInvariantMassDistributionCoherentZNCanyZNAzeroShiftMinusTwoH(0),
      fInvariantMassDistributionCoherentZNCanyZNAzeroShiftMinusOneH(0),
      fInvariantMassDistributionCoherentZNCanyZNAzeroShiftPlusOneH(0),
      fInvariantMassDistributionCoherentZNCanyZNAzeroShiftPlusTwoH(0),
      fInvariantMassDistributionCoherentZNCanyZNAzeroHv2(0),
      fInvariantMassDistributionCoherentZNCanyZNAzeroRapidityHv2{0, 0, 0},
      fInvariantMassDistributionCoherentZNCanyZNAanyH(0),
      fInvariantMassDistributionCoherentZNCanyZNAanyShiftMinusTwoH(0),
      fInvariantMassDistributionCoherentZNCanyZNAanyShiftMinusOneH(0),
      fInvariantMassDistributionCoherentZNCanyZNAanyShiftPlusOneH(0),
      fInvariantMassDistributionCoherentZNCanyZNAanyShiftPlusTwoH(0),
      fInvariantMassDistributionCoherentZNCanyZNAanyHv2(0),
      fInvariantMassDistributionCoherentZNCanyZNAanyRapidityHv2{0, 0, 0},
      fInvariantMassDistributionIncoherentZNCzeroZNAzeroH(0),
      fInvariantMassDistributionIncoherentZNCzeroZNAzeroShiftMinusTwoH(0),
      fInvariantMassDistributionIncoherentZNCzeroZNAzeroShiftMinusOneH(0),
      fInvariantMassDistributionIncoherentZNCzeroZNAzeroShiftPlusOneH(0),
      fInvariantMassDistributionIncoherentZNCzeroZNAzeroShiftPlusTwoH(0),
      fInvariantMassDistributionIncoherentZNCzeroZNAzeroHv2(0),
      fInvariantMassDistributionIncoherentZNCzeroZNAzeroRapidityHv2{0, 0, 0},
      fInvariantMassDistributionIncoherentZNCzeroZNAanyH(0),
      fInvariantMassDistributionIncoherentZNCzeroZNAanyShiftMinusTwoH(0),
      fInvariantMassDistributionIncoherentZNCzeroZNAanyShiftMinusOneH(0),
      fInvariantMassDistributionIncoherentZNCzeroZNAanyShiftPlusOneH(0),
      fInvariantMassDistributionIncoherentZNCzeroZNAanyShiftPlusTwoH(0),
      fInvariantMassDistributionIncoherentZNCzeroZNAanyHv2(0),
      fInvariantMassDistributionIncoherentZNCzeroZNAanyRapidityHv2{0, 0, 0},
      fInvariantMassDistributionIncoherentZNCanyZNAzeroH(0),
      fInvariantMassDistributionIncoherentZNCanyZNAzeroShiftMinusTwoH(0),
      fInvariantMassDistributionIncoherentZNCanyZNAzeroShiftMinusOneH(0),
      fInvariantMassDistributionIncoherentZNCanyZNAzeroShiftPlusOneH(0),
      fInvariantMassDistributionIncoherentZNCanyZNAzeroShiftPlusTwoH(0),
      fInvariantMassDistributionIncoherentZNCanyZNAzeroHv2(0),
      fInvariantMassDistributionIncoherentZNCanyZNAzeroRapidityHv2{0, 0, 0},
      fInvariantMassDistributionIncoherentZNCanyZNAanyH(0),
      fInvariantMassDistributionIncoherentZNCanyZNAanyShiftMinusTwoH(0),
      fInvariantMassDistributionIncoherentZNCanyZNAanyShiftMinusOneH(0),
      fInvariantMassDistributionIncoherentZNCanyZNAanyShiftPlusOneH(0),
      fInvariantMassDistributionIncoherentZNCanyZNAanyShiftPlusTwoH(0),
      fInvariantMassDistributionIncoherentZNCanyZNAanyHv2(0),
      fInvariantMassDistributionIncoherentZNCanyZNAanyRapidityHv2{0, 0, 0},
      fDimuonPtDistributionZNCzeroZNAzeroH(0),
      fDimuonPtDistributionZNCzeroZNAanyH(0),
      fDimuonPtDistributionZNCanyZNAzeroH(0),
      fDimuonPtDistributionZNCanyZNAanyH(0),
      fDimuonPtDistributionZNCzeroZNAzeroShiftPlusOneH(0),
      fDimuonPtDistributionZNCzeroZNAanyShiftPlusOneH(0),
      fDimuonPtDistributionZNCanyZNAzeroShiftPlusOneH(0),
      fDimuonPtDistributionZNCanyZNAanyShiftPlusOneH(0),
      fDimuonPtDistributionZNCzeroZNAzeroHv2(0),
      fDimuonPtDistributionZNCzeroZNAanyHv2(0),
      fDimuonPtDistributionZNCanyZNAzeroHv2(0),
      fDimuonPtDistributionZNCanyZNAanyHv2(0),
      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv2{0, 0, 0},
      fDimuonPtDistributionZNCzeroZNAanyRapidityHv2{0, 0, 0},
      fDimuonPtDistributionZNCanyZNAzeroRapidityHv2{0, 0, 0},
      fDimuonPtDistributionZNCanyZNAanyRapidityHv2{0, 0, 0},
      fDimuonPtDistributionCoherentZNCzeroZNAzeroH(0),
      fDimuonPtDistributionCoherentZNCzeroZNAanyH(0),
      fDimuonPtDistributionCoherentZNCanyZNAzeroH(0),
      fDimuonPtDistributionCoherentZNCanyZNAanyH(0),
      fDimuonPtDistributionIncoherentZNCzeroZNAzeroH(0),
      fDimuonPtDistributionIncoherentZNCzeroZNAanyH(0),
      fDimuonPtDistributionIncoherentZNCanyZNAzeroH(0),
      fDimuonPtDistributionIncoherentZNCanyZNAanyH(0),
      fAngularDistribOfPositiveMuonRestFrameJPsiH(0),
      fAngularDistribOfNegativeMuonRestFrameJPsiH(0),
      fCheckHelicityRestFrameJPsiH(0),
      fThetaDistribOfPositiveMuonRestFrameJPsiRapidityBinH{ 0, 0, 0, 0, 0, 0, 0, 0},
      fCosThetaHelicityFrameJPsiH(0),
      fCosThetaCollinsSoperFrameJPsiH(0),
      fPhiHelicityFrameJPsiH(0),
      fPhiCollinsSoperFrameJPsiH(0),
      fCosThetaHelicityFrameJPsiRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0},
      fCosThetaCollinsSoperFrameJPsiRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0},
      fPhiHelicityFrameJPsiRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0},
      fPhiCollinsSoperFrameJPsiRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0},
      fCosThetaHelicityFrameJPsiTenRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      fCosThetaCollinsSoperFrameJPsiTenRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      fPhiHelicityFrameJPsiTenRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      fPhiCollinsSoperFrameJPsiTenRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      fInvariantMassDistributionInBinsOfCosThetaHelicityFrameH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      fInvariantMassDistributionBinsOfCosThetaAndPhiHelicityFrameH(0),
      fCosThetaAndPhiHelicityFrameInclusivePeopleBinningH(0),
      fInvariantMassDistributionForSignalExtractionHelicityFrameH(0),
      fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameTryingH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameTryingH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                                0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                                0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                                0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                                0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fInvariantMassDistributionForSignalExtractionHelicityFrameMyBinningH(0),
      fInvariantMassDistributionForSignalExtractionCsFrameMyBinningH(0),
      fInvariantMassDistributionForSignalExtractionHelicityFrameMyBinning0N0NH(0),
      fInvariantMassDistributionForSignalExtractionCsFrameMyBinning0N0NH(0),
      fInvariantMassDistributionStrictPtH(0),
      fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameMyVariableBinningH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                                                0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                                                0, 0, 0, 0, 0, 0 },
      fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameMyVariableBinningH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameMyStrictVariableBinningH{ 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameMyStrictVariableBinningH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameMySeventeenBinsVariableBinningH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fCosThetaHelicityFrameJPsiAlreadyCorrectedH(0),
      fDcaAgainstPtOfVectorMesonH(0),
      fInvariantMassDistributionStrictPtStrictDcaH(0),
      fInvariantMassDistributionVsPtStrictDcaH(0),
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
      fInvariantMassDistributionOnlyPhiHeFrameTwentyfiveBins0N0NH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                   0, 0, 0, 0, 0 },
      fInvariantMassDistributionOnlyCosThetaHeFrameTwentyfiveBins0N0NH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                        0, 0, 0, 0, 0 },
      fInvariantMassDistributionOnlyTildePhiHeFrameTwentyfiveBins0N0NH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                        0, 0, 0, 0, 0 },
      fInvariantMassDistributionOnlyPhiCsFrameTwentyfiveBins0N0NH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                   0, 0, 0, 0, 0 },
      fInvariantMassDistributionOnlyCosThetaCsFrameTwentyfiveBins0N0NH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                        0, 0, 0, 0, 0 },
      fInvariantMassDistributionOnlyTildePhiCsFrameTwentyfiveBins0N0NH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                        0, 0, 0, 0, 0 }
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}

//_____________________________________________________________________________
AliAnalysisTaskUPCforward::AliAnalysisTaskUPCforward( const char* name, Int_t _fSetSingleMuonPt = 0 )
    : AliAnalysisTaskSE(name),
      fAOD(0),
      fOutputList(0),
      fSetSingleMuonPt(_fSetSingleMuonPt),
      fPtSingleMuonH(0),
      fNumberMuonsH(0),
      fCounterH(0),
      fEtaMuonH(0),
      fRAbsMuonH(0),
      fInvariantMassDistributionH(0),
      fInvariantMassDistributionAtDcaH(0),
      fEntriesAgainstRunNumberH(0),
      fEntriesAgainstRunNumberProperlyH(0),
      fRunNumberTriggerCMUP11ClassH(0),
      fRunNumberTriggerCMUP11ClassProperlyH(0),
      fRunNumberTriggerCMUP26ClassH(0),
      fRunNumberTriggerCMUP26ClassProperlyH(0),
      fRunNumberTriggerCMUP6ClassH(0),
      fRunNumberTriggerCMUP6ClassProperlyH(0),
      fRunNumberTriggerCMUP10ClassH(0),
      fRunNumberTriggerCMUP10ClassProperlyH(0),
      fRunNumberTriggerCMUP13ClassH(0),
      fRunNumberTriggerCMUP13ClassProperlyH(0),
      fTriggersVsRunH(0),
      fInvariantMassDistributionCoherentH(0),
      fInvariantMassDistributionCoherentRapidityBinsH{ 0, 0, 0, 0, 0, 0 },
      fInvariantMassDistributionCoherentShiftMinusTwoH(0),
      fInvariantMassDistributionCoherentShiftMinusOneH(0),
      fInvariantMassDistributionCoherentShiftPlusOneH(0),
      fInvariantMassDistributionCoherentShiftPlusTwoH(0),
      fInvariantMassDistributionIncoherentH(0),
      fInvariantMassDistributionIncoherentShiftMinusTwoH(0),
      fInvariantMassDistributionIncoherentShiftMinusOneH(0),
      fInvariantMassDistributionIncoherentShiftPlusOneH(0),
      fInvariantMassDistributionIncoherentShiftPlusTwoH(0),
      fDimuonPtDistributionH(0),
      fZNCEnergyAgainstEntriesH(0),
      fZNAEnergyAgainstEntriesH(0),
      fZNCEnergyBeforeTimingSelectionH(0),
      fZNAEnergyBeforeTimingSelectionH(0),
      fZNCEnergyAgainstEntriesExtendedH(0),
      fZNAEnergyAgainstEntriesExtendedH(0),
      fZNCEnergyAgainstEntriesExtendedHv2(0),
      fZNAEnergyAgainstEntriesExtendedHv2(0),
      fZNCEnergyBeforeTimingSelectionExtendedH(0),
      fZNAEnergyBeforeTimingSelectionExtendedH(0),
      fZNCEnergyCalibratedH(0),
      fZNAEnergyCalibratedH(0),
      fZNCEnergyUncalibratedH(0),
      fZNAEnergyUncalibratedH(0),
      fZNCEnergyCalibratedHigherGainH(0),
      fZNAEnergyCalibratedHigherGainH(0),
      fZNCTimeAgainstEntriesH(0),
      fZNATimeAgainstEntriesH(0),
      fZNCTimeStrictTimeWindowH(0),
      fZNATimeStrictTimeWindowH(0),
      fZNCTimeWithoutTimingH{0, 0, 0, 0},
      fZNATimeWithoutTimingH{0, 0, 0, 0},
      fZNCTime4FillingH(0),
      fZNATime4FillingH(0),
      fZNCminusZNAtimeVsZNCplusZNAtimeH{0, 0, 0, 0},
      fZNCminusZNAtimeVsZNCplusZNAtime4FillingH(0),
      fCounterZNCH(0),
      fCounterZNAH(0),
      fInvariantMassDistributionNoNeutronsH(0),
      fInvariantMassDistributionOneNeutronH(0),
      fInvariantMassDistributionAtLeastOneNeutronH(0),
      fInvariantMassDistributionCoherentNoNeutronsH(0),
      fInvariantMassDistributionCoherentOneNeutronH(0),
      fInvariantMassDistributionCoherentAtLeastOneNeutronH(0),
      fInvariantMassDistributionIncoherentNoNeutronsH(0),
      fInvariantMassDistributionIncoherentOneNeutronH(0),
      fInvariantMassDistributionIncoherentAtLeastOneNeutronH(0),
      fDcaAgainstInvariantMassH(0),
      fInvariantMassDistributionExtendedH(0),
      fInvariantMassDistributionCoherentExtendedH(0),
      fInvariantMassDistributionIncoherentExtendedH(0),
      fMuonTrackCuts(0x0),
      fRunNum(0),
      fTracklets(0),
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
      fInvariantMassDistributionCoherentZNCzeroZNAzeroH(0),
      fInvariantMassDistributionCoherentZNCzeroZNAzeroShiftMinusTwoH(0),
      fInvariantMassDistributionCoherentZNCzeroZNAzeroShiftMinusOneH(0),
      fInvariantMassDistributionCoherentZNCzeroZNAzeroShiftPlusOneH(0),
      fInvariantMassDistributionCoherentZNCzeroZNAzeroShiftPlusTwoH(0),
      fInvariantMassDistributionCoherentZNCzeroZNAzeroHv2(0),
      fInvariantMassDistributionCoherentZNCzeroZNAzeroRapidityHv2{0, 0, 0},
      fInvariantMassDistributionCoherentZNCzeroZNAanyH(0),
      fInvariantMassDistributionCoherentZNCzeroZNAanyShiftMinusTwoH(0),
      fInvariantMassDistributionCoherentZNCzeroZNAanyShiftMinusOneH(0),
      fInvariantMassDistributionCoherentZNCzeroZNAanyShiftPlusOneH(0),
      fInvariantMassDistributionCoherentZNCzeroZNAanyShiftPlusTwoH(0),
      fInvariantMassDistributionCoherentZNCzeroZNAanyHv2(0),
      fInvariantMassDistributionCoherentZNCzeroZNAanyRapidityHv2{0, 0, 0},
      fInvariantMassDistributionCoherentZNCanyZNAzeroH(0),
      fInvariantMassDistributionCoherentZNCanyZNAzeroShiftMinusTwoH(0),
      fInvariantMassDistributionCoherentZNCanyZNAzeroShiftMinusOneH(0),
      fInvariantMassDistributionCoherentZNCanyZNAzeroShiftPlusOneH(0),
      fInvariantMassDistributionCoherentZNCanyZNAzeroShiftPlusTwoH(0),
      fInvariantMassDistributionCoherentZNCanyZNAzeroHv2(0),
      fInvariantMassDistributionCoherentZNCanyZNAzeroRapidityHv2{0, 0, 0},
      fInvariantMassDistributionCoherentZNCanyZNAanyH(0),
      fInvariantMassDistributionCoherentZNCanyZNAanyShiftMinusTwoH(0),
      fInvariantMassDistributionCoherentZNCanyZNAanyShiftMinusOneH(0),
      fInvariantMassDistributionCoherentZNCanyZNAanyShiftPlusOneH(0),
      fInvariantMassDistributionCoherentZNCanyZNAanyShiftPlusTwoH(0),
      fInvariantMassDistributionCoherentZNCanyZNAanyHv2(0),
      fInvariantMassDistributionCoherentZNCanyZNAanyRapidityHv2{0, 0, 0},
      fInvariantMassDistributionIncoherentZNCzeroZNAzeroH(0),
      fInvariantMassDistributionIncoherentZNCzeroZNAzeroShiftMinusTwoH(0),
      fInvariantMassDistributionIncoherentZNCzeroZNAzeroShiftMinusOneH(0),
      fInvariantMassDistributionIncoherentZNCzeroZNAzeroShiftPlusOneH(0),
      fInvariantMassDistributionIncoherentZNCzeroZNAzeroShiftPlusTwoH(0),
      fInvariantMassDistributionIncoherentZNCzeroZNAzeroHv2(0),
      fInvariantMassDistributionIncoherentZNCzeroZNAzeroRapidityHv2{0, 0, 0},
      fInvariantMassDistributionIncoherentZNCzeroZNAanyH(0),
      fInvariantMassDistributionIncoherentZNCzeroZNAanyShiftMinusTwoH(0),
      fInvariantMassDistributionIncoherentZNCzeroZNAanyShiftMinusOneH(0),
      fInvariantMassDistributionIncoherentZNCzeroZNAanyShiftPlusOneH(0),
      fInvariantMassDistributionIncoherentZNCzeroZNAanyShiftPlusTwoH(0),
      fInvariantMassDistributionIncoherentZNCzeroZNAanyHv2(0),
      fInvariantMassDistributionIncoherentZNCzeroZNAanyRapidityHv2{0, 0, 0},
      fInvariantMassDistributionIncoherentZNCanyZNAzeroH(0),
      fInvariantMassDistributionIncoherentZNCanyZNAzeroShiftMinusTwoH(0),
      fInvariantMassDistributionIncoherentZNCanyZNAzeroShiftMinusOneH(0),
      fInvariantMassDistributionIncoherentZNCanyZNAzeroShiftPlusOneH(0),
      fInvariantMassDistributionIncoherentZNCanyZNAzeroShiftPlusTwoH(0),
      fInvariantMassDistributionIncoherentZNCanyZNAzeroHv2(0),
      fInvariantMassDistributionIncoherentZNCanyZNAzeroRapidityHv2{0, 0, 0},
      fInvariantMassDistributionIncoherentZNCanyZNAanyH(0),
      fInvariantMassDistributionIncoherentZNCanyZNAanyShiftMinusTwoH(0),
      fInvariantMassDistributionIncoherentZNCanyZNAanyShiftMinusOneH(0),
      fInvariantMassDistributionIncoherentZNCanyZNAanyShiftPlusOneH(0),
      fInvariantMassDistributionIncoherentZNCanyZNAanyShiftPlusTwoH(0),
      fInvariantMassDistributionIncoherentZNCanyZNAanyHv2(0),
      fInvariantMassDistributionIncoherentZNCanyZNAanyRapidityHv2{0, 0, 0},
      fDimuonPtDistributionZNCzeroZNAzeroH(0),
      fDimuonPtDistributionZNCzeroZNAanyH(0),
      fDimuonPtDistributionZNCanyZNAzeroH(0),
      fDimuonPtDistributionZNCanyZNAanyH(0),
      fDimuonPtDistributionZNCzeroZNAzeroShiftPlusOneH(0),
      fDimuonPtDistributionZNCzeroZNAanyShiftPlusOneH(0),
      fDimuonPtDistributionZNCanyZNAzeroShiftPlusOneH(0),
      fDimuonPtDistributionZNCanyZNAanyShiftPlusOneH(0),
      fDimuonPtDistributionZNCzeroZNAzeroHv2(0),
      fDimuonPtDistributionZNCzeroZNAanyHv2(0),
      fDimuonPtDistributionZNCanyZNAzeroHv2(0),
      fDimuonPtDistributionZNCanyZNAanyHv2(0),
      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv2{0, 0, 0},
      fDimuonPtDistributionZNCzeroZNAanyRapidityHv2{0, 0, 0},
      fDimuonPtDistributionZNCanyZNAzeroRapidityHv2{0, 0, 0},
      fDimuonPtDistributionZNCanyZNAanyRapidityHv2{0, 0, 0},
      fDimuonPtDistributionCoherentZNCzeroZNAzeroH(0),
      fDimuonPtDistributionCoherentZNCzeroZNAanyH(0),
      fDimuonPtDistributionCoherentZNCanyZNAzeroH(0),
      fDimuonPtDistributionCoherentZNCanyZNAanyH(0),
      fDimuonPtDistributionIncoherentZNCzeroZNAzeroH(0),
      fDimuonPtDistributionIncoherentZNCzeroZNAanyH(0),
      fDimuonPtDistributionIncoherentZNCanyZNAzeroH(0),
      fDimuonPtDistributionIncoherentZNCanyZNAanyH(0),
      fAngularDistribOfPositiveMuonRestFrameJPsiH(0),
      fAngularDistribOfNegativeMuonRestFrameJPsiH(0),
      fCheckHelicityRestFrameJPsiH(0),
      fThetaDistribOfPositiveMuonRestFrameJPsiRapidityBinH{ 0, 0, 0, 0, 0, 0, 0, 0},
      fCosThetaHelicityFrameJPsiH(0),
      fCosThetaCollinsSoperFrameJPsiH(0),
      fPhiHelicityFrameJPsiH(0),
      fPhiCollinsSoperFrameJPsiH(0),
      fCosThetaHelicityFrameJPsiRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0},
      fCosThetaCollinsSoperFrameJPsiRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0},
      fPhiHelicityFrameJPsiRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0},
      fPhiCollinsSoperFrameJPsiRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0},
      fCosThetaHelicityFrameJPsiTenRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      fCosThetaCollinsSoperFrameJPsiTenRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      fPhiHelicityFrameJPsiTenRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      fPhiCollinsSoperFrameJPsiTenRapidityBinsH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      fInvariantMassDistributionInBinsOfCosThetaHelicityFrameH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      fInvariantMassDistributionBinsOfCosThetaAndPhiHelicityFrameH(0),
      fCosThetaAndPhiHelicityFrameInclusivePeopleBinningH(0),
      fInvariantMassDistributionForSignalExtractionHelicityFrameH(0),
      fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameTryingH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameTryingH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                                0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                                0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                                0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                                0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fInvariantMassDistributionForSignalExtractionHelicityFrameMyBinningH(0),
      fInvariantMassDistributionForSignalExtractionCsFrameMyBinningH(0),
      fInvariantMassDistributionForSignalExtractionHelicityFrameMyBinning0N0NH(0),
      fInvariantMassDistributionForSignalExtractionCsFrameMyBinning0N0NH(0),
      fInvariantMassDistributionStrictPtH(0),
      fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameMyVariableBinningH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                                                0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                                                0, 0, 0, 0, 0, 0 },
      fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameMyVariableBinningH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameMyStrictVariableBinningH{ 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameMyStrictVariableBinningH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameMySeventeenBinsVariableBinningH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fCosThetaHelicityFrameJPsiAlreadyCorrectedH(0),
      fDcaAgainstPtOfVectorMesonH(0),
      fInvariantMassDistributionStrictPtStrictDcaH(0),
      fInvariantMassDistributionVsPtStrictDcaH(0),
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
      fInvariantMassDistributionOnlyPhiHeFrameTwentyfiveBins0N0NH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                   0, 0, 0, 0, 0 },
      fInvariantMassDistributionOnlyCosThetaHeFrameTwentyfiveBins0N0NH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                        0, 0, 0, 0, 0 },
      fInvariantMassDistributionOnlyTildePhiHeFrameTwentyfiveBins0N0NH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                        0, 0, 0, 0, 0 },
      fInvariantMassDistributionOnlyPhiCsFrameTwentyfiveBins0N0NH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                   0, 0, 0, 0, 0 },
      fInvariantMassDistributionOnlyCosThetaCsFrameTwentyfiveBins0N0NH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                        0, 0, 0, 0, 0 },
      fInvariantMassDistributionOnlyTildePhiCsFrameTwentyfiveBins0N0NH{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                        0, 0, 0, 0, 0 }
{
    // FillGoodRunVector(fVectorGoodRunNumbers);

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
AliAnalysisTaskUPCforward::~AliAnalysisTaskUPCforward()
{
    // destructor
    if(fOutputList)    {delete fOutputList;}     	// at the end of your task, it is deleted
    if(fMuonTrackCuts) {delete fMuonTrackCuts;}   // from memory by calling this function
}
//_____________________________________________________________________________
// void AliAnalysisTaskUPCforward::FillGoodRunVector(std::vector<Int_t> &fVectorGoodRunNumbers)
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
//                                          297590, 297595/*, 297623, 297624*/ };
//   /* - This good run number list has been taken from the analysis
//      - note of Kay's talk for DIS 2017, see:
//      - https://alice-notes.web.cern.ch/system/files/notes/analysis/596/2017-Feb-08-analysis_note-2017-Feb-08-analysis-note.pdf
//      -
//    */
//   Int_t listOfGoodRunNumbersLHC15o[] = { 244918, 244980, 244982, 244983, 245064, 245066, 245068, 245145, 245146, 245151,
//                                          245152, 245231, 245232, 245233, 245253, 245259, 245343, 245345, 245346, 245347,
//                                          245353, 245401, 245407, 245409, 245410, 245446, 245450, 245496, 245501, 245504,
//                                          245505, 245507, 245535, 245540, 245542, 245543, 245554, 245683, 245692, 245700,
//                                          245705, 245729, 245731, 245738, 245752, 245759, 245766, 245775, 245785, 245793,
//                                          245829, 245831, 245833, 245949, 245952, 245954, 245963, 245996, 246001, 246003,
//                                          246012, 246036, 246037, 246042, 246048, 246049, 246053, 246087, 246089, 246113,
//                                          246115, 246148, 246151, 246152, 246153, 246178, 246181, 246182, 246217, 246220,
//                                          246222, 246225, 246272, 246275, 246276, 246390, 246391, 246392, 246424, 246428,
//                                          246431, 246433, 246434, 246487, 246488, 246493, 246495, 246675, 246676, 246750,
//                                          246751, 246755, 246757, 246758, 246759, 246760, 246763, 246765, 246804, 246805,
//                                          246806, 246807, 246808, 246809, 246844, 246845, 246846, 246847, 246851, 246855,
//                                          246859, 246864, 246865, 246867, 246871, 246930, 246937, 246942, 246945, 246948,
//                                          246949, 246980, 246982, 246984, 246989, 246991, 246994
//                                        };
//   Int_t sizeOfLHC18q = 0;
//   Int_t sizeOfLHC18r = 0;
//   Int_t sizeOfLHC15o = 0;
//   for ( Int_t GoodRunNumberLHC18q : listOfGoodRunNumbersLHC18q ) {
//         fVectorGoodRunNumbers.push_back(GoodRunNumberLHC18q);
//         sizeOfLHC18q++;
//   }
//   for ( Int_t GoodRunNumberLHC18r : listOfGoodRunNumbersLHC18r ) {
//         fVectorGoodRunNumbers.push_back(GoodRunNumberLHC18r);
//         sizeOfLHC18r++;
//   }
//   for ( Int_t GoodRunNumberLHC15o : listOfGoodRunNumbersLHC15o ) {
//         fVectorGoodRunNumbers.push_back(GoodRunNumberLHC15o);
//         sizeOfLHC15o++;
//   }
//   cout << std::endl << "LHC18q GOOD RUNS:  " << std::endl;
//   for ( Int_t i = 0; i < sizeOfLHC18q; i++ ) {
//         cout << fVectorGoodRunNumbers.at(i) << ",   number: " << i << std::endl;
//   }
//   cout << std::endl << "LHC18r GOOD RUNS:  " << std::endl;
//   for ( Int_t i = sizeOfLHC18q; i < sizeOfLHC18q + sizeOfLHC18r; i++ ) {
//         cout << fVectorGoodRunNumbers.at(i) << ",   number: " << (i-sizeOfLHC18q) << std::endl;
//   }
//   cout << std::endl << "LHC15o GOOD RUNS:  " << std::endl;
//   for ( Int_t i = sizeOfLHC18q + sizeOfLHC18r; i < sizeOfLHC18q + sizeOfLHC18r + sizeOfLHC15o; i++ ) {
//         cout << fVectorGoodRunNumbers.at(i) << ",   number: " << (i-sizeOfLHC18q-sizeOfLHC18r) << std::endl;
//   }
// }
//_____________________________________________________________________________
void AliAnalysisTaskUPCforward::UserCreateOutputObjects()
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

  fPtSingleMuonH = new TH1F("fPtSingleMuonH", "fPtSingleMuonH", 1000, 0., 10.);
  fOutputList->Add(fPtSingleMuonH);

  fRAbsMuonH = new TH1F("fRAbsMuonH", "fRAbsMuonH", 100, 0, 100);
  fOutputList->Add(fRAbsMuonH);

  fInvariantMassDistributionH = new TH1F("fInvariantMassDistributionH", "fInvariantMassDistributionH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionH);

  fInvariantMassDistributionAtDcaH = new TH1F("fInvariantMassDistributionAtDcaH", "fInvariantMassDistributionAtDcaH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionAtDcaH);

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
  // fEntriesAgainstRunNumberProperlyH->SetCanExtend(TH1::kAllAxes);
  fEntriesAgainstRunNumberProperlyH->LabelsDeflate();
  fOutputList->Add(fEntriesAgainstRunNumberProperlyH);

  fRunNumberTriggerCMUP11ClassH = new TH1F("fRunNumberTriggerCMUP11ClassH", "fRunNumberTriggerCMUP11ClassH", 10000, 290000, 300000);
  fOutputList->Add(fRunNumberTriggerCMUP11ClassH);

  fRunNumberTriggerCMUP11ClassProperlyH = new TH1F("fRunNumberTriggerCMUP11ClassProperlyH", "fRunNumberTriggerCMUP11ClassProperlyH", 10000, 290000, 300000);
  fRunNumberTriggerCMUP11ClassProperlyH->SetStats(0);
  fRunNumberTriggerCMUP11ClassProperlyH->SetFillColor(38);
  // fRunNumberTriggerCMUP11ClassProperlyH->SetCanExtend(TH1::kAllAxes);
  fRunNumberTriggerCMUP11ClassProperlyH->LabelsDeflate();
  fOutputList->Add(fRunNumberTriggerCMUP11ClassProperlyH);

  fRunNumberTriggerCMUP26ClassH = new TH1F("fRunNumberTriggerCMUP26ClassH", "fRunNumberTriggerCMUP26ClassH", 10000, 290000, 300000);
  fOutputList->Add(fRunNumberTriggerCMUP26ClassH);

  fRunNumberTriggerCMUP26ClassProperlyH = new TH1F("fRunNumberTriggerCMUP26ClassProperlyH", "fRunNumberTriggerCMUP26ClassProperlyH", 10000, 290000, 300000);
  fRunNumberTriggerCMUP26ClassProperlyH->SetStats(0);
  fRunNumberTriggerCMUP26ClassProperlyH->SetFillColor(38);
  // fRunNumberTriggerCMUP26ClassProperlyH->SetCanExtend(TH1::kAllAxes);
  fRunNumberTriggerCMUP26ClassProperlyH->LabelsDeflate();
  fOutputList->Add(fRunNumberTriggerCMUP26ClassProperlyH);

  fRunNumberTriggerCMUP6ClassH = new TH1F("fRunNumberTriggerCMUP6ClassH", "fRunNumberTriggerCMUP6ClassH", 10000, 290000, 300000);
  fOutputList->Add(fRunNumberTriggerCMUP6ClassH);

  fRunNumberTriggerCMUP6ClassProperlyH = new TH1F("fRunNumberTriggerCMUP6ClassProperlyH", "fRunNumberTriggerCMUP6ClassProperlyH", 10000, 290000, 300000);
  fRunNumberTriggerCMUP6ClassProperlyH->SetStats(0);
  fRunNumberTriggerCMUP6ClassProperlyH->SetFillColor(38);
  // fRunNumberTriggerCMUP6ClassProperlyH->SetCanExtend(TH1::kAllAxes);
  fRunNumberTriggerCMUP6ClassProperlyH->LabelsDeflate();
  fOutputList->Add(fRunNumberTriggerCMUP6ClassProperlyH);

  fRunNumberTriggerCMUP10ClassH = new TH1F("fRunNumberTriggerCMUP10ClassH", "fRunNumberTriggerCMUP10ClassH", 10000, 290000, 300000);
  fOutputList->Add(fRunNumberTriggerCMUP10ClassH);

  fRunNumberTriggerCMUP10ClassProperlyH = new TH1F("fRunNumberTriggerCMUP10ClassProperlyH", "fRunNumberTriggerCMUP10ClassProperlyH", 10000, 290000, 300000);
  fRunNumberTriggerCMUP10ClassProperlyH->SetStats(0);
  fRunNumberTriggerCMUP10ClassProperlyH->SetFillColor(38);
  // fRunNumberTriggerCMUP10ClassProperlyH->SetCanExtend(TH1::kAllAxes);
  fRunNumberTriggerCMUP10ClassProperlyH->LabelsDeflate();
  fOutputList->Add(fRunNumberTriggerCMUP10ClassProperlyH);

  fRunNumberTriggerCMUP13ClassH = new TH1F("fRunNumberTriggerCMUP13ClassH", "fRunNumberTriggerCMUP13ClassH", 10000, 290000, 300000);
  fOutputList->Add(fRunNumberTriggerCMUP13ClassH);

  fRunNumberTriggerCMUP13ClassProperlyH = new TH1F("fRunNumberTriggerCMUP13ClassProperlyH", "fRunNumberTriggerCMUP13ClassProperlyH", 10000, 290000, 300000);
  fRunNumberTriggerCMUP13ClassProperlyH->SetStats(0);
  fRunNumberTriggerCMUP13ClassProperlyH->SetFillColor(38);
  // fRunNumberTriggerCMUP13ClassProperlyH->SetCanExtend(TH1::kAllAxes);
  fRunNumberTriggerCMUP13ClassProperlyH->LabelsDeflate();
  fOutputList->Add(fRunNumberTriggerCMUP13ClassProperlyH);

  fTriggersVsRunH = new TH2F("fTriggersVsRunH","",5,0,5,60000,240000,300000);
  fOutputList->Add(fTriggersVsRunH);

  fInvariantMassDistributionCoherentH = new TH1F("fInvariantMassDistributionCoherentH", "fInvariantMassDistributionCoherentH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionCoherentH);

  for( Int_t iRapidity = 0; iRapidity < 6; iRapidity++ ) {
    fInvariantMassDistributionCoherentRapidityBinsH[iRapidity]
            = new TH1F( Form("fInvariantMassDistributionCoherentRapidityBinsH_%d", iRapidity),
                        Form("fInvariantMassDistributionCoherentRapidityBinsH_%d", iRapidity),
                        2000, 0, 20);
    fOutputList->Add(fInvariantMassDistributionCoherentRapidityBinsH[iRapidity]);
  }


  fInvariantMassDistributionCoherentShiftMinusTwoH = new TH1F("fInvariantMassDistributionCoherentShiftMinusTwoH", "fInvariantMassDistributionCoherentShiftMinusTwoH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionCoherentShiftMinusTwoH);

  fInvariantMassDistributionCoherentShiftMinusOneH = new TH1F("fInvariantMassDistributionCoherentShiftMinusOneH", "fInvariantMassDistributionCoherentShiftMinusOneH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionCoherentShiftMinusOneH);

  fInvariantMassDistributionCoherentShiftPlusOneH = new TH1F("fInvariantMassDistributionCoherentShiftPlusOneH", "fInvariantMassDistributionCoherentShiftPlusOneH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionCoherentShiftPlusOneH);

  fInvariantMassDistributionCoherentShiftPlusTwoH = new TH1F("fInvariantMassDistributionCoherentShiftPlusTwoH", "fInvariantMassDistributionCoherentShiftPlusTwoH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionCoherentShiftPlusTwoH);

  fInvariantMassDistributionIncoherentH = new TH1F("fInvariantMassDistributionIncoherentH", "fInvariantMassDistributionIncoherentH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionIncoherentH);

  fInvariantMassDistributionIncoherentShiftMinusTwoH = new TH1F("fInvariantMassDistributionIncoherentShiftMinusTwoH", "fInvariantMassDistributionIncoherentShiftMinusTwoH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionIncoherentShiftMinusTwoH);

  fInvariantMassDistributionIncoherentShiftMinusOneH = new TH1F("fInvariantMassDistributionIncoherentShiftMinusOneH", "fInvariantMassDistributionIncoherentShiftMinusOneH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionIncoherentShiftMinusOneH);

  fInvariantMassDistributionIncoherentShiftPlusOneH = new TH1F("fInvariantMassDistributionIncoherentShiftPlusOneH", "fInvariantMassDistributionIncoherentShiftPlusOneH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionIncoherentShiftPlusOneH);

  fInvariantMassDistributionIncoherentShiftPlusTwoH = new TH1F("fInvariantMassDistributionIncoherentShiftPlusTwoH", "fInvariantMassDistributionIncoherentShiftPlusTwoH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionIncoherentShiftPlusTwoH);

  fDimuonPtDistributionH = new TH1F("fDimuonPtDistributionH", "fDimuonPtDistributionH", 4000, 0, 20);
  fOutputList->Add(fDimuonPtDistributionH);

  fDimuonPtDistributionShiftPlusOneH = new TH1F("fDimuonPtDistributionShiftPlusOneH", "fDimuonPtDistributionShiftPlusOneH", 4000, 0.02, 20.02);
  fOutputList->Add(fDimuonPtDistributionShiftPlusOneH);

  fZNCEnergyAgainstEntriesH = new TH1F("fZNCEnergyAgainstEntriesH", "fZNCEnergyAgainstEntriesH", 20000, -10000, 40000);
  fOutputList->Add(fZNCEnergyAgainstEntriesH);

  fZNAEnergyAgainstEntriesH = new TH1F("fZNAEnergyAgainstEntriesH", "fZNAEnergyAgainstEntriesH", 20000, -10000, 40000);
  fOutputList->Add(fZNAEnergyAgainstEntriesH);

  fZNCEnergyBeforeTimingSelectionH = new TH1F("fZNCEnergyBeforeTimingSelectionH", "fZNCEnergyBeforeTimingSelectionH", 20000, -10000, 40000);
  fOutputList->Add(fZNCEnergyBeforeTimingSelectionH);

  fZNAEnergyBeforeTimingSelectionH = new TH1F("fZNAEnergyBeforeTimingSelectionH", "fZNAEnergyBeforeTimingSelectionH", 20000, -10000, 40000);
  fOutputList->Add(fZNAEnergyBeforeTimingSelectionH);

  fZNCEnergyAgainstEntriesExtendedH = new TH1F("fZNCEnergyAgainstEntriesExtendedH", "fZNCEnergyAgainstEntriesExtendedH", 20000, -10000, 400000);
  fOutputList->Add(fZNCEnergyAgainstEntriesExtendedH);

  fZNAEnergyAgainstEntriesExtendedH = new TH1F("fZNAEnergyAgainstEntriesExtendedH", "fZNAEnergyAgainstEntriesExtendedH", 20000, -10000, 400000);
  fOutputList->Add(fZNAEnergyAgainstEntriesExtendedH);

  fZNCEnergyAgainstEntriesExtendedHv2 = new TH1F("fZNCEnergyAgainstEntriesExtendedHv2", "fZNCEnergyAgainstEntriesExtendedHv2", 20000, -10000, 400000);
  fOutputList->Add(fZNCEnergyAgainstEntriesExtendedHv2);

  fZNAEnergyAgainstEntriesExtendedHv2 = new TH1F("fZNAEnergyAgainstEntriesExtendedHv2", "fZNAEnergyAgainstEntriesExtendedHv2", 20000, -10000, 400000);
  fOutputList->Add(fZNAEnergyAgainstEntriesExtendedHv2);

  fZNCEnergyBeforeTimingSelectionExtendedH = new TH1F("fZNCEnergyBeforeTimingSelectionExtendedH", "fZNCEnergyBeforeTimingSelectionExtendedH", 20000, -10000, 400000);
  fOutputList->Add(fZNCEnergyBeforeTimingSelectionExtendedH);

  fZNAEnergyBeforeTimingSelectionExtendedH = new TH1F("fZNAEnergyBeforeTimingSelectionExtendedH", "fZNAEnergyBeforeTimingSelectionExtendedH", 20000, -10000, 400000);
  fOutputList->Add(fZNAEnergyBeforeTimingSelectionExtendedH);

  fZNCEnergyCalibratedH = new TH1F("fZNCEnergyCalibratedH", "fZNCEnergyCalibratedH", 20000, -10000, 40000);
  fOutputList->Add(fZNCEnergyCalibratedH);

  fZNAEnergyCalibratedH = new TH1F("fZNAEnergyCalibratedH", "fZNAEnergyCalibratedH", 20000, -10000, 40000);
  fOutputList->Add(fZNAEnergyCalibratedH);

  fZNCEnergyUncalibratedH = new TH1F("fZNCEnergyUncalibratedH", "fZNCEnergyUncalibratedH", 20000, -10000, 40000);
  fOutputList->Add(fZNCEnergyUncalibratedH);

  fZNAEnergyUncalibratedH = new TH1F("fZNAEnergyUncalibratedH", "fZNAEnergyUncalibratedH", 20000, -10000, 40000);
  fOutputList->Add(fZNAEnergyUncalibratedH);

  fZNCEnergyCalibratedHigherGainH = new TH1F("fZNCEnergyCalibratedHigherGainH", "fZNCEnergyCalibratedHigherGainH", 20000, -80000, 320000);
  fOutputList->Add(fZNCEnergyCalibratedHigherGainH);

  fZNAEnergyCalibratedHigherGainH = new TH1F("fZNAEnergyCalibratedHigherGainH", "fZNAEnergyCalibratedHigherGainH", 20000, -80000, 320000);
  fOutputList->Add(fZNAEnergyCalibratedHigherGainH);

  fZNCTimeAgainstEntriesH = new TH1F("fZNCTimeAgainstEntriesH", "fZNCTimeAgainstEntriesH", 6000, -1500, 1500);
  fOutputList->Add(fZNCTimeAgainstEntriesH);

  fZNATimeAgainstEntriesH = new TH1F("fZNATimeAgainstEntriesH", "fZNATimeAgainstEntriesH", 6000, -1500, 1500);
  fOutputList->Add(fZNATimeAgainstEntriesH);

  fZNCTimeStrictTimeWindowH = new TH1F("fZNCTimeStrictTimeWindowH", "fZNCTimeStrictTimeWindowH", 6000, -1500, 1500);
  fOutputList->Add(fZNCTimeStrictTimeWindowH);

  fZNATimeStrictTimeWindowH = new TH1F("fZNATimeStrictTimeWindowH", "fZNATimeStrictTimeWindowH", 6000, -1500, 1500);
  fOutputList->Add(fZNATimeStrictTimeWindowH);

  for(int iTiming = 0; iTiming < 4; iTiming++) {
    fZNCTimeWithoutTimingH[iTiming] = new TH1F( Form("fZNCTimeWithoutTimingH_%d", iTiming),
                                                Form("fZNCTimeWithoutTimingH_%d", iTiming),
                                                6000, -1500, 1500
                                               );
    fOutputList->Add(fZNCTimeWithoutTimingH[iTiming]);
  }

  for(int iTiming = 0; iTiming < 4; iTiming++) {
    fZNATimeWithoutTimingH[iTiming] = new TH1F( Form("fZNATimeWithoutTimingH_%d", iTiming),
                                                Form("fZNATimeWithoutTimingH_%d", iTiming),
                                                6000, -1500, 1500
                                               );
    fOutputList->Add(fZNATimeWithoutTimingH[iTiming]);
  }

  fZNCTime4FillingH = new TH1F("fZNCTime4FillingH", "fZNCTime4FillingH", 6000, -1500, 1500);
  fOutputList->Add(fZNCTime4FillingH);

  fZNATime4FillingH = new TH1F("fZNATime4FillingH", "fZNATime4FillingH", 6000, -1500, 1500);
  fOutputList->Add(fZNATime4FillingH);

  for(int iTiming = 0; iTiming < 4; iTiming++) {
    fZNCminusZNAtimeVsZNCplusZNAtimeH[iTiming] = new TH2F( Form("fZNCminusZNAtimeVsZNCplusZNAtimeH_%d", iTiming),
                                                           Form("fZNCminusZNAtimeVsZNCplusZNAtimeH_%d", iTiming),
                                                           1200, -300, 300, 1200, -300, 300
                                                           );
    fOutputList->Add(fZNCminusZNAtimeVsZNCplusZNAtimeH[iTiming]);
  }

  fZNCminusZNAtimeVsZNCplusZNAtime4FillingH = new TH2F("fZNCminusZNAtimeVsZNCplusZNAtime4FillingH", "fZNCminusZNAtimeVsZNCplusZNAtime4FillingH", 1200, -300, 300, 1200, -300, 300);
  fOutputList->Add(fZNCminusZNAtimeVsZNCplusZNAtime4FillingH);

  fCounterZNCH = new TH1F("fCounterZNCH", "fCounterZNCH", 6, -0.5, 5.5);
  fOutputList->Add(fCounterZNCH);

  fCounterZNAH = new TH1F("fCounterZNAH", "fCounterZNAH", 6, -0.5, 5.5);
  fOutputList->Add(fCounterZNAH);

  fInvariantMassDistributionNoNeutronsH = new TH1F("fInvariantMassDistributionNoNeutronsH", "fInvariantMassDistributionNoNeutronsH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionNoNeutronsH);

  fInvariantMassDistributionOneNeutronH = new TH1F("fInvariantMassDistributionOneNeutronH", "fInvariantMassDistributionOneNeutronH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionOneNeutronH);

  fInvariantMassDistributionAtLeastOneNeutronH = new TH1F("fInvariantMassDistributionAtLeastOneNeutronH", "fInvariantMassDistributionAtLeastOneNeutronH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionAtLeastOneNeutronH);

  fInvariantMassDistributionCoherentNoNeutronsH = new TH1F("fInvariantMassDistributionCoherentNoNeutronsH", "fInvariantMassDistributionCoherentNoNeutronsH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionCoherentNoNeutronsH);

  fInvariantMassDistributionCoherentOneNeutronH = new TH1F("fInvariantMassDistributionCoherentOneNeutronH", "fInvariantMassDistributionCoherentOneNeutronH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionCoherentOneNeutronH);

  fInvariantMassDistributionCoherentAtLeastOneNeutronH = new TH1F("fInvariantMassDistributionCoherentAtLeastOneNeutronH", "fInvariantMassDistributionCoherentAtLeastOneNeutronH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionCoherentAtLeastOneNeutronH);

  fInvariantMassDistributionIncoherentNoNeutronsH = new TH1F("fInvariantMassDistributionIncoherentNoNeutronsH", "fInvariantMassDistributionIncoherentNoNeutronsH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionIncoherentNoNeutronsH);

  fInvariantMassDistributionIncoherentOneNeutronH = new TH1F("fInvariantMassDistributionIncoherentOneNeutronH", "fInvariantMassDistributionIncoherentOneNeutronH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionIncoherentOneNeutronH);

  fInvariantMassDistributionIncoherentAtLeastOneNeutronH = new TH1F("fInvariantMassDistributionIncoherentAtLeastOneNeutronH", "fInvariantMassDistributionIncoherentAtLeastOneNeutronH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionIncoherentAtLeastOneNeutronH);

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


  /* - These histograms pertain the differential neutron emission analysis.
     -
   */
  fInvariantMassDistributionCoherentZNCzeroZNAzeroH = new TH1F("fInvariantMassDistributionCoherentZNCzeroZNAzeroH", "fInvariantMassDistributionCoherentZNCzeroZNAzeroH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionCoherentZNCzeroZNAzeroH);

  fInvariantMassDistributionCoherentZNCzeroZNAzeroShiftMinusTwoH = new TH1F("fInvariantMassDistributionCoherentZNCzeroZNAzeroShiftMinusTwoH", "fInvariantMassDistributionCoherentZNCzeroZNAzeroShiftMinusTwoH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionCoherentZNCzeroZNAzeroShiftMinusTwoH);

  fInvariantMassDistributionCoherentZNCzeroZNAzeroShiftMinusOneH = new TH1F("fInvariantMassDistributionCoherentZNCzeroZNAzeroShiftMinusOneH", "fInvariantMassDistributionCoherentZNCzeroZNAzeroShiftMinusOneH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionCoherentZNCzeroZNAzeroShiftMinusOneH);

  fInvariantMassDistributionCoherentZNCzeroZNAzeroShiftPlusOneH = new TH1F("fInvariantMassDistributionCoherentZNCzeroZNAzeroShiftPlusOneH", "fInvariantMassDistributionCoherentZNCzeroZNAzeroShiftPlusOneH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionCoherentZNCzeroZNAzeroShiftPlusOneH);

  fInvariantMassDistributionCoherentZNCzeroZNAzeroShiftPlusTwoH = new TH1F("fInvariantMassDistributionCoherentZNCzeroZNAzeroShiftPlusTwoH", "fInvariantMassDistributionCoherentZNCzeroZNAzeroShiftPlusTwoH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionCoherentZNCzeroZNAzeroShiftPlusTwoH);

  fInvariantMassDistributionCoherentZNCzeroZNAzeroHv2 = new TH1F("fInvariantMassDistributionCoherentZNCzeroZNAzeroHv2", "fInvariantMassDistributionCoherentZNCzeroZNAzeroHv2", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionCoherentZNCzeroZNAzeroHv2);

  for( Int_t iRapidity = 0; iRapidity < 3; iRapidity++ ){
    fInvariantMassDistributionCoherentZNCzeroZNAzeroRapidityHv2[iRapidity] = new TH1F(
              Form("fInvariantMassDistributionCoherentZNCzeroZNAzeroRapidityHv2_%d", iRapidity),
              Form("fInvariantMassDistributionCoherentZNCzeroZNAzeroRapidityHv2_%d", iRapidity),
              2000, 0, 20
              );
    fOutputList->Add(fInvariantMassDistributionCoherentZNCzeroZNAzeroRapidityHv2[iRapidity]);
  }

  fInvariantMassDistributionCoherentZNCzeroZNAanyH = new TH1F("fInvariantMassDistributionCoherentZNCzeroZNAanyH", "fInvariantMassDistributionCoherentZNCzeroZNAanyH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionCoherentZNCzeroZNAanyH);

  fInvariantMassDistributionCoherentZNCzeroZNAanyShiftMinusTwoH = new TH1F("fInvariantMassDistributionCoherentZNCzeroZNAanyShiftMinusTwoH", "fInvariantMassDistributionCoherentZNCzeroZNAanyShiftMinusTwoH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionCoherentZNCzeroZNAanyShiftMinusTwoH);

  fInvariantMassDistributionCoherentZNCzeroZNAanyShiftMinusOneH = new TH1F("fInvariantMassDistributionCoherentZNCzeroZNAanyShiftMinusOneH", "fInvariantMassDistributionCoherentZNCzeroZNAanyShiftMinusOneH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionCoherentZNCzeroZNAanyShiftMinusOneH);

  fInvariantMassDistributionCoherentZNCzeroZNAanyShiftPlusOneH = new TH1F("fInvariantMassDistributionCoherentZNCzeroZNAanyShiftPlusOneH", "fInvariantMassDistributionCoherentZNCzeroZNAanyShiftPlusOneH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionCoherentZNCzeroZNAanyShiftPlusOneH);

  fInvariantMassDistributionCoherentZNCzeroZNAanyShiftPlusTwoH = new TH1F("fInvariantMassDistributionCoherentZNCzeroZNAanyShiftPlusTwoH", "fInvariantMassDistributionCoherentZNCzeroZNAanyShiftPlusTwoH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionCoherentZNCzeroZNAanyShiftPlusTwoH);

  fInvariantMassDistributionCoherentZNCzeroZNAanyHv2 = new TH1F("fInvariantMassDistributionCoherentZNCzeroZNAanyHv2", "fInvariantMassDistributionCoherentZNCzeroZNAanyHv2", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionCoherentZNCzeroZNAanyHv2);

  for( Int_t iRapidity = 0; iRapidity < 3; iRapidity++ ){
    fInvariantMassDistributionCoherentZNCzeroZNAanyRapidityHv2[iRapidity] = new TH1F(
              Form("fInvariantMassDistributionCoherentZNCzeroZNAanyRapidityHv2_%d", iRapidity),
              Form("fInvariantMassDistributionCoherentZNCzeroZNAanyRapidityHv2_%d", iRapidity),
              2000, 0, 20
              );
    fOutputList->Add(fInvariantMassDistributionCoherentZNCzeroZNAanyRapidityHv2[iRapidity]);
  }

  fInvariantMassDistributionCoherentZNCanyZNAzeroH = new TH1F("fInvariantMassDistributionCoherentZNCanyZNAzeroH", "fInvariantMassDistributionCoherentZNCanyZNAzeroH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionCoherentZNCanyZNAzeroH);

  fInvariantMassDistributionCoherentZNCanyZNAzeroShiftMinusTwoH = new TH1F("fInvariantMassDistributionCoherentZNCanyZNAzeroShiftMinusTwoH", "fInvariantMassDistributionCoherentZNCanyZNAzeroShiftMinusTwoH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionCoherentZNCanyZNAzeroShiftMinusTwoH);

  fInvariantMassDistributionCoherentZNCanyZNAzeroShiftMinusOneH = new TH1F("fInvariantMassDistributionCoherentZNCanyZNAzeroShiftMinusOneH", "fInvariantMassDistributionCoherentZNCanyZNAzeroShiftMinusOneH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionCoherentZNCanyZNAzeroShiftMinusOneH);

  fInvariantMassDistributionCoherentZNCanyZNAzeroShiftPlusOneH = new TH1F("fInvariantMassDistributionCoherentZNCanyZNAzeroShiftPlusOneH", "fInvariantMassDistributionCoherentZNCanyZNAzeroShiftPlusOneH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionCoherentZNCanyZNAzeroShiftPlusOneH);

  fInvariantMassDistributionCoherentZNCanyZNAzeroShiftPlusTwoH = new TH1F("fInvariantMassDistributionCoherentZNCanyZNAzeroShiftPlusTwoH", "fInvariantMassDistributionCoherentZNCanyZNAzeroShiftPlusTwoH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionCoherentZNCanyZNAzeroShiftPlusTwoH);

  fInvariantMassDistributionCoherentZNCanyZNAzeroHv2 = new TH1F("fInvariantMassDistributionCoherentZNCanyZNAzeroHv2", "fInvariantMassDistributionCoherentZNCanyZNAzeroHv2", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionCoherentZNCanyZNAzeroHv2);

  for( Int_t iRapidity = 0; iRapidity < 3; iRapidity++ ){
    fInvariantMassDistributionCoherentZNCanyZNAzeroRapidityHv2[iRapidity] = new TH1F(
              Form("fInvariantMassDistributionCoherentZNCanyZNAzeroRapidityHv2_%d", iRapidity),
              Form("fInvariantMassDistributionCoherentZNCanyZNAzeroRapidityHv2_%d", iRapidity),
              2000, 0, 20
              );
    fOutputList->Add(fInvariantMassDistributionCoherentZNCanyZNAzeroRapidityHv2[iRapidity]);
  }

  fInvariantMassDistributionCoherentZNCanyZNAanyH = new TH1F("fInvariantMassDistributionCoherentZNCanyZNAanyH", "fInvariantMassDistributionCoherentZNCanyZNAanyH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionCoherentZNCanyZNAanyH);

  fInvariantMassDistributionCoherentZNCanyZNAanyShiftMinusTwoH = new TH1F("fInvariantMassDistributionCoherentZNCanyZNAanyShiftMinusTwoH", "fInvariantMassDistributionCoherentZNCanyZNAanyShiftMinusTwoH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionCoherentZNCanyZNAanyShiftMinusTwoH);

  fInvariantMassDistributionCoherentZNCanyZNAanyShiftMinusOneH = new TH1F("fInvariantMassDistributionCoherentZNCanyZNAanyShiftMinusOneH", "fInvariantMassDistributionCoherentZNCanyZNAanyShiftMinusOneH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionCoherentZNCanyZNAanyShiftMinusOneH);

  fInvariantMassDistributionCoherentZNCanyZNAanyShiftPlusOneH = new TH1F("fInvariantMassDistributionCoherentZNCanyZNAanyShiftPlusOneH", "fInvariantMassDistributionCoherentZNCanyZNAanyShiftPlusOneH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionCoherentZNCanyZNAanyShiftPlusOneH);

  fInvariantMassDistributionCoherentZNCanyZNAanyShiftPlusTwoH = new TH1F("fInvariantMassDistributionCoherentZNCanyZNAanyShiftPlusTwoH", "fInvariantMassDistributionCoherentZNCanyZNAanyShiftPlusTwoH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionCoherentZNCanyZNAanyShiftPlusTwoH);

  fInvariantMassDistributionCoherentZNCanyZNAanyHv2 = new TH1F("fInvariantMassDistributionCoherentZNCanyZNAanyHv2", "fInvariantMassDistributionCoherentZNCanyZNAanyHv2", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionCoherentZNCanyZNAanyHv2);

  for( Int_t iRapidity = 0; iRapidity < 3; iRapidity++ ){
    fInvariantMassDistributionCoherentZNCanyZNAanyRapidityHv2[iRapidity] = new TH1F(
              Form("fInvariantMassDistributionCoherentZNCanyZNAanyRapidityHv2_%d", iRapidity),
              Form("fInvariantMassDistributionCoherentZNCanyZNAanyRapidityHv2_%d", iRapidity),
              2000, 0, 20
              );
    fOutputList->Add(fInvariantMassDistributionCoherentZNCanyZNAanyRapidityHv2[iRapidity]);
  }

  fInvariantMassDistributionIncoherentZNCzeroZNAzeroH = new TH1F("fInvariantMassDistributionIncoherentZNCzeroZNAzeroH", "fInvariantMassDistributionIncoherentZNCzeroZNAzeroH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionIncoherentZNCzeroZNAzeroH);

  fInvariantMassDistributionIncoherentZNCzeroZNAzeroShiftMinusTwoH = new TH1F("fInvariantMassDistributionIncoherentZNCzeroZNAzeroShiftMinusTwoH", "fInvariantMassDistributionIncoherentZNCzeroZNAzeroShiftMinusTwoH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionIncoherentZNCzeroZNAzeroShiftMinusTwoH);

  fInvariantMassDistributionIncoherentZNCzeroZNAzeroShiftMinusOneH = new TH1F("fInvariantMassDistributionIncoherentZNCzeroZNAzeroShiftMinusOneH", "fInvariantMassDistributionIncoherentZNCzeroZNAzeroShiftMinusOneH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionIncoherentZNCzeroZNAzeroShiftMinusOneH);

  fInvariantMassDistributionIncoherentZNCzeroZNAzeroShiftPlusOneH = new TH1F("fInvariantMassDistributionIncoherentZNCzeroZNAzeroShiftPlusOneH", "fInvariantMassDistributionIncoherentZNCzeroZNAzeroShiftPlusOneH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionIncoherentZNCzeroZNAzeroShiftPlusOneH);

  fInvariantMassDistributionIncoherentZNCzeroZNAzeroShiftPlusTwoH = new TH1F("fInvariantMassDistributionIncoherentZNCzeroZNAzeroShiftPlusTwoH", "fInvariantMassDistributionIncoherentZNCzeroZNAzeroShiftPlusTwoH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionIncoherentZNCzeroZNAzeroShiftPlusTwoH);

  fInvariantMassDistributionIncoherentZNCzeroZNAzeroHv2 = new TH1F("fInvariantMassDistributionIncoherentZNCzeroZNAzeroHv2", "fInvariantMassDistributionIncoherentZNCzeroZNAzeroHv2", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionIncoherentZNCzeroZNAzeroHv2);

  for( Int_t iRapidity = 0; iRapidity < 3; iRapidity++ ){
    fInvariantMassDistributionIncoherentZNCzeroZNAzeroRapidityHv2[iRapidity] = new TH1F(
              Form("fInvariantMassDistributionIncoherentZNCzeroZNAzeroRapidityHv2_%d", iRapidity),
              Form("fInvariantMassDistributionIncoherentZNCzeroZNAzeroRapidityHv2_%d", iRapidity),
              2000, 0, 20
              );
    fOutputList->Add(fInvariantMassDistributionIncoherentZNCzeroZNAzeroRapidityHv2[iRapidity]);
  }

  fInvariantMassDistributionIncoherentZNCzeroZNAanyH = new TH1F("fInvariantMassDistributionIncoherentZNCzeroZNAanyH", "fInvariantMassDistributionIncoherentZNCzeroZNAanyH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionIncoherentZNCzeroZNAanyH);

  fInvariantMassDistributionIncoherentZNCzeroZNAanyShiftMinusTwoH = new TH1F("fInvariantMassDistributionIncoherentZNCzeroZNAanyShiftMinusTwoH", "fInvariantMassDistributionIncoherentZNCzeroZNAanyShiftMinusTwoH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionIncoherentZNCzeroZNAanyShiftMinusTwoH);

  fInvariantMassDistributionIncoherentZNCzeroZNAanyShiftMinusOneH = new TH1F("fInvariantMassDistributionIncoherentZNCzeroZNAanyShiftMinusOneH", "fInvariantMassDistributionIncoherentZNCzeroZNAanyShiftMinusOneH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionIncoherentZNCzeroZNAanyShiftMinusOneH);

  fInvariantMassDistributionIncoherentZNCzeroZNAanyShiftPlusOneH = new TH1F("fInvariantMassDistributionIncoherentZNCzeroZNAanyShiftPlusOneH", "fInvariantMassDistributionIncoherentZNCzeroZNAanyShiftPlusOneH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionIncoherentZNCzeroZNAanyShiftPlusOneH);

  fInvariantMassDistributionIncoherentZNCzeroZNAanyShiftPlusTwoH = new TH1F("fInvariantMassDistributionIncoherentZNCzeroZNAanyShiftPlusTwoH", "fInvariantMassDistributionIncoherentZNCzeroZNAanyShiftPlusTwoH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionIncoherentZNCzeroZNAanyShiftPlusTwoH);

  fInvariantMassDistributionIncoherentZNCzeroZNAanyHv2 = new TH1F("fInvariantMassDistributionIncoherentZNCzeroZNAanyHv2", "fInvariantMassDistributionIncoherentZNCzeroZNAanyHv2", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionIncoherentZNCzeroZNAanyHv2);

  for( Int_t iRapidity = 0; iRapidity < 3; iRapidity++ ){
    fInvariantMassDistributionIncoherentZNCzeroZNAanyRapidityHv2[iRapidity] = new TH1F(
              Form("fInvariantMassDistributionIncoherentZNCzeroZNAanyRapidityHv2_%d", iRapidity),
              Form("fInvariantMassDistributionIncoherentZNCzeroZNAanyRapidityHv2_%d", iRapidity),
              2000, 0, 20
              );
    fOutputList->Add(fInvariantMassDistributionIncoherentZNCzeroZNAanyRapidityHv2[iRapidity]);
  }

  fInvariantMassDistributionIncoherentZNCanyZNAzeroH = new TH1F("fInvariantMassDistributionIncoherentZNCanyZNAzeroH", "fInvariantMassDistributionIncoherentZNCanyZNAzeroH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionIncoherentZNCanyZNAzeroH);

  fInvariantMassDistributionIncoherentZNCanyZNAzeroShiftMinusTwoH = new TH1F("fInvariantMassDistributionIncoherentZNCanyZNAzeroShiftMinusTwoH", "fInvariantMassDistributionIncoherentZNCanyZNAzeroShiftMinusTwoH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionIncoherentZNCanyZNAzeroShiftMinusTwoH);

  fInvariantMassDistributionIncoherentZNCanyZNAzeroShiftMinusOneH = new TH1F("fInvariantMassDistributionIncoherentZNCanyZNAzeroShiftMinusOneH", "fInvariantMassDistributionIncoherentZNCanyZNAzeroShiftMinusOneH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionIncoherentZNCanyZNAzeroShiftMinusOneH);

  fInvariantMassDistributionIncoherentZNCanyZNAzeroShiftPlusOneH = new TH1F("fInvariantMassDistributionIncoherentZNCanyZNAzeroShiftPlusOneH", "fInvariantMassDistributionIncoherentZNCanyZNAzeroShiftPlusOneH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionIncoherentZNCanyZNAzeroShiftPlusOneH);

  fInvariantMassDistributionIncoherentZNCanyZNAzeroShiftPlusTwoH = new TH1F("fInvariantMassDistributionIncoherentZNCanyZNAzeroShiftPlusTwoH", "fInvariantMassDistributionIncoherentZNCanyZNAzeroShiftPlusTwoH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionIncoherentZNCanyZNAzeroShiftPlusTwoH);

  fInvariantMassDistributionIncoherentZNCanyZNAzeroHv2 = new TH1F("fInvariantMassDistributionIncoherentZNCanyZNAzeroHv2", "fInvariantMassDistributionIncoherentZNCanyZNAzeroHv2", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionIncoherentZNCanyZNAzeroHv2);

  for( Int_t iRapidity = 0; iRapidity < 3; iRapidity++ ){
    fInvariantMassDistributionIncoherentZNCanyZNAzeroRapidityHv2[iRapidity] = new TH1F(
              Form("fInvariantMassDistributionIncoherentZNCanyZNAzeroRapidityHv2_%d", iRapidity),
              Form("fInvariantMassDistributionIncoherentZNCanyZNAzeroRapidityHv2_%d", iRapidity),
              2000, 0, 20
              );
    fOutputList->Add(fInvariantMassDistributionIncoherentZNCanyZNAzeroRapidityHv2[iRapidity]);
  }

  fInvariantMassDistributionIncoherentZNCanyZNAanyH = new TH1F("fInvariantMassDistributionIncoherentZNCanyZNAanyH", "fInvariantMassDistributionIncoherentZNCanyZNAanyH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionIncoherentZNCanyZNAanyH);

  fInvariantMassDistributionIncoherentZNCanyZNAanyShiftMinusTwoH = new TH1F("fInvariantMassDistributionIncoherentZNCanyZNAanyShiftMinusTwoH", "fInvariantMassDistributionIncoherentZNCanyZNAanyShiftMinusTwoH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionIncoherentZNCanyZNAanyShiftMinusTwoH);

  fInvariantMassDistributionIncoherentZNCanyZNAanyShiftMinusOneH = new TH1F("fInvariantMassDistributionIncoherentZNCanyZNAanyShiftMinusOneH", "fInvariantMassDistributionIncoherentZNCanyZNAanyShiftMinusOneH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionIncoherentZNCanyZNAanyShiftMinusOneH);

  fInvariantMassDistributionIncoherentZNCanyZNAanyShiftPlusOneH = new TH1F("fInvariantMassDistributionIncoherentZNCanyZNAanyShiftPlusOneH", "fInvariantMassDistributionIncoherentZNCanyZNAanyShiftPlusOneH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionIncoherentZNCanyZNAanyShiftPlusOneH);

  fInvariantMassDistributionIncoherentZNCanyZNAanyShiftPlusTwoH = new TH1F("fInvariantMassDistributionIncoherentZNCanyZNAanyShiftPlusTwoH", "fInvariantMassDistributionIncoherentZNCanyZNAanyShiftPlusTwoH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionIncoherentZNCanyZNAanyShiftPlusTwoH);

  fInvariantMassDistributionIncoherentZNCanyZNAanyHv2 = new TH1F("fInvariantMassDistributionIncoherentZNCanyZNAanyHv2", "fInvariantMassDistributionIncoherentZNCanyZNAanyHv2", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionIncoherentZNCanyZNAanyHv2);

  for( Int_t iRapidity = 0; iRapidity < 3; iRapidity++ ){
    fInvariantMassDistributionIncoherentZNCanyZNAanyRapidityHv2[iRapidity] = new TH1F(
              Form("fInvariantMassDistributionIncoherentZNCanyZNAanyRapidityHv2_%d", iRapidity),
              Form("fInvariantMassDistributionIncoherentZNCanyZNAanyRapidityHv2_%d", iRapidity),
              2000, 0, 20
              );
    fOutputList->Add(fInvariantMassDistributionIncoherentZNCanyZNAanyRapidityHv2[iRapidity]);
  }

  fDimuonPtDistributionZNCzeroZNAzeroH = new TH1F("fDimuonPtDistributionZNCzeroZNAzeroH", "fDimuonPtDistributionZNCzeroZNAzeroH", 4000, 0, 20);
  fOutputList->Add(fDimuonPtDistributionZNCzeroZNAzeroH);

  fDimuonPtDistributionZNCzeroZNAanyH = new TH1F("fDimuonPtDistributionZNCzeroZNAanyH", "fDimuonPtDistributionZNCzeroZNAanyH", 4000, 0, 20);
  fOutputList->Add(fDimuonPtDistributionZNCzeroZNAanyH);

  fDimuonPtDistributionZNCanyZNAzeroH = new TH1F("fDimuonPtDistributionZNCanyZNAzeroH", "fDimuonPtDistributionZNCanyZNAzeroH", 4000, 0, 20);
  fOutputList->Add(fDimuonPtDistributionZNCanyZNAzeroH);

  fDimuonPtDistributionZNCanyZNAanyH = new TH1F("fDimuonPtDistributionZNCanyZNAanyH", "fDimuonPtDistributionZNCanyZNAanyH", 4000, 0, 20);
  fOutputList->Add(fDimuonPtDistributionZNCanyZNAanyH);

  fDimuonPtDistributionZNCzeroZNAzeroShiftPlusOneH = new TH1F("fDimuonPtDistributionZNCzeroZNAzeroShiftPlusOneH", "fDimuonPtDistributionZNCzeroZNAzeroShiftPlusOneH", 4000, 0.02, 20.02);
  fOutputList->Add(fDimuonPtDistributionZNCzeroZNAzeroShiftPlusOneH);

  fDimuonPtDistributionZNCzeroZNAanyShiftPlusOneH = new TH1F("fDimuonPtDistributionZNCzeroZNAanyShiftPlusOneH", "fDimuonPtDistributionZNCzeroZNAanyShiftPlusOneH", 4000, 0.02, 20.02);
  fOutputList->Add(fDimuonPtDistributionZNCzeroZNAanyShiftPlusOneH);

  fDimuonPtDistributionZNCanyZNAzeroShiftPlusOneH = new TH1F("fDimuonPtDistributionZNCanyZNAzeroShiftPlusOneH", "fDimuonPtDistributionZNCanyZNAzeroShiftPlusOneH", 4000, 0.02, 20.02);
  fOutputList->Add(fDimuonPtDistributionZNCanyZNAzeroShiftPlusOneH);

  fDimuonPtDistributionZNCanyZNAanyShiftPlusOneH = new TH1F("fDimuonPtDistributionZNCanyZNAanyShiftPlusOneH", "fDimuonPtDistributionZNCanyZNAanyShiftPlusOneH", 4000, 0.02, 20.02);
  fOutputList->Add(fDimuonPtDistributionZNCanyZNAanyShiftPlusOneH);

  fDimuonPtDistributionZNCzeroZNAzeroHv2 = new TH1F("fDimuonPtDistributionZNCzeroZNAzeroHv2", "fDimuonPtDistributionZNCzeroZNAzeroHv2", 4000, 0, 20);
  fOutputList->Add(fDimuonPtDistributionZNCzeroZNAzeroHv2);

  fDimuonPtDistributionZNCzeroZNAanyHv2 = new TH1F("fDimuonPtDistributionZNCzeroZNAanyHv2", "fDimuonPtDistributionZNCzeroZNAanyHv2", 4000, 0, 20);
  fOutputList->Add(fDimuonPtDistributionZNCzeroZNAanyHv2);

  fDimuonPtDistributionZNCanyZNAzeroHv2 = new TH1F("fDimuonPtDistributionZNCanyZNAzeroHv2", "fDimuonPtDistributionZNCanyZNAzeroHv2", 4000, 0, 20);
  fOutputList->Add(fDimuonPtDistributionZNCanyZNAzeroHv2);

  fDimuonPtDistributionZNCanyZNAanyHv2 = new TH1F("fDimuonPtDistributionZNCanyZNAanyHv2", "fDimuonPtDistributionZNCanyZNAanyHv2", 4000, 0, 20);
  fOutputList->Add(fDimuonPtDistributionZNCanyZNAanyHv2);

  for( Int_t iRapidity = 0; iRapidity < 3; iRapidity++ ){
    fDimuonPtDistributionZNCzeroZNAzeroRapidityHv2[iRapidity] = new TH1F(
              Form("fDimuonPtDistributionZNCzeroZNAzeroRapidityHv2_%d", iRapidity),
              Form("fDimuonPtDistributionZNCzeroZNAzeroRapidityHv2_%d", iRapidity),
              4000, 0, 20
              );
    fOutputList->Add(fDimuonPtDistributionZNCzeroZNAzeroRapidityHv2[iRapidity]);
  }

  for( Int_t iRapidity = 0; iRapidity < 3; iRapidity++ ){
    fDimuonPtDistributionZNCzeroZNAanyRapidityHv2[iRapidity] = new TH1F(
              Form("fDimuonPtDistributionZNCzeroZNAanyRapidityHv2_%d", iRapidity),
              Form("fDimuonPtDistributionZNCzeroZNAanyRapidityHv2_%d", iRapidity),
              4000, 0, 20
              );
    fOutputList->Add(fDimuonPtDistributionZNCzeroZNAanyRapidityHv2[iRapidity]);
  }

  for( Int_t iRapidity = 0; iRapidity < 3; iRapidity++ ){
    fDimuonPtDistributionZNCanyZNAzeroRapidityHv2[iRapidity] = new TH1F(
              Form("fDimuonPtDistributionZNCanyZNAzeroRapidityHv2_%d", iRapidity),
              Form("fDimuonPtDistributionZNCanyZNAzeroRapidityHv2_%d", iRapidity),
              4000, 0, 20
              );
    fOutputList->Add(fDimuonPtDistributionZNCanyZNAzeroRapidityHv2[iRapidity]);
  }

  for( Int_t iRapidity = 0; iRapidity < 3; iRapidity++ ){
    fDimuonPtDistributionZNCanyZNAanyRapidityHv2[iRapidity] = new TH1F(
              Form("fDimuonPtDistributionZNCanyZNAanyRapidityHv2_%d", iRapidity),
              Form("fDimuonPtDistributionZNCanyZNAanyRapidityHv2_%d", iRapidity),
              4000, 0, 20
              );
    fOutputList->Add(fDimuonPtDistributionZNCanyZNAanyRapidityHv2[iRapidity]);
  }


  fDimuonPtDistributionCoherentZNCzeroZNAzeroH = new TH1F("fDimuonPtDistributionCoherentZNCzeroZNAzeroH", "fDimuonPtDistributionCoherentZNCzeroZNAzeroH", 4000, 0, 20);
  fOutputList->Add(fDimuonPtDistributionCoherentZNCzeroZNAzeroH);

  fDimuonPtDistributionCoherentZNCzeroZNAanyH = new TH1F("fDimuonPtDistributionCoherentZNCzeroZNAanyH", "fDimuonPtDistributionCoherentZNCzeroZNAanyH", 4000, 0, 20);
  fOutputList->Add(fDimuonPtDistributionCoherentZNCzeroZNAanyH);

  fDimuonPtDistributionCoherentZNCanyZNAzeroH = new TH1F("fDimuonPtDistributionCoherentZNCanyZNAzeroH", "fDimuonPtDistributionCoherentZNCanyZNAzeroH", 4000, 0, 20);
  fOutputList->Add(fDimuonPtDistributionCoherentZNCanyZNAzeroH);

  fDimuonPtDistributionCoherentZNCanyZNAanyH = new TH1F("fDimuonPtDistributionCoherentZNCanyZNAanyH", "fDimuonPtDistributionCoherentZNCanyZNAanyH", 4000, 0, 20);
  fOutputList->Add(fDimuonPtDistributionCoherentZNCanyZNAanyH);

  fDimuonPtDistributionIncoherentZNCzeroZNAzeroH = new TH1F("fDimuonPtDistributionIncoherentZNCzeroZNAzeroH", "fDimuonPtDistributionIncoherentZNCzeroZNAzeroH", 4000, 0, 20);
  fOutputList->Add(fDimuonPtDistributionIncoherentZNCzeroZNAzeroH);

  fDimuonPtDistributionIncoherentZNCzeroZNAanyH = new TH1F("fDimuonPtDistributionIncoherentZNCzeroZNAanyH", "fDimuonPtDistributionIncoherentZNCzeroZNAanyH", 4000, 0, 20);
  fOutputList->Add(fDimuonPtDistributionIncoherentZNCzeroZNAanyH);

  fDimuonPtDistributionIncoherentZNCanyZNAzeroH = new TH1F("fDimuonPtDistributionIncoherentZNCanyZNAzeroH", "fDimuonPtDistributionIncoherentZNCanyZNAzeroH", 4000, 0, 20);
  fOutputList->Add(fDimuonPtDistributionIncoherentZNCanyZNAzeroH);

  fDimuonPtDistributionIncoherentZNCanyZNAanyH = new TH1F("fDimuonPtDistributionIncoherentZNCanyZNAanyH", "fDimuonPtDistributionIncoherentZNCanyZNAanyH", 4000, 0, 20);
  fOutputList->Add(fDimuonPtDistributionIncoherentZNCanyZNAanyH);



  /* - Here starts the list of histograms needed for the analysis of the J/Psi's
     - polarization.
     -
   */
  fAngularDistribOfPositiveMuonRestFrameJPsiH = new TH1F("fAngularDistribOfPositiveMuonRestFrameJPsiH", "fAngularDistribOfPositiveMuonRestFrameJPsiH", 1000, -1., 1.);
  fOutputList->Add(fAngularDistribOfPositiveMuonRestFrameJPsiH);

  fAngularDistribOfNegativeMuonRestFrameJPsiH = new TH1F("fAngularDistribOfNegativeMuonRestFrameJPsiH", "fAngularDistribOfNegativeMuonRestFrameJPsiH", 1000, -1., 1.);
  fOutputList->Add(fAngularDistribOfNegativeMuonRestFrameJPsiH);

  fCheckHelicityRestFrameJPsiH = new TH1F("fCheckHelicityRestFrameJPsiH", "fCheckHelicityRestFrameJPsiH", 100000, -50., 50.);
  fOutputList->Add(fCheckHelicityRestFrameJPsiH);

  for(Int_t iRapidityBin = 0; iRapidityBin < 8; iRapidityBin++ ){
    fThetaDistribOfPositiveMuonRestFrameJPsiRapidityBinH[iRapidityBin] = new TH1F(
                Form("fThetaDistribOfPositiveMuonRestFrameJPsiRapidityBinH_%d", iRapidityBin),
                Form("fThetaDistribOfPositiveMuonRestFrameJPsiRapidityBinH_%d", iRapidityBin),
                1000, -1., 1.
              );
    fOutputList->Add(fThetaDistribOfPositiveMuonRestFrameJPsiRapidityBinH[iRapidityBin]);
  }

  fCosThetaHelicityFrameJPsiH = new TH1F("fCosThetaHelicityFrameJPsiH", "fCosThetaHelicityFrameJPsiH", 1000, -1., 1.);
  fOutputList->Add(fCosThetaHelicityFrameJPsiH);

  fCosThetaCollinsSoperFrameJPsiH = new TH1F("fCosThetaCollinsSoperFrameJPsiH", "fCosThetaCollinsSoperFrameJPsiH", 1000, -1., 1.);
  fOutputList->Add(fCosThetaCollinsSoperFrameJPsiH);

  fPhiHelicityFrameJPsiH = new TH1F("fPhiHelicityFrameJPsiH", "fPhiHelicityFrameJPsiH", 4000, -4., 4.);
  fOutputList->Add(fPhiHelicityFrameJPsiH);

  fPhiCollinsSoperFrameJPsiH = new TH1F("fPhiCollinsSoperFrameJPsiH", "fPhiCollinsSoperFrameJPsiH", 4000, -4., 4.);
  fOutputList->Add(fPhiCollinsSoperFrameJPsiH);

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

  for(Int_t iRapidityBin = 0; iRapidityBin < 10; iRapidityBin++ ){
    fCosThetaHelicityFrameJPsiTenRapidityBinsH[iRapidityBin] = new TH1F(
                Form("fCosThetaHelicityFrameJPsiTenRapidityBinsH_%d", iRapidityBin),
                Form("fCosThetaHelicityFrameJPsiTenRapidityBinsH_%d", iRapidityBin),
                1000, -1., 1.
              );
    fOutputList->Add(fCosThetaHelicityFrameJPsiTenRapidityBinsH[iRapidityBin]);
  }

  for(Int_t iRapidityBin = 0; iRapidityBin < 10; iRapidityBin++ ){
    fCosThetaCollinsSoperFrameJPsiTenRapidityBinsH[iRapidityBin] = new TH1F(
                Form("fCosThetaCollinsSoperFrameJPsiTenRapidityBinsH_%d", iRapidityBin),
                Form("fCosThetaCollinsSoperFrameJPsiTenRapidityBinsH_%d", iRapidityBin),
                1000, -1., 1.
              );
    fOutputList->Add(fCosThetaCollinsSoperFrameJPsiTenRapidityBinsH[iRapidityBin]);
  }

  for(Int_t iRapidityBin = 0; iRapidityBin < 10; iRapidityBin++ ){
    fPhiHelicityFrameJPsiTenRapidityBinsH[iRapidityBin] = new TH1F(
                Form("fPhiHelicityFrameJPsiTenRapidityBinsH_%d", iRapidityBin),
                Form("fPhiHelicityFrameJPsiTenRapidityBinsH_%d", iRapidityBin),
                4000, -4., 4.
              );
    fOutputList->Add(fPhiHelicityFrameJPsiTenRapidityBinsH[iRapidityBin]);
  }

  for(Int_t iRapidityBin = 0; iRapidityBin < 10; iRapidityBin++ ){
    fPhiCollinsSoperFrameJPsiTenRapidityBinsH[iRapidityBin] = new TH1F(
                Form("fPhiCollinsSoperFrameJPsiTenRapidityBinsH_%d", iRapidityBin),
                Form("fPhiCollinsSoperFrameJPsiTenRapidityBinsH_%d", iRapidityBin),
                4000, -4., 4.
              );
    fOutputList->Add(fPhiCollinsSoperFrameJPsiTenRapidityBinsH[iRapidityBin]);
  }

  for(Int_t iCosThetaBins = 0; iCosThetaBins < 10; iCosThetaBins++ ){
    fInvariantMassDistributionInBinsOfCosThetaHelicityFrameH[iCosThetaBins] = new TH1F(
                Form("fInvariantMassDistributionInBinsOfCosThetaHelicityFrameH_%d", iCosThetaBins),
                Form("fInvariantMassDistributionInBinsOfCosThetaHelicityFrameH_%d", iCosThetaBins),
                2000, 0, 20
                );
    fOutputList->Add(fInvariantMassDistributionInBinsOfCosThetaHelicityFrameH[iCosThetaBins]);
  }

  fInvariantMassDistributionBinsOfCosThetaAndPhiHelicityFrameH =
        new TH2F( "fInvariantMassDistributionBinsOfCosThetaAndPhiHelicityFrameH",
                  "fInvariantMassDistributionBinsOfCosThetaAndPhiHelicityFrameH",
                  80, -1, 1, 80, -4, 4
                  );
  fOutputList->Add(fInvariantMassDistributionBinsOfCosThetaAndPhiHelicityFrameH);

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

  /* - Invariant mass distributions for signal extraction for POLARISATION.
     - The usage will be:    histo[CosTheta][Phi];
     -
   */
  fInvariantMassDistributionForSignalExtractionHelicityFrameH = new TH1F**[10];
  for( Int_t iCosTheta = 0; iCosTheta < 10; iCosTheta++ ){
    fInvariantMassDistributionForSignalExtractionHelicityFrameH[iCosTheta] = new TH1F*[10];
    for( Int_t iPhi = 0; iPhi < 10; iPhi++ ){
      fInvariantMassDistributionForSignalExtractionHelicityFrameH[iCosTheta][iPhi] =
          new TH1F( Form("fInvariantMassDistributionForSignalExtractionHelicityFrameH_%d_%d", iCosTheta, iPhi),
                    Form("fInvariantMassDistributionForSignalExtractionHelicityFrameH_%d_%d", iCosTheta, iPhi),
                    2000, 0, 20
                    );
      fOutputList->Add(fInvariantMassDistributionForSignalExtractionHelicityFrameH[iCosTheta][iPhi]);
    }
  }

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

  for(Int_t iCosThetaBins = 0; iCosThetaBins < 40; iCosThetaBins++ ){
    fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameTryingH[iCosThetaBins] = new TH1F(
                Form("fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameTryingH_%d", iCosThetaBins),
                Form("fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameTryingH_%d", iCosThetaBins),
                2000, 0, 20
                );
    fOutputList->Add(fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameTryingH[iCosThetaBins]);
  }

  for(Int_t iPhiBins = 0; iPhiBins < 50; iPhiBins++ ){
    fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameTryingH[iPhiBins] = new TH1F(
                Form("fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameTryingH_%d", iPhiBins),
                Form("fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameTryingH_%d", iPhiBins),
                2000, 0, 20
                );
    fOutputList->Add(fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameTryingH[iPhiBins]);
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

  fInvariantMassDistributionForSignalExtractionHelicityFrameMyBinning0N0NH = new TH1F**[7];
  for( Int_t iCosTheta = 0; iCosTheta < 7; iCosTheta++ ){
    fInvariantMassDistributionForSignalExtractionHelicityFrameMyBinning0N0NH[iCosTheta] = new TH1F*[20];
    for( Int_t iPhi = 0; iPhi < 20; iPhi++ ){
      fInvariantMassDistributionForSignalExtractionHelicityFrameMyBinning0N0NH[iCosTheta][iPhi] =
          new TH1F( Form("fInvariantMassDistributionForSignalExtractionHelicityFrameMyBinning0N0NH_%d_%d", iCosTheta, iPhi),
                    Form("fInvariantMassDistributionForSignalExtractionHelicityFrameMyBinning0N0NH_%d_%d", iCosTheta, iPhi),
                    2000, 0, 20
                    );
      fOutputList->Add(fInvariantMassDistributionForSignalExtractionHelicityFrameMyBinning0N0NH[iCosTheta][iPhi]);
    }
  }

  fInvariantMassDistributionForSignalExtractionCsFrameMyBinning0N0NH = new TH1F**[7];
  for( Int_t iCosTheta = 0; iCosTheta < 7; iCosTheta++ ){
    fInvariantMassDistributionForSignalExtractionCsFrameMyBinning0N0NH[iCosTheta] = new TH1F*[20];
    for( Int_t iPhi = 0; iPhi < 20; iPhi++ ){
      fInvariantMassDistributionForSignalExtractionCsFrameMyBinning0N0NH[iCosTheta][iPhi] =
          new TH1F( Form("fInvariantMassDistributionForSignalExtractionCsFrameMyBinning0N0NH_%d_%d", iCosTheta, iPhi),
                    Form("fInvariantMassDistributionForSignalExtractionCsFrameMyBinning0N0NH_%d_%d", iCosTheta, iPhi),
                    2000, 0, 20
                    );
      fOutputList->Add(fInvariantMassDistributionForSignalExtractionCsFrameMyBinning0N0NH[iCosTheta][iPhi]);
    }
  }

  fInvariantMassDistributionStrictPtH = new TH1F("fInvariantMassDistributionStrictPtH", "fInvariantMassDistributionStrictPtH", 8000, 0, 40);
  fOutputList->Add(fInvariantMassDistributionStrictPtH);

  /* - My variable binning
   * - for 1D analysis.
   */
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

  for(Int_t iCosThetaBins = 0; iCosThetaBins < 9; iCosThetaBins++ ){
    fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameMyStrictVariableBinningH[iCosThetaBins] = new TH1F(
                Form("fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameMyStrictVariableBinningH_%d", iCosThetaBins),
                Form("fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameMyStrictVariableBinningH_%d", iCosThetaBins),
                2000, 0, 20
                );
    fOutputList->Add(fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameMyStrictVariableBinningH[iCosThetaBins]);
  }

  for(Int_t iPhiBins = 0; iPhiBins < 15; iPhiBins++ ){
    fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameMyStrictVariableBinningH[iPhiBins] = new TH1F(
                Form("fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameMyStrictVariableBinningH_%d", iPhiBins),
                Form("fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameMyStrictVariableBinningH_%d", iPhiBins),
                2000, 0, 20
                );
    fOutputList->Add(fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameMyStrictVariableBinningH[iPhiBins]);
  }

  for(Int_t iCosThetaBins = 0; iCosThetaBins < 17; iCosThetaBins++ ){
    fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameMySeventeenBinsVariableBinningH[iCosThetaBins] = new TH1F(
                Form("fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameMySeventeenBinsVariableBinningH_%d", iCosThetaBins),
                Form("fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameMySeventeenBinsVariableBinningH_%d", iCosThetaBins),
                2000, 0, 20
                );
    fOutputList->Add(fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameMySeventeenBinsVariableBinningH[iCosThetaBins]);
  }

  fCosThetaHelicityFrameJPsiAlreadyCorrectedH = new TH1F("fCosThetaHelicityFrameJPsiAlreadyCorrectedH", "fCosThetaHelicityFrameJPsiAlreadyCorrectedH", 100, -1., 1.);
  fOutputList->Add(fCosThetaHelicityFrameJPsiAlreadyCorrectedH);

  for(Int_t iCosThetaBins = 0; iCosThetaBins < 50; iCosThetaBins++ ){
    fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameForAlreadyCorrectedFiftyH[iCosThetaBins] = new TH1F(
                Form("fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameForAlreadyCorrectedFiftyH_%d", iCosThetaBins),
                Form("fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameForAlreadyCorrectedFiftyH_%d", iCosThetaBins),
                2000, 0, 20
                );
    fOutputList->Add(fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameForAlreadyCorrectedFiftyH[iCosThetaBins]);
  }

  for(Int_t iCosThetaBins = 0; iCosThetaBins < 100; iCosThetaBins++ ){
    fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameForAlreadyCorrectedHundredH[iCosThetaBins] = new TH1F(
                Form("fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameForAlreadyCorrectedHundredH_%d", iCosThetaBins),
                Form("fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameForAlreadyCorrectedHundredH_%d", iCosThetaBins),
                2000, 0, 20
                );
    fOutputList->Add(fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameForAlreadyCorrectedHundredH[iCosThetaBins]);
  }

  fDcaAgainstPtOfVectorMesonH = new TH2F("fDcaAgainstPtOfVectorMesonH", "fDcaAgainstPtOfVectorMesonH", 1000, 0, 10, 4000, 0, 40);
  fOutputList->Add(fDcaAgainstPtOfVectorMesonH);

  fInvariantMassDistributionStrictPtStrictDcaH = new TH1F("fInvariantMassDistributionStrictPtStrictDcaH", "fInvariantMassDistributionStrictPtStrictDcaH", 8000, 0, 40);
  fOutputList->Add(fInvariantMassDistributionStrictPtStrictDcaH);

  fInvariantMassDistributionVsPtStrictDcaH = new TH2F("fInvariantMassDistributionVsPtStrictDcaH", "fInvariantMassDistributionVsPtStrictDcaH", 2000, 0, 20, 1000, 0, 10);
  fOutputList->Add(fInvariantMassDistributionVsPtStrictDcaH);

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

  for(Int_t iCosThetaBins = 0; iCosThetaBins < 25; iCosThetaBins++ ){
    fInvariantMassDistributionOnlyCosThetaHeFrameTwentyfiveBins0N0NH[iCosThetaBins] = new TH1F(
                Form("fInvariantMassDistributionOnlyCosThetaHeFrameTwentyfiveBins0N0NH_%d", iCosThetaBins),
                Form("fInvariantMassDistributionOnlyCosThetaHeFrameTwentyfiveBins0N0NH_%d", iCosThetaBins),
                2000, 0, 20
                );
    fOutputList->Add(fInvariantMassDistributionOnlyCosThetaHeFrameTwentyfiveBins0N0NH[iCosThetaBins]);
  }

  for(Int_t iPhiBins = 0; iPhiBins < 25; iPhiBins++ ){
    fInvariantMassDistributionOnlyPhiHeFrameTwentyfiveBins0N0NH[iPhiBins] = new TH1F(
                Form("fInvariantMassDistributionOnlyPhiHeFrameTwentyfiveBins0N0NH_%d", iPhiBins),
                Form("fInvariantMassDistributionOnlyPhiHeFrameTwentyfiveBins0N0NH_%d", iPhiBins),
                2000, 0, 20
                );
    fOutputList->Add(fInvariantMassDistributionOnlyPhiHeFrameTwentyfiveBins0N0NH[iPhiBins]);
  }

  for(Int_t iPhiBins = 0; iPhiBins < 25; iPhiBins++ ){
    fInvariantMassDistributionOnlyTildePhiHeFrameTwentyfiveBins0N0NH[iPhiBins] = new TH1F(
                Form("fInvariantMassDistributionOnlyTildePhiHeFrameTwentyfiveBins0N0NH_%d", iPhiBins),
                Form("fInvariantMassDistributionOnlyTildePhiHeFrameTwentyfiveBins0N0NH_%d", iPhiBins),
                2000, 0, 20
                );
    fOutputList->Add(fInvariantMassDistributionOnlyTildePhiHeFrameTwentyfiveBins0N0NH[iPhiBins]);
  }

  for(Int_t iCosThetaBins = 0; iCosThetaBins < 25; iCosThetaBins++ ){
    fInvariantMassDistributionOnlyCosThetaCsFrameTwentyfiveBins0N0NH[iCosThetaBins] = new TH1F(
                Form("fInvariantMassDistributionOnlyCosThetaCsFrameTwentyfiveBins0N0NH_%d", iCosThetaBins),
                Form("fInvariantMassDistributionOnlyCosThetaCsFrameTwentyfiveBins0N0NH_%d", iCosThetaBins),
                2000, 0, 20
                );
    fOutputList->Add(fInvariantMassDistributionOnlyCosThetaCsFrameTwentyfiveBins0N0NH[iCosThetaBins]);
  }

  for(Int_t iPhiBins = 0; iPhiBins < 25; iPhiBins++ ){
    fInvariantMassDistributionOnlyPhiCsFrameTwentyfiveBins0N0NH[iPhiBins] = new TH1F(
                Form("fInvariantMassDistributionOnlyPhiCsFrameTwentyfiveBins0N0NH_%d", iPhiBins),
                Form("fInvariantMassDistributionOnlyPhiCsFrameTwentyfiveBins0N0NH_%d", iPhiBins),
                2000, 0, 20
                );
    fOutputList->Add(fInvariantMassDistributionOnlyPhiCsFrameTwentyfiveBins0N0NH[iPhiBins]);
  }

  for(Int_t iPhiBins = 0; iPhiBins < 25; iPhiBins++ ){
    fInvariantMassDistributionOnlyTildePhiCsFrameTwentyfiveBins0N0NH[iPhiBins] = new TH1F(
                Form("fInvariantMassDistributionOnlyTildePhiCsFrameTwentyfiveBins0N0NH_%d", iPhiBins),
                Form("fInvariantMassDistributionOnlyTildePhiCsFrameTwentyfiveBins0N0NH_%d", iPhiBins),
                2000, 0, 20
                );
    fOutputList->Add(fInvariantMassDistributionOnlyTildePhiCsFrameTwentyfiveBins0N0NH[iPhiBins]);
  }


  //_______________________________
  // - End of the function
  PostData(1, fOutputList);           // postdata will notify the analysis manager of changes / updates to the
  // fOutputList object. the manager will in the end take care of writing your output to file
  // so it needs to know what's in the output
}
//_____________________________________________________________________________
void AliAnalysisTaskUPCforward::NotifyRun()
{
  /// Set run number for cuts
  fMuonTrackCuts->SetRun(fInputHandler);
}
//_____________________________________________________________________________
void AliAnalysisTaskUPCforward::UserExec(Option_t *)
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
  fCounterH->Fill(iSelectionCounter); // AOD event found
  iSelectionCounter++;

  /* - Is it the right trigger?
     - In 2018 there were the following CMUP triggers:
     - CMUP11-B-NOPF-MUFAST,
     - CMUP26-B-NOPF-MUFAST,
     - CMUP6-B-NOPF-MUFAST.
     - At this 2nd step of fCounterH, all events are in fact proper events,
     - and with the correct trigger as well.
     -
     -
     - NEW: trigger class for the 2015 data.
     - The available classes in 2015 data were:
     - CMUP10-B-NOPF-MUFAST,
     - CMUP11-B-NOPF-MUFAST,
     - CMUP13-B-NOPF-MUFAST,
     - CTEST63-B-NOPF-MUFAST,
     - CTEST64-B-NOPF-MUFAST.
     - However, the CTESTs contain little to no data, so we
     - can just skip on them...
     - One thing I didn't know before, was that the CMUP11 trigger
     - class is in common with the 2018 dataset!!
     - This means that I could still obtain something
     - even with this string request for the 2015 data...
     -
   */
  TString trigger = fAOD->GetFiredTriggerClasses();
  if (    !(trigger.Contains("CMUP11-B-NOPF-MUFAST") ||
	          trigger.Contains("CMUP26-B-NOPF-MUFAST") ||
	          trigger.Contains("CMUP6-B-NOPF-MUFAST")  ||
            trigger.Contains("CMUP10-B-NOPF-MUFAST") ||
            trigger.Contains("CMUP13-B-NOPF-MUFAST")  )
          )  {
                    PostData(1, fOutputList);
                    return;
  }
  fCounterH->Fill(iSelectionCounter); // right trigger found
  iSelectionCounter++;


  /* - The following lines concern the LUMI computation.
   * - What is being done is that we fill the histograms
   * - with the number of events which pass the relative
   * - trigger conditions.
   * -
   */
  // Int_t  counterForTrigger = 0;
  fRunNum    = fAOD->GetRunNumber();
  if ( trigger.Contains("CMUP11-B-NOPF-MUFAST") )  {
    fTriggersVsRunH->Fill( 0.5, fRunNum );
    fRunNumberTriggerCMUP11ClassH        ->Fill(fRunNum);
    fRunNumberTriggerCMUP11ClassProperlyH->Fill( Form("%d", fRunNum) , 1 );
    // counterForTrigger++;
  }
  if ( trigger.Contains("CMUP26-B-NOPF-MUFAST") )  {
    fTriggersVsRunH->Fill( 1.5, fRunNum );
    fRunNumberTriggerCMUP26ClassH        ->Fill(fRunNum);
    fRunNumberTriggerCMUP26ClassProperlyH->Fill( Form("%d", fRunNum) , 1 );
    // counterForTrigger++;
  }
  if ( trigger.Contains("CMUP6-B-NOPF-MUFAST") )  {
    fTriggersVsRunH->Fill( 2.5, fRunNum );
    fRunNumberTriggerCMUP6ClassH        ->Fill(fRunNum);
    fRunNumberTriggerCMUP6ClassProperlyH->Fill( Form("%d", fRunNum) , 1 );
    // counterForTrigger++;
  }
  if ( trigger.Contains("CMUP10-B-NOPF-MUFAST") )  {
    fTriggersVsRunH->Fill( 3.5, fRunNum );
    fRunNumberTriggerCMUP10ClassH        ->Fill(fRunNum);
    fRunNumberTriggerCMUP10ClassProperlyH->Fill( Form("%d", fRunNum) , 1 );
    // counterForTrigger++;
  }
  if ( trigger.Contains("CMUP13-B-NOPF-MUFAST") )  {
    fTriggersVsRunH->Fill( 4.5, fRunNum );
    fRunNumberTriggerCMUP13ClassH        ->Fill(fRunNum);
    fRunNumberTriggerCMUP13ClassProperlyH->Fill( Form("%d", fRunNum) , 1 );
    // counterForTrigger++;
  }

  // if( counterForTrigger > 1 ) cout << "counterForTrigger = " << counterForTrigger << endl;

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
  // fRunNum    = fAOD->GetRunNumber();
  fTracklets = fAOD->GetTracklets()->GetNumberOfTracklets();

  /* - Trigger Inputs:
     - L0: ..... ;
     - L1: ..... .
   */
  fL0inputs = fAOD->GetHeader()->GetL0TriggerInputs();
  fL1inputs = fAOD->GetHeader()->GetL1TriggerInputs();

  /* - Past-future protection maps:
     - IR1: .... ;
     - IR2: .... .
   */
  fIR1Map = fAOD->GetHeader()->GetIRInt1InteractionMap();
  fIR2Map = fAOD->GetHeader()->GetIRInt2InteractionMap();

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

  /* - Reset Event information.
   * -
   */
  fZNAEnergy  = -8999;
  fZNCEnergy  = -8999;
  fZPAEnergy  = -8999;
  fZPCEnergy  = -8999;

  fZNAEnergy  = dataZDC->GetZNATowerEnergy()[0];
  fZNCEnergy  = dataZDC->GetZNCTowerEnergy()[0];
  fZPAEnergy  = dataZDC->GetZPATowerEnergy()[0];
  fZPCEnergy  = dataZDC->GetZPCTowerEnergy()[0];

  fZNATime    = dataZDC->GetZNATime();
  fZNCTime    = dataZDC->GetZNCTime();

  /* - Reset Event information.
   * -
   */
  for (Int_t i=0;i<4;i++) fZNATDC[i] = -999;
  for (Int_t i=0;i<4;i++) fZNCTDC[i] = -999;
  for (Int_t i=0;i<4;i++) fZPATDC[i] = -999;
  for (Int_t i=0;i<4;i++) fZPCTDC[i] = -999;

  for (Int_t i=0;i<4;i++) fZNATDC[i] = dataZDC->GetZNATDCm(i);
  for (Int_t i=0;i<4;i++) fZNCTDC[i] = dataZDC->GetZNCTDCm(i);
  for (Int_t i=0;i<4;i++) fZPATDC[i] = dataZDC->GetZPATDCm(i);
  for (Int_t i=0;i<4;i++) fZPCTDC[i] = dataZDC->GetZPCTDCm(i);

  /* - These lines are the calibration for the ZDC as provided by Evgeny Kryshen.
     -
   */
  Bool_t calibrated = 0;
  if ( fRunNum <= 245068 ) {
    calibrated = 1;
  } else if ( ( fRunNum > 245068 ) && ( fRunNum <  246995 ) ){
    calibrated = 0;
  } else {
    calibrated = 1;
  }
  // if ( fRunNum == 296509 ) calibrated = 1;
  // if ( fRunNum >  296689 ) calibrated = 1;
  // if ( fRunNum >  296695 ) calibrated = 0;
  // if ( fRunNum == 297219 ) calibrated = 1;
  // if ( fRunNum == 297221 ) calibrated = 1;
  // if ( fRunNum == 297415 ) calibrated = 1;

  if ( calibrated == 0 ) {
    if( fRunNum <= 246994 ) {
      fZNAEnergy *= (2500./250.);
      fZNCEnergy *= (2500./250.);
    }
  //   if( fRunNum >  246994 ) {
  //     fZNAEnergy *= (2500./190.);
  //     fZNCEnergy *= (2500./190.);
  //   }
  }

  /* - V0: we try to find the V0 object data in the nano-AOD. If we cannot,
     - we return, because there would be no way to actually select the events
     - otherwise! We are here, so we could even check if there is a discrepancy
     - between good events with and without V0's information. Or at least, this
     - is my impression of it (filling fCounterH). V0 information:
     - fV0ADecision: ..... ;
     - fV0CDecision: ..... .
     -
     -
     -
     - Plot the V0 variables to try to understand whether it is cells we are
     - talking about or boolean responses or something else altogether.
  */
  AliVVZERO *dataVZERO = dynamic_cast<AliVVZERO*>(fAOD->GetVZEROData());
  if(!dataVZERO) {
        PostData(1, fOutputList);
        return;
  }
  fCounterH->Fill(iSelectionCounter);
  iSelectionCounter++;
  fCounterH->Fill(12);


  fV0ADecision = dataVZERO->GetV0ADecision();
  fCounterH->Fill(13);
  fV0CDecision = dataVZERO->GetV0CDecision();
  fCounterH->Fill(14);


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
  //     fCounterH->Fill(15);
  // } else {
  //     // std::cout << "fVectorGoodRunNumbers does not contain: " << fRunNum << std::endl;
  //     fCounterH->Fill(16);
  //     PostData(1, fOutputList);
  //     return;
  // }
  // fCounterH->Fill(17);


  //_______________________________
  /* - NOTE: total run selection.
   * -
   */
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
  Bool_t checkIfGoodRun = kFALSE;
  for( Int_t iRunLHC18q = 0; iRunLHC18q < 128; iRunLHC18q++){
  // for( Int_t iRunLHC18q = 0; iRunLHC18q < 129; iRunLHC18q++){
  // for( Int_t iRunLHC18q = 0; iRunLHC18q < 125; iRunLHC18q++){
    if( fRunNum == listOfGoodRunNumbersLHC18q[iRunLHC18q] ) checkIfGoodRun = kTRUE;
  }
  for( Int_t iRunLHC18r = 0; iRunLHC18r <  97; iRunLHC18r++){
  // for( Int_t iRunLHC18r = 0; iRunLHC18r <  98; iRunLHC18r++){
  // for( Int_t iRunLHC18r = 0; iRunLHC18r <  82; iRunLHC18r++){
    if( fRunNum == listOfGoodRunNumbersLHC18r[iRunLHC18r] ) checkIfGoodRun = kTRUE;
  }
  for( Int_t iRunLHC15o = 0; iRunLHC15o < 136/*137*/; iRunLHC15o++){
  // for( Int_t iRunLHC15o = 0; iRunLHC15o < 134; iRunLHC15o++){
    if( fRunNum == listOfGoodRunNumbersLHC15o[iRunLHC15o] ) checkIfGoodRun = kTRUE;
  }
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
     -
     - Weird fact: this doesn't seem to work... I have changed it so that if
     - the single cell has recorded a signal (kTRUE) then it adds up to the
     - total number of cells. Hope for the best.
     -
     - I am an idiot!!!!!! I have to reset the variable everytime!!!!
   */
  fV0TotalNCells = 0;
  for(Int_t iV0Hits = 0; iV0Hits < 64; iV0Hits++) {
        fV0Hits[iV0Hits] = dataVZERO->GetBBFlag(iV0Hits);
        if(fV0Hits[iV0Hits] == kTRUE) {
              // if(iV0Hits < 32) fV0TotalNCells += fV0Hits[iV0Hits];
              if(iV0Hits < 32) fV0TotalNCells += 1;
        }
        // std::cout << "fV0Hits[iV0Hits = " << iV0Hits << ", fRunNum=" << fRunNum << "] = " << fV0Hits[iV0Hits] << endl;
        // std::cout << "fV0TotalNCells (fRunNum = " << fRunNum << ") = " << fV0TotalNCells << endl;
  }
  fCounterH->Fill(18);

  /* - AD: we try to find the AD object data in the nano-AOD. If we cannot,
     - we return, because there would be no way to actually select the events
     - otherwise! We are here, so we could even check if there is a discrepancy
     - between good events with and without AD's information. Or at least, this
     - is my impression of it (filling fCounterH). AD information:
     - fADADecision: small detector in ALICE, ADA and ADC at large distances;
     - fADCDecision: again, maybe check whether it is cells or boolean, same as V0.
  */
  // AD
  AliVAD *dataAD = dynamic_cast<AliVAD*>(fAOD->GetADData());
  fCounterH->Fill(19);
  if(dataAD) {
        fCounterH->Fill(iSelectionCounter);
        iSelectionCounter++;
        fCounterH->Fill(20);

        fADADecision = dataAD->GetADADecision();
        fADCDecision = dataAD->GetADCDecision();
        fCounterH->Fill(21);
  }
  fCounterH->Fill(22);

  // END EVENT DATA EXTRACTION
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
  /* - CMUP11-B triggers: I have to check with my supervisor, but this requirement
     - may have already been satisfied with the requirements for the trigger info
   */
  /* - Maximum 2 V0C cells fired:
     -
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
  /* - Empty V0C decision
   * - or at least in beam timing.
   */
  if( !(fV0CDecision == 0 || fV0CDecision == 1) ) {
       PostData(1, fOutputList);
       return;
  }
  /* - 0 tracklets in SPD
     - Is it like this?? Not too sure what fTracklets was!
   */
  // if(fTracklets != 0) {
  //      PostData(1, fOutputList);
  //      return;
  // }
  /* - Maximum 2 V0C cells fired.
     -
     - Trying a more readable and immediate approach.
   */
  // if( !(fV0TotalNCells < 2) ) {
  //      PostData(1, fOutputList);
  //      return;
  // }
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
    // get track
    // AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(iTrack));
    // if(!track) return;
    //
    // // is it a good muon track?
    // if(!track->IsMuonTrack()) continue;
    // if(!fMuonTrackCuts->IsSelected(track)) continue;

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

    /* -
     * - Systematics:
     * -
     * - Trigger Efficiency
     * -
     */
    if        ( fSetSingleMuonPt == 1 ) {
      if ( !(track[nGoodMuons]->Pt() > 0.85) ) {
      // track[nGoodMuons] = 0x0;
      continue;
      }
    } else if ( fSetSingleMuonPt == 2 ) {
      if ( !(track[nGoodMuons]->Pt() > 0.90) ) {
      // track[nGoodMuons] = 0x0;
      continue;
      }
    } else if ( fSetSingleMuonPt == 3 ) {
      if ( !(track[nGoodMuons]->Pt() > 0.95) ) {
      // track[nGoodMuons] = 0x0;
      continue;
      }
    } else if ( fSetSingleMuonPt == 4 ) {
      if ( !(track[nGoodMuons]->Pt() > 1.00) ) {
      // track[nGoodMuons] = 0x0;
      continue;
      }
    } else if ( fSetSingleMuonPt == 5 ) {
      if ( !(track[nGoodMuons]->Pt() > 1.05) ) {
      // track[nGoodMuons] = 0x0;
      continue;
      }
    } else if ( fSetSingleMuonPt == 6 ) {
      if ( !(track[nGoodMuons]->Pt() > 1.10) ) {
      // track[nGoodMuons] = 0x0;
      continue;
      }
    } else if ( fSetSingleMuonPt == 7 ) {
      if ( !(track[nGoodMuons]->Pt() > 1.15) ) {
      // track[nGoodMuons] = 0x0;
      continue;
      }
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

    // fill muon info
    // fEtaMuonH->Fill(track->Eta());
    // fRAbsMuonH->Fill(track->GetRAtAbsorberEnd());

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
        fEtaMuonH     ->Fill(track[iFilling]->Eta());
        fRAbsMuonH    ->Fill(track[iFilling]->GetRAtAbsorberEnd());
        fPtSingleMuonH->Fill(track[iFilling]->Pt());
  }
  // store muons
  fNumberMuonsH->Fill(nGoodMuons);
  fEntriesAgainstRunNumberH->Fill(fRunNum);
  /* - This is the last part of my try to obtain a proper RunNumbers histogram...
     -
   */
  fEntriesAgainstRunNumberProperlyH->Fill( Form("%d", fRunNum) , 1 );
  if (nGoodMuons>0) fCounterH->Fill(iSelectionCounter); // At least one good muon
  iSelectionCounter++;



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
  fInvariantMassDistributionExtendedH->Fill(possibleJPsi.Mag());

  /* - This is a TH2F histogram filled with DCA against the invariant mass of
     - the dimuon pair. This should be plotted as LEGO or as a heat map...
     - Is this the correct way of doing it?? I have to ask to my supervisor!!
     -
   */
  for(int indexMuonForDCA = 0; indexMuonForDCA < 2; indexMuonForDCA++) {
    fDcaAgainstInvariantMassH  ->Fill(possibleJPsi.Mag(), track[indexMuonForDCA]->DCA());
    fDcaAgainstPtOfVectorMesonH->Fill(possibleJPsi.Pt(),  track[indexMuonForDCA]->DCA());
  }
  if( (track[0]->DCA() < 10) && (track[1]->DCA() < 10) ) {
    if( (possibleJPsi.Pt() > 0.3) && (possibleJPsi.Pt() < 0.7) ){
      fInvariantMassDistributionStrictPtStrictDcaH->Fill( possibleJPsi.Mag() );
    }
    fInvariantMassDistributionVsPtStrictDcaH->Fill( possibleJPsi.Mag(), possibleJPsi.Pt() );
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
  }
  // fDimuonPtDistributionH->Fill(ptOfTheDimuonPair);
  if ( (possibleJPsi.Mag() > 2.8) && (possibleJPsi.Mag() < 3.3) ) {
    fDimuonPtDistributionH            ->Fill(ptOfTheDimuonPair);
    fDimuonPtDistributionShiftPlusOneH->Fill(ptOfTheDimuonPair);
  }

  if( ptOfTheDimuonPair < 0.200 ) {
    fInvariantMassDistributionCoherentShiftMinusTwoH->Fill(possibleJPsi.Mag());
  } else {
    fInvariantMassDistributionIncoherentShiftMinusTwoH->Fill(possibleJPsi.Mag());
  }
  if( ptOfTheDimuonPair < 0.225 ) {
    fInvariantMassDistributionCoherentShiftMinusOneH->Fill(possibleJPsi.Mag());
  } else {
    fInvariantMassDistributionIncoherentShiftMinusOneH->Fill(possibleJPsi.Mag());
  }
  if( ptOfTheDimuonPair < 0.275 ) {
    fInvariantMassDistributionCoherentShiftPlusOneH ->Fill(possibleJPsi.Mag());
  } else {
    fInvariantMassDistributionIncoherentShiftPlusOneH ->Fill(possibleJPsi.Mag());
  }
  if( ptOfTheDimuonPair < 0.300 ) {
    fInvariantMassDistributionCoherentShiftPlusTwoH ->Fill(possibleJPsi.Mag());
  } else {
    fInvariantMassDistributionIncoherentShiftPlusTwoH ->Fill(possibleJPsi.Mag());
  }




  //_____________________________________
  // STARTING THE CHECK FOR THE UPC GROUP
  /* - Here we will be using the information about the fMomentumAtDCA, to try
     - and see if there are any differences compared to before. Hopefully this
     - will be doing something. Maybe not. However there is already a major
     - difference with the way they are implemented. Here,
     - fMomentumAtDCA = [px, py, pz], while
     - fMomentum      = [pt, eta,phi]
     -
     -
     - RESULTS: the information here is too imprecise to be properly used.
     -          Switching to the usual fMomentum only. With the information
     -          at the DCA it is impossible to resolve the PsiPrime peak.
     -
   */
  TLorentzVector muonsAtDCA[2];
  TLorentzVector possibleJPsiAtDCA;
  for(int indexMuonAtDCA = 0; indexMuonAtDCA < 2; indexMuonAtDCA++) {
        /* - If dummy values then we cannot compute the J/Psi...
           - The dummy values are -999. but I do not trust the threshold for the
           - comparison... Better safe than sorry.
         */
        if( (track[indexMuonAtDCA]->PxAtDCA()) < -998. ) {
              PostData(1, fOutputList);
              return;
        }
        if( (track[indexMuonAtDCA]->PyAtDCA()) < -998. ) {
              PostData(1, fOutputList);
              return;
        }
        if( (track[indexMuonAtDCA]->PzAtDCA()) < -998. ) {
              PostData(1, fOutputList);
              return;
        }
        Double_t EnergyOfTheTrack = TMath::Sqrt(  track[indexMuonAtDCA]->PxAtDCA()*track[indexMuonAtDCA]->PxAtDCA() +
                                                  track[indexMuonAtDCA]->PyAtDCA()*track[indexMuonAtDCA]->PyAtDCA() +
                                                  track[indexMuonAtDCA]->PzAtDCA()*track[indexMuonAtDCA]->PzAtDCA() +
                                                  TDatabasePDG::Instance()->GetParticle(13)->Mass()*TDatabasePDG::Instance()->GetParticle(13)->Mass()/1000000
                                                  );
        muonsAtDCA[indexMuonAtDCA].SetPxPyPzE(   track[indexMuonAtDCA]->PxAtDCA(),
                                                 track[indexMuonAtDCA]->PyAtDCA(),
                                                 track[indexMuonAtDCA]->PzAtDCA(),
                                                 EnergyOfTheTrack
                                                 );
        possibleJPsiAtDCA += muonsAtDCA[indexMuonAtDCA];
  }
  fInvariantMassDistributionAtDcaH->Fill(possibleJPsiAtDCA.Mag());
  // END THE CHECK
  //_____________________________________




  /* - Now this is a critical part of  the analysis. What happens next is a
     - differential analysis in terms of the energy perceived by the neutron ZDC.
     - What it means is that now we may cut on those sensible values to select
     - only those J/Psi candidates falling under a certain peak of the neutron
     - ZNC energy distribution. It will be seen that the fZNCEnergyAgainstEntriesH
     - plot will present many gaussian-like peaks. Each peak represent an
     - increasingly large number of neutrons seen by the ZNC.
     -
     - Starting from the first peak, 0n, then 1n, hopefully 2n, but anything
     - else is more like a guess. If my understanding is good enough, even the
     - 2n peak requires user input to facilitate the minimizer's job.
     -
     - So, first thing first, Guillermo Contreras has suggested the preliminary
     - cut on the ZDC time, quoting:
     - "The energy value makes sense only when the time information is not
     - -999... You can choose times |t|<5 ns to plot the energy distributions
     - in the neutron ZDC".
     -
     - This happens with the request |fZNCTime|<5 if I have understood correctly.
     - After this we can fill whatever histogram we want to.
     -
     -
     -
     - NEW: after UPC meeting 5/3/2019
     - On ZDC timing. Usually we use time information from TDCs corresponding to
     - the common PMT (reads all four ZN sectors) on both sides. Each AOD event
     - contains information on up to four consecutive timing hits from these
     - TDCs within +/-12 bcs around the trigger bunch crossing. These hits are
     - stored in fZNATDCm and fZNCTDCm arrays:
     - https://github.com/alisw/AliRoot/blob/master/STEER/AOD/AliAODZDC.h#L153
     - and can be accessed as in:
     -
     - AliAODZDC* aodZDC = aod->GetZDCData();
     - for (Int_t i=0;i<4;i++) fZNATDC[i] = aodZDC->GetZNATDCm(i);
     - for (Int_t i=0;i<4;i++) fZNCTDC[i] = aodZDC->GetZNCTDCm(i);
     -
     - These hits may come from hadronic or EMD processes in neighbouring bcs.
     - In Pb-Pb we usually have 0-2 hits within +/-12 bcs mainly due to EMD.
     - Unused timing slots in these arrays are filled with large negative value
     - (-999). In order to check if there was a timing hit in the trigger bc,
     - you have to check if at least one timing hit out of four is within +/-2
     - ns around 0.
     -
     - Regarding these getters GetZNATime() and GetZNCTime(), defined here:
     - https://github.com/alisw/AliRoot/blob/master/STEER/AOD/AliAODZDC.h#L51
     - They are outdated because, as mentioned here, they return timing
     - information from the first slot in those arrays (fZNATDCm[0], fZNCTDCm[0]):
     - https://github.com/alisw/AliRoot/blob/master/STEER/AOD/AliAODZDC.h#L145
     - The first hit often corresponds to previous bunch crossings (e.g. EMD),
     - while interesting hit around 0 may be stored in the next slots.
     -
     -
   */
  Bool_t isZNAfired = kFALSE;
  Bool_t isZNCfired = kFALSE;
  Bool_t isZNAfiredStrict = kFALSE;
  Bool_t isZNCfiredStrict = kFALSE;
  Int_t  counterZNA = 0;
  Int_t  counterZNC = 0;
  /* - Note that in C++ the && and || operators "short-circuit". That means that
     - they only evaluate a parameter if required. If the first parameter to &&
     - is false, or the first to || is true, the rest will not be evaluated.
     - That means that writing:
     - if ( (isZNAfired == 0) && (...) )
     - should mean effectively
     - if ( isZNAfired != 0 ) continue;
     - hence it should be *at least* one hit!!!
     -
   */
  for(Int_t iZDC = 0; iZDC < 4 ; iZDC++) {
    if ( (isZNAfired == 0) && (fZNATDC[iZDC] > -2.) && (fZNATDC[iZDC] < 2.) ) {
      isZNAfired = kTRUE;
      /* - After mail with Chiara Oppedisano, it seems like the best way
         - to proceed is to firstly call the IsZNAfired() and then filling...
         -
         - If this doesn't appear in later pulls it is because this
         - doesn't seem to suit my case...
         -
       */
      if( dataZDC->IsZNAfired() ) {
        if ( (possibleJPsi.Mag() > 2.85) && (possibleJPsi.Mag() < 3.35) ){
          fZNATimeAgainstEntriesH->Fill(fZNATDC[iZDC]);
        }
      }
      fCounterZNAH->Fill(counterZNA);
    }
    if ( (isZNCfired == 0) && (fZNCTDC[iZDC] > -2.) && (fZNCTDC[iZDC] < 2.) ) {
      isZNCfired = kTRUE;
      if( dataZDC->IsZNCfired() ) {
        if ( (possibleJPsi.Mag() > 2.85) && (possibleJPsi.Mag() < 3.35) ){
          fZNCTimeAgainstEntriesH->Fill(fZNCTDC[iZDC]);
        }
      }
      fCounterZNCH->Fill(counterZNC);
    }
    counterZNA++;
    counterZNC++;
  }

  if ( isZNCfired != 0 ) {
    fZNCEnergyAgainstEntriesH        ->Fill(fZNCEnergy);
    fZNCEnergyAgainstEntriesExtendedH->Fill(fZNCEnergy);
    if ( calibrated == 0 ) fZNCEnergyUncalibratedH->Fill(fZNCEnergy);
    if ( calibrated == 1 ) {
      fZNCEnergyCalibratedH          ->Fill( fZNCEnergy );
      fZNCEnergyCalibratedHigherGainH->Fill( dataZDC->GetZNCTowerEnergyLR()[0] );
    }
  }
  fZNCEnergyBeforeTimingSelectionH        ->Fill(fZNCEnergy);
  fZNCEnergyBeforeTimingSelectionExtendedH->Fill(fZNCEnergy);
  if ( dataZDC->IsZNCfired() && ( isZNCfired != 0 ) ) {
    fZNCEnergyAgainstEntriesExtendedHv2->Fill(fZNCEnergy);
  }
  if ( dataZDC->IsZNAfired() && ( isZNAfired != 0 ) ) {
    fZNAEnergyAgainstEntriesExtendedHv2->Fill(fZNAEnergy);
  }
  if ( isZNAfired != 0 ) {
    fZNAEnergyAgainstEntriesH        ->Fill(fZNAEnergy);
    fZNAEnergyAgainstEntriesExtendedH->Fill(fZNAEnergy);
    if ( calibrated == 0 ) fZNAEnergyUncalibratedH->Fill(fZNAEnergy);
    if ( calibrated == 1 ) {
      fZNAEnergyCalibratedH          ->Fill( fZNAEnergy );
      fZNAEnergyCalibratedHigherGainH->Fill( dataZDC->GetZNATowerEnergyLR()[0] );
    }
  }
  fZNAEnergyBeforeTimingSelectionH        ->Fill(fZNAEnergy);
  fZNAEnergyBeforeTimingSelectionExtendedH->Fill(fZNAEnergy);

  /* - CHECKS for the timing:
     - Stricter timing window AND without timing window at all!
     -
   */
  for(Int_t iZDC = 0; iZDC < 4 ; iZDC++) {
    if ( (isZNAfiredStrict == 0) && (fZNATDC[iZDC] > -1.) && (fZNATDC[iZDC] < 1.) ) {
      isZNAfiredStrict = kTRUE;
      if( dataZDC->IsZNAfired() ) fZNATimeStrictTimeWindowH->Fill(fZNATDC[iZDC]);
    }
    if ( (isZNCfiredStrict == 0) && (fZNCTDC[iZDC] > -1.) && (fZNCTDC[iZDC] < 1.) ) {
      isZNCfiredStrict = kTRUE;
      if( dataZDC->IsZNCfired() ) fZNCTimeStrictTimeWindowH->Fill(fZNCTDC[iZDC]);
    }
    fZNATimeWithoutTimingH[iZDC]             ->Fill(fZNATDC[iZDC]);
    fZNCTimeWithoutTimingH[iZDC]             ->Fill(fZNCTDC[iZDC]);
    fZNCTime4FillingH                        ->Fill(fZNCTDC[iZDC]);
    fZNATime4FillingH                        ->Fill(fZNATDC[iZDC]);
    fZNCminusZNAtimeVsZNCplusZNAtimeH[iZDC]  ->Fill(fZNCTDC[iZDC]-fZNATDC[iZDC], fZNCTDC[iZDC]+fZNATDC[iZDC]);
    fZNCminusZNAtimeVsZNCplusZNAtime4FillingH->Fill(fZNCTDC[iZDC]-fZNATDC[iZDC], fZNCTDC[iZDC]+fZNATDC[iZDC]);
  }

  /* - Filling the v2 histogram only if the
   * - ZNC or the ZNA have detected any activity at all...
   */
  if( isZNCfired == 0 ) {
        if( isZNAfired == 0 ) {
              if ( (possibleJPsi.Mag() > 2.85) && (possibleJPsi.Mag() < 3.35) ) {
                  fDimuonPtDistributionZNCzeroZNAzeroHv2            ->Fill(ptOfTheDimuonPair);
                  // fDimuonPtDistributionZNCzeroZNAzeroShiftPlusOneHv2->Fill(ptOfTheDimuonPair);
                  if (        possibleJPsi.Rapidity() > -4.0  && possibleJPsi.Rapidity() <= -3.50 ) {
                    fDimuonPtDistributionZNCzeroZNAzeroRapidityHv2[0]->Fill(ptOfTheDimuonPair);
                  } else if ( possibleJPsi.Rapidity() > -3.50 && possibleJPsi.Rapidity() <= -3.00 ) {
                    fDimuonPtDistributionZNCzeroZNAzeroRapidityHv2[1]->Fill(ptOfTheDimuonPair);
                  } else if ( possibleJPsi.Rapidity() > -3.00 && possibleJPsi.Rapidity() <= -2.50 ) {
                    fDimuonPtDistributionZNCzeroZNAzeroRapidityHv2[2]->Fill(ptOfTheDimuonPair);
                  }
              }
              if( ptOfTheDimuonPair < 0.25 ) {
                  fInvariantMassDistributionCoherentZNCzeroZNAzeroHv2->Fill(possibleJPsi.Mag());
                  if (        possibleJPsi.Rapidity() > -4.0  && possibleJPsi.Rapidity() <= -3.50 ) {
                    fInvariantMassDistributionCoherentZNCzeroZNAzeroRapidityHv2[0]->Fill(possibleJPsi.Mag());
                  } else if ( possibleJPsi.Rapidity() > -3.50 && possibleJPsi.Rapidity() <= -3.00 ) {
                    fInvariantMassDistributionCoherentZNCzeroZNAzeroRapidityHv2[1]->Fill(possibleJPsi.Mag());
                  } else if ( possibleJPsi.Rapidity() > -3.00 && possibleJPsi.Rapidity() <= -2.50 ) {
                    fInvariantMassDistributionCoherentZNCzeroZNAzeroRapidityHv2[2]->Fill(possibleJPsi.Mag());
                  }
              } else {
                  fInvariantMassDistributionIncoherentZNCzeroZNAzeroHv2->Fill(possibleJPsi.Mag());
                  if (        possibleJPsi.Rapidity() > -4.0  && possibleJPsi.Rapidity() <= -3.50 ) {
                    fInvariantMassDistributionIncoherentZNCzeroZNAzeroRapidityHv2[0]->Fill(possibleJPsi.Mag());
                  } else if ( possibleJPsi.Rapidity() > -3.50 && possibleJPsi.Rapidity() <= -3.00 ) {
                    fInvariantMassDistributionIncoherentZNCzeroZNAzeroRapidityHv2[1]->Fill(possibleJPsi.Mag());
                  } else if ( possibleJPsi.Rapidity() > -3.00 && possibleJPsi.Rapidity() <= -2.50 ) {
                    fInvariantMassDistributionIncoherentZNCzeroZNAzeroRapidityHv2[2]->Fill(possibleJPsi.Mag());
                  }
              }
              // if( ptOfTheDimuonPair < 0.200 ) {
              //     fInvariantMassDistributionCoherentZNCzeroZNAzeroShiftMinusTwoHv2->Fill(possibleJPsi.Mag());
              // } else {
              //     fInvariantMassDistributionIncoherentZNCzeroZNAzeroShiftMinusTwoHv2->Fill(possibleJPsi.Mag());
              // }
              // if( ptOfTheDimuonPair < 0.225 ) {
              //     fInvariantMassDistributionCoherentZNCzeroZNAzeroShiftMinusOneHv2->Fill(possibleJPsi.Mag());
              // } else {
              //     fInvariantMassDistributionIncoherentZNCzeroZNAzeroShiftMinusOneHv2->Fill(possibleJPsi.Mag());
              // }
              // if( ptOfTheDimuonPair < 0.275 ) {
              //     fInvariantMassDistributionCoherentZNCzeroZNAzeroShiftPlusOneHv2->Fill(possibleJPsi.Mag());
              // } else {
              //     fInvariantMassDistributionIncoherentZNCzeroZNAzeroShiftPlusOneHv2->Fill(possibleJPsi.Mag());
              // }
              // if( ptOfTheDimuonPair < 0.300 ) {
              //     fInvariantMassDistributionCoherentZNCzeroZNAzeroShiftPlusTwoHv2->Fill(possibleJPsi.Mag());
              // } else {
              //     fInvariantMassDistributionIncoherentZNCzeroZNAzeroShiftPlusTwoHv2->Fill(possibleJPsi.Mag());
              // }
        } else {
              if ( (possibleJPsi.Mag() > 2.85) && (possibleJPsi.Mag() < 3.35) ) {
                  fDimuonPtDistributionZNCzeroZNAanyHv2            ->Fill(ptOfTheDimuonPair);
                  // fDimuonPtDistributionZNCzeroZNAanyShiftPlusOneHv2->Fill(ptOfTheDimuonPair);
                  if (        possibleJPsi.Rapidity() > -4.0  && possibleJPsi.Rapidity() <= -3.50 ) {
                    fDimuonPtDistributionZNCzeroZNAanyRapidityHv2[0]->Fill(ptOfTheDimuonPair);
                  } else if ( possibleJPsi.Rapidity() > -3.50 && possibleJPsi.Rapidity() <= -3.00 ) {
                    fDimuonPtDistributionZNCzeroZNAanyRapidityHv2[1]->Fill(ptOfTheDimuonPair);
                  } else if ( possibleJPsi.Rapidity() > -3.00 && possibleJPsi.Rapidity() <= -2.50 ) {
                    fDimuonPtDistributionZNCzeroZNAanyRapidityHv2[2]->Fill(ptOfTheDimuonPair);
                  }
              }
              if( ptOfTheDimuonPair < 0.25 ) {
                  fInvariantMassDistributionCoherentZNCzeroZNAanyHv2->Fill(possibleJPsi.Mag());
                  // fDimuonPtDistributionCoherentZNCzeroZNAanyHv2->Fill(ptOfTheDimuonPair);
                  if (        possibleJPsi.Rapidity() > -4.0  && possibleJPsi.Rapidity() <= -3.50 ) {
                    fInvariantMassDistributionCoherentZNCzeroZNAanyRapidityHv2[0]->Fill(possibleJPsi.Mag());
                  } else if ( possibleJPsi.Rapidity() > -3.50 && possibleJPsi.Rapidity() <= -3.00 ) {
                    fInvariantMassDistributionCoherentZNCzeroZNAanyRapidityHv2[1]->Fill(possibleJPsi.Mag());
                  } else if ( possibleJPsi.Rapidity() > -3.00 && possibleJPsi.Rapidity() <= -2.50 ) {
                    fInvariantMassDistributionCoherentZNCzeroZNAanyRapidityHv2[2]->Fill(possibleJPsi.Mag());
                  }
              } else {
                  fInvariantMassDistributionIncoherentZNCzeroZNAanyHv2->Fill(possibleJPsi.Mag());
                  // fDimuonPtDistributionIncoherentZNCzeroZNAanyHv2->Fill(ptOfTheDimuonPair);
                  if (        possibleJPsi.Rapidity() > -4.0  && possibleJPsi.Rapidity() <= -3.50 ) {
                    fInvariantMassDistributionIncoherentZNCzeroZNAanyRapidityHv2[0]->Fill(possibleJPsi.Mag());
                  } else if ( possibleJPsi.Rapidity() > -3.50 && possibleJPsi.Rapidity() <= -3.00 ) {
                    fInvariantMassDistributionIncoherentZNCzeroZNAanyRapidityHv2[1]->Fill(possibleJPsi.Mag());
                  } else if ( possibleJPsi.Rapidity() > -3.00 && possibleJPsi.Rapidity() <= -2.50 ) {
                    fInvariantMassDistributionIncoherentZNCzeroZNAanyRapidityHv2[2]->Fill(possibleJPsi.Mag());
                  }
              }
              // if( ptOfTheDimuonPair < 0.200 ) {
              //     fInvariantMassDistributionCoherentZNCzeroZNAanyShiftMinusTwoHv2->Fill(possibleJPsi.Mag());
              // } else {
              //     fInvariantMassDistributionIncoherentZNCzeroZNAanyShiftMinusTwoHv2->Fill(possibleJPsi.Mag());
              // }
              // if( ptOfTheDimuonPair < 0.225 ) {
              //     fInvariantMassDistributionCoherentZNCzeroZNAanyShiftMinusOneHv2->Fill(possibleJPsi.Mag());
              // } else {
              //     fInvariantMassDistributionIncoherentZNCzeroZNAanyShiftMinusOneHv2->Fill(possibleJPsi.Mag());
              // }
              // if( ptOfTheDimuonPair < 0.275 ) {
              //     fInvariantMassDistributionCoherentZNCzeroZNAanyShiftPlusOneHv2->Fill(possibleJPsi.Mag());
              // } else {
              //     fInvariantMassDistributionIncoherentZNCzeroZNAanyShiftPlusOneHv2->Fill(possibleJPsi.Mag());
              // }
              // if( ptOfTheDimuonPair < 0.300 ) {
              //     fInvariantMassDistributionCoherentZNCzeroZNAanyShiftPlusTwoHv2->Fill(possibleJPsi.Mag());
              // } else {
              //     fInvariantMassDistributionIncoherentZNCzeroZNAanyShiftPlusTwoHv2->Fill(possibleJPsi.Mag());
              // }
        }
  } else {
        if( isZNAfired == 0 ) {
              if ( (possibleJPsi.Mag() > 2.85) && (possibleJPsi.Mag() < 3.35) ) {
                  fDimuonPtDistributionZNCanyZNAzeroHv2            ->Fill(ptOfTheDimuonPair);
                  // fDimuonPtDistributionZNCanyZNAzeroShiftPlusOneHv2->Fill(ptOfTheDimuonPair);
                  if (        possibleJPsi.Rapidity() > -4.0  && possibleJPsi.Rapidity() <= -3.50 ) {
                    fDimuonPtDistributionZNCanyZNAzeroRapidityHv2[0]->Fill(ptOfTheDimuonPair);
                  } else if ( possibleJPsi.Rapidity() > -3.50 && possibleJPsi.Rapidity() <= -3.00 ) {
                    fDimuonPtDistributionZNCanyZNAzeroRapidityHv2[1]->Fill(ptOfTheDimuonPair);
                  } else if ( possibleJPsi.Rapidity() > -3.00 && possibleJPsi.Rapidity() <= -2.50 ) {
                    fDimuonPtDistributionZNCanyZNAzeroRapidityHv2[2]->Fill(ptOfTheDimuonPair);
                  }
              }
              if( ptOfTheDimuonPair < 0.25 ) {
                  fInvariantMassDistributionCoherentZNCanyZNAzeroHv2->Fill(possibleJPsi.Mag());
                  // fDimuonPtDistributionCoherentZNCanyZNAzeroHv2->Fill(ptOfTheDimuonPair);
                  if (        possibleJPsi.Rapidity() > -4.0  && possibleJPsi.Rapidity() <= -3.50 ) {
                    fInvariantMassDistributionCoherentZNCanyZNAzeroRapidityHv2[0]->Fill(possibleJPsi.Mag());
                  } else if ( possibleJPsi.Rapidity() > -3.50 && possibleJPsi.Rapidity() <= -3.00 ) {
                    fInvariantMassDistributionCoherentZNCanyZNAzeroRapidityHv2[1]->Fill(possibleJPsi.Mag());
                  } else if ( possibleJPsi.Rapidity() > -3.00 && possibleJPsi.Rapidity() <= -2.50 ) {
                    fInvariantMassDistributionCoherentZNCanyZNAzeroRapidityHv2[2]->Fill(possibleJPsi.Mag());
                  }
              } else {
                  fInvariantMassDistributionIncoherentZNCanyZNAzeroHv2->Fill(possibleJPsi.Mag());
                  // fDimuonPtDistributionIncoherentZNCanyZNAzeroHv2->Fill(ptOfTheDimuonPair);
                  if (        possibleJPsi.Rapidity() > -4.0  && possibleJPsi.Rapidity() <= -3.50 ) {
                    fInvariantMassDistributionIncoherentZNCanyZNAzeroRapidityHv2[0]->Fill(possibleJPsi.Mag());
                  } else if ( possibleJPsi.Rapidity() > -3.50 && possibleJPsi.Rapidity() <= -3.00 ) {
                    fInvariantMassDistributionIncoherentZNCanyZNAzeroRapidityHv2[1]->Fill(possibleJPsi.Mag());
                  } else if ( possibleJPsi.Rapidity() > -3.00 && possibleJPsi.Rapidity() <= -2.50 ) {
                    fInvariantMassDistributionIncoherentZNCanyZNAzeroRapidityHv2[2]->Fill(possibleJPsi.Mag());
                  }
              }
              // if( ptOfTheDimuonPair < 0.200 ) {
              //     fInvariantMassDistributionCoherentZNCanyZNAzeroShiftMinusTwoHv2->Fill(possibleJPsi.Mag());
              // } else {
              //     fInvariantMassDistributionIncoherentZNCanyZNAzeroShiftMinusTwoHv2->Fill(possibleJPsi.Mag());
              // }
              // if( ptOfTheDimuonPair < 0.225 ) {
              //     fInvariantMassDistributionCoherentZNCanyZNAzeroShiftMinusOneHv2->Fill(possibleJPsi.Mag());
              // } else {
              //     fInvariantMassDistributionIncoherentZNCanyZNAzeroShiftMinusOneHv2->Fill(possibleJPsi.Mag());
              // }
              // if( ptOfTheDimuonPair < 0.275 ) {
              //     fInvariantMassDistributionCoherentZNCanyZNAzeroShiftPlusOneHv2->Fill(possibleJPsi.Mag());
              // } else {
              //     fInvariantMassDistributionIncoherentZNCanyZNAzeroShiftPlusOneHv2->Fill(possibleJPsi.Mag());
              // }
              // if( ptOfTheDimuonPair < 0.300 ) {
              //     fInvariantMassDistributionCoherentZNCanyZNAzeroShiftPlusTwoHv2->Fill(possibleJPsi.Mag());
              // } else {
              //     fInvariantMassDistributionIncoherentZNCanyZNAzeroShiftPlusTwoHv2->Fill(possibleJPsi.Mag());
              // }
        } else {
              if ( (possibleJPsi.Mag() > 2.85) && (possibleJPsi.Mag() < 3.35) ) {
                  fDimuonPtDistributionZNCanyZNAanyHv2            ->Fill(ptOfTheDimuonPair);
                  // fDimuonPtDistributionZNCanyZNAanyShiftPlusOneHv2->Fill(ptOfTheDimuonPair);
                  if (        possibleJPsi.Rapidity() > -4.0  && possibleJPsi.Rapidity() <= -3.50 ) {
                    fDimuonPtDistributionZNCanyZNAanyRapidityHv2[0]->Fill(ptOfTheDimuonPair);
                  } else if ( possibleJPsi.Rapidity() > -3.50 && possibleJPsi.Rapidity() <= -3.00 ) {
                    fDimuonPtDistributionZNCanyZNAanyRapidityHv2[1]->Fill(ptOfTheDimuonPair);
                  } else if ( possibleJPsi.Rapidity() > -3.00 && possibleJPsi.Rapidity() <= -2.50 ) {
                    fDimuonPtDistributionZNCanyZNAanyRapidityHv2[2]->Fill(ptOfTheDimuonPair);
                  }
              }
              if( ptOfTheDimuonPair < 0.25 ) {
                  fInvariantMassDistributionCoherentZNCanyZNAanyHv2->Fill(possibleJPsi.Mag());
                  // fDimuonPtDistributionCoherentZNCanyZNAanyHv2->Fill(ptOfTheDimuonPair);
                  if (        possibleJPsi.Rapidity() > -4.0  && possibleJPsi.Rapidity() <= -3.50 ) {
                    fInvariantMassDistributionCoherentZNCanyZNAanyRapidityHv2[0]->Fill(possibleJPsi.Mag());
                  } else if ( possibleJPsi.Rapidity() > -3.50 && possibleJPsi.Rapidity() <= -3.00 ) {
                    fInvariantMassDistributionCoherentZNCanyZNAanyRapidityHv2[1]->Fill(possibleJPsi.Mag());
                  } else if ( possibleJPsi.Rapidity() > -3.00 && possibleJPsi.Rapidity() <= -2.50 ) {
                    fInvariantMassDistributionCoherentZNCanyZNAanyRapidityHv2[2]->Fill(possibleJPsi.Mag());
                  }
              } else {
                  fInvariantMassDistributionIncoherentZNCanyZNAanyHv2->Fill(possibleJPsi.Mag());
                  // fDimuonPtDistributionIncoherentZNCanyZNAanyHv2->Fill(ptOfTheDimuonPair);
                  if (        possibleJPsi.Rapidity() > -4.0  && possibleJPsi.Rapidity() <= -3.50 ) {
                    fInvariantMassDistributionIncoherentZNCanyZNAanyRapidityHv2[0]->Fill(possibleJPsi.Mag());
                  } else if ( possibleJPsi.Rapidity() > -3.50 && possibleJPsi.Rapidity() <= -3.00 ) {
                    fInvariantMassDistributionIncoherentZNCanyZNAanyRapidityHv2[1]->Fill(possibleJPsi.Mag());
                  } else if ( possibleJPsi.Rapidity() > -3.00 && possibleJPsi.Rapidity() <= -2.50 ) {
                    fInvariantMassDistributionIncoherentZNCanyZNAanyRapidityHv2[2]->Fill(possibleJPsi.Mag());
                  }
              }
              // if( ptOfTheDimuonPair < 0.200 ) {
              //     fInvariantMassDistributionCoherentZNCanyZNAanyShiftMinusTwoHv2->Fill(possibleJPsi.Mag());
              // } else {
              //     fInvariantMassDistributionIncoherentZNCanyZNAanyShiftMinusTwoHv2->Fill(possibleJPsi.Mag());
              // }
              // if( ptOfTheDimuonPair < 0.225 ) {
              //     fInvariantMassDistributionCoherentZNCanyZNAanyShiftMinusOneHv2->Fill(possibleJPsi.Mag());
              // } else {
              //     fInvariantMassDistributionIncoherentZNCanyZNAanyShiftMinusOneHv2->Fill(possibleJPsi.Mag());
              // }
              // if( ptOfTheDimuonPair < 0.275 ) {
              //     fInvariantMassDistributionCoherentZNCanyZNAanyShiftPlusOneHv2->Fill(possibleJPsi.Mag());
              // } else {
              //     fInvariantMassDistributionIncoherentZNCanyZNAanyShiftPlusOneHv2->Fill(possibleJPsi.Mag());
              // }
              // if( ptOfTheDimuonPair < 0.300 ) {
              //     fInvariantMassDistributionCoherentZNCanyZNAanyShiftPlusTwoHv2->Fill(possibleJPsi.Mag());
              // } else {
              //     fInvariantMassDistributionIncoherentZNCanyZNAanyShiftPlusTwoHv2->Fill(possibleJPsi.Mag());
              // }
        }
  }


  /*
  fZNCTimeAgainstEntriesH->Fill(fZNCTime);
  if( fZNCTime > -5.0 ) {
    if( fZNCTime < 5.0 ) {
          At any levels, this means |fZNCTime| < 5.
          fZNCEnergyAgainstEntriesH->Fill(fZNCEnergy);
          fZNAEnergyAgainstEntriesH->Fill(fZNAEnergy);
          if ( calibrated == 0 ) fZNAEnergyUncalibratedH->Fill(fZNAEnergy);
          if ( calibrated == 1 ) fZNAEnergyCalibratedH  ->Fill(fZNAEnergy);
          if ( calibrated == 0 ) fZNCEnergyUncalibratedH->Fill(fZNCEnergy);
          if ( calibrated == 1 ) fZNCEnergyCalibratedH  ->Fill(fZNCEnergy);

             - Now this offers the oppurtunity to do differential mass studies.
             - This can be seen here. When we try to do everything while cutting
             - on the ZNC energy.
             -
             - I don't if by the time you will be reading these lines of mine
             - the ZDC calibration will be done or not. For now I am Implementing
             - the cut based on Evgeny Kryshen's plot. Then we will see in the
             - future.
             -
             - NB: this is wrong and outdated. See next cycle for new code!

          if( fZNCEnergy > -300 ) {
                    if( fZNCEnergy < 125 ) {
                            fInvariantMassDistributionNoNeutronsH->Fill(possibleJPsi.Mag());
                            if( ptOfTheDimuonPair < 0.25) {
                                  fInvariantMassDistributionCoherentNoNeutronsH->Fill(possibleJPsi.Mag());
                            } else {
                                  fInvariantMassDistributionIncoherentNoNeutronsH->Fill(possibleJPsi.Mag());
                            }
                    } else if( fZNCEnergy < 375 ) {
                            fInvariantMassDistributionOneNeutronH->Fill(possibleJPsi.Mag());
                            fInvariantMassDistributionAtLeastOneNeutronH->Fill(possibleJPsi.Mag());
                            if( ptOfTheDimuonPair < 0.25) {
                                  fInvariantMassDistributionCoherentOneNeutronH->Fill(possibleJPsi.Mag());
                                  fInvariantMassDistributionCoherentAtLeastOneNeutronH->Fill(possibleJPsi.Mag());
                            } else {
                                  fInvariantMassDistributionIncoherentOneNeutronH->Fill(possibleJPsi.Mag());
                                  fInvariantMassDistributionIncoherentAtLeastOneNeutronH->Fill(possibleJPsi.Mag());
                            }
                    } else  {
                            fInvariantMassDistributionAtLeastOneNeutronH->Fill(possibleJPsi.Mag());
                            if( ptOfTheDimuonPair < 0.25) {
                                  fInvariantMassDistributionCoherentAtLeastOneNeutronH->Fill(possibleJPsi.Mag());
                            } else {
                                  fInvariantMassDistributionIncoherentAtLeastOneNeutronH->Fill(possibleJPsi.Mag());
                            }
                    }
          }
    }
  }
  */

  //_____________________________________
  // DIFFERENTIAL ANALYSIS NEUTRON EMISSION
  /* - This if should be really wrong...
     - But now I think I can do the same implementing th ZNA timing information
     - too by simply requesting:
     -  (1)   isZNAfired == kTRUE;
     -  (2)   isZNCfired == kTRUE;
     -
   */
  // if( fZNCTime > -5.0 ) {
  //   if( fZNCTime < 5.0 ) {
  // if ( isZNAfired != 0 ) {
  //   if ( isZNCfired != 0 ) {
          /* At any levels, this means |fZNCTime| < 2. */
          if( fZNCEnergy > -5000 ) {
                      if( fZNCEnergy < 1250 ) {
                                  if( fZNAEnergy > -5000 ) {
                                              if( fZNAEnergy < 1250 ) {
                                                          fDimuonPtDistributionZNCzeroZNAzeroH            ->Fill(ptOfTheDimuonPair);
                                                          fDimuonPtDistributionZNCzeroZNAzeroShiftPlusOneH->Fill(ptOfTheDimuonPair);
                                                          if( ptOfTheDimuonPair < 0.25 ) {
                                                                    fInvariantMassDistributionCoherentZNCzeroZNAzeroH->Fill(possibleJPsi.Mag());
                                                                    fDimuonPtDistributionCoherentZNCzeroZNAzeroH->Fill(ptOfTheDimuonPair);
                                                          } else {
                                                                    fInvariantMassDistributionIncoherentZNCzeroZNAzeroH->Fill(possibleJPsi.Mag());
                                                                    fDimuonPtDistributionIncoherentZNCzeroZNAzeroH->Fill(ptOfTheDimuonPair);
                                                          }
                                                          if( ptOfTheDimuonPair < 0.200 ) {
                                                            fInvariantMassDistributionCoherentZNCzeroZNAzeroShiftMinusTwoH->Fill(possibleJPsi.Mag());
                                                          } else {
                                                            fInvariantMassDistributionIncoherentZNCzeroZNAzeroShiftMinusTwoH->Fill(possibleJPsi.Mag());
                                                          }
                                                          if( ptOfTheDimuonPair < 0.225 ) {
                                                            fInvariantMassDistributionCoherentZNCzeroZNAzeroShiftMinusOneH->Fill(possibleJPsi.Mag());
                                                          } else {
                                                            fInvariantMassDistributionIncoherentZNCzeroZNAzeroShiftMinusOneH->Fill(possibleJPsi.Mag());
                                                          }
                                                          if( ptOfTheDimuonPair < 0.275 ) {
                                                            fInvariantMassDistributionCoherentZNCzeroZNAzeroShiftPlusOneH->Fill(possibleJPsi.Mag());
                                                          } else {
                                                            fInvariantMassDistributionIncoherentZNCzeroZNAzeroShiftPlusOneH->Fill(possibleJPsi.Mag());
                                                          }
                                                          if( ptOfTheDimuonPair < 0.300 ) {
                                                            fInvariantMassDistributionCoherentZNCzeroZNAzeroShiftPlusTwoH->Fill(possibleJPsi.Mag());
                                                          } else {
                                                            fInvariantMassDistributionIncoherentZNCzeroZNAzeroShiftPlusTwoH->Fill(possibleJPsi.Mag());
                                                          }
                                              } else {
                                                          fDimuonPtDistributionZNCzeroZNAanyH            ->Fill(ptOfTheDimuonPair);
                                                          fDimuonPtDistributionZNCzeroZNAanyShiftPlusOneH->Fill(ptOfTheDimuonPair);
                                                          if( ptOfTheDimuonPair < 0.25 ) {
                                                                    fInvariantMassDistributionCoherentZNCzeroZNAanyH->Fill(possibleJPsi.Mag());
                                                                    fDimuonPtDistributionCoherentZNCzeroZNAanyH->Fill(ptOfTheDimuonPair);
                                                          } else {
                                                                    fInvariantMassDistributionIncoherentZNCzeroZNAanyH->Fill(possibleJPsi.Mag());
                                                                    fDimuonPtDistributionIncoherentZNCzeroZNAanyH->Fill(ptOfTheDimuonPair);
                                                          }
                                                          if( ptOfTheDimuonPair < 0.200 ) {
                                                            fInvariantMassDistributionCoherentZNCzeroZNAanyShiftMinusTwoH->Fill(possibleJPsi.Mag());
                                                          } else {
                                                            fInvariantMassDistributionIncoherentZNCzeroZNAanyShiftMinusTwoH->Fill(possibleJPsi.Mag());
                                                          }
                                                          if( ptOfTheDimuonPair < 0.225 ) {
                                                            fInvariantMassDistributionCoherentZNCzeroZNAanyShiftMinusOneH->Fill(possibleJPsi.Mag());
                                                          } else {
                                                            fInvariantMassDistributionIncoherentZNCzeroZNAanyShiftMinusOneH->Fill(possibleJPsi.Mag());
                                                          }
                                                          if( ptOfTheDimuonPair < 0.275 ) {
                                                            fInvariantMassDistributionCoherentZNCzeroZNAanyShiftPlusOneH->Fill(possibleJPsi.Mag());
                                                          } else {
                                                            fInvariantMassDistributionIncoherentZNCzeroZNAanyShiftPlusOneH->Fill(possibleJPsi.Mag());
                                                          }
                                                          if( ptOfTheDimuonPair < 0.300 ) {
                                                            fInvariantMassDistributionCoherentZNCzeroZNAanyShiftPlusTwoH->Fill(possibleJPsi.Mag());
                                                          } else {
                                                            fInvariantMassDistributionIncoherentZNCzeroZNAanyShiftPlusTwoH->Fill(possibleJPsi.Mag());
                                                          }
                                              }
                                  }
                      } else {
                                  if( fZNAEnergy > -5000 ) {
                                              if( fZNAEnergy < 1250 ) {
                                                          fDimuonPtDistributionZNCanyZNAzeroH            ->Fill(ptOfTheDimuonPair);
                                                          fDimuonPtDistributionZNCanyZNAzeroShiftPlusOneH->Fill(ptOfTheDimuonPair);
                                                          if( ptOfTheDimuonPair < 0.25 ) {
                                                                    fInvariantMassDistributionCoherentZNCanyZNAzeroH->Fill(possibleJPsi.Mag());
                                                                    fDimuonPtDistributionCoherentZNCanyZNAzeroH->Fill(ptOfTheDimuonPair);
                                                          } else {
                                                                    fInvariantMassDistributionIncoherentZNCanyZNAzeroH->Fill(possibleJPsi.Mag());
                                                                    fDimuonPtDistributionIncoherentZNCanyZNAzeroH->Fill(ptOfTheDimuonPair);
                                                          }
                                                          if( ptOfTheDimuonPair < 0.200 ) {
                                                            fInvariantMassDistributionCoherentZNCanyZNAzeroShiftMinusTwoH->Fill(possibleJPsi.Mag());
                                                          } else {
                                                            fInvariantMassDistributionIncoherentZNCanyZNAzeroShiftMinusTwoH->Fill(possibleJPsi.Mag());
                                                          }
                                                          if( ptOfTheDimuonPair < 0.225 ) {
                                                            fInvariantMassDistributionCoherentZNCanyZNAzeroShiftMinusOneH->Fill(possibleJPsi.Mag());
                                                          } else {
                                                            fInvariantMassDistributionIncoherentZNCanyZNAzeroShiftMinusOneH->Fill(possibleJPsi.Mag());
                                                          }
                                                          if( ptOfTheDimuonPair < 0.275 ) {
                                                            fInvariantMassDistributionCoherentZNCanyZNAzeroShiftPlusOneH->Fill(possibleJPsi.Mag());
                                                          } else {
                                                            fInvariantMassDistributionIncoherentZNCanyZNAzeroShiftPlusOneH->Fill(possibleJPsi.Mag());
                                                          }
                                                          if( ptOfTheDimuonPair < 0.300 ) {
                                                            fInvariantMassDistributionCoherentZNCanyZNAzeroShiftPlusTwoH->Fill(possibleJPsi.Mag());
                                                          } else {
                                                            fInvariantMassDistributionIncoherentZNCanyZNAzeroShiftPlusTwoH->Fill(possibleJPsi.Mag());
                                                          }
                                              } else {
                                                          fDimuonPtDistributionZNCanyZNAanyH            ->Fill(ptOfTheDimuonPair);
                                                          fDimuonPtDistributionZNCanyZNAanyShiftPlusOneH->Fill(ptOfTheDimuonPair);
                                                          if( ptOfTheDimuonPair < 0.25 ) {
                                                                    fInvariantMassDistributionCoherentZNCanyZNAanyH->Fill(possibleJPsi.Mag());
                                                                    fDimuonPtDistributionCoherentZNCanyZNAanyH->Fill(ptOfTheDimuonPair);
                                                          } else {
                                                                    fInvariantMassDistributionIncoherentZNCanyZNAanyH->Fill(possibleJPsi.Mag());
                                                                    fDimuonPtDistributionIncoherentZNCanyZNAanyH->Fill(ptOfTheDimuonPair);
                                                          }
                                                          if( ptOfTheDimuonPair < 0.200 ) {
                                                            fInvariantMassDistributionCoherentZNCanyZNAanyShiftMinusTwoH->Fill(possibleJPsi.Mag());
                                                          } else {
                                                            fInvariantMassDistributionIncoherentZNCanyZNAanyShiftMinusTwoH->Fill(possibleJPsi.Mag());
                                                          }
                                                          if( ptOfTheDimuonPair < 0.225 ) {
                                                            fInvariantMassDistributionCoherentZNCanyZNAanyShiftMinusOneH->Fill(possibleJPsi.Mag());
                                                          } else {
                                                            fInvariantMassDistributionIncoherentZNCanyZNAanyShiftMinusOneH->Fill(possibleJPsi.Mag());
                                                          }
                                                          if( ptOfTheDimuonPair < 0.275 ) {
                                                            fInvariantMassDistributionCoherentZNCanyZNAanyShiftPlusOneH->Fill(possibleJPsi.Mag());
                                                          } else {
                                                            fInvariantMassDistributionIncoherentZNCanyZNAanyShiftPlusOneH->Fill(possibleJPsi.Mag());
                                                          }
                                                          if( ptOfTheDimuonPair < 0.300 ) {
                                                            fInvariantMassDistributionCoherentZNCanyZNAanyShiftPlusTwoH->Fill(possibleJPsi.Mag());
                                                          } else {
                                                            fInvariantMassDistributionIncoherentZNCanyZNAanyShiftPlusTwoH->Fill(possibleJPsi.Mag());
                                                          }
                                              }
                                  }

                      }
          }
  //   }
  // }


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
  muonsCopy2[0] = muonsCopy[0];
  muonsCopy2[1] = muonsCopy[1];
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
  Double_t possibleJPsiCopyMag =possibleJPsiCopy.Mag();


  /* - NEW:
     -
   */
  Bool_t controlFlag = 0;
  if ( possibleJPsiCopy.Pt() < 0.25 ) {
        Double_t CosThetaHelicityFrameValue = CosThetaHelicityFrame( muonsCopy2[0], muonsCopy2[1], possibleJPsiCopy );
        Double_t PhiHelicityFrameValue      =   CosPhiHelicityFrame( muonsCopy2[0], muonsCopy2[1], possibleJPsiCopy );
        for(Int_t iCosThetaBins = 0; iCosThetaBins < 10; iCosThetaBins++) {
          if( controlFlag == 1) break;
          if( (CosThetaHelicityFrameValue + 1.) < 2.*((Double_t)iCosThetaBins + 1.)/10. ){
            for(Int_t iPhiBins = 0; iPhiBins < 10; iPhiBins++) {
              if( controlFlag == 1) break;
              if( (PhiHelicityFrameValue + 3.14) < 6.28*((Double_t)iPhiBins + 1.)/10. ){
                  fInvariantMassDistributionForSignalExtractionHelicityFrameH[iCosThetaBins][iPhiBins]->Fill(possibleJPsiCopyMag);
                  controlFlag = 1;
              }
            }
          }
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
              fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameH[iCosThetaBins]->Fill(possibleJPsiCopyMag);
              controlFlag2 = 1;
          }
        }
        for(Int_t iPhiBins = 0; iPhiBins < 50; iPhiBins++) {
          if( controlFlag3 == 1) break;
          if( (PhiHelicityFrameValue + 3.14) < 6.28*((Double_t)iPhiBins + 1.)/50. ){
              fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameH[iPhiBins]->Fill(possibleJPsiCopyMag);
              controlFlag3 = 1;
          }
        }

  }

  /* - NEW:
     -
   */
  if ( possibleJPsiCopy.Pt() < 0.25 ) {
        Double_t CosThetaHelicityFrameValue2 = CosThetaHelicityFrame( muonsCopy2[0], muonsCopy2[1], possibleJPsiCopy );
        Double_t PhiHelicityFrameValue2      =   CosPhiHelicityFrame( muonsCopy2[0], muonsCopy2[1], possibleJPsiCopy );
        if(      CosThetaHelicityFrameValue2 < -0.95 ){fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameTryingH[0]->Fill(possibleJPsiCopyMag);}
        else if( CosThetaHelicityFrameValue2 < -0.90 ){fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameTryingH[1]->Fill(possibleJPsiCopyMag);}
        else if( CosThetaHelicityFrameValue2 < -0.85 ){fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameTryingH[2]->Fill(possibleJPsiCopyMag);}
        else if( CosThetaHelicityFrameValue2 < -0.80 ){fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameTryingH[3]->Fill(possibleJPsiCopyMag);}
        else if( CosThetaHelicityFrameValue2 < -0.75 ){fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameTryingH[4]->Fill(possibleJPsiCopyMag);}
        else if( CosThetaHelicityFrameValue2 < -0.70 ){fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameTryingH[5]->Fill(possibleJPsiCopyMag);}
        else if( CosThetaHelicityFrameValue2 < -0.65 ){fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameTryingH[6]->Fill(possibleJPsiCopyMag);}
        else if( CosThetaHelicityFrameValue2 < -0.60 ){fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameTryingH[7]->Fill(possibleJPsiCopyMag);}
        else if( CosThetaHelicityFrameValue2 < -0.55 ){fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameTryingH[8]->Fill(possibleJPsiCopyMag);}
        else if( CosThetaHelicityFrameValue2 < -0.50 ){fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameTryingH[9]->Fill(possibleJPsiCopyMag);}
        else if( CosThetaHelicityFrameValue2 < -0.45 ){fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameTryingH[10]->Fill(possibleJPsiCopyMag);}
        else if( CosThetaHelicityFrameValue2 < -0.40 ){fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameTryingH[11]->Fill(possibleJPsiCopyMag);}
        else if( CosThetaHelicityFrameValue2 < -0.35 ){fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameTryingH[12]->Fill(possibleJPsiCopyMag);}
        else if( CosThetaHelicityFrameValue2 < -0.30 ){fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameTryingH[13]->Fill(possibleJPsiCopyMag);}
        else if( CosThetaHelicityFrameValue2 < -0.25 ){fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameTryingH[14]->Fill(possibleJPsiCopyMag);}
        else if( CosThetaHelicityFrameValue2 < -0.20 ){fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameTryingH[15]->Fill(possibleJPsiCopyMag);}
        else if( CosThetaHelicityFrameValue2 < -0.15 ){fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameTryingH[16]->Fill(possibleJPsiCopyMag);}
        else if( CosThetaHelicityFrameValue2 < -0.10 ){fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameTryingH[17]->Fill(possibleJPsiCopyMag);}
        else if( CosThetaHelicityFrameValue2 < -0.05 ){fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameTryingH[18]->Fill(possibleJPsiCopyMag);}
        else if( CosThetaHelicityFrameValue2 <  0.00 ){fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameTryingH[19]->Fill(possibleJPsiCopyMag);}
        else if( CosThetaHelicityFrameValue2 < +0.05 ){fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameTryingH[20]->Fill(possibleJPsiCopyMag);}
        else if( CosThetaHelicityFrameValue2 < +0.10 ){fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameTryingH[21]->Fill(possibleJPsiCopyMag);}
        else if( CosThetaHelicityFrameValue2 < +0.15 ){fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameTryingH[22]->Fill(possibleJPsiCopyMag);}
        else if( CosThetaHelicityFrameValue2 < +0.20 ){fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameTryingH[23]->Fill(possibleJPsiCopyMag);}
        else if( CosThetaHelicityFrameValue2 < +0.25 ){fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameTryingH[24]->Fill(possibleJPsiCopyMag);}
        else if( CosThetaHelicityFrameValue2 < +0.30 ){fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameTryingH[25]->Fill(possibleJPsiCopyMag);}
        else if( CosThetaHelicityFrameValue2 < +0.35 ){fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameTryingH[26]->Fill(possibleJPsiCopyMag);}
        else if( CosThetaHelicityFrameValue2 < +0.40 ){fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameTryingH[27]->Fill(possibleJPsiCopyMag);}
        else if( CosThetaHelicityFrameValue2 < +0.45 ){fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameTryingH[28]->Fill(possibleJPsiCopyMag);}
        else if( CosThetaHelicityFrameValue2 < +0.50 ){fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameTryingH[29]->Fill(possibleJPsiCopyMag);}
        else if( CosThetaHelicityFrameValue2 < +0.55 ){fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameTryingH[30]->Fill(possibleJPsiCopyMag);}
        else if( CosThetaHelicityFrameValue2 < +0.60 ){fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameTryingH[31]->Fill(possibleJPsiCopyMag);}
        else if( CosThetaHelicityFrameValue2 < +0.65 ){fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameTryingH[32]->Fill(possibleJPsiCopyMag);}
        else if( CosThetaHelicityFrameValue2 < +0.70 ){fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameTryingH[33]->Fill(possibleJPsiCopyMag);}
        else if( CosThetaHelicityFrameValue2 < +0.75 ){fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameTryingH[34]->Fill(possibleJPsiCopyMag);}
        else if( CosThetaHelicityFrameValue2 < +0.80 ){fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameTryingH[35]->Fill(possibleJPsiCopyMag);}
        else if( CosThetaHelicityFrameValue2 < +0.85 ){fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameTryingH[36]->Fill(possibleJPsiCopyMag);}
        else if( CosThetaHelicityFrameValue2 < +0.90 ){fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameTryingH[37]->Fill(possibleJPsiCopyMag);}
        else if( CosThetaHelicityFrameValue2 < +0.95 ){fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameTryingH[38]->Fill(possibleJPsiCopyMag);}
        else if( CosThetaHelicityFrameValue2 < +1.00 ){fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameTryingH[39]->Fill(possibleJPsiCopyMag);}

        if( PhiHelicityFrameValue2 > -3.14 && PhiHelicityFrameValue2 < -3.0144 ){fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameTryingH[0]->Fill(possibleJPsiCopyMag);}
        else if(                              PhiHelicityFrameValue2 < -2.8888 ){fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameTryingH[1]->Fill(possibleJPsiCopyMag);}
        else if(                              PhiHelicityFrameValue2 < -2.7632 ){fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameTryingH[2]->Fill(possibleJPsiCopyMag);}
        else if(                              PhiHelicityFrameValue2 < -2.6376 ){fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameTryingH[3]->Fill(possibleJPsiCopyMag);}
        else if(                              PhiHelicityFrameValue2 < -2.5120 ){fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameTryingH[4]->Fill(possibleJPsiCopyMag);}
        else if(                              PhiHelicityFrameValue2 < -2.3864 ){fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameTryingH[5]->Fill(possibleJPsiCopyMag);}
        else if(                              PhiHelicityFrameValue2 < -2.2608 ){fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameTryingH[6]->Fill(possibleJPsiCopyMag);}
        else if(                              PhiHelicityFrameValue2 < -2.1352 ){fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameTryingH[7]->Fill(possibleJPsiCopyMag);}
        else if(                              PhiHelicityFrameValue2 < -2.0096 ){fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameTryingH[8]->Fill(possibleJPsiCopyMag);}
        else if(                              PhiHelicityFrameValue2 < -1.8840 ){fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameTryingH[9]->Fill(possibleJPsiCopyMag);}
        else if(                              PhiHelicityFrameValue2 < -1.7584 ){fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameTryingH[10]->Fill(possibleJPsiCopyMag);}
        else if(                              PhiHelicityFrameValue2 < -1.6328 ){fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameTryingH[11]->Fill(possibleJPsiCopyMag);}
        else if(                              PhiHelicityFrameValue2 < -1.5072 ){fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameTryingH[12]->Fill(possibleJPsiCopyMag);}
        else if(                              PhiHelicityFrameValue2 < -1.3816 ){fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameTryingH[13]->Fill(possibleJPsiCopyMag);}
        else if(                              PhiHelicityFrameValue2 < -1.2560 ){fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameTryingH[14]->Fill(possibleJPsiCopyMag);}
        else if(                              PhiHelicityFrameValue2 < -1.1304 ){fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameTryingH[15]->Fill(possibleJPsiCopyMag);}
        else if(                              PhiHelicityFrameValue2 < -1.0048 ){fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameTryingH[16]->Fill(possibleJPsiCopyMag);}
        else if(                              PhiHelicityFrameValue2 < -0.8792 ){fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameTryingH[17]->Fill(possibleJPsiCopyMag);}
        else if(                              PhiHelicityFrameValue2 < -0.7536 ){fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameTryingH[18]->Fill(possibleJPsiCopyMag);}
        else if(                              PhiHelicityFrameValue2 < -0.6280 ){fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameTryingH[19]->Fill(possibleJPsiCopyMag);}
        else if(                              PhiHelicityFrameValue2 < -0.5024 ){fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameTryingH[20]->Fill(possibleJPsiCopyMag);}
        else if(                              PhiHelicityFrameValue2 < -0.3768 ){fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameTryingH[21]->Fill(possibleJPsiCopyMag);}
        else if(                              PhiHelicityFrameValue2 < -0.2512 ){fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameTryingH[22]->Fill(possibleJPsiCopyMag);}
        else if(                              PhiHelicityFrameValue2 < -0.1256 ){fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameTryingH[23]->Fill(possibleJPsiCopyMag);}
        else if(                              PhiHelicityFrameValue2 <  0.0000 ){fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameTryingH[24]->Fill(possibleJPsiCopyMag);}
        else if(                              PhiHelicityFrameValue2 < +0.1256 ){fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameTryingH[25]->Fill(possibleJPsiCopyMag);}
        else if(                              PhiHelicityFrameValue2 < +0.2512 ){fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameTryingH[26]->Fill(possibleJPsiCopyMag);}
        else if(                              PhiHelicityFrameValue2 < +0.3768 ){fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameTryingH[27]->Fill(possibleJPsiCopyMag);}
        else if(                              PhiHelicityFrameValue2 < +0.5024 ){fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameTryingH[28]->Fill(possibleJPsiCopyMag);}
        else if(                              PhiHelicityFrameValue2 < +0.6280 ){fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameTryingH[29]->Fill(possibleJPsiCopyMag);}
        else if(                              PhiHelicityFrameValue2 < +0.7536 ){fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameTryingH[30]->Fill(possibleJPsiCopyMag);}
        else if(                              PhiHelicityFrameValue2 < +0.8792 ){fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameTryingH[31]->Fill(possibleJPsiCopyMag);}
        else if(                              PhiHelicityFrameValue2 < +1.0048 ){fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameTryingH[32]->Fill(possibleJPsiCopyMag);}
        else if(                              PhiHelicityFrameValue2 < +1.1304 ){fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameTryingH[33]->Fill(possibleJPsiCopyMag);}
        else if(                              PhiHelicityFrameValue2 < +1.2560 ){fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameTryingH[34]->Fill(possibleJPsiCopyMag);}
        else if(                              PhiHelicityFrameValue2 < +1.3816 ){fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameTryingH[35]->Fill(possibleJPsiCopyMag);}
        else if(                              PhiHelicityFrameValue2 < +1.5072 ){fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameTryingH[36]->Fill(possibleJPsiCopyMag);}
        else if(                              PhiHelicityFrameValue2 < +1.6328 ){fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameTryingH[37]->Fill(possibleJPsiCopyMag);}
        else if(                              PhiHelicityFrameValue2 < +1.7584 ){fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameTryingH[38]->Fill(possibleJPsiCopyMag);}
        else if(                              PhiHelicityFrameValue2 < +1.8840 ){fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameTryingH[39]->Fill(possibleJPsiCopyMag);}
        else if(                              PhiHelicityFrameValue2 < +2.0096 ){fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameTryingH[40]->Fill(possibleJPsiCopyMag);}
        else if(                              PhiHelicityFrameValue2 < +2.1352 ){fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameTryingH[41]->Fill(possibleJPsiCopyMag);}
        else if(                              PhiHelicityFrameValue2 < +2.2608 ){fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameTryingH[42]->Fill(possibleJPsiCopyMag);}
        else if(                              PhiHelicityFrameValue2 < +2.3864 ){fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameTryingH[43]->Fill(possibleJPsiCopyMag);}
        else if(                              PhiHelicityFrameValue2 < +2.5120 ){fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameTryingH[44]->Fill(possibleJPsiCopyMag);}
        else if(                              PhiHelicityFrameValue2 < +2.6376 ){fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameTryingH[45]->Fill(possibleJPsiCopyMag);}
        else if(                              PhiHelicityFrameValue2 < +2.7632 ){fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameTryingH[46]->Fill(possibleJPsiCopyMag);}
        else if(                              PhiHelicityFrameValue2 < +2.8888 ){fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameTryingH[47]->Fill(possibleJPsiCopyMag);}
        else if(                              PhiHelicityFrameValue2 < +3.0144 ){fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameTryingH[48]->Fill(possibleJPsiCopyMag);}
        else if(                              PhiHelicityFrameValue2 < +3.14   ){fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameTryingH[49]->Fill(possibleJPsiCopyMag);}

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
                  fInvariantMassDistributionForSignalExtractionHelicityFrameMyBinningH[iCosThetaBins][iPhiBins]->Fill(possibleJPsiCopyMag);
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
                  fInvariantMassDistributionForSignalExtractionCsFrameMyBinningH[iCosThetaBins][iPhiBins]->Fill(possibleJPsiCopyMag);
                  controlFlag4_CS = 1;
              }
            }
          }
        }
  }

  /* - NEW:
   * - My variable binning for 1D analysis
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
        for(Int_t iCosThetaBins = 0; iCosThetaBins < 26; iCosThetaBins++) {
          if( controlFlag5 == 1) break;
          if( CosThetaHelicityFrameValue4 < MyVariableCosThetaBinning1D[iCosThetaBins + 1] ){
            fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameMyVariableBinningH[iCosThetaBins]->Fill(possibleJPsiCopyMag);
            controlFlag5 = 1;
          }
        }
        for(Int_t iPhiBins = 0; iPhiBins < 30; iPhiBins++) {
          if( controlFlag6 == 1) break;
          if( PhiHelicityFrameValue4  < MyVariablePhiBinning1D[iPhiBins + 1] ){
            fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameMyVariableBinningH[iPhiBins]->Fill(possibleJPsiCopyMag);
            controlFlag6 = 1;
          }
        }
  }

  /* - NEW:
   * - My STRICT variable binning for 1D analysis
   */
  Bool_t controlFlag7 = 0;
  Bool_t controlFlag8 = 0;
  Double_t MyVariableCosThetaBinning1Dv2[] = { -0.65, -0.45, -0.3, -0.15,  -0.05,
                                                0.05,  0.15,  0.3,  0.45,   0.65 };
  Double_t MyVariablePhiBinning1Dv2[] = { -3.14*1,    -3.14*13/15,  -3.14*11/15,     -3.14*9/15,
                                          -3.14*7/15, -3.14*5/15,   -3.14*3/15,      -3.14*1/15,
                                          +3.14*1/15, +3.14*3/15,   +3.14*5/15,      +3.14*7/15,
                                          +3.14*9/15, +3.14*11/15,  +3.14*13/15,     +3.14*1 };
  if ( possibleJPsiCopy.Pt() < 0.25 ) {
        Double_t CosThetaHelicityFrameValue5 = CosThetaHelicityFrame( muonsCopy2[0], muonsCopy2[1], possibleJPsiCopy );
        Double_t PhiHelicityFrameValue5      =   CosPhiHelicityFrame( muonsCopy2[0], muonsCopy2[1], possibleJPsiCopy );
        for(Int_t iCosThetaBins = 0; iCosThetaBins < 9; iCosThetaBins++) {
          if( controlFlag7 == 1) break;
          if( CosThetaHelicityFrameValue5 < MyVariableCosThetaBinning1Dv2[iCosThetaBins + 1] ){
            fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameMyStrictVariableBinningH[iCosThetaBins]->Fill(possibleJPsiCopyMag);
            controlFlag7 = 1;
          }
        }
        for(Int_t iPhiBins = 0; iPhiBins < 15; iPhiBins++) {
          if( controlFlag8 == 1) break;
          if( PhiHelicityFrameValue5  < MyVariablePhiBinning1Dv2[iPhiBins + 1] ){
            fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameMyStrictVariableBinningH[iPhiBins]->Fill(possibleJPsiCopyMag);
            controlFlag8 = 1;
          }
        }
  }

  /* - NEW:
   * - My SEVENTEEN BINS variable binning for 1D analysis
   */
  Bool_t controlFlag9  = 0;
  Bool_t controlFlag10 = 0;
  Bool_t controlFlag11 = 0;
  Bool_t controlFlag12 = 0;
  Double_t MyVariableCosThetaBinning1Dv3[] = { -0.65,  -0.5,  -0.425, -0.35,
                                               -0.275, -0.2,  -0.125, -0.075,
                                               -0.025,  0.025, 0.075,  0.125,
                                                0.2,    0.275, 0.35,   0.425,
                                                0.5,    0.65 };
  if ( possibleJPsiCopy.Pt() < 0.25 ) {
        Double_t CosThetaHelicityFrameValue6 = CosThetaHelicityFrame( muonsCopy2[0], muonsCopy2[1], possibleJPsiCopy );
        Double_t PhiHelicityFrameValue6      =   CosPhiHelicityFrame( muonsCopy2[0], muonsCopy2[1], possibleJPsiCopy );
        Double_t WeightFromAcceptance        =      AccEffCorrection( CosThetaHelicityFrameValue6 );
        if ( ( possibleJPsiCopyMag > 2.85 ) && ( possibleJPsiCopyMag < 3.35 ) ) {
          fCosThetaHelicityFrameJPsiAlreadyCorrectedH->Fill( CosThetaHelicityFrameValue6, 1/WeightFromAcceptance );
        }
        for(Int_t iCosThetaBins = 0; iCosThetaBins < 17; iCosThetaBins++) {
          if( controlFlag9 == 1) break;
          if( CosThetaHelicityFrameValue6 < MyVariableCosThetaBinning1Dv3[iCosThetaBins + 1] ){
            fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameMySeventeenBinsVariableBinningH[iCosThetaBins]->Fill(possibleJPsiCopyMag);
            controlFlag9 = 1;
          }
        }
        for(Int_t iCosThetaBins = 0; iCosThetaBins < 50; iCosThetaBins++) {
          if( controlFlag11 == 1) break;
          if( (CosThetaHelicityFrameValue6 + 1.) < 2.*((Double_t)iCosThetaBins + 1.)/50. ) {
            fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameForAlreadyCorrectedFiftyH[iCosThetaBins]->Fill(possibleJPsiCopyMag);
            controlFlag11 = 1;
          }
        }
        for(Int_t iCosThetaBins = 0; iCosThetaBins < 100; iCosThetaBins++) {
          if( controlFlag12 == 1) break;
          if( (CosThetaHelicityFrameValue6 + 1.) < 2.*((Double_t)iCosThetaBins + 1.)/100. ) {
            fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameForAlreadyCorrectedHundredH[iCosThetaBins]->Fill(possibleJPsiCopyMag);
            controlFlag12 = 1;
          }
        }
  }

  /* - NEW:
     -
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
            fInvariantMassDistributionOnlyCosThetaHeFrameTwentyfiveBinsH[iCosThetaBins]->Fill(possibleJPsiCopyMag);
            controlFlag13 = 1;
          }
        }
        for(Int_t iPhiBins = 0; iPhiBins < 25; iPhiBins++) {
          if( controlFlag14 == 1) break;
          if( (PhiHelicityFrameValue7 + 3.14) < 6.28*((Double_t)iPhiBins + 1.)/25. ){
            fInvariantMassDistributionOnlyPhiHeFrameTwentyfiveBinsH[iPhiBins]->Fill(possibleJPsiCopyMag);
            controlFlag14 = 1;
          }
        }
        if( (CosThetaHelicityFrameValue7 > 0) || (CosThetaHelicityFrameValue7 == 0) ){
          for(Int_t iTildePhiBins = 0; iTildePhiBins < 25; iTildePhiBins++) {
            if( controlFlag15 == 1) break;
            if( (TildePhiPositiveCosTheta) < 6.28*((Double_t)iTildePhiBins + 1.)/25. ){
              fInvariantMassDistributionOnlyTildePhiHeFrameTwentyfiveBinsH[iTildePhiBins]->Fill(possibleJPsiCopyMag);
              controlFlag15 = 1;
            }
          }
        } else if ( CosThetaHelicityFrameValue7 < 0 ){
          for(Int_t iTildePhiBins = 0; iTildePhiBins < 25; iTildePhiBins++) {
            if( controlFlag15 == 1) break;
            if( (TildePhiNegativeCosTheta) < 6.28*((Double_t)iTildePhiBins + 1.)/25. ){
              fInvariantMassDistributionOnlyTildePhiHeFrameTwentyfiveBinsH[iTildePhiBins]->Fill(possibleJPsiCopyMag);
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
            fInvariantMassDistributionOnlyCosThetaCsFrameTwentyfiveBinsH[iCosThetaBins]->Fill(possibleJPsiCopyMag);
            controlFlag16 = 1;
          }
        }
        for(Int_t iPhiBins = 0; iPhiBins < 25; iPhiBins++) {
          if( controlFlag17 == 1) break;
          if( (PhiCollinsSoperValue + 3.14) < 6.28*((Double_t)iPhiBins + 1.)/25. ){
            fInvariantMassDistributionOnlyPhiCsFrameTwentyfiveBinsH[iPhiBins]->Fill(possibleJPsiCopyMag);
            controlFlag17 = 1;
          }
        }
        if( (CosThetaCollinsSoperValue > 0) || (CosThetaCollinsSoperValue == 0) ){
          for(Int_t iTildePhiBins = 0; iTildePhiBins < 25; iTildePhiBins++) {
            if( controlFlag18 == 1) break;
            if( (TildePhiPositiveCosThetaCS) < 6.28*((Double_t)iTildePhiBins + 1.)/25. ){
              fInvariantMassDistributionOnlyTildePhiCsFrameTwentyfiveBinsH[iTildePhiBins]->Fill(possibleJPsiCopyMag);
              controlFlag18 = 1;
            }
          }
        } else if ( CosThetaCollinsSoperValue < 0 ){
          for(Int_t iTildePhiBins = 0; iTildePhiBins < 25; iTildePhiBins++) {
            if( controlFlag18 == 1) break;
            if( (TildePhiNegativeCosThetaCS) < 6.28*((Double_t)iTildePhiBins + 1.)/25. ){
              fInvariantMassDistributionOnlyTildePhiCsFrameTwentyfiveBinsH[iTildePhiBins]->Fill(possibleJPsiCopyMag);
              controlFlag18 = 1;
            }
          }
        }
  }



  /* -
   * - Final incarnation of the polarisation analysis
   * - but with the 0N0N neutron emission class.
   * - This way there is no incoherent and feed-down
   * - straightaway!!
   */
  Bool_t controlFlag19 = 0;
  Bool_t controlFlag20 = 0;
  Bool_t controlFlag21 = 0;
  Bool_t controlFlag22 = 0;
  Bool_t controlFlag23 = 0;
  Bool_t controlFlag24 = 0;
  if( isZNCfired == 0 ) {
        if( isZNAfired == 0 ) {
              if ( possibleJPsiCopy.Pt() < 0.25 ) {
                    Double_t CosThetaHelicityFrameValue8 = CosThetaHelicityFrame( muonsCopy2[0], muonsCopy2[1], possibleJPsiCopy );
                    Double_t PhiHelicityFrameValue8      =   CosPhiHelicityFrame( muonsCopy2[0], muonsCopy2[1], possibleJPsiCopy );
                    // Double_t TildePhiPositiveCosTheta    = CosThetaHelicityFrameValue7 - 0.25 * TMath::Pi();
                    // Double_t TildePhiNegativeCosTheta    = CosThetaHelicityFrameValue7 - 0.75 * TMath::Pi();
                    Double_t TildePhiPositiveCosTheta2   = PhiHelicityFrameValue8 - 0.25 * 3.14;
                    Double_t TildePhiNegativeCosTheta2   = PhiHelicityFrameValue8 - 0.75 * 3.14;

                    Double_t CosThetaCollinsSoperValue2  = CosThetaCollinsSoper(  muonsCopy2[0], muonsCopy2[1], possibleJPsiCopy );
                    Double_t PhiCollinsSoperValue2       =   CosPhiCollinsSoper(  muonsCopy2[0], muonsCopy2[1], possibleJPsiCopy );
                    Double_t TildePhiPositiveCosThetaCS2 = PhiCollinsSoperValue2 - 0.25 * 3.14;
                    Double_t TildePhiNegativeCosThetaCS2 = PhiCollinsSoperValue2 - 0.75 * 3.14;

                    if( TildePhiPositiveCosTheta2 < 0. ) {
                      TildePhiPositiveCosTheta2 += 2. * TMath::Pi();
                    }
                    if( TildePhiNegativeCosTheta2 < 0. ) {
                      TildePhiNegativeCosTheta2 += 2. * TMath::Pi();
                    }

                    if( TildePhiPositiveCosThetaCS2 < 0. ) {
                      TildePhiPositiveCosThetaCS2 += 2. * TMath::Pi();
                    }
                    if( TildePhiNegativeCosThetaCS2 < 0. ) {
                      TildePhiNegativeCosThetaCS2 += 2. * TMath::Pi();
                    }

                    /* - HELICITY FRAME ANALYSIS
                     * -
                     */
                    for(Int_t iCosThetaBins = 0; iCosThetaBins < 25; iCosThetaBins++) {
                      if( controlFlag19 == 1) break;
                      if( (CosThetaHelicityFrameValue8 + 1.) < 2.*((Double_t)iCosThetaBins + 1.)/25. ) {
                        fInvariantMassDistributionOnlyCosThetaHeFrameTwentyfiveBins0N0NH[iCosThetaBins]->Fill(possibleJPsiCopyMag);
                        controlFlag19 = 1;
                      }
                    }
                    for(Int_t iPhiBins = 0; iPhiBins < 25; iPhiBins++) {
                      if( controlFlag20 == 1) break;
                      if( (PhiHelicityFrameValue8 + 3.14) < 6.28*((Double_t)iPhiBins + 1.)/25. ){
                        fInvariantMassDistributionOnlyPhiHeFrameTwentyfiveBins0N0NH[iPhiBins]->Fill(possibleJPsiCopyMag);
                        controlFlag20 = 1;
                      }
                    }
                    if( (CosThetaHelicityFrameValue8 > 0) || (CosThetaHelicityFrameValue8 == 0) ){
                      for(Int_t iTildePhiBins = 0; iTildePhiBins < 25; iTildePhiBins++) {
                        if( controlFlag21 == 1) break;
                        if( (TildePhiPositiveCosTheta2) < 6.28*((Double_t)iTildePhiBins + 1.)/25. ){
                          fInvariantMassDistributionOnlyTildePhiHeFrameTwentyfiveBins0N0NH[iTildePhiBins]->Fill(possibleJPsiCopyMag);
                          controlFlag21 = 1;
                        }
                      }
                    } else if ( CosThetaHelicityFrameValue8 < 0 ){
                      for(Int_t iTildePhiBins = 0; iTildePhiBins < 25; iTildePhiBins++) {
                        if( controlFlag21 == 1) break;
                        if( (TildePhiNegativeCosTheta2) < 6.28*((Double_t)iTildePhiBins + 1.)/25. ){
                          fInvariantMassDistributionOnlyTildePhiHeFrameTwentyfiveBins0N0NH[iTildePhiBins]->Fill(possibleJPsiCopyMag);
                          controlFlag21 = 1;
                        }
                      }
                    }

                    /* - COLLINS SOPER ANALYSIS
                     * -
                     */
                    for(Int_t iCosThetaBins = 0; iCosThetaBins < 25; iCosThetaBins++) {
                      if( controlFlag22 == 1) break;
                      if( (CosThetaCollinsSoperValue2 + 1.) < 2.*((Double_t)iCosThetaBins + 1.)/25. ) {
                        fInvariantMassDistributionOnlyCosThetaCsFrameTwentyfiveBins0N0NH[iCosThetaBins]->Fill(possibleJPsiCopyMag);
                        controlFlag22 = 1;
                      }
                    }
                    for(Int_t iPhiBins = 0; iPhiBins < 25; iPhiBins++) {
                      if( controlFlag23 == 1) break;
                      if( (PhiCollinsSoperValue2 + 3.14) < 6.28*((Double_t)iPhiBins + 1.)/25. ){
                        fInvariantMassDistributionOnlyPhiCsFrameTwentyfiveBins0N0NH[iPhiBins]->Fill(possibleJPsiCopyMag);
                        controlFlag23 = 1;
                      }
                    }
                    if( (CosThetaCollinsSoperValue2 > 0) || (CosThetaCollinsSoperValue2 == 0) ){
                      for(Int_t iTildePhiBins = 0; iTildePhiBins < 25; iTildePhiBins++) {
                        if( controlFlag24 == 1) break;
                        if( (TildePhiPositiveCosThetaCS2) < 6.28*((Double_t)iTildePhiBins + 1.)/25. ){
                          fInvariantMassDistributionOnlyTildePhiCsFrameTwentyfiveBins0N0NH[iTildePhiBins]->Fill(possibleJPsiCopyMag);
                          controlFlag24 = 1;
                        }
                      }
                    } else if ( CosThetaCollinsSoperValue2 < 0 ){
                      for(Int_t iTildePhiBins = 0; iTildePhiBins < 25; iTildePhiBins++) {
                        if( controlFlag24 == 1) break;
                        if( (TildePhiNegativeCosThetaCS2) < 6.28*((Double_t)iTildePhiBins + 1.)/25. ){
                          fInvariantMassDistributionOnlyTildePhiCsFrameTwentyfiveBins0N0NH[iTildePhiBins]->Fill(possibleJPsiCopyMag);
                          controlFlag24 = 1;
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
              Bool_t controlFlag25    = 0;
              Bool_t controlFlag25_CS = 0;
              Double_t MyVariableCosThetaBinning2[] = { -0.65, -0.35, -0.15, -0.05,
                                                         0.05,  0.15,  0.35,  0.65 };
              Double_t MyVariablePhiBinning2[] = { -3.14*1,       -3.14*19*0.05, -3.14*18*0.05, -3.14*17*0.05,
                                                   -3.14*13*0.05, -3.14*9*0.05,  -3.14*6*0.05,  -3.14*4*0.05,
                                                   -3.14*2*0.05,  -3.14*1*0.05,   0,            +3.14*1*0.05,
                                                   +3.14*2*0.05,  +3.14*4*0.05,  +3.14*6*0.05,  +3.14*9*0.05,
                                                   +3.14*13*0.05, +3.14*17*0.05, +3.14*18*0.05, +3.14*19*0.05,
                                                   +3.14*1 };
              if ( possibleJPsiCopy.Pt() < 0.25 ) {
                    Double_t CosThetaHelicityFrameValue9 = CosThetaHelicityFrame( muonsCopy2[0], muonsCopy2[1], possibleJPsiCopy );
                    Double_t PhiHelicityFrameValue9      =   CosPhiHelicityFrame( muonsCopy2[0], muonsCopy2[1], possibleJPsiCopy );
                    Double_t CosThetaCollinsSoperValue2  =  CosThetaCollinsSoper( muonsCopy2[0], muonsCopy2[1], possibleJPsiCopy );
                    Double_t PhiCollinsSoperValue2       =    CosPhiCollinsSoper( muonsCopy2[0], muonsCopy2[1], possibleJPsiCopy );
                    for(Int_t iCosThetaBins = 0; iCosThetaBins < 7; iCosThetaBins++) {
                      if( controlFlag25 == 1) break;
                      if( CosThetaHelicityFrameValue9 < MyVariableCosThetaBinning2[iCosThetaBins + 1] ){
                        for(Int_t iPhiBins = 0; iPhiBins < 20; iPhiBins++) {
                          if( controlFlag25 == 1) break;
                          if( PhiHelicityFrameValue9  < MyVariablePhiBinning2[iPhiBins + 1] ){
                              fInvariantMassDistributionForSignalExtractionHelicityFrameMyBinning0N0NH[iCosThetaBins][iPhiBins]->Fill(possibleJPsiCopyMag);
                              controlFlag25 = 1;
                          }
                        }
                      }
                    }
                    for(Int_t iCosThetaBins = 0; iCosThetaBins < 7; iCosThetaBins++) {
                      if( controlFlag25_CS == 1) break;
                      if( CosThetaCollinsSoperValue2 < MyVariableCosThetaBinning2[iCosThetaBins + 1] ){
                        for(Int_t iPhiBins = 0; iPhiBins < 20; iPhiBins++) {
                          if( controlFlag25_CS == 1) break;
                          if( PhiCollinsSoperValue2  < MyVariablePhiBinning2[iPhiBins + 1] ){
                              fInvariantMassDistributionForSignalExtractionCsFrameMyBinning0N0NH[iCosThetaBins][iPhiBins]->Fill(possibleJPsiCopyMag);
                              controlFlag25_CS = 1;
                          }
                        }
                      }
                    }
              }

        }
  }






  /* - Strict binning in Pt.
   * - This is my hope to see the bottomonium states...
   */
  if ( (possibleJPsiCopy.Pt() > 0.4) && (possibleJPsiCopy.Pt() < 0.6) ) {
    fInvariantMassDistributionStrictPtH->Fill( possibleJPsiCopy.Mag() );
  }

  if ( (possibleJPsiCopy.Mag() > 2.8) && (possibleJPsiCopy.Mag() < 3.3) && (possibleJPsiCopy.Pt() < 0.25) ) {
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

    /* - The following is a new addition.
       - The plots before were divided in 8 rapidity bins.
       - This means that the rapidity boundaries were not round numbers...
       - After having discussed with my supervisor, we agreed that
       - 10 rapidity bins are reasonable:
       - 1) it is a round number to divide 1.5 rapidity units into;
       - 2) it is not too much bigger than 8;
       - 3) the point above implies that each rapidity bin should not be
       -    that much more depleted than when divided in 8...
       - For example, 4 rapidity bins are nowhere enough,
       - 8 would be ok, but the edges of the bins are weird,
       - 15 would be perfect, but there are too many rapidity bins,
       - and this would imply scarcely populated bins!
       -
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

    /* - The next few lines of code are better explained in the header.
       - What happens here is that we fill the invariant mass distributions
       - in terms of CosTheta bins. In this way we can compute the relative
       - abundance of J/Psi and GammaGamma to the decay angular distribution.
       -
     */
    Double_t steeringVariable = CosThetaHelicityFrame(muonsCopy2[0],muonsCopy2[1],possibleJPsiCopy);
    for(Int_t iCosThetaBins = 0; iCosThetaBins < 10; iCosThetaBins++) {
      if( (steeringVariable + 1.) < 2.*((Double_t)iCosThetaBins + 1.)/10. ){
        // cout << steeringVariable << endl;
        fInvariantMassDistributionInBinsOfCosThetaHelicityFrameH[iCosThetaBins]->Fill(possibleJPsiCopyMag);
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

  }

  // post the data
  PostData(1, fOutputList);
}
//_____________________________________________________________________________
/* - The following are code snippets adapted from the AliAODDimuon class.
   - The problem is that that class was adapted specifically for the
   - inclusive people's analysis, hence it is not fit for the UPC...
   -
 */
Double_t AliAnalysisTaskUPCforward::CosThetaCollinsSoper( TLorentzVector muonPositive,
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
Double_t AliAnalysisTaskUPCforward::CosThetaHelicityFrame( TLorentzVector muonPositive,
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
Double_t AliAnalysisTaskUPCforward::CosPhiCollinsSoper( TLorentzVector muonPositive,
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
Double_t AliAnalysisTaskUPCforward::CosPhiHelicityFrame(  TLorentzVector muonPositive,
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
Double_t AliAnalysisTaskUPCforward::AccEffCorrection( Double_t CosThetaToBeWeighted )
{
  Double_t val = 0;
  Double_t par[13] = { 1.84270e-04, 1.35000e-02, 2.52572e+00, 4.49481e-01, 7.71374e-01,
                       4.10537e-01,-1.91331e-03,-4.33903e+00, 4.44258e-01,-7.64919e-01,
                       4.32631e-04, 1.48228e-02, 2.75000e+00 };
  if (        CosThetaToBeWeighted < -0.550 ) {
    val = par[0] + par[1] * ( CosThetaToBeWeighted + 0.650 ) + par[2] * ( CosThetaToBeWeighted + 0.650 ) * ( CosThetaToBeWeighted + 0.650 );
  } else if ( CosThetaToBeWeighted < -0.125 ) {
    val = par[3] + par[4] * CosThetaToBeWeighted;
  } else if ( CosThetaToBeWeighted <  0.125 ) {
    val = par[5] + par[6] * CosThetaToBeWeighted + par[7] * CosThetaToBeWeighted * CosThetaToBeWeighted;
  } else if ( CosThetaToBeWeighted <  0.550 ) {
    val = par[8] + par[9] * CosThetaToBeWeighted;
  } else {
    val = par[10] + par[11] * ( CosThetaToBeWeighted - 0.650 ) + par[12] * ( CosThetaToBeWeighted - 0.650 ) * ( CosThetaToBeWeighted - 0.650 );
  }
  return val;

}
//_____________________________________________________________________________
void AliAnalysisTaskUPCforward::Terminate(Option_t *)
{
    cout << endl;
    // terminate
    // called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________
