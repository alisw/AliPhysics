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
#include <bitset>

// my headers
#include "AliAnalysisTaskADin2018.h"



class AliAnalysisTaskADin2018;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

typedef std::bitset<32> IntBits;
// typedef std::bitset<sizeof(UShort_t)> IntBits;

ClassImp(AliAnalysisTaskADin2018) // classimp: necessary for root

//_____________________________________________________________________________
AliAnalysisTaskADin2018::AliAnalysisTaskADin2018()
    : AliAnalysisTaskSE(),
      fAOD(0),
      fOutputList(0),
      fADcheck(0),
      fNumberMuonsH(0),
      fCounterH(0),
      fEtaMuonH(0),
      fRAbsMuonH(0),
      fInvariantMassDistributionH(0),
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
      fDimuonPtDistributionRapidityHv3{0, 0, 0, 0, 0, 0},
      fDimuonPtDistributionRapidityH{0, 0, 0, 0, 0, 0},
      fInvariantMassDistributionExtendedH(0),
      fInvariantMassDistributionCoherentExtendedH(0),
      fInvariantMassDistributionIncoherentExtendedH(0),


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
      fDimuonPtDistributionZNCzeroZNAzeroHv3(0),
      fDimuonPtDistributionZNCzeroZNAanyHv3(0),
      fDimuonPtDistributionZNCanyZNAzeroHv3(0),
      fDimuonPtDistributionZNCanyZNAanyHv3(0),
      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3{0, 0, 0},
      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3{0, 0, 0},
      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3{0, 0, 0},
      fDimuonPtDistributionZNCanyZNAanyRapidityHv3{0, 0, 0},
      //_______________________________
      /* -
       * - SIDEBANDS
       */
      fDimuonPtDistributionZNCzeroZNAzeroHv2LowerSide(0),
      fDimuonPtDistributionZNCzeroZNAanyHv2LowerSide(0),
      fDimuonPtDistributionZNCanyZNAzeroHv2LowerSide(0),
      fDimuonPtDistributionZNCanyZNAanyHv2LowerSide(0),
      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv2LowerSide{0, 0, 0},
      fDimuonPtDistributionZNCzeroZNAanyRapidityHv2LowerSide{0, 0, 0},
      fDimuonPtDistributionZNCanyZNAzeroRapidityHv2LowerSide{0, 0, 0},
      fDimuonPtDistributionZNCanyZNAanyRapidityHv2LowerSide{0, 0, 0},
      fDimuonPtDistributionZNCzeroZNAzeroHv3LowerSide(0),
      fDimuonPtDistributionZNCzeroZNAanyHv3LowerSide(0),
      fDimuonPtDistributionZNCanyZNAzeroHv3LowerSide(0),
      fDimuonPtDistributionZNCanyZNAanyHv3LowerSide(0),
      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3LowerSide{0, 0, 0},
      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3LowerSide{0, 0, 0},
      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3LowerSide{0, 0, 0},
      fDimuonPtDistributionZNCanyZNAanyRapidityHv3LowerSide{0, 0, 0},

      fDimuonPtDistributionZNCzeroZNAzeroHv2HigherSide(0),
      fDimuonPtDistributionZNCzeroZNAanyHv2HigherSide(0),
      fDimuonPtDistributionZNCanyZNAzeroHv2HigherSide(0),
      fDimuonPtDistributionZNCanyZNAanyHv2HigherSide(0),
      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv2HigherSide{0, 0, 0},
      fDimuonPtDistributionZNCzeroZNAanyRapidityHv2HigherSide{0, 0, 0},
      fDimuonPtDistributionZNCanyZNAzeroRapidityHv2HigherSide{0, 0, 0},
      fDimuonPtDistributionZNCanyZNAanyRapidityHv2HigherSide{0, 0, 0},
      fDimuonPtDistributionZNCzeroZNAzeroHv3HigherSide(0),
      fDimuonPtDistributionZNCzeroZNAanyHv3HigherSide(0),
      fDimuonPtDistributionZNCanyZNAzeroHv3HigherSide(0),
      fDimuonPtDistributionZNCanyZNAanyHv3HigherSide(0),
      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3HigherSide{0, 0, 0},
      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3HigherSide{0, 0, 0},
      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3HigherSide{0, 0, 0},
      fDimuonPtDistributionZNCanyZNAanyRapidityHv3HigherSide{0, 0, 0},
      //_______________________________
      fDimuonPtDistributionCoherentZNCzeroZNAzeroH(0),
      fDimuonPtDistributionCoherentZNCzeroZNAanyH(0),
      fDimuonPtDistributionCoherentZNCanyZNAzeroH(0),
      fDimuonPtDistributionCoherentZNCanyZNAanyH(0),
      fDimuonPtDistributionIncoherentZNCzeroZNAzeroH(0),
      fDimuonPtDistributionIncoherentZNCzeroZNAanyH(0),
      fDimuonPtDistributionIncoherentZNCanyZNAzeroH(0),
      fDimuonPtDistributionIncoherentZNCanyZNAanyH(0),
      fADmultiplicityH{ 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0 },
      fADmultiplicity0N0NclassH{ 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0 },
      fADmultiplicity0NXNclassH{ 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0 },
      fADmultiplicityXN0NclassH{ 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0 },
      fADmultiplicityXNXNclassH{ 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0 },
      fADmultiplicity0N0NclassRapidityH{ 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
                                         0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
                                         0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0 },
      fADmultiplicity0NXNclassRapidityH{ 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
                                         0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
                                         0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0 },
      fADmultiplicityXN0NclassRapidityH{ 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
                                         0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
                                         0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0 },
      fADmultiplicityXNXNclassRapidityH{ 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
                                         0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
                                         0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0 },
      fADmultiplicityTotalH(0),
      fADmultiplicity0N0NclassTotalH(0),
      fADmultiplicity0NXNclassTotalH(0),
      fADmultiplicityXN0NclassTotalH(0),
      fADmultiplicityXNXNclassTotalH(0),
      fADAmultiplicityTotalH(0),
      fADAmultiplicity0N0NclassTotalH(0),
      fADAmultiplicity0NXNclassTotalH(0),
      fADAmultiplicityXN0NclassTotalH(0),
      fADAmultiplicityXNXNclassTotalH(0),
      fADCmultiplicityTotalH(0),
      fADCmultiplicity0N0NclassTotalH(0),
      fADCmultiplicity0NXNclassTotalH(0),
      fADCmultiplicityXN0NclassTotalH(0),
      fADCmultiplicityXNXNclassTotalH(0),
      fADAmultiplicityTotalVsZNAenergyH(0),
      fADCmultiplicityTotalVsZNCenergyH(0),

      fVZEROmultiplicityTotalH(0),
      fVZEROmultiplicity0N0NclassTotalH(0),
      fVZEROmultiplicity0NXNclassTotalH(0),
      fVZEROmultiplicityXN0NclassTotalH(0),
      fVZEROmultiplicityXNXNclassTotalH(0),
      fVZEROAmultiplicityTotalH(0),
      fVZEROAmultiplicity0N0NclassTotalH(0),
      fVZEROAmultiplicity0NXNclassTotalH(0),
      fVZEROAmultiplicityXN0NclassTotalH(0),
      fVZEROAmultiplicityXNXNclassTotalH(0),
      fVZEROCmultiplicityTotalH(0),
      fVZEROCmultiplicity0N0NclassTotalH(0),
      fVZEROCmultiplicity0NXNclassTotalH(0),
      fVZEROCmultiplicityXN0NclassTotalH(0),
      fVZEROCmultiplicityXNXNclassTotalH(0),
      fVZEROAmultiplicityTotalVsZNAenergyH(0),
      fVZEROCmultiplicityTotalVsZNCenergyH(0),

      fADAmultiplicityVsVZEROAmultiplicityH(0),
      fADCmultiplicityVsVZEROCmultiplicityH(0),


      fADAmultiplicityTotalVsZNAenergyH_ADAno(0),
      fADCmultiplicityTotalVsZNCenergyH_ADCno(0),
      fVZEROAmultiplicityTotalVsZNAenergyH_VZEROAno(0),
      fVZEROCmultiplicityTotalVsZNCenergyH_VZEROCno(0),

      fADAmultiplicityVsVZEROAmultiplicityH_VZEROAno(0),
      fADAmultiplicityVsVZEROAmultiplicityH_ADAno_VZEROAno(0),
      fADCmultiplicityVsVZEROCmultiplicityH_VZEROCno(0),
      fADCmultiplicityVsVZEROCmultiplicityH_ADCno_VZEROCno(0),


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
      fV0TotalNCells(0)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}

//_____________________________________________________________________________
AliAnalysisTaskADin2018::AliAnalysisTaskADin2018(const char* name)
    : AliAnalysisTaskSE(name),
      fAOD(0),
      fOutputList(0),
      fADcheck(0),
      fNumberMuonsH(0),
      fCounterH(0),
      fEtaMuonH(0),
      fRAbsMuonH(0),
      fInvariantMassDistributionH(0),
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
      fDimuonPtDistributionRapidityHv3{0, 0, 0, 0, 0, 0},
      fDimuonPtDistributionRapidityH{0, 0, 0, 0, 0, 0},
      fInvariantMassDistributionExtendedH(0),
      fInvariantMassDistributionCoherentExtendedH(0),
      fInvariantMassDistributionIncoherentExtendedH(0),



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
      fDimuonPtDistributionZNCzeroZNAzeroHv3(0),
      fDimuonPtDistributionZNCzeroZNAanyHv3(0),
      fDimuonPtDistributionZNCanyZNAzeroHv3(0),
      fDimuonPtDistributionZNCanyZNAanyHv3(0),
      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3{0, 0, 0},
      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3{0, 0, 0},
      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3{0, 0, 0},
      fDimuonPtDistributionZNCanyZNAanyRapidityHv3{0, 0, 0},
      //_______________________________
      /* -
       * - SIDEBANDS
       */
      fDimuonPtDistributionZNCzeroZNAzeroHv2LowerSide(0),
      fDimuonPtDistributionZNCzeroZNAanyHv2LowerSide(0),
      fDimuonPtDistributionZNCanyZNAzeroHv2LowerSide(0),
      fDimuonPtDistributionZNCanyZNAanyHv2LowerSide(0),
      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv2LowerSide{0, 0, 0},
      fDimuonPtDistributionZNCzeroZNAanyRapidityHv2LowerSide{0, 0, 0},
      fDimuonPtDistributionZNCanyZNAzeroRapidityHv2LowerSide{0, 0, 0},
      fDimuonPtDistributionZNCanyZNAanyRapidityHv2LowerSide{0, 0, 0},
      fDimuonPtDistributionZNCzeroZNAzeroHv3LowerSide(0),
      fDimuonPtDistributionZNCzeroZNAanyHv3LowerSide(0),
      fDimuonPtDistributionZNCanyZNAzeroHv3LowerSide(0),
      fDimuonPtDistributionZNCanyZNAanyHv3LowerSide(0),
      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3LowerSide{0, 0, 0},
      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3LowerSide{0, 0, 0},
      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3LowerSide{0, 0, 0},
      fDimuonPtDistributionZNCanyZNAanyRapidityHv3LowerSide{0, 0, 0},

      fDimuonPtDistributionZNCzeroZNAzeroHv2HigherSide(0),
      fDimuonPtDistributionZNCzeroZNAanyHv2HigherSide(0),
      fDimuonPtDistributionZNCanyZNAzeroHv2HigherSide(0),
      fDimuonPtDistributionZNCanyZNAanyHv2HigherSide(0),
      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv2HigherSide{0, 0, 0},
      fDimuonPtDistributionZNCzeroZNAanyRapidityHv2HigherSide{0, 0, 0},
      fDimuonPtDistributionZNCanyZNAzeroRapidityHv2HigherSide{0, 0, 0},
      fDimuonPtDistributionZNCanyZNAanyRapidityHv2HigherSide{0, 0, 0},
      fDimuonPtDistributionZNCzeroZNAzeroHv3HigherSide(0),
      fDimuonPtDistributionZNCzeroZNAanyHv3HigherSide(0),
      fDimuonPtDistributionZNCanyZNAzeroHv3HigherSide(0),
      fDimuonPtDistributionZNCanyZNAanyHv3HigherSide(0),
      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3HigherSide{0, 0, 0},
      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3HigherSide{0, 0, 0},
      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3HigherSide{0, 0, 0},
      fDimuonPtDistributionZNCanyZNAanyRapidityHv3HigherSide{0, 0, 0},
      //_______________________________
      fDimuonPtDistributionCoherentZNCzeroZNAzeroH(0),
      fDimuonPtDistributionCoherentZNCzeroZNAanyH(0),
      fDimuonPtDistributionCoherentZNCanyZNAzeroH(0),
      fDimuonPtDistributionCoherentZNCanyZNAanyH(0),
      fDimuonPtDistributionIncoherentZNCzeroZNAzeroH(0),
      fDimuonPtDistributionIncoherentZNCzeroZNAanyH(0),
      fDimuonPtDistributionIncoherentZNCanyZNAzeroH(0),
      fDimuonPtDistributionIncoherentZNCanyZNAanyH(0),
      fADmultiplicityH{ 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0 },
      fADmultiplicity0N0NclassH{ 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0 },
      fADmultiplicity0NXNclassH{ 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0 },
      fADmultiplicityXN0NclassH{ 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0 },
      fADmultiplicityXNXNclassH{ 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0 },
      fADmultiplicity0N0NclassRapidityH{ 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
                                         0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
                                         0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0 },
      fADmultiplicity0NXNclassRapidityH{ 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
                                         0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
                                         0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0 },
      fADmultiplicityXN0NclassRapidityH{ 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
                                         0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
                                         0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0 },
      fADmultiplicityXNXNclassRapidityH{ 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
                                         0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
                                         0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0 },
      fADmultiplicityTotalH(0),
      fADmultiplicity0N0NclassTotalH(0),
      fADmultiplicity0NXNclassTotalH(0),
      fADmultiplicityXN0NclassTotalH(0),
      fADmultiplicityXNXNclassTotalH(0),
      fADAmultiplicityTotalH(0),
      fADAmultiplicity0N0NclassTotalH(0),
      fADAmultiplicity0NXNclassTotalH(0),
      fADAmultiplicityXN0NclassTotalH(0),
      fADAmultiplicityXNXNclassTotalH(0),
      fADCmultiplicityTotalH(0),
      fADCmultiplicity0N0NclassTotalH(0),
      fADCmultiplicity0NXNclassTotalH(0),
      fADCmultiplicityXN0NclassTotalH(0),
      fADCmultiplicityXNXNclassTotalH(0),
      fADAmultiplicityTotalVsZNAenergyH(0),
      fADCmultiplicityTotalVsZNCenergyH(0),

      fVZEROmultiplicityTotalH(0),
      fVZEROmultiplicity0N0NclassTotalH(0),
      fVZEROmultiplicity0NXNclassTotalH(0),
      fVZEROmultiplicityXN0NclassTotalH(0),
      fVZEROmultiplicityXNXNclassTotalH(0),
      fVZEROAmultiplicityTotalH(0),
      fVZEROAmultiplicity0N0NclassTotalH(0),
      fVZEROAmultiplicity0NXNclassTotalH(0),
      fVZEROAmultiplicityXN0NclassTotalH(0),
      fVZEROAmultiplicityXNXNclassTotalH(0),
      fVZEROCmultiplicityTotalH(0),
      fVZEROCmultiplicity0N0NclassTotalH(0),
      fVZEROCmultiplicity0NXNclassTotalH(0),
      fVZEROCmultiplicityXN0NclassTotalH(0),
      fVZEROCmultiplicityXNXNclassTotalH(0),
      fVZEROAmultiplicityTotalVsZNAenergyH(0),
      fVZEROCmultiplicityTotalVsZNCenergyH(0),

      fADAmultiplicityVsVZEROAmultiplicityH(0),
      fADCmultiplicityVsVZEROCmultiplicityH(0),

      fADAmultiplicityTotalVsZNAenergyH_ADAno(0),
      fADCmultiplicityTotalVsZNCenergyH_ADCno(0),
      fVZEROAmultiplicityTotalVsZNAenergyH_VZEROAno(0),
      fVZEROCmultiplicityTotalVsZNCenergyH_VZEROCno(0),

      fADAmultiplicityVsVZEROAmultiplicityH_VZEROAno(0),
      fADAmultiplicityVsVZEROAmultiplicityH_ADAno_VZEROAno(0),
      fADCmultiplicityVsVZEROCmultiplicityH_VZEROCno(0),
      fADCmultiplicityVsVZEROCmultiplicityH_ADCno_VZEROCno(0),

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
      fV0TotalNCells(0)
{

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
AliAnalysisTaskADin2018::~AliAnalysisTaskADin2018()
{
    // destructor
    if(fOutputList)    {delete fOutputList;}     	// at the end of your task, it is deleted
    if(fMuonTrackCuts) {delete fMuonTrackCuts;}   // from memory by calling this function
}
//_____________________________________________________________________________
void AliAnalysisTaskADin2018::UserCreateOutputObjects()
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


  /* - These histograms have an EXTENDED range (0,20)->(0,40)
     -
   */
  fInvariantMassDistributionExtendedH = new TH1F("fInvariantMassDistributionExtendedH", "fInvariantMassDistributionExtendedH", 4000, 0, 40);
  fOutputList->Add(fInvariantMassDistributionExtendedH);

  fInvariantMassDistributionCoherentExtendedH = new TH1F("fInvariantMassDistributionCoherentExtendedH", "fInvariantMassDistributionCoherentExtendedH", 4000, 0, 40);
  fOutputList->Add(fInvariantMassDistributionCoherentExtendedH);

  fInvariantMassDistributionIncoherentExtendedH = new TH1F("fInvariantMassDistributionIncoherentExtendedH", "fInvariantMassDistributionIncoherentExtendedH", 4000, 0, 40);
  fOutputList->Add(fInvariantMassDistributionIncoherentExtendedH);




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

  fADmultiplicityTotalH = new TH1F("fADmultiplicityTotalH","fADmultiplicityTotalH",1000000, 0, 100000);
  fOutputList->Add(fADmultiplicityTotalH);

  fADmultiplicity0N0NclassTotalH = new TH1F("fADmultiplicity0N0NclassTotalH","fADmultiplicity0N0NclassTotalH",1000000, 0, 100000);
  fOutputList->Add(fADmultiplicity0N0NclassTotalH);

  fADmultiplicity0NXNclassTotalH = new TH1F("fADmultiplicity0NXNclassTotalH","fADmultiplicity0NXNclassTotalH",1000000, 0, 100000);
  fOutputList->Add(fADmultiplicity0NXNclassTotalH);

  fADmultiplicityXN0NclassTotalH = new TH1F("fADmultiplicityXN0NclassTotalH","fADmultiplicityXN0NclassTotalH",1000000, 0, 100000);
  fOutputList->Add(fADmultiplicityXN0NclassTotalH);

  fADmultiplicityXNXNclassTotalH = new TH1F("fADmultiplicityXNXNclassTotalH","fADmultiplicityXNXNclassTotalH",1000000, 0, 100000);
  fOutputList->Add(fADmultiplicityXNXNclassTotalH);

  fADAmultiplicityTotalH = new TH1F("fADAmultiplicityTotalH","fADAmultiplicityTotalH",1000000, 0, 100000);
  fOutputList->Add(fADAmultiplicityTotalH);

  fADAmultiplicity0N0NclassTotalH = new TH1F("fADAmultiplicity0N0NclassTotalH","fADAmultiplicity0N0NclassTotalH",1000000, 0, 100000);
  fOutputList->Add(fADAmultiplicity0N0NclassTotalH);

  fADAmultiplicity0NXNclassTotalH = new TH1F("fADAmultiplicity0NXNclassTotalH","fADAmultiplicity0NXNclassTotalH",1000000, 0, 100000);
  fOutputList->Add(fADAmultiplicity0NXNclassTotalH);

  fADAmultiplicityXN0NclassTotalH = new TH1F("fADAmultiplicityXN0NclassTotalH","fADAmultiplicityXN0NclassTotalH",1000000, 0, 100000);
  fOutputList->Add(fADAmultiplicityXN0NclassTotalH);

  fADAmultiplicityXNXNclassTotalH = new TH1F("fADAmultiplicityXNXNclassTotalH","fADAmultiplicityXNXNclassTotalH",1000000, 0, 100000);
  fOutputList->Add(fADAmultiplicityXNXNclassTotalH);

  fADCmultiplicityTotalH = new TH1F("fADCmultiplicityTotalH","fADCmultiplicityTotalH",1000000, 0, 100000);
  fOutputList->Add(fADCmultiplicityTotalH);

  fADCmultiplicity0N0NclassTotalH = new TH1F("fADCmultiplicity0N0NclassTotalH","fADCmultiplicity0N0NclassTotalH",1000000, 0, 100000);
  fOutputList->Add(fADCmultiplicity0N0NclassTotalH);

  fADCmultiplicity0NXNclassTotalH = new TH1F("fADCmultiplicity0NXNclassTotalH","fADCmultiplicity0NXNclassTotalH",1000000, 0, 100000);
  fOutputList->Add(fADCmultiplicity0NXNclassTotalH);

  fADCmultiplicityXN0NclassTotalH = new TH1F("fADCmultiplicityXN0NclassTotalH","fADCmultiplicityXN0NclassTotalH",1000000, 0, 100000);
  fOutputList->Add(fADCmultiplicityXN0NclassTotalH);

  fADCmultiplicityXNXNclassTotalH = new TH1F("fADCmultiplicityXNXNclassTotalH","fADCmultiplicityXNXNclassTotalH",1000000, 0, 100000);
  fOutputList->Add(fADCmultiplicityXNXNclassTotalH);

  fADCmultiplicityTotalVsZNCenergyH = new TH2F("fADCmultiplicityTotalVsZNCenergyH","", 500, 0, 50, 20000, -10000, 40000);
  fOutputList->Add(fADCmultiplicityTotalVsZNCenergyH);

  fADAmultiplicityTotalVsZNAenergyH = new TH2F("fADAmultiplicityTotalVsZNAenergyH","", 500, 0, 50, 20000, -10000, 40000);
  fOutputList->Add(fADAmultiplicityTotalVsZNAenergyH);

  //_______________________________
  /*
   * VZERO plots.
   */
  fVZEROmultiplicityTotalH = new TH1F("fVZEROmultiplicityTotalH","fVZEROmultiplicityTotalH",1000000, 0, 100000);
  fOutputList->Add(fVZEROmultiplicityTotalH);

  fVZEROmultiplicity0N0NclassTotalH = new TH1F("fVZEROmultiplicity0N0NclassTotalH","fVZEROmultiplicity0N0NclassTotalH",1000000, 0, 100000);
  fOutputList->Add(fVZEROmultiplicity0N0NclassTotalH);

  fVZEROmultiplicity0NXNclassTotalH = new TH1F("fVZEROmultiplicity0NXNclassTotalH","fVZEROmultiplicity0NXNclassTotalH",1000000, 0, 100000);
  fOutputList->Add(fVZEROmultiplicity0NXNclassTotalH);

  fVZEROmultiplicityXN0NclassTotalH = new TH1F("fVZEROmultiplicityXN0NclassTotalH","fVZEROmultiplicityXN0NclassTotalH",1000000, 0, 100000);
  fOutputList->Add(fVZEROmultiplicityXN0NclassTotalH);

  fVZEROmultiplicityXNXNclassTotalH = new TH1F("fVZEROmultiplicityXNXNclassTotalH","fVZEROmultiplicityXNXNclassTotalH",1000000, 0, 100000);
  fOutputList->Add(fVZEROmultiplicityXNXNclassTotalH);

  fVZEROAmultiplicityTotalH = new TH1F("fVZEROAmultiplicityTotalH","fVZEROAmultiplicityTotalH",1000000, 0, 100000);
  fOutputList->Add(fVZEROAmultiplicityTotalH);

  fVZEROAmultiplicity0N0NclassTotalH = new TH1F("fVZEROAmultiplicity0N0NclassTotalH","fVZEROAmultiplicity0N0NclassTotalH",1000000, 0, 100000);
  fOutputList->Add(fVZEROAmultiplicity0N0NclassTotalH);

  fVZEROAmultiplicity0NXNclassTotalH = new TH1F("fVZEROAmultiplicity0NXNclassTotalH","fVZEROAmultiplicity0NXNclassTotalH",1000000, 0, 100000);
  fOutputList->Add(fVZEROAmultiplicity0NXNclassTotalH);

  fVZEROAmultiplicityXN0NclassTotalH = new TH1F("fVZEROAmultiplicityXN0NclassTotalH","fVZEROAmultiplicityXN0NclassTotalH",1000000, 0, 100000);
  fOutputList->Add(fVZEROAmultiplicityXN0NclassTotalH);

  fVZEROAmultiplicityXNXNclassTotalH = new TH1F("fVZEROAmultiplicityXNXNclassTotalH","fVZEROAmultiplicityXNXNclassTotalH",1000000, 0, 100000);
  fOutputList->Add(fVZEROAmultiplicityXNXNclassTotalH);

  fVZEROCmultiplicityTotalH = new TH1F("fVZEROCmultiplicityTotalH","fVZEROCmultiplicityTotalH",1000000, 0, 100000);
  fOutputList->Add(fVZEROCmultiplicityTotalH);

  fVZEROCmultiplicity0N0NclassTotalH = new TH1F("fVZEROCmultiplicity0N0NclassTotalH","fVZEROCmultiplicity0N0NclassTotalH",1000000, 0, 100000);
  fOutputList->Add(fVZEROCmultiplicity0N0NclassTotalH);

  fVZEROCmultiplicity0NXNclassTotalH = new TH1F("fVZEROCmultiplicity0NXNclassTotalH","fVZEROCmultiplicity0NXNclassTotalH",1000000, 0, 100000);
  fOutputList->Add(fVZEROCmultiplicity0NXNclassTotalH);

  fVZEROCmultiplicityXN0NclassTotalH = new TH1F("fVZEROCmultiplicityXN0NclassTotalH","fVZEROCmultiplicityXN0NclassTotalH",1000000, 0, 100000);
  fOutputList->Add(fVZEROCmultiplicityXN0NclassTotalH);

  fVZEROCmultiplicityXNXNclassTotalH = new TH1F("fVZEROCmultiplicityXNXNclassTotalH","fVZEROCmultiplicityXNXNclassTotalH",1000000, 0, 100000);
  fOutputList->Add(fVZEROCmultiplicityXNXNclassTotalH);

  fVZEROCmultiplicityTotalVsZNCenergyH = new TH2F("fVZEROCmultiplicityTotalVsZNCenergyH","fVZEROCmultiplicityTotalVsZNCenergyH", 500, 0, 50, 1000, -10000, 40000);
  fOutputList->Add(fVZEROCmultiplicityTotalVsZNCenergyH);

  fVZEROAmultiplicityTotalVsZNAenergyH = new TH2F("fVZEROAmultiplicityTotalVsZNAenergyH","fVZEROAmultiplicityTotalVsZNAenergyH", 500, 0, 50, 1000, -10000, 40000);
  fOutputList->Add(fVZEROAmultiplicityTotalVsZNAenergyH);

  fADCmultiplicityVsVZEROCmultiplicityH = new TH2F("fADCmultiplicityVsVZEROCmultiplicityH","fADCmultiplicityVsVZEROCmultiplicityH", 500, 0, 50, 500, 0, 50);
  fOutputList->Add(fADCmultiplicityVsVZEROCmultiplicityH);

  fADAmultiplicityVsVZEROAmultiplicityH = new TH2F("fADAmultiplicityVsVZEROAmultiplicityH","fADAmultiplicityVsVZEROAmultiplicityH", 500, 0, 50, 500, 0, 50);
  fOutputList->Add(fADAmultiplicityVsVZEROAmultiplicityH);

  //_______________________________
  /* -
   * - TRIGGER BITS check
   * -
   */
  fADCmultiplicityTotalVsZNCenergyH_ADCno = new TH2F("fADCmultiplicityTotalVsZNCenergyH_ADCno","", 500, 0, 50, 20000, -10000, 40000);
  fOutputList->Add(fADCmultiplicityTotalVsZNCenergyH_ADCno);

  fADAmultiplicityTotalVsZNAenergyH_ADAno = new TH2F("fADAmultiplicityTotalVsZNAenergyH_ADAno","", 500, 0, 50, 20000, -10000, 40000);
  fOutputList->Add(fADAmultiplicityTotalVsZNAenergyH_ADAno);

  fVZEROCmultiplicityTotalVsZNCenergyH_VZEROCno = new TH2F("fVZEROCmultiplicityTotalVsZNCenergyH_VZEROCno","fVZEROCmultiplicityTotalVsZNCenergyH_VZEROCno", 500, 0, 50, 1000, -10000, 40000);
  fOutputList->Add(fVZEROCmultiplicityTotalVsZNCenergyH_VZEROCno);

  fVZEROAmultiplicityTotalVsZNAenergyH_VZEROAno = new TH2F("fVZEROAmultiplicityTotalVsZNAenergyH_VZEROAno","fVZEROAmultiplicityTotalVsZNAenergyH_VZEROAno", 500, 0, 50, 1000, -10000, 40000);
  fOutputList->Add(fVZEROAmultiplicityTotalVsZNAenergyH_VZEROAno);

  fADCmultiplicityVsVZEROCmultiplicityH_VZEROCno = new TH2F("fADCmultiplicityVsVZEROCmultiplicityH_VZEROCno","fADCmultiplicityVsVZEROCmultiplicityH_VZEROCno", 500, 0, 50, 500, 0, 50);
  fOutputList->Add(fADCmultiplicityVsVZEROCmultiplicityH_VZEROCno);

  fADCmultiplicityVsVZEROCmultiplicityH_ADCno_VZEROCno = new TH2F("fADCmultiplicityVsVZEROCmultiplicityH_ADCno_VZEROCno","fADCmultiplicityVsVZEROCmultiplicityH_ADCno_VZEROCno", 500, 0, 50, 500, 0, 50);
  fOutputList->Add(fADCmultiplicityVsVZEROCmultiplicityH_ADCno_VZEROCno);

  fADAmultiplicityVsVZEROAmultiplicityH_VZEROAno = new TH2F("fADAmultiplicityVsVZEROAmultiplicityH_VZEROAno","fADAmultiplicityVsVZEROAmultiplicityH_VZEROAno", 500, 0, 50, 500, 0, 50);
  fOutputList->Add(fADAmultiplicityVsVZEROAmultiplicityH_VZEROAno);

  fADAmultiplicityVsVZEROAmultiplicityH_ADAno_VZEROAno = new TH2F("fADAmultiplicityVsVZEROAmultiplicityH_ADAno_VZEROAno","fADAmultiplicityVsVZEROAmultiplicityH_ADAno_VZEROAno", 500, 0, 50, 500, 0, 50);
  fOutputList->Add(fADAmultiplicityVsVZEROAmultiplicityH_ADAno_VZEROAno);
  //_______________________________


  for( Int_t iChannel = 0; iChannel < 16; iChannel++ ){
    fADmultiplicityH[iChannel] = new TH1F(
              Form("fADmultiplicityH_%d", iChannel),
              Form("fADmultiplicityH_%d", iChannel),
              100000, 0, 10000
              );
    fOutputList->Add(fADmultiplicityH[iChannel]);

    fADmultiplicity0N0NclassH[iChannel] = new TH1F(
              Form("fADmultiplicity0N0NclassH_%d", iChannel),
              Form("fADmultiplicity0N0NclassH_%d", iChannel),
              100000, 0, 10000
              );
    fOutputList->Add(fADmultiplicity0N0NclassH[iChannel]);

    fADmultiplicity0NXNclassH[iChannel] = new TH1F(
              Form("fADmultiplicity0NXNclassH_%d", iChannel),
              Form("fADmultiplicity0NXNclassH_%d", iChannel),
              100000, 0, 10000
              );
    fOutputList->Add(fADmultiplicity0NXNclassH[iChannel]);

    fADmultiplicityXN0NclassH[iChannel] = new TH1F(
              Form("fADmultiplicityXN0NclassH_%d", iChannel),
              Form("fADmultiplicityXN0NclassH_%d", iChannel),
              100000, 0, 10000
              );
    fOutputList->Add(fADmultiplicityXN0NclassH[iChannel]);

    fADmultiplicityXNXNclassH[iChannel] = new TH1F(
              Form("fADmultiplicityXNXNclassH_%d", iChannel),
              Form("fADmultiplicityXNXNclassH_%d", iChannel),
              100000, 0, 10000
              );
    fOutputList->Add(fADmultiplicityXNXNclassH[iChannel]);

    for( Int_t iRapidity = 0; iRapidity < 3; iRapidity++ )
    {
      fADmultiplicity0N0NclassRapidityH[iChannel + 16*iRapidity] = new TH1F(
                Form("fADmultiplicity0N0NclassRapidityH_%d", iChannel + 16*iRapidity),
                Form("fADmultiplicity0N0NclassRapidityH_%d", iChannel + 16*iRapidity),
                100000, 0, 10000
                );
      fOutputList->Add(fADmultiplicity0N0NclassRapidityH[iChannel + 16*iRapidity]);

      fADmultiplicity0NXNclassRapidityH[iChannel + 16*iRapidity] = new TH1F(
                Form("fADmultiplicity0NXNclassRapidityH_%d", iChannel + 16*iRapidity),
                Form("fADmultiplicity0NXNclassRapidityH_%d", iChannel + 16*iRapidity),
                100000, 0, 10000
                );
      fOutputList->Add(fADmultiplicity0NXNclassRapidityH[iChannel + 16*iRapidity]);

      fADmultiplicityXN0NclassRapidityH[iChannel + 16*iRapidity] = new TH1F(
                Form("fADmultiplicityXN0NclassRapidityH_%d", iChannel + 16*iRapidity),
                Form("fADmultiplicityXN0NclassRapidityH_%d", iChannel + 16*iRapidity),
                100000, 0, 10000
                );
      fOutputList->Add(fADmultiplicityXN0NclassRapidityH[iChannel + 16*iRapidity]);

      fADmultiplicityXNXNclassRapidityH[iChannel + 16*iRapidity] = new TH1F(
                Form("fADmultiplicityXNXNclassRapidityH_%d", iChannel + 16*iRapidity),
                Form("fADmultiplicityXNXNclassRapidityH_%d", iChannel + 16*iRapidity),
                100000, 0, 10000
                );
      fOutputList->Add(fADmultiplicityXNXNclassRapidityH[iChannel + 16*iRapidity]);
    }
  }


  // Float_t PtBins[]    = { 0.000, 0.025, 0.050, 0.075, 0.100, 0.125, 0.150, 0.175,
  //                         0.200, 0.225, 0.250, 0.275, 0.300, 0.325, 0.350, 0.375,
  //                         0.400, 0.425, 0.450, 0.475, 0.500, 0.525, 0.550, 0.575,
  //                         0.600, 0.625, 0.650, 0.675, 0.700, 0.725, 0.750, 0.775,
  //                         0.800, 0.825, 0.850, 0.875, 0.900, 0.925, 0.950, 0.975,
  //                         1.000,
  //                       };
  // Float_t PtBins[]    = { 0.000, 0.025, 0.050, 0.075, 0.100, 0.125, 0.150, 0.175,
  //                         0.200, 0.225, 0.250, 0.275, 0.300, 0.325, 0.350, 0.375,
  //                         0.400, 0.425, 0.450, 0.475, 0.500, 0.575, 0.650, 0.725,
  //                         0.800, 0.875, 0.950, 1.100, 1.250, 1.400, 1.600, 1.800,
  //                         2.000, 2.500, 3.000, 3.500, 4.000, 5.000
  //                       };
  Float_t PtBins[]    = { 0.000, 0.025, 0.050, 0.075, 0.100, 0.125, 0.150, 0.175,
                          0.200, 0.225, 0.250, 0.275, 0.350, 0.425, 0.500, 0.575,
                          0.650, 0.725,
                          0.800, 0.875, 0.950, 1.100, 1.250, 1.400, 1.600, 1.800,
                          2.000, 2.500, 3.000, 3.500, 4.000, 5.000
                        };
  Int_t   PtBinNumber = sizeof(PtBins)/sizeof(Float_t) - 1; // or just = 9

  for( Int_t iRapidity = 0; iRapidity < 6; iRapidity++ ){
    fDimuonPtDistributionRapidityHv3[iRapidity] = new TH1F(
              Form("fDimuonPtDistributionRapidityHv3_%d", iRapidity),
              Form("fDimuonPtDistributionRapidityHv3_%d", iRapidity),
              PtBinNumber, PtBins
              );
    fOutputList->Add(fDimuonPtDistributionRapidityHv3[iRapidity]);
  }

  for( Int_t iRapidity = 0; iRapidity < 6; iRapidity++ ){
    fDimuonPtDistributionRapidityH[iRapidity] = new TH1F(
              Form("fDimuonPtDistributionRapidityH_%d", iRapidity),
              Form("fDimuonPtDistributionRapidityH_%d", iRapidity),
              4000, 0, 20
              );
    fOutputList->Add(fDimuonPtDistributionRapidityH[iRapidity]);
  }


  fDimuonPtDistributionZNCzeroZNAzeroHv3 = new TH1F("fDimuonPtDistributionZNCzeroZNAzeroHv3", "fDimuonPtDistributionZNCzeroZNAzeroHv3", PtBinNumber, PtBins);
  fOutputList->Add(fDimuonPtDistributionZNCzeroZNAzeroHv3);

  fDimuonPtDistributionZNCzeroZNAanyHv3 = new TH1F("fDimuonPtDistributionZNCzeroZNAanyHv3", "fDimuonPtDistributionZNCzeroZNAanyHv3", PtBinNumber, PtBins);
  fOutputList->Add(fDimuonPtDistributionZNCzeroZNAanyHv3);

  fDimuonPtDistributionZNCanyZNAzeroHv3 = new TH1F("fDimuonPtDistributionZNCanyZNAzeroHv3", "fDimuonPtDistributionZNCanyZNAzeroHv3", PtBinNumber, PtBins);
  fOutputList->Add(fDimuonPtDistributionZNCanyZNAzeroHv3);

  fDimuonPtDistributionZNCanyZNAanyHv3 = new TH1F("fDimuonPtDistributionZNCanyZNAanyHv3", "fDimuonPtDistributionZNCanyZNAanyHv3", PtBinNumber, PtBins);
  fOutputList->Add(fDimuonPtDistributionZNCanyZNAanyHv3);

  for( Int_t iRapidity = 0; iRapidity < 3; iRapidity++ ){
    fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3[iRapidity] = new TH1F(
              Form("fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3_%d", iRapidity),
              Form("fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3_%d", iRapidity),
              PtBinNumber, PtBins
              );
    fOutputList->Add(fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3[iRapidity]);
  }

  for( Int_t iRapidity = 0; iRapidity < 3; iRapidity++ ){
    fDimuonPtDistributionZNCzeroZNAanyRapidityHv3[iRapidity] = new TH1F(
              Form("fDimuonPtDistributionZNCzeroZNAanyRapidityHv3_%d", iRapidity),
              Form("fDimuonPtDistributionZNCzeroZNAanyRapidityHv3_%d", iRapidity),
              PtBinNumber, PtBins
              );
    fOutputList->Add(fDimuonPtDistributionZNCzeroZNAanyRapidityHv3[iRapidity]);
  }

  for( Int_t iRapidity = 0; iRapidity < 3; iRapidity++ ){
    fDimuonPtDistributionZNCanyZNAzeroRapidityHv3[iRapidity] = new TH1F(
              Form("fDimuonPtDistributionZNCanyZNAzeroRapidityHv3_%d", iRapidity),
              Form("fDimuonPtDistributionZNCanyZNAzeroRapidityHv3_%d", iRapidity),
              PtBinNumber, PtBins
              );
    fOutputList->Add(fDimuonPtDistributionZNCanyZNAzeroRapidityHv3[iRapidity]);
  }

  for( Int_t iRapidity = 0; iRapidity < 3; iRapidity++ ){
    fDimuonPtDistributionZNCanyZNAanyRapidityHv3[iRapidity] = new TH1F(
              Form("fDimuonPtDistributionZNCanyZNAanyRapidityHv3_%d", iRapidity),
              Form("fDimuonPtDistributionZNCanyZNAanyRapidityHv3_%d", iRapidity),
              PtBinNumber, PtBins
              );
    fOutputList->Add(fDimuonPtDistributionZNCanyZNAanyRapidityHv3[iRapidity]);
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





  //_______________________________
  /* -
   * - SIDEBANDS
   * -
   */
  //
  // LOWER SIDE
  //
  fDimuonPtDistributionZNCzeroZNAzeroHv2LowerSide = new TH1F("fDimuonPtDistributionZNCzeroZNAzeroHv2LowerSide", "fDimuonPtDistributionZNCzeroZNAzeroHv2LowerSide", 4000, 0, 20);
  fOutputList->Add(fDimuonPtDistributionZNCzeroZNAzeroHv2LowerSide);

  fDimuonPtDistributionZNCzeroZNAanyHv2LowerSide = new TH1F("fDimuonPtDistributionZNCzeroZNAanyHv2LowerSide", "fDimuonPtDistributionZNCzeroZNAanyHv2LowerSide", 4000, 0, 20);
  fOutputList->Add(fDimuonPtDistributionZNCzeroZNAanyHv2LowerSide);

  fDimuonPtDistributionZNCanyZNAzeroHv2LowerSide = new TH1F("fDimuonPtDistributionZNCanyZNAzeroHv2LowerSide", "fDimuonPtDistributionZNCanyZNAzeroHv2LowerSide", 4000, 0, 20);
  fOutputList->Add(fDimuonPtDistributionZNCanyZNAzeroHv2LowerSide);

  fDimuonPtDistributionZNCanyZNAanyHv2LowerSide = new TH1F("fDimuonPtDistributionZNCanyZNAanyHv2LowerSide", "fDimuonPtDistributionZNCanyZNAanyHv2LowerSide", 4000, 0, 20);
  fOutputList->Add(fDimuonPtDistributionZNCanyZNAanyHv2LowerSide);

  for( Int_t iRapidity = 0; iRapidity < 3; iRapidity++ ){
    fDimuonPtDistributionZNCzeroZNAzeroRapidityHv2LowerSide[iRapidity] = new TH1F(
              Form("fDimuonPtDistributionZNCzeroZNAzeroRapidityHv2LowerSide_%d", iRapidity),
              Form("fDimuonPtDistributionZNCzeroZNAzeroRapidityHv2LowerSide_%d", iRapidity),
              4000, 0, 20
              );
    fOutputList->Add(fDimuonPtDistributionZNCzeroZNAzeroRapidityHv2LowerSide[iRapidity]);
  }

  for( Int_t iRapidity = 0; iRapidity < 3; iRapidity++ ){
    fDimuonPtDistributionZNCzeroZNAanyRapidityHv2LowerSide[iRapidity] = new TH1F(
              Form("fDimuonPtDistributionZNCzeroZNAanyRapidityHv2LowerSide_%d", iRapidity),
              Form("fDimuonPtDistributionZNCzeroZNAanyRapidityHv2LowerSide_%d", iRapidity),
              4000, 0, 20
              );
    fOutputList->Add(fDimuonPtDistributionZNCzeroZNAanyRapidityHv2LowerSide[iRapidity]);
  }

  for( Int_t iRapidity = 0; iRapidity < 3; iRapidity++ ){
    fDimuonPtDistributionZNCanyZNAzeroRapidityHv2LowerSide[iRapidity] = new TH1F(
              Form("fDimuonPtDistributionZNCanyZNAzeroRapidityHv2LowerSide_%d", iRapidity),
              Form("fDimuonPtDistributionZNCanyZNAzeroRapidityHv2LowerSide_%d", iRapidity),
              4000, 0, 20
              );
    fOutputList->Add(fDimuonPtDistributionZNCanyZNAzeroRapidityHv2LowerSide[iRapidity]);
  }

  for( Int_t iRapidity = 0; iRapidity < 3; iRapidity++ ){
    fDimuonPtDistributionZNCanyZNAanyRapidityHv2LowerSide[iRapidity] = new TH1F(
              Form("fDimuonPtDistributionZNCanyZNAanyRapidityHv2LowerSide_%d", iRapidity),
              Form("fDimuonPtDistributionZNCanyZNAanyRapidityHv2LowerSide_%d", iRapidity),
              4000, 0, 20
              );
    fOutputList->Add(fDimuonPtDistributionZNCanyZNAanyRapidityHv2LowerSide[iRapidity]);
  }

  fDimuonPtDistributionZNCzeroZNAzeroHv3LowerSide = new TH1F("fDimuonPtDistributionZNCzeroZNAzeroHv3LowerSide", "fDimuonPtDistributionZNCzeroZNAzeroHv3LowerSide", PtBinNumber, PtBins);
  fOutputList->Add(fDimuonPtDistributionZNCzeroZNAzeroHv3LowerSide);

  fDimuonPtDistributionZNCzeroZNAanyHv3LowerSide = new TH1F("fDimuonPtDistributionZNCzeroZNAanyHv3LowerSide", "fDimuonPtDistributionZNCzeroZNAanyHv3LowerSide", PtBinNumber, PtBins);
  fOutputList->Add(fDimuonPtDistributionZNCzeroZNAanyHv3LowerSide);

  fDimuonPtDistributionZNCanyZNAzeroHv3LowerSide = new TH1F("fDimuonPtDistributionZNCanyZNAzeroHv3LowerSide", "fDimuonPtDistributionZNCanyZNAzeroHv3LowerSide", PtBinNumber, PtBins);
  fOutputList->Add(fDimuonPtDistributionZNCanyZNAzeroHv3LowerSide);

  fDimuonPtDistributionZNCanyZNAanyHv3LowerSide = new TH1F("fDimuonPtDistributionZNCanyZNAanyHv3LowerSide", "fDimuonPtDistributionZNCanyZNAanyHv3LowerSide", PtBinNumber, PtBins);
  fOutputList->Add(fDimuonPtDistributionZNCanyZNAanyHv3LowerSide);

  for( Int_t iRapidity = 0; iRapidity < 3; iRapidity++ ){
    fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3LowerSide[iRapidity] = new TH1F(
              Form("fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3LowerSide_%d", iRapidity),
              Form("fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3LowerSide_%d", iRapidity),
              PtBinNumber, PtBins
              );
    fOutputList->Add(fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3LowerSide[iRapidity]);
  }

  for( Int_t iRapidity = 0; iRapidity < 3; iRapidity++ ){
    fDimuonPtDistributionZNCzeroZNAanyRapidityHv3LowerSide[iRapidity] = new TH1F(
              Form("fDimuonPtDistributionZNCzeroZNAanyRapidityHv3LowerSide_%d", iRapidity),
              Form("fDimuonPtDistributionZNCzeroZNAanyRapidityHv3LowerSide_%d", iRapidity),
              PtBinNumber, PtBins
              );
    fOutputList->Add(fDimuonPtDistributionZNCzeroZNAanyRapidityHv3LowerSide[iRapidity]);
  }

  for( Int_t iRapidity = 0; iRapidity < 3; iRapidity++ ){
    fDimuonPtDistributionZNCanyZNAzeroRapidityHv3LowerSide[iRapidity] = new TH1F(
              Form("fDimuonPtDistributionZNCanyZNAzeroRapidityHv3LowerSide_%d", iRapidity),
              Form("fDimuonPtDistributionZNCanyZNAzeroRapidityHv3LowerSide_%d", iRapidity),
              PtBinNumber, PtBins
              );
    fOutputList->Add(fDimuonPtDistributionZNCanyZNAzeroRapidityHv3LowerSide[iRapidity]);
  }

  for( Int_t iRapidity = 0; iRapidity < 3; iRapidity++ ){
    fDimuonPtDistributionZNCanyZNAanyRapidityHv3LowerSide[iRapidity] = new TH1F(
              Form("fDimuonPtDistributionZNCanyZNAanyRapidityHv3LowerSide_%d", iRapidity),
              Form("fDimuonPtDistributionZNCanyZNAanyRapidityHv3LowerSide_%d", iRapidity),
              PtBinNumber, PtBins
              );
    fOutputList->Add(fDimuonPtDistributionZNCanyZNAanyRapidityHv3LowerSide[iRapidity]);
  }

  //
  // HIGHER SIDE
  //
  fDimuonPtDistributionZNCzeroZNAzeroHv2HigherSide = new TH1F("fDimuonPtDistributionZNCzeroZNAzeroHv2HigherSide", "fDimuonPtDistributionZNCzeroZNAzeroHv2HigherSide", 4000, 0, 20);
  fOutputList->Add(fDimuonPtDistributionZNCzeroZNAzeroHv2HigherSide);

  fDimuonPtDistributionZNCzeroZNAanyHv2HigherSide = new TH1F("fDimuonPtDistributionZNCzeroZNAanyHv2HigherSide", "fDimuonPtDistributionZNCzeroZNAanyHv2HigherSide", 4000, 0, 20);
  fOutputList->Add(fDimuonPtDistributionZNCzeroZNAanyHv2HigherSide);

  fDimuonPtDistributionZNCanyZNAzeroHv2HigherSide = new TH1F("fDimuonPtDistributionZNCanyZNAzeroHv2HigherSide", "fDimuonPtDistributionZNCanyZNAzeroHv2HigherSide", 4000, 0, 20);
  fOutputList->Add(fDimuonPtDistributionZNCanyZNAzeroHv2HigherSide);

  fDimuonPtDistributionZNCanyZNAanyHv2HigherSide = new TH1F("fDimuonPtDistributionZNCanyZNAanyHv2HigherSide", "fDimuonPtDistributionZNCanyZNAanyHv2HigherSide", 4000, 0, 20);
  fOutputList->Add(fDimuonPtDistributionZNCanyZNAanyHv2HigherSide);

  for( Int_t iRapidity = 0; iRapidity < 3; iRapidity++ ){
    fDimuonPtDistributionZNCzeroZNAzeroRapidityHv2HigherSide[iRapidity] = new TH1F(
              Form("fDimuonPtDistributionZNCzeroZNAzeroRapidityHv2HigherSide_%d", iRapidity),
              Form("fDimuonPtDistributionZNCzeroZNAzeroRapidityHv2HigherSide_%d", iRapidity),
              4000, 0, 20
              );
    fOutputList->Add(fDimuonPtDistributionZNCzeroZNAzeroRapidityHv2HigherSide[iRapidity]);
  }

  for( Int_t iRapidity = 0; iRapidity < 3; iRapidity++ ){
    fDimuonPtDistributionZNCzeroZNAanyRapidityHv2HigherSide[iRapidity] = new TH1F(
              Form("fDimuonPtDistributionZNCzeroZNAanyRapidityHv2HigherSide_%d", iRapidity),
              Form("fDimuonPtDistributionZNCzeroZNAanyRapidityHv2HigherSide_%d", iRapidity),
              4000, 0, 20
              );
    fOutputList->Add(fDimuonPtDistributionZNCzeroZNAanyRapidityHv2HigherSide[iRapidity]);
  }

  for( Int_t iRapidity = 0; iRapidity < 3; iRapidity++ ){
    fDimuonPtDistributionZNCanyZNAzeroRapidityHv2HigherSide[iRapidity] = new TH1F(
              Form("fDimuonPtDistributionZNCanyZNAzeroRapidityHv2HigherSide_%d", iRapidity),
              Form("fDimuonPtDistributionZNCanyZNAzeroRapidityHv2HigherSide_%d", iRapidity),
              4000, 0, 20
              );
    fOutputList->Add(fDimuonPtDistributionZNCanyZNAzeroRapidityHv2HigherSide[iRapidity]);
  }

  for( Int_t iRapidity = 0; iRapidity < 3; iRapidity++ ){
    fDimuonPtDistributionZNCanyZNAanyRapidityHv2HigherSide[iRapidity] = new TH1F(
              Form("fDimuonPtDistributionZNCanyZNAanyRapidityHv2HigherSide_%d", iRapidity),
              Form("fDimuonPtDistributionZNCanyZNAanyRapidityHv2HigherSide_%d", iRapidity),
              4000, 0, 20
              );
    fOutputList->Add(fDimuonPtDistributionZNCanyZNAanyRapidityHv2HigherSide[iRapidity]);
  }

  fDimuonPtDistributionZNCzeroZNAzeroHv3HigherSide = new TH1F("fDimuonPtDistributionZNCzeroZNAzeroHv3HigherSide", "fDimuonPtDistributionZNCzeroZNAzeroHv3HigherSide", PtBinNumber, PtBins);
  fOutputList->Add(fDimuonPtDistributionZNCzeroZNAzeroHv3HigherSide);

  fDimuonPtDistributionZNCzeroZNAanyHv3HigherSide = new TH1F("fDimuonPtDistributionZNCzeroZNAanyHv3HigherSide", "fDimuonPtDistributionZNCzeroZNAanyHv3HigherSide", PtBinNumber, PtBins);
  fOutputList->Add(fDimuonPtDistributionZNCzeroZNAanyHv3HigherSide);

  fDimuonPtDistributionZNCanyZNAzeroHv3HigherSide = new TH1F("fDimuonPtDistributionZNCanyZNAzeroHv3HigherSide", "fDimuonPtDistributionZNCanyZNAzeroHv3HigherSide", PtBinNumber, PtBins);
  fOutputList->Add(fDimuonPtDistributionZNCanyZNAzeroHv3HigherSide);

  fDimuonPtDistributionZNCanyZNAanyHv3HigherSide = new TH1F("fDimuonPtDistributionZNCanyZNAanyHv3HigherSide", "fDimuonPtDistributionZNCanyZNAanyHv3HigherSide", PtBinNumber, PtBins);
  fOutputList->Add(fDimuonPtDistributionZNCanyZNAanyHv3HigherSide);

  for( Int_t iRapidity = 0; iRapidity < 3; iRapidity++ ){
    fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3HigherSide[iRapidity] = new TH1F(
              Form("fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3HigherSide_%d", iRapidity),
              Form("fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3HigherSide_%d", iRapidity),
              PtBinNumber, PtBins
              );
    fOutputList->Add(fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3HigherSide[iRapidity]);
  }

  for( Int_t iRapidity = 0; iRapidity < 3; iRapidity++ ){
    fDimuonPtDistributionZNCzeroZNAanyRapidityHv3HigherSide[iRapidity] = new TH1F(
              Form("fDimuonPtDistributionZNCzeroZNAanyRapidityHv3HigherSide_%d", iRapidity),
              Form("fDimuonPtDistributionZNCzeroZNAanyRapidityHv3HigherSide_%d", iRapidity),
              PtBinNumber, PtBins
              );
    fOutputList->Add(fDimuonPtDistributionZNCzeroZNAanyRapidityHv3HigherSide[iRapidity]);
  }

  for( Int_t iRapidity = 0; iRapidity < 3; iRapidity++ ){
    fDimuonPtDistributionZNCanyZNAzeroRapidityHv3HigherSide[iRapidity] = new TH1F(
              Form("fDimuonPtDistributionZNCanyZNAzeroRapidityHv3HigherSide_%d", iRapidity),
              Form("fDimuonPtDistributionZNCanyZNAzeroRapidityHv3HigherSide_%d", iRapidity),
              PtBinNumber, PtBins
              );
    fOutputList->Add(fDimuonPtDistributionZNCanyZNAzeroRapidityHv3HigherSide[iRapidity]);
  }

  for( Int_t iRapidity = 0; iRapidity < 3; iRapidity++ ){
    fDimuonPtDistributionZNCanyZNAanyRapidityHv3HigherSide[iRapidity] = new TH1F(
              Form("fDimuonPtDistributionZNCanyZNAanyRapidityHv3HigherSide_%d", iRapidity),
              Form("fDimuonPtDistributionZNCanyZNAanyRapidityHv3HigherSide_%d", iRapidity),
              PtBinNumber, PtBins
              );
    fOutputList->Add(fDimuonPtDistributionZNCanyZNAanyRapidityHv3HigherSide[iRapidity]);
  }
  //
  // END SIDEBANDS
  //_______________________________

  //_______________________________
  // - End of the function
  PostData(1, fOutputList);           // postdata will notify the analysis manager of changes / updates to the
  // fOutputList object. the manager will in the end take care of writing your output to file
  // so it needs to know what's in the output
}
//_____________________________________________________________________________
void AliAnalysisTaskADin2018::NotifyRun()
{
  /// Set run number for cuts
  fMuonTrackCuts->SetRun(fInputHandler);
}
//_____________________________________________________________________________
void AliAnalysisTaskADin2018::UserExec(Option_t *)
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
  // fADcheck = 1;

  // get AOD event
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if(!fAOD) {
      PostData(1, fOutputList);
      return;
  }
  fCounterH->Fill(iSelectionCounter); // AOD event found
  iSelectionCounter++;

  /* - Is it the right trigger?
   * -
   * -
   */
  TString trigger = fAOD->GetFiredTriggerClasses();
  if ( !( trigger.Contains("CMUP6-B-NOPF-MUFAST") )  )  {
          PostData(1, fOutputList);
          return;
  }


  //_______________________________
  //
  // COHERENT check
  //
  // if (    !(trigger.Contains("CMUP11-B-NOPF-MUFAST") ||
  //         trigger.Contains("CMUP26-B-NOPF-MUFAST") ||
  //         trigger.Contains("CMUP6-B-NOPF-MUFAST")  ||
  //         trigger.Contains("CMUP10-B-NOPF-MUFAST") ||
  //         trigger.Contains("CMUP13-B-NOPF-MUFAST")  )
  //       )  {
  //                 PostData(1, fOutputList);
  //                 return;
  // }

  //_______________________________
  /* -
   * - Requiring both triggers to do all
   * - sort of checks.
   * -
   * - NB: CMUP11 = *VZEROA *ADA *ADC DIMUON_U
   * -
   */
  // if ( !( trigger.Contains("CMUP11-B-NOPF-MUFAST") )  )  {
  //         PostData(1, fOutputList);
  //         return;
  // }
  //_______________________________
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
  if ( trigger.Contains("CMUP6-B-NOPF-MUFAST") )  {
    fTriggersVsRunH->Fill( 2.5, fRunNum );
    fRunNumberTriggerCMUP6ClassH        ->Fill(fRunNum);
    fRunNumberTriggerCMUP6ClassProperlyH->Fill( Form("%d", fRunNum) , 1 );
    // counterForTrigger++;
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
  // if ( fRunNum <= 245068 ) {
  //   calibrated = 1;
  // } else if ( ( fRunNum > 245068 ) && ( fRunNum <  246995 ) ){
  //   calibrated = 0;
  // } else {
  //   calibrated = 1;
  // }

  // if ( calibrated == 0 ) {
  //   if( fRunNum <= 246994 ) {
  //     fZNAEnergy *= (2500./250.);
  //     fZNCEnergy *= (2500./250.);
  //   }
  // }

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

  Int_t is_VZEROA_set = -9;
  Int_t is_VZEROC_set = -9;
  is_VZEROA_set = IntBits( dataVZERO->GetTriggerBits() ).test(12);
  is_VZEROC_set = IntBits( dataVZERO->GetTriggerBits() ).test(13);


  Int_t VZEROAPastFutureBeamBeamFlags[32][21];
  Int_t VZEROCPastFutureBeamBeamFlags[32][21];
  for(   Int_t iChannel = 0; iChannel < 32; iChannel++ ){
    for( Int_t iClock   = 0; iClock   < 21; iClock++   ){
      VZEROAPastFutureBeamBeamFlags[iChannel][iClock] = 0;
      VZEROCPastFutureBeamBeamFlags[iChannel][iClock] = 0;
    }
  }

  for(   Int_t iChannel = 0; iChannel < 32; iChannel++ ){
    for( Int_t iClock   = 0; iClock   < 21; iClock++   ){
      VZEROAPastFutureBeamBeamFlags[iChannel][iClock] = dataVZERO->GetPFBBFlag(iChannel + 32, iClock);
      VZEROCPastFutureBeamBeamFlags[iChannel][iClock] = dataVZERO->GetPFBBFlag(iChannel, iClock);
    }
  }

  Int_t VZEROAPastFutureBoolean = 0;
  Int_t VZEROCPastFutureBoolean = 0;
  for(   Int_t iChannel = 0; iChannel < 32; iChannel++ ){
    for( Int_t iClock   = 0; iClock   < 21; iClock++   ){
      if( dataVZERO->GetPFBBFlag(iChannel + 32, iClock) != 0 ) VZEROAPastFutureBoolean = 1;
      if( dataVZERO->GetPFBBFlag(iChannel, iClock)      != 0 ) VZEROCPastFutureBoolean = 1;
    }
  }

  Double_t VZEROmultiplicities[64]   = { -1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1,
                                         -1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1,
                                         -1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1,
                                         -1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1 };
  Double_t VZEROmultiplicitiesTotal  = 0;
  Double_t VZEROAmultiplicitiesTotal = 0;
  Double_t VZEROCmultiplicitiesTotal = 0;

  for( Int_t iChannel = 0; iChannel < 64; iChannel++ ){
    VZEROmultiplicities[iChannel] = dataVZERO->GetMultiplicity(iChannel);
    VZEROmultiplicitiesTotal     += dataVZERO->GetMultiplicity(iChannel);
    if ( iChannel < 32 ) {
      VZEROCmultiplicitiesTotal  += dataVZERO->GetMultiplicity(iChannel);
    } else {
      VZEROAmultiplicitiesTotal  += dataVZERO->GetMultiplicity(iChannel);
    }
  }


  //_____________________________________
  // RUN SELECTION
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
  Int_t is_ADA_set = -9;
  Int_t is_ADC_set = -9;
  Double_t ADmultiplicities[16]   = { -1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1 };
  Double_t ADmultiplicitiesTotal  = 0;
  Double_t ADAmultiplicitiesTotal = 0;
  Double_t ADCmultiplicitiesTotal = 0;

  Int_t ADAPastFutureBeamBeamFlags[8][21];
  Int_t ADCPastFutureBeamBeamFlags[8][21];


  Int_t ADAPastFutureBoolean = 0;
  Int_t ADCPastFutureBoolean = 0;

  if(dataAD) {
        fCounterH->Fill(iSelectionCounter);
        iSelectionCounter++;
        fCounterH->Fill(20);

        fADADecision = dataAD->GetADADecision();
        fADCDecision = dataAD->GetADCDecision();
        fCounterH->Fill(21);

        is_ADA_set = IntBits( dataAD->GetTriggerBits() ).test(12);
        is_ADC_set = IntBits( dataAD->GetTriggerBits() ).test(13);
        // cout << "is_ADA_set = " << is_ADA_set << endl;
        // cout << "is_ADC_set = " << is_ADC_set << endl;
        // cout << "is_ADA_set = " << IntBits( dataAD->GetTriggerBits() ) << endl;
        // cout << "is_ADC_set = " << dataAD->GetTriggerBits() << endl;
        for( Int_t iChannel = 0; iChannel < 16; iChannel++ ){
          ADmultiplicities[iChannel] = dataAD->GetMultiplicity(iChannel);
          ADmultiplicitiesTotal     += dataAD->GetMultiplicity(iChannel);
          if ( iChannel < 8 ) {
            ADCmultiplicitiesTotal  += dataAD->GetMultiplicity(iChannel);
          } else {
            ADAmultiplicitiesTotal  += dataAD->GetMultiplicity(iChannel);
          }
        }


        for(   Int_t iChannel = 0; iChannel < 8; iChannel++ ){
          for( Int_t iClock   = 0; iClock   < 21; iClock++   ){
            ADAPastFutureBeamBeamFlags[iChannel][iClock] = 0;
            ADCPastFutureBeamBeamFlags[iChannel][iClock] = 0;
          }
        }

        for(   Int_t iChannel = 0; iChannel < 8; iChannel++ ){
          for( Int_t iClock   = 0; iClock   < 21; iClock++   ){
            ADAPastFutureBeamBeamFlags[iChannel][iClock] = dataAD->GetPFBBFlag(iChannel + 8, iClock);
            ADCPastFutureBeamBeamFlags[iChannel][iClock] = dataAD->GetPFBBFlag(iChannel, iClock);
          }
        }

        for(   Int_t iChannel = 0; iChannel < 8;  iChannel++ ){
          for( Int_t iClock   = 0; iClock   < 21; iClock++   ){
            if( dataAD->GetPFBBFlag(iChannel + 8, iClock) != 0 ) ADAPastFutureBoolean = 1;
            if( dataAD->GetPFBBFlag(iChannel, iClock)     != 0 ) ADCPastFutureBoolean = 1;
          }
        }


  }
  fCounterH->Fill(22);

  // if( ADAPastFutureBoolean == 0 ) cout << "ADAPastFutureBoolean" << endl;

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
  /**
   * Check with no AD veto at all...
   */
  if(fADADecision != 0) {
       PostData(1, fOutputList);
       return;
  }
  if(fADCDecision != 0) {
       PostData(1, fOutputList);
       return;
  }

  //
  // Int_t fADactivity = -9;
  // if ( fADADecision == 1) fADactivity = 1;
  // if ( fADCDecision == 1) fADactivity = 1;
  // if ( fADactivity != 1) {
  //     /* -
  //      * - Remember that this is akin
  //      * - to asking Beam-Beam
  //      * - activities.
  //      */
  //     PostData(1, fOutputList);
  //     return;
  // }


  /**
   * - This is the AD check.
   * - This is selected at the
   * - level of train set up...
   * - Needed for both the neutron
   * - emission analysis and the
   * - p-Pb analysis!!
   * -
   */
  // if( fADcheck != 0){
      // if(fADCDecision != 0) {
      //     // fEntriesAgainstRunNumberProperlyH->Fill( Form("%d", fRunNum) , 1 );
      //     PostData(1, fOutputList);
      //     return;
      // }
  // }
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
  // if( fV0TotalNCells < 3 ) {
  //      PostData(1, fOutputList);
  //      return;
  // }



  //
  //
  //
  // /******************************
  //  *       HW + SW check        *
  //  ******************************/
  // /**
  //  * Check with no AD veto at all...
  //  */
  // if ( is_ADA_set == 0 ){
  //     /* -
  //      * - This means that we don't have
  //      * - the hardware trigger.
  //      * -
  //      */
  //     fEntriesAgainstRunNumberProperlyH          ->Fill( "ADA0" , 1 );
  //     if ( is_ADC_set == 0 ){
  //         fEntriesAgainstRunNumberProperlyH      ->Fill( "ADA0_ADC0" , 1 );
  //         if (fADADecision != 0) {
  //             fEntriesAgainstRunNumberProperlyH  ->Fill( "ADA0_ADC0_Adec1" , 1 );
  //             if(fADCDecision != 0) {
  //               fEntriesAgainstRunNumberProperlyH->Fill( "ADA0_ADC0_Adec1_Cdec1" , 1 );
  //             } else {
  //               fEntriesAgainstRunNumberProperlyH->Fill( "ADA0_ADC0_Adec1_Cdec0" , 1 );
  //             }
  //         } else  {
  //             fEntriesAgainstRunNumberProperlyH  ->Fill( "ADA0_ADC0_Adec0" , 1 );
  //             if(fADCDecision != 0) {
  //               fEntriesAgainstRunNumberProperlyH->Fill( "ADA0_ADC0_Adec0_Cdec1" , 1 );
  //             } else {
  //               fEntriesAgainstRunNumberProperlyH->Fill( "ADA0_ADC0_Adec0_Cdec0" , 1 );
  //             }
  //         }
  //     } else {
  //         fEntriesAgainstRunNumberProperlyH      ->Fill( "ADA0_ADC1" , 1 );
  //         if (fADADecision != 0) {
  //             fEntriesAgainstRunNumberProperlyH  ->Fill( "ADA0_ADC1_Adec1" , 1 );
  //             if(fADCDecision != 0) {
  //               fEntriesAgainstRunNumberProperlyH->Fill( "ADA0_ADC1_Adec1_Cdec1" , 1 );
  //             } else {
  //               fEntriesAgainstRunNumberProperlyH->Fill( "ADA0_ADC1_Adec1_Cdec0" , 1 );
  //             }
  //         } else  {
  //             fEntriesAgainstRunNumberProperlyH  ->Fill( "ADA0_ADC1_Adec0" , 1 );
  //             if(fADCDecision != 0) {
  //               fEntriesAgainstRunNumberProperlyH->Fill( "ADA0_ADC1_Adec0_Cdec1" , 1 );
  //             } else {
  //               fEntriesAgainstRunNumberProperlyH->Fill( "ADA0_ADC1_Adec0_Cdec0" , 1 );
  //             }
  //         }
  //     }
  // } else {
  //     /* -
  //      * - This means that we have
  //      * - the hardware trigger.
  //      * -
  //      */
  //     fEntriesAgainstRunNumberProperlyH          ->Fill( "ADA1" , 1 );
  //     if ( is_ADC_set == 0 ){
  //         fEntriesAgainstRunNumberProperlyH      ->Fill( "ADA1_ADC0" , 1 );
  //         if (fADADecision != 0) {
  //             fEntriesAgainstRunNumberProperlyH  ->Fill( "ADA1_ADC0_Adec1" , 1 );
  //             if(fADCDecision != 0) {
  //               fEntriesAgainstRunNumberProperlyH->Fill( "ADA1_ADC0_Adec1_Cdec1" , 1 );
  //             } else {
  //               fEntriesAgainstRunNumberProperlyH->Fill( "ADA1_ADC0_Adec1_Cdec0" , 1 );
  //             }
  //         } else  {
  //             fEntriesAgainstRunNumberProperlyH  ->Fill( "ADA1_ADC0_Adec0" , 1 );
  //             if(fADCDecision != 0) {
  //               fEntriesAgainstRunNumberProperlyH->Fill( "ADA1_ADC0_Adec0_Cdec1" , 1 );
  //             } else {
  //               fEntriesAgainstRunNumberProperlyH->Fill( "ADA1_ADC0_Adec0_Cdec0" , 1 );
  //             }
  //         }
  //     } else {
  //         fEntriesAgainstRunNumberProperlyH      ->Fill( "ADA1_ADC1" , 1 );
  //         if (fADADecision != 0) {
  //             fEntriesAgainstRunNumberProperlyH  ->Fill( "ADA1_ADC1_Adec1" , 1 );
  //             if(fADCDecision != 0) {
  //               fEntriesAgainstRunNumberProperlyH->Fill( "ADA1_ADC1_Adec1_Cdec1" , 1 );
  //             } else {
  //               fEntriesAgainstRunNumberProperlyH->Fill( "ADA1_ADC1_Adec1_Cdec0" , 1 );
  //             }
  //         } else  {
  //             fEntriesAgainstRunNumberProperlyH  ->Fill( "ADA1_ADC1_Adec0" , 1 );
  //             if(fADCDecision != 0) {
  //               fEntriesAgainstRunNumberProperlyH->Fill( "ADA1_ADC1_Adec0_Cdec1" , 1 );
  //             } else {
  //               fEntriesAgainstRunNumberProperlyH->Fill( "ADA1_ADC1_Adec0_Cdec0" , 1 );
  //             }
  //         }
  //     }
  // }
  //
  //
  //
  // /**
  //  * Check with no AD veto at all...
  //  */
  // if(fADADecision != 0) {
  //      PostData(1, fOutputList);
  //      return;
  // }
  //
  // /**
  //  * - This is the AD check.
  //  * - This is selected at the
  //  * - level of train set up...
  //  * - Needed for both the neutron
  //  * - emission analysis and the
  //  * - p-Pb analysis!!
  //  * -
  //  */
  // if( fADcheck != 0){
  //     if(fADCDecision != 0) {
  //         // fEntriesAgainstRunNumberProperlyH->Fill( Form("%d", fRunNum) , 1 );
  //         PostData(1, fOutputList);
  //         return;
  //     }
  // }






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
        fEtaMuonH ->Fill(track[iFilling]->Eta());
        fRAbsMuonH->Fill(track[iFilling]->GetRAtAbsorberEnd());
  }
  // store muons
  fNumberMuonsH->Fill(nGoodMuons);
  fEntriesAgainstRunNumberH->Fill(fRunNum);
  /* - This is the last part of my try to obtain a proper RunNumbers histogram...
     -
   */
  // fEntriesAgainstRunNumberProperlyH->Fill( Form("%d", fRunNum) , 1 );
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
  if ( (possibleJPsi.Mag() > 2.85) && (possibleJPsi.Mag() < 3.35) ) {
    fDimuonPtDistributionH            ->Fill(ptOfTheDimuonPair);
    fDimuonPtDistributionShiftPlusOneH->Fill(ptOfTheDimuonPair);
    if (        possibleJPsi.Rapidity() > -4.0  && possibleJPsi.Rapidity() <= -3.75 ) {
      fDimuonPtDistributionRapidityH[0]->Fill(ptOfTheDimuonPair);
      if (        ptOfTheDimuonPair < 0.275 ) {
        fDimuonPtDistributionRapidityHv3[0]->Fill(ptOfTheDimuonPair);
      } else if ( ptOfTheDimuonPair < 0.950 ) {
        fDimuonPtDistributionRapidityHv3[0]->Fill( ptOfTheDimuonPair, 0.33333333333 );
      } else if ( ptOfTheDimuonPair < 1.400 ) {
        fDimuonPtDistributionRapidityHv3[0]->Fill( ptOfTheDimuonPair, 0.16666666666 );
      } else if ( ptOfTheDimuonPair < 2.000 ) {
        fDimuonPtDistributionRapidityHv3[0]->Fill( ptOfTheDimuonPair, 0.125 );
      } else if ( ptOfTheDimuonPair < 4.000 ) {
        fDimuonPtDistributionRapidityHv3[0]->Fill( ptOfTheDimuonPair, 0.050 );
      } else if ( ptOfTheDimuonPair < 5.000 ) {
        fDimuonPtDistributionRapidityHv3[0]->Fill( ptOfTheDimuonPair, 0.025 );
      }
    } else if ( possibleJPsi.Rapidity() > -3.75 && possibleJPsi.Rapidity() <= -3.50 ) {
      fDimuonPtDistributionRapidityH[1]->Fill(ptOfTheDimuonPair);
      if (        ptOfTheDimuonPair < 0.275 ) {
        fDimuonPtDistributionRapidityHv3[1]->Fill(ptOfTheDimuonPair);
      } else if ( ptOfTheDimuonPair < 0.950 ) {
        fDimuonPtDistributionRapidityHv3[1]->Fill( ptOfTheDimuonPair, 0.33333333333 );
      } else if ( ptOfTheDimuonPair < 1.400 ) {
        fDimuonPtDistributionRapidityHv3[1]->Fill( ptOfTheDimuonPair, 0.16666666666 );
      } else if ( ptOfTheDimuonPair < 2.000 ) {
        fDimuonPtDistributionRapidityHv3[1]->Fill( ptOfTheDimuonPair, 0.125 );
      } else if ( ptOfTheDimuonPair < 4.000 ) {
        fDimuonPtDistributionRapidityHv3[1]->Fill( ptOfTheDimuonPair, 0.050 );
      } else if ( ptOfTheDimuonPair < 5.000 ) {
        fDimuonPtDistributionRapidityHv3[1]->Fill( ptOfTheDimuonPair, 0.025 );
      }
    } else if ( possibleJPsi.Rapidity() > -3.50 && possibleJPsi.Rapidity() <= -3.25 ) {
      fDimuonPtDistributionRapidityH[2]->Fill(ptOfTheDimuonPair);
      if (        ptOfTheDimuonPair < 0.275 ) {
        fDimuonPtDistributionRapidityHv3[2]->Fill(ptOfTheDimuonPair);
      } else if ( ptOfTheDimuonPair < 0.950 ) {
        fDimuonPtDistributionRapidityHv3[2]->Fill( ptOfTheDimuonPair, 0.33333333333 );
      } else if ( ptOfTheDimuonPair < 1.400 ) {
        fDimuonPtDistributionRapidityHv3[2]->Fill( ptOfTheDimuonPair, 0.16666666666 );
      } else if ( ptOfTheDimuonPair < 2.000 ) {
        fDimuonPtDistributionRapidityHv3[2]->Fill( ptOfTheDimuonPair, 0.125 );
      } else if ( ptOfTheDimuonPair < 4.000 ) {
        fDimuonPtDistributionRapidityHv3[2]->Fill( ptOfTheDimuonPair, 0.050 );
      } else if ( ptOfTheDimuonPair < 5.000 ) {
        fDimuonPtDistributionRapidityHv3[2]->Fill( ptOfTheDimuonPair, 0.025 );
      }
    } else if ( possibleJPsi.Rapidity() > -3.25 && possibleJPsi.Rapidity() <= -3.00 ) {
      fDimuonPtDistributionRapidityH[3]->Fill(ptOfTheDimuonPair);
      if (        ptOfTheDimuonPair < 0.275 ) {
        fDimuonPtDistributionRapidityHv3[3]->Fill(ptOfTheDimuonPair);
      } else if ( ptOfTheDimuonPair < 0.950 ) {
        fDimuonPtDistributionRapidityHv3[3]->Fill( ptOfTheDimuonPair, 0.33333333333 );
      } else if ( ptOfTheDimuonPair < 1.400 ) {
        fDimuonPtDistributionRapidityHv3[3]->Fill( ptOfTheDimuonPair, 0.16666666666 );
      } else if ( ptOfTheDimuonPair < 2.000 ) {
        fDimuonPtDistributionRapidityHv3[3]->Fill( ptOfTheDimuonPair, 0.125 );
      } else if ( ptOfTheDimuonPair < 4.000 ) {
        fDimuonPtDistributionRapidityHv3[3]->Fill( ptOfTheDimuonPair, 0.050 );
      } else if ( ptOfTheDimuonPair < 5.000 ) {
        fDimuonPtDistributionRapidityHv3[3]->Fill( ptOfTheDimuonPair, 0.025 );
      }
    } else if ( possibleJPsi.Rapidity() > -3.00 && possibleJPsi.Rapidity() <= -2.75 ) {
      fDimuonPtDistributionRapidityH[4]->Fill(ptOfTheDimuonPair);
      if (        ptOfTheDimuonPair < 0.275 ) {
        fDimuonPtDistributionRapidityHv3[4]->Fill(ptOfTheDimuonPair);
      } else if ( ptOfTheDimuonPair < 0.950 ) {
        fDimuonPtDistributionRapidityHv3[4]->Fill( ptOfTheDimuonPair, 0.33333333333 );
      } else if ( ptOfTheDimuonPair < 1.400 ) {
        fDimuonPtDistributionRapidityHv3[4]->Fill( ptOfTheDimuonPair, 0.16666666666 );
      } else if ( ptOfTheDimuonPair < 2.000 ) {
        fDimuonPtDistributionRapidityHv3[4]->Fill( ptOfTheDimuonPair, 0.125 );
      } else if ( ptOfTheDimuonPair < 4.000 ) {
        fDimuonPtDistributionRapidityHv3[4]->Fill( ptOfTheDimuonPair, 0.050 );
      } else if ( ptOfTheDimuonPair < 5.000 ) {
        fDimuonPtDistributionRapidityHv3[4]->Fill( ptOfTheDimuonPair, 0.025 );
      }
    } else if ( possibleJPsi.Rapidity() > -2.75 && possibleJPsi.Rapidity() <= -2.50 ) {
      fDimuonPtDistributionRapidityH[5]->Fill(ptOfTheDimuonPair);
      if (        ptOfTheDimuonPair < 0.275 ) {
        fDimuonPtDistributionRapidityHv3[5]->Fill(ptOfTheDimuonPair);
      } else if ( ptOfTheDimuonPair < 0.950 ) {
        fDimuonPtDistributionRapidityHv3[5]->Fill( ptOfTheDimuonPair, 0.33333333333 );
      } else if ( ptOfTheDimuonPair < 1.400 ) {
        fDimuonPtDistributionRapidityHv3[5]->Fill( ptOfTheDimuonPair, 0.16666666666 );
      } else if ( ptOfTheDimuonPair < 2.000 ) {
        fDimuonPtDistributionRapidityHv3[5]->Fill( ptOfTheDimuonPair, 0.125 );
      } else if ( ptOfTheDimuonPair < 4.000 ) {
        fDimuonPtDistributionRapidityHv3[5]->Fill( ptOfTheDimuonPair, 0.050 );
      } else if ( ptOfTheDimuonPair < 5.000 ) {
        fDimuonPtDistributionRapidityHv3[5]->Fill( ptOfTheDimuonPair, 0.025 );
      }
    }
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
    if ( (possibleJPsi.Mag() > 2.85) && (possibleJPsi.Mag() < 3.35) && (possibleJPsi.Pt() < 0.25)){
      fZNCEnergyAgainstEntriesH        ->Fill(fZNCEnergy);
    }
    fZNCEnergyAgainstEntriesExtendedH->Fill(fZNCEnergy);
    if ( calibrated == 0 ) fZNCEnergyUncalibratedH->Fill(fZNCEnergy);
    if ( calibrated == 1 ) {
      fZNCEnergyCalibratedH          ->Fill( fZNCEnergy );
      fZNCEnergyCalibratedHigherGainH->Fill( dataZDC->GetZNCTowerEnergyLR()[0] );
    }
  }
  if ( (possibleJPsi.Mag() > 2.85) && (possibleJPsi.Mag() < 3.35) && (possibleJPsi.Pt() < 0.25)){
    fZNCEnergyBeforeTimingSelectionH               ->Fill(fZNCEnergy);
    if ( ADCPastFutureBoolean    == 0 ) fADCmultiplicityTotalVsZNCenergyH   ->Fill(ADCmultiplicitiesTotal,    fZNCEnergy);
    if ( VZEROCPastFutureBoolean == 0 ) fVZEROCmultiplicityTotalVsZNCenergyH->Fill(VZEROCmultiplicitiesTotal, fZNCEnergy);
    if ( is_ADC_set == 0 )    {
      if ( ADCPastFutureBoolean == 0 ) fADCmultiplicityTotalVsZNCenergyH_ADCno      ->Fill(ADCmultiplicitiesTotal,    fZNCEnergy);
    }
    if ( is_VZEROC_set == 0 ) {
      if ( VZEROCPastFutureBoolean == 0 ) fVZEROCmultiplicityTotalVsZNCenergyH_VZEROCno->Fill(VZEROCmultiplicitiesTotal, fZNCEnergy);
    }

  }
  fZNCEnergyBeforeTimingSelectionExtendedH->Fill(fZNCEnergy);
  if ( dataZDC->IsZNCfired() && ( isZNCfired != 0 ) ) {
    fZNCEnergyAgainstEntriesExtendedHv2->Fill(fZNCEnergy);
  }
  if ( dataZDC->IsZNAfired() && ( isZNAfired != 0 ) ) {
    fZNAEnergyAgainstEntriesExtendedHv2->Fill(fZNAEnergy);
  }
  if ( isZNAfired != 0 ) {
    if ( (possibleJPsi.Mag() > 2.85) && (possibleJPsi.Mag() < 3.35) && (possibleJPsi.Pt() < 0.25)){
      fZNAEnergyAgainstEntriesH        ->Fill(fZNAEnergy);
    }
    fZNAEnergyAgainstEntriesExtendedH->Fill(fZNAEnergy);
    if ( calibrated == 0 ) fZNAEnergyUncalibratedH->Fill(fZNAEnergy);
    if ( calibrated == 1 ) {
      fZNAEnergyCalibratedH          ->Fill( fZNAEnergy );
      fZNAEnergyCalibratedHigherGainH->Fill( dataZDC->GetZNATowerEnergyLR()[0] );
    }
  }
  if ( (possibleJPsi.Mag() > 2.85) && (possibleJPsi.Mag() < 3.35) && (possibleJPsi.Pt() < 0.25)){
    fZNAEnergyBeforeTimingSelectionH               ->Fill(fZNAEnergy);
    if ( ADAPastFutureBoolean    == 0 ) fADAmultiplicityTotalVsZNAenergyH   ->Fill(ADAmultiplicitiesTotal,    fZNAEnergy);
    if ( VZEROAPastFutureBoolean == 0 ) fVZEROAmultiplicityTotalVsZNAenergyH->Fill(VZEROAmultiplicitiesTotal, fZNAEnergy);
    if ( is_ADA_set == 0 )    {
      if ( ADAPastFutureBoolean == 0 ) fADAmultiplicityTotalVsZNAenergyH_ADAno      ->Fill(ADAmultiplicitiesTotal,    fZNAEnergy);
    }
    if ( is_VZEROA_set == 0 ) {
      if ( VZEROAPastFutureBoolean == 0 ) fVZEROAmultiplicityTotalVsZNAenergyH_VZEROAno->Fill(VZEROAmultiplicitiesTotal, fZNAEnergy);
    }
  }
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
  for( Int_t iChannel = 0; iChannel < 16; iChannel++ ){
    fADmultiplicityH[iChannel]->Fill(ADmultiplicities[iChannel]);
  }


  if ( ADAPastFutureBoolean == 0 && ADCPastFutureBoolean == 0 && VZEROAPastFutureBoolean == 0 && VZEROCPastFutureBoolean == 0 ){
  fADmultiplicityTotalH ->Fill(ADmultiplicitiesTotal);
  fADAmultiplicityTotalH->Fill(ADAmultiplicitiesTotal);
  fADCmultiplicityTotalH->Fill(ADCmultiplicitiesTotal);

  fVZEROmultiplicityTotalH ->Fill(VZEROmultiplicitiesTotal);
  fVZEROAmultiplicityTotalH->Fill(VZEROAmultiplicitiesTotal);
  fVZEROCmultiplicityTotalH->Fill(VZEROCmultiplicitiesTotal);

  fADAmultiplicityVsVZEROAmultiplicityH->Fill(ADAmultiplicitiesTotal, VZEROAmultiplicitiesTotal);
  fADCmultiplicityVsVZEROCmultiplicityH->Fill(ADCmultiplicitiesTotal, VZEROCmultiplicitiesTotal);

  if ( is_VZEROA_set == 0 ) {
    fADAmultiplicityVsVZEROAmultiplicityH_VZEROAno        ->Fill(ADAmultiplicitiesTotal, VZEROAmultiplicitiesTotal);
    if ( is_ADA_set == 0 )  {
      fADAmultiplicityVsVZEROAmultiplicityH_ADAno_VZEROAno->Fill(ADAmultiplicitiesTotal, VZEROAmultiplicitiesTotal);
    }
  }
  if ( is_VZEROC_set == 0 ) {
    fADCmultiplicityVsVZEROCmultiplicityH_VZEROCno        ->Fill(ADCmultiplicitiesTotal, VZEROCmultiplicitiesTotal);
    if ( is_ADC_set == 0 )  {
      fADCmultiplicityVsVZEROCmultiplicityH_ADCno_VZEROCno->Fill(ADCmultiplicitiesTotal, VZEROCmultiplicitiesTotal);
    }
  }
  }


  /* - Filling the v2 histogram only if the
   * - ZNC or the ZNA have detected any activity at all...
   */
  if( isZNCfired == 0 ) {
        if( isZNAfired == 0 ) {
              //_______________________________
              /* -
               * - CREATING SIDEBANDS
               * - (LOWER side)
               * - for templates
               * -
               */
              if ( (possibleJPsi.Mag() > 2.4) && (possibleJPsi.Mag() < 2.8) ) {
                  fDimuonPtDistributionZNCzeroZNAzeroHv2LowerSide             ->Fill(ptOfTheDimuonPair);
                  /* -
                   * - Variable pt-binning.
                   * -
                   */
                  // if (        ptOfTheDimuonPair < 0.500 ) {
                  if (        ptOfTheDimuonPair < 0.275 ) {
                    fDimuonPtDistributionZNCzeroZNAzeroHv3LowerSide           ->Fill(ptOfTheDimuonPair);
                  } else if ( ptOfTheDimuonPair < 0.950 ) {
                    fDimuonPtDistributionZNCzeroZNAzeroHv3LowerSide           ->Fill( ptOfTheDimuonPair, 0.33333333333 );
                  } else if ( ptOfTheDimuonPair < 1.400 ) {
                    fDimuonPtDistributionZNCzeroZNAzeroHv3LowerSide           ->Fill( ptOfTheDimuonPair, 0.16666666666 );
                  } else if ( ptOfTheDimuonPair < 2.000 ) {
                    fDimuonPtDistributionZNCzeroZNAzeroHv3LowerSide           ->Fill( ptOfTheDimuonPair, 0.125 );
                  } else if ( ptOfTheDimuonPair < 4.000 ) {
                    fDimuonPtDistributionZNCzeroZNAzeroHv3LowerSide           ->Fill( ptOfTheDimuonPair, 0.050 );
                  } else if ( ptOfTheDimuonPair < 5.000 ) {
                    fDimuonPtDistributionZNCzeroZNAzeroHv3LowerSide           ->Fill( ptOfTheDimuonPair, 0.025 );
                  }
                  // fDimuonPtDistributionZNCzeroZNAzeroShiftPlusOneHv2->Fill(ptOfTheDimuonPair);
                  if (        possibleJPsi.Rapidity() > -4.0  && possibleJPsi.Rapidity() <= -3.50 ) {
                    fDimuonPtDistributionZNCzeroZNAzeroRapidityHv2LowerSide[0]->Fill(ptOfTheDimuonPair);
                    /* -
                     * - Variable pt-binning.
                     * -
                     */
                    // if (        ptOfTheDimuonPair < 0.500 ) {
                    if (        ptOfTheDimuonPair < 0.275 ) {
                      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3LowerSide[0]->Fill(ptOfTheDimuonPair);
                    } else if ( ptOfTheDimuonPair < 0.950 ) {
                      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3LowerSide[0]->Fill(ptOfTheDimuonPair, 0.33333333333);
                    } else if ( ptOfTheDimuonPair < 1.400 ) {
                      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3LowerSide[0]->Fill(ptOfTheDimuonPair, 0.16666666666);
                    } else if ( ptOfTheDimuonPair < 2.000 ) {
                      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3LowerSide[0]->Fill(ptOfTheDimuonPair, 0.125);
                    } else if ( ptOfTheDimuonPair < 4.000 ) {
                      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3LowerSide[0]->Fill(ptOfTheDimuonPair, 0.050);
                    } else if ( ptOfTheDimuonPair < 5.000 ) {
                      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3LowerSide[0]->Fill(ptOfTheDimuonPair, 0.025);
                    }
                  } else if ( possibleJPsi.Rapidity() > -3.50 && possibleJPsi.Rapidity() <= -3.00 ) {
                    fDimuonPtDistributionZNCzeroZNAzeroRapidityHv2LowerSide[1]->Fill(ptOfTheDimuonPair);
                    /* -
                     * - Variable pt-binning.
                     * -
                     */
                    // if (        ptOfTheDimuonPair < 0.500 ) {
                    if (        ptOfTheDimuonPair < 0.275 ) {
                      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3LowerSide[1]->Fill(ptOfTheDimuonPair);
                    } else if ( ptOfTheDimuonPair < 0.950 ) {
                      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3LowerSide[1]->Fill(ptOfTheDimuonPair, 0.33333333333);
                    } else if ( ptOfTheDimuonPair < 1.400 ) {
                      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3LowerSide[1]->Fill(ptOfTheDimuonPair, 0.16666666666);
                    } else if ( ptOfTheDimuonPair < 2.000 ) {
                      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3LowerSide[1]->Fill(ptOfTheDimuonPair, 0.125);
                    } else if ( ptOfTheDimuonPair < 4.000 ) {
                      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3LowerSide[1]->Fill(ptOfTheDimuonPair, 0.050);
                    } else if ( ptOfTheDimuonPair < 5.000 ) {
                      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3LowerSide[1]->Fill(ptOfTheDimuonPair, 0.025);
                    }
                  } else if ( possibleJPsi.Rapidity() > -3.00 && possibleJPsi.Rapidity() <= -2.50 ) {
                    fDimuonPtDistributionZNCzeroZNAzeroRapidityHv2LowerSide[2]->Fill(ptOfTheDimuonPair);
                    /* -
                     * - Variable pt-binning.
                     * -
                     */
                    // if (        ptOfTheDimuonPair < 0.500 ) {
                    if (        ptOfTheDimuonPair < 0.275 ) {
                      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3LowerSide[2]->Fill(ptOfTheDimuonPair);
                    } else if ( ptOfTheDimuonPair < 0.950 ) {
                      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3LowerSide[2]->Fill(ptOfTheDimuonPair, 0.33333333333);
                    } else if ( ptOfTheDimuonPair < 1.400 ) {
                      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3LowerSide[2]->Fill(ptOfTheDimuonPair, 0.16666666666);
                    } else if ( ptOfTheDimuonPair < 2.000 ) {
                      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3LowerSide[2]->Fill(ptOfTheDimuonPair, 0.125);
                    } else if ( ptOfTheDimuonPair < 4.000 ) {
                      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3LowerSide[2]->Fill(ptOfTheDimuonPair, 0.050);
                    } else if ( ptOfTheDimuonPair < 5.000 ) {
                      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3LowerSide[2]->Fill(ptOfTheDimuonPair, 0.025);
                    }
                  }
              }
              //
              // END SIDEBANDS
              //_______________________________
              /* -
               * - CREATING SIDEBANDS
               * - (HIGHER side)
               * - for templates
               * -
               */
              if ( (possibleJPsi.Mag() > 4.0) && (possibleJPsi.Mag() < 5.5) ) {
                  fDimuonPtDistributionZNCzeroZNAzeroHv2HigherSide             ->Fill(ptOfTheDimuonPair);
                  /* -
                   * - Variable pt-binning.
                   * -
                   */
                  // if (        ptOfTheDimuonPair < 0.500 ) {
                  if (        ptOfTheDimuonPair < 0.275 ) {
                    fDimuonPtDistributionZNCzeroZNAzeroHv3HigherSide           ->Fill(ptOfTheDimuonPair);
                  } else if ( ptOfTheDimuonPair < 0.950 ) {
                    fDimuonPtDistributionZNCzeroZNAzeroHv3HigherSide           ->Fill( ptOfTheDimuonPair, 0.33333333333 );
                  } else if ( ptOfTheDimuonPair < 1.400 ) {
                    fDimuonPtDistributionZNCzeroZNAzeroHv3HigherSide           ->Fill( ptOfTheDimuonPair, 0.16666666666 );
                  } else if ( ptOfTheDimuonPair < 2.000 ) {
                    fDimuonPtDistributionZNCzeroZNAzeroHv3HigherSide           ->Fill( ptOfTheDimuonPair, 0.125 );
                  } else if ( ptOfTheDimuonPair < 4.000 ) {
                    fDimuonPtDistributionZNCzeroZNAzeroHv3HigherSide           ->Fill( ptOfTheDimuonPair, 0.050 );
                  } else if ( ptOfTheDimuonPair < 5.000 ) {
                    fDimuonPtDistributionZNCzeroZNAzeroHv3HigherSide           ->Fill( ptOfTheDimuonPair, 0.025 );
                  }
                  // fDimuonPtDistributionZNCzeroZNAzeroShiftPlusOneHv2->Fill(ptOfTheDimuonPair);
                  if (        possibleJPsi.Rapidity() > -4.0  && possibleJPsi.Rapidity() <= -3.50 ) {
                    fDimuonPtDistributionZNCzeroZNAzeroRapidityHv2HigherSide[0]->Fill(ptOfTheDimuonPair);
                    /* -
                     * - Variable pt-binning.
                     * -
                     */
                    // if (        ptOfTheDimuonPair < 0.500 ) {
                    if (        ptOfTheDimuonPair < 0.275 ) {
                      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3HigherSide[0]->Fill(ptOfTheDimuonPair);
                    } else if ( ptOfTheDimuonPair < 0.950 ) {
                      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3HigherSide[0]->Fill(ptOfTheDimuonPair, 0.33333333333);
                    } else if ( ptOfTheDimuonPair < 1.400 ) {
                      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3HigherSide[0]->Fill(ptOfTheDimuonPair, 0.16666666666);
                    } else if ( ptOfTheDimuonPair < 2.000 ) {
                      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3HigherSide[0]->Fill(ptOfTheDimuonPair, 0.125);
                    } else if ( ptOfTheDimuonPair < 4.000 ) {
                      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3HigherSide[0]->Fill(ptOfTheDimuonPair, 0.050);
                    } else if ( ptOfTheDimuonPair < 5.000 ) {
                      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3HigherSide[0]->Fill(ptOfTheDimuonPair, 0.025);
                    }
                  } else if ( possibleJPsi.Rapidity() > -3.50 && possibleJPsi.Rapidity() <= -3.00 ) {
                    fDimuonPtDistributionZNCzeroZNAzeroRapidityHv2HigherSide[1]->Fill(ptOfTheDimuonPair);
                    /* -
                     * - Variable pt-binning.
                     * -
                     */
                    // if (        ptOfTheDimuonPair < 0.500 ) {
                    if (        ptOfTheDimuonPair < 0.275 ) {
                      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3HigherSide[1]->Fill(ptOfTheDimuonPair);
                    } else if ( ptOfTheDimuonPair < 0.950 ) {
                      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3HigherSide[1]->Fill(ptOfTheDimuonPair, 0.33333333333);
                    } else if ( ptOfTheDimuonPair < 1.400 ) {
                      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3HigherSide[1]->Fill(ptOfTheDimuonPair, 0.16666666666);
                    } else if ( ptOfTheDimuonPair < 2.000 ) {
                      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3HigherSide[1]->Fill(ptOfTheDimuonPair, 0.125);
                    } else if ( ptOfTheDimuonPair < 4.000 ) {
                      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3HigherSide[1]->Fill(ptOfTheDimuonPair, 0.050);
                    } else if ( ptOfTheDimuonPair < 5.000 ) {
                      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3HigherSide[1]->Fill(ptOfTheDimuonPair, 0.025);
                    }
                  } else if ( possibleJPsi.Rapidity() > -3.00 && possibleJPsi.Rapidity() <= -2.50 ) {
                    fDimuonPtDistributionZNCzeroZNAzeroRapidityHv2HigherSide[2]->Fill(ptOfTheDimuonPair);
                    /* -
                     * - Variable pt-binning.
                     * -
                     */
                    // if (        ptOfTheDimuonPair < 0.500 ) {
                    if (        ptOfTheDimuonPair < 0.275 ) {
                      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3HigherSide[2]->Fill(ptOfTheDimuonPair);
                    } else if ( ptOfTheDimuonPair < 0.950 ) {
                      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3HigherSide[2]->Fill(ptOfTheDimuonPair, 0.33333333333);
                    } else if ( ptOfTheDimuonPair < 1.400 ) {
                      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3HigherSide[2]->Fill(ptOfTheDimuonPair, 0.16666666666);
                    } else if ( ptOfTheDimuonPair < 2.000 ) {
                      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3HigherSide[2]->Fill(ptOfTheDimuonPair, 0.125);
                    } else if ( ptOfTheDimuonPair < 4.000 ) {
                      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3HigherSide[2]->Fill(ptOfTheDimuonPair, 0.050);
                    } else if ( ptOfTheDimuonPair < 5.000 ) {
                      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3HigherSide[2]->Fill(ptOfTheDimuonPair, 0.025);
                    }
                  }
              }
              //
              // END SIDEBANDS
              //_______________________________
              /* -
               * - NORMAL analysis
               */
              if ( (possibleJPsi.Mag() > 2.85) && (possibleJPsi.Mag() < 3.35) ) {
                  fDimuonPtDistributionZNCzeroZNAzeroHv2             ->Fill(ptOfTheDimuonPair);
                  /* -
                   * - Variable pt-binning.
                   * -
                   */
                  // if (        ptOfTheDimuonPair < 0.500 ) {
                  if (        ptOfTheDimuonPair < 0.275 ) {
                    fDimuonPtDistributionZNCzeroZNAzeroHv3           ->Fill(ptOfTheDimuonPair);
                  } else if ( ptOfTheDimuonPair < 0.950 ) {
                    fDimuonPtDistributionZNCzeroZNAzeroHv3           ->Fill( ptOfTheDimuonPair, 0.33333333333 );
                  } else if ( ptOfTheDimuonPair < 1.400 ) {
                    fDimuonPtDistributionZNCzeroZNAzeroHv3           ->Fill( ptOfTheDimuonPair, 0.16666666666 );
                  } else if ( ptOfTheDimuonPair < 2.000 ) {
                    fDimuonPtDistributionZNCzeroZNAzeroHv3           ->Fill( ptOfTheDimuonPair, 0.125 );
                  } else if ( ptOfTheDimuonPair < 4.000 ) {
                    fDimuonPtDistributionZNCzeroZNAzeroHv3           ->Fill( ptOfTheDimuonPair, 0.050 );
                  } else if ( ptOfTheDimuonPair < 5.000 ) {
                    fDimuonPtDistributionZNCzeroZNAzeroHv3           ->Fill( ptOfTheDimuonPair, 0.025 );
                  }
                  // fDimuonPtDistributionZNCzeroZNAzeroShiftPlusOneHv2->Fill(ptOfTheDimuonPair);
                  if (        possibleJPsi.Rapidity() > -4.0  && possibleJPsi.Rapidity() <= -3.50 ) {
                    fDimuonPtDistributionZNCzeroZNAzeroRapidityHv2[0]->Fill(ptOfTheDimuonPair);
                    /* -
                     * - Variable pt-binning.
                     * -
                     */
                    // if (        ptOfTheDimuonPair < 0.500 ) {
                    if (        ptOfTheDimuonPair < 0.275 ) {
                      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3[0]->Fill(ptOfTheDimuonPair);
                    } else if ( ptOfTheDimuonPair < 0.950 ) {
                      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3[0]->Fill(ptOfTheDimuonPair, 0.33333333333);
                    } else if ( ptOfTheDimuonPair < 1.400 ) {
                      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3[0]->Fill(ptOfTheDimuonPair, 0.16666666666);
                    } else if ( ptOfTheDimuonPair < 2.000 ) {
                      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3[0]->Fill(ptOfTheDimuonPair, 0.125);
                    } else if ( ptOfTheDimuonPair < 4.000 ) {
                      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3[0]->Fill(ptOfTheDimuonPair, 0.050);
                    } else if ( ptOfTheDimuonPair < 5.000 ) {
                      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3[0]->Fill(ptOfTheDimuonPair, 0.025);
                    }
                  } else if ( possibleJPsi.Rapidity() > -3.50 && possibleJPsi.Rapidity() <= -3.00 ) {
                    fDimuonPtDistributionZNCzeroZNAzeroRapidityHv2[1]->Fill(ptOfTheDimuonPair);
                    /* -
                     * - Variable pt-binning.
                     * -
                     */
                    // if (        ptOfTheDimuonPair < 0.500 ) {
                    if (        ptOfTheDimuonPair < 0.275 ) {
                      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3[1]->Fill(ptOfTheDimuonPair);
                    } else if ( ptOfTheDimuonPair < 0.950 ) {
                      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3[1]->Fill(ptOfTheDimuonPair, 0.33333333333);
                    } else if ( ptOfTheDimuonPair < 1.400 ) {
                      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3[1]->Fill(ptOfTheDimuonPair, 0.16666666666);
                    } else if ( ptOfTheDimuonPair < 2.000 ) {
                      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3[1]->Fill(ptOfTheDimuonPair, 0.125);
                    } else if ( ptOfTheDimuonPair < 4.000 ) {
                      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3[1]->Fill(ptOfTheDimuonPair, 0.050);
                    } else if ( ptOfTheDimuonPair < 5.000 ) {
                      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3[1]->Fill(ptOfTheDimuonPair, 0.025);
                    }
                  } else if ( possibleJPsi.Rapidity() > -3.00 && possibleJPsi.Rapidity() <= -2.50 ) {
                    fDimuonPtDistributionZNCzeroZNAzeroRapidityHv2[2]->Fill(ptOfTheDimuonPair);
                    /* -
                     * - Variable pt-binning.
                     * -
                     */
                    // if (        ptOfTheDimuonPair < 0.500 ) {
                    if (        ptOfTheDimuonPair < 0.275 ) {
                      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3[2]->Fill(ptOfTheDimuonPair);
                    } else if ( ptOfTheDimuonPair < 0.950 ) {
                      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3[2]->Fill(ptOfTheDimuonPair, 0.33333333333);
                    } else if ( ptOfTheDimuonPair < 1.400 ) {
                      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3[2]->Fill(ptOfTheDimuonPair, 0.16666666666);
                    } else if ( ptOfTheDimuonPair < 2.000 ) {
                      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3[2]->Fill(ptOfTheDimuonPair, 0.125);
                    } else if ( ptOfTheDimuonPair < 4.000 ) {
                      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3[2]->Fill(ptOfTheDimuonPair, 0.050);
                    } else if ( ptOfTheDimuonPair < 5.000 ) {
                      fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3[2]->Fill(ptOfTheDimuonPair, 0.025);
                    }
                  }
              }
              if( ptOfTheDimuonPair < 0.25 ) {
              // if( ptOfTheDimuonPair < 0.10 ) {
                  fInvariantMassDistributionCoherentZNCzeroZNAzeroHv2->Fill(possibleJPsi.Mag());
                  if (        possibleJPsi.Rapidity() > -4.0  && possibleJPsi.Rapidity() <= -3.50 ) {
                    fInvariantMassDistributionCoherentZNCzeroZNAzeroRapidityHv2[0]->Fill(possibleJPsi.Mag());
                    if ( (possibleJPsi.Mag() > 2.85) && (possibleJPsi.Mag() < 3.35) ){
                      for( Int_t iChannel = 0; iChannel < 16; iChannel++ ){
                        fADmultiplicity0N0NclassRapidityH[iChannel]->Fill(ADmultiplicities[iChannel]);
                      }
                    }
                  } else if ( possibleJPsi.Rapidity() > -3.50 && possibleJPsi.Rapidity() <= -3.00 ) {
                    fInvariantMassDistributionCoherentZNCzeroZNAzeroRapidityHv2[1]->Fill(possibleJPsi.Mag());
                    if ( (possibleJPsi.Mag() > 2.85) && (possibleJPsi.Mag() < 3.35) ){
                      for( Int_t iChannel = 0; iChannel < 16; iChannel++ ){
                        fADmultiplicity0N0NclassRapidityH[iChannel + 16]->Fill(ADmultiplicities[iChannel]);
                      }
                    }
                  } else if ( possibleJPsi.Rapidity() > -3.00 && possibleJPsi.Rapidity() <= -2.50 ) {
                    fInvariantMassDistributionCoherentZNCzeroZNAzeroRapidityHv2[2]->Fill(possibleJPsi.Mag());
                    if ( (possibleJPsi.Mag() > 2.85) && (possibleJPsi.Mag() < 3.35) ){
                      for( Int_t iChannel = 0; iChannel < 16; iChannel++ ){
                        fADmultiplicity0N0NclassRapidityH[iChannel + 32]->Fill(ADmultiplicities[iChannel]);
                      }
                    }
                  }
                  if ( (possibleJPsi.Mag() > 2.85) && (possibleJPsi.Mag() < 3.35) ){
                    if ( ADAPastFutureBoolean == 0 && ADCPastFutureBoolean == 0 && VZEROAPastFutureBoolean == 0 && VZEROCPastFutureBoolean == 0 ){
                    for( Int_t iChannel = 0; iChannel < 16; iChannel++ ){
                      fADmultiplicity0N0NclassH[iChannel]->Fill(ADmultiplicities[iChannel]);
                    }
                    fADmultiplicity0N0NclassTotalH ->Fill(ADmultiplicitiesTotal);
                    fADAmultiplicity0N0NclassTotalH->Fill(ADAmultiplicitiesTotal);
                    fADCmultiplicity0N0NclassTotalH->Fill(ADCmultiplicitiesTotal);

                    fVZEROmultiplicity0N0NclassTotalH ->Fill(VZEROmultiplicitiesTotal);
                    fVZEROAmultiplicity0N0NclassTotalH->Fill(VZEROAmultiplicitiesTotal);
                    fVZEROCmultiplicity0N0NclassTotalH->Fill(VZEROCmultiplicitiesTotal);
                    }
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
              //_______________________________
              /* -
               * - CREATING SIDEBANDS
               * - (LOWER side)
               * - for templates
               * -
               */
              if ( (possibleJPsi.Mag() > 2.4) && (possibleJPsi.Mag() < 2.8) ) {
                  fDimuonPtDistributionZNCzeroZNAanyHv2LowerSide            ->Fill(ptOfTheDimuonPair);
                  /* -
                   * - Variable pt-binning.
                   * -
                   */
                  // if (        ptOfTheDimuonPair < 0.500 ) {
                  if (        ptOfTheDimuonPair < 0.275 ) {
                    fDimuonPtDistributionZNCzeroZNAanyHv3LowerSide           ->Fill(ptOfTheDimuonPair);
                  } else if ( ptOfTheDimuonPair < 0.950 ) {
                    fDimuonPtDistributionZNCzeroZNAanyHv3LowerSide           ->Fill( ptOfTheDimuonPair, 0.33333333333 );
                  } else if ( ptOfTheDimuonPair < 1.400 ) {
                    fDimuonPtDistributionZNCzeroZNAanyHv3LowerSide           ->Fill( ptOfTheDimuonPair, 0.16666666666 );
                  } else if ( ptOfTheDimuonPair < 2.000 ) {
                    fDimuonPtDistributionZNCzeroZNAanyHv3LowerSide           ->Fill( ptOfTheDimuonPair, 0.125 );
                  } else if ( ptOfTheDimuonPair < 4.000 ) {
                    fDimuonPtDistributionZNCzeroZNAanyHv3LowerSide           ->Fill( ptOfTheDimuonPair, 0.050 );
                  } else if ( ptOfTheDimuonPair < 5.000 ) {
                    fDimuonPtDistributionZNCzeroZNAanyHv3LowerSide           ->Fill( ptOfTheDimuonPair, 0.025 );
                  }
                  // fDimuonPtDistributionZNCzeroZNAanyShiftPlusOneHv2->Fill(ptOfTheDimuonPair);
                  if (        possibleJPsi.Rapidity() > -4.0  && possibleJPsi.Rapidity() <= -3.50 ) {
                    fDimuonPtDistributionZNCzeroZNAanyRapidityHv2LowerSide[0]->Fill(ptOfTheDimuonPair);
                    /* -
                     * - Variable pt-binning.
                     * -
                     */
                    // if (        ptOfTheDimuonPair < 0.500 ) {
                    if (        ptOfTheDimuonPair < 0.275 ) {
                      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3LowerSide[0]->Fill(ptOfTheDimuonPair);
                    } else if ( ptOfTheDimuonPair < 0.950 ) {
                      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3LowerSide[0]->Fill(ptOfTheDimuonPair, 0.33333333333);
                    } else if ( ptOfTheDimuonPair < 1.400 ) {
                      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3LowerSide[0]->Fill(ptOfTheDimuonPair, 0.16666666666);
                    } else if ( ptOfTheDimuonPair < 2.000 ) {
                      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3LowerSide[0]->Fill(ptOfTheDimuonPair, 0.125);
                    } else if ( ptOfTheDimuonPair < 4.000 ) {
                      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3LowerSide[0]->Fill(ptOfTheDimuonPair, 0.050);
                    } else if ( ptOfTheDimuonPair < 5.000 ) {
                      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3LowerSide[0]->Fill(ptOfTheDimuonPair, 0.025);
                    }
                  } else if ( possibleJPsi.Rapidity() > -3.50 && possibleJPsi.Rapidity() <= -3.00 ) {
                    fDimuonPtDistributionZNCzeroZNAanyRapidityHv2LowerSide[1]->Fill(ptOfTheDimuonPair);
                    /* -
                     * - Variable pt-binning.
                     * -
                     */
                    // if (        ptOfTheDimuonPair < 0.500 ) {
                    if (        ptOfTheDimuonPair < 0.275 ) {
                      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3LowerSide[1]->Fill(ptOfTheDimuonPair);
                    } else if ( ptOfTheDimuonPair < 0.950 ) {
                      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3LowerSide[1]->Fill(ptOfTheDimuonPair, 0.33333333333);
                    } else if ( ptOfTheDimuonPair < 1.400 ) {
                      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3LowerSide[1]->Fill(ptOfTheDimuonPair, 0.16666666666);
                    } else if ( ptOfTheDimuonPair < 2.000 ) {
                      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3LowerSide[1]->Fill(ptOfTheDimuonPair, 0.125);
                    } else if ( ptOfTheDimuonPair < 4.000 ) {
                      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3LowerSide[1]->Fill(ptOfTheDimuonPair, 0.050);
                    } else if ( ptOfTheDimuonPair < 5.000 ) {
                      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3LowerSide[1]->Fill(ptOfTheDimuonPair, 0.025);
                    }
                  } else if ( possibleJPsi.Rapidity() > -3.00 && possibleJPsi.Rapidity() <= -2.50 ) {
                    fDimuonPtDistributionZNCzeroZNAanyRapidityHv2LowerSide[2]->Fill(ptOfTheDimuonPair);
                    /* -
                     * - Variable pt-binning.
                     * -
                     */
                    // if (        ptOfTheDimuonPair < 0.500 ) {
                    if (        ptOfTheDimuonPair < 0.275 ) {
                      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3LowerSide[2]->Fill(ptOfTheDimuonPair);
                    } else if ( ptOfTheDimuonPair < 0.950 ) {
                      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3LowerSide[2]->Fill(ptOfTheDimuonPair, 0.33333333333);
                    } else if ( ptOfTheDimuonPair < 1.400 ) {
                      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3LowerSide[2]->Fill(ptOfTheDimuonPair, 0.16666666666);
                    } else if ( ptOfTheDimuonPair < 2.000 ) {
                      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3LowerSide[2]->Fill(ptOfTheDimuonPair, 0.125);
                    } else if ( ptOfTheDimuonPair < 4.000 ) {
                      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3LowerSide[2]->Fill(ptOfTheDimuonPair, 0.050);
                    } else if ( ptOfTheDimuonPair < 5.000 ) {
                      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3LowerSide[2]->Fill(ptOfTheDimuonPair, 0.025);
                    }
                  }
              }
              // END SIDEBANDS
              //_______________________________
              /* -
               * - CREATING SIDEBANDS
               * - (HIGHER side)
               * - for templates
               * -
               */
              if ( (possibleJPsi.Mag() > 4.0) && (possibleJPsi.Mag() < 5.5) ) {
                  fDimuonPtDistributionZNCzeroZNAanyHv2HigherSide            ->Fill(ptOfTheDimuonPair);
                  /* -
                   * - Variable pt-binning.
                   * -
                   */
                  // if (        ptOfTheDimuonPair < 0.500 ) {
                  if (        ptOfTheDimuonPair < 0.275 ) {
                    fDimuonPtDistributionZNCzeroZNAanyHv3HigherSide           ->Fill(ptOfTheDimuonPair);
                  } else if ( ptOfTheDimuonPair < 0.950 ) {
                    fDimuonPtDistributionZNCzeroZNAanyHv3HigherSide           ->Fill( ptOfTheDimuonPair, 0.33333333333 );
                  } else if ( ptOfTheDimuonPair < 1.400 ) {
                    fDimuonPtDistributionZNCzeroZNAanyHv3HigherSide           ->Fill( ptOfTheDimuonPair, 0.16666666666 );
                  } else if ( ptOfTheDimuonPair < 2.000 ) {
                    fDimuonPtDistributionZNCzeroZNAanyHv3HigherSide           ->Fill( ptOfTheDimuonPair, 0.125 );
                  } else if ( ptOfTheDimuonPair < 4.000 ) {
                    fDimuonPtDistributionZNCzeroZNAanyHv3HigherSide           ->Fill( ptOfTheDimuonPair, 0.050 );
                  } else if ( ptOfTheDimuonPair < 5.000 ) {
                    fDimuonPtDistributionZNCzeroZNAanyHv3HigherSide           ->Fill( ptOfTheDimuonPair, 0.025 );
                  }
                  // fDimuonPtDistributionZNCzeroZNAanyShiftPlusOneHv2->Fill(ptOfTheDimuonPair);
                  if (        possibleJPsi.Rapidity() > -4.0  && possibleJPsi.Rapidity() <= -3.50 ) {
                    fDimuonPtDistributionZNCzeroZNAanyRapidityHv2HigherSide[0]->Fill(ptOfTheDimuonPair);
                    /* -
                     * - Variable pt-binning.
                     * -
                     */
                    // if (        ptOfTheDimuonPair < 0.500 ) {
                    if (        ptOfTheDimuonPair < 0.275 ) {
                      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3HigherSide[0]->Fill(ptOfTheDimuonPair);
                    } else if ( ptOfTheDimuonPair < 0.950 ) {
                      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3HigherSide[0]->Fill(ptOfTheDimuonPair, 0.33333333333);
                    } else if ( ptOfTheDimuonPair < 1.400 ) {
                      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3HigherSide[0]->Fill(ptOfTheDimuonPair, 0.16666666666);
                    } else if ( ptOfTheDimuonPair < 2.000 ) {
                      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3HigherSide[0]->Fill(ptOfTheDimuonPair, 0.125);
                    } else if ( ptOfTheDimuonPair < 4.000 ) {
                      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3HigherSide[0]->Fill(ptOfTheDimuonPair, 0.050);
                    } else if ( ptOfTheDimuonPair < 5.000 ) {
                      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3HigherSide[0]->Fill(ptOfTheDimuonPair, 0.025);
                    }
                  } else if ( possibleJPsi.Rapidity() > -3.50 && possibleJPsi.Rapidity() <= -3.00 ) {
                    fDimuonPtDistributionZNCzeroZNAanyRapidityHv2HigherSide[1]->Fill(ptOfTheDimuonPair);
                    /* -
                     * - Variable pt-binning.
                     * -
                     */
                    // if (        ptOfTheDimuonPair < 0.500 ) {
                    if (        ptOfTheDimuonPair < 0.275 ) {
                      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3HigherSide[1]->Fill(ptOfTheDimuonPair);
                    } else if ( ptOfTheDimuonPair < 0.950 ) {
                      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3HigherSide[1]->Fill(ptOfTheDimuonPair, 0.33333333333);
                    } else if ( ptOfTheDimuonPair < 1.400 ) {
                      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3HigherSide[1]->Fill(ptOfTheDimuonPair, 0.16666666666);
                    } else if ( ptOfTheDimuonPair < 2.000 ) {
                      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3HigherSide[1]->Fill(ptOfTheDimuonPair, 0.125);
                    } else if ( ptOfTheDimuonPair < 4.000 ) {
                      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3HigherSide[1]->Fill(ptOfTheDimuonPair, 0.050);
                    } else if ( ptOfTheDimuonPair < 5.000 ) {
                      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3HigherSide[1]->Fill(ptOfTheDimuonPair, 0.025);
                    }
                  } else if ( possibleJPsi.Rapidity() > -3.00 && possibleJPsi.Rapidity() <= -2.50 ) {
                    fDimuonPtDistributionZNCzeroZNAanyRapidityHv2HigherSide[2]->Fill(ptOfTheDimuonPair);
                    /* -
                     * - Variable pt-binning.
                     * -
                     */
                    // if (        ptOfTheDimuonPair < 0.500 ) {
                    if (        ptOfTheDimuonPair < 0.275 ) {
                      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3HigherSide[2]->Fill(ptOfTheDimuonPair);
                    } else if ( ptOfTheDimuonPair < 0.950 ) {
                      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3HigherSide[2]->Fill(ptOfTheDimuonPair, 0.33333333333);
                    } else if ( ptOfTheDimuonPair < 1.400 ) {
                      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3HigherSide[2]->Fill(ptOfTheDimuonPair, 0.16666666666);
                    } else if ( ptOfTheDimuonPair < 2.000 ) {
                      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3HigherSide[2]->Fill(ptOfTheDimuonPair, 0.125);
                    } else if ( ptOfTheDimuonPair < 4.000 ) {
                      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3HigherSide[2]->Fill(ptOfTheDimuonPair, 0.050);
                    } else if ( ptOfTheDimuonPair < 5.000 ) {
                      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3HigherSide[2]->Fill(ptOfTheDimuonPair, 0.025);
                    }
                  }
              }
              // END SIDEBANDS
              //_______________________________
              /* -
               * - NORMAL analysis
               */
              if ( (possibleJPsi.Mag() > 2.85) && (possibleJPsi.Mag() < 3.35) ) {
                  fDimuonPtDistributionZNCzeroZNAanyHv2            ->Fill(ptOfTheDimuonPair);
                  /* -
                   * - Variable pt-binning.
                   * -
                   */
                  // if (        ptOfTheDimuonPair < 0.500 ) {
                  if (        ptOfTheDimuonPair < 0.275 ) {
                    fDimuonPtDistributionZNCzeroZNAanyHv3           ->Fill(ptOfTheDimuonPair);
                  } else if ( ptOfTheDimuonPair < 0.950 ) {
                    fDimuonPtDistributionZNCzeroZNAanyHv3           ->Fill( ptOfTheDimuonPair, 0.33333333333 );
                  } else if ( ptOfTheDimuonPair < 1.400 ) {
                    fDimuonPtDistributionZNCzeroZNAanyHv3           ->Fill( ptOfTheDimuonPair, 0.16666666666 );
                  } else if ( ptOfTheDimuonPair < 2.000 ) {
                    fDimuonPtDistributionZNCzeroZNAanyHv3           ->Fill( ptOfTheDimuonPair, 0.125 );
                  } else if ( ptOfTheDimuonPair < 4.000 ) {
                    fDimuonPtDistributionZNCzeroZNAanyHv3           ->Fill( ptOfTheDimuonPair, 0.050 );
                  } else if ( ptOfTheDimuonPair < 5.000 ) {
                    fDimuonPtDistributionZNCzeroZNAanyHv3           ->Fill( ptOfTheDimuonPair, 0.025 );
                  }
                  // fDimuonPtDistributionZNCzeroZNAanyShiftPlusOneHv2->Fill(ptOfTheDimuonPair);
                  if (        possibleJPsi.Rapidity() > -4.0  && possibleJPsi.Rapidity() <= -3.50 ) {
                    fDimuonPtDistributionZNCzeroZNAanyRapidityHv2[0]->Fill(ptOfTheDimuonPair);
                    /* -
                     * - Variable pt-binning.
                     * -
                     */
                    // if (        ptOfTheDimuonPair < 0.500 ) {
                    if (        ptOfTheDimuonPair < 0.275 ) {
                      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3[0]->Fill(ptOfTheDimuonPair);
                    } else if ( ptOfTheDimuonPair < 0.950 ) {
                      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3[0]->Fill(ptOfTheDimuonPair, 0.33333333333);
                    } else if ( ptOfTheDimuonPair < 1.400 ) {
                      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3[0]->Fill(ptOfTheDimuonPair, 0.16666666666);
                    } else if ( ptOfTheDimuonPair < 2.000 ) {
                      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3[0]->Fill(ptOfTheDimuonPair, 0.125);
                    } else if ( ptOfTheDimuonPair < 4.000 ) {
                      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3[0]->Fill(ptOfTheDimuonPair, 0.050);
                    } else if ( ptOfTheDimuonPair < 5.000 ) {
                      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3[0]->Fill(ptOfTheDimuonPair, 0.025);
                    }
                  } else if ( possibleJPsi.Rapidity() > -3.50 && possibleJPsi.Rapidity() <= -3.00 ) {
                    fDimuonPtDistributionZNCzeroZNAanyRapidityHv2[1]->Fill(ptOfTheDimuonPair);
                    /* -
                     * - Variable pt-binning.
                     * -
                     */
                    // if (        ptOfTheDimuonPair < 0.500 ) {
                    if (        ptOfTheDimuonPair < 0.275 ) {
                      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3[1]->Fill(ptOfTheDimuonPair);
                    } else if ( ptOfTheDimuonPair < 0.950 ) {
                      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3[1]->Fill(ptOfTheDimuonPair, 0.33333333333);
                    } else if ( ptOfTheDimuonPair < 1.400 ) {
                      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3[1]->Fill(ptOfTheDimuonPair, 0.16666666666);
                    } else if ( ptOfTheDimuonPair < 2.000 ) {
                      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3[1]->Fill(ptOfTheDimuonPair, 0.125);
                    } else if ( ptOfTheDimuonPair < 4.000 ) {
                      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3[1]->Fill(ptOfTheDimuonPair, 0.050);
                    } else if ( ptOfTheDimuonPair < 5.000 ) {
                      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3[1]->Fill(ptOfTheDimuonPair, 0.025);
                    }
                  } else if ( possibleJPsi.Rapidity() > -3.00 && possibleJPsi.Rapidity() <= -2.50 ) {
                    fDimuonPtDistributionZNCzeroZNAanyRapidityHv2[2]->Fill(ptOfTheDimuonPair);
                    /* -
                     * - Variable pt-binning.
                     * -
                     */
                    // if (        ptOfTheDimuonPair < 0.500 ) {
                    if (        ptOfTheDimuonPair < 0.275 ) {
                      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3[2]->Fill(ptOfTheDimuonPair);
                    } else if ( ptOfTheDimuonPair < 0.950 ) {
                      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3[2]->Fill(ptOfTheDimuonPair, 0.33333333333);
                    } else if ( ptOfTheDimuonPair < 1.400 ) {
                      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3[2]->Fill(ptOfTheDimuonPair, 0.16666666666);
                    } else if ( ptOfTheDimuonPair < 2.000 ) {
                      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3[2]->Fill(ptOfTheDimuonPair, 0.125);
                    } else if ( ptOfTheDimuonPair < 4.000 ) {
                      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3[2]->Fill(ptOfTheDimuonPair, 0.050);
                    } else if ( ptOfTheDimuonPair < 5.000 ) {
                      fDimuonPtDistributionZNCzeroZNAanyRapidityHv3[2]->Fill(ptOfTheDimuonPair, 0.025);
                    }
                  }
              }
              if( ptOfTheDimuonPair < 0.25 ) {
              // if( ptOfTheDimuonPair < 0.10 ) {
                  fInvariantMassDistributionCoherentZNCzeroZNAanyHv2->Fill(possibleJPsi.Mag());
                  // fDimuonPtDistributionCoherentZNCzeroZNAanyHv2->Fill(ptOfTheDimuonPair);
                  if (        possibleJPsi.Rapidity() > -4.0  && possibleJPsi.Rapidity() <= -3.50 ) {
                    fInvariantMassDistributionCoherentZNCzeroZNAanyRapidityHv2[0]->Fill(possibleJPsi.Mag());
                    if ( (possibleJPsi.Mag() > 2.85) && (possibleJPsi.Mag() < 3.35) ){
                      for( Int_t iChannel = 0; iChannel < 16; iChannel++ ){
                        fADmultiplicity0NXNclassRapidityH[iChannel]->Fill(ADmultiplicities[iChannel]);
                      }
                    }
                  } else if ( possibleJPsi.Rapidity() > -3.50 && possibleJPsi.Rapidity() <= -3.00 ) {
                    fInvariantMassDistributionCoherentZNCzeroZNAanyRapidityHv2[1]->Fill(possibleJPsi.Mag());
                    if ( (possibleJPsi.Mag() > 2.85) && (possibleJPsi.Mag() < 3.35) ){
                      for( Int_t iChannel = 0; iChannel < 16; iChannel++ ){
                        fADmultiplicity0NXNclassRapidityH[iChannel + 16]->Fill(ADmultiplicities[iChannel]);
                      }
                    }
                  } else if ( possibleJPsi.Rapidity() > -3.00 && possibleJPsi.Rapidity() <= -2.50 ) {
                    fInvariantMassDistributionCoherentZNCzeroZNAanyRapidityHv2[2]->Fill(possibleJPsi.Mag());
                    if ( (possibleJPsi.Mag() > 2.85) && (possibleJPsi.Mag() < 3.35) ){
                      for( Int_t iChannel = 0; iChannel < 16; iChannel++ ){
                        fADmultiplicity0NXNclassRapidityH[iChannel + 32]->Fill(ADmultiplicities[iChannel]);
                      }
                    }
                  }
                  if ( (possibleJPsi.Mag() > 2.85) && (possibleJPsi.Mag() < 3.35) ){
                    if ( ADAPastFutureBoolean == 0 && ADCPastFutureBoolean == 0 && VZEROAPastFutureBoolean == 0 && VZEROCPastFutureBoolean == 0 ){
                    for( Int_t iChannel = 0; iChannel < 16; iChannel++ ){
                      fADmultiplicity0NXNclassH[iChannel]->Fill(ADmultiplicities[iChannel]);
                    }
                    fADmultiplicity0NXNclassTotalH ->Fill(ADmultiplicitiesTotal);
                    fADAmultiplicity0NXNclassTotalH->Fill(ADAmultiplicitiesTotal);
                    fADCmultiplicity0NXNclassTotalH->Fill(ADCmultiplicitiesTotal);

                    fVZEROmultiplicity0NXNclassTotalH ->Fill(VZEROmultiplicitiesTotal);
                    fVZEROAmultiplicity0NXNclassTotalH->Fill(VZEROAmultiplicitiesTotal);
                    fVZEROCmultiplicity0NXNclassTotalH->Fill(VZEROCmultiplicitiesTotal);
                    }
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
              //_______________________________
              /* -
               * - CREATING SIDEBANDS
               * - (LOWER side)
               * - for templates
               * -
               */
              if ( (possibleJPsi.Mag() > 2.4) && (possibleJPsi.Mag() < 2.8) ) {
                  fDimuonPtDistributionZNCanyZNAzeroHv2LowerSide            ->Fill(ptOfTheDimuonPair);
                  /* -
                   * - Variable pt-binning.
                   * -
                   */
                  // if (        ptOfTheDimuonPair < 0.500 ) {
                  if (        ptOfTheDimuonPair < 0.275 ) {
                    fDimuonPtDistributionZNCanyZNAzeroHv3LowerSide           ->Fill(ptOfTheDimuonPair);
                  } else if ( ptOfTheDimuonPair < 0.950 ) {
                    fDimuonPtDistributionZNCanyZNAzeroHv3LowerSide           ->Fill( ptOfTheDimuonPair, 0.33333333333 );
                  } else if ( ptOfTheDimuonPair < 1.400 ) {
                    fDimuonPtDistributionZNCanyZNAzeroHv3LowerSide           ->Fill( ptOfTheDimuonPair, 0.16666666666 );
                  } else if ( ptOfTheDimuonPair < 2.000 ) {
                    fDimuonPtDistributionZNCanyZNAzeroHv3LowerSide           ->Fill( ptOfTheDimuonPair, 0.125 );
                  } else if ( ptOfTheDimuonPair < 4.000 ) {
                    fDimuonPtDistributionZNCanyZNAzeroHv3LowerSide           ->Fill( ptOfTheDimuonPair, 0.050 );
                  } else if ( ptOfTheDimuonPair < 5.000 ) {
                    fDimuonPtDistributionZNCanyZNAzeroHv3LowerSide           ->Fill( ptOfTheDimuonPair, 0.025 );
                  }
                  // fDimuonPtDistributionZNCanyZNAzeroShiftPlusOneHv2->Fill(ptOfTheDimuonPair);
                  if (        possibleJPsi.Rapidity() > -4.0  && possibleJPsi.Rapidity() <= -3.50 ) {
                    fDimuonPtDistributionZNCanyZNAzeroRapidityHv2LowerSide[0]->Fill(ptOfTheDimuonPair);
                    /* -
                     * - Variable pt-binning.
                     * -
                     */
                    // if (        ptOfTheDimuonPair < 0.500 ) {
                    if (        ptOfTheDimuonPair < 0.275 ) {
                      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3LowerSide[0]->Fill(ptOfTheDimuonPair);
                    } else if ( ptOfTheDimuonPair < 0.950 ) {
                      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3LowerSide[0]->Fill(ptOfTheDimuonPair, 0.33333333333);
                    } else if ( ptOfTheDimuonPair < 1.400 ) {
                      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3LowerSide[0]->Fill(ptOfTheDimuonPair, 0.16666666666);
                    } else if ( ptOfTheDimuonPair < 2.000 ) {
                      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3LowerSide[0]->Fill(ptOfTheDimuonPair, 0.125);
                    } else if ( ptOfTheDimuonPair < 4.000 ) {
                      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3LowerSide[0]->Fill(ptOfTheDimuonPair, 0.050);
                    } else if ( ptOfTheDimuonPair < 5.000 ) {
                      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3LowerSide[0]->Fill(ptOfTheDimuonPair, 0.025);
                    }
                  } else if ( possibleJPsi.Rapidity() > -3.50 && possibleJPsi.Rapidity() <= -3.00 ) {
                    fDimuonPtDistributionZNCanyZNAzeroRapidityHv2LowerSide[1]->Fill(ptOfTheDimuonPair);
                    /* -
                     * - Variable pt-binning.
                     * -
                     */
                    // if (        ptOfTheDimuonPair < 0.500 ) {
                    if (        ptOfTheDimuonPair < 0.275 ) {
                      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3LowerSide[1]->Fill(ptOfTheDimuonPair);
                    } else if ( ptOfTheDimuonPair < 0.950 ) {
                      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3LowerSide[1]->Fill(ptOfTheDimuonPair, 0.33333333333);
                    } else if ( ptOfTheDimuonPair < 1.400 ) {
                      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3LowerSide[1]->Fill(ptOfTheDimuonPair, 0.16666666666);
                    } else if ( ptOfTheDimuonPair < 2.000 ) {
                      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3LowerSide[1]->Fill(ptOfTheDimuonPair, 0.125);
                    } else if ( ptOfTheDimuonPair < 4.000 ) {
                      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3LowerSide[1]->Fill(ptOfTheDimuonPair, 0.050);
                    } else if ( ptOfTheDimuonPair < 5.000 ) {
                      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3LowerSide[1]->Fill(ptOfTheDimuonPair, 0.025);
                    }
                  } else if ( possibleJPsi.Rapidity() > -3.00 && possibleJPsi.Rapidity() <= -2.50 ) {
                    fDimuonPtDistributionZNCanyZNAzeroRapidityHv2LowerSide[2]->Fill(ptOfTheDimuonPair);
                    /* -
                     * - Variable pt-binning.
                     * -
                     */
                    // if (        ptOfTheDimuonPair < 0.500 ) {
                    if (        ptOfTheDimuonPair < 0.275 ) {
                      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3LowerSide[2]->Fill(ptOfTheDimuonPair);
                    } else if ( ptOfTheDimuonPair < 0.950 ) {
                      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3LowerSide[2]->Fill(ptOfTheDimuonPair, 0.33333333333);
                    } else if ( ptOfTheDimuonPair < 1.400 ) {
                      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3LowerSide[2]->Fill(ptOfTheDimuonPair, 0.16666666666);
                    } else if ( ptOfTheDimuonPair < 2.000 ) {
                      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3LowerSide[2]->Fill(ptOfTheDimuonPair, 0.125);
                    } else if ( ptOfTheDimuonPair < 4.000 ) {
                      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3LowerSide[2]->Fill(ptOfTheDimuonPair, 0.050);
                    } else if ( ptOfTheDimuonPair < 5.000 ) {
                      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3LowerSide[2]->Fill(ptOfTheDimuonPair, 0.025);
                    }
                  }
              }
              // END SIDEBANDS
              //_______________________________
              /* -
               * - CREATING SIDEBANDS
               * - (HIGHER side)
               * - for templates
               * -
               */
              if ( (possibleJPsi.Mag() > 4.0) && (possibleJPsi.Mag() < 5.5) ) {
                  fDimuonPtDistributionZNCanyZNAzeroHv2HigherSide            ->Fill(ptOfTheDimuonPair);
                  /* -
                   * - Variable pt-binning.
                   * -
                   */
                  // if (        ptOfTheDimuonPair < 0.500 ) {
                  if (        ptOfTheDimuonPair < 0.275 ) {
                    fDimuonPtDistributionZNCanyZNAzeroHv3HigherSide           ->Fill(ptOfTheDimuonPair);
                  } else if ( ptOfTheDimuonPair < 0.950 ) {
                    fDimuonPtDistributionZNCanyZNAzeroHv3HigherSide           ->Fill( ptOfTheDimuonPair, 0.33333333333 );
                  } else if ( ptOfTheDimuonPair < 1.400 ) {
                    fDimuonPtDistributionZNCanyZNAzeroHv3HigherSide           ->Fill( ptOfTheDimuonPair, 0.16666666666 );
                  } else if ( ptOfTheDimuonPair < 2.000 ) {
                    fDimuonPtDistributionZNCanyZNAzeroHv3HigherSide           ->Fill( ptOfTheDimuonPair, 0.125 );
                  } else if ( ptOfTheDimuonPair < 4.000 ) {
                    fDimuonPtDistributionZNCanyZNAzeroHv3HigherSide           ->Fill( ptOfTheDimuonPair, 0.050 );
                  } else if ( ptOfTheDimuonPair < 5.000 ) {
                    fDimuonPtDistributionZNCanyZNAzeroHv3HigherSide           ->Fill( ptOfTheDimuonPair, 0.025 );
                  }
                  // fDimuonPtDistributionZNCanyZNAzeroShiftPlusOneHv2->Fill(ptOfTheDimuonPair);
                  if (        possibleJPsi.Rapidity() > -4.0  && possibleJPsi.Rapidity() <= -3.50 ) {
                    fDimuonPtDistributionZNCanyZNAzeroRapidityHv2HigherSide[0]->Fill(ptOfTheDimuonPair);
                    /* -
                     * - Variable pt-binning.
                     * -
                     */
                    // if (        ptOfTheDimuonPair < 0.500 ) {
                    if (        ptOfTheDimuonPair < 0.275 ) {
                      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3HigherSide[0]->Fill(ptOfTheDimuonPair);
                    } else if ( ptOfTheDimuonPair < 0.950 ) {
                      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3HigherSide[0]->Fill(ptOfTheDimuonPair, 0.33333333333);
                    } else if ( ptOfTheDimuonPair < 1.400 ) {
                      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3HigherSide[0]->Fill(ptOfTheDimuonPair, 0.16666666666);
                    } else if ( ptOfTheDimuonPair < 2.000 ) {
                      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3HigherSide[0]->Fill(ptOfTheDimuonPair, 0.125);
                    } else if ( ptOfTheDimuonPair < 4.000 ) {
                      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3HigherSide[0]->Fill(ptOfTheDimuonPair, 0.050);
                    } else if ( ptOfTheDimuonPair < 5.000 ) {
                      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3HigherSide[0]->Fill(ptOfTheDimuonPair, 0.025);
                    }
                  } else if ( possibleJPsi.Rapidity() > -3.50 && possibleJPsi.Rapidity() <= -3.00 ) {
                    fDimuonPtDistributionZNCanyZNAzeroRapidityHv2HigherSide[1]->Fill(ptOfTheDimuonPair);
                    /* -
                     * - Variable pt-binning.
                     * -
                     */
                    // if (        ptOfTheDimuonPair < 0.500 ) {
                    if (        ptOfTheDimuonPair < 0.275 ) {
                      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3HigherSide[1]->Fill(ptOfTheDimuonPair);
                    } else if ( ptOfTheDimuonPair < 0.950 ) {
                      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3HigherSide[1]->Fill(ptOfTheDimuonPair, 0.33333333333);
                    } else if ( ptOfTheDimuonPair < 1.400 ) {
                      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3HigherSide[1]->Fill(ptOfTheDimuonPair, 0.16666666666);
                    } else if ( ptOfTheDimuonPair < 2.000 ) {
                      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3HigherSide[1]->Fill(ptOfTheDimuonPair, 0.125);
                    } else if ( ptOfTheDimuonPair < 4.000 ) {
                      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3HigherSide[1]->Fill(ptOfTheDimuonPair, 0.050);
                    } else if ( ptOfTheDimuonPair < 5.000 ) {
                      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3HigherSide[1]->Fill(ptOfTheDimuonPair, 0.025);
                    }
                  } else if ( possibleJPsi.Rapidity() > -3.00 && possibleJPsi.Rapidity() <= -2.50 ) {
                    fDimuonPtDistributionZNCanyZNAzeroRapidityHv2HigherSide[2]->Fill(ptOfTheDimuonPair);
                    /* -
                     * - Variable pt-binning.
                     * -
                     */
                    // if (        ptOfTheDimuonPair < 0.500 ) {
                    if (        ptOfTheDimuonPair < 0.275 ) {
                      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3HigherSide[2]->Fill(ptOfTheDimuonPair);
                    } else if ( ptOfTheDimuonPair < 0.950 ) {
                      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3HigherSide[2]->Fill(ptOfTheDimuonPair, 0.33333333333);
                    } else if ( ptOfTheDimuonPair < 1.400 ) {
                      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3HigherSide[2]->Fill(ptOfTheDimuonPair, 0.16666666666);
                    } else if ( ptOfTheDimuonPair < 2.000 ) {
                      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3HigherSide[2]->Fill(ptOfTheDimuonPair, 0.125);
                    } else if ( ptOfTheDimuonPair < 4.000 ) {
                      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3HigherSide[2]->Fill(ptOfTheDimuonPair, 0.050);
                    } else if ( ptOfTheDimuonPair < 5.000 ) {
                      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3HigherSide[2]->Fill(ptOfTheDimuonPair, 0.025);
                    }
                  }
              }
              // END SIDEBANDS
              //_______________________________
              /* -
               * - NORMAL analysis
               */
              if ( (possibleJPsi.Mag() > 2.85) && (possibleJPsi.Mag() < 3.35) ) {
                  fDimuonPtDistributionZNCanyZNAzeroHv2            ->Fill(ptOfTheDimuonPair);
                  /* -
                   * - Variable pt-binning.
                   * -
                   */
                  // if (        ptOfTheDimuonPair < 0.500 ) {
                  if (        ptOfTheDimuonPair < 0.275 ) {
                    fDimuonPtDistributionZNCanyZNAzeroHv3           ->Fill(ptOfTheDimuonPair);
                  } else if ( ptOfTheDimuonPair < 0.950 ) {
                    fDimuonPtDistributionZNCanyZNAzeroHv3           ->Fill( ptOfTheDimuonPair, 0.33333333333 );
                  } else if ( ptOfTheDimuonPair < 1.400 ) {
                    fDimuonPtDistributionZNCanyZNAzeroHv3           ->Fill( ptOfTheDimuonPair, 0.16666666666 );
                  } else if ( ptOfTheDimuonPair < 2.000 ) {
                    fDimuonPtDistributionZNCanyZNAzeroHv3           ->Fill( ptOfTheDimuonPair, 0.125 );
                  } else if ( ptOfTheDimuonPair < 4.000 ) {
                    fDimuonPtDistributionZNCanyZNAzeroHv3           ->Fill( ptOfTheDimuonPair, 0.050 );
                  } else if ( ptOfTheDimuonPair < 5.000 ) {
                    fDimuonPtDistributionZNCanyZNAzeroHv3           ->Fill( ptOfTheDimuonPair, 0.025 );
                  }
                  // fDimuonPtDistributionZNCanyZNAzeroShiftPlusOneHv2->Fill(ptOfTheDimuonPair);
                  if (        possibleJPsi.Rapidity() > -4.0  && possibleJPsi.Rapidity() <= -3.50 ) {
                    fDimuonPtDistributionZNCanyZNAzeroRapidityHv2[0]->Fill(ptOfTheDimuonPair);
                    /* -
                     * - Variable pt-binning.
                     * -
                     */
                    // if (        ptOfTheDimuonPair < 0.500 ) {
                    if (        ptOfTheDimuonPair < 0.275 ) {
                      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3[0]->Fill(ptOfTheDimuonPair);
                    } else if ( ptOfTheDimuonPair < 0.950 ) {
                      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3[0]->Fill(ptOfTheDimuonPair, 0.33333333333);
                    } else if ( ptOfTheDimuonPair < 1.400 ) {
                      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3[0]->Fill(ptOfTheDimuonPair, 0.16666666666);
                    } else if ( ptOfTheDimuonPair < 2.000 ) {
                      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3[0]->Fill(ptOfTheDimuonPair, 0.125);
                    } else if ( ptOfTheDimuonPair < 4.000 ) {
                      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3[0]->Fill(ptOfTheDimuonPair, 0.050);
                    } else if ( ptOfTheDimuonPair < 5.000 ) {
                      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3[0]->Fill(ptOfTheDimuonPair, 0.025);
                    }
                  } else if ( possibleJPsi.Rapidity() > -3.50 && possibleJPsi.Rapidity() <= -3.00 ) {
                    fDimuonPtDistributionZNCanyZNAzeroRapidityHv2[1]->Fill(ptOfTheDimuonPair);
                    /* -
                     * - Variable pt-binning.
                     * -
                     */
                    // if (        ptOfTheDimuonPair < 0.500 ) {
                    if (        ptOfTheDimuonPair < 0.275 ) {
                      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3[1]->Fill(ptOfTheDimuonPair);
                    } else if ( ptOfTheDimuonPair < 0.950 ) {
                      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3[1]->Fill(ptOfTheDimuonPair, 0.33333333333);
                    } else if ( ptOfTheDimuonPair < 1.400 ) {
                      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3[1]->Fill(ptOfTheDimuonPair, 0.16666666666);
                    } else if ( ptOfTheDimuonPair < 2.000 ) {
                      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3[1]->Fill(ptOfTheDimuonPair, 0.125);
                    } else if ( ptOfTheDimuonPair < 4.000 ) {
                      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3[1]->Fill(ptOfTheDimuonPair, 0.050);
                    } else if ( ptOfTheDimuonPair < 5.000 ) {
                      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3[1]->Fill(ptOfTheDimuonPair, 0.025);
                    }
                  } else if ( possibleJPsi.Rapidity() > -3.00 && possibleJPsi.Rapidity() <= -2.50 ) {
                    fDimuonPtDistributionZNCanyZNAzeroRapidityHv2[2]->Fill(ptOfTheDimuonPair);
                    /* -
                     * - Variable pt-binning.
                     * -
                     */
                    // if (        ptOfTheDimuonPair < 0.500 ) {
                    if (        ptOfTheDimuonPair < 0.275 ) {
                      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3[2]->Fill(ptOfTheDimuonPair);
                    } else if ( ptOfTheDimuonPair < 0.950 ) {
                      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3[2]->Fill(ptOfTheDimuonPair, 0.33333333333);
                    } else if ( ptOfTheDimuonPair < 1.400 ) {
                      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3[2]->Fill(ptOfTheDimuonPair, 0.16666666666);
                    } else if ( ptOfTheDimuonPair < 2.000 ) {
                      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3[2]->Fill(ptOfTheDimuonPair, 0.125);
                    } else if ( ptOfTheDimuonPair < 4.000 ) {
                      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3[2]->Fill(ptOfTheDimuonPair, 0.050);
                    } else if ( ptOfTheDimuonPair < 5.000 ) {
                      fDimuonPtDistributionZNCanyZNAzeroRapidityHv3[2]->Fill(ptOfTheDimuonPair, 0.025);
                    }
                  }
              }
              if( ptOfTheDimuonPair < 0.25 ) {
              // if( ptOfTheDimuonPair < 0.10 ) {
                  fInvariantMassDistributionCoherentZNCanyZNAzeroHv2->Fill(possibleJPsi.Mag());
                  // fDimuonPtDistributionCoherentZNCanyZNAzeroHv2->Fill(ptOfTheDimuonPair);
                  if (        possibleJPsi.Rapidity() > -4.0  && possibleJPsi.Rapidity() <= -3.50 ) {
                    fInvariantMassDistributionCoherentZNCanyZNAzeroRapidityHv2[0]->Fill(possibleJPsi.Mag());
                    if ( (possibleJPsi.Mag() > 2.85) && (possibleJPsi.Mag() < 3.35) ){
                      for( Int_t iChannel = 0; iChannel < 16; iChannel++ ){
                        fADmultiplicityXN0NclassRapidityH[iChannel]->Fill(ADmultiplicities[iChannel]);
                      }
                    }
                  } else if ( possibleJPsi.Rapidity() > -3.50 && possibleJPsi.Rapidity() <= -3.00 ) {
                    fInvariantMassDistributionCoherentZNCanyZNAzeroRapidityHv2[1]->Fill(possibleJPsi.Mag());
                    if ( (possibleJPsi.Mag() > 2.85) && (possibleJPsi.Mag() < 3.35) ){
                      for( Int_t iChannel = 0; iChannel < 16; iChannel++ ){
                        fADmultiplicityXN0NclassRapidityH[iChannel + 16]->Fill(ADmultiplicities[iChannel]);
                      }
                    }
                  } else if ( possibleJPsi.Rapidity() > -3.00 && possibleJPsi.Rapidity() <= -2.50 ) {
                    fInvariantMassDistributionCoherentZNCanyZNAzeroRapidityHv2[2]->Fill(possibleJPsi.Mag());
                    if ( (possibleJPsi.Mag() > 2.85) && (possibleJPsi.Mag() < 3.35) ){
                      for( Int_t iChannel = 0; iChannel < 16; iChannel++ ){
                        fADmultiplicityXN0NclassRapidityH[iChannel + 32]->Fill(ADmultiplicities[iChannel]);
                      }
                    }
                  }
                  if ( (possibleJPsi.Mag() > 2.85) && (possibleJPsi.Mag() < 3.35) ){
                    if ( ADAPastFutureBoolean == 0 && ADCPastFutureBoolean == 0 && VZEROAPastFutureBoolean == 0 && VZEROCPastFutureBoolean == 0 ){
                    for( Int_t iChannel = 0; iChannel < 16; iChannel++ ){
                      fADmultiplicityXN0NclassH[iChannel]->Fill(ADmultiplicities[iChannel]);
                    }
                    fADmultiplicityXN0NclassTotalH ->Fill(ADmultiplicitiesTotal);
                    fADAmultiplicityXN0NclassTotalH->Fill(ADAmultiplicitiesTotal);
                    fADCmultiplicityXN0NclassTotalH->Fill(ADCmultiplicitiesTotal);

                    fVZEROmultiplicityXN0NclassTotalH ->Fill(VZEROmultiplicitiesTotal);
                    fVZEROAmultiplicityXN0NclassTotalH->Fill(VZEROAmultiplicitiesTotal);
                    fVZEROCmultiplicityXN0NclassTotalH->Fill(VZEROCmultiplicitiesTotal);
                    }
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
              /* -
               * - CREATING SIDEBANDS
               * - (LOWER side)
               * - for templates
               * -
               */
              if ( (possibleJPsi.Mag() > 2.4) && (possibleJPsi.Mag() < 2.8) ) {
                  fDimuonPtDistributionZNCanyZNAanyHv2LowerSide            ->Fill(ptOfTheDimuonPair);
                  /* -
                   * - Variable pt-binning.
                   * -
                   */
                  // if (        ptOfTheDimuonPair < 0.500 ) {
                  if (        ptOfTheDimuonPair < 0.275 ) {
                    fDimuonPtDistributionZNCanyZNAanyHv3LowerSide           ->Fill(ptOfTheDimuonPair);
                  } else if ( ptOfTheDimuonPair < 0.950 ) {
                    fDimuonPtDistributionZNCanyZNAanyHv3LowerSide           ->Fill( ptOfTheDimuonPair, 0.33333333333 );
                  } else if ( ptOfTheDimuonPair < 1.400 ) {
                    fDimuonPtDistributionZNCanyZNAanyHv3LowerSide           ->Fill( ptOfTheDimuonPair, 0.16666666666 );
                  } else if ( ptOfTheDimuonPair < 2.000 ) {
                    fDimuonPtDistributionZNCanyZNAanyHv3LowerSide           ->Fill( ptOfTheDimuonPair, 0.125 );
                  } else if ( ptOfTheDimuonPair < 4.000 ) {
                    fDimuonPtDistributionZNCanyZNAanyHv3LowerSide           ->Fill( ptOfTheDimuonPair, 0.050 );
                  } else if ( ptOfTheDimuonPair < 5.000 ) {
                    fDimuonPtDistributionZNCanyZNAanyHv3LowerSide           ->Fill( ptOfTheDimuonPair, 0.025 );
                  }
                  // fDimuonPtDistributionZNCanyZNAanyShiftPlusOneHv2->Fill(ptOfTheDimuonPair);
                  if (        possibleJPsi.Rapidity() > -4.0  && possibleJPsi.Rapidity() <= -3.50 ) {
                    fDimuonPtDistributionZNCanyZNAanyRapidityHv2LowerSide[0]->Fill(ptOfTheDimuonPair);
                    /* -
                     * - Variable pt-binning.
                     * -
                     */
                    // if (        ptOfTheDimuonPair < 0.500 ) {
                    if (        ptOfTheDimuonPair < 0.275 ) {
                      fDimuonPtDistributionZNCanyZNAanyRapidityHv3LowerSide[0]->Fill(ptOfTheDimuonPair);
                    } else if ( ptOfTheDimuonPair < 0.950 ) {
                      fDimuonPtDistributionZNCanyZNAanyRapidityHv3LowerSide[0]->Fill( ptOfTheDimuonPair, 0.33333333333 );
                    } else if ( ptOfTheDimuonPair < 1.400 ) {
                      fDimuonPtDistributionZNCanyZNAanyRapidityHv3LowerSide[0]->Fill( ptOfTheDimuonPair, 0.16666666666 );
                    } else if ( ptOfTheDimuonPair < 2.000 ) {
                      fDimuonPtDistributionZNCanyZNAanyRapidityHv3LowerSide[0]->Fill( ptOfTheDimuonPair, 0.125 );
                    } else if ( ptOfTheDimuonPair < 4.000 ) {
                      fDimuonPtDistributionZNCanyZNAanyRapidityHv3LowerSide[0]->Fill( ptOfTheDimuonPair, 0.050 );
                    } else if ( ptOfTheDimuonPair < 5.000 ) {
                      fDimuonPtDistributionZNCanyZNAanyRapidityHv3LowerSide[0]->Fill( ptOfTheDimuonPair, 0.025 );
                    }
                  } else if ( possibleJPsi.Rapidity() > -3.50 && possibleJPsi.Rapidity() <= -3.00 ) {
                    fDimuonPtDistributionZNCanyZNAanyRapidityHv2LowerSide[1]->Fill(ptOfTheDimuonPair);
                    /* -
                     * - Variable pt-binning.
                     * -
                     */
                    // if (        ptOfTheDimuonPair < 0.500 ) {
                    if (        ptOfTheDimuonPair < 0.275 ) {
                      fDimuonPtDistributionZNCanyZNAanyRapidityHv3LowerSide[1]->Fill(ptOfTheDimuonPair);
                    } else if ( ptOfTheDimuonPair < 0.950 ) {
                      fDimuonPtDistributionZNCanyZNAanyRapidityHv3LowerSide[1]->Fill( ptOfTheDimuonPair, 0.33333333333 );
                    } else if ( ptOfTheDimuonPair < 1.400 ) {
                      fDimuonPtDistributionZNCanyZNAanyRapidityHv3LowerSide[1]->Fill( ptOfTheDimuonPair, 0.16666666666 );
                    } else if ( ptOfTheDimuonPair < 2.000 ) {
                      fDimuonPtDistributionZNCanyZNAanyRapidityHv3LowerSide[1]->Fill( ptOfTheDimuonPair, 0.125 );
                    } else if ( ptOfTheDimuonPair < 4.000 ) {
                      fDimuonPtDistributionZNCanyZNAanyRapidityHv3LowerSide[1]->Fill( ptOfTheDimuonPair, 0.050 );
                    } else if ( ptOfTheDimuonPair < 5.000 ) {
                      fDimuonPtDistributionZNCanyZNAanyRapidityHv3LowerSide[1]->Fill( ptOfTheDimuonPair, 0.025 );
                    }
                  } else if ( possibleJPsi.Rapidity() > -3.00 && possibleJPsi.Rapidity() <= -2.50 ) {
                    fDimuonPtDistributionZNCanyZNAanyRapidityHv2LowerSide[2]->Fill(ptOfTheDimuonPair);
                    /* -
                     * - Variable pt-binning.
                     * -
                     */
                    // if (        ptOfTheDimuonPair < 0.500 ) {
                    if (        ptOfTheDimuonPair < 0.275 ) {
                      fDimuonPtDistributionZNCanyZNAanyRapidityHv3LowerSide[2]->Fill(ptOfTheDimuonPair);
                    } else if ( ptOfTheDimuonPair < 0.950 ) {
                      fDimuonPtDistributionZNCanyZNAanyRapidityHv3LowerSide[2]->Fill( ptOfTheDimuonPair, 0.33333333333 );
                    } else if ( ptOfTheDimuonPair < 1.400 ) {
                      fDimuonPtDistributionZNCanyZNAanyRapidityHv3LowerSide[2]->Fill( ptOfTheDimuonPair, 0.16666666666 );
                    } else if ( ptOfTheDimuonPair < 2.000 ) {
                      fDimuonPtDistributionZNCanyZNAanyRapidityHv3LowerSide[2]->Fill( ptOfTheDimuonPair, 0.125 );
                    } else if ( ptOfTheDimuonPair < 4.000 ) {
                      fDimuonPtDistributionZNCanyZNAanyRapidityHv3LowerSide[2]->Fill( ptOfTheDimuonPair, 0.050 );
                    } else if ( ptOfTheDimuonPair < 5.000 ) {
                      fDimuonPtDistributionZNCanyZNAanyRapidityHv3LowerSide[2]->Fill( ptOfTheDimuonPair, 0.025 );
                    }
                  }
              }
              // END SIDEBANDS
              //_______________________________
              /* -
               * - CREATING SIDEBANDS
               * - (HIGHER side)
               * - for templates
               * -
               */
              if ( (possibleJPsi.Mag() > 4.0) && (possibleJPsi.Mag() < 5.5) ) {
                  fDimuonPtDistributionZNCanyZNAanyHv2HigherSide            ->Fill(ptOfTheDimuonPair);
                  /* -
                   * - Variable pt-binning.
                   * -
                   */
                  // if (        ptOfTheDimuonPair < 0.500 ) {
                  if (        ptOfTheDimuonPair < 0.275 ) {
                    fDimuonPtDistributionZNCanyZNAanyHv3HigherSide           ->Fill(ptOfTheDimuonPair);
                  } else if ( ptOfTheDimuonPair < 0.950 ) {
                    fDimuonPtDistributionZNCanyZNAanyHv3HigherSide           ->Fill( ptOfTheDimuonPair, 0.33333333333 );
                  } else if ( ptOfTheDimuonPair < 1.400 ) {
                    fDimuonPtDistributionZNCanyZNAanyHv3HigherSide           ->Fill( ptOfTheDimuonPair, 0.16666666666 );
                  } else if ( ptOfTheDimuonPair < 2.000 ) {
                    fDimuonPtDistributionZNCanyZNAanyHv3HigherSide           ->Fill( ptOfTheDimuonPair, 0.125 );
                  } else if ( ptOfTheDimuonPair < 4.000 ) {
                    fDimuonPtDistributionZNCanyZNAanyHv3HigherSide           ->Fill( ptOfTheDimuonPair, 0.050 );
                  } else if ( ptOfTheDimuonPair < 5.000 ) {
                    fDimuonPtDistributionZNCanyZNAanyHv3HigherSide           ->Fill( ptOfTheDimuonPair, 0.025 );
                  }
                  // fDimuonPtDistributionZNCanyZNAanyShiftPlusOneHv2->Fill(ptOfTheDimuonPair);
                  if (        possibleJPsi.Rapidity() > -4.0  && possibleJPsi.Rapidity() <= -3.50 ) {
                    fDimuonPtDistributionZNCanyZNAanyRapidityHv2HigherSide[0]->Fill(ptOfTheDimuonPair);
                    /* -
                     * - Variable pt-binning.
                     * -
                     */
                    // if (        ptOfTheDimuonPair < 0.500 ) {
                    if (        ptOfTheDimuonPair < 0.275 ) {
                      fDimuonPtDistributionZNCanyZNAanyRapidityHv3HigherSide[0]->Fill(ptOfTheDimuonPair);
                    } else if ( ptOfTheDimuonPair < 0.950 ) {
                      fDimuonPtDistributionZNCanyZNAanyRapidityHv3HigherSide[0]->Fill( ptOfTheDimuonPair, 0.33333333333 );
                    } else if ( ptOfTheDimuonPair < 1.400 ) {
                      fDimuonPtDistributionZNCanyZNAanyRapidityHv3HigherSide[0]->Fill( ptOfTheDimuonPair, 0.16666666666 );
                    } else if ( ptOfTheDimuonPair < 2.000 ) {
                      fDimuonPtDistributionZNCanyZNAanyRapidityHv3HigherSide[0]->Fill( ptOfTheDimuonPair, 0.125 );
                    } else if ( ptOfTheDimuonPair < 4.000 ) {
                      fDimuonPtDistributionZNCanyZNAanyRapidityHv3HigherSide[0]->Fill( ptOfTheDimuonPair, 0.050 );
                    } else if ( ptOfTheDimuonPair < 5.000 ) {
                      fDimuonPtDistributionZNCanyZNAanyRapidityHv3HigherSide[0]->Fill( ptOfTheDimuonPair, 0.025 );
                    }
                  } else if ( possibleJPsi.Rapidity() > -3.50 && possibleJPsi.Rapidity() <= -3.00 ) {
                    fDimuonPtDistributionZNCanyZNAanyRapidityHv2HigherSide[1]->Fill(ptOfTheDimuonPair);
                    /* -
                     * - Variable pt-binning.
                     * -
                     */
                    // if (        ptOfTheDimuonPair < 0.500 ) {
                    if (        ptOfTheDimuonPair < 0.275 ) {
                      fDimuonPtDistributionZNCanyZNAanyRapidityHv3HigherSide[1]->Fill(ptOfTheDimuonPair);
                    } else if ( ptOfTheDimuonPair < 0.950 ) {
                      fDimuonPtDistributionZNCanyZNAanyRapidityHv3HigherSide[1]->Fill( ptOfTheDimuonPair, 0.33333333333 );
                    } else if ( ptOfTheDimuonPair < 1.400 ) {
                      fDimuonPtDistributionZNCanyZNAanyRapidityHv3HigherSide[1]->Fill( ptOfTheDimuonPair, 0.16666666666 );
                    } else if ( ptOfTheDimuonPair < 2.000 ) {
                      fDimuonPtDistributionZNCanyZNAanyRapidityHv3HigherSide[1]->Fill( ptOfTheDimuonPair, 0.125 );
                    } else if ( ptOfTheDimuonPair < 4.000 ) {
                      fDimuonPtDistributionZNCanyZNAanyRapidityHv3HigherSide[1]->Fill( ptOfTheDimuonPair, 0.050 );
                    } else if ( ptOfTheDimuonPair < 5.000 ) {
                      fDimuonPtDistributionZNCanyZNAanyRapidityHv3HigherSide[1]->Fill( ptOfTheDimuonPair, 0.025 );
                    }
                  } else if ( possibleJPsi.Rapidity() > -3.00 && possibleJPsi.Rapidity() <= -2.50 ) {
                    fDimuonPtDistributionZNCanyZNAanyRapidityHv2HigherSide[2]->Fill(ptOfTheDimuonPair);
                    /* -
                     * - Variable pt-binning.
                     * -
                     */
                    // if (        ptOfTheDimuonPair < 0.500 ) {
                    if (        ptOfTheDimuonPair < 0.275 ) {
                      fDimuonPtDistributionZNCanyZNAanyRapidityHv3HigherSide[2]->Fill(ptOfTheDimuonPair);
                    } else if ( ptOfTheDimuonPair < 0.950 ) {
                      fDimuonPtDistributionZNCanyZNAanyRapidityHv3HigherSide[2]->Fill( ptOfTheDimuonPair, 0.33333333333 );
                    } else if ( ptOfTheDimuonPair < 1.400 ) {
                      fDimuonPtDistributionZNCanyZNAanyRapidityHv3HigherSide[2]->Fill( ptOfTheDimuonPair, 0.16666666666 );
                    } else if ( ptOfTheDimuonPair < 2.000 ) {
                      fDimuonPtDistributionZNCanyZNAanyRapidityHv3HigherSide[2]->Fill( ptOfTheDimuonPair, 0.125 );
                    } else if ( ptOfTheDimuonPair < 4.000 ) {
                      fDimuonPtDistributionZNCanyZNAanyRapidityHv3HigherSide[2]->Fill( ptOfTheDimuonPair, 0.050 );
                    } else if ( ptOfTheDimuonPair < 5.000 ) {
                      fDimuonPtDistributionZNCanyZNAanyRapidityHv3HigherSide[2]->Fill( ptOfTheDimuonPair, 0.025 );
                    }
                  }
              }
              // END SIDEBANDS
              //_______________________________
              /* -
               * - NORMAL analysis
               */
              if ( (possibleJPsi.Mag() > 2.85) && (possibleJPsi.Mag() < 3.35) ) {
                  fDimuonPtDistributionZNCanyZNAanyHv2            ->Fill(ptOfTheDimuonPair);
                  /* -
                   * - Variable pt-binning.
                   * -
                   */
                  // if (        ptOfTheDimuonPair < 0.500 ) {
                  if (        ptOfTheDimuonPair < 0.275 ) {
                    fDimuonPtDistributionZNCanyZNAanyHv3           ->Fill(ptOfTheDimuonPair);
                  } else if ( ptOfTheDimuonPair < 0.950 ) {
                    fDimuonPtDistributionZNCanyZNAanyHv3           ->Fill( ptOfTheDimuonPair, 0.33333333333 );
                  } else if ( ptOfTheDimuonPair < 1.400 ) {
                    fDimuonPtDistributionZNCanyZNAanyHv3           ->Fill( ptOfTheDimuonPair, 0.16666666666 );
                  } else if ( ptOfTheDimuonPair < 2.000 ) {
                    fDimuonPtDistributionZNCanyZNAanyHv3           ->Fill( ptOfTheDimuonPair, 0.125 );
                  } else if ( ptOfTheDimuonPair < 4.000 ) {
                    fDimuonPtDistributionZNCanyZNAanyHv3           ->Fill( ptOfTheDimuonPair, 0.050 );
                  } else if ( ptOfTheDimuonPair < 5.000 ) {
                    fDimuonPtDistributionZNCanyZNAanyHv3           ->Fill( ptOfTheDimuonPair, 0.025 );
                  }
                  // fDimuonPtDistributionZNCanyZNAanyShiftPlusOneHv2->Fill(ptOfTheDimuonPair);
                  if (        possibleJPsi.Rapidity() > -4.0  && possibleJPsi.Rapidity() <= -3.50 ) {
                    fDimuonPtDistributionZNCanyZNAanyRapidityHv2[0]->Fill(ptOfTheDimuonPair);
                    /* -
                     * - Variable pt-binning.
                     * -
                     */
                    // if (        ptOfTheDimuonPair < 0.500 ) {
                    if (        ptOfTheDimuonPair < 0.275 ) {
                      fDimuonPtDistributionZNCanyZNAanyRapidityHv3[0]->Fill(ptOfTheDimuonPair);
                    } else if ( ptOfTheDimuonPair < 0.950 ) {
                      fDimuonPtDistributionZNCanyZNAanyRapidityHv3[0]->Fill( ptOfTheDimuonPair, 0.33333333333 );
                    } else if ( ptOfTheDimuonPair < 1.400 ) {
                      fDimuonPtDistributionZNCanyZNAanyRapidityHv3[0]->Fill( ptOfTheDimuonPair, 0.16666666666 );
                    } else if ( ptOfTheDimuonPair < 2.000 ) {
                      fDimuonPtDistributionZNCanyZNAanyRapidityHv3[0]->Fill( ptOfTheDimuonPair, 0.125 );
                    } else if ( ptOfTheDimuonPair < 4.000 ) {
                      fDimuonPtDistributionZNCanyZNAanyRapidityHv3[0]->Fill( ptOfTheDimuonPair, 0.050 );
                    } else if ( ptOfTheDimuonPair < 5.000 ) {
                      fDimuonPtDistributionZNCanyZNAanyRapidityHv3[0]->Fill( ptOfTheDimuonPair, 0.025 );
                    }
                  } else if ( possibleJPsi.Rapidity() > -3.50 && possibleJPsi.Rapidity() <= -3.00 ) {
                    fDimuonPtDistributionZNCanyZNAanyRapidityHv2[1]->Fill(ptOfTheDimuonPair);
                    /* -
                     * - Variable pt-binning.
                     * -
                     */
                    // if (        ptOfTheDimuonPair < 0.500 ) {
                    if (        ptOfTheDimuonPair < 0.275 ) {
                      fDimuonPtDistributionZNCanyZNAanyRapidityHv3[1]->Fill(ptOfTheDimuonPair);
                    } else if ( ptOfTheDimuonPair < 0.950 ) {
                      fDimuonPtDistributionZNCanyZNAanyRapidityHv3[1]->Fill( ptOfTheDimuonPair, 0.33333333333 );
                    } else if ( ptOfTheDimuonPair < 1.400 ) {
                      fDimuonPtDistributionZNCanyZNAanyRapidityHv3[1]->Fill( ptOfTheDimuonPair, 0.16666666666 );
                    } else if ( ptOfTheDimuonPair < 2.000 ) {
                      fDimuonPtDistributionZNCanyZNAanyRapidityHv3[1]->Fill( ptOfTheDimuonPair, 0.125 );
                    } else if ( ptOfTheDimuonPair < 4.000 ) {
                      fDimuonPtDistributionZNCanyZNAanyRapidityHv3[1]->Fill( ptOfTheDimuonPair, 0.050 );
                    } else if ( ptOfTheDimuonPair < 5.000 ) {
                      fDimuonPtDistributionZNCanyZNAanyRapidityHv3[1]->Fill( ptOfTheDimuonPair, 0.025 );
                    }
                  } else if ( possibleJPsi.Rapidity() > -3.00 && possibleJPsi.Rapidity() <= -2.50 ) {
                    fDimuonPtDistributionZNCanyZNAanyRapidityHv2[2]->Fill(ptOfTheDimuonPair);
                    /* -
                     * - Variable pt-binning.
                     * -
                     */
                    // if (        ptOfTheDimuonPair < 0.500 ) {
                    if (        ptOfTheDimuonPair < 0.275 ) {
                      fDimuonPtDistributionZNCanyZNAanyRapidityHv3[2]->Fill(ptOfTheDimuonPair);
                    } else if ( ptOfTheDimuonPair < 0.950 ) {
                      fDimuonPtDistributionZNCanyZNAanyRapidityHv3[2]->Fill( ptOfTheDimuonPair, 0.33333333333 );
                    } else if ( ptOfTheDimuonPair < 1.400 ) {
                      fDimuonPtDistributionZNCanyZNAanyRapidityHv3[2]->Fill( ptOfTheDimuonPair, 0.16666666666 );
                    } else if ( ptOfTheDimuonPair < 2.000 ) {
                      fDimuonPtDistributionZNCanyZNAanyRapidityHv3[2]->Fill( ptOfTheDimuonPair, 0.125 );
                    } else if ( ptOfTheDimuonPair < 4.000 ) {
                      fDimuonPtDistributionZNCanyZNAanyRapidityHv3[2]->Fill( ptOfTheDimuonPair, 0.050 );
                    } else if ( ptOfTheDimuonPair < 5.000 ) {
                      fDimuonPtDistributionZNCanyZNAanyRapidityHv3[2]->Fill( ptOfTheDimuonPair, 0.025 );
                    }
                  }
              }
              if( ptOfTheDimuonPair < 0.25 ) {
              // if( ptOfTheDimuonPair < 0.10 ) {
                  fInvariantMassDistributionCoherentZNCanyZNAanyHv2->Fill(possibleJPsi.Mag());
                  // fDimuonPtDistributionCoherentZNCanyZNAanyHv2->Fill(ptOfTheDimuonPair);
                  if (        possibleJPsi.Rapidity() > -4.0  && possibleJPsi.Rapidity() <= -3.50 ) {
                    fInvariantMassDistributionCoherentZNCanyZNAanyRapidityHv2[0]->Fill(possibleJPsi.Mag());
                    if ( (possibleJPsi.Mag() > 2.85) && (possibleJPsi.Mag() < 3.35) ){
                      for( Int_t iChannel = 0; iChannel < 16; iChannel++ ){
                        fADmultiplicityXNXNclassRapidityH[iChannel]->Fill(ADmultiplicities[iChannel]);
                      }
                    }
                  } else if ( possibleJPsi.Rapidity() > -3.50 && possibleJPsi.Rapidity() <= -3.00 ) {
                    fInvariantMassDistributionCoherentZNCanyZNAanyRapidityHv2[1]->Fill(possibleJPsi.Mag());
                    if ( (possibleJPsi.Mag() > 2.85) && (possibleJPsi.Mag() < 3.35) ){
                      for( Int_t iChannel = 0; iChannel < 16; iChannel++ ){
                        fADmultiplicityXNXNclassRapidityH[iChannel + 16]->Fill(ADmultiplicities[iChannel]);
                      }
                    }
                  } else if ( possibleJPsi.Rapidity() > -3.00 && possibleJPsi.Rapidity() <= -2.50 ) {
                    fInvariantMassDistributionCoherentZNCanyZNAanyRapidityHv2[2]->Fill(possibleJPsi.Mag());
                    if ( (possibleJPsi.Mag() > 2.85) && (possibleJPsi.Mag() < 3.35) ){
                      for( Int_t iChannel = 0; iChannel < 16; iChannel++ ){
                        fADmultiplicityXNXNclassRapidityH[iChannel + 32]->Fill(ADmultiplicities[iChannel]);
                      }
                    }
                  }
                  if ( (possibleJPsi.Mag() > 2.85) && (possibleJPsi.Mag() < 3.35) ){
                    if ( ADAPastFutureBoolean == 0 && ADCPastFutureBoolean == 0 && VZEROAPastFutureBoolean == 0 && VZEROCPastFutureBoolean == 0 ){
                    for( Int_t iChannel = 0; iChannel < 16; iChannel++ ){
                      fADmultiplicityXNXNclassH[iChannel]->Fill(ADmultiplicities[iChannel]);
                    }
                    fADmultiplicityXNXNclassTotalH ->Fill(ADmultiplicitiesTotal);
                    fADAmultiplicityXNXNclassTotalH->Fill(ADAmultiplicitiesTotal);
                    fADCmultiplicityXNXNclassTotalH->Fill(ADCmultiplicitiesTotal);

                    fVZEROmultiplicityXNXNclassTotalH ->Fill(VZEROmultiplicitiesTotal);
                    fVZEROAmultiplicityXNXNclassTotalH->Fill(VZEROAmultiplicitiesTotal);
                    fVZEROCmultiplicityXNXNclassTotalH->Fill(VZEROCmultiplicitiesTotal);
                    }
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













  /* -
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
  /* - If we are in the J/Psi peak, hence 2.8 < M < 3.3 GeV/c, AND if we are
     - in the coherent regime, so if the Pt < 0.25 GeV/c, we fill the plots.
     -
     - In the following note that the rapidity is well computed, so we are
     - dealing with negative values... -4.0 < Y < -2.5 !!!
     -
   */
  Double_t possibleJPsiCopyMag =possibleJPsiCopy.Mag();





  // post the data
  PostData(1, fOutputList);
}
//_____________________________________________________________________________
void AliAnalysisTaskADin2018::Terminate(Option_t *)
{
    cout << endl;
    // terminate
    // called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________
