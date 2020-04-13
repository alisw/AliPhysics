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

///////////////////////////////////////////////////////////////////////////
//                Dielectron MC signal description container             //
//                                                                       //
//                                                                       //
/*
 * A container to describe the decay of a two body process
 *
 *
 *
 *
 *
 */
//                                                                       //
///////////////////////////////////////////////////////////////////////////


#include "AliDielectronSignalMC.h"

ClassImp(AliDielectronSignalMC)

const char *AliDielectronSignalMC::fgkJpsiSignals[AliDielectronSignalMC::kEnd] = {"none",
  "inclusiveJpsi", "beautyJpsi", "promptJpsi", "promptJpsiRad", "promptJpsiNonRad", "directJpsi", "gammaConversion", "gammaConversionDiffMother", "electrons", "directElec", "elecPrim", "fromBgTest"};

//_________________________________________________________________________
AliDielectronSignalMC::AliDielectronSignalMC() :
  TNamed("AliDielectronSignalMC", "AliDielectronSignalMC"),
  fCheckUnlikeSign(kTRUE),
  fCheckLikeSignPP(kFALSE),
  fCheckLikeSignMM(kFALSE),
  fLeg1(0),
  fLeg2(0),
  fMother1(0),
  fMother2(0),
  fGrandMother1(0),
  fGrandMother2(0),
  fStackPDG(0),
  fLeg1Exclude(kFALSE),
  fLeg2Exclude(kFALSE),
  fMother1Exclude(kFALSE),
  fMother2Exclude(kFALSE),
  fGrandMother1Exclude(kFALSE),
  fGrandMother2Exclude(kFALSE),
  fLeg1Source(kDontCare),
  fLeg2Source(kDontCare),
  fMother1Source(kDontCare),
  fMother2Source(kDontCare),
  fGrandMother1Source(kDontCare),
  fGrandMother2Source(kDontCare),
  fCheckBothChargesLeg1(kFALSE),
  fCheckBothChargesLeg2(kFALSE),
  fCheckBothChargesMother1(kFALSE),
  fCheckBothChargesMother2(kFALSE),
  fCheckBothChargesGrandMother1(kFALSE),
  fCheckBothChargesGrandMother2(kFALSE),
  fCheckGEANTProcess(kFALSE),
  fCheckMotherGrandmother(kFALSE),
  fCheckMotherGrandmotherDiffPair(kFALSE),
  fMotherIsGrandmother(kFALSE),
  fMotherIsGrandmotherDiffPair(kFALSE),
  fMothersRelation(kUndefined),
  fGrandMothersRelation(kUndefined),
  fGEANTProcess(kPPrimary),
  fJpsiRadiative(kAll),
  fCheckCorrelatedHF(kFALSE),
  fCheckStackForPDG(kFALSE),
  fFillPureMCStep(kFALSE)
{
  //
  // Default constructor
  //
}


//_________________________________________________________________________
AliDielectronSignalMC::AliDielectronSignalMC(const Char_t* name, const Char_t* title) :
  TNamed(name, title),
  fCheckUnlikeSign(kTRUE),
  fCheckLikeSignPP(kFALSE),
  fCheckLikeSignMM(kFALSE),
  fLeg1(0),
  fLeg2(0),
  fMother1(0),
  fMother2(0),
  fGrandMother1(0),
  fGrandMother2(0),
  fStackPDG(0),
  fLeg1Exclude(kFALSE),
  fLeg2Exclude(kFALSE),
  fMother1Exclude(kFALSE),
  fMother2Exclude(kFALSE),
  fGrandMother1Exclude(kFALSE),
  fGrandMother2Exclude(kFALSE),
  fLeg1Source(kDontCare),
  fLeg2Source(kDontCare),
  fMother1Source(kDontCare),
  fMother2Source(kDontCare),
  fGrandMother1Source(kDontCare),
  fGrandMother2Source(kDontCare),
  fCheckBothChargesLeg1(kFALSE),
  fCheckBothChargesLeg2(kFALSE),
  fCheckBothChargesMother1(kFALSE),
  fCheckBothChargesMother2(kFALSE),
  fCheckBothChargesGrandMother1(kFALSE),
  fCheckBothChargesGrandMother2(kFALSE),
  fCheckGEANTProcess(kFALSE),
  fCheckMotherGrandmother(kFALSE),
  fCheckMotherGrandmotherDiffPair(kFALSE),
  fMotherIsGrandmother(kFALSE),
  fMotherIsGrandmotherDiffPair(kFALSE),
  fMothersRelation(kUndefined),
  fGrandMothersRelation(kUndefined),
  fGEANTProcess(kPPrimary),
  fJpsiRadiative(kAll),
  fCheckCorrelatedHF(kFALSE),
  fCheckStackForPDG(kFALSE),
  fFillPureMCStep(kFALSE)
{
  //
  // Named constructor
  //
}



//_________________________________________________________________________
AliDielectronSignalMC::~AliDielectronSignalMC() {
  //
  //  Destructor
  //
}

//_________________________________________________________________________
AliDielectronSignalMC* AliDielectronSignalMC::GetJpsiMCsignalDef(EJpsiSignals kSignal)
{
  AliDielectronSignalMC *mcSignal = new AliDielectronSignalMC();
  mcSignal->SetName(fgkJpsiSignals[kSignal]);
  mcSignal->SetTitle(fgkJpsiSignals[kSignal]);
  switch (kSignal) {
    case kBegin:
      printf("No AliDielectronSignalMC defined for kBegin returning NULL");
      return 0x0;
    case kInclusiveJpsi:
      // Inclusive Jpsi
      // All jpsi available decaying into dielectron pairs
      mcSignal->SetLegPDGs(11,-11);
      mcSignal->SetMotherPDGs(443,443);
      mcSignal->SetMothersRelation(AliDielectronSignalMC::kSame);
      mcSignal->SetFillPureMCStep(kTRUE);
      mcSignal->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
      mcSignal->SetCheckBothChargesLegs(kTRUE,kTRUE);
      mcSignal->SetCheckBothChargesMothers(kTRUE,kTRUE);
      return mcSignal;
    case kBeautyJpsi:
      // Jpsi from beauty decays
      // Only b-Mesons, b-Baryons decay to fast to measure them separately from promptJpsi in real data (anyhow small branching ratio b-Baryons->Jpsi)
      mcSignal->SetLegPDGs(11,-11);
      mcSignal->SetMotherPDGs(443,443);
      mcSignal->SetMothersRelation(AliDielectronSignalMC::kSame);
      mcSignal->SetGrandMotherPDGs(500,500);
      mcSignal->SetFillPureMCStep(kTRUE);
      mcSignal->SetCheckBothChargesLegs(kTRUE,kTRUE);
      mcSignal->SetCheckBothChargesMothers(kTRUE,kTRUE);
      mcSignal->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
      return mcSignal;
    case kPromptJpsi:
      // Inclusive Jpsi without Jpsi from beauty decays
      mcSignal->SetLegPDGs(11,-11);
      mcSignal->SetMotherPDGs(443,443);
      mcSignal->SetGrandMotherPDGs(503,503,kTRUE,kTRUE);
      mcSignal->SetMothersRelation(AliDielectronSignalMC::kSame);
      mcSignal->SetFillPureMCStep(kTRUE);
      mcSignal->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
      mcSignal->SetCheckBothChargesLegs(kTRUE,kTRUE);
      mcSignal->SetCheckBothChargesMothers(kTRUE,kTRUE);
      mcSignal->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
      return mcSignal;
    case kPromptRadJpsi:
      mcSignal->SetLegPDGs(11,-11);
      mcSignal->SetMotherPDGs(443,443);
      mcSignal->SetGrandMotherPDGs(503,503,kTRUE,kTRUE);   // not from beauty hadrons
      mcSignal->SetMothersRelation(AliDielectronSignalMC::kSame);
      mcSignal->SetFillPureMCStep(kTRUE);
      mcSignal->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
      mcSignal->SetCheckBothChargesLegs(kTRUE,kTRUE);
      mcSignal->SetCheckBothChargesMothers(kTRUE,kTRUE);
      mcSignal->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
      mcSignal->SetJpsiRadiative(AliDielectronSignalMC::kIsRadiative);
      return mcSignal;
    case kPromptNonRadJpsi:
      mcSignal->SetLegPDGs(11,-11);
      mcSignal->SetMotherPDGs(443,443);
      mcSignal->SetGrandMotherPDGs(503,503,kTRUE,kTRUE);   // not from beauty hadrons
      mcSignal->SetMothersRelation(AliDielectronSignalMC::kSame);
      mcSignal->SetFillPureMCStep(kTRUE);
      mcSignal->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
      mcSignal->SetCheckBothChargesLegs(kTRUE,kTRUE);
      mcSignal->SetCheckBothChargesMothers(kTRUE,kTRUE);
      mcSignal->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
      mcSignal->SetJpsiRadiative(AliDielectronSignalMC::kIsNotRadiative);
      return mcSignal;
    case kDirectJpsi:
      mcSignal->SetLegPDGs(11,-11);
      mcSignal->SetMotherPDGs(443,443);
      mcSignal->SetMothersRelation(AliDielectronSignalMC::kSame);
      mcSignal->SetFillPureMCStep(kTRUE);
      mcSignal->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
      mcSignal->SetMotherSources(AliDielectronSignalMC::kDirect, AliDielectronSignalMC::kDirect);
      mcSignal->SetCheckBothChargesLegs(kTRUE,kTRUE);
      mcSignal->SetCheckBothChargesMothers(kTRUE,kTRUE);
      return mcSignal;
    case kGammaConv:
      mcSignal->SetLegPDGs(11,-11);
      mcSignal->SetCheckBothChargesLegs(kTRUE,kTRUE);
      mcSignal->SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
      mcSignal->SetMotherPDGs(22,22);
      mcSignal->SetMothersRelation(AliDielectronSignalMC::kSame);
      return mcSignal;
    case kGammaConvDiffMother:
      mcSignal->SetLegPDGs(11,-11);
      mcSignal->SetCheckBothChargesLegs(kTRUE,kTRUE);
      mcSignal->SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
      mcSignal->SetMotherPDGs(22,22);
      mcSignal->SetMothersRelation(AliDielectronSignalMC::kDifferent);
      return mcSignal;
    case kElectrons:
      // All electrons
      mcSignal->SetLegPDGs(11,1);  //dummy second leg (never MCtrue)
      mcSignal->SetCheckBothChargesLegs(kTRUE,kTRUE);
      // mcSignal->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
      mcSignal->SetGrandMotherPDGs(500,500,kTRUE,kTRUE); // exclude non-prompt jpsi eletrons
      mcSignal->SetFillPureMCStep(kTRUE);
      // mcSignal->SetMothersRelation(AliDielectronSignalMC::kSame);
      return mcSignal;
    case kDirectElectrons:
      mcSignal->SetLegPDGs(11,1); //NEW
      mcSignal->SetMothersRelation(AliDielectronSignalMC::kSame);
      mcSignal->SetGrandMotherPDGs(-1103,-1103);
      mcSignal->SetFillPureMCStep(kTRUE);
      mcSignal->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
      mcSignal->SetMotherSources(AliDielectronSignalMC::kDirect, AliDielectronSignalMC::kDirect);
      mcSignal->SetCheckBothChargesLegs(kTRUE,kTRUE);
      mcSignal->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
      mcSignal->SetGrandMotherPDGs(902,902,kTRUE,kTRUE); // exclude open charm,beauty hadrons
      return mcSignal;
    case kPrimaryElectrons:
      mcSignal->SetLegPDGs(11,1);  //dummy second leg (never MCtrue)
      mcSignal->SetCheckBothChargesLegs(kTRUE,kTRUE);
      mcSignal->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
      mcSignal->SetCheckBothChargesMothers(kTRUE,kTRUE);
      mcSignal->SetMotherPDGs(902,902,kTRUE,kTRUE); // exclude open charm,beauty hadrons
      mcSignal->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
      mcSignal->SetGrandMotherPDGs(902,902,kTRUE,kTRUE); // exclude open charm,beauty hadrons
      mcSignal->SetFillPureMCStep(kTRUE);
      return mcSignal;
    case kFromBgTest:
      mcSignal->SetLegPDGs(11,1);  //dummy second leg (never MCtrue)
      mcSignal->SetCheckBothChargesLegs(kTRUE,kTRUE);
      mcSignal->SetLegSources(AliDielectronSignalMC::kFromBGEvent, AliDielectronSignalMC::kFinalState);
      mcSignal->SetFillPureMCStep(kTRUE);
      return mcSignal;
    case kEnd:
      printf("No AliDielectronSignalMC defined for kEnd returning NULL");
      return 0x0;
  }
}
