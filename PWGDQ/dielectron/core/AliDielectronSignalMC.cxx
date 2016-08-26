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

//_________________________________________________________________________
AliDielectronSignalMC::AliDielectronSignalMC() :
  TNamed("AliDielectronSignalMC", "AliDielectronSignalMC"),
  fCheckLikeSignPP(kFALSE),
  fCheckLikeSignMM(kFALSE),
  fLeg1(0),
  fLeg2(0),
  fMother1(0),
  fMother2(0),
  fGrandMother1(0),
  fGrandMother2(0),
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
  fMothersRelation(kUndefined),
  fGrandMothersRelation(kUndefined),
  fGEANTProcess(kPPrimary),
  fJpsiRadiative(kAll),
  fStackPDG(0),
  fCheckMotherGrandmother(kFALSE),
  fMotherIsGrandmother(kFALSE),
  fCheckStackForPDG(kFALSE),
  fFillPureMCStep(kFALSE) {

  //
  // Default constructor
  //
}


//_________________________________________________________________________
AliDielectronSignalMC::AliDielectronSignalMC(const Char_t* name, const Char_t* title) :
  TNamed(name, title),
  fCheckLikeSignPP(kFALSE),
  fCheckLikeSignMM(kFALSE),
  fLeg1(0),
  fLeg2(0),
  fMother1(0),
  fMother2(0),
  fGrandMother1(0),
  fGrandMother2(0),
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
  fMothersRelation(kUndefined),
  fGrandMothersRelation(kUndefined),
  fGEANTProcess(kPPrimary),
  fJpsiRadiative(kAll),
  fStackPDG(0),
  fCheckMotherGrandmother(kFALSE),
  fMotherIsGrandmother(kFALSE),
  fCheckStackForPDG(kFALSE),
  fFillPureMCStep(kFALSE){

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
