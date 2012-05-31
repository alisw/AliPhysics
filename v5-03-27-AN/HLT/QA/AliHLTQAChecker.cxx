// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/** @file   AliHLTQAChecker.cxx
    @author Matthias Richter
    @date   2009-11-24
    @brief  HLT QA checker instance
*/
#include "AliHLTQAChecker.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTQAChecker)

AliHLTQAChecker::AliHLTQAChecker() 
  : AliQACheckerBase("HLT","TPC Quality Assurance Checker")
{
}

AliHLTQAChecker::~AliHLTQAChecker()
{
}

Double_t * AliHLTQAChecker::Check(AliQAv1::ALITASK_t /*task*/, TObjArray ** /*pTarget*/, const AliDetectorRecoParam * /*recoParam*/)
{
  return NULL;
}

void AliHLTQAChecker::Init(const AliQAv1::DETECTORINDEX_t /*det*/)
{
}

void AliHLTQAChecker::SetQA(AliQAv1::ALITASK_t /*index*/, Double_t * /*value*/) const
{
}

Double_t AliHLTQAChecker::CheckRAW(Int_t /*specie*/, TObjArray* /*list*/)
{
  return 0.0;
}

Double_t AliHLTQAChecker::CheckREC(Int_t /*specie*/, TObjArray* /*list*/)
{
  return 0.0;
}

Double_t AliHLTQAChecker::CheckESD(Int_t /*specie*/, TObjArray* /*list*/)
{
  return 0.0;
}
