// $Id$
//***************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Gaute Ovrebekk <st05886@alf.uib.no>                   *
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

/// @file   AliHLTD0Candidate.cxx
/// @author Gaute Ovrebekk
/// @date   2010-11-19
/// @brief  Class for storing the D0 candidates

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTD0Candidate.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTD0Candidate)

AliHLTD0Candidate::AliHLTD0Candidate()
  : fInvD0(0)
  , fInvD0bar(0)
  , fpt(0)
  , flabelPos(0)
  , flabelNeg(0)
  , fPtPos(0)
  , fPtNeg(0)
{
}

AliHLTD0Candidate::AliHLTD0Candidate(Double_t InvD0,Double_t InvD0bar,Double_t pt,
				      Int_t lPos, Int_t lNeg, Double_t PtPos, Double_t PtNeg)
				     : fInvD0(InvD0)
				     , fInvD0bar(InvD0bar)
				     , fpt(pt)
				     , flabelPos(lPos)
				     , flabelNeg(lNeg)
				     , fPtPos(PtPos)
				     , fPtNeg(PtNeg)
{
}

AliHLTD0Candidate::~AliHLTD0Candidate()
{
}
