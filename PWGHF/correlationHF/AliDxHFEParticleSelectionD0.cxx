// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE Project            * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                  Sedat Altinpinar <Sedat.Altinpinar@cern.ch>           *
//*                  Hege Erdal       <hege.erdal@gmail.com>               *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/// @file   AliDxHFEParticleSelectionD0.cxx
/// @author Sedat Altinpinar, Hege Erdal, Matthias Richter
/// @date   2012-03-19
/// @brief  D0 selection for D0-HFE correlation
///

#include "AliDxHFEParticleSelectionD0.h"
#include "AliVParticle.h"
#include "TObjArray.h"

/// ROOT macro for the implementation of ROOT specific class methods
ClassImp(AliDxHFEParticleSelectionD0)

AliDxHFEParticleSelectionD0::AliDxHFEParticleSelectionD0(const char* opt)
  : AliDxHFEParticleSelection(opt)
{
  // constructor
  // 
  // 
  // 
  // 
}

AliDxHFEParticleSelectionD0::~AliDxHFEParticleSelectionD0()
{
  // destructor
}

bool AliDxHFEParticleSelectionD0::IsSelected(AliVParticle* /*p*/)
{
  /// TODO: implement specific selection of D0 candidates
  return true;
}
