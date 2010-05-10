//-*- Mode: C++ -*-
// $Id: AliHLTEMCALConstants.cxx $
//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Svein Lindal <slindal@fys.uio.no>                     *
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

/// @file   AliHLTEMCALConstants.cxx
/// @author Svein Lindal
/// @date   2009-11-12
/// @brief  Class containing constants for EMCAL and EMCAL
///         loaded libraries


//#include "AliHLTCaloConstants.h"
#include "AliHLTEMCALConstants.h"

ClassImp(AliHLTEMCALConstants)

AliHLTEMCALConstants::AliHLTEMCALConstants() : AliHLTCaloConstants()
{
  InitConstants();
  //Default constructor
}



AliHLTEMCALConstants::~AliHLTEMCALConstants()
{
  //Default destructor
}


void 
AliHLTEMCALConstants::InitConstants() 
{
  fkDETNAME = "EMCAL";
  fkCELLSTEP = 6.0;
  fkDDLOFFSET =  4608;
}
				   
