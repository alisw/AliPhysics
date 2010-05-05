//-*- Mode: C++ -*-
// $Id: AliHLTPHOSConstants.cxx $
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

/// @file   AliHLTPHOSConstants.cxx
/// @author Svein Lindal
/// @date   2009-11-12
/// @brief  Class containing constants for PHOS
///         loaded libraries

#include "AliHLTCaloConstants.h"
#include "AliHLTPHOSConstants.h"

ClassImp(AliHLTPHOSConstants);

AliHLTPHOSConstants::AliHLTPHOSConstants() :
  AliHLTCaloConstants()
{
  //Default constructor

  
  fkNZROWSRCU = 56;
  fkNXCOLUMNSRCU = 16;
  fkNZROWSMOD = 56;
  fkNXCOLUMNSMOD = 64;
  fkNMODULES= 5;
  fkNRCUS = 4; 
  fkNRCUSPERMODULE = 4;
  fkNRCUSPERTOTAL = fkNMODULES*fkNRCUSPERMODULE;
  fkNFEECS = 14;

  fkDETNAME = "PHOS";
  fkCELLSTEP = 2.255;
  fkDDLOFFSET = 1792;
}

AliHLTPHOSConstants::~AliHLTPHOSConstants()
{
  //Default destructor
}

