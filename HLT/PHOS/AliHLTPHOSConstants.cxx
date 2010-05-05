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
  AliHLTCaloConstants(),
  fkMAXHOSTS(20),
  fkDEFAULTEVENTPORT(42001),
  fkNZROWSRCU(56),
  fkNXCOLUMNSRCU(16),
  fkNZROWSMOD(56),
  fkNXCOLUMNSMOD(64),
  fkNDATATYPES(10),
  fkPFMAXPATHLENGTH(256),
  fkPFDEFAULTNSAMPLES(70),
  fkPFDEFAULTSTARTINDEX(0),
  fkDEFAULTTAU(2.),
  fkDEFAULTFS(10),
  fkMODULE0(0),
  fkMODULE1(1),
  fkMODULE2(2),
  fkMODULE3(3),
  fkMODULE4(4),
  fkCSPSPERFEE(32),
  fkRCU0(0),
  fkRCU1(1),
  fkRCU2(2),
  fkRCU3(3),
  fkZ0(0),
  fkZ1(1),
  fkX0(0),
  fkX1(1),
  fkNMODULES(5),
  fkNRCUS(4),
  fkNRCUSPERMODULE(4),
  fkNRCUSPERTOTAL(fkNMODULES*fkNRCUSPERMODULE),
  fkNFEECS(14),
  fkNALTROS(4),
  fkNALTROCHANNELS(16),
  fkNBRANCHES(2),
  fkCELLSTEP(2.255),
  fkNRCUSPERSECTOR(-9999),
  fkDDLOFFSET(1792),
  fkDETNAME("PHOS")
{
  //Default constructor
}

AliHLTPHOSConstants::~AliHLTPHOSConstants()
{
  //Default destructor
}

