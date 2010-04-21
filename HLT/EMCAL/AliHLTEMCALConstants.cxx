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

#include "AliHLTCaloConstants.h"
#include "AliHLTEMCALConstants.h"

ClassImp(AliHLTEMCALConstants);

AliHLTEMCALConstants::AliHLTEMCALConstants() :
  AliHLTCaloConstants(),
  fkMAXHOSTS(20),
  fkDEFAULTEVENTPORT(42001),
  fkMAXBINVALUE(1023),
  fkHIGHGAIN(1),
  fkLOWGAIN(0),
  fkALTROMAXSAMPLES(1008),
  fkALTROMAXPRESAMPLES(15),
  fkNZROWSRCU(56),
  fkNXCOLUMNSRCU(16),
  fkNZROWSMOD(48),
  fkNXCOLUMNSMOD(24),
  fkNGAINS(2),
  fkNDATATYPES(10),
  fkPFMAXPATHLENGTH(256),
  fkPFDEFAULTNSAMPLES(70),
  fkPFDEFAULTSTARTINDEX(0),
  fkDEFAULTTAU(0.2),
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
  fkNMODULES(10),
  fkNRCUS(4),
  fkNRCUSPERMODULE(2),
  fkNRCUSPERTOTAL(fkNMODULES*fkNRCUSPERMODULE),
  fkNFEECS(9),
  fkNALTROS(4),
  fkNALTROCHANNELS(16),
  fkNBRANCHES(2), 
  fkCELLSTEP(6.0),
  fkMAXCELLSTEPETA(6.32), 	// FR: max tower dimension along eta
  fkMINCELLSTEPETA(5.99), 	// FR: min tower dimension along eta
  fkCELLSTEPPHI(6.04667), 	// FR: tower dimension along phi
  fkCELLHEIGHT(27.74),  	// FR: tower height
  fkCELLANGLE(1.50),  		// FR: tower tapeiring angle (DEG)
  fkRADLENGTH(1.23),
  fkCRITICENERGY(8),
  fkCJ(0.5),
  fkNRCUSPERSECTOR(4),
  fkDDLOFFSET(4608),
  fkDETNAME("EMCAL")
{
  //Default constructor
}

AliHLTEMCALConstants::~AliHLTEMCALConstants()
{
  //Default destructor
}

