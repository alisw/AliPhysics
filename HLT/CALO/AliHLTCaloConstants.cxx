//-*- Mode: C++ -*-
// $Id: AliHLTCaloConstants.cxx $
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

/// @file   AliHLTCaloConstants.cxx
/// @author Svein Lindal
/// @date   2009-11-12
/// @brief  Class containing constants for EMCAL and PHOS
///         loaded libraries

#include "AliHLTCaloConstants.h"



//ClassImp(AliHLTCaloConstants);

AliHLTCaloConstants::AliHLTCaloConstants() :   fkCELLSTEP(-1),
					       fkMAXCELLSTEPETA(-1),
					       fkMINCELLSTEPETA(-1),
					       fkCELLSTEPPHI(-1),
					       fkCELLHEIGHT(-1),
					       fkCELLANGLE(-1),
					       fkRADLENGTH(-1),
					       fkCRITICENERGY(-1),
					       fkCJ(-1),
					       fkDDLOFFSET(-1),
					       fkNZROWSRCU(-1),
					       fkNXCOLUMNSRCU(-1),
					       fkNZROWSMOD(-1),
					       fkNXCOLUMNSMOD(-1),
					       fkNMODULES(-1),  
					       fkNRCUS(-1),  
					       fkNRCUSPERMODULE(-1),
					       fkNRCUSPERTOTAL(-1), 
					       fkNFEECS(-1), 
					       fkNALTROS(4),
					       fkNALTROCHANNELS(16),
					       fkNBRANCHES(2),
					       fkCSPSPERFEE(32),
					       fkNDATATYPES(10),
					       fkMAXHOSTS(20),
					       fkDEFAULTEVENTPORT(42001)
					       
					
{
  //Default constructor
}


AliHLTCaloConstants::~AliHLTCaloConstants()
{
  //Default destructor
}

