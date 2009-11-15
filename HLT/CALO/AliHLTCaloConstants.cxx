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
#include "TString.h"

AliHLTCaloConstants::AliHLTCaloConstants(TString det) :
  fkMAXHOSTS(20),
  fkDEFAULTEVENTPORT(42001),
  fkMAXBINVALUE(1023),
  fkHIGHGAIN(0),
  fkLOWGAIN(1),
  fkALTROMAXSAMPLES(1008), 
  fkALTROMAXPRESAMPLES(15),        
  fkNZROWSRCU(56),         
  fkNXCOLUMNSRCU(16), 
  fkNZROWSMOD((det == "EMCAL") ? 48 : 56),            
  fkNXCOLUMNSMOD((det == "EMCAL") ? 24 : 64 ),            
  fkNGAINS(2),                           
  fkNDATATYPES(10),    
  fkPFMAXPATHLENGTH(256),
  fkPFDEFAULTNSAMPLES(70),
  fkPFDEFAULTSTARTINDEX(0),
  fkDEFAULTTAU((det == "EMCAL") ? 0.2 : 2),
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
  fkNMODULES((det == "EMCAL") ? 13 : 5),   
  fkNRCUS(4),       
  fkNRCUSPERMODULE((det == "EMCAL") ? 2 : 4 ),                 
  fkNRCUSPERTOTAL(fkNMODULES*fkNRCUSPERMODULE),
  fkNFEECS((det == "EMCAL") ? 9 : 14),                            
  fkNALTROS(4),                           
  fkNALTROCHANNELS(16),
  fkNBRANCHES(2),    
  fkCELLSTEP((det == "EMCAL") ? -9999999.9 : 2.2 ),      ///XXX right value?
  fkNRCUSPERSECTOR(4)
{
  //Default constructor
}

AliHLTCaloConstants::~AliHLTCaloConstants()
{
  //Default destructor
}

