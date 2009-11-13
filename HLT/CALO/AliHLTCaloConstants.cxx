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
  fMAXHOSTS(20),
  fDEFAULTEVENTPORT(42001),
  fMAXBINVALUE(1023),
  fHIGHGAIN(0),
  fLOWGAIN(1),
  fALTROMAXSAMPLES(1008), 
  fALTROMAXPRESAMPLES(15),        
  fNZROWSRCU(56),         
  fNXCOLUMNSRCU(16), 
  fNZROWSMOD((det == "EMCAL") ? 48 : 56),            
  fNXCOLUMNSMOD((det == "EMCAL") ? 24 : 64 ),            
  fNGAINS(2),                           
  fNDATATYPES(10),    
  fPFMAXPATHLENGTH(256),
  fPFDEFAULTNSAMPLES(70),
  fPFDEFAULTSTARTINDEX(0),
  fDEFAULTTAU((det == "EMCAL") ? 0.2 : 2),
  fDEFAULTFS(10),  
  fMODULE0(0),
  fMODULE1(1),
  fMODULE2(2),
  fMODULE3(3),
  fMODULE4(4),
  fCSPSPERFEE(32),
  fRCU0(0),
  fRCU1(1),
  fRCU2(2),
  fRCU3(3),
  fZ0(0),
  fZ1(1),
  fX0(0),
  fX1(1),
  fNMODULES((det == "EMCAL") ? 13 : 5),   
  fNRCUS(4),       
  fNRCUSPERMODULE((det == "EMCAL") ? 2 : 4 ),                 
  fNRCUSPERTOTAL(fNMODULES*fNRCUSPERMODULE),
  fNFEECS((det == "EMCAL") ? 9 : 14),                            
  fNALTROS(4),                           
  fNALTROCHANNELS(16),
  fNBRANCHES(2),    
  fCELLSTEP((det == "EMCAL") ? -9999999.9 : 2.2 ),      ///XXX right value?
  fNRCUSPERSECTOR(4)
{
  //Default constructor
}

AliHLTCaloConstants::~AliHLTCaloConstants()
{
  //Default destructor
}

