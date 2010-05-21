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
  fkMAXCELLSTEPETA = 6.32; 	// FR: max tower dimension along eta
  fkMINCELLSTEPETA = 5.99; 	// FR: min tower dimension along eta
  fkCELLSTEPPHI = 6.04667; 	// FR: tower dimension along phi
  fkCELLHEIGHT = 27.74;  	// FR: tower height
  fkCELLANGLE = 1.50;  		// FR: tower tapeiring angle (DEG)
  fkRADLENGTH = 1.23;
  fkCRITICENERGY = 8.0;
  fkCJ = 0.5 ;

}
				   
