// $Id: AliHLTCaloConstantsHandler.cxx 34223 2009-08-12 07:55:37Z richterm $

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Svein Lindal <slindal@fys.uio.no>
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

/// @file   AliHLTCaloConstantsHandler.h
/// @author Svein Lindal
/// @date   2009-10-21
/// @brief  Handler class that helps create an instance of the right 
///         AliHLTCaloConstants child class 
///         (e.g. AliHLTPHOSConstants or AliHLTEMCALConstants)


#include "TString.h"
#include "AliHLTCaloConstantsHandler.h"
#include "AliHLTCaloConstants.h"
#include "AliHLTMisc.h"
#include "AliHLTLogging.h"


ClassImp(AliHLTCaloConstantsHandler)


AliHLTCaloConstantsHandler::AliHLTCaloConstantsHandler(TString det):
  fCaloConstants(0)
{
  if (det.CompareTo("PHOS") == 0) {
    fCaloConstants = AliHLTMisc::LoadInstance( ( AliHLTCaloConstants*  ) NULL, 
					       "AliHLTPHOSConstants", "libAliHLTPHOS.so");
  } 
  else 
    {
      fCaloConstants = AliHLTMisc::LoadInstance( ( AliHLTCaloConstants* ) NULL, 
						 "AliHLTEMCALConstants" , "libAliHLTEMCAL.so");
    }
  
  if( fCaloConstants == 0 )
    {
      AliHLTLogging *log = new AliHLTLogging();
      log->Logging(kHLTLogFatal, __FILE__, "fCaloConstants == ZERO ",  "fCaloConstants == ZERO " ); 
      delete log;
    }
}


AliHLTCaloConstantsHandler::~AliHLTCaloConstantsHandler()
{
  //Default destructor
}
