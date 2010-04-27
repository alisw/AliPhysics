//-*- Mode: C++ -*-
// $Id: AliHLTCALOConstantsHandler.h 34622 2009-09-04 13:22:01Z odjuvsla $


//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               */

/// @file   AliHLTCaloConstantsHandler.h
/// @author Svein Lindal
/// @date   
/// @brief  Handler class that helps create an instance of the right 
///         AliHLTCaloConstants child class 
///         (e.g. AliHLTPHOSConstants or AliHLTEMCALConstants)


#ifndef ALIHLTCALOCONSTANTSHANDLER_H
#define ALIHLTCALOCONSTANTSHANDLER_H

#include "AliHLTCaloConstants.h"
#include "TString.h"

class AliHLTCaloConstantsHandler
{
public:
  AliHLTCaloConstantsHandler(TString det);
  virtual ~AliHLTCaloConstantsHandler();


protected:
  AliHLTCaloConstants* fCaloConstants;
  
private:
  

  /** Keep the standard constructor private, since we must alway initialize by specific calorimeter**/
  AliHLTCaloConstantsHandler();
  
  /** Keep the copy constructor private since it should not be used */
  AliHLTCaloConstantsHandler(const AliHLTCaloConstantsHandler & );
  
  /** Keep the assignement operator private since it should not be used */
  AliHLTCaloConstantsHandler & operator = (const AliHLTCaloConstantsHandler &);
  
  ClassDef(AliHLTCaloConstantsHandler, 1);

};

#endif
