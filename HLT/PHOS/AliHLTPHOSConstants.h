//-*- Mode: C++ -*-
// $Id$

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               */

/// @file   AliHLTPHOSConstants.h
/// @author Svein Lindal
/// @date   
/// @brief  Class containing constants for PHOS libraries.

#ifndef ALIHLTPHOSCONSTANTS_H
#define ALIHLTPHOSCONSTANTS_H


#include <TString.h>
#include "Rtypes.h"


#include "AliHLTCaloConstants.h"


class AliHLTPHOSConstants : public AliHLTCaloConstants
{
public:
  AliHLTPHOSConstants();
  virtual ~AliHLTPHOSConstants();
  virtual void InitConstants(); 
  virtual Int_t GetNZROWSMOD() const      { return PHOS::NZROWSMOD;} 
  virtual Int_t GetNXCOLUMNSMOD() const   { return PHOS::NXCOLUMNSMOD;}; 
  virtual Int_t GetNMODULES() const       { return PHOS::NMODULES; }; 
  virtual Int_t GetNRCUSPERMODULE() const { return PHOS::NRCUSPERMODULE;};
  virtual Int_t GetNFEECS() const         { return PHOS::NFEECS; } ;
  
  ClassDef(AliHLTPHOSConstants, 1)
};

#endif

