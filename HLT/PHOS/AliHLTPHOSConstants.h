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

class AliHLTCaloConstants;

class AliHLTPHOSConstants : public AliHLTCaloConstants
{

public:
  AliHLTPHOSConstants();
  ~AliHLTPHOSConstants();
  Int_t GetNZROWSRCU() const { return fkNZROWSRCU;}				
  Int_t GetNXCOLUMNSRCU() const { return fkNXCOLUMNSRCU;} 			
  Int_t GetNZROWSMOD() const { return fkNZROWSMOD;} 				
  Int_t GetNXCOLUMNSMOD() const { return fkNXCOLUMNSMOD;} 			
  Int_t GetNMODULES() const { return fkNMODULES;} 				
  Int_t GetNRCUS() const { return fkNRCUS;} 					
  Int_t GetNRCUSPERMODULE() const { return fkNRCUSPERMODULE;} 			
  Int_t GetNRCUSPERTOTAL() const { return fkNRCUSPERTOTAL;} 			
  Int_t GetNFEECS() const { return fkNFEECS;} 					
  //  Float_t GetCELLSTEP() const { return fkCELLSTEP; } 				
  
private:
  /** Constant members */
  const Int_t fkNZROWSRCU; /**<Number of rows per module*/ 
  const Int_t fkNXCOLUMNSRCU;          //Constant
  const Int_t fkNZROWSMOD;  /**<Number of rows per module*/ 
  const Int_t fkNXCOLUMNSMOD;  /**<Number of columns per module*/ 
  const Int_t fkNMODULES;   /**<Number of modules of the PHOS detector*/
  const Int_t fkNRCUS;   /**<Number of RCUs per Module*/
  const Int_t fkNRCUSPERMODULE;   /**<Number of RCUs per Module*/
  const Int_t fkNRCUSPERTOTAL; /**<Total number of RCUs for PHOS*/
  const Int_t fkNFEECS;  /**<Number of Frontend cards per branch*/
  
  ClassDef(AliHLTPHOSConstants, 1);

};

#endif
