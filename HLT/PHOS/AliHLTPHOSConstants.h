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
  Int_t GetMAXHOSTS() const { return fkMAXHOSTS;}				
  Int_t GetDEFAULTEVENTPORT() const { return fkDEFAULTEVENTPORT; }		
  Int_t GetNZROWSRCU() const { return fkNZROWSRCU;}				
  Int_t GetNXCOLUMNSRCU() const { return fkNXCOLUMNSRCU;} 			
  Int_t GetNZROWSMOD() const { return fkNZROWSMOD;} 				
  Int_t GetNXCOLUMNSMOD() const { return fkNXCOLUMNSMOD;} 			
  Int_t GetNDATATYPES() const { return fkNDATATYPES;} 				
  Int_t GetCSPSPERFEE() const { return fkCSPSPERFEE;} 				
  Int_t GetNMODULES() const { return fkNMODULES;} 				
  Int_t GetNRCUS() const { return fkNRCUS;} 					
  Int_t GetNRCUSPERMODULE() const { return fkNRCUSPERMODULE;} 			
  Int_t GetNRCUSPERTOTAL() const { return fkNRCUSPERTOTAL;} 			
  Int_t GetNFEECS() const { return fkNFEECS;} 					
  Int_t GetNALTROS() const { return fkNALTROS;} 				
  Int_t GetNALTROCHANNELS() const { return fkNALTROCHANNELS;} 			
  Int_t GetNBRANCHES() const { return fkNBRANCHES;} 				
  Float_t GetCELLSTEP() const { return fkCELLSTEP; } 				//   Int_t GetNRCUSPERSECTOR() const { return fkNRCUSPERSECTOR; } 			  Int_t GetDDLOFFSET() const { return fkDDLOFFSET; }
  TString GetDETNAME() const { return fkDETNAME; }
  
private:
  /** Constant members */
  const Int_t fkMAXHOSTS;                  //Constant
  const Int_t fkDEFAULTEVENTPORT;          //Constant
  const Int_t fkNZROWSRCU; /**<Number of rows per module*/ 
  const Int_t fkNXCOLUMNSRCU;          //Constant
  const Int_t fkNZROWSMOD;  /**<Number of rows per module*/ 
  const Int_t fkNXCOLUMNSMOD;  /**<Number of columns per module*/ 
  const Int_t fkNDATATYPES;          //Constant
  const Int_t fkCSPSPERFEE; //Constant
  
  const Int_t fkNMODULES;   /**<Number of modules of the PHOS detector*/
  const Int_t fkNRCUS;   /**<Number of RCUs per Module*/
 
  const Int_t fkNRCUSPERMODULE;   /**<Number of RCUs per Module*/
  const Int_t fkNRCUSPERTOTAL; /**<Total number of RCUs for PHOS*/
  const Int_t fkNFEECS;  /**<Number of Frontend cards per branch*/
  const Int_t fkNALTROS;  /**<Number of ALTROs per frontend card*/
  const Int_t fkNALTROCHANNELS;  //Constant
  const Int_t fkNBRANCHES; //Constant

  const Float_t fkCELLSTEP;  //Constant
  //  const Int_t fkNRCUSPERSECTOR;  //Constant
  const Int_t fkDDLOFFSET;  //Constant
  const TString fkDETNAME;  //Constant
   
  ClassDef(AliHLTPHOSConstants, 1);

};

#endif
