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
  //  Int_t GetMAXBINVALUE() const { return fkMAXBINVALUE; } 			

  //  Int_t GetHIGHGAIN() const { return fkHIGHGAIN;} 				
  //  Int_t GetLOWGAIN() const { return fkLOWGAIN;} 				
  										
  //  Int_t GetALTROMAXSAMPLES() const {return fkALTROMAXSAMPLES;} 			
  //  Int_t GetALTROMAXPRESAMPLES() const { return fkALTROMAXPRESAMPLES;} 		
 										
  Int_t GetNZROWSRCU() const { return fkNZROWSRCU;}				
  Int_t GetNXCOLUMNSRCU() const { return fkNXCOLUMNSRCU;} 			
  Int_t GetNZROWSMOD() const { return fkNZROWSMOD;} 				
  Int_t GetNXCOLUMNSMOD() const { return fkNXCOLUMNSMOD;} 			
  //  Int_t GetNGAINS() const { return fkNGAINS;} 					
  Int_t GetNDATATYPES() const { return fkNDATATYPES;} 				
 										
  Int_t GetPFMAXPATHLENGTH() const { return fkPFMAXPATHLENGTH;} 		
 										
  Int_t GetPFDEFAULTNSAMPLES() const{ return fkPFDEFAULTNSAMPLES;} 		
  Int_t GetPFDEFAULTSTARTINDEX() const { return fkPFDEFAULTSTARTINDEX;} 	
 										
  Double_t GetDEFAULTTAU() const { return fkDEFAULTTAU;} 			
  Int_t GetDEFAULTFS() const { return fkDEFAULTFS;} 				
 										
  Int_t GetMODULE0() const { return fkMODULE0;} 				
  Int_t GetMODULE1() const { return fkMODULE1;} 				
  Int_t GetMODULE2() const { return fkMODULE2;} 				
  Int_t GetMODULE3() const { return fkMODULE3;} 				
  Int_t GetMODULE4() const { return fkMODULE4;} 				
 										
  Int_t GetCSPSPERFEE() const { return fkCSPSPERFEE;} 				
  Int_t GetRCU0() const { return fkRCU0;} 					
  Int_t GetRCU1() const { return fkRCU1;} 					
  Int_t GetRCU2() const { return fkRCU2;} 					
  Int_t GetRCU3() const { return fkRCU3;} 					
 										
  Int_t GetZ0() const { return fkZ0;} 						
  Int_t GetZ1() const { return fkZ1;} 						
  Int_t GetX0() const { return fkX0;} 						
  Int_t GetX1() const { return fkX1;} 						
 										
  Int_t GetNMODULES() const { return fkNMODULES;} 				
  Int_t GetNRCUS() const { return fkNRCUS;} 					
										
  Int_t GetNRCUSPERMODULE() const { return fkNRCUSPERMODULE;} 			
  Int_t GetNRCUSPERTOTAL() const { return fkNRCUSPERTOTAL;} 			
  Int_t GetNFEECS() const { return fkNFEECS;} 					
  Int_t GetNALTROS() const { return fkNALTROS;} 				
  Int_t GetNALTROCHANNELS() const { return fkNALTROCHANNELS;} 			
  Int_t GetNBRANCHES() const { return fkNBRANCHES;} 				
										
  Float_t GetCELLSTEP() const { return fkCELLSTEP; } 					
  Int_t GetNRCUSPERSECTOR() const { return fkNRCUSPERSECTOR; } 					
  
  Int_t GetDDLOFFSET() const { return fkDDLOFFSET; }
  
  TString GetDETNAME() const { return fkDETNAME; }

private:

  /** Constant members */
  const Int_t fkMAXHOSTS;                  //Constant
  const Int_t fkDEFAULTEVENTPORT;          //Constant
  //   const Int_t fkMAXBINVALUE;          //Constant
  //  const Int_t fkHIGHGAIN;          //Constant
  //  const Int_t fkLOWGAIN;          //Constant

  //  const Int_t fkALTROMAXSAMPLES; /**<The maximum number of samples of the ALTRO*/
  //   const Int_t fkALTROMAXPRESAMPLES;           //Constant

  const Int_t fkNZROWSRCU; /**<Number of rows per module*/ 
  const Int_t fkNXCOLUMNSRCU;          //Constant
  const Int_t fkNZROWSMOD;  /**<Number of rows per module*/ 
  const Int_t fkNXCOLUMNSMOD;  /**<Number of columns per module*/ 
  //   const Int_t fkNGAINS;  /**<Number of gains per ALTRO channel*/
  const Int_t fkNDATATYPES;          //Constant

  const Int_t fkPFMAXPATHLENGTH;          //Constant

  const Int_t fkPFDEFAULTNSAMPLES;          //Constant
  const Int_t fkPFDEFAULTSTARTINDEX;          //Constant

  const Double_t fkDEFAULTTAU; /**<Assume that the signal rise time of the altrp pulses is 2 us (nominal value of the electronics)*/
  const Int_t fkDEFAULTFS; /**<Assume that the signal is samples with 10 MHZ samle rate*/

  const Int_t fkMODULE0;          //Constant
  const Int_t fkMODULE1;          //Constant
  const Int_t fkMODULE2;          //Constant
  const Int_t fkMODULE3;          //Constant
  const Int_t fkMODULE4;          //Constant
         
  const Int_t fkCSPSPERFEE; //Constant
  const Int_t fkRCU0; //Constant
  const Int_t fkRCU1; //Constant
  const Int_t fkRCU2; //Constant
  const Int_t fkRCU3; //Constant

  const Int_t fkZ0; //Constant
  const Int_t fkZ1; //Constant
  const Int_t fkX0; //Constant
  const Int_t fkX1; //Constant

  const Int_t fkNMODULES;   /**<Number of modules of the PHOS detector*/
  const Int_t fkNRCUS;   /**<Number of RCUs per Module*/
 
  const Int_t fkNRCUSPERMODULE;   /**<Number of RCUs per Module*/
  const Int_t fkNRCUSPERTOTAL; /**<Total number of RCUs for PHOS*/
  const Int_t fkNFEECS;  /**<Number of Frontend cards per branch*/
  const Int_t fkNALTROS;  /**<Number of ALTROs per frontend card*/
  const Int_t fkNALTROCHANNELS;  //Constant
  const Int_t fkNBRANCHES; //Constant

  const Float_t fkCELLSTEP;  //Constant
  const Int_t fkNRCUSPERSECTOR;  //Constant
  
  const Int_t fkDDLOFFSET;  //Constant

  const TString fkDETNAME;  //Constant
   
  ClassDef(AliHLTPHOSConstants, 1);

};

#endif
