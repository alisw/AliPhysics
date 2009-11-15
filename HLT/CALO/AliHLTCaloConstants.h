//-*- Mode: C++ -*-
// $Id: AliHLTCALOConstants.h $

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               */

/// @file   AliHLCaloConstants.h
/// @author Svein Lindal
/// @date   
/// @brief  Class containing constants for PHOS and EMCAL
///         loaded libraries

#include "Rtypes.h"

#ifndef ALIHLTCALOCONSTANTS_H
#define ALIHLTCALOCONSTANTS_H

class TString;

class AliHLTCaloConstants
{
public:
  
  AliHLTCaloConstants(TString det);
  virtual ~AliHLTCaloConstants();
  
  Int_t GetMAXHOSTS() const {return fkMAXHOSTS; } 
  Int_t GetDEFAULTEVENTPORT() const { return fkDEFAULTEVENTPORT; }  
  Int_t GetMAXBINVALUE() const { return fkMAXBINVALUE; }	   	 
  Int_t GetHIGHGAIN() const { return fkHIGHGAIN; }		 
  Int_t GetLOWGAIN() const { return fkLOWGAIN; }		 
  Int_t GetALTROMAXSAMPLES() const { return fkALTROMAXSAMPLES; }	 
  Int_t GetALTROMAXPRESAMPLES() const { return fkALTROMAXPRESAMPLES; }   
  Int_t GetNZROWSRCU() const { return fkNZROWSRCU; }	   	 
  Int_t GetNXCOLUMNSRCU() const { return fkNXCOLUMNSRCU; }	 
  Int_t GetNZROWSMOD() const { return fkNZROWSMOD; }	   	 
  Int_t GetNXCOLUMNSMOD() const { return fkNXCOLUMNSMOD; }	 
  Int_t GetNGAINS() const{ return fkNGAINS; }		 
  Int_t GetNDATATYPES() const { return fkNDATATYPES; }	   	 
  Int_t GetPFMAXPATHLENGTH() const { return fkPFMAXPATHLENGTH; }	 
  Int_t GetPFDEFAULTNSAMPLES() const { return fkPFDEFAULTNSAMPLES; } 
  Int_t GetPFDEFAULTSTARTINDEX() const { return fkPFDEFAULTSTARTINDEX; }  
  Double_t GetDEFAULTTAU() const { return fkDEFAULTTAU; }	   	 
  Int_t GetDEFAULTFS() const { return fkDEFAULTFS; }	   	 
  Int_t GetMODULE0() const { return fkMODULE0; }		 
  Int_t GetMODULE1() const { return fkMODULE1; }		 
  Int_t GetMODULE2() const { return fkMODULE2; }  	 
  Int_t GetMODULE3() const { return fkMODULE3; }		 
  Int_t GetMODULE4() const { return fkMODULE4; }		 
  Int_t GetCSPSPERFEE() const { return fkCSPSPERFEE; }	   	 
  Int_t GetRCU0() const { return fkRCU0; }		 
  Int_t GetRCU1() const { return fkRCU1; }		 
  Int_t GetRCU2() const { return fkRCU2; }		 
  Int_t GetRCU3() const { return fkRCU3; }		 
  Int_t GetZ0() const { return fkZ0; }		   	 
  Int_t GetZ1() const { return fkZ1; }		   	 
  Int_t GetX0() const { return fkX0; }		   	 
  Int_t GetX1() const { return fkX1; }		   	 
  Int_t GetNMODULES() const { return fkNMODULES; }		 
  Int_t GetNRCUS() const { return fkNRCUS; }		 
  Int_t GetNRCUSPERMODULE() const { return fkNRCUSPERMODULE; }	 
  Int_t GetNRCUSPERTOTAL() const { return fkNRCUSPERTOTAL; }	 
  Int_t GetNFEECS() const { return fkNFEECS; }		 
  Int_t GetNALTROS() const { return fkNALTROS; }		 
  Int_t GetNALTROCHANNELS() const { return fkNALTROCHANNELS; }	 
  Int_t GetNBRANCHES() const { return fkNBRANCHES; }	   	
  Float_t GetCELLSTEP() const {return fkCELLSTEP; }
  Int_t GetNRCUSPERSECTOR() const {return fkNRCUSPERSECTOR; }


private:

  AliHLTCaloConstants();

  const Int_t fkMAXHOSTS; // soon to be obsolete
  const Int_t fkDEFAULTEVENTPORT; // soon to be obsolete
  const Int_t fkMAXBINVALUE; // 1023;
  const Int_t fkHIGHGAIN; //   1;
  const Int_t fkLOWGAIN; //   0;
  const Int_t fkALTROMAXSAMPLES; // 1008;                           /**<The maximum number of samples of the ALTRO*/
  const Int_t fkALTROMAXPRESAMPLES; // 15;        
  const Int_t fkNZROWSRCU; //   56;                    /**<Number of rows per module*/       
  const Int_t fkNXCOLUMNSRCU; //   16; 
  const Int_t fkNZROWSMOD; //  48;            /**<Number of rows per module*/       
  const Int_t fkNXCOLUMNSMOD; //  24;            /**<Number of columns per module*/ 
  const Int_t fkNGAINS; //   2;                             /**<Number of gains per ALTRO channel*/
  const Int_t fkNDATATYPES; //   10;    
  const Int_t fkPFMAXPATHLENGTH; // 256;
  const Int_t fkPFDEFAULTNSAMPLES; // 70;
  const Int_t fkPFDEFAULTSTARTINDEX; // 0;
  const Double_t fkDEFAULTTAU; // 0.2;                /**<Assume that the signal rise time of the altrp pulses is 2 us (nominal value of the electronics)*/
  const Int_t fkDEFAULTFS; // 10;   /**<Assume that the signal is samples with 10 MHZ samle rate*/
  const Int_t fkMODULE0; // 0;
  const Int_t fkMODULE1; // 1;
  const Int_t fkMODULE2; // 2;
  const Int_t fkMODULE3; // 3;
  const Int_t fkMODULE4; // 4;
  const Int_t fkCSPSPERFEE; // 32;
  const Int_t fkRCU0; // 0;
  const Int_t fkRCU1; // 1;
  const Int_t fkRCU2; // 2;
  const Int_t fkRCU3; // 3;
  const Int_t fkZ0; // 0;
  const Int_t fkZ1; // 1;
  const Int_t fkX0; // 0;
  const Int_t fkX1; // 1;
  const Int_t fkNMODULES; //      13;                            /**<Number of modules of the EMCAL detector*/
  const Int_t fkNRCUS; //      4;                             /**<Number of RCUs per Module*/
  const Int_t fkNRCUSPERMODULE; //  2;                            /**<Number of RCUs per Module*/
  const Int_t fkNRCUSPERTOTAL; //  NMODULES*NRCUSPERMODULE;        /**<Total number of RCUs for EMCAL*/
  const Int_t fkNFEECS; //  9;                             /**<Number of frontend cards per branch*/
  const Int_t fkNALTROS; //   4;                            /**<Number of ALTROs per fkrontend card*/
  const Int_t fkNALTROCHANNELS; //  16;
  const Int_t fkNBRANCHES; //   2;      
  const Float_t fkCELLSTEP;   // Obsolete variable? Called in Calomapper!
  const Int_t fkNRCUSPERSECTOR; // 4;

};
#endif
