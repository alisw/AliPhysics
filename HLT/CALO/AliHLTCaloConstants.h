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

#ifndef ALIHLTCALOCONSTANTS_H
#define ALIHLTCALOCONSTANTS_H

class TString;

class AliHLTCaloConstants
{
public:
  
  AliHLTCaloConstants(TString det);
  virtual ~AliHLTCaloConstants();
    
  inline int GetMAXHOSTS() {return fMAXHOSTS; }
  inline int GetDEFAULTEVENTPORT() { return fDEFAULTEVENTPORT; }  
  inline int GetMAXBINVALUE() { return fMAXBINVALUE; }	   	 
  inline int GetHIGHGAIN() { return fHIGHGAIN; }		 
  inline int GetLOWGAIN() { return fLOWGAIN; }		 
  inline int GetALTROMAXSAMPLES() { return fALTROMAXSAMPLES; }	 
  inline int GetALTROMAXPRESAMPLES() { return fALTROMAXPRESAMPLES; }   
  inline int GetNZROWSRCU() { return fNZROWSRCU; }	   	 
  inline int GetNXCOLUMNSRCU() { return fNXCOLUMNSRCU; }	 
  inline int GetNZROWSMOD() { return fNZROWSMOD; }	   	 
  inline int GetNXCOLUMNSMOD() { return fNXCOLUMNSMOD; }	 
  inline int GetNGAINS() { return fNGAINS; }		 
  inline int GetNDATATYPES() { return fNDATATYPES; }	   	 
  inline int GetPFMAXPATHLENGTH() { return fPFMAXPATHLENGTH; }	 
  inline int GetPFDEFAULTNSAMPLES() { return fPFDEFAULTNSAMPLES; } 
  inline int GetPFDEFAULTSTARTINDEX() { return fPFDEFAULTSTARTINDEX; }  
  inline double GetDEFAULTTAU() { return fDEFAULTTAU; }	   	 
  inline int GetDEFAULTFS() { return fDEFAULTFS; }	   	 
  inline int GetMODULE0() { return fMODULE0; }		 
  inline int GetMODULE1() { return fMODULE1; }		 
  inline int GetMODULE2() { return fMODULE2; }  	 
  inline int GetMODULE3() { return fMODULE3; }		 
  inline int GetMODULE4() { return fMODULE4; }		 
  inline int GetCSPSPERFEE() { return fCSPSPERFEE; }	   	 
  inline int GetRCU0() { return fRCU0; }		 
  inline int GetRCU1() { return fRCU1; }		 
  inline int GetRCU2() { return fRCU2; }		 
  inline int GetRCU3() { return fRCU3; }		 
  inline int GetZ0() { return fZ0; }		   	 
  inline int GetZ1() { return fZ1; }		   	 
  inline int GetX0() { return fX0; }		   	 
  inline int GetX1() { return fX1; }		   	 
  inline int GetNMODULES() { return fNMODULES; }		 
  inline int GetNRCUS() { return fNRCUS; }		 
  inline int GetNRCUSPERMODULE() { return fNRCUSPERMODULE; }	 
  inline int GetNRCUSPERTOTAL() { return fNRCUSPERTOTAL; }	 
  inline int GetNFEECS() { return fNFEECS; }		 
  inline int GetNALTROS() { return fNALTROS; }		 
  inline int GetNALTROCHANNELS() { return fNALTROCHANNELS; }	 
  inline int GetNBRANCHES() { return fNBRANCHES; }	   	
  inline float GetCELLSTEP() {return fCELLSTEP; }
  inline int GetNRCUSPERSECTOR() {return fNRCUSPERSECTOR; }


private:

  AliHLTCaloConstants();

  const int fMAXHOSTS;
  const int fDEFAULTEVENTPORT;
  const int fMAXBINVALUE; // 1023;
  const int fHIGHGAIN; //   1;
  const int fLOWGAIN; //   0;
  const int fALTROMAXSAMPLES; // 1008;                           /**<The maximum number of samples of the ALTRO*/
  const int fALTROMAXPRESAMPLES; // 15;        
  const int fNZROWSRCU; //   56;                    /**<Number of rows per module*/       
  const int fNXCOLUMNSRCU; //   16; 
  const int fNZROWSMOD; //  48;            /**<Number of rows per module*/       
  const int fNXCOLUMNSMOD; //  24;            /**<Number of columns per module*/ 
  const int fNGAINS; //   2;                             /**<Number of gains per ALTRO channel*/
  const int fNDATATYPES; //   10;    
  const int fPFMAXPATHLENGTH; // 256;
  const int fPFDEFAULTNSAMPLES; // 70;
  const int fPFDEFAULTSTARTINDEX; // 0;
  const double fDEFAULTTAU; // 0.2;                /**<Assume that the signal rise time of the altrp pulses is 2 us (nominal value of the electronics)*/
  const int fDEFAULTFS; // 10;   /**<Assume that the signal is samples with 10 MHZ samle rate*/
  const int fMODULE0; // 0;
  const int fMODULE1; // 1;
  const int fMODULE2; // 2;
  const int fMODULE3; // 3;
  const int fMODULE4; // 4;
  const int fCSPSPERFEE; // 32;
  const int fRCU0; // 0;
  const int fRCU1; // 1;
  const int fRCU2; // 2;
  const int fRCU3; // 3;
  const int fZ0; // 0;
  const int fZ1; // 1;
  const int fX0; // 0;
  const int fX1; // 1;
  const int fNMODULES; //      13;                            /**<Number of modules of the EMCAL detector*/
  const int fNRCUS; //      4;                             /**<Number of RCUs per Module*/
  const int fNRCUSPERMODULE; //  2;                            /**<Number of RCUs per Module*/
  const int fNRCUSPERTOTAL; //  NMODULES*NRCUSPERMODULE;        /**<Total number of RCUs for EMCAL*/
  const int fNFEECS; //  9;                             /**<Number of Frontend cards per branch*/
  const int fNALTROS; //   4;                            /**<Number of ALTROs per frontend card*/
  const int fNALTROCHANNELS; //  16;
  const int fNBRANCHES; //   2;      
  const float fCELLSTEP;   // Obsolete variable? Called in Calomapper!
  const int fNRCUSPERSECTOR; // 4;

};
#endif
