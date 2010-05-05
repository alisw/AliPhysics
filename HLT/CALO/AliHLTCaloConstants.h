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
#include "TString.h"

#ifndef ALIHLTCALOCONSTANTS_H
#define ALIHLTCALOCONSTANTS_H



#define  fgkALTROMAXSAMPLES   1008  /**<The maximum number of samples of the ALTRO*/ 
#define  fgkALTROMAXPRESAMPLES  15  //Constant 
#define  fgkNGAINS  2
#define  fgkHIGHGAIN  1
#define  fgkLOWGAIN  0
#define  fkMAXBINVALUE  1023 //Constant



class AliHLTCaloConstants
{

public:

  AliHLTCaloConstants();
  virtual ~AliHLTCaloConstants();
  virtual Int_t GetMAXHOSTS() const  { return fkMAXHOSTS; } ; 				       	
  virtual Int_t GetDEFAULTEVENTPORT()  const { return fkDEFAULTEVENTPORT; }; 			       
  static  Int_t GetALTROMAXSAMPLES()    { return fgkALTROMAXSAMPLES; }; 
  static  Int_t GetALTROMAXPRESAMPLES()  { return fgkALTROMAXPRESAMPLES; } ; 
  
  static  Int_t GetNGAINS()     {  return fgkNGAINS; }; 	
  static  Int_t GetHIGHGAIN()   { return fgkHIGHGAIN; };
  static  Int_t GetLOWGAIN()    { return fgkLOWGAIN; }; 
  static  Int_t GetMAXBINVALUE()  {return fkMAXBINVALUE; }; 	
  
  
  // virtual Int_t GetNZROWSMOD() const { return fgkNZROWSMOD;} 
  virtual Int_t GetNZROWSMOD() const = 0;

  virtual Int_t GetNZROWSRCU() const = 0;
  virtual Int_t GetNXCOLUMNSRCU() const = 0;
  virtual Int_t GetNXCOLUMNSMOD() const = 0; 			     	
  //  virtual Int_t GetNGAINS() const = 0; 				
  virtual Int_t GetNDATATYPES() const { return fkNDATATYPES; } ; 				       	
  virtual Int_t GetPFMAXPATHLENGTH() const = 0;
  virtual Int_t GetPFDEFAULTNSAMPLES() const = 0; 			       	
  virtual Int_t GetPFDEFAULTSTARTINDEX() const = 0; 				
  virtual Double_t GetDEFAULTTAU() const = 0; 					
  virtual Int_t GetDEFAULTFS() const = 0; 				       
  Int_t GetCSPSPERFEE() const {  return fkCSPSPERFEE; }; 			       		
 
  virtual Int_t GetNMODULES() const = 0; 					
  virtual Int_t GetNRCUS() const = 0; 						
  virtual Int_t GetNRCUSPERMODULE() const = 0; 				       	
  virtual Int_t GetNRCUSPERTOTAL() const = 0; 					
  virtual Int_t GetNFEECS() const = 0; 						
  Int_t GetNALTROS() {  return fkNALTROS; }; 					
  Int_t GetNALTROCHANNELS() const { return fkNALTROCHANNELS; }; 					
  Int_t GetNBRANCHES() const { return fkNBRANCHES; }; 					
 
   // EMCAL specific
  Float_t GetCELLSTEP() const { return fkCELLSTEP; }
  Float_t GetMAXCELLSTEPETA() const { return fkMAXCELLSTEPETA; }  //FR
  Float_t GetMINCELLSTEPETA() const { return fkMINCELLSTEPETA; }  //FR
  Float_t GetCELLSTEPPHI() const { return fkCELLSTEPPHI; }        //FR
  Float_t GetCELLHEIGHT() const { return fkCELLHEIGHT; }        //FR
  Float_t GetCELLANGLE() const { return fkCELLANGLE; }        //FR
  Float_t GetRADLENGTH() const { return fkRADLENGTH; }        //FR
  Float_t GetCRITICENERGY() const { return fkCRITICENERGY; }        //FR
  Float_t GetCJ() const { return fkCJ;} //FR
  virtual Int_t GetDDLOFFSET() const = 0;
  TString GetDETNAME() { return fkDETNAME; };
  
  
protected:
  TString fkDETNAME;
  Float_t fkCELLSTEP; //Constant
  Float_t fkMAXCELLSTEPETA;
  Float_t fkMINCELLSTEPETA;
  Float_t fkCELLSTEPPHI;
  Float_t fkCELLHEIGHT;
  Float_t fkCELLANGLE;
  Float_t fkRADLENGTH;
  Float_t fkCRITICENERGY;
  Float_t fkCJ;
  Int_t fkDDLOFFSET;   //Constant
  
private:
  const Int_t fkNALTROS;  /**<Number of ALTROs per frontend card*/
  const Int_t fkNALTROCHANNELS;  //Constant
  const Int_t fkNBRANCHES; //Constant
  const Int_t fkCSPSPERFEE; //Constant
  const Int_t fkNDATATYPES;          //Constant
  const Int_t fkMAXHOSTS;                  //Constant
  const Int_t fkDEFAULTEVENTPORT;          //Constant
  ClassDef(AliHLTCaloConstants, 1);
};



#endif
