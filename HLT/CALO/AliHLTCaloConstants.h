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

#include "AliCaloConstants.h"


class AliHLTCaloConstants
{
public:
  AliHLTCaloConstants();
  virtual ~AliHLTCaloConstants();
  virtual void InitConstants() = 0; 

  // Common PHOS / EMCAL stuff
  static  Int_t GetALTROMAXSAMPLES()    { return ALTRO::ALTROMAXSAMPLES; };
  static  Int_t GetNGAINS()             { return ALTRO::NGAINS; };
  static  Int_t GetHIGHGAIN()           { return ALTRO::HIGHGAIN; };
  static  Int_t GetLOWGAIN()          { return ALTRO::LOWGAIN; };
  static  Int_t GetHGLGFACTOR()       { return CALO::HGLGFACTOR;}; //FR
  static  Int_t GetMAXBINVALUE()      { return ALTRO::MAXBINVALUE; };
  static  Int_t GetCSPSPERFEE()       { return CALO::CSPSPERFEE; };
  static  Int_t GetNALTROS()           { return ALTRO::NALTROS; };
  static  Int_t GetNALTROCHANNELS()    { return ALTRO::NALTROCHANNELS; };
  static  Int_t GetNBRANCHES()         { return CALO::NBRANCHES; }; 	
  //static  Int_t GetMAXHWADDRESSES() { return CALO::MAXHWADDRESSES; }
  static  Int_t GetMAXHWADDRESSES() { return PHOS::MAXHWADDR; }
  // Detector specific stuff
  // PHOS Only, bad move somewher else, PTH
  virtual Int_t GetNZROWSRCU() const  { return   56 ; } ; 
  virtual Int_t GetNXCOLUMNSRCU() const { return 16; } ;
  // END PHOS Only

  virtual Int_t GetNZROWSMOD() const      = 0; 
  virtual Int_t GetNXCOLUMNSMOD() const   = 0; 
  virtual Int_t GetNMODULES() const       = 0; 
  virtual Int_t GetNRCUSPERMODULE() const = 0;
  virtual Int_t GetNFEECS() const         = 0;
  
  Int_t GetDDLOFFSET() const { return fkDDLOFFSET; }
  Float_t GetCELLSTEP() const { return fkCELLSTEP; }
  
  //EMCAL specific, !! Move somewhere else, PTH
  Float_t GetCELLSTEPPHI() const { return fkCELLSTEPPHI; }        //FR
  Float_t GetCELLHEIGHT() const { return fkCELLHEIGHT; }        //FR
  Float_t GetCELLANGLE() const { return fkCELLANGLE; }        //FR
  Float_t GetRADLENGTH() const { return fkRADLENGTH; }        //FR
  Float_t GetCRITICENERGY() const { return fkCRITICENERGY; }        //FR
  Float_t GetCJ() const { return fkCJ;} //FR
  TString GetDETNAME() { return fkDETNAME; };

protected:
  //EMCAL specific, !! Move somewhere else, PTH
  // @todo: These variables should be declared constant, doesnt work right now
  // because the default copy contructor is called somewhere.
  
  Float_t fkCELLSTEPPHI;
  Float_t fkCELLHEIGHT;
  Float_t fkCELLANGLE;
  Float_t fkRADLENGTH;
  Float_t fkCRITICENERGY;
 
  Float_t fkCJ;
  Float_t fkCELLSTEP; //Constant

  Int_t fkDDLOFFSET;   //Constant
  TString fkDETNAME;

private:
  AliHLTCaloConstants(const AliHLTCaloConstants & );
  AliHLTCaloConstants & operator = (const AliHLTCaloConstants &); 

  ClassDef(AliHLTCaloConstants, 1)
};



#endif
