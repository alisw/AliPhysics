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


namespace CaloHLTConst
{
  //Constants related to the ALTRO chip (Common to EMCAL / PHOS )
  const int MAXBINVALUE = 1023;
  const int NGAINS         =   2;    
  const int HIGHGAIN    =   1;
  const int LOWGAIN     =   0;
  const int ALTROMAXSAMPLES = 1008;    /**<The maximum number of samples of the ALTRO*/
  const int ALTROMAXPRESAMPLES = 15;        
  const int NALTROS        =   4;      /**<Number of ALTROs per frontend card*/
  const int NALTROCHANNELS =  16;

  //FEE constants common to PHOS EMCAL
  const int CSPSPERFEE    = 32;
  const int NBRANCHES      =   2;   
  
  namespace EmcalHLTConst
  {
    const int NZROWSMOD      =  48;   /**<Number of rows per module*/       
    const int NXCOLUMNSMOD   =  24;   /**<Number of columns per module*/ 
    const int NRCUSPERSECTOR = 4;
    const int NMODULES    =    10;    /**<Number of modules of the EMCAL detector*/
    const int NRCUSPERMODULE =  2 ;   /**<Number of RCUs per Module*/
    const int NFEECS         =  9; 
  };

  namespace PhosHLTConst
  {
    const int NZROWSMOD      =  56;   /**<Number of rows per module*/       
    const int NXCOLUMNSMOD   =  64;   /**<Number of columns per module*/ 
    const int NMODULES    =    5;     /**<Number of modules of the EMCAL detector*/
    const int NRCUSPERMODULE =  4 ;   /**<Number of RCUs per Module*/
    const int NFEECS         =  14;   /**<Number of Frontend cards per branch*/
  };
};


namespace CALO  =  CaloHLTConst; // just for easier notation
namespace EMCAL =  CaloHLTConst::EmcalHLTConst;
namespace PHOS  =  CaloHLTConst::PhosHLTConst;


class AliHLTCaloConstants
{
public:
  AliHLTCaloConstants();
  virtual ~AliHLTCaloConstants();
  virtual void InitConstants() = 0; 

  // Common PHOS / EMCAL stuff
  static  Int_t GetALTROMAXSAMPLES()    { return CALO::ALTROMAXSAMPLES; }; 
  static  Int_t GetNGAINS()             { return CALO::NGAINS; }; 	
  static  Int_t GetHIGHGAIN()           { return CALO::HIGHGAIN; };
  static  Int_t GetLOWGAIN()          { return CALO::LOWGAIN; }; 
  static  Int_t GetMAXBINVALUE()      { return CALO::MAXBINVALUE; }; 	
  static  Int_t GetCSPSPERFEE()       { return CALO::CSPSPERFEE; }; 			       		
  static  Int_t GetNALTROS()           { return CALO::NALTROS; }; 					
  static  Int_t GetNALTROCHANNELS()    { return CALO::NALTROCHANNELS; }; 					
  static  Int_t GetNBRANCHES()         { return CALO::NBRANCHES; }; 	

  // Detector specific stuff

 
  //  virtual Int_t GetNZROWSRCU() const      = 0; 
  // virtual Int_t GetNXCOLUMNSRCU() const   = 0;
  
  // PHOS Only, bad move somewher else, PTH
  virtual Int_t GetNZROWSRCU() const  { return   56 ; } ; 
  virtual Int_t GetNXCOLUMNSRCU() const { return 16; } ;
  // END PHOS Only

  virtual Int_t GetNZROWSMOD() const      = 0; 
  virtual Int_t GetNXCOLUMNSMOD() const   = 0; 
  virtual Int_t GetNMODULES() const       = 0; 
  //  virtual Int_t GetNRCUS() const          = 0;
  virtual Int_t GetNRCUSPERMODULE() const = 0;
  //  virtual Int_t GetNRCUSPERTOTAL() const  = 0;
  virtual Int_t GetNFEECS() const         = 0;
  
  Int_t GetDDLOFFSET() const { return fkDDLOFFSET; }
  Float_t GetCELLSTEP() const { return fkCELLSTEP; }
  
  //EMCAL specific, !! Move somewhere else, PTH
  Float_t GetMAXCELLSTEPETA() const { return fkMAXCELLSTEPETA; }  //FR
  Float_t GetMINCELLSTEPETA() const { return fkMINCELLSTEPETA; }  //FR
  Float_t GetCELLSTEPPHI() const { return fkCELLSTEPPHI; }        //FR
  Float_t GetCELLHEIGHT() const { return fkCELLHEIGHT; }        //FR
  Float_t GetCELLANGLE() const { return fkCELLANGLE; }        //FR
  Float_t GetRADLENGTH() const { return fkRADLENGTH; }        //FR
  Float_t GetCRITICENERGY() const { return fkCRITICENERGY; }        //FR
  Float_t GetCJ() const { return fkCJ;} //FR
  TString GetDETNAME() { return fkDETNAME; };

protected:
  //  TString fkDETNAME;
  //  Float_t fkCELLSTEP; //Constant
  //  Int_t fkDDLOFFSET;   //Constant
  
  //EMCAL specific, !! Move somewhere else, PTH
  // @todo: These variables should be declared constant, doesnt work right now
  // because the default copy contructor is called somewhere.
  Float_t fkMAXCELLSTEPETA;
  Float_t fkMINCELLSTEPETA;
  Float_t fkCELLSTEPPHI;
  Float_t fkCELLHEIGHT;
  Float_t fkCELLANGLE;
  Float_t fkRADLENGTH;
  Float_t fkCRITICENERGY;
  Float_t fkCELLSTEP; //Constant
  Float_t fkCJ;
  Int_t fkDDLOFFSET;   //Constant
  TString fkDETNAME;
private:
  
  ClassDef(AliHLTCaloConstants, 1)
};



#endif
