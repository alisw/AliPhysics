//-*- Mode: C++ -*-
// $Id: AliHLTEMCALConstants.h 35357 2009-10-08 13:24:38Z phille $

/**************************************************************************
 * This file is property of and copyright by the Experimental Nuclear     *
 * Physics Group, Dep. of Physics                                         *
 * University of Oslo, Norway, 2006                                       *
 *                                                                        * 
 * Author: Svein Lindal slindal@fys.uio.no for the ALICE DCS Project.  *
 * Contributors are mentioned in the code where appropriate.              *
 * Please report bugs to slindal@fys.uio.no                                * 
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#ifndef ALIHLTEMCALCONSTANTS_H
#define ALIHLTEMCALCONSTANTS_H

class AliHLTCaloConstants;


class AliHLTEMCALConstants : public AliHLTCaloConstants
{


public:
  AliHLTEMCALConstants();
  ~AliHLTEMCALConstants();
 
  // Int_t GetNZROWSRCU() const { return fkNZROWSRCU;}
  // Int_t GetNXCOLUMNSRCU() const { return fkNXCOLUMNSRCU;} 
  // Int_t GetNZROWSMOD() const { return fkNZROWSMOD;} 
  // Int_t GetNXCOLUMNSMOD() const { return fkNXCOLUMNSMOD;} 
  // Int_t GetNMODULES() const { return fkNMODULES;} 
  // Int_t GetNRCUS() const { return fkNRCUS;} 
  // Int_t GetNRCUSPERMODULE() const { return fkNRCUSPERMODULE;} 
  // Int_t GetNRCUSPERTOTAL() const { return fkNRCUSPERTOTAL;} 
  // Int_t GetNFEECS() const { return fkNFEECS;} 
  
  // Float_t GetCELLSTEP() const { return fkCELLSTEP; }
  // Float_t GetMAXCELLSTEPETA() const { return fkMAXCELLSTEPETA; }  //FR
  // Float_t GetMINCELLSTEPETA() const { return fkMINCELLSTEPETA; }  //FR
  // Float_t GetCELLSTEPPHI() const { return fkCELLSTEPPHI; }        //FR
  // Float_t GetCELLHEIGHT() const { return fkCELLHEIGHT; }        //FR
  // Float_t GetCELLANGLE() const { return fkCELLANGLE; }        //FR
  // Float_t GetRADLENGTH() const { return fkRADLENGTH; }        //FR
  // Float_t GetCRITICENERGY() const { return fkCRITICENERGY; }        //FR
  // Float_t GetCJ() const { return fkCJ;} //FR
  // Int_t GetDDLOFFSET() const { return fkDDLOFFSET; }

  //private:

  /** Constant members */
 
  // const Int_t fkNZROWSRCU; /**<Number of rows per module*/ 
  // const Int_t fkNXCOLUMNSRCU;//Constant
  // const Int_t fkNZROWSMOD;  /**<Number of rows per module*/ 
  // const Int_t fkNXCOLUMNSMOD;  /**<Number of columns per module*/ 
  // const Int_t fkNMODULES;   /**<Number of modules of the EMCAL detector*/
  // const Int_t fkNRCUS;   /**<Number of RCUs per Module*/
  // const Int_t fkNRCUSPERMODULE;   /**<Number of RCUs per Module*/
  // const Int_t fkNRCUSPERTOTAL; /**<Total number of RCUs for EMCAL*/
  // const Int_t fkNFEECS;  /**<Number of Frontend cards per branch*/

  ClassDef(AliHLTEMCALConstants, 1)

};
#endif
