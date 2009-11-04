//-*- Mode: C++ -*-
// $Id: AliHLTCALOConstant.h 34622 2009-09-04 13:22:01Z odjuvsla $

/**************************************************************************
 * This file is property of and copyright by the Experimental Nuclear     *
 * Physics Group, Dep. of Physics                                         *
 * University of Oslo, Norway, 2009                                       *
 *                                                                        * 
 * Author: Svein Lindal, slindal@fys.uio.no for the ALICE PHOS Project.  *
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

#ifndef ALIHLTCALOCONSTANTS_H
#define ALIHLTCALOCONSTANTS_H


class AliHLTCaloConstants
{
public:
  AliHLTCaloConstants();
  virtual ~AliHLTCaloConstants();

  virtual int GetMAXHOSTS() = 0;
  virtual int GetDEFAULTEVENTPORT() = 0; 
  virtual int GetMAXBINVALUE() = 0;
  virtual int GetHIGHGAIN() = 0;
  virtual int GetLOWGAIN() = 0;
  virtual int GetALTROMAXSAMPLES() = 0;
  virtual int GetALTROMAXPRESAMPLES() = 0;
  virtual int GetNZROWSRCU() = 0;
  virtual int GetNXCOLUMNSRCU() = 0;
  virtual int GetNZROWSMOD() = 0;
  virtual int GetNXCOLUMNSMOD() = 0;
  virtual int GetNGAINS() = 0;
  virtual int GetNDATATYPES() = 0;
  virtual int GetPFMAXPATHLENGTH() = 0;
  virtual int GetPFDEFAULTNSAMPLES() = 0;
  virtual int GetPFDEFAULTSTARTINDEX() = 0;
  virtual int GetDEFAULTTAU() = 0;
  virtual int GetDEFAULTFS() = 0;
  virtual int GetMODULE0() = 0;
  virtual int GetMODULE1() = 0;
  virtual int GetMODULE2() = 0;     
  virtual int GetMODULE3() = 0;
  virtual int GetMODULE4() = 0;
  virtual int GetCSPSPERFEE() = 0;
  virtual int GetRCU0() = 0;
  virtual int GetRCU1() = 0;
  virtual int GetRCU2() = 0;
  virtual int GetRCU3() = 0;
  virtual int GetZ0() = 0;
  virtual int GetZ1() = 0;
  virtual int GetX0() = 0;
  virtual int GetX1() = 0;
  virtual int GetNMODULES() = 0;
  virtual int GetNRCUS() = 0;
  virtual int GetNRCUSPERMODULE() = 0;
  virtual int GetNRCUSPERTOTAL() = 0;
  virtual int GetNFEECS() = 0;
  virtual int GetNALTROS() = 0;
  virtual int GetNALTROCHANNELS() = 0;
  virtual int GetNBRANCHES() = 0;
  virtual float GetCELLSTEP() = 0;

// #ifndef __CINT__
//   const unsigned char PFVECTORDIR[] = "/HLT/PHOS/PFVectors";
// #endif




};

#endif
