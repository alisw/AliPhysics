#ifndef ALIHLTPHOSRCUPROCESSOR_H
#define ALIHLTPHOSRCUPROCESSOR_H

/**************************************************************************
 * This file is property of and copyright by the Experimental Nuclear     *
 * Physics Group, Dep. of Physics                                         *
 * University of Oslo, Norway, 2007                                       *
 *                                                                        *
 * Author: Per Thomas Hille <perthi@fys.uio.no> for the ALICE HLT Project.*
 * Contributors are mentioned in the code where appropriate.              *
 * Please report bugs to perthi@fys.uio.no                                *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliHLTPHOSProcessor.h"


class  AliHLTPHOSRcuProcessor : public AliHLTPHOSProcessor
{
 public:
  AliHLTPHOSRcuProcessor();
  virtual ~AliHLTPHOSRcuProcessor();
  const AliHLTUInt16_t  GetEquippmentID() const;
  int ScanArguments(int argc, const char** argv);
  void SetEquippmentID(AliHLTUInt16_t id);
  void SetCoordinates(AliHLTUInt16_t equippmentID);
  const AliHLTUInt16_t fkEquippmentID;  /**<Equippment ID as defined by ALICE*/
  AliHLTUInt8_t  fRcuX;                 /**<X position of RCU the data from this Equippment comes from (0 or 1)*/
  AliHLTUInt8_t  fRcuZ;                 /**<Z position of RCU the data from this Equippment comes from (0 or 1)*/
  AliHLTUInt8_t  fRcuZOffset;           /**<offset in therms of towers in the Z direction relative to the module*/ 
  AliHLTUInt8_t  fRcuXOffset;           /**<offset in therms of towers in the X direction relative to the module*/
  Bool_t fIsSetEquippmentID;            /**<wether or not the EquippmentID is set*/
};

#endif


