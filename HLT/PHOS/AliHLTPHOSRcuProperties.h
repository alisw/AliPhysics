#ifndef ALIHLTPHOSRCUPROPERTIES_HPP
#define ALIHLTPHOSRCUPROPERTIES_HPP

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

#include "AliHLTDataTypes.h"
#include "AliHLTLogging.h"
#include "AliHLTPHOSConstants.h"

using namespace PhosHLTConst;

//class  AliHLTPHOSRcuProperties : public AliHLTLogging
class  AliHLTPHOSRcuProperties
{
public:
  AliHLTPHOSRcuProperties();
  virtual ~AliHLTPHOSRcuProperties();

  const AliHLTUInt16_t  GetEquippmentID() const;
  const int GetRCUID() const;
  virtual int ScanArguments(int argc, const char** argv);
  void SetEquippmentID(AliHLTUInt16_t id);
  void SetCoordinates(AliHLTUInt16_t equippmentID);

  
  int fTest;

 protected:
  const AliHLTUInt16_t fkEquippmentID;  /**<Equippment ID as defined by ALICE*/
  AliHLTUInt8_t  fModID; /**<ID of the module this component read data from (0-4)*/
  int fRcuID;
  AliHLTUInt8_t  fRcuX;                 /**<X position of RCU the data from this Equippment comes from (0 or 1)*/
  AliHLTUInt8_t  fRcuZ;                 /**<Z position of RCU the data from this Equippment comes from (0 or 1)*/
  AliHLTUInt8_t  fRcuZOffset;           /**<offset in therms of towers in the Z direction relative to the module*/ 
  AliHLTUInt8_t  fRcuXOffset;           /**<offset in therms of towers in the X direction relative to the module*/
  Bool_t fPrintInfo;                    /**<wether or not to print debugg info to std out*/
  int fPrintInfoFrequncy;               /**<Defines the update frequency for information printet to std out*/

  Bool_t fIsSetEquippmentID;            /**<wether or not the EquippmentID is set*/
  AliHLTLogging fLog;

  bool fIsInitialized; 
 
  

};

#endif
