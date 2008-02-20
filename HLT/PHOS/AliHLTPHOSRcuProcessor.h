#ifndef ALIHLTPHOSRCUPROCESSOR_H
#define ALIHLTPHOSRCUPROCESSOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


#include "AliHLTPHOSProcessor.h"


class  AliHLTPHOSRcuProcessor : public AliHLTPHOSProcessor
{
 public:
  AliHLTPHOSRcuProcessor();
  virtual ~AliHLTPHOSRcuProcessor();
  const AliHLTUInt16_t  GetEquippmentID() const;
  virtual int ScanArguments(int argc, const char** argv);
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


