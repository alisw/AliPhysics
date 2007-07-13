#ifndef ALIHLTPHOSPROCESSOR_H
#define ALIHLTPHOSPROCESSOR_H

#include "AliHLTProcessor.h"
#include "AliHLTPHOSConstants.h"
#include "AliHLTPHOSCommonDefs.h"
#include "TString.h"
#include "AliHLTPHOSDefinitions.h"

using namespace PhosHLTConst;

class AliHLTPHOSProcessor:public AliHLTProcessor
{
 public:
  AliHLTPHOSProcessor();
  virtual ~AliHLTPHOSProcessor();
  AliHLTPHOSProcessor(const AliHLTPHOSProcessor & );
  AliHLTPHOSProcessor & operator = (const AliHLTPHOSProcessor &)
   {
      return *this;
   };

  virtual int DoInit(int argc, const char** argv) = 0;
  virtual int Deinit() = 0;
  virtual const char* GetComponentID() = 0;
  const AliHLTUInt16_t  GetEquippmentID() const;
  virtual void GetInputDataTypes( std::vector <AliHLTComponentDataType>& list) =0;
  virtual AliHLTComponentDataType GetOutputDataType() =0;
  virtual void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier) =0;
  virtual AliHLTComponent* Spawn() = 0; 
 protected:
  int fPhosEventCount;                  /**<Global event counter for this component*/
  const AliHLTUInt16_t fkEquippmentID;  /**<Equippment ID as defined by ALICE*/
  AliHLTUInt8_t  fModuleID;             /**<ID of the module this component read data from (0-4)*/
  AliHLTUInt8_t  fRcuX;                 /**<X position of RCU the data from this Equippment comes from (0 or 1)*/
  AliHLTUInt8_t  fRcuZ;                 /**<Z position of RCU the data from this Equippment comes from (0 or 1)*/
  AliHLTUInt8_t  fRcuZOffset;           /**<offset in therms of towers in the Z direction relative to the module*/ 
  AliHLTUInt8_t  fRcuXOffset;           /**<offset in therms of towers in the X direction relative to the module*/ 
  Bool_t fPrintInfo;                    /**<wether or not to print debugg info to std out*/
  Bool_t fIsSetEquippmentID;            /**<wether or not the EquippmentID is set*/ 
  int fPrintInfoFrequncy;               /**<Defines the update frequency for information printet to std out*/
  void SetEquippmentID(AliHLTUInt16_t id);
  void SetCoordinates(AliHLTUInt16_t equippmentID);
  int ScanArguments(int argc, const char** argv);
  static const AliHLTComponentDataType fgkInputDataTypes[]; /**<List of  datatypes that can be given to this component*/
 private:
  

};


#endif
