#ifndef ALIHLTPHOSRAWANALYZERCOMPONENT_H
#define ALIHLTPHOSRAWANALYZERCOMPONENT_H

/* Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */

//
//Base class for PHOS HLT raw data analysis components
// see cxx file for more details

#include "AliHLTProcessor.h"
#include "AliHLTPHOSDefinitions.h"
#include "AliHLTPHOSCommonDefs.h"


class AliRawReaderMemory;
class AliCaloRawStream;
class AliHLTPHOSRawAnalyzer;
class AliHLTPHOSRcuCellEnergyDataStruct;
class AliHLTPHOSRcuChannelDataStruct;

class AliHLTPHOSRawAnalyzerComponent: public AliHLTProcessor
{
 public:
  AliHLTPHOSRawAnalyzerComponent();
  virtual ~AliHLTPHOSRawAnalyzerComponent();
  AliHLTPHOSRawAnalyzerComponent(const AliHLTPHOSRawAnalyzerComponent & );
  AliHLTPHOSRawAnalyzerComponent & operator = (const AliHLTPHOSRawAnalyzerComponent &)
   {
      return *this;
   };
  virtual int DoInit(int argc =0, const char** argv  = 0);
  virtual int Deinit();
  virtual int DoDeinit();
  void DumpData(int gain =0) const;
  void DumpChannelData(Double_t *data =0) const; 
  void SetEquippmentID(AliHLTUInt16_t id =0);
  const AliHLTUInt16_t  GetEquippmentID() const;
  void SetCoordinates(AliHLTUInt16_t equippmentID =0);
  virtual const char* GetComponentID() = 0;
  //  virtual void GetInputDataTypes(std::vector<AliHLTComponentDataType, std::allocator<AliHLTComponentDataType> >&);
  virtual void GetInputDataTypes( std::vector <AliHLTComponentDataType>& list);

  virtual AliHLTComponentDataType GetOutputDataType();
  virtual void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);
  virtual AliHLTComponent* Spawn() = 0; 
 protected:
  AliHLTPHOSRawAnalyzer *fAnalyzerPtr;  /**<Pointer to an analyzer object used for raw data anlysis*/ 
 private:
  virtual int DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
		     AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
		     AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks ); 
  void Reset();
  void ResetDataPtr(int startindex = 0, int sampleCnt = 0);
  static int fgEventCount;       /**<Global event counter for this component*/
  const AliHLTUInt16_t fkEquippmentID;  /**<Equippment ID as defined by ALICE*/
  AliHLTUInt8_t  fModuleID;      /**<ID of the module this component read data from (0-4)*/
  AliHLTUInt8_t  fRcuX;          /**<X position of RCU the data from this Equippment comes from (0 or 1)*/
  AliHLTUInt8_t  fRcuZ;          /**<Z position of RCU the data from this Equippment comes from (0 or 1)*/
  AliHLTUInt8_t  fRcuZOffset;    /**<offset in therms of towers in the Z direction relative to the module*/ 
  AliHLTUInt8_t  fRcuXOffset;    /**<offset in therms of towers in the X direction relative to the module*/ 
  Bool_t fPrintInfo;             /**<wether or not to print debugg info to std out*/
  Bool_t fSendChannelData;       /**<wether or not to send raw data from the component into shared memory*/
  Double_t fTmpChannelData[ALTRO_MAX_SAMPLES];                        /**<temporary variable to store raw samples from a single altro channel*/
  Double_t fMaxValues[N_MODULES][N_ROWS_MOD][N_COLUMNS_MOD][N_GAINS]; /**<array to store cell energies*/
  int fPrintInfoFrequncy;                             /**<Defines the update frequency for information printet to std out*/
  AliCaloRawStream *fPHOSRawStream;                   /**<Streamer for PHOS raw data, used by fPHOSRawMemory reader*/ 
  AliRawReaderMemory *fRawMemoryReader;               /**<Decoder to read PHOS raw data on the altro format*/  
  AliHLTPHOSRcuCellEnergyDataStruct* fOutPtr;         /**<Pointer to outputbuffer to write results from the component into shared memory*/
  static const AliHLTComponentDataType fgkInputDataTypes[]; /**<List of  datatypes that can be given to this component*/
};
#endif

