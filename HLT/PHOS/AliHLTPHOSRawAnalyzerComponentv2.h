#ifndef ALIHLTPHOSRAWANALYZERCOMPONENTV2_H
#define ALIHLTPHOSRAWANALYZERCOMPONENTV2_H

/* Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */
#include "AliHLTPHOSRcuProcessor.h"


class AliHLTPHOSRawAnalyzer;
class AliHLTPHOSRcuCellEnergyDataStruct;
class AliHLTPHOSMapper;
class AliHLTPHOSSanityInspector;
class AliAltroDecoder;      // decoder for altro payload
class AliAltroData;         // container for altro payload
class AliAltroBunch;        // container for altro bunches
class AliHLTPHOSDigitMaker;
class AliHLTPHOSDigitContainerDataStruct;


class AliHLTPHOSRawAnalyzerComponentv2 : public AliHLTPHOSRcuProcessor
{
 public:
  AliHLTPHOSRawAnalyzerComponentv2();
  virtual ~AliHLTPHOSRawAnalyzerComponentv2();
  virtual int DoInit(int argc =0, const char** argv  = 0);
  virtual int Deinit();
  virtual const char* GetComponentID() = 0;
  virtual void GetInputDataTypes( vector <AliHLTComponentDataType>& list);
  virtual AliHLTComponentDataType GetOutputDataType();
  virtual void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);
  virtual AliHLTComponent* Spawn() = 0; 

 protected:
  AliHLTPHOSRawAnalyzer *fAnalyzerPtr;  /**<Pointer to an analyzer object used for raw data anlysis*/ 

  using AliHLTPHOSRcuProcessor::DoEvent;
  virtual int DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
		     AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
		     AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks ); 

  virtual Int_t DoIt(const AliHLTComponentBlockData* iter, AliHLTUInt8_t* outputPtr, const AliHLTUInt32_t size, UInt_t& totSize);

 private:
  AliHLTPHOSRawAnalyzerComponentv2(const AliHLTPHOSRawAnalyzerComponentv2 & );
  AliHLTPHOSRawAnalyzerComponentv2 & operator = (const AliHLTPHOSRawAnalyzerComponentv2 &);

  AliHLTPHOSMapper *fMapperPtr; //Mapping from harware address to geometrical address
  AliHLTPHOSSanityInspector *fSanityInspectorPtr; //comment

  AliAltroDecoder *fDecoderPtr;           // decoder for altro payload
  AliAltroData    *fAltroDataPtr;         // container for altro payload
  AliAltroBunch   *fAltroBunchPtr;        // container for altro bunches

  unsigned long fNCorruptedBlocks;
  unsigned long fNOKBlocks;

  Short_t fAlgorithm;

};
#endif

