//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTPHOSRAWANALYZERCOMPONENT_H
#define ALIHLTPHOSRAWANALYZERCOMPONENT_H

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


class AliHLTPHOSRawAnalyzerComponent: public AliHLTPHOSRcuProcessor
{
 public:
  AliHLTPHOSRawAnalyzerComponent();
  virtual ~AliHLTPHOSRawAnalyzerComponent();
  virtual int DoInit(int argc =0, const char** argv  = 0);
  virtual int Deinit();
  virtual const char* GetComponentID() = 0;
  virtual void GetInputDataTypes( vector <AliHLTComponentDataType>& list);
  virtual AliHLTComponentDataType GetOutputDataType();
  virtual int GetOutputDataTypes(AliHLTComponentDataTypeList& list);
  virtual void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);
  virtual AliHLTComponent* Spawn() = 0; 
 protected:
  AliHLTPHOSRawAnalyzer *fAnalyzerPtr;  /**<Pointer to an analyzer object used for raw data anlysis*/ 

  using AliHLTPHOSRcuProcessor::DoEvent;
  virtual int DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
		     AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
		     AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks ); 

 private:
  AliHLTPHOSRawAnalyzerComponent(const AliHLTPHOSRawAnalyzerComponent & );
  AliHLTPHOSRawAnalyzerComponent & operator = (const AliHLTPHOSRawAnalyzerComponent &);
  
  //AliAltroBunch* GetFirstBunch();
  void GetFirstBunch(AliAltroData *altrodata,  AliAltroBunch *altrobunch);
  void FillDataArray(UInt_t *data, const AliAltroData *altrodataptr, const int channel = -1);
  //  bool CheckBuffer();

  void Reset(AliHLTPHOSRcuCellEnergyDataStruct* cellDataPtr);
  void ResetDataPtr(int startindex = 0, int sampleCnt = 0);
  void SetBaselines(const char* baselineFile); 
  void SetSelectiveReadOutThresholds(const char* filepath, Int_t nSignams);
  
  Bool_t fSendChannelData;       /**<wether or not to send raw data from the component into shared memory*/


  //  Double_t fTmpChannelData[ALTROMAXSAMPLES];                        /**<temporary variable to store raw samples from a single altro channel*/
  UInt_t fTmpChannelData[ALTROMAXSAMPLES];                        /**<temporary variable to store raw samples from a single altro channel*/

  Double_t fMaxValues[NMODULES][NZROWSMOD][NXCOLUMNSMOD][NGAINS]; /**<array to store cell energies*/

  AliHLTPHOSRcuCellEnergyDataStruct* fOutPtr;  //comment
  AliHLTPHOSMapper *fMapperPtr; //Mapping from harware address to geometrical address
  AliHLTPHOSSanityInspector *fSanityInspectorPtr; //comment
  Bool_t fUseBaselineSubtraction; //comment
  Float_t fBaselines[NXCOLUMNSRCU][NZROWSRCU][NGAINS]; //comment

  AliAltroDecoder *fDecoderPtr;           // decoder for altro payload
  AliAltroData    *fAltroDataPtr;         // container for altro payload
  AliAltroBunch   *fAltroBunchPtr;        // container for altro bunches

  /** Are we pushing the raw data */
  Bool_t fDoPushRawData; //Added by OD

  /** Are we making digits? */
  Bool_t fDoMakeDigits;  //Added by OD

  /** The digit maker */
  AliHLTPHOSDigitMaker *fDigitMakerPtr;                    //! transient Added by OD

  /** The digit containerto use for digit making */
  AliHLTPHOSDigitContainerDataStruct *fDigitContainerPtr;  //! transient Added by OD

  /** Are we doing selective read out? */
  Bool_t fDoSelectiveReadOut;  //Added by OD

  /** Threshold for selective readout ( zero suppression threshold) */
  Float_t fSelectiveReadOutThresholds[NXCOLUMNSMOD][NZROWSMOD][NGAINS]; //Added by OD

  /** The selected HW addresses */
  AliHLTUInt16_t *fSelectedChannelsList;  //! transient Added by OD

  /** Should we check data size? */
  Bool_t fDoCheckDataSize;      //Added by OD

  unsigned long fNCorruptedBlocks;
  unsigned long fNOKBlocks;

};
#endif

