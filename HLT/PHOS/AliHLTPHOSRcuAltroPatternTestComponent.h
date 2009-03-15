//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTPHOSRCUALTRPATTERNTESTCOMPONENT_H
#define ALIHLTPHOSRCUALTRPATTERNTESTCOMPONENT_H 

/* Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice  */ 

#include "AliHLTPHOSRcuProcessor.h"





//
// Class for validation of PHOS
// Frontend electronics using playback
// of pedestal patterns from the altro
//


class AliHLTPHOSSharedMemoryInterface;
class AliHLTPHOSRcuAltroPatternTest;

class AliHLTPHOSRcuAltroPatternTestComponent:public AliHLTPHOSRcuProcessor
{
 public:
  AliHLTPHOSRcuAltroPatternTestComponent();
  virtual ~AliHLTPHOSRcuAltroPatternTestComponent();
  virtual int DoInit( int argc, const char** argv );
  virtual int Deinit();
  virtual int DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
		     AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
		     AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks );
  virtual void GetInputDataTypes( vector <AliHLTComponentDataType>&);
  virtual AliHLTComponentDataType GetOutputDataType();
  virtual void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);
  virtual AliHLTComponent* Spawn();
  virtual const char* GetComponentID();

 protected:
  using AliHLTPHOSRcuProcessor::DoEvent;

 private:
  AliHLTPHOSRcuAltroPatternTestComponent(const AliHLTPHOSRcuAltroPatternTestComponent &);
  AliHLTPHOSRcuAltroPatternTestComponent & operator = (const AliHLTPHOSRcuAltroPatternTestComponent &);
  AliHLTPHOSRcuAltroPatternTest *fPatternTestPtr;
  void ScanPatternFromFile(const char *filename, int *pattern, const int lengt = ALTROMAXSAMPLES) const;
  AliHLTPHOSSharedMemoryInterface *fShmPtr; // Interface to read altro channel data from shared memory
  unsigned long fNTotalPatterns;   /**<The total number of patterns analyzed*/
  unsigned long fNWrongPatterns;   /**<The total number of incorrect patterns found*/
  unsigned long fNTotalSamples;    /**<The total number of samples analyzed*/
  unsigned long fNWrongSamples;    /**<The total number of incorrect samples found*/
};

#endif
