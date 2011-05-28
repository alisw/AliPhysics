#ifndef ALIHLTTPCHWCFMERGERUNIT_H
#define ALIHLTTPCHWCFMERGERUNIT_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *


#include "AliHLTTPCHWCFDataTypes.h"


//  @class   AliHLTTPCHWCFMergerUnit
//  @author Sergey Gorbunov <sergey.gorbunov@fias.uni-frankfurt.de>
//  @author Torsten Alt <talt@cern.ch> 
//  @brief  Channel Merger unit of FPGA ClusterFinder Emulator for TPC
//  @brief  ( see AliHLTTPCHWCFEmulator class )
//  @note
//
class AliHLTTPCHWCFMergerUnit
{
 public:  

  /** standard constructor */
  AliHLTTPCHWCFMergerUnit();
  
  /** destructor */
  ~AliHLTTPCHWCFMergerUnit();

  /** set debug level */
  void SetDebugLevel( int val ){ fDebug = val; }

  /** do cluster deconvolution in pad direction */
  void SetDeconvolution( bool val ){ fDeconvolute = val; }

  /** bypass the merger */
  void SetByPassMerger( bool val  ){ fByPassMerger = val; }

  /** set max distance in timebins for matching fragments */
  void SetMatchDistance( unsigned int val ){ fMatchDistance = val; }
  
  /** update MeanTime value when merging fragments */
  void SetMatchTimeFollow( bool val ){ fMatchTimeFollow = val; }

 /** initialise */
  int Init();
  
  /** input stream of data */
  int InputStream( const AliHLTTPCHWCFClusterFragment *fragment );

  /** output stream of data */
  const AliHLTTPCHWCFClusterFragment *OutputStream();

 private: 

  /** copy constructor prohibited */
  AliHLTTPCHWCFMergerUnit(const AliHLTTPCHWCFMergerUnit&);
  /** assignment operator prohibited */
  AliHLTTPCHWCFMergerUnit& operator=(const AliHLTTPCHWCFMergerUnit&);  
  
  int  fDebug; // debug level
  unsigned int fMatchDistance; // max distance in timebins for matching fragments
  bool fMatchTimeFollow;    // update MeanTime value when merging fragments 
  bool fDeconvolute; // do cluster deconvolution in pad direction
  bool fByPassMerger;// bypass the merger 
  AliHLTTPCHWCFClusterFragment fInput; // current input
  AliHLTTPCHWCFClusterFragment fMemory[2][AliHLTTPCHWCFDefinitions::kMaxNTimeBins*2]; // memory for 2 channels
  AliHLTTPCHWCFClusterFragment *fSearchRange[2]; // search range array
  AliHLTTPCHWCFClusterFragment *fInsertRange[2]; // insert range array
  int fSearchStart[2];  // index of the first candidate in SR
  int fSearchEnd[2];    // index of end of SR
  int fInsertEnd[2];    // index of end of IR
  int fInsertRow[2];    // current row number in IR
  int fInsertPad[2];    // current pad number in IR
};

#endif
