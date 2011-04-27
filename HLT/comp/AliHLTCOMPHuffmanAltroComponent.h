//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTCOMPHUFFMANALTROCOMPONENT_H
#define ALIHLTCOMPHUFFMANALTROCOMPONENT_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTCOMPHuffmanAltroComponent.h
/// @author Jenny Wagner, maintained by Matthias.Richter@cern.ch
/// @date   20-11-2007
/// @brief  The Huffman compressor component.
///

#include "AliHLTProcessor.h"

class AliHLTCOMPHuffmanAltro;
class AliHLTCOMPHuffmanData;

/**
 * @class AliHLTCOMPHuffmanAltroComponent
 * Implementation of the Huffman compressor component.
 * The component implements the interface methods of the @ref AliHLTProcessor.
 * The actual Huffman compression and decompression algorithms are implemented
 * in AliHLTCOMPHuffmanAltro.
 * The component can handle compressed and decompressed data of different formats.
 * Two components are registered, the HuffmanCompressor and the HuffmanDecompressor.
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b COMPHuffmanCompressor, \b COMPHuffmanDecompressor   <br>
 * Library: \b libAliHLTComp.so                                         <br>
 * Input Data Types:                                                    <br>
 * \li  COMPHuffmanDecompressor: AliHLTCompDefinitions::fgkDDLEncodedHuffmanAltroDataType
 * \li  COMPHuffmanCompressor: ::kAliHLTDataTypeDDLRaw                  <br>
 *
 * Output Data Types:                                                   <br>
 * \li  COMPHuffmanCompressor: AliHLTCompDefinitions::fgkDDLEncodedHuffmanAltroDataType
 * \li  COMPHuffmanDecompressor: ::kAliHLTDataTypeDDLRaw                <br>
 * 
 * <h2>Mandatory arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * \li -origin <i> detector </i> <br>
 *  set origin of data to specify Huffman code table to be read in (parameter transient)
 * \li -runnumber <i> decimal number </i> <br>
 *  set runnumber to specify Huffman code table to be read in (parameter transient)
 * \li -dataspec <i> 0xYYXXaabb </i> <br>
 *  set usual HLT dataspec (last slice, first slice, last patch, first patch)_Hexadezimal to specify Huffman code table to be read in
 * \li -trailerwords <i> decimal number </i> <br>
 *  set number of trailerwords of incoming data (ranging from 1 to 3)
 *
 * <h2>Optional arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * \li -tablepath <i> path to Huffman code table </i> <br>
 *  set path to read Huffman code table from root file, if no path is given, path is set to current path (parameter transient)
 *
 * @ingroup alihlt_comp_components
 */
class AliHLTCOMPHuffmanAltroComponent : public AliHLTProcessor {
 public:
  /**
   * constructor 
   * @param compression    whether to use the compressor or decompressor 
   **/
  AliHLTCOMPHuffmanAltroComponent(bool compression);
  /** destructor */
  virtual ~AliHLTCOMPHuffmanAltroComponent();
  
  // Public functions to implement AliHLTComponent's interface.
  // These functions are required for the registration process
  
  const char* GetComponentID();
  void GetInputDataTypes( vector<AliHLTComponentDataType>& list);
  AliHLTComponentDataType GetOutputDataType();
  virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );
  AliHLTComponent* Spawn();
  
  
 protected:
  
  // Protected functions to implement AliHLTComponent's interface.
  // These functions provide initialization as well as the actual processing
  // capabilities of the component. 
  
  int DoInit( int argc, const char** argv );
  int DoDeinit();
  int DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
	       AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
	       AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks );
  
 private:
  /** copy constructor prohibited */
  AliHLTCOMPHuffmanAltroComponent(const AliHLTCOMPHuffmanAltroComponent&);
  /** assignment operator prohibited */
  AliHLTCOMPHuffmanAltroComponent& operator=(const AliHLTCOMPHuffmanAltroComponent&);

  /** the entropy encoder object */
  AliHLTCOMPHuffmanAltro* fHuffmanCompressor;        // instance of Huffman compressor component   
  /** bool to decide whether to compress (TRUE) or to decompress (FALSE)  */                            
  Bool_t fCompressionSwitch;                         // mode choice for compressor instance from input line arguments
  /** bool to decide wheter to calibrate (create a new code table) (TRUE) or to (de)compress data (FALSE) */
  Bool_t fTrainingMode;                              // mode choice for training instance from input line arguments
  /** specification of origin of input data to load correct Huffman code table */
  TString fOrigin;   // gets input line argument to specify the origin of the Huffman table
  /** specification of run number to load correct Huffman code table */
  AliHLTUInt64_t fRunNumber;                           // gets input line argument to specify run type of Huffman Table
  /** data specifications (0xYYXXaabb) to load correct Huffman code table */
  AliHLTUInt64_t fDataSpec;                            // gets input line argument to specify the data of the Huffman table
  /** path to load Huffman code table from (if not explicitly given, table is taken from current path) */
  TString fTablePath;                                  // gets explicit path for Huffman table from command line
  /** number of NRCU trailer words of incoming data */
  Int_t fNrcuTrailerwords;                           // number of rcu trailer words
  /** pointer to Huffman code table for read in */
  AliHLTCOMPHuffmanData* fHuffmanData;               // instance of Huffman Data containing the code table
  
  //AliHLTUInt8_t fSlice;  // determine slice number from input event block
  //AliHLTUInt8_t fPatch;  // determine patch number from input event bloc
  
  ClassDef(AliHLTCOMPHuffmanAltroComponent, 0)
    };
#endif
