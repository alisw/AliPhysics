//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTONLINECONFIGURATION_H
#define ALIHLTONLINECONFIGURATION_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

//  @file   AliHLTOnlineConfiguration.h
//  @author 
//  @date   
//  @brief  
//  @note   

#include "TObject.h"
#include "TArrayC.h"

/**
 * @class AliHLTOnlineConfiguration
 * @brief Description of the HLT online configuration
 * This class wraps an online HLT configuration file for storage in a
 * CDB object.
 *
 * Ideas:
 * The class provides easy access to the XML tags of the configuration through
 * an XML parser.
 *
 * The online configuration is translated into a tree like configuration
 * structure. The Draw function can be used to create a graph.
 *
 * The xml configuration is loaded into an internal TArrayC buffer, which
 * is automatically compressed depending on the compression level used for
 * root file storage. Later extension will be the implementation of a custom
 * Streamer function implementing the most efficient compression.
 */
class AliHLTOnlineConfiguration : public TObject {
 public:
  /// standard constructor
  AliHLTOnlineConfiguration();
  /// destructor
  virtual ~AliHLTOnlineConfiguration();

  /// load configuration from file
  int LoadConfiguration(const char* filename);

  /// compress the xml buffer
  int Compress();

  /// compress the xml buffer
  int Uncompress();

  /// overloaded from TObject, print info
  virtual void        Print(const char* options) const;

  /// overloaded from TObject, more crude data dump
  virtual void        Dump() const;

  /// overloaded from TObject, clear object
  virtual void        Clear(Option_t * option="");

  /// overloaded from TObject, clone object
  virtual TObject    *Clone(const char *newname="") const;

  /// overloaded from TObject, copy object
  virtual void        Copy(TObject &object) const;

  /// overloaded from TObject, draw graph of the configuration
  virtual void        Draw(Option_t *option="");

  /// custom status bits of TObject
  /// bit 14 to 23 can be freely used
  /// use functions SetBit, ResetBit, TestBit
  enum {
    kLoaded          = BIT(14),   // xml buffer is loaded
    kCompressed      = BIT(15),   // xml buffer is compressed
    kParsed          = BIT(16),   // already parsed
  };

 private:
  TArrayC fXMLBuffer; // buffer for XML configuration
  UInt_t fXMLSize;

  ClassDef(AliHLTOnlineConfiguration, 1); // description of HLT online configuration
};

#endif
