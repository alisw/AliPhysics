//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTONLINECONFIGURATION_H
#define ALIHLTONLINECONFIGURATION_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

///  @file   AliHLTOnlineConfiguration.h
///  @author Matthias Richter
///  @author Lars Christian Raae
///  @date   
///  @brief  Description of the HLT online configuration

#include <vector>

#include "AliHLTLogging.h"
#include <TObject.h>
#include <TArrayC.h>
#include <TList.h>
#include <TXMLNode.h>

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
class AliHLTOnlineConfiguration : public TObject, public AliHLTLogging {
 public:
  /// standard constructor
  AliHLTOnlineConfiguration();
  /// destructor
  virtual ~AliHLTOnlineConfiguration();

  /// load configuration from file
  int LoadConfiguration(const char* filename);

  /// compress the xml buffer
  int Compress();

  /// uncompress the xml buffer
  int Uncompress();
  
  /// parse the xml buffer
  int Parse();
  
  /// get default chains (sources of HLTOutFormatter)
  const char* GetDefaultChains() const {return fDefaultChains.Data();}
  
  /// get component libraries
  TString GetComponentLibraries();

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
  /// buffer for XML configuration
  TArrayC fXMLBuffer;
  /// size of XML buffer
  UInt_t fXMLSize;
  /// list of parsed configuration entries
  TList fConfEntryList;
  /// default chains (sources of HLTOutFormatter)
  TString fDefaultChains;
  
  /**
   * Parse XML configuration.
   * @param node       XML root node of configuration
   * @return
   *   -EINVAL if unparsable or empty configuration
   *   -EPROTO if no configuration is loaded
   *   0 if any elements were successfully parsed
   */
  int ParseConfiguration(TXMLNode* node);

  /**
   * Parse XML configuration entry.
   * @param node       XML root node of entry
   * @param id         online component ID
   * @param type       online component type
   * @return
   *   -EINVAL if parsing error
   *   0 if entry was successfully parsed
   */
  int ParseEntry(TXMLNode* node, const char* id, const char* type);
  
  /**
   * Parse standard component configuration.
   * @param id         online component ID
   * @param type       online component type
   * @param cmd        online command
   * @param sources    component sources
   * @param nodes      online computing nodes
   * @return
   *   -EINVAL if parsing error
   *   0 if entry was successfully parsed
   */  
  int ParseStandardComponent(const char* id, const char* type, const char* cmd,
    TString& sources, TString& nodes);
  
  /**
   * Parse RORCPublisher configuration.
   * @param id         online component ID
   * @param type       online component type
   * @param cmd        online command
   * @param sources    component sources
   * @param nodes      online computing nodes
   * @return
   *   -EINVAL if parsing error
   *   0 if entry was successfully parsed
   */  
  int ParseRORCPublisher(const char* id, const char* type, const char* cmd,
    TString& sources, TString& nodes);

  /**
   * Parse TCPDumpSubscriber configuration.
   * @param id         online component ID
   * @param type       online component type
   * @param cmd        online command
   * @param sources    component sources
   * @param nodes      online computing nodes
   * @return
   *   -EINVAL if parsing error
   *   0 if entry was successfully parsed
   */  
  int ParseTCPDumpSubscriber(const char* id, const char* type, const char* cmd,
    TString& sources, TString& nodes);

  /**
   * Parse Relay configuration.
   * @param id         online component ID
   * @param type       online component type
   * @param cmd        online command
   * @param sources    component sources
   * @param nodes      online computing nodes
   * @return
   *   -EINVAL if parsing error
   *   0 if entry was successfully parsed
   */  
  int ParseRelay(const char* id, const char* type, const char* cmd,
    TString& sources, TString& nodes);

  /**
   * Parse HLTOutFormatter configuration.
   * @param id         online component ID
   * @param type       online component type
   * @param cmd        online command
   * @param sources    component sources
   * @param nodes      online computing nodes
   * @return
   *   -EINVAL if parsing error
   *   0 if entry was successfully parsed
   */  
  int ParseHLTOutFormatter(const char* id, const char* type, const char* cmd,
    TString& sources, TString& nodes);

  /**
   * Parse HLTOutWriterSubscriber configuration.
   * @param id         online component ID
   * @param type       online component type
   * @param cmd        online command
   * @param sources    component sources
   * @param nodes      online computing nodes
   * @return
   *   -EINVAL if parsing error
   *   0 if entry was successfully parsed
   */  
  int ParseHLTOutWriterSubscriber(const char* id, const char* type,
    const char* cmd, TString& sources, TString& nodes);

  ClassDef(AliHLTOnlineConfiguration, 1); // description of HLT online configuration
};

#endif
