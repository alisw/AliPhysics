#ifndef _AliHLTTriggerCTPRE_h_
#define _AliHLTTriggerCTPRE_h_
/// @file   AliHLTTriggerCTPRE.h
/// @author Mikolaj Krzewicki, mikolaj.krzewicki@cern.ch
/// @date   2015
/// @brief  HLT trigger CTP truggers based on a glob or regex

#include "AliHLTTrigger.h"
#include "TString.h"

class TPRegexp;

/**
 * @class  AliHLTTriggerCTPRE
 */
class AliHLTTriggerCTPRE : public AliHLTTrigger
{
 public:
  AliHLTTriggerCTPRE();
  virtual ~AliHLTTriggerCTPRE();

  /// inherited from AliHLTTrigger: name of this trigger
  virtual const char* GetTriggerName() const;
  /// inherited from AliHLTComponent: create an instance
  virtual AliHLTComponent* Spawn();

 protected:
  /// inherited from AliHLTComponent: handle the initialization
  int DoInit(int argc, const char** argv);

  /// inherited from AliHLTComponent: handle cleanup
  int DoDeinit();

  /// inherited from AliHLTComponent: handle re-configuration event
  int Reconfigure(const char* cdbEntry, const char* chainId);

  /// inherited from AliHLTComponent: handle dcs update event
  int ReadPreprocessorValues(const char* modules);

  /// Configure from CDB object, checking if AliHLTESDTrackCuts or TObjString
  int ConfigureFromCDBObject(TString cdbPath);

  /// inherited from AliHLTComponent, scan one argument and
  /// its parameters
  int ScanConfigurationArgument(int argc, const char** argv);

 private:
  /// copy constructor prohibited 
  AliHLTTriggerCTPRE (const AliHLTTriggerCTPRE&);
  
  /// assignment operator prohibited
  AliHLTTriggerCTPRE& operator=(const AliHLTTriggerCTPRE&);

  /// inherited from AliHLTTrigger: calculate the trigger
  virtual int DoTrigger();

  //helper for simple globbing
  Bool_t Globncmp(const char* triggerName, const char* glob, int triggerNameSize, int globSize );

  /// Name of the trigger
  TString              fName;               //! transient
  TString              fGlob;               //! glob for simple trig selection
  TPRegexp*            fRegexp;             //! a perl comp. regexp

  /// the default configuration entry for this component
  static const char*   fgkDefaultOCDBEntry; //!transient

  ClassDef(AliHLTTriggerCTPRE, 0)
};
#endif
