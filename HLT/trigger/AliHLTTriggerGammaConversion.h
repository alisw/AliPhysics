//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTTRIGGERGAMMACONVERSION_H
#define ALIHLTTRIGGERGAMMACONVERSION_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTTriggerGammaConversion.h
/// @author Kenneth Aamodt
/// @date   2009-11-01
/// @brief  HLT trigger component for gamma conversions.
///         

#include "AliHLTTrigger.h"

class AliESDtrack;

/**
 * @class  AliHLTTriggerGammaConversion
 * HLT trigger component for gamma conversions
 * 
 * 
 * Triggers on gamma conversions wich satisfy cuts on invariant mass, radius, dca and Pt.
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b GammaConversionTrigger                             <br>
 * Library: \b libAliHLTTrigger.so                                        <br>
 * Input Data Types:  kAliHLTDataTypeESDObject                            <br>
 * Output Data Types: ::kAliHLTAnyDataType                                <br>
 *
 * <h2>Mandatory arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * <h2>Optional arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * <h2>Configuration:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * \li -max-invmass    <i> mass  </i> <br>
 *      invariant mass of the two gamma daughters
 * \li -minpt    <i> pt  </i> <br>
 *      required minimum pt for a trigger
 * \li -maxpt    <i> pt  </i> <br>
 *      required maximum pt for a trigger
 * \li -max-dca    <i> distance  </i> <br>
 *      dca between the two gamma daughters
  * \li -max-radius    <i> r  </i> <br>
 *      maximum radius from the collision point in xy-plane 

 *
 * By default, configuration is loaded from OCDB, can be overridden by
 * component arguments.
 *
 * <h2>Default CDB entries:</h2>
 * HLT/ConfigHLT/GammaConversionTrigger: TObjString storing the arguments
 *
 * <h2>Performance:</h2>
 * 
 *
 * <h2>Memory consumption:</h2>
 * 
 *
 * <h2>Output size:</h2>
 * 
 *
 * \ingroup alihlt_trigger_components
 */
class AliHLTTriggerGammaConversion : public AliHLTTrigger
{
 public:
  AliHLTTriggerGammaConversion();
  ~AliHLTTriggerGammaConversion();

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

  /// inherited from AliHLTComponent, scan one argument and
  /// its parameters
  int ScanConfigurationArgument(int argc, const char** argv);

 private:

  /// inherited from AliHLTTrigger: calculate the trigger
  virtual int DoTrigger();

  /// mass cut, maximum
  double fMaxInvMass; //! transient
  /// pt cut, maximum
  double fPtMax; //! transient
  /// pt cut, minimum
  double fPtMin; //! transient
  /// maximum dca to qualify as a gamma conversion
  double fMaxDca; //! transient
  /// maximum radius from collision point in xy-plane
  double fMaxR;

  /// number of reconstructed gammas
  Int_t fNReconstructedGammas; //! transient

  /// the default configuration entry for this component
  static const char* fgkOCDBEntry; //!transient

  ClassDef(AliHLTTriggerGammaConversion, 0)
};
#endif //ALIHLTTRIGGERGAMMACONVERSION_H
