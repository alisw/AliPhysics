#ifndef ALIHLTTRIGGERDETECTORGEOMRECTANGLE_H
#define ALIHLTTRIGGERDETECTORGEOMRECTANGLE_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTTriggerDetectorGeomRectangle.h
/// @author Oystein Djuvsland
/// @date   2009-10-08
/// @brief  HLT class describing simple rectangular geometry of (sub-)detectors.
///         Used for the AliHLTTriggerBarrelGeomMultiplicity classes

#include "AliHLTTriggerDetectorGeom.h"

/**
 * @class  AliHLTTriggerDetectorGeomRectangle
 * HLT class describing simple rectangular geometry cuts of (sub-)detectors.
 * Used for the AliHLTTriggerBarrelGeomMultiplicity classes
 *
 * \ingroup alihlt_trigger
 */

#include "AliHLTTriggerDetectorGeom.h"

class AliHLTTriggerDetectorGeomRectangle : public AliHLTTriggerDetectorGeom
{
public: 

  /** Default constructor */
  AliHLTTriggerDetectorGeomRectangle();
  
  /** Default destructor */
  virtual ~AliHLTTriggerDetectorGeomRectangle();

  /**
   * Check if a point is in the detector geometry. 
   * @param point is the point in global coordinates (cm)
   * @return true if the point is in the geometry
   */
  Bool_t IsInDetector(Double_t point[3]);
  
  ClassDef(AliHLTTriggerDetectorGeomRectangle, 1);
};

#endif
