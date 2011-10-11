#ifndef ALIHLTTRIGGERDETECTORGEOM_H
#define ALIHLTTRIGGERDETECTORGEOM_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTTriggerDetectorGeom.h
/// @author Oystein Djuvsland
/// @date   2009-10-08
/// @brief  HLT class describing simple geometry of (sub-)detectors.
///         Used for the AliHLTTriggerBarrelGeomMultiplicity classes


#include "TString.h"
#include <ostream>
#include "TObject.h"
#include "TString.h"
#include <ostream>
#include "TObject.h"

/**
 * @class  AliHLTTriggerDetectorGeom
 * HLT class describing simple geometry cuts of (sub-)detectors.
 * Used for the AliHLTTriggerBarrelGeomMultiplicity class
 *
 * \ingroup alihlt_trigger
 */

class AliHLTTriggerDetectorGeom : public TObject
{
public: 
  
  /** Default constructor */
  AliHLTTriggerDetectorGeom();   // Default constructor

  /** Destructor */
  virtual ~AliHLTTriggerDetectorGeom(); // Destructor

  /**
   * Set the minimum in eta
   * @param etaMin is the minumum
   */
  void SetEtaMin(Double_t etaMin) { fEtaMin = etaMin; } // Set min in eta

  /**
   * Set the maximum in eta
   * @param etaMax is the maxumum
   */
  void SetEtaMax(Double_t etaMax) { fEtaMax = etaMax; } // Set max in eta

  /**
   * Set the minimum in phi
   * @param phiMin is the minumum
   */
  void SetPhiMin(Double_t phiMin) { fPhiMin = phiMin; } // Set min in phi

  /**
   * Set the maximum in phi
   * @param phiMax is the maxumum
   */
  void SetPhiMax(Double_t phiMax) { fPhiMax = phiMax; } // Set max in phi

  /** 
   * Set the initial point describing the plane
   * @param point is the point
   */
  void SetInitialPoint(Double_t *point); // Set initial point

  /** 
   * Set the vector describing the plane
   * @param nVector is a the vector
   */
  void SetNormVector(Double_t *nVector); // Set normal vector

  /**
   * Set the name of the (sub-)detector
   * @param name is the name
   */
  void SetDetName(TString &name) { fName = name; } // Set name

  /** 
   * Get the minimum in eta (should be moved...)
   */
  Double_t EtaMin() { return fEtaMin; }

  /** 
   * Get the maximum in eta (should be moved...)
   */
  Double_t EtaMax() { return fEtaMax; }

  /** 
   * Get the minimum in phi (should be moved...)
   */
  Double_t PhiMin() { return fPhiMin; }

  /** 
   * Get the maximum in phi (should be moved...)
   */
  Double_t PhiMax() { return fPhiMax; }

  /** 
   * Get the initial point
   */
  void GetInitialPoint(Double_t *point);

  /**
   * Get the normal vector
   */
  void GetNormVector(Double_t *vec);
  
  /** 
   * Get the detector name
   */
  TString& DetName() { return fName; }

  /** 
   * Print the geometry
   */
  void PrintDetectorGeom(std::ostream &out);

  /** 
   * Abstract method to check if a point is in the 
   * acceptance of the geometry 
   */
  virtual Bool_t IsInDetector(Double_t point[3]) = 0;

protected:

  /** The minimum in eta */
  Double_t fEtaMin;   // The minimum in eta 

  /** The maximum in eta */
  Double_t fEtaMax;   // The maximum in eta 

  /** The minimum in phi */
  Double_t fPhiMin;   // The minimum in phi 

  /** The maximum in phi */
  Double_t fPhiMax;   // The maximum in phi 

  /** Name of the (sub-)detector */
  TString fName;      // Name of the (sub-)detector 
  
private: 

  /** 
   * The point which together with a normal vector 
   * defines the plane of the (sub-)detector
   */
  Double_t fInitalPoint[3]; // Point representing the plane

  /** 
   * The normal vector which together with a point
   * defines the plane of the (sub-)detector
   */
  Double_t fNormVector[3]; // Normal vector representing the plane

  ClassDef(AliHLTTriggerDetectorGeom, 1);

};

#endif
