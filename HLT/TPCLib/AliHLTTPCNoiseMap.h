// XEmacs -*-C++-*-
// $Id$

#ifndef AliHLTTPCNOISEMAP_H
#define AliHLTTPCNOISEMAP_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTTPCNoiseMap.h
    @author Kalliopi Kanaki
    @date   
    @brief Class for reading the noise map from HCDB.
*/

#include "AliHLTLogging.h"
class AliTPCCalPad;

/** 
 * @class AliHLTTPCNoiseMap
 *
 * This singleton class enables the reading 
 * of the noise map from the HCDB. There will only 
 * be a single instance called. It returns the noise
 * map, as it is created by the offline code, i.e.
 * in global coordinates. 
 *
 * The use is simple:
 * AliHLTTPCNoiseMap *nm = AliHLTTPCNoiseMap::Instance();
 * AliTPCCalPad *noisePad = nm->ReadNoiseMap();
 *
 * @ingroup alihlt_tpc
 */  
 
  
class AliHLTTPCNoiseMap : public AliHLTLogging {
public:
  
  /** returns a pointer to the sole instance */
  static AliHLTTPCNoiseMap* Instance();
  
  /** method to retrieve the noise map from HCDB */
  AliTPCCalPad* ReadNoiseMap(Int_t runNo);

private:
  /** standard constructor prohibited */
  AliHLTTPCNoiseMap();
  /** copy constructor prohibited */
  AliHLTTPCNoiseMap(const AliHLTTPCNoiseMap&);
  /** assignment operator prohibited */
  AliHLTTPCNoiseMap& operator=(const AliHLTTPCNoiseMap&);
  
  /** pointer to sole instance */
  static AliHLTTPCNoiseMap *pNoiseMapInstance;

 ClassDef(AliHLTTPCNoiseMap, 0)
};
#endif // AliHLTTPCNOISEMAP_H
