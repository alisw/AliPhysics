//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTTPCQADATAMAKER_H
#define ALIHLTTPCQADATAMAKER_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTTPCQADataMaker.h
    @author Zhongbao Yin, Matthias Richter
    @date   2009-05-14
    @brief  Container for the HLT offline QA
*/

#include "AliHLTQADataMakerBase.h"

/**
 * @class AliHLTTPCQADataMaker
 * Steering class for HLT QA for reconstruction
 */
class AliHLTTPCQADataMaker: public AliHLTQADataMakerBase {

public:

  enum HESDsType_t {kMultiplicity=0, kMultiplicityFired, kNCls, 
		    kNClsFired, kPHLT, kPOffline, kPRatio, 
		    kPHLTFired, kPOfflineFired, kPRatioFired,
		    kPtHLT, kPtOffline, 
		    kPtHLTFired, kPtOfflineFired, 
		    kNClsPerTrkHLT, kNClsPerTrkOffline, 
		    kNClsPerTrkHLTFired, kNClsPerTrkOfflineFired, 
		    kPhiHLT, kPhiOffline,
		    kPhiHLTFired, kPhiOfflineFired,
		    kEtaHLT, kEtaOffline,
		    kEtaHLTFired, kEtaOfflineFired};  
  
  AliHLTTPCQADataMaker();
  virtual ~AliHLTTPCQADataMaker();

private:
  /** copy constructor prohibited */
  AliHLTTPCQADataMaker(const AliHLTTPCQADataMaker&);   
  /** assignment operator prohibited */
  AliHLTTPCQADataMaker& operator = (const AliHLTTPCQADataMaker&);

  virtual void StartOfDetectorCycle();
  virtual void EndOfDetectorCycle(AliQAv1::TASKINDEX_t, TObjArray** list);
  virtual void MakeRaws(AliRawReader * rawReader);
  virtual void InitESDs();
  virtual void MakeESDs(AliESDEvent * esd, AliESDEvent* hltesd);

  ClassDef(AliHLTTPCQADataMaker,0)  // HLT Quality Assurance Data Maker for reconstruction
};

#endif // ALIHLTTPCQADATAMAKER_H
