//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTT0CALIBOBJECT_H
#define ALIHLTT0CALIBOBJECT_H

//* This file is property of and copyright by the                          * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               */

//  @file   AliHLTT0CalibrationComponent.h
//  @author Alla
//  @date   2014-06-20
//  @brief  A sample calibration component for the HLT.
//  
//all T0 calibration params in 1 vector: 0-23 mean CFD
                          //                                       24-47 diff CFD
                          //                                       48 T0AC shift
                          //                                       49 T0A shift  
                          //                                       50 T0C shift
#include "TObject.h"

class AliHLTT0CalibObject : public TObject {
public:
  AliHLTT0CalibObject();
  //virtual ~AliHLTT0CalibObject();
  AliHLTT0CalibObject (const AliHLTT0CalibObject  &o);   
   Float_t  T0CalibParamsCFDMean(int ipmt) {return fT0CalibParams[ipmt];}  
  Float_t  T0CalibParamsCFDDiff(int ipmt) {return fT0CalibParams[ipmt+24];}  
  Float_t  T0CalibParamsT0shift(int icase) {return fT0CalibParams[icase+48];}
  void SetT0calibParams(int i, Float_t par) { fT0CalibParams[i]=par;}
private:

  Float_t  fT0CalibParams[52]; 
  AliHLTT0CalibObject& operator= (const AliHLTT0CalibObject &) { return *this;}

 
  ClassDef(AliHLTT0CalibObject,1)  // HLT T0 AliT0CalibObject (Header)
};

#endif 

