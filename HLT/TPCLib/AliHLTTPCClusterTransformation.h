// -*- Mode: C++ -*-
// $Id: AliHLTTPCClusterTransformation.h 40939 2010-05-04 15:35:58Z kkanaki $

#ifndef ALIHLTTPCCLUSTERTRANSFORMATION_H
#define ALIHLTTPCCLUSTERTRANSFORMATION_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTTPCClusterTransformation.h
    @author Kalliopi Kanaki, Sergey Gorbunov
    @date   
    @brief
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include"Rtypes.h"

class AliTPCParam;
class AliTPCTransform;

/**
 * @class AliHLTTPCClusterTransformation
 *
 * The class transforms internal TPC coordinates (pad,time) to XYZ.
 * Allnecessary calibration and alignment corrections are applied
 * 
 * @ingroup alihlt_tpc_components
 */

class AliHLTTPCClusterTransformation{
    
 public:

  /** standard constructor */    
  AliHLTTPCClusterTransformation();           
  /** destructor */
  virtual ~AliHLTTPCClusterTransformation();

  int  Init( double FieldBz, UInt_t TimeStamp );
  void SetCurrentTimeStamp( UInt_t TimeStamp );
  int  Transform( int Slice, int Row, float Pad, float Time, float XYZ[] );

  int  ReverseAlignment( float XYZ[], int slice, int padrow);
  void SetRotationMatrix(const Double_t *rot=NULL, bool bCalcAdjugate=false);
  bool CalcAdjugateRotation(bool bCheck=false);

  void Print(const char* option=NULL) const;

 protected:

  AliTPCTransform * fOfflineTransform;                             //! transient
  AliTPCParam     * fOfflineTPCParam;                                 //! transient
  Int_t fLastSector; // last sector
  Double_t fAliT[3]; // alignment - translation
  Double_t fAliR[9]; // alignment - rotation
  Double_t fAdjR[9]; // alignment - inverse rotation (adjugate)

 private:

  /** copy constructor prohibited */
  AliHLTTPCClusterTransformation(const AliHLTTPCClusterTransformation&);
  /** assignment operator prohibited */
  AliHLTTPCClusterTransformation& operator=(const AliHLTTPCClusterTransformation&);

  ClassDef(AliHLTTPCClusterTransformation, 0)
};

#endif
