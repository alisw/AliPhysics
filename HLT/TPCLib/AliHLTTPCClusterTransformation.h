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

#include "Rtypes.h"
#include "TString.h"
#include "AliHLTTPCGeometry.h"
#include "AliHLTTPCFastTransform.h"

class AliTPCParam;
class AliRecoParam;
class AliHLTTPCReverseTransformInfoV1;

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

  /** Initialisation  */
  Int_t  Init( double FieldBz, Long_t TimeStamp, bool isMC );
 
  /** Initialisation  */
  Int_t  Init( const AliHLTTPCFastTransformObject &obj );

  /** Initialised flag */
   Bool_t IsInitialised() const;

  /** Deinitialisation  */
   void DeInit();

  /** Setting the current time stamp  */
  Int_t SetCurrentTimeStamp( Long_t TimeStamp );
 
  /** Returns the current time stamp  */
  Long_t GetCurrentTimeStamp() const { return fFastTransform.GetCurrentTimeStamp(); }

  /** Transformation: calibration and alignment*/
  Int_t  Transform( int Slice, int Row, float Pad, float Time, float XYZ[] );

  /** Applying reverse alignment */
  int  ReverseAlignment( float XYZ[], int slice, int padrow);

  /** Last error message */
  const char* GetLastError() const { return fError.Data(); }

  /** Printout */
  void Print(const char* option=NULL) const;

  /** total size of the object*/
  Int_t GetSize() const ;

  const AliHLTTPCFastTransform& GetFastTransform(){ return fFastTransform; }
  
  AliHLTTPCFastTransform& GetFastTransformNonConst(){ return fFastTransform; }
  
  void SetInitSec(Int_t min, Int_t max) {fFastTransform.SetInitSec(min, max);}
  
  const AliHLTTPCReverseTransformInfoV1* GetReverseTransformInfo() {return fFastTransform.GetReverseTransformInfo();}

 private:

  /** copy constructor prohibited */
  AliHLTTPCClusterTransformation(const AliHLTTPCClusterTransformation&);
  /** assignment operator prohibited */
  AliHLTTPCClusterTransformation& operator=(const AliHLTTPCClusterTransformation&);

 /** Set error string */
  Int_t Error(Int_t code, const char *msg);

  static AliRecoParam    fOfflineRecoParam;  //! static container for TPC Reco Param

  TString fError; // Last error message

  AliHLTTPCFastTransform fFastTransform;// fast transformation object
  
  bool fIsMC; //Do we process MC?

  ClassDef(AliHLTTPCClusterTransformation, 1)
};

inline Int_t AliHLTTPCClusterTransformation::Error(Int_t code, const char *msg)
{
  // Set error 
  fError = msg;
  return code;
}

inline Int_t  AliHLTTPCClusterTransformation::Transform( int Slice, int Row, float Pad, float Time, float XYZ[] )
{
  // Convert row, pad, time to X Y Z   	   
  Int_t sector=-99, thisrow=-99;  
  AliHLTTPCGeometry::Slice2Sector( Slice, Row, sector, thisrow);
  int err = fFastTransform.Transform(sector, thisrow, Pad, Time, XYZ);
  if( err!=0 ) return Error(-1,Form( "AliHLTTPCClusterTransformation::Transform: Fast Transformation failed with error %d :%s",err,fFastTransform.GetLastError()) );
  return 0;
}

inline Int_t  AliHLTTPCClusterTransformation::ReverseAlignment( float XYZ[], int slice, int padrow)
{
  // reverse the alignment correction
  Int_t sector=-99, thisrow=-99;
  AliHLTTPCGeometry::Slice2Sector( slice, padrow, sector, thisrow);
  int err = fFastTransform.ReverseAlignment(sector, XYZ);
  if( err!=0 ) return Error(-1,Form( "AliHLTTPCClusterTransformation::ReverseAlignment: Fast Transformation failed with error %d :%s",err,fFastTransform.GetLastError()) );
  return 0;
}
#endif
