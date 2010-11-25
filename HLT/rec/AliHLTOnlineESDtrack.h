//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTONLINEESDTRACK_H
#define ALIHLTONLINEESDTRACK_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTOnlineESDtrack.h
/// @author Matthias Richter
/// @date   2010-10-29
/// @brief  A streamlined container class for AliESDtrack.
/// @note   

#include "AliExternalTrackParam.h"
#include "AliESDtrack.h"

/**
 * @class AliHLTOnlineESDtrack
 * @brief Container for AliESDtrack information relevant for HLT.
 *
 * The class implements a reduced set of member variables of AliESDtrack
 * and corresponding transformation functions and operators.
 * It is used by the custom streamer of AliHLTESDEvent.
 */
class AliHLTOnlineESDtrack : public AliExternalTrackParam {
 public:
  /// standard constructor
  AliHLTOnlineESDtrack();
  /// copy constructor
  AliHLTOnlineESDtrack(const AliHLTOnlineESDtrack& t);
  /// destructor
  virtual ~AliHLTOnlineESDtrack();

  AliHLTOnlineESDtrack& operator=(const AliHLTOnlineESDtrack& t);
  AliHLTOnlineESDtrack& operator=(const AliESDtrack& t);

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

  static void CopyInternalParam(AliExternalTrackParam* &internalParam, const AliExternalTrackParam* pSrc);

  ///////////////////////////////////////////////////////////////////////////
  // access methods as in AliESDtrack
  //
  const AliExternalTrackParam * GetConstrainedParam() const {return fCp;}
  const AliExternalTrackParam * GetInnerParam() const { return fIp;}
  const AliExternalTrackParam * GetTPCInnerParam() const {return fTPCInner;}
  const AliExternalTrackParam * GetOuterParam() const { return fOp;}

  ULong_t   GetStatus() const {return fFlags;}
  Int_t     GetID() const { return fID;}
  Int_t     GetLabel() const {return fLabel;}
  Int_t     GetTPCLabel() const {return fTPCLabel;}
  Int_t     GetITSLabel() const {return fITSLabel;}
  Int_t     GetTRDLabel() const {return fTRDLabel;}
  Float_t   GetIntegratedLength() const {return fTrackLength;}
  void      GetImpactParametersTPC(Float_t p[2], Float_t cov[3]) const {
    p[0]=fdTPC; p[1]=fzTPC; cov[0]=fCddTPC; cov[1]=fCdzTPC; cov[2]=fCzzTPC;
  }
  Float_t   GetConstrainedChi2TPC() const {return fCchi2TPC;}
  void      GetImpactParameters(Float_t p[2], Float_t cov[3]) const {
    p[0]=fD; p[1]=fZ; cov[0]=fCdd; cov[1]=fCdz; cov[2]=fCzz;
  }
  Float_t   GetConstrainedChi2() const {return fCchi2;}

  Float_t   GetITSchi2() const {return fITSchi2;}
  Float_t   GetTPCchi2() const {return fTPCchi2;}
  Float_t   GetTPCchi2Iter1() const {return fTPCchi2Iter1;}
  UShort_t  GetTPCNcls() const { return fTPCncls;}
  UShort_t  GetTPCNclsF() const { return fTPCnclsF;}
  UShort_t  GetTPCNclsIter1() const { return fTPCnclsIter1;}
  UShort_t  GetTPCNclsFIter1() const { return fTPCnclsFIter1;}
  UShort_t  GetITSNcls() const { return fITSncls;}
  UChar_t   GetTRDncls() const {return fTRDncls;}
  UChar_t   GetTRDncls0() const {return fTRDncls0;}

private:
  // reduced set of parameters from AliESDtrack
  AliExternalTrackParam *fCp; // Track parameters constrained to the primary vertex
  AliExternalTrackParam *fIp; // Track parameters estimated at the inner wall of TPC
  AliExternalTrackParam *fTPCInner; // Track parameters estimated at the inner wall of TPC using the TPC stand-alone 
  AliExternalTrackParam *fOp; // Track parameters estimated at the point of maximal radial coordinate reached during the tracking

  ULong_t   fFlags;          // Reconstruction status flags 
  Int_t     fID;             // Unique ID of the track
  Int_t     fLabel;          // Track label
  Int_t     fITSLabel;       // label according ITS
  Int_t     fTPCLabel;       // label according TPC
  Int_t     fTRDLabel;       // label according TRD

  Float_t   fTrackLength;   // Track length

  Float_t   fdTPC;          // TPC-only impact parameter in XY plane
  Float_t   fzTPC;          // TPC-only impact parameter in Z
  Float_t   fCddTPC,fCdzTPC,fCzzTPC; // Covariance matrix of the TPC-only impact parameters 
  Float_t   fCchi2TPC;      // [0.,0.,8] TPC-only chi2 at the primary vertex

  Float_t   fD;             // Impact parameter in XY plane
  Float_t   fZ;             // Impact parameter in Z
  Float_t   fCdd,fCdz,fCzz; // Covariance matrix of the impact parameters 
  Float_t   fCchi2;          // [0.,0.,8] chi2 at the primary vertex

  Float_t   fITSchi2;        // [0.,0.,8] chi2 in the ITS
  Float_t   fTPCchi2;        // [0.,0.,8] chi2 in the TPC
  Float_t   fTPCchi2Iter1;  // [0.,0.,8] chi2 in the TPC

  UShort_t  fTPCncls;       // number of clusters assigned in the TPC
  UShort_t  fTPCnclsF;      // number of findable clusters in the TPC
  UShort_t  fTPCnclsIter1;  // number of clusters assigned in the TPC - iteration 1
  UShort_t  fTPCnclsFIter1; // number of findable clusters in the TPC - iteration 1

  Char_t    fITSncls;        // number of clusters assigned in the ITS
  UChar_t   fITSClusterMap;  // map of clusters, one bit per a layer
  UChar_t   fTRDncls;        // number of clusters assigned in the TRD
  UChar_t   fTRDncls0;       // number of clusters assigned in the TRD before first material cross

  ClassDef(AliHLTOnlineESDtrack, 1); // AliESDtrack instance optimized for HLT
};

#endif
