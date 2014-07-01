/* $Id$ */

/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/**
 * >> interface class to a track/particle (normal/flat) <<
 *
 * To be used in the online and offline calibration schema.
 *
 * Primary Authors : Mikolaj Krzewicki (mkrzewic@cern.ch)
 *
 **************************************************************************/

#include "AliVVtrack.h"
  AliFlatExternalTrackParam* AliVVtrack::GetTrackParamRefitted() { return NULL; } 
  AliFlatExternalTrackParam* AliVVtrack::GetTrackParamIp() { return NULL; } 
  AliFlatExternalTrackParam* AliVVtrack::GetTrackParamTPCInner() { return NULL; } 
  AliFlatExternalTrackParam* AliVVtrack::GetTrackParamOp() { return NULL; } 
  AliFlatExternalTrackParam* AliVVtrack::GetTrackParamCp() { return NULL; } 
  AliFlatExternalTrackParam* AliVVtrack::GetTrackParamITSOut() { return NULL; } 

  // --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  
  Int_t AliVVtrack::GetNumberOfTPCClusters() { return 0; } 
  AliFlatTPCCluster* AliVVtrack::GetTPCClusters() { return NULL; } 
  AliFlatTPCCluster* AliVVtrack::GetTPCCluster(Int_t /*ind*/) { return NULL; } 
  Int_t AliVVtrack::GetNumberOfITSClusters() { return 0; }
  AliVVtrack *AliVVtrack::GetNextTrack() {return NULL; }

  Bool_t AliVVtrack::GetXYZ(Double_t* /*p*/) const {return kFALSE;}
  void  AliVVtrack::GetXYZ(Double_t& /*x*/, Double_t& /*y*/, Double_t& /*z*/) const {}
  Double_t AliVVtrack::GetTgl()  const {return 0.;}
  UShort_t AliVVtrack::GetTPCNclsF() const { return 0;}

  Double_t AliVVtrack::GetTOFsignalDz() const {return 0.;}
  void AliVVtrack::GetImpactParameters(Float_t& /*xy*/,Float_t& /*z*/) const {}
  //TODO:
  void AliVVtrack::GetDZ(Double_t /*x*/,Double_t /*y*/,Double_t /*z*/,Double_t /*b*/, Float_t dz[2]) const {if (dz[0]==0) return;}

  Float_t AliVVtrack::GetTPCClusterInfo(Int_t nNeighbours, Int_t type, Int_t row0, Int_t row1, Int_t bitType ) const {return 0.*nNeighbours*type*row0*row1*bitType;}
  UShort_t AliVVtrack::GetTPCncls(Int_t row0,Int_t row1) const {return 0*row0*row1;}
  Bool_t AliVVtrack::IsOn(Int_t /*mask*/) const {return kFALSE;}
  void AliVVtrack::GetDirection(Double_t d[3]) const {if (d[0]==0) return;}
  const Double_t *AliVVtrack::GetParameter() const {return 0;}
  void AliVVtrack::GetImpactParametersTPC(Float_t& /*xy*/,Float_t& /*z*/) const {}
  Int_t AliVVtrack::GetNumberOfClusters() const {return 0;} 
  const AliVVtrack* AliVVtrack::GetTPCInnerParam() const {return NULL;}
  Double_t AliVVtrack::Pt() const {return 0.;}
  Double_t AliVVtrack::GetP() const {return 0.;}
  Double_t AliVVtrack::GetTPCmomentum() const {return 0.;}
  ULong_t AliVVtrack::GetStatus() const {return 0;}
  const AliVVtrack * AliVVtrack::GetOuterParam() const { return NULL;}
  const AliVVtrack * AliVVtrack::GetInnerParam() const { return NULL;}
  Int_t AliVVtrack::GetKinkIndex(Int_t /*i*/) const { return 0;}
  Double_t AliVVtrack::Eta() const {return 0.;}
  Double_t AliVVtrack::GetY() const {return 0.;}
  
  Double_t AliVVtrack::GetX() const {return 0.;}
  Double_t AliVVtrack::GetZ() const {return 0.;}
  Int_t AliVVtrack::GetNcls(Int_t /*idet*/) const {return 0;}
  void AliVVtrack::GetIntegratedTimes(Double_t* /*times*/, Int_t nspec) const {if (nspec<0) return;}
  Double_t AliVVtrack::GetSigned1Pt()  const {return 0.;}  
  Double_t AliVVtrack::GetLinearD(Double_t /*xv*/, Double_t /*yv*/) const {return 0.;}
  const AliVVtrack *AliVVtrack::GetConstrainedParam() const {return NULL;}
  Double_t AliVVtrack::GetAlpha() const {return 0.;}
  Char_t AliVVtrack::GetITSclusters(Int_t* /*idx*/) const {return 0;}
  Double_t AliVVtrack::GetSign() const {return 0.;}
  UShort_t AliVVtrack::GetTPCNcls() const { return 0;}
  Float_t AliVVtrack::GetTPCCrossedRows() const {return 0.;}
  Double_t AliVVtrack::GetTPCsignal() const {return 0.;}
  Double_t AliVVtrack::GetTOFsignal() const {return 0.;}
  UChar_t AliVVtrack::GetTRDclusters(Int_t* /*idx*/) const {return 0;}
  
  //AliTPCtrack
  Int_t AliVVtrack::GetNFoundable() const {return 0;} 
  Double_t AliVVtrack::GetdEdx()  const {return 0.;}
