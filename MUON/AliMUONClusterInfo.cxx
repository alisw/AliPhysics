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

// $Id$

//-----------------------------------------------------------------------------
/// \class AliMUONClusterInfo
///
/// Class to summarize ESD data at cluster
///
/// \author Philippe Pillot, Subatech
//-----------------------------------------------------------------------------

#include "AliMUONClusterInfo.h"

#include "AliLog.h"

#include <Riostream.h>

/// \cond CLASSIMP
ClassImp(AliMUONClusterInfo)
/// \endcond

//_____________________________________________________________________________
AliMUONClusterInfo::AliMUONClusterInfo()
: TObject(),
  fEventId(0),
  fZ(0.),
  fClusterId(0),
  fClusterX(0.),
  fClusterY(0.),
  fClusterXErr(0.),
  fClusterYErr(0.),
  fClusterChi2(0.),
  fClusterCharge(0.),
  fTrackId(0),
  fTrackX(0.),
  fTrackY(0.),
  fTrackThetaX(0.),
  fTrackThetaY(0.),
  fTrackP(0.),
  fTrackXErr(0.),
  fTrackYErr(0.),
  fTrackChi2(0.),
  fTrackCharge(0),
  fNPads(0),
  fPads(new TClonesArray("AliMUONPadInfo",10))
{
  /// default constructor
}

//_____________________________________________________________________________
AliMUONClusterInfo::AliMUONClusterInfo (const AliMUONClusterInfo& clusterInfo)
: TObject(clusterInfo),
  fEventId(clusterInfo.fEventId),
  fZ(clusterInfo.fZ),
  fClusterId(clusterInfo.fClusterId),
  fClusterX(clusterInfo.fClusterX),
  fClusterY(clusterInfo.fClusterY),
  fClusterXErr(clusterInfo.fClusterXErr),
  fClusterYErr(clusterInfo.fClusterYErr),
  fClusterChi2(clusterInfo.fClusterChi2),
  fClusterCharge(clusterInfo.fClusterCharge),
  fTrackId(clusterInfo.fTrackId),
  fTrackX(clusterInfo.fTrackX),
  fTrackY(clusterInfo.fTrackY),
  fTrackThetaX(clusterInfo.fTrackThetaX),
  fTrackThetaY(clusterInfo.fTrackThetaY),
  fTrackP(clusterInfo.fTrackP),
  fTrackXErr(clusterInfo.fTrackXErr),
  fTrackYErr(clusterInfo.fTrackYErr),
  fTrackChi2(clusterInfo.fTrackChi2),
  fTrackCharge(clusterInfo.fTrackCharge),
  fNPads(clusterInfo.fNPads),
  fPads(new TClonesArray("AliMUONPadInfo",clusterInfo.fNPads))
{
  /// Copy constructor
  AliMUONPadInfo *pad = (AliMUONPadInfo*) clusterInfo.fPads->First();
  while (pad) {
    new ((*fPads)[fPads->GetEntriesFast()]) AliMUONPadInfo(*pad);
    pad = (AliMUONPadInfo*) clusterInfo.fPads->After(pad);
  }
}

//_____________________________________________________________________________
AliMUONClusterInfo& AliMUONClusterInfo::operator=(const AliMUONClusterInfo& clusterInfo)
{
  /// Equal operator
  if (this == &clusterInfo) return *this;
  
  TObject::operator=(clusterInfo); // don't forget to invoke the base class' assignment operator
  
  fEventId = clusterInfo.fEventId;
  fZ = clusterInfo.fZ;
  fClusterId = clusterInfo.fClusterId;
  fClusterX = clusterInfo.fClusterX;
  fClusterY = clusterInfo.fClusterY;
  fClusterXErr = clusterInfo.fClusterXErr;
  fClusterYErr = clusterInfo.fClusterYErr;
  fClusterChi2 = clusterInfo.fClusterChi2;
  fClusterCharge = clusterInfo.fClusterCharge;
  fTrackId = clusterInfo.fTrackId;
  fTrackX = clusterInfo.fTrackX;
  fTrackY = clusterInfo.fTrackY;
  fTrackThetaX = clusterInfo.fTrackThetaX;
  fTrackThetaY = clusterInfo.fTrackThetaY;
  fTrackP = clusterInfo.fTrackP;
  fTrackXErr = clusterInfo.fTrackXErr;
  fTrackYErr = clusterInfo.fTrackYErr;
  fTrackChi2 = clusterInfo.fTrackChi2;
  fTrackCharge = clusterInfo.fTrackCharge;
  fNPads = clusterInfo.fNPads;
  
  fPads->Clear("C");
  AliMUONPadInfo *pad = (AliMUONPadInfo*) clusterInfo.fPads->First();
  while (pad) {
    new ((*fPads)[fPads->GetEntriesFast()]) AliMUONPadInfo(*pad);
    pad = (AliMUONPadInfo*) clusterInfo.fPads->After(pad);
  }
  
  return *this;
}

//__________________________________________________________________________
AliMUONClusterInfo::~AliMUONClusterInfo()
{
  /// Destructor
  delete fPads;
}

//__________________________________________________________________________
void AliMUONClusterInfo::Clear(Option_t* opt)
{
  /// Clear arrays
  fPads->Clear(opt);
  fNPads = 0;
}

//_____________________________________________________________________________
void AliMUONClusterInfo::Print(Option_t* option) const
{
  /// print cluster info content
  /// print also pad info if option=FULL
  
  // global info
  cout<<Form("eventID=%d", GetEventId())<<endl;
  
  // cluster info
  cout<<Form("- clusterID=%u (ch=%d, det=%d, index=%d)",
	     GetClusterId(), GetChamberId(), GetDetElemId(), GetClusterIndex())<<endl;
  
  cout<<Form("    position=(%5.2f, %5.2f, %5.2f), sigma=(%8.5f, %8.5f, 0.0)",
	     GetClusterX(), GetClusterY(), GetZ(), GetClusterXErr(), GetClusterYErr())<<endl;
  
  cout<<Form("    charge=%5.2f, chi2=%5.2f", GetClusterCharge(), GetClusterChi2())<<endl;
  
  // track info
  cout<<Form("- trackID=%u", GetTrackId())<<endl;
  
  cout<<Form("    position=(%5.2f, %5.2f, %5.2f), angles=(%5.2f, %5.2f), momentum=%5.2f",
	     GetTrackX(), GetTrackY(), GetZ(), GetTrackThetaX(), GetTrackThetaY(), GetTrackP())<<endl;
  
  cout<<Form("    sigma_XY=(%8.5f, %8.5f), charge=%d, chi2=%5.2f",
	     GetTrackXErr(), GetTrackYErr(), GetTrackCharge(), GetTrackChi2())<<endl;
  
  // pad info
  if (strstr(option,"FULL")) {
    AliMUONPadInfo *pad = (AliMUONPadInfo*) fPads->First();
    while (pad) {
      pad->Print("FULL");
      pad = (AliMUONPadInfo*) fPads->After(pad);
    }
  }
  
}

