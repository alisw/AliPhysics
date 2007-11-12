/**************************************************************************
* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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


//-----------------------------------------------------------------------------
/// \class AliMUONRecoParam
///
/// Class with MUON reconstruction parameters
///
///  \author Philippe Pillot
//-----------------------------------------------------------------------------



#include "AliMUONRecoParam.h"

#include "AliLog.h"

#include <Riostream.h>

ClassImp(AliMUONRecoParam)


//_____________________________________________________________________________
AliMUONRecoParam::AliMUONRecoParam()
: TObject(),
  fClusteringMode("MLEM"),
  fTrackingMode("KALMAN"),
  fMinBendingMomentum(0.),
  fMaxBendingMomentum(0.),
  fNonBendingVertexDispersion(0.),
  fBendingVertexDispersion(0.),
  fMaxNonBendingDistanceToTrack(0.),
  fMaxBendingDistanceToTrack(0.),
  fSigmaCutForTracking(0.),
  fSigmaCutForImprovement(0.),
  fMaxNormChi2MatchTrigger(0.),
  fTrackAllTracks(kFALSE),
  fRecoverTracks(kFALSE),
  fMakeTrackCandidatesFast(kFALSE),
  fComplementTracks(kFALSE),
  fImproveTracks(kFALSE),
  fUseSmoother(kFALSE)
{
  /// Constructor
  
  // use the default parameters for low flux environment
  SetLowFluxParam();
}

//_____________________________________________________________________________
AliMUONRecoParam::~AliMUONRecoParam() 
{
  /// Destructor
}

//_____________________________________________________________________________
AliMUONRecoParam *AliMUONRecoParam::GetLowFluxParam() 
{
  /// Return default reconstruction parameters for low flux environment
  
  AliMUONRecoParam *param = new AliMUONRecoParam();
  param->SetLowFluxParam();
  
  return param;
}

//_____________________________________________________________________________
AliMUONRecoParam *AliMUONRecoParam::GetHighFluxParam() 
{
  /// Return default reconstruction parameters for high flux environment
  
  AliMUONRecoParam *param = new AliMUONRecoParam();
  param->SetHighFluxParam();
 
  return param;
}

//_____________________________________________________________________________
void AliMUONRecoParam::SetLowFluxParam() 
{
  /// Set reconstruction parameters for low flux environment
  
  fMinBendingMomentum = 0.5;
  fMaxBendingMomentum = 3000.;
  fNonBendingVertexDispersion = 10.;
  fBendingVertexDispersion = 10.;
  fMaxNonBendingDistanceToTrack = 2.;
  fMaxBendingDistanceToTrack = 2.;
  fSigmaCutForTracking = 6.;
  fSigmaCutForImprovement = 5.;
  fMaxNormChi2MatchTrigger = 16.;
  fTrackAllTracks = kTRUE;
  fRecoverTracks = kTRUE;
  fMakeTrackCandidatesFast = kFALSE;
  fComplementTracks = kTRUE;
  fImproveTracks = kTRUE;
  fUseSmoother = kTRUE;
  
}

//_____________________________________________________________________________
void AliMUONRecoParam::SetHighFluxParam() 
{
  /// Set reconstruction parameters for high flux environment
  
  fMinBendingMomentum = 0.5;
  fMaxBendingMomentum = 3000.;
  fNonBendingVertexDispersion = 10.;
  fBendingVertexDispersion = 10.;
  fMaxNonBendingDistanceToTrack = 2.;
  fMaxBendingDistanceToTrack = 2.;
  fSigmaCutForTracking = 6.;
  fSigmaCutForImprovement = 5.;
  fMaxNormChi2MatchTrigger = 16.;
  fTrackAllTracks = kTRUE;
  fRecoverTracks = kTRUE;
  fMakeTrackCandidatesFast = kFALSE;
  fComplementTracks = kTRUE;
  fImproveTracks = kTRUE;
  fUseSmoother = kTRUE;
  
}

//_____________________________________________________________________________
void AliMUONRecoParam::Print(Option_t *option) const
{
  /// print reconstruction parameters
  /// if option = FULL then print also unused parameters
  
  cout<<endl<<"\t------Reconstruction parameters------"<<endl;
  
  cout<<Form("Clustering mode = %s",fClusteringMode.Data())<<endl;
  
  cout<<Form("Tracking mode = %s",fTrackingMode.Data())<<endl;
  
  cout<<Form("Bending momentum range = [%5.2f,%5.2f]",fMinBendingMomentum,fMaxBendingMomentum)<<endl;
  
  if (strstr(fTrackingMode,"ORIGINAL"))
    cout<<Form("Vertex dispertion = (%5.2f,%5.2f)",fNonBendingVertexDispersion,fBendingVertexDispersion)<<endl;
  else if (strstr(option,"FULL"))
    cout<<Form("Vertex dispertion (used for original tracking only) = (%5.2f,%5.2f)",fNonBendingVertexDispersion,fBendingVertexDispersion)<<endl;
  
  cout<<Form("Maximum distance to track = (%5.2f,%5.2f)",fMaxNonBendingDistanceToTrack,fMaxBendingDistanceToTrack)<<endl;
  
  cout<<Form("Sigma cut for tracking = %5.2f",fSigmaCutForTracking)<<endl;
  
  if (fTrackAllTracks) cout<<"Track all the possible candidates"<<endl;
  else cout<<"Track only the best candidates"<<endl;
  
  if (strstr(option,"FULL")) {
    cout<<"Make track candidates assuming linear propagation between stations 4 and 5: ";
    if (fMakeTrackCandidatesFast) cout<<"ON"<<endl;
    else cout<<"OFF"<<endl;
  } else if (fMakeTrackCandidatesFast)
    cout<<"Make track candidates assuming linear propagation between stations 4 and 5"<<endl;
  
  if (strstr(option,"FULL")) {
    cout<<"Try to recover tracks getting lost during tracking: ";
    if (fRecoverTracks) cout<<"ON"<<endl;
    else cout<<"OFF"<<endl;
  } else if (fRecoverTracks)
    cout<<"Try to recover tracks getting lost during tracking"<<endl;
  
  if (strstr(option,"FULL")) {
    cout<<"Try to complete the reconstructed tracks by adding missing clusters: ";
    if (fComplementTracks) cout<<"ON"<<endl;
    else cout<<"OFF"<<endl;
  } else if (fComplementTracks)
    cout<<"Try to complete the reconstructed tracks by adding missing clusters"<<endl;
  
  if (strstr(option,"FULL")) {
    cout<<"Try to improve the reconstructed tracks by removing bad clusters: ";
    if (fImproveTracks) cout<<Form("ON (sigma cut = %5.2f)",fSigmaCutForImprovement)<<endl;
    else cout<<"OFF"<<endl;
  } else if (fImproveTracks)
    cout<<Form("Try to improve the reconstructed tracks by removing bad clusters (sigma cut = %5.2f)",fSigmaCutForImprovement)<<endl;
  
  if (strstr(option,"FULL")) {
    cout<<"Use smoother to compute final track parameters, etc, at each cluster (used for Kalman tracking only): ";
    if (fUseSmoother) cout<<"ON"<<endl;
    else cout<<"OFF"<<endl;
  } else if (fUseSmoother)
    cout<<"Use smoother to compute final track parameters, etc, at each cluster"<<endl;
  
  cout<<Form("Maximum normalized chi2 of tracking/trigger track matching = %5.2f",fMaxNormChi2MatchTrigger)<<endl;
  
  cout<<"\t-------------------------------------"<<endl<<endl;
  
}

