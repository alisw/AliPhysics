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
: AliDetectorRecoParam(),
  fClusteringMode("MLEM"),
  fTrackingMode("KALMAN"),
  fMostProbBendingMomentum(0.),
  fMinBendingMomentum(0.),
  fMaxBendingMomentum(0.),
  fMaxNonBendingSlope(0.),
  fNonBendingVertexDispersion(0.),
  fBendingVertexDispersion(0.),
  fMaxNonBendingDistanceToTrack(0.),
  fMaxBendingDistanceToTrack(0.),
  fSigmaCutForTracking(0.),
  fSigmaCutForImprovement(0.),
  fSigmaCutForTrigger(0.),
  fMaxNormChi2MatchTrigger(0.),
  fPercentOfFullClusterInESD(10.),
  fCombinedClusterTrackReco(kFALSE),
  fTrackAllTracks(kFALSE),
  fRecoverTracks(kFALSE),
  fMakeTrackCandidatesFast(kFALSE),
  fMakeMoreTrackCandidates(kFALSE),
  fComplementTracks(kFALSE),
  fImproveTracks(kFALSE),
  fUseSmoother(kFALSE),
  fSaveFullClusterInESD(kTRUE),
  fCalibrationMode("NOGAIN"),
  fBypassSt45(kFALSE)
{
  /// Constructor
  SetNameTitle("MUON","MUON");
  
  // use the default parameters for low flux environment
  SetLowFluxParam();
}

//_____________________________________________________________________________
AliMUONRecoParam::~AliMUONRecoParam() 
{
  /// Destructor
}

//_____________________________________________________________________________
Option_t*
AliMUONRecoParam::GetCalibrationMode() const
{
  /// Return the calibration mode. Can be : 
  /// NOGAIN : only do pedestal subtraction
  /// GAIN : do pedestal subtraction, and apply gain correction, but with a
  ///        single capacitance value for all channels
  /// GAINCONSTANTCAPA : as GAIN, but with a channel-dependent capacitance value
  
  return fCalibrationMode.Data();
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
AliMUONRecoParam *AliMUONRecoParam::GetCosmicParam() 
{
  /// Return default reconstruction parameters for high flux environment
  
  AliMUONRecoParam *param = new AliMUONRecoParam();
  param->SetCosmicParam();
  
  return param;
}

//_____________________________________________________________________________
void AliMUONRecoParam::SetLowFluxParam() 
{
  /// Set reconstruction parameters for low flux environment
  
  fMostProbBendingMomentum = 2.;
  fMinBendingMomentum = 1.;
  fMaxBendingMomentum = 3000.;
  fMaxNonBendingSlope = 0.3;
  fNonBendingVertexDispersion = 10.;
  fBendingVertexDispersion = 10.;
  fMaxNonBendingDistanceToTrack = 1.;
  fMaxBendingDistanceToTrack = 1.;
  fSigmaCutForTracking = 6.;
  fSigmaCutForImprovement = 5.;
  fSigmaCutForTrigger = 8.;
  fMaxNormChi2MatchTrigger = 16.;
  fCombinedClusterTrackReco = kFALSE;
  fTrackAllTracks = kTRUE;
  fRecoverTracks = kTRUE;
  fMakeTrackCandidatesFast = kFALSE;
  fMakeMoreTrackCandidates = kFALSE;
  fComplementTracks = kTRUE;
  fImproveTracks = kTRUE;
  fUseSmoother = kTRUE;
  for (Int_t iCh = 0; iCh < 10; iCh++) fUseChamber[iCh] = kTRUE;
  for (Int_t iSt = 0; iSt < 5; iSt++) fRequestStation[iSt] = kTRUE;
  fBypassSt45 = kFALSE;
  
}

//_____________________________________________________________________________
void AliMUONRecoParam::SetHighFluxParam() 
{
  /// Set reconstruction parameters for high flux environment
  
  fMostProbBendingMomentum = 2.;
  fMinBendingMomentum = 1.;
  fMaxBendingMomentum = 3000.;
  fMaxNonBendingSlope = 0.3;
  fNonBendingVertexDispersion = 10.;
  fBendingVertexDispersion = 10.;
  fMaxNonBendingDistanceToTrack = 1.;
  fMaxBendingDistanceToTrack = 1.;
  fSigmaCutForTracking = 6.;
  fSigmaCutForImprovement = 5.;
  fSigmaCutForTrigger = 8.;
  fMaxNormChi2MatchTrigger = 16.;
  fCombinedClusterTrackReco = kFALSE;
  fTrackAllTracks = kTRUE;
  fRecoverTracks = kTRUE;
  fMakeTrackCandidatesFast = kFALSE;
  fMakeMoreTrackCandidates = kFALSE;
  fComplementTracks = kTRUE;
  fImproveTracks = kTRUE;
  fUseSmoother = kTRUE;
  for (Int_t iCh = 0; iCh < 10; iCh++) fUseChamber[iCh] = kTRUE;
  for (Int_t iSt = 0; iSt < 5; iSt++) fRequestStation[iSt] = kTRUE;
  fBypassSt45 = kFALSE;
  
}

//_____________________________________________________________________________
void AliMUONRecoParam::SetCosmicParam() 
{
  /// Set reconstruction parameters for high flux environment
  
  fMostProbBendingMomentum = 2.;
  fMinBendingMomentum = 1.;
  fMaxBendingMomentum = 10000000.;
  fMaxNonBendingSlope = 0.3;
  fNonBendingVertexDispersion = 10.;
  fBendingVertexDispersion = 10.;
  fMaxNonBendingDistanceToTrack = 10.;
  fMaxBendingDistanceToTrack = 10.;
  fSigmaCutForTracking = 20.;
  fSigmaCutForImprovement = 20.;
  fSigmaCutForTrigger = 8.;
  fMaxNormChi2MatchTrigger = 16.;
  fPercentOfFullClusterInESD = 100.;
  fCombinedClusterTrackReco = kFALSE;
  fTrackAllTracks = kTRUE;
  fRecoverTracks = kTRUE;
  fMakeTrackCandidatesFast = kFALSE;
  fMakeMoreTrackCandidates = kFALSE;
  fComplementTracks = kTRUE;
  fImproveTracks = kTRUE;
  fUseSmoother = kTRUE;
  fSaveFullClusterInESD = kTRUE;
  for (Int_t iCh = 0; iCh < 10; iCh++) fUseChamber[iCh] = kTRUE;
  for (Int_t iSt = 0; iSt < 5; iSt++) fRequestStation[iSt] = kTRUE;
  fBypassSt45 = kFALSE;
  
}

//_____________________________________________________________________________
void AliMUONRecoParam::Print(Option_t *option) const
{
  /// print reconstruction parameters
  /// if option = FULL then print also unused parameters
  
  cout<<endl<<"\t------Reconstruction parameters------"<<endl;
  
  cout<<Form("Calibration mode = %s",fCalibrationMode.Data())<<endl;
  cout<<Form("Clustering mode = %s",fClusteringMode.Data())<<endl;
  cout<<Form("Tracking mode = %s",fTrackingMode.Data())<<endl;

  if (BypassSt45()) cout << "Will bypass St45 (replacing their clusters by generated ones from trigger tracks)" << endl;
  
  if (fCombinedClusterTrackReco) cout<<"Combined cluster/track reconstruction: ON"<<endl;
  else cout<<"Combined cluster/track reconstruction: OFF"<<endl;
  
  if (fSaveFullClusterInESD) cout<<Form("Save all cluster info in ESD for %5.2f %% of events",fPercentOfFullClusterInESD)<<endl;
  else cout<<"Save partial cluster info in ESD"<<endl;
  
  cout<<Form("Most probable bending momentum (used only if B=0) = %5.2f",fMostProbBendingMomentum)<<endl;
  
  cout<<Form("Bending momentum range = [%5.2f,%5.2f]",fMinBendingMomentum,fMaxBendingMomentum)<<endl;
  
  cout<<Form("Maximum non bending slope = %5.2f",fMaxNonBendingSlope)<<endl;
  
  if (strstr(fTrackingMode,"ORIGINAL"))
    cout<<Form("Vertex dispertion = (%5.2f,%5.2f)",fNonBendingVertexDispersion,fBendingVertexDispersion)<<endl;
  else if (strstr(option,"FULL"))
    cout<<Form("Vertex dispertion (used for original tracking only) = (%5.2f,%5.2f)",fNonBendingVertexDispersion,fBendingVertexDispersion)<<endl;
  
  cout<<Form("Maximum distance to track = (%5.2f,%5.2f)",fMaxNonBendingDistanceToTrack,fMaxBendingDistanceToTrack)<<endl;
  
  cout<<Form("Sigma cut for tracking = %5.2f",fSigmaCutForTracking)<<endl;

  cout<<Form("Sigma cut for trigger hit pattern = %5.2f",fSigmaCutForTrigger)<<endl;
  
  if (fTrackAllTracks) cout<<"Track all the possible candidates"<<endl;
  else cout<<"Track only the best candidates"<<endl;
  
  if (strstr(option,"FULL")) {
    cout<<"Make track candidates assuming linear propagation between stations 4 and 5: ";
    if (fMakeTrackCandidatesFast) cout<<"ON"<<endl;
    else cout<<"OFF"<<endl;
  } else if (fMakeTrackCandidatesFast)
    cout<<"Make track candidates assuming linear propagation between stations 4 and 5"<<endl;
  
  if (strstr(option,"FULL")) {
    cout<<"Make track candidates starting from 1 cluster in each of the stations 4 and 5: ";
    if (fMakeMoreTrackCandidates) cout<<"ON"<<endl;
    else cout<<"OFF"<<endl;
  } else if (fMakeMoreTrackCandidates)
    cout<<"Make track candidates starting from 1 cluster in each of the stations 4 and 5"<<endl;
  
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
  
  Bool_t discardedCh = kFALSE;
  Int_t ch = 0;
  do {
    if (!UseChamber(ch)) {
      if (!discardedCh) {
	cout<<"Discarded chambers(1..): "<<ch+1;
	discardedCh = kTRUE;
      }
      else cout<<" "<<ch+1;
    }
  } while (++ch < 10);
  if (discardedCh) cout<<endl;
  
  Bool_t discardedSt = kFALSE;
  Int_t st = 0;
  do {
    if (!RequestStation(st)) {
      if (!discardedSt) {
	cout<<"Not requested stations(1..): "<<st+1;
	discardedSt = kTRUE;
      }
      else cout<<" "<<st+1;
    }
  } while (++st < 5);
  if (discardedSt) cout<<endl;
  
  cout<<"\t-------------------------------------"<<endl<<endl;
  
}

