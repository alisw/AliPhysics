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

// $Id$

//-----------------------------------------------------------------------------
/// \class AliMUONRecoParam
///
/// Class with MUON reconstruction parameters
///
///  \author Philippe Pillot
//-----------------------------------------------------------------------------



#include "AliMUONRecoParam.h"
#include "AliMUONPadStatusMaker.h"

#include "AliRecoParam.h"
#include "AliLog.h"

#include <Riostream.h>

ClassImp(AliMUONRecoParam)


//_____________________________________________________________________________
AliMUONRecoParam::AliMUONRecoParam()
: AliDetectorRecoParam(),
  fClusteringMode("MLEM"),
  fTrackingMode("KALMAN"),
  fMinBendingMomentum(0.),
  fMaxBendingMomentum(0.),
  fMaxNonBendingSlope(0.),
  fMaxBendingSlope(0.),
  fNonBendingVertexDispersion(0.),
  fBendingVertexDispersion(0.),
  fMaxNonBendingDistanceToTrack(0.),
  fMaxBendingDistanceToTrack(0.),
  fSigmaCutForTracking(0.),
  fSigmaCutForImprovement(0.),
  fSigmaCutForTrigger(0.),
  fStripCutForTrigger(0.),
  fMaxStripAreaForTrigger(0.),
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
  fBypassSt45(0),
  fPadGoodnessMask(0),
  fChargeSigmaCut(4.0),
  fRemoveConnectedTracksInSt12(kFALSE)
{
  /// Constructor
  
  SetNameTitle("Dummy","Dummy");
  SetDefaultLimits();
}

//_____________________________________________________________________________
AliMUONRecoParam::~AliMUONRecoParam() 
{
  /// Destructor
}

//_____________________________________________________________________________
void
AliMUONRecoParam::BypassSt45(Bool_t st4, Bool_t st5)
{
	/// Set the bypass status
	
	if ( st4 && st5 ) fBypassSt45 = 45;
	else if ( st4 ) fBypassSt45 = 4;
	else if ( st5 ) fBypassSt45 = 5;
	else fBypassSt45 = 0;
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
  
  SetNameTitle("Low Flux","Low Flux");
  SetEventSpecie(AliRecoParam::kLowMult);
  fMinBendingMomentum = 1.;
  fMaxBendingMomentum = 3000.;
  fMaxNonBendingSlope = 0.3;
  fMaxBendingSlope = 0.4;
  fNonBendingVertexDispersion = 10.;
  fBendingVertexDispersion = 10.;
  fMaxNonBendingDistanceToTrack = 1.;
  fMaxBendingDistanceToTrack = 1.;
  fSigmaCutForTracking = 6.;
  fSigmaCutForImprovement = 5.;
  fSigmaCutForTrigger = 8.;
  fStripCutForTrigger = 1.;
  fMaxStripAreaForTrigger = 3.;
  fMaxNormChi2MatchTrigger = 16.;
  fCombinedClusterTrackReco = kFALSE;
  fTrackAllTracks = kTRUE;
  fRecoverTracks = kTRUE;
  fMakeTrackCandidatesFast = kFALSE;
  fMakeMoreTrackCandidates = kFALSE;
  fComplementTracks = kTRUE;
  fImproveTracks = kTRUE;
  fRemoveConnectedTracksInSt12 = kTRUE;
  fUseSmoother = kTRUE;
  for (Int_t iCh = 0; iCh < 10; iCh++) {
    fUseChamber[iCh] = kTRUE;
    fDefaultNonBendingReso[iCh] = 0.144;
    fDefaultBendingReso[iCh] = 0.01;
  }
  for (Int_t iSt = 0; iSt < 5; iSt++) fRequestStation[iSt] = kTRUE;
  fBypassSt45 = 0;
  
}

//_____________________________________________________________________________
void AliMUONRecoParam::SetHighFluxParam() 
{
  /// Set reconstruction parameters for high flux environment
  
  SetNameTitle("High Flux","High Flux");
  SetEventSpecie(AliRecoParam::kHighMult);
  fMinBendingMomentum = 1.;
  fMaxBendingMomentum = 3000.;
  fMaxNonBendingSlope = 0.3;
  fMaxBendingSlope = 0.4;
  fNonBendingVertexDispersion = 10.;
  fBendingVertexDispersion = 10.;
  fMaxNonBendingDistanceToTrack = 1.;
  fMaxBendingDistanceToTrack = 1.;
  fSigmaCutForTracking = 6.;
  fSigmaCutForImprovement = 5.;
  fSigmaCutForTrigger = 8.;
  fStripCutForTrigger = 1.;
  fMaxStripAreaForTrigger = 3.;
  fMaxNormChi2MatchTrigger = 16.;
  fCombinedClusterTrackReco = kFALSE;
  fTrackAllTracks = kTRUE;
  fRecoverTracks = kTRUE;
  fMakeTrackCandidatesFast = kFALSE;
  fMakeMoreTrackCandidates = kFALSE;
  fComplementTracks = kTRUE;
  fImproveTracks = kTRUE;
  fRemoveConnectedTracksInSt12 = kFALSE;
  fUseSmoother = kTRUE;
  for (Int_t iCh = 0; iCh < 10; iCh++) {
    fUseChamber[iCh] = kTRUE;
    fDefaultNonBendingReso[iCh] = 0.144;
    fDefaultBendingReso[iCh] = 0.01;
  }
  for (Int_t iSt = 0; iSt < 5; iSt++) fRequestStation[iSt] = kTRUE;
  fBypassSt45 = 0;
  
}

//_____________________________________________________________________________
void AliMUONRecoParam::SetCosmicParam() 
{
  /// Set reconstruction parameters for high flux environment
  
  SetNameTitle("Cosmic","Cosmic");
  SetEventSpecie(AliRecoParam::kCosmic);
  fMinBendingMomentum = 1.;
  fMaxBendingMomentum = 10000000.;
  fMaxNonBendingSlope = 0.3;
  fMaxBendingSlope = 0.4;
  fNonBendingVertexDispersion = 10.;
  fBendingVertexDispersion = 10.;
  fMaxNonBendingDistanceToTrack = 10.;
  fMaxBendingDistanceToTrack = 10.;
  fSigmaCutForTracking = 20.;
  fSigmaCutForImprovement = 20.;
  fSigmaCutForTrigger = 8.;
  fStripCutForTrigger = 1.5;
  fMaxStripAreaForTrigger = 3.;
  fMaxNormChi2MatchTrigger = 16.;
  fPercentOfFullClusterInESD = 100.;
  fCombinedClusterTrackReco = kFALSE;
  fTrackAllTracks = kTRUE;
  fRecoverTracks = kTRUE;
  fMakeTrackCandidatesFast = kFALSE;
  fMakeMoreTrackCandidates = kFALSE;
  fComplementTracks = kTRUE;
  fImproveTracks = kTRUE;
  fRemoveConnectedTracksInSt12 = kTRUE;
  fUseSmoother = kTRUE;
  fSaveFullClusterInESD = kTRUE;
  for (Int_t iCh = 0; iCh < 10; iCh++) {
    fUseChamber[iCh] = kTRUE;
    fDefaultNonBendingReso[iCh] = 0.144;
    fDefaultBendingReso[iCh] = 0.01;
  }
  for (Int_t iSt = 0; iSt < 5; iSt++) fRequestStation[iSt] = kTRUE;
  fBypassSt45 = 0;
  
}

//_____________________________________________________________________________
UInt_t
AliMUONRecoParam::RequestedStationMask() const
{
  /// Get the mask of the requested station, i.e. an integer where 
  /// bit n is set to one if the station n was requested
  
  UInt_t m(0);
  
  for ( Int_t i = 0; i < 5; ++i ) 
  {
    if ( RequestStation(i) ) m |= ( 1 << i );
  }
  return m;
}

//_____________________________________________________________________________
void AliMUONRecoParam::Print(Option_t *option) const
{
  /// print reconstruction parameters
  /// if option = FULL then print also unused parameters
  
  cout<<endl<<"\t------MUON Reconstruction parameters ("<<GetName()<<")------"<<endl;
  
  if (IsDefault()) cout<<"\t\t*** Parameters used by default ***"<<endl;
  
  cout<<Form("Calibration mode = %s",fCalibrationMode.Data())<<endl;
  cout<<Form("Clustering mode = %s",fClusteringMode.Data())<<endl;
  cout<<Form("Tracking mode = %s",fTrackingMode.Data())<<endl;

	TString bypass;
	
	if ( BypassSt45() )
	{
		bypass = "stations 4 and 5";
	}
	else if ( BypassSt4() ) 
	{
		bypass = "station 4";
	}
	else if ( BypassSt5() ) 
	{
		bypass = "station 5";
	}
	
  if (bypass.Length()) cout << "Will bypass " << bypass.Data() << " (replacing real clusters by generated ones from trigger tracks)" << endl;
  
  if (fCombinedClusterTrackReco) cout<<"Combined cluster/track reconstruction: ON"<<endl;
  else cout<<"Combined cluster/track reconstruction: OFF"<<endl;
  
  if (fSaveFullClusterInESD) cout<<Form("Save all cluster info in ESD for %5.2f %% of events",fPercentOfFullClusterInESD)<<endl;
  else cout<<"Save partial cluster info in ESD"<<endl;
    
  cout<<Form("Bending momentum range = [%5.2f,%5.2f]",fMinBendingMomentum,fMaxBendingMomentum)<<endl;
  
  cout<<Form("Maximum non bending slope = %5.2f",fMaxNonBendingSlope)<<endl;
  
  cout<<Form("Maximum bending slope (used only if B=0) = %5.2f",fMaxBendingSlope)<<endl;
  
  if (strstr(fTrackingMode,"ORIGINAL"))
    cout<<Form("Vertex dispertion = (%5.2f,%5.2f)",fNonBendingVertexDispersion,fBendingVertexDispersion)<<endl;
  else if (strstr(option,"FULL"))
    cout<<Form("Vertex dispertion (used for original tracking only) = (%5.2f,%5.2f)",fNonBendingVertexDispersion,fBendingVertexDispersion)<<endl;
  
  cout<<Form("Maximum distance to track = (%5.2f,%5.2f)",fMaxNonBendingDistanceToTrack,fMaxBendingDistanceToTrack)<<endl;
  
  cout<<Form("Sigma cut for tracking = %5.2f",fSigmaCutForTracking)<<endl;

  cout<<Form("Sigma cut for trigger hit pattern = %5.2f",fSigmaCutForTrigger)<<endl;

  cout<<Form("Cut in strips for trigger chamber efficiency = %5.2f",fStripCutForTrigger)<<endl;

  cout<<Form("Max search area in strips for trigger chamber efficiency = %5.2f",fMaxStripAreaForTrigger)<<endl;

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
  
  if (fRemoveConnectedTracksInSt12) cout<<"Remove tracks sharing one cluster or more in any station"<<endl;
  else cout<<"Remove tracks sharing one cluster or more in stations 3, 4 and 5"<<endl;
  
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
  
  cout << Form("Pad goodness policy mask is 0x%x",PadGoodnessMask()) << endl;
  cout << "Which means we reject pads having the condition = " <<
  AliMUONPadStatusMaker::AsCondition(PadGoodnessMask()).Data() << endl;
  
  cout << "The pad limits we are using are :" << endl;
  
  cout << Form("%5.0f <= HVSt12 <= %5.0f Volts",HVSt12LowLimit(),HVSt12HighLimit()) << endl;
  cout << Form("%5.0f <= HVSt345 <= %5.0f Volts",HVSt345LowLimit(),HVSt345HighLimit()) << endl;
  cout << Form("%7.2f <= Pedestal mean <= %7.2f",PedMeanLowLimit(),PedMeanHighLimit()) << endl;
  cout << Form("%7.2f <= Pedestal sigma <= %7.2f",PedSigmaLowLimit(),PedSigmaHighLimit()) << endl;
  cout << Form("%e <= Gain linear term <= %e",GainA1LowLimit(),GainA1HighLimit()) << endl;
  cout << Form("%e <= Gain quadratic term <= %e",GainA2LowLimit(),GainA2HighLimit()) << endl;
  cout << Form("%5.0f <= Gain threshold term <= %5.0f",GainThresLowLimit(),GainThresHighLimit()) << endl;
  
  cout << Form("And we cut on charge >= %7.2f x ( pedestal sigma ) ",ChargeSigmaCut()) << endl;
  
  cout << "chamber non bending resolution = |";
  for (Int_t iCh = 0; iCh < 10; iCh++) cout << Form(" %6.3f |",fDefaultNonBendingReso[iCh]);
  cout << endl;
  cout << "chamber bending resolution = |";
  for (Int_t iCh = 0; iCh < 10; iCh++) cout << Form(" %6.3f |",fDefaultBendingReso[iCh]);
  cout << endl;
  
  cout<<"\t-----------------------------------------------------"<<endl<<endl;
  
}

//_____________________________________________________________________________
void
AliMUONRecoParam::SetDefaultLimits()
{
	/// Set the default limits and pad goodness policy

	fHVSt12Limits[0]=1500;
	fHVSt12Limits[1]=2000;

	fHVSt345Limits[0]=1500;
	fHVSt345Limits[1]=2000;

	fPedMeanLimits[0] = 50;
	fPedMeanLimits[1] = 1024;
	
	fPedSigmaLimits[0] = 0.1;
	fPedSigmaLimits[1] = 100;

	fGainA1Limits[0] = 0.1;
	fGainA1Limits[1] = 10;

	fGainA2Limits[0] = -1E30;
	fGainA2Limits[1] = 1E30;
	
	fGainThresLimits[0] = 0;
	fGainThresLimits[1] = 4095;
	
	fPadGoodnessMask = 0x8080;
  
  fChargeSigmaCut = 4.0;
}

