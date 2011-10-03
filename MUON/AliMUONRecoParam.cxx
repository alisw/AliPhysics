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

#include "AliCDBManager.h"
#include "AliCDBEntry.h"

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
  fRemoveConnectedTracksInSt12(kFALSE),
  fMaxTriggerTracks(0),
  fMaxTrackCandidates(0),
  fSelectTrackOnSlope(kFALSE),
  fMissingPadFractionLimit(-1),
  fFractionOfBuspatchOutsideOccupancyLimit(0),
  fAverageNoisePadCharge(0.22875),
  fClusterChargeCut(2.0),
  fEventSizeSoftLimit(35.0),
  fEventSizeHardLimit(45.0),
  fTokenLostLimit(0.0),
  fTryRecover(kFALSE)
{  
  /// Constructor
  
  SetNameTitle("Dummy","Dummy");
  for (Int_t iCh = 0; iCh < 10; iCh++) {
    fUseChamber[iCh] = kTRUE;
    fDefaultNonBendingReso[iCh] = 0.;
    fDefaultBendingReso[iCh] = 0.;
  }
  for (Int_t iSt = 0; iSt < 5; iSt++) fRequestStation[iSt] = kTRUE;
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
  /// INJECTIONGAIN : as GAIN, but with gain values taken as EMELEC factory values
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
AliMUONRecoParam *AliMUONRecoParam::GetCalibrationParam() 
{
  /// Return default (dummy) reconstruction parameters for calibration environment
  
  AliMUONRecoParam *param = new AliMUONRecoParam();
  param->SetCalibrationParam();
  
  return param;
}


//_____________________________________________________________________________
void AliMUONRecoParam::SetLowFluxParam() 
{
  /// Set reconstruction parameters for low flux environment
  
  SetNameTitle("Low Flux","Low Flux");
  SetEventSpecie(AliRecoParam::kLowMult);
  fMinBendingMomentum = 0.8;
  fMaxBendingMomentum = 1.e10;
  fMaxNonBendingSlope = 0.3;
  fMaxBendingSlope = 0.4;
  fSelectTrackOnSlope = kFALSE;
  fNonBendingVertexDispersion = 70.;
  fBendingVertexDispersion = 70.;
  fMaxNonBendingDistanceToTrack = 1.;
  fMaxBendingDistanceToTrack = 1.;
  fSigmaCutForTracking = 6.;
  fSigmaCutForImprovement = 5.;
  fSigmaCutForTrigger = 4.;
  fStripCutForTrigger = 1.;
  fMaxStripAreaForTrigger = 3.;
  fMaxNormChi2MatchTrigger = fSigmaCutForTrigger * fSigmaCutForTrigger;
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
  fMaxTriggerTracks = 100;
  fMaxTrackCandidates = 10000;  
}

//_____________________________________________________________________________
void AliMUONRecoParam::SetHighFluxParam() 
{
  /// Set reconstruction parameters for high flux environment
  
  SetNameTitle("High Flux","High Flux");
  SetEventSpecie(AliRecoParam::kHighMult);
  fMinBendingMomentum = 0.8;
  fMaxBendingMomentum = 1.e10;
  fMaxNonBendingSlope = 0.3;
  fMaxBendingSlope = 0.4;
  fSelectTrackOnSlope = kFALSE;
  fNonBendingVertexDispersion = 70.;
  fBendingVertexDispersion = 70.;
  fMaxNonBendingDistanceToTrack = 1.;
  fMaxBendingDistanceToTrack = 1.;
  fSigmaCutForTracking = 6.;
  fSigmaCutForImprovement = 5.;
  fSigmaCutForTrigger = 4.;
  fStripCutForTrigger = 1.;
  fMaxStripAreaForTrigger = 3.;
  fMaxNormChi2MatchTrigger = fSigmaCutForTrigger * fSigmaCutForTrigger;
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
  fMaxTriggerTracks = 100;
  fMaxTrackCandidates = 10000;
}

//_____________________________________________________________________________
void AliMUONRecoParam::SetCosmicParam() 
{
  /// Set reconstruction parameters for high flux environment
  
  SetNameTitle("Cosmic","Cosmic");
  SetEventSpecie(AliRecoParam::kCosmic);
  fMinBendingMomentum = 0.8;
  fMaxBendingMomentum = 1.e10;
  fMaxNonBendingSlope = 0.3;
  fMaxBendingSlope = 0.4;
  fSelectTrackOnSlope = kTRUE;
  fNonBendingVertexDispersion = 170.;
  fBendingVertexDispersion = 170.;
  fMaxNonBendingDistanceToTrack = 1.;
  fMaxBendingDistanceToTrack = 1.;
  fSigmaCutForTracking = 7.;
  fSigmaCutForImprovement = 6.;
  fSigmaCutForTrigger = 4.;
  fStripCutForTrigger = 1.5;
  fMaxStripAreaForTrigger = 3.;
  fMaxNormChi2MatchTrigger = fSigmaCutForTrigger * fSigmaCutForTrigger;
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
    fDefaultNonBendingReso[iCh] = 0.4;
    fDefaultBendingReso[iCh] = 0.4;
  }
  fRequestStation[0] = kTRUE;
  fRequestStation[1] = kTRUE;
  fRequestStation[2] = kTRUE;
  fRequestStation[3] = kTRUE;
  fRequestStation[4] = kTRUE;
  fBypassSt45 = 0;
  fPadGoodnessMask = 0x400BE80; // Ped Mean is Zero | Ped Mean Too Low | Ped Mean Too High | Ped Sigma Too Low | Ped Sigma Too High | Ped is missing | HV is missing | manu occupancy too high
  fMaxTriggerTracks = 100;
  fMaxTrackCandidates = 10000;
  
  SetPedMeanLimits(20, 700);
  SetManuOccupancyLimits(-1.,0.01); // reject manu above occ=1%

  SetBuspatchOccupancyLimits(-1,0.05);  
  SetFractionOfBuspatchOutsideOccupancyLimit(0.10); // 10 %
}


//_____________________________________________________________________________
void AliMUONRecoParam::SetCalibrationParam() 
{
  /// Set (dummy) reconstruction parameters for calibration environment
  
  SetNameTitle("Calibration","Calibration");
  SetEventSpecie(AliRecoParam::kCalib);

  fPedMeanLimits[0] = 5000;
  fPedMeanLimits[1] = 0;

  fPadGoodnessMask = 0x8C00; // Pedestal is missing | is too low | too high

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
  
  cout << "Event Specie=" << GetEventSpecie() << endl;
  
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
    
  cout<<"Selection of track candidates:"<<endl;
  if (fSelectTrackOnSlope) cout<<Form("\t- Non-bending slope < %5.2f",fMaxNonBendingSlope)<<endl;
  else cout<<"\t- Impact parameter < 3 * vertex dispersion in the non-bending direction"<<endl;
  cout<<Form("\t- if B!=0: Bending momentum > %5.2f",fMinBendingMomentum)<<endl;
  if (fSelectTrackOnSlope) cout<<Form("\t  if B==0: Bending slope < %5.2f",fMaxBendingSlope)<<endl;
  else cout<<"\t  if B==0: Impact parameter < 3 * vertex dispersion in the bending direction"<<endl;
  
  cout<<Form("Vertex dispersion (used to estimate initial bending momentum resolution) = (%5.2f,%5.2f)",fNonBendingVertexDispersion,fBendingVertexDispersion)<<endl;
  
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
  
  for ( int ichamber = 0; ichamber < 10; ++ichamber ) 
  {
    cout << Form("HV Ch %d must be >= %5.2f",ichamber,HVLimit(ichamber)) << endl;
  }

  cout << Form("%7.2f <= Pedestal mean <= %7.2f",PedMeanLowLimit(),PedMeanHighLimit()) << endl;
  cout << Form("%7.2f <= Pedestal sigma <= %7.2f",PedSigmaLowLimit(),PedSigmaHighLimit()) << endl;
  cout << Form("%e <= Gain linear term <= %e",GainA1LowLimit(),GainA1HighLimit()) << endl;
  cout << Form("%e <= Gain quadratic term <= %e",GainA2LowLimit(),GainA2HighLimit()) << endl;
  cout << Form("%5.0f <= Gain threshold term <= %5.0f",GainThresLowLimit(),GainThresHighLimit()) << endl;
    
  cout << Form("And we cut on charge >= %7.2f x ( pedestal sigma ) ",ChargeSigmaCut()) << endl;
  
  cout << "Occupancy limits are :" << endl;
  
  cout << Form("%e <= Manu occupancy < %7.2f",ManuOccupancyLowLimit(),ManuOccupancyHighLimit()) << endl;
  cout << Form("%e <= Buspatch occupancy < %7.2f",BuspatchOccupancyLowLimit(),BuspatchOccupancyHighLimit()) << endl;
  cout << Form("%e <= DE occupancy < %7.2f",DEOccupancyLowLimit(),DEOccupancyHighLimit()) << endl;
  
  cout << "'QAChecker' limits" << endl;  
  cout << Form("FractionOfBuspatchOutsideOccupancyLimit = %5.2f %%",FractionOfBuspatchOutsideOccupancyLimit()*100.0) << endl;
  cout << Form("Event size limit = %5.2f KB/event (soft) and %5.2f KB/event (hard)",fEventSizeSoftLimit,fEventSizeHardLimit) << endl;
  if ( fTokenLostLimit > 0 )
  {
    cout << Form("We tolerate up to %5.2f %% token lost errors per event",fTokenLostLimit) << endl;
  }
  else
  {
    cout << "We dot not tolerate any token lost error !" << endl;
  }
  
  cout << "chamber non bending resolution = |";
  for (Int_t iCh = 0; iCh < 10; iCh++) cout << Form(" %6.3f |",fDefaultNonBendingReso[iCh]);
  cout << endl;
  cout << "chamber bending resolution = |";
  for (Int_t iCh = 0; iCh < 10; iCh++) cout << Form(" %6.3f |",fDefaultBendingReso[iCh]);
  cout << endl;
  cout<<Form("maximum number of trigger tracks above which the tracking is cancelled = %d",fMaxTriggerTracks)<<endl;
  cout<<Form("maximum number of track candidates above which the tracking is abandonned = %d",fMaxTrackCandidates)<<endl;

  cout << Form("The average noise pad charge is assumed to be %7.2f fC",AverageNoisePadCharge()) << endl;
  cout << Form("and clusters below %5.2f times this noise charge (i.e. %7.2f fC) are discarded",
               ClusterChargeCut(),ClusterChargeCut()*AverageNoisePadCharge()) << endl;
  cout << Form("Note that LowestPadCharge is then %7.2f fC",LowestPadCharge()) << endl;
  
  if (TryRecover())
  {
    cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    cout << "!!! WILL TRY TO RECOVER CORRUPTED RAW DATA !!!" << endl;
    cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;    
  }
  cout<<"\t-----------------------------------------------------"<<endl<<endl;
  
}

//_____________________________________________________________________________
void
AliMUONRecoParam::SetHVLimit(Int_t chamberId, Double_t value)
{
  /// Set the HV limit for a given chamber (or all chambers 
  /// if chamberId==-1
  
  if ( chamberId == -1 ) 
  {
    for ( Int_t i = 0; i < 10; ++i ) 
    {
      fHVLimit[i] = value;
    }
  }
  else if ( chamberId >= 0 && chamberId < 10 ) 
  {
    fHVLimit[chamberId]=value;
  }
  else
  {
    AliError(Form("chamberId = %d is not a valid chamberId",chamberId));
  }
}

//_____________________________________________________________________________
Double_t AliMUONRecoParam::HVLimit(Int_t chamberId) const
{
  /// Get the HV limit for a given chamber
  if ( chamberId >= 0 && chamberId < 10 )
  {
    return fHVLimit[chamberId];
  }
  AliError(Form("chamberId = %d is not a valid chamberId",chamberId));

  return 0.0;
}

//_____________________________________________________________________________
void
AliMUONRecoParam::SetDefaultLimits()
{
	/// Set the default limits and pad goodness policy

  fHVSt12Limits[0]=1500; // kept for backward compatibility only
	fHVSt12Limits[1]=2000; // kept for backward compatibility only
	fHVSt345Limits[0]=1500; // kept for backward compatibility only
	fHVSt345Limits[1]=2000; // kept for backward compatibility only
  
  SetHVLimit(-1,1600); // this one is the real HV limit used now
  
	fPedMeanLimits[0] = 20;
	fPedMeanLimits[1] = 1024;
	
	fPedSigmaLimits[0] = 0.6;
	fPedSigmaLimits[1] = 100;

	fGainA1Limits[0] = 0.1;
	fGainA1Limits[1] = 10;

	fGainA2Limits[0] = -1E30;
	fGainA2Limits[1] = 1E30;
	
	fGainThresLimits[0] = 0;
	fGainThresLimits[1] = 4095;
	
	fPadGoodnessMask = 0x8080; // Ped is missing | HV is missing

  fManuOccupancyLimits[0] = -1.0; 
  fManuOccupancyLimits[1] = 1.0;

  fBuspatchOccupancyLimits[0] = 1E-6; 
  fBuspatchOccupancyLimits[1] = 1.0;

  fDEOccupancyLimits[0] = -1.0; 
  fDEOccupancyLimits[1] = 1.0;

  fMissingPadFractionLimit = -1; // DEPRECATED
  fFractionOfBuspatchOutsideOccupancyLimit = 0.05; // 5 % 

  ChargeSigmaCut(4.0); // pad with charge < 4.0 x sigma will be removed (where sigma is the actual noise of that very pad, i.e. not the average)
  
  AverageNoisePadCharge(0.22875); // 0.22875 coulombs ~ 1.5 ADC channels

  ClusterChargeCut(2.0); // will cut cluster below 2.0 x LowestPadCharge()
  
  SetEventSizeLimits(35.0,45.0);
  
  SetTokenLostLimit(0.0);
  
  fTryRecover = kFALSE;
}


//-----------------------------------------------------------------------
TObjArray* 
AliMUONRecoParam::Create(const char* settings)
{
  /// Create pre-defined recoparam array, according to settings.
  /// settings is case-insensitive.
  ///
  /// Currently defined are :
  ///
  /// "cosmics" :
  ///      Cosmic (default)
  ///      Calibration
  /// "ppideal"
  ///      LowFlux (default)
  ///      Calibration
  /// "ppreal"
  ///      LowFlux (modified to reconstruct real p-p data)
  ///      Calibration
  /// "pprealsim"
  ///      LowFlux (modified to reconstruct realistic p-p simulation)
  ///      Calibration
  
  AliMUONRecoParam* param(0x0);
  
  AliRecoParam::EventSpecie_t defaultParam = AliRecoParam::kLowMult;
  
  TString stype(settings);
  stype.ToLower();
  
  if ( stype == "cosmics" )
  {
    // set parameters for cosmic runs
    param = AliMUONRecoParam::GetCosmicParam();
    defaultParam = AliRecoParam::kCosmic;
  }
  else if ( stype == "ppideal" ) 
  {
    // set default lowFlux parameters
    param = AliMUONRecoParam::GetLowFluxParam();
  }
  else if ( stype == "ppreal" || stype == "pprealsim" || stype == "pprealnofield" ) 
  {      
    // common parameters for p-p data and realistic p-p simu
    param = AliMUONRecoParam::GetLowFluxParam();
    param->SaveFullClusterInESD(kTRUE, 100.);
    for (Int_t iCh=0; iCh<10; iCh++) 
    {
      param->SetDefaultNonBendingReso(iCh,0.4);
      param->SetDefaultBendingReso(iCh,0.4);
    }
    param->SetSigmaCutForTracking(7.);
    param->SetStripCutForTrigger(1.5);
    param->SetSigmaCutForTrigger(6.);
    param->ImproveTracks(kTRUE, 6.);
    param->SetPedMeanLimits(20, 700);
    param->SetManuOccupancyLimits(-1.,0.01);
    param->SetBuspatchOccupancyLimits(-1.,0.01);  
    param->SetFractionOfBuspatchOutsideOccupancyLimit(0.05); // 5 %
    param->SetEventSizeLimits(45., 65.);
    
    // specific parameters for p-p data or realistic p-p simu
    if ( stype == "ppreal" || stype == "pprealnofield" )
    {
      param->SetPadGoodnessMask(0x400BE9B);
    }
    else
    {
      param->SetPadGoodnessMask(0x8080);      
    }
    
    if ( stype == "pprealnofield" )
    {
      param->TryRecover(kTRUE);
    }
  }
  else if ( stype == "pbpbreal" || stype == "pbpbrealsim" ) 
  {      
    // common parameters for Pb-Pb data and realistic Pb-Pb simu
    param = AliMUONRecoParam::GetHighFluxParam();
    defaultParam = AliRecoParam::kHighMult;
    param->SaveFullClusterInESD(kTRUE, 100.);
    for (Int_t iCh=0; iCh<10; iCh++) 
    {
      param->SetDefaultNonBendingReso(iCh,0.2);
      param->SetDefaultBendingReso(iCh,0.2);
    }
    param->SetSigmaCutForTracking(5.);
    param->SetStripCutForTrigger(1.5);
    param->SetSigmaCutForTrigger(4.);
    param->ImproveTracks(kTRUE, 4.);
    param->SetPedMeanLimits(20, 700);
    param->SetManuOccupancyLimits(-1.,0.01);
    param->SetBuspatchOccupancyLimits(-1.,0.01);  
    param->SetFractionOfBuspatchOutsideOccupancyLimit(0.05); // 5 %
    param->SetEventSizeLimits(100., 150.);
    
    // specific parameters for Pb-Pb data or realistic Pb-Pb simu
    if ( stype == "pbpbreal" )
    {
      param->SetPadGoodnessMask(0x400BE9B);
    }
    else
    {
      param->SetPadGoodnessMask(0x8080);      
    }
  }
  else
  {
    AliErrorClass("Unknown settings !");
    return 0x0;
  }

  TObjArray* recoParams = new TObjArray;

  recoParams->AddLast(param);
  
  // set (dummy) parameters for calibration runs
  param = AliMUONRecoParam::GetCalibrationParam();
  recoParams->AddLast(param);
  
  // set parameters for Pb-Pb runs
  // param = AliMUONRecoParam::GetHighFluxParam();
  // recoParams.AddLast(param);
  
  // identify default parameters (exit if identification failed)
  Bool_t defaultIsSet = kFALSE;
  TIter next(recoParams);
  while ( (param = static_cast<AliMUONRecoParam*>(next())) ) 
  {
    if (param->GetEventSpecie() == defaultParam) 
    {
      param->SetAsDefault();
      defaultIsSet = kTRUE;
    }
    param->Print("FULL");
  }
  
  if (!defaultIsSet) 
  {
    AliErrorClass("The default reconstruction parameters are not set! Exiting...");
    return 0x0;
  }  
  
  return recoParams;
}

//______________________________________________________________________________
void 
AliMUONRecoParam::Show(Int_t runNumber, const char* ocdb)
{
  /// Show what we have in the designated OCDB for that run, as far as RecoParams are concerned
  
  AliCDBManager::Instance()->SetDefaultStorage(ocdb);
  AliCDBManager::Instance()->SetRun(runNumber);
  
  AliCDBEntry* entry = AliCDBManager::Instance()->Get("MUON/Calib/RecoParam");
  
  if (!entry) return;
  
  TObject* o = entry->GetObject();
  
  if ( o->IsA() == TObjArray::Class() ) 
  {
    TObjArray* array = static_cast<TObjArray*>(o);
    for ( Int_t i = 0; i <= array->GetLast(); ++i ) 
    {
      AliDetectorRecoParam* p = static_cast<AliDetectorRecoParam*>(array->At(i));
      cout << Form("array[%d]=%s %s %s",i,
                   p ? p->ClassName() : "",
                   p ? AliRecoParam::GetEventSpecieName(AliRecoParam::Convert(p->GetEventSpecie())) :"",
                   p ? ( p->IsDefault() ? "default" : "") : "" ) << endl;
    }
    cout << "=========== dumps below ====== " << endl;
    
    for ( Int_t i = 0; i <= array->GetLast(); ++i ) 
    {
      AliDetectorRecoParam* p = static_cast<AliDetectorRecoParam*>(array->At(i));
      if ( p ) p->Print("");
    }
  }
  else
  {
    o->Print();
  }
}
