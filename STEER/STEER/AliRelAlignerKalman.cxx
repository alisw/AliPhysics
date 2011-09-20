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

///////////////////////////////////////////////////////////////////////////////
//
//    Kalman filter based aligner:
//    Finds alignement constants for  two tracking volumes (by default ITS
//    and TPC)
//    Determines the inverse transformation of the second volume (TPC)
//    with respect to the first (ITS) (how to realign TPC to ITS)
//    by measuring the residual between the 2 tracks.
//    Additionally calculates some callibration parameters for TPC
//    Fit parameters are:
//    - 3 shifts, x,y,z
//    - 3 Cardan angles, psi, theta, phi (see definition in alignment docs),
//    - TPC drift velocity correction,
//    - TPC time offset correction.
//
//    Basic usage:
//    When aligning two volumes, at any given time a single instance of
//    the class should be active. The fit of the parameters is updated
//    by adding new data using one of the Add.... methods:
//
//    In collision events add an ESD event to update the fit (adds all tracks):
//
//        Bool_t AddESDevent( AliESDevent* pTrack );
//    
//    or add each individual track
//
//        AddESDtrack( AliESDtrack* pTrack );
//
//    For cosmic data, the assumption is that the tracking is done twice:
//    once global and once only ITS and the tracklets are saved inside
//    one AliESDEvent. The method
//
//        Bool_t AddCosmicEvent( AliESDEvent* pEvent );
//
//    then searches the event for matching tracklets and upon succes it updates.
//    One cosmic ideally triggers two updates: for the upper and lower half of
//    the cosmic (upper ITS tracklet+upper TPC tracklet, idem dito for lower)
//
//    by default give misalignment parameters for TPC as they appear to be.
//    TPC calibration parameters are always given as correction to values used in reco.
//
//    _________________________________________________________________________
//    Expert options:
//    look at AddESDevent() and AddCosmicEvent() to get the idea of how the
//    aligner works, it's safe to repeat the needed steps outside of the class,
//    only public methods are used.
//
//    Origin: Mikolaj Krzewicki, Nikhef, Mikolaj.Krzewicki@cern.ch
//
//////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <TObject.h>
#include <TMath.h>
#include <TMatrix.h>
#include <TVector.h>
#include <TVector3.h>
#include <TDecompLU.h>
#include <TArrayI.h>
#include <TObjArray.h>
#include <TH1D.h>
#include <TF1.h>

#include "AliESDtrack.h"
#include "AliESDEvent.h"
#include "AliExternalTrackParam.h"

#include "AliRelAlignerKalman.h"

ClassImp(AliRelAlignerKalman)

//______________________________________________________________________________
AliRelAlignerKalman::AliRelAlignerKalman():
    TObject(),
    fPTrackParam1(new AliExternalTrackParam()),
    fPTrackParam2(new AliExternalTrackParam()),
    fMagField(0.),
    fNMeasurementParams(4),
    fPX(new TVectorD( fgkNSystemParams )),
    fPXcov(new TMatrixDSym( fgkNSystemParams )),
    fPH(new TMatrixD( fNMeasurementParams, fgkNSystemParams )),
    fQ(1.e-15),
    fPMeasurement(new TVectorD( fNMeasurementParams )),
    fPMeasurementCov(new TMatrixDSym( fNMeasurementParams )),
    fPMeasurementPrediction(new TVectorD( fNMeasurementParams )),
    fOutRejSigmas(1.),
    fOutRejSigma2Median(5.),
    fYZOnly(kFALSE),
    fNumericalParanoia(kTRUE),
    fRejectOutliers(kTRUE),
    fRejectOutliersSigma2Median(kFALSE),
    fRequireMatchInTPC(kFALSE),
    fCuts(kFALSE),
    fMinPointsVol1(3),
    fMinPointsVol2(50),
    fMinPt(0.),
    fMaxPt(1.e100),
    fMaxMatchingAngle(0.1),
    fMaxMatchingDistance(10.),  //in cm
    fCorrectionMode(kFALSE),
    fNTracks(0),
    fNUpdates(0),
    fNOutliers(0),
    fNOutliersSigma2Median(0),
    fNMatchedCosmics(0),
    fNMatchedTPCtracklets(0),
    fNProcessedEvents(0),
    fTimeStamp(0),
    fRunNumber(0),
    fNMerges(0),
    fNMergesFailed(0),
    fTPCvd(2.64),
    fTPCZLengthA(2.4972500e02),
    fTPCZLengthC(2.4969799e02)
{
  //Default constructor
  for (Int_t i=0;i<fgkNSystemParams;i++) fDelta[i] = 1.e-6;
  for (Int_t i=0; i<4;i++){fResArrSigma2Median[i]=NULL;}
  Reset();
}

//______________________________________________________________________________
AliRelAlignerKalman::AliRelAlignerKalman(const AliRelAlignerKalman& a):
    TObject(static_cast<TObject>(a)),
    fPTrackParam1(new AliExternalTrackParam()),
    fPTrackParam2(new AliExternalTrackParam()),
    fMagField(a.fMagField),
    fNMeasurementParams(a.fNMeasurementParams),
    fPX(new TVectorD( *a.fPX )),
    fPXcov(new TMatrixDSym( *a.fPXcov )),
    fPH(new TMatrixD( fNMeasurementParams, fgkNSystemParams )),
    fQ(a.fQ),
    fPMeasurement(new TVectorD( fNMeasurementParams )),
    fPMeasurementCov(new TMatrixDSym( fNMeasurementParams )),
    fPMeasurementPrediction(new TVectorD( fNMeasurementParams )),
    fOutRejSigmas(a.fOutRejSigmas),
    fOutRejSigma2Median(a.fOutRejSigma2Median),
    fYZOnly(a.fYZOnly),
    fNumericalParanoia(a.fNumericalParanoia),
    fRejectOutliers(a.fRejectOutliers),
    fRejectOutliersSigma2Median(a.fRejectOutliersSigma2Median),
    fRequireMatchInTPC(a.fRequireMatchInTPC),
    fCuts(a.fCuts),
    fMinPointsVol1(a.fMinPointsVol1),
    fMinPointsVol2(a.fMinPointsVol2),
    fMinPt(a.fMinPt),
    fMaxPt(a.fMaxPt),
    fMaxMatchingAngle(a.fMaxMatchingAngle),
    fMaxMatchingDistance(a.fMaxMatchingDistance),  //in cm
    fCorrectionMode(a.fCorrectionMode),
    fNTracks(a.fNTracks),
    fNUpdates(a.fNUpdates),
    fNOutliers(a.fNOutliers),
    fNOutliersSigma2Median(a.fNOutliersSigma2Median),
    fNMatchedCosmics(a.fNMatchedCosmics),
    fNMatchedTPCtracklets(a.fNMatchedTPCtracklets),
    fNProcessedEvents(a.fNProcessedEvents),
    fTimeStamp(a.fTimeStamp),
    fRunNumber(a.fRunNumber),
    fNMerges(a.fNMerges),
    fNMergesFailed(a.fNMergesFailed),
    fTPCvd(a.fTPCvd),
    fTPCZLengthA(a.fTPCZLengthA),
    fTPCZLengthC(a.fTPCZLengthC)
{
  //copy constructor
  memcpy(fDelta,a.fDelta,fgkNSystemParams*sizeof(Double_t));

  //copy contents of the residuals array for sigma2median scheme
  for (Int_t i=0;i<4;i++)
  {
    if ((a.fResArrSigma2Median)[i]) 
    {
      fResArrSigma2Median[i] = new Double_t[fgkNtracksSigma2Median];
      memcpy(fResArrSigma2Median[i],(a.fResArrSigma2Median)[i],
             fgkNtracksSigma2Median*sizeof(Double_t));
    }
    else
      fResArrSigma2Median[i] = NULL;
  }
}

//______________________________________________________________________________
AliRelAlignerKalman& AliRelAlignerKalman::operator=(const AliRelAlignerKalman& a)
{
  //assignment operator
  fMagField=a.fMagField,
  fNMeasurementParams=a.fNMeasurementParams;
  *fPX = *a.fPX;
  *fPXcov = *a.fPXcov;
  fQ=a.fQ;
  fOutRejSigmas=a.fOutRejSigmas;
  fOutRejSigma2Median=a.fOutRejSigma2Median;
  memcpy(fDelta,a.fDelta,fgkNSystemParams*sizeof(Double_t));
  fYZOnly=a.fYZOnly;
  fNumericalParanoia=a.fNumericalParanoia;
  fRejectOutliers=a.fRejectOutliers;
  fRejectOutliersSigma2Median=a.fRejectOutliersSigma2Median;
  fRequireMatchInTPC=a.fRequireMatchInTPC;
  fCuts=a.fCuts;
  fMinPointsVol1=a.fMinPointsVol1;
  fMinPointsVol2=a.fMinPointsVol2;
  fMinPt=a.fMinPt;
  fMaxPt=a.fMaxPt;
  fMaxMatchingAngle=a.fMaxMatchingAngle;
  fMaxMatchingDistance=a.fMaxMatchingDistance;  //in c;
  fCorrectionMode=a.fCorrectionMode;
  fNTracks=a.fNTracks;
  fNUpdates=a.fNUpdates;
  fNOutliers=a.fNOutliers;
  fNOutliersSigma2Median=a.fNOutliersSigma2Median;
  fNMatchedCosmics=a.fNMatchedCosmics;
  fNMatchedTPCtracklets=a.fNMatchedTPCtracklets;
  fNProcessedEvents=a.fNProcessedEvents;
  fTimeStamp=a.fTimeStamp;
  fRunNumber=a.fRunNumber;
  fNMerges=a.fNMerges;
  fTPCvd=a.fTPCvd;
  fTPCZLengthA=a.fTPCZLengthA;
  fTPCZLengthC=a.fTPCZLengthC;

  //copy contents of the residuals array for sigma2median scheme
  for (Int_t i=0;i<4;i++)
  {
    if ((a.fResArrSigma2Median)[i]) 
    {
      if (!(fResArrSigma2Median[i])) fResArrSigma2Median[i] = 
                                     new Double_t[fgkNtracksSigma2Median];
      memcpy(fResArrSigma2Median[i],(a.fResArrSigma2Median)[i],
             fgkNtracksSigma2Median*sizeof(Double_t));
    }
    else
      fResArrSigma2Median[i] = NULL;
  }
  return *this;
}

//______________________________________________________________________________
AliRelAlignerKalman::~AliRelAlignerKalman()
{
  //destructor
  delete fPTrackParam1;
  delete fPTrackParam2;
  delete fPX;
  delete fPXcov;
  delete fPH;
  delete fPMeasurement;
  delete fPMeasurementCov;
  for (Int_t i=0;i<4;i++) 
  {
    delete [] (fResArrSigma2Median[i]);
  }
}

//______________________________________________________________________________
Bool_t AliRelAlignerKalman::AddESDevent( const AliESDEvent* pEvent )
{
  //Add all tracks in an ESD event

  fNProcessedEvents++; //update the counter
  
  Bool_t success=kFALSE;
  SetMagField( pEvent->GetMagneticField() );
  AliESDtrack* track=NULL;
  
  for (Int_t i=0; i<pEvent->GetNumberOfTracks(); i++)
  {
    track = pEvent->GetTrack(i);
    if (!track) continue;
    if ( ((track->GetStatus()&AliESDtrack::kTPCin)>0)&&
         ((track->GetStatus()&AliESDtrack::kITSrefit)>0)&&
         (track->GetNcls(0)>=fMinPointsVol1)&&
         (track->GetNcls(1)>=fMinPointsVol2) )
    { 
      success = ( AddESDtrack( track ) || success );
    }
  }
  if (success)
  {
    fTimeStamp = pEvent->GetTimeStamp();
    fRunNumber = pEvent->GetRunNumber();
  }
  return success;
}

//______________________________________________________________________________
Bool_t AliRelAlignerKalman::AddESDtrack( const AliESDtrack* pTrack )
{
  //Adds a full track, returns true if results in a new estimate
  //  gets the inner TPC parameters from AliESDTrack::GetInnerParam()
  //  gets the outer ITS parameters from AliESDfriendTrack::GetITSout()

  const AliExternalTrackParam* pconstparamsITS = pTrack->GetOuterParam();
  if (!pconstparamsITS) return kFALSE;
  const AliExternalTrackParam* pconstparamsTPC = pTrack->GetInnerParam();
  if (!pconstparamsTPC) return kFALSE;
  
  //TPC part
  AliExternalTrackParam paramsTPC = (*pconstparamsTPC);
  paramsTPC.Rotate(pconstparamsITS->GetAlpha());
  paramsTPC.PropagateTo(pconstparamsITS->GetX(), fMagField);

  return (AddTrackParams(pconstparamsITS, &paramsTPC));
}

//______________________________________________________________________________
Bool_t AliRelAlignerKalman::AddTrackParams( const AliExternalTrackParam* p1, const AliExternalTrackParam* p2 )
{
  //Update the estimate using new matching tracklets

  if (!SetTrackParams(p1, p2)) return kFALSE;
  return Update();
}

//______________________________________________________________________________
Bool_t AliRelAlignerKalman::AddCosmicEvent( const AliESDEvent* pEvent )
{
  //Add an cosmic with separately tracked ITS and TPC parts, do trackmatching

  fNProcessedEvents++; //update the counter

  Bool_t success=kFALSE;
  TArrayI trackTArrITS(1);
  TArrayI trackTArrTPC(1);
  if (!FindCosmicTrackletNumbersInEvent( trackTArrITS, trackTArrTPC, pEvent )) return kFALSE;
  SetMagField( pEvent->GetMagneticField() );
  AliESDtrack* ptrack=NULL;
  const AliExternalTrackParam* pconstparams1;
  const AliExternalTrackParam* pconstparams2;
  AliExternalTrackParam params1;
  AliExternalTrackParam params2;
  
  ////////////////////////////////
  for (Int_t i=0;i<trackTArrITS.GetSize();i++)
  {
    //ITS track
    ptrack = pEvent->GetTrack(trackTArrITS[i]);
    pconstparams1 = ptrack->GetOuterParam();
    if (!pconstparams1) continue;
    params1 = *pconstparams1; //make copy to be safe
    
    //TPC track
    ptrack = pEvent->GetTrack(trackTArrTPC[i]);
    pconstparams2 = ptrack->GetInnerParam();
    if (!pconstparams2) continue;
    params2 = *pconstparams2; //make copy
    params2.Rotate(params1.GetAlpha());
    params2.PropagateTo( params1.GetX(), fMagField );

    if (!SetTrackParams( &params1, &params2 )) continue;
    
    //do some accounting and update
    if (Update())
      success = kTRUE;
    else
      continue;
  }
  fTimeStamp=pEvent->GetTimeStamp(); //always update timestamp even when no update performed
  fRunNumber=pEvent->GetRunNumber();
  return success;
}

//______________________________________________________________________________
void AliRelAlignerKalman::SetPoint2Track( Bool_t set )
{
  fNMeasurementParams = (set)?2:4;
  delete fPH;
  fPH = new TMatrixD( fNMeasurementParams, fgkNSystemParams );
  delete fPMeasurement;
  fPMeasurement = new TVectorD( fNMeasurementParams );
  delete fPMeasurementCov;
  fPMeasurementCov = new TMatrixDSym( fNMeasurementParams );
  delete fPMeasurementPrediction;
  fPMeasurementPrediction = new TVectorD( fNMeasurementParams );
  fYZOnly = set;
}

//______________________________________________________________________________
void AliRelAlignerKalman::Print(Option_t*) const
{
  //Print some useful info
  Double_t rad2deg = 180./TMath::Pi();
  printf("\nAliRelAlignerKalman\n");
  if (fCorrectionMode) printf("(Correction mode)\n");
  printf("  run: %i, timestamp: %i, magfield: %.3f\n", fRunNumber, fTimeStamp, fMagField);
  printf("  %i(-%i) inputs, %i(-%i) updates, %i(-%i) merges\n", fNTracks, fNOutliersSigma2Median, fNUpdates, fNOutliers, fNMerges, fNMergesFailed );
  printf("  psi(x):           % .3f ± (%.2f) mrad  |  % .3f ± (%.2f) deg\n",1e3*(*fPX)(0), 1e3*TMath::Sqrt((*fPXcov)(0,0)),(*fPX)(0)*rad2deg,TMath::Sqrt((*fPXcov)(0,0))*rad2deg);
  printf("  theta(y):         % .3f ± (%.2f) mrad  |  % .3f ± (%.2f) deg\n",1e3*(*fPX)(1), 1e3*TMath::Sqrt((*fPXcov)(1,1)),(*fPX)(1)*rad2deg,TMath::Sqrt((*fPXcov)(1,1))*rad2deg);
  printf("  phi(z):           % .3f ± (%.2f) mrad  |  % .3f ± (%.2f) deg\n",1e3*(*fPX)(2), 1e3*TMath::Sqrt((*fPXcov)(2,2)),(*fPX)(2)*rad2deg,TMath::Sqrt((*fPXcov)(2,2))*rad2deg);
  printf("  x:                % .3f ± (%.2f) micron\n", 1e4*(*fPX)(3), 1e4*TMath::Sqrt((*fPXcov)(3,3)));
  printf("  y:                % .3f ± (%.2f) micron\n", 1e4*(*fPX)(4), 1e4*TMath::Sqrt((*fPXcov)(4,4)));
  printf("  z:                % .3f ± (%.2f) micron\n", 1e4*(*fPX)(5), 1e4*TMath::Sqrt((*fPXcov)(5,5)));
  if (fgkNSystemParams>6) printf("  vd corr           % .5g ± (%.2g)    [ vd should be %.4g (was %.4g in reco) ]\n", (*fPX)(6), TMath::Sqrt((*fPXcov)(6,6)), (*fPX)(6)*fTPCvd, fTPCvd);
  if (fgkNSystemParams>7) printf("  t0                % .5g ± (%.2g) us  |  %.4g ± (%.2g) cm     [ t0_real = t0_rec+t0 ]\n",(*fPX)(7), TMath::Sqrt((*fPXcov)(7,7)), fTPCvd*(*fPX)(7), fTPCvd*TMath::Sqrt((*fPXcov)(7,7)));
  if (fgkNSystemParams>8) printf("  vd/dy             % .5f ± (%.2f) (cm/us)/m\n", (*fPX)(8), TMath::Sqrt((*fPXcov)(8,8)));
  printf("\n");
  return;
}

//______________________________________________________________________________
void AliRelAlignerKalman::PrintSystemMatrix()
{
  //Print the system matrix for this measurement
  printf("Kalman system matrix:\n");
  for ( Int_t i=0; i<fNMeasurementParams; i++ )
  {
    for ( Int_t j=0; j<fgkNSystemParams; j++ )
    {
      printf("% -2.2f  ", (*fPH)(i,j) );
    }//for i
    printf("\n");
  }//for j
  printf("\n");
  return;
}

//______________________________________________________________________________
Bool_t AliRelAlignerKalman::SetTrackParams( const AliExternalTrackParam* exparam1, const AliExternalTrackParam* exparam2 )
{
  //Set the parameters, exparam1 will normally be ITS and exparam 2 tht TPC
  fNTracks++; //count added input sets

  //INPUT OUTLIER REJECTION
  if (fRejectOutliersSigma2Median && IsOutlierSigma2Median(exparam1,exparam2) ) return kFALSE;

  *fPTrackParam1 = *exparam1;
  *fPTrackParam2 = *exparam2;
  
  return kTRUE;
}

//______________________________________________________________________________
Bool_t AliRelAlignerKalman::Update()
{
  //perform the update
  
  //if (fCalibrationMode) return UpdateCalibration();
  //if (fFillHistograms)
  //{
  //  if (!UpdateEstimateKalman()) return kFALSE;
  //  return UpdateCalibration(); //Update histograms only when update ok.
  //}
  //else return UpdateEstimateKalman();
  if (!PrepareMeasurement()) return kFALSE;
  if (!PrepareSystemMatrix()) return kFALSE;
  if (!PreparePrediction()) return kFALSE;
  return UpdateEstimateKalman();
}

//______________________________________________________________________________
void AliRelAlignerKalman::RotMat( TMatrixD &R, const TVectorD& angles )
{
  //Get Rotation matrix R given the Cardan angles psi, theta, phi (around x, y, z).
  Double_t sinpsi = TMath::Sin(angles(0));
  Double_t sintheta = TMath::Sin(angles(1));
  Double_t sinphi = TMath::Sin(angles(2));
  Double_t cospsi = TMath::Cos(angles(0));
  Double_t costheta = TMath::Cos(angles(1));
  Double_t cosphi = TMath::Cos(angles(2));

  R(0,0) = costheta*cosphi;
  R(0,1) = -costheta*sinphi;
  R(0,2) = sintheta;
  R(1,0) = sinpsi*sintheta*cosphi + cospsi*sinphi;
  R(1,1) = -sinpsi*sintheta*sinphi + cospsi*cosphi;
  R(1,2) = -costheta*sinpsi;
  R(2,0) = -cospsi*sintheta*cosphi + sinpsi*sinphi;
  R(2,1) = cospsi*sintheta*sinphi + sinpsi*cosphi;
  R(2,2) = costheta*cospsi;
}

//______________________________________________________________________________
Bool_t AliRelAlignerKalman::PrepareMeasurement()
{
  //Calculate the residuals and their covariance matrix
  
    const Double_t* pararr1 = fPTrackParam1->GetParameter();
    const Double_t* pararr2 = fPTrackParam2->GetParameter();

    //Take the track parameters and calculate the input to the Kalman filter
    (*fPMeasurement)(0) = pararr2[0]-pararr1[0];
    (*fPMeasurement)(1) = pararr2[1]-pararr1[1];
    if (!fYZOnly)
    {
      (*fPMeasurement)(2) = pararr2[2]-pararr1[2];
      (*fPMeasurement)(3) = pararr2[3]-pararr1[3];
    }

    //the covariance
    const Double_t* parcovarr1 = fPTrackParam1->GetCovariance();
    const Double_t* parcovarr2 = fPTrackParam2->GetCovariance();
    (*fPMeasurementCov)(0,0)=parcovarr1[0];
    (*fPMeasurementCov)(0,1)=parcovarr1[1];
    (*fPMeasurementCov)(1,0)=parcovarr1[1];
    (*fPMeasurementCov)(1,1)=parcovarr1[2];
    (*fPMeasurementCov)(0,0)+=parcovarr2[0];
    (*fPMeasurementCov)(0,1)+=parcovarr2[1];
    (*fPMeasurementCov)(1,0)+=parcovarr2[1];
    (*fPMeasurementCov)(1,1)+=parcovarr2[2];
    if (!fYZOnly)
    {
      (*fPMeasurementCov)(0,2)=parcovarr1[3];
      (*fPMeasurementCov)(0,3)=parcovarr1[6];
      (*fPMeasurementCov)(1,2)=parcovarr1[4];
      (*fPMeasurementCov)(1,3)=parcovarr1[7];
      (*fPMeasurementCov)(2,0)=parcovarr1[3];
      (*fPMeasurementCov)(2,1)=parcovarr1[4];
      (*fPMeasurementCov)(2,2)=parcovarr1[5];
      (*fPMeasurementCov)(2,3)=parcovarr1[8];
      (*fPMeasurementCov)(3,0)=parcovarr1[6];
      (*fPMeasurementCov)(3,1)=parcovarr1[7];
      (*fPMeasurementCov)(3,2)=parcovarr1[8];
      (*fPMeasurementCov)(3,3)=parcovarr1[9];
      (*fPMeasurementCov)(0,2)+=parcovarr2[3];
      (*fPMeasurementCov)(0,3)+=parcovarr2[6];
      (*fPMeasurementCov)(1,2)+=parcovarr2[4];
      (*fPMeasurementCov)(1,3)+=parcovarr2[7];
      (*fPMeasurementCov)(2,0)+=parcovarr2[3];
      (*fPMeasurementCov)(2,1)+=parcovarr2[4];
      (*fPMeasurementCov)(2,2)+=parcovarr2[5];
      (*fPMeasurementCov)(2,3)+=parcovarr2[8];
      (*fPMeasurementCov)(3,0)+=parcovarr2[6];
      (*fPMeasurementCov)(3,1)+=parcovarr2[7];
      (*fPMeasurementCov)(3,2)+=parcovarr2[8];
      (*fPMeasurementCov)(3,3)+=parcovarr2[9];
    }
    
  //if (fApplyCovarianceCorrection)
  //  *fPMeasurementCov += *fPMeasurementCovCorr;
  return kTRUE;
}

//______________________________________________________________________________
Bool_t AliRelAlignerKalman::PrepareSystemMatrix()
{
  //Calculate the system matrix for the Kalman filter
  //approximate the system using as reference the track in the first volume

  TVectorD z1( fNMeasurementParams );
  TVectorD z2( fNMeasurementParams );
  TVectorD x1( fgkNSystemParams );
  TVectorD x2( fgkNSystemParams );
  //get the derivatives
  for ( Int_t i=0; i<fgkNSystemParams; i++ )
  {
    x1 = *fPX;
    x2 = *fPX;
    x1(i) = x1(i) - fDelta[i]/(2.0);
    x2(i) = x2(i) + fDelta[i]/(2.0);
    if (!PredictMeasurement( z1, x1 )) return kFALSE;
    if (!PredictMeasurement( z2, x2 )) return kFALSE;
    for (Int_t j=0; j<fNMeasurementParams; j++ )
    {
      (*fPH)(j,i) = ( z2(j)-z1(j) ) / fDelta[i];
    }
  }
  return kTRUE;
}

//______________________________________________________________________________
Bool_t AliRelAlignerKalman::PreparePrediction()
{
  //Prepare the prediction of the measurement using state vector
  return PredictMeasurement( (*fPMeasurementPrediction), (*fPX) );
}

//______________________________________________________________________________
Bool_t AliRelAlignerKalman::PredictMeasurement( TVectorD& pred, const TVectorD& state )
{
  // Implements a system model for the Kalman fit
  // pred is [dy,dz,dsinphi,dtgl]
  // state is [psi,theta,phi,x,y,z,driftTPC,offsetTPC]
  // note: the measurement is in a local frame, so the prediction also has to be
  // note: state is the misalignment in global reference system

  if (fCorrectionMode)
  {
      AliExternalTrackParam track(*fPTrackParam2); //make a copy track
      if (!CorrectTrack( &track, state )) return kFALSE; //predict what the ideal track would be by applying correction
      
      const Double_t* oldparam = fPTrackParam2->GetParameter();
      const Double_t* newparam = track.GetParameter();

      //calculate the predicted residual
      pred(0) = oldparam[0] - newparam[0];
      pred(1) = oldparam[1] - newparam[1];
      if (!fYZOnly)
      {
        pred(2) = oldparam[2] - newparam[2];
        pred(3) = oldparam[3] - newparam[3];
      }
      return kTRUE;
  }
  else
  {
      AliExternalTrackParam track(*fPTrackParam1); //make a copy track
      if (!MisalignTrack( &track, state )) return kFALSE; //predict what the measured track would be by applying misalignment

      const Double_t* oldparam = fPTrackParam1->GetParameter();
      const Double_t* newparam = track.GetParameter();

      //calculate the predicted residual
      pred(0) = newparam[0] - oldparam[0];
      pred(1) = newparam[1] - oldparam[1];
      if (!fYZOnly)
      {
        pred(2) = newparam[2] - oldparam[2];
        pred(3) = newparam[3] - oldparam[3];
      }
      return kTRUE;
    }
  return kFALSE;
}

//______________________________________________________________________________
Bool_t AliRelAlignerKalman::UpdateEstimateKalman()
{
  //Kalman estimation of noisy constants: in the model A=1
  //The arguments are (following the usual convention):
  //  fPX - the state vector (parameters)
  //  fPXcov - the state covariance matrix (parameter errors)
  //  fPMeasurement - measurement vector
  //  fPMeasurementCov - measurement covariance matrix
  //  fPH - measurement model matrix ( fPMeasurement = Hx + v ) v being measurement noise (error fR)

  TMatrixDSym identity(TMatrixDSym::kUnit, (*fPXcov));            //unit matrix

  //predict the state
  //(*fPXcov) = (*fPXcov) + fQ*identity;  //add some process noise (diagonal)

  // update prediction with measurement
  // calculate Kalman gain
  // K = PH/(HPH+fPMeasurementCov)
  TMatrixD pht( (*fPXcov), TMatrixD::kMultTranspose, (*fPH) );  //common factor (used twice)
  TMatrixD hpht( (*fPH), TMatrixD::kMult, pht );
  hpht += (*fPMeasurementCov);

  //shit happens so protect yourself!
//  if (fNumericalParanoia)
//  {
//    TDecompLU lu(hpht);
//    if (lu.Condition() > 1e12) return kFALSE;
//    lu.Invert(hpht);
//  }
//  else
//  {
    Double_t det=0.0;
    hpht.Invert(&det); //since the matrix is small...
    if (det < 2e-99) return kFALSE; //we need some sort of protection even in this case....
//  }
  //printf("KalmanUpdate: det(hpht): %.4g\n",det);

  TMatrixD k(pht, TMatrixD::kMult, hpht ); //compute K (hpht is already inverted)

  // update the state and its covariance matrix
  TVectorD xupdate(fgkNSystemParams);
  xupdate = k*((*fPMeasurement)-(*fPMeasurementPrediction));

  //SIMPLE OUTLIER REJECTION
  if ( IsOutlier( xupdate, (*fPXcov) ) && fRejectOutliers )
  {
    fNOutliers++;
    //printf("AliRelAlignerKalman: outlier\n");
    return kFALSE;
  }

  TMatrixD kh( k, TMatrixD::kMult, (*fPH) );
  TMatrixD ikh(fgkNSystemParams,fgkNSystemParams); //this is because for some reason TMatrixD::kAdd didn't work
  ikh = identity - kh;
  TMatrixD ikhp( ikh, TMatrixD::kMult, (*fPXcov) ); // (identity-KH)fPXcov
  if (!IsPositiveDefinite(ikhp)) return kFALSE;

  (*fPX) += xupdate;
  TMatrixDSymFromTMatrixD( (*fPXcov), ikhp ); //make the matrix completely symetrical

  fNUpdates++;

  return kTRUE;
}

//______________________________________________________________________________
Bool_t AliRelAlignerKalman::IsOutlier( const TVectorD& update, const TMatrixDSym& covmatrix )
{
  //check whether an update is an outlier given the covariance matrix of the fit

  Bool_t is=kFALSE;
  for (Int_t i=0;i<fgkNSystemParams;i++)
  {
    if (covmatrix(i,i)<0.) return kTRUE; //if cov matrix has neg diagonals something went wrong
    is = (is) || (TMath::Abs(update(i)) > fOutRejSigmas*TMath::Sqrt((covmatrix)(i,i)));
  }
  return is;
}

//______________________________________________________________________________
Bool_t AliRelAlignerKalman::IsOutlierSigma2Median( const AliExternalTrackParam* pITS, 
                                                   const AliExternalTrackParam* pTPC )
{
  //check if the input residuals are not too far off their median
  TVectorD vecDelta(4),vecMedian(4), vecRMS(4);
  TVectorD vecDeltaN(5);
  Double_t sign=(pITS->GetParameter()[1]>0)? 1.:-1.;
  vecDeltaN[4]=0;
  for (Int_t i=0;i<4;i++){
    vecDelta[i]=(pITS->GetParameter()[i]-pTPC->GetParameter()[i])*sign;
    (fResArrSigma2Median[i])[(fNTracks-1)%fgkNtracksSigma2Median]=vecDelta[i];
  }
  Int_t entries=(fNTracks<fgkNtracksSigma2Median)?fNTracks:fgkNtracksSigma2Median;
  for (Int_t i=0;i<fNMeasurementParams;i++){       //in point2trackmode just take the first 2 params (zy)
    vecMedian[i] = TMath::Median(entries,fResArrSigma2Median[i]);
    vecRMS[i]    = TMath::RMS(entries,fResArrSigma2Median[i]);
    vecDeltaN[i] = 0;
    if (vecRMS[i]>0.){
      vecDeltaN[i] = (vecDelta[i]-vecMedian[i])/vecRMS[i];
      vecDeltaN[4]+= TMath::Abs(vecDeltaN[i]);  //sum of abs residuals
    }
  }
  Bool_t outlier=kFALSE;
  if (fNTracks<3)  outlier=kTRUE;   //median and RMS still to be defined
  if ( vecDeltaN[4]/fNMeasurementParams>fOutRejSigma2Median) outlier=kTRUE;
  if (outlier) fNOutliersSigma2Median++;
  return outlier;
}

//______________________________________________________________________________
Bool_t AliRelAlignerKalman::IsPositiveDefinite( const TMatrixD& mat ) const
{
  //check for positive definiteness

  for (Int_t i=0; i<mat.GetNcols(); i++)
  {
    if (mat(i,i)<=0.) return kFALSE;
  }

  if (!fNumericalParanoia) return kTRUE;

  TDecompLU lu(mat);
  return (lu.Decompose());
}

//______________________________________________________________________________
void AliRelAlignerKalman::TMatrixDSymFromTMatrixD( TMatrixDSym& matsym, const TMatrixD& mat )
{
  //Produce a valid symmetric matrix out of an almost symmetric TMatrixD

  for (Int_t i=0; i<mat.GetNcols(); i++)
  {
    matsym(i,i) = mat(i,i); //copy diagonal
    for (Int_t j=i+1; j<mat.GetNcols(); j++)
    {
      //copy the rest
      Double_t average = (mat(i,j)+mat(j,i))/2.;
      matsym(i,j)=average;
      matsym(j,i)=average;
    }
  }
  matsym.MakeValid();
  return;
}

//______________________________________________________________________________
void AliRelAlignerKalman::Angles( TVectorD &angles, const TMatrixD &rotmat )
{
  //Calculate the Cardan angles (psi,theta,phi) from rotation matrix
  //b = R*a
  angles(0) = TMath::ATan2( -rotmat(1,2), rotmat(2,2) );
  angles(1) = TMath::ASin( rotmat(0,2) );
  angles(2) = TMath::ATan2( -rotmat(0,1), rotmat(0,0) );
  return;
}

//______________________________________________________________________________
void AliRelAlignerKalman::PrintCorrelationMatrix()
{
  //Print the correlation matrix for the fitted parameters
  printf("Correlation matrix for system parameters:\n");
  for ( Int_t i=0; i<fgkNSystemParams; i++ )
  {
    for ( Int_t j=0; j<i+1; j++ )
    {
      if ((*fPXcov)(i,i)==0. || (*fPXcov)(j,j)==0.) printf("   NaN  ");
      else
        printf("% -1.3f  ", (*fPXcov)(i,j)/TMath::Sqrt( (*fPXcov)(i,i) * (*fPXcov)(j,j) ) );
    }//for j
    printf("\n");
  }//for i
  printf("\n");
  return;
}

//______________________________________________________________________________
Bool_t AliRelAlignerKalman::FindCosmicTrackletNumbersInEvent( TArrayI& outITSindexTArr, TArrayI& outTPCindexTArr, const AliESDEvent* pEvent )
{
  //Find matching track segments in an event with tracks in TPC and ITS(standalone)

  //Sanity cuts on tracks + check which tracks are ITS which are TPC
  Int_t ntracks = pEvent->GetNumberOfTracks(); ////printf("number of tracks in event: %i\n", ntracks);
  fMagField = pEvent->GetMagneticField();
  if (ntracks<2)
  {
    //printf("TrackFinder: less than 2 tracks!\n");
    return kFALSE;
  }
  Float_t* phiArr = new Float_t[ntracks];
  Float_t* thetaArr = new Float_t[ntracks];
  Int_t* goodtracksArr = new Int_t[ntracks];
  Int_t* candidateTPCtracksArr = new Int_t[ntracks];
  Int_t* matchedITStracksArr = new Int_t[ntracks];
  Int_t* matchedTPCtracksArr = new Int_t[ntracks];
  Int_t* tracksArrITS = new Int_t[ntracks];
  Int_t* tracksArrTPC = new Int_t[ntracks];
  Int_t* nPointsArr = new Int_t[ntracks];
  Int_t nITStracks = 0;
  Int_t nTPCtracks = 0;
  Int_t nGoodTracks = 0;
  Int_t nCandidateTPCtracks = 0;
  Int_t nMatchedITStracks = 0;
  AliESDtrack* pTrack = NULL;
  Bool_t foundMatchTPC = kFALSE;

  //select and clasify tracks
  for (Int_t itrack=0; itrack < ntracks; itrack++)
  {
    pTrack = pEvent->GetTrack(itrack);
    if (!pTrack)
    {
      //std::cout<<"no track!"<<std::endl;
      continue;
    }
    if (fCuts)
    {
      if (pTrack->GetP()<fMinPt || pTrack->GetP()>fMaxPt) continue;
    }
    goodtracksArr[nGoodTracks]=itrack;
    Float_t phi = pTrack->GetAlpha()+TMath::ASin(pTrack->GetSnp());
    Float_t theta = 0.5*TMath::Pi()-TMath::ATan(pTrack->GetTgl());
    phiArr[nGoodTracks]=phi;
    thetaArr[nGoodTracks]=theta;

    //check if track is ITS
    Int_t nClsITS = pTrack->GetNcls(0);
    Int_t nClsTPC = pTrack->GetNcls(1);
    if ( ((pTrack->GetStatus()&AliESDtrack::kITSout)>0)&&
         !((pTrack->GetStatus()&AliESDtrack::kTPCin)>0)&&
         !(nClsITS<fMinPointsVol1) )  //enough points
    {
      tracksArrITS[nITStracks] = nGoodTracks;
      nITStracks++;
      nGoodTracks++;
      continue;
    }

    //check if track is TPC
    if ( ((pTrack->GetStatus()&AliESDtrack::kTPCin)>0)&&
         !(nClsTPC<fMinPointsVol2) )  //enough points
    {
      tracksArrTPC[nTPCtracks] = nGoodTracks;
      nTPCtracks++;
      nGoodTracks++;
      //printf("tracksArrTPC[%d]=%d, goodtracksArr[%d]=%d\n",nTPCtracks-1,tracksArrTPC[nTPCtracks-1],nGoodTracks-1,goodtracksArr[nGoodTracks-1]);
      continue;
    }
  }//for itrack   -   selection fo tracks

  //printf("TrackFinder: %d ITS | %d TPC out of %d tracks in event\n", nITStracks,nTPCtracks,ntracks);

  if ( nITStracks < 1 || nTPCtracks < 1 )
  {
    delete [] phiArr;
    delete [] thetaArr;
    delete [] goodtracksArr;
    delete [] matchedITStracksArr;
    delete [] candidateTPCtracksArr;
    delete [] matchedTPCtracksArr;
    delete [] tracksArrITS;
    delete [] tracksArrTPC;
    delete [] nPointsArr;
    return kFALSE;
  }

  //find matching in TPC
  if (nTPCtracks>1)  //if there is something to be matched, try and match it
  {
    Float_t min = 10000000.;
    for (Int_t itr1=0; itr1<nTPCtracks; itr1++)
    {
      for (Int_t itr2=itr1+1; itr2<nTPCtracks; itr2++)
      {
        Float_t deltatheta = TMath::Abs(TMath::Pi()-thetaArr[tracksArrTPC[itr1]]-thetaArr[tracksArrTPC[itr2]]);
        if (deltatheta > fMaxMatchingAngle) continue;
        Float_t deltaphi = TMath::Abs(TMath::Abs(phiArr[tracksArrTPC[itr1]]-phiArr[tracksArrTPC[itr2]])-TMath::Pi());
        if (deltaphi > fMaxMatchingAngle) continue;
        if (deltatheta+deltaphi<min) //only the best matching pair
        {
          min=deltatheta+deltaphi;
          candidateTPCtracksArr[0] = tracksArrTPC[itr1];  //store the index of track in goodtracksArr[]
          candidateTPCtracksArr[1] = tracksArrTPC[itr2];
          nCandidateTPCtracks = 2;
          foundMatchTPC = kTRUE;
          //printf("TrackFinder: Matching TPC tracks candidates:\n");
          //printf("TrackFinder: candidateTPCtracksArr[0]=%d\n",tracksArrTPC[itr1]);
          //printf("TrackFinder: candidateTPCtracksArr[1]=%d\n",tracksArrTPC[itr2]);
        }
      }
    }
  }//if nTPCtracks>1
  else //if nTPCtracks==1 - if nothing to match, take the only one we've got
  {
    candidateTPCtracksArr[0] = tracksArrTPC[0];
    nCandidateTPCtracks = 1;
    foundMatchTPC = kFALSE;
  }
  if (foundMatchTPC) fNMatchedTPCtracklets++;
  //if no match but the requirement is set return kFALSE
  if (fRequireMatchInTPC && !foundMatchTPC)
  {
    delete [] phiArr;
    delete [] thetaArr;
    delete [] goodtracksArr;
    delete [] candidateTPCtracksArr;
    delete [] matchedITStracksArr;
    delete [] matchedTPCtracksArr;
    delete [] tracksArrITS;
    delete [] tracksArrTPC;
    delete [] nPointsArr;
    //printf("TrackFinder: no match in TPC && required\n");
    return kFALSE;
  }

  //if no match and more than one track take all TPC tracks
  if (!fRequireMatchInTPC && !foundMatchTPC)
  {
    for (Int_t i=0;i<nTPCtracks;i++)
    {
      candidateTPCtracksArr[i] = tracksArrTPC[i];
    }
    nCandidateTPCtracks = nTPCtracks;
  }
  //printf("TrackFinder: nCandidateTPCtracks: %i\n", nCandidateTPCtracks);

  Double_t* minDifferenceArr = new Double_t[nCandidateTPCtracks];

  //find ITS matches for good TPC tracks
  Bool_t matchedITStracks=kFALSE;
  for (Int_t itpc=0;itpc<nCandidateTPCtracks;itpc++)
  {
    minDifferenceArr[nMatchedITStracks] = 10000000.;
    matchedITStracks=kFALSE;
    for (Int_t iits=0; iits<nITStracks; iits++)
    {
      AliESDtrack* itstrack = pEvent->GetTrack(goodtracksArr[tracksArrITS[iits]]);
      const AliExternalTrackParam* parits = itstrack->GetOuterParam();
      AliESDtrack* tpctrack = pEvent->GetTrack(goodtracksArr[candidateTPCtracksArr[itpc]]);
      const AliExternalTrackParam* tmp = tpctrack->GetInnerParam();
      AliExternalTrackParam partpc(*tmp);  //make a copy to avoid tampering with original params
      partpc.Rotate(parits->GetAlpha());
      partpc.PropagateTo(parits->GetX(),fMagField);
      Float_t dtgl = TMath::Abs(partpc.GetTgl()-parits->GetTgl());
      if (dtgl > fMaxMatchingAngle) continue;
      Float_t dsnp = TMath::Abs(partpc.GetSnp()-parits->GetSnp());
      if (dsnp > fMaxMatchingAngle) continue;
      Float_t dy = TMath::Abs(partpc.GetY()-parits->GetY());
      Float_t dz = TMath::Abs(partpc.GetZ()-parits->GetZ());
      if (TMath::Sqrt(dy*dy+dz*dz) > fMaxMatchingDistance) continue;
      if (dtgl+dsnp<minDifferenceArr[nMatchedITStracks]) //only the best matching pair
      {
        minDifferenceArr[nMatchedITStracks]=dtgl+dsnp;
        matchedITStracksArr[nMatchedITStracks] = tracksArrITS[iits];
        matchedTPCtracksArr[nMatchedITStracks] = candidateTPCtracksArr[itpc]; //this tells us minDifferenceArrwhich TPC track this ITS track belongs to
        //printf("TrackFinder: Matching ITS to TPC:\n");
        //printf("TrackFinder: minDifferenceArr[%i]=%.2f\n",nMatchedITStracks,minDifferenceArr[nMatchedITStracks]);
        //printf("TrackFinder: matchedITStracksArr[%i]=%i\n",nMatchedITStracks,matchedITStracksArr[nMatchedITStracks]);
        //printf("TrackFinder: matchedTPCtracksArr[%i]=%i\n",nMatchedITStracks,matchedTPCtracksArr[nMatchedITStracks]);
        matchedITStracks=kTRUE;;
      }
    }
    if (matchedITStracks) nMatchedITStracks++;
  }

  if (nMatchedITStracks==0) //no match in ITS
  {
    delete [] phiArr;
    delete [] thetaArr;
    delete [] minDifferenceArr;
    delete [] goodtracksArr;
    delete [] matchedITStracksArr;
    delete [] candidateTPCtracksArr;
    delete [] matchedTPCtracksArr;
    delete [] tracksArrITS;
    delete [] tracksArrTPC;
    delete [] nPointsArr;
    //printf("TrackFinder: No match in ITS\n");
    return kFALSE;
  }

  //printf("TrackFinder: nMatchedITStracks: %i\n",nMatchedITStracks);
  //we found a cosmic
  fNMatchedCosmics++;

  //Now we may have ended up with more matches than we want in the case there was
  //no TPC match and there were many TPC tracks
  //a cosmic in a scenario like this will only have produced 1 pair, the rest is garbage
  //so take only the best pair.
  //same applies when there are more matches than ITS tracks - means that one ITS track
  //matches more TPC tracks.
  if ((nMatchedITStracks>2 && !foundMatchTPC) || nMatchedITStracks>nITStracks)
  {
    Int_t imin = TMath::LocMin(nMatchedITStracks,minDifferenceArr);
    matchedITStracksArr[0] = matchedITStracksArr[imin];
    matchedTPCtracksArr[0] = matchedTPCtracksArr[imin];
    nMatchedITStracks = 1;
    //printf("TrackFinder: too many matches - take only the best one\n");
    //printf("TrackFinder: LocMin in matched its tracks: %d\n",imin);
    //printf("TrackFinder: matchedITStracksArr[0]=%d\n",matchedITStracksArr[0]);
    //printf("TrackFinder: matchedTPCtracksArr[0]=%d\n",matchedTPCtracksArr[0]);
  }

  ///////////////////////////////////////////////////////////////////////////
  outITSindexTArr.Set(nMatchedITStracks);
  outTPCindexTArr.Set(nMatchedITStracks);
  for (Int_t i=0;i<nMatchedITStracks;i++)
  {
    outITSindexTArr.AddAt( goodtracksArr[matchedITStracksArr[i]], i );
    outTPCindexTArr.AddAt( goodtracksArr[matchedTPCtracksArr[i]], i );
    //printf("TrackFinder: Fill the output\n");
    //printf("TrackFinder: matchedITStracksArr[%d]=%d\n",i,matchedITStracksArr[i]);
    //printf("TrackFinder: matchedTPCtracksArr[%d]=%d\n",i,matchedTPCtracksArr[i]);
  }
  //printf("TrackFinder: Size of outputarrays: %d, %d\n", outITSindexTArr.GetSize(), outTPCindexTArr.GetSize());
  ///////////////////////////////////////////////////////////////////////////

  delete [] phiArr;
  delete [] thetaArr;
  delete [] minDifferenceArr;
  delete [] goodtracksArr;
  delete [] candidateTPCtracksArr;
  delete [] matchedITStracksArr;
  delete [] matchedTPCtracksArr;
  delete [] tracksArrITS;
  delete [] tracksArrTPC;
  delete [] nPointsArr;
  return kTRUE;
}

//______________________________________________________________________________
Bool_t AliRelAlignerKalman::CorrectTrack( AliExternalTrackParam* tr, const TVectorD& misal ) const
{
  //implements the system model -
  //applies correction for misalignment and calibration to track
  //track needs to be already propagated to the global reference plane

  Double_t x = tr->GetX();
  Double_t alpha = tr->GetAlpha();
  Double_t point[3],dir[3];
  tr->GetXYZ(point);
  tr->GetDirection(dir);
  TVector3 Point(point);
  TVector3 Dir(dir);
  
  //Apply corrections to track

  //Shift
  Point(0) -= misal(3); //add shift in x
  Point(1) -= misal(4); //add shift in y
  Point(2) -= misal(5); //add shift in z
  //Rotation
  TMatrixD rotmat(3,3);
  RotMat( rotmat, misal );
  Point = rotmat.T() * Point;
  Dir = rotmat * Dir;
  
  //TPC vdrift and T0 corrections
  TVector3 Point2; //second point of the track
  Point2 = Point + Dir;
  Double_t vdCorr = 1./misal(6);
  Double_t t0 = misal(7);
  Double_t vdY = 0.0;
  if (fgkNSystemParams>8) vdY = misal(8)/100.; //change over 100cm.

  //my model
  if (Point(2)>0)
  {
    //A-Side
    Point(2) = Point(2)   - (fTPCZLengthA-Point(2))  * (vdCorr-1.+vdY*Point(1)/fTPCvd)  - (fTPCvd*vdCorr+vdY*Point(1))*t0;
    Point2(2) = Point2(2) - (fTPCZLengthA-Point2(2)) * (vdCorr-1.+vdY*Point2(1)/fTPCvd) - (fTPCvd*vdCorr+vdY*Point2(1))*t0;
  }
  else
  {
    //C-side
    Point(2) = Point(2)   - (fTPCZLengthC+Point(2))  * (1.-vdCorr-vdY*Point(1)/fTPCvd)  + (fTPCvd*vdCorr+vdY*Point(1))*t0;
    Point2(2) = Point2(2) - (fTPCZLengthC+Point2(2)) * (1.-vdCorr-vdY*Point2(1)/fTPCvd) + (fTPCvd*vdCorr+vdY*Point2(1))*t0;
  }

  //Stefan's model
  //if (Point(2)>0)
  //{
  //  //A-Side
  //  Point(2) = Point(2)   - (fTPCZLengthA-Point(2))  * (1.-vdCorr+vdY*Point(1)/fTPCvd)  - (fTPCvd*vdCorr+vdY*Point(1))*t0;
  //  Point2(2) = Point2(2) - (fTPCZLengthA-Point2(2)) * (1.-vdCorr+vdY*Point2(1)/fTPCvd) - (fTPCvd*vdCorr+vdY*Point2(1))*t0;
  //}
  //else
  //{
  //  //C-side
  //  Point(2) = Point(2)   + (fTPCZLengthC+Point(2))  * (1.-vdCorr+vdY*Point(1)/fTPCvd)  + (fTPCvd*vdCorr+vdY*Point(1))*t0;
  //  Point2(2) = Point2(2) + (fTPCZLengthC+Point2(2)) * (1.-vdCorr+vdY*Point2(1)/fTPCvd) + (fTPCvd*vdCorr+vdY*Point2(1))*t0;
  //}

  Dir = Point2-Point;
  Dir=Dir.Unit(); //keep unit length

  //Turn back to local system
  Dir.GetXYZ(dir);
  Point.GetXYZ(point);
  tr->Global2LocalPosition(point,alpha);
  tr->Global2LocalPosition(dir,alpha);

  //Calculate new intersection point with ref plane
  Double_t p[5],pcov[15];
  if (dir[0]==0.0) return kFALSE;
  Double_t s=(x-point[0])/dir[0];
  p[0] = point[1]+s*dir[1];
  p[1] = point[2]+s*dir[2];
  Double_t pt = TMath::Sqrt(dir[0]*dir[0]+dir[1]*dir[1]);
  if (pt==0.0) return kFALSE;
  p[2] = dir[1]/pt;
  p[3] = dir[2]/pt;
  //insert everything back into track
  const Double_t* pcovtmp = tr->GetCovariance();
  p[4] = tr->GetSigned1Pt(); //copy the momentum
  memcpy(pcov,pcovtmp,15*sizeof(Double_t));
  tr->Set(x,alpha,p,pcov);
  return kTRUE;

  ////put params back into track and propagate to ref
  //Double_t p[5],pcov[15];
  //p[0] = point[1];
  //p[1] = point[2];
  //Double_t xnew = point[0];
  //Double_t pt = TMath::Sqrt(dir[0]*dir[0]+dir[1]*dir[1]);
  //if (pt==0.0) return kFALSE;
  //p[2] = dir[1]/pt;
  //p[3] = dir[2]/pt;
  //p[4] = tr->GetSigned1Pt(); //copy the momentum
  //const Double_t* pcovtmp = tr->GetCovariance();
  //memcpy(pcov,pcovtmp,15*sizeof(Double_t));
  //tr->Set(xnew,alpha,p,pcov);
  //return tr->PropagateTo(x,fMagField);
}

//______________________________________________________________________________
Bool_t AliRelAlignerKalman::MisalignTrack( AliExternalTrackParam* tr, const TVectorD& misal ) const
{
  //implements the system model -
  //applies misalignment and miscalibration to reference track
  //trackparams have to be at the global reference plane

  Double_t x = tr->GetX();
  Double_t alpha = tr->GetAlpha();
  Double_t point[3],dir[3];
  tr->GetXYZ(point);
  tr->GetDirection(dir);
  TVector3 Point(point);
  TVector3 Dir(dir);
  
  //Apply misalignment to track
  
  //TPC vdrift and T0 corrections
  TVector3 Point2; //second point of the track
  Point2 = Point + Dir;
  Double_t vdCorr = 1./misal(6);
  Double_t t0 = misal(7);
  Double_t vdY = 0.0;
  if (fgkNSystemParams>8) vdY = misal(8)/100.; //change over 100cm.

  if (Point(2)>0)
  {
    //A-Side
    Point(2) = Point(2)   + ((fTPCZLengthA-Point(2))/(vdCorr*fTPCvd+vdY*Point(1)))
                          * (fTPCvd*(vdCorr-1.)+vdY*Point(1)) + fTPCvd*t0;
    Point2(2) = Point2(2) + ((fTPCZLengthA-Point2(2))/(vdCorr*fTPCvd+vdY*Point2(1)))
                          * (fTPCvd*(vdCorr-1.)+vdY*Point2(1)) + fTPCvd*t0;
  }
  else
  {
    //C-side
    Point(2) = Point(2)   + (fTPCZLengthC+Point(2))/(vdCorr*fTPCvd+vdY*Point(1))
                          * (fTPCvd*(1.-vdCorr)-vdY*Point(1)) - fTPCvd*t0;
    Point2(2) = Point2(2) + (fTPCZLengthC+Point2(2))/(vdCorr*fTPCvd+vdY*Point2(1))
                          * (fTPCvd*(1.-vdCorr)-vdY*Point2(1)) - fTPCvd*t0;
  }
  Dir = Point2-Point;
  Dir=Dir.Unit(); //keep unit length

  //Rotation
  TMatrixD rotmat(3,3);
  RotMat( rotmat, misal );
  Point = rotmat * Point;
  Dir = rotmat * Dir;
  //Shift
  Point(0) += misal(3); //add shift in x
  Point(1) += misal(4); //add shift in y
  Point(2) += misal(5); //add shift in z

  //Turn back to local system
  Dir.GetXYZ(dir);
  Point.GetXYZ(point);
  tr->Global2LocalPosition(point,alpha);
  tr->Global2LocalPosition(dir,alpha);

  //Calculate new intersection point with ref plane
  Double_t p[5],pcov[15];
  if (dir[0]==0.0) return kFALSE;
  Double_t s=(x-point[0])/dir[0];
  p[0] = point[1]+s*dir[1];
  p[1] = point[2]+s*dir[2];
  Double_t pt = TMath::Sqrt(dir[0]*dir[0]+dir[1]*dir[1]);
  if (pt==0.0) return kFALSE;
  p[2] = dir[1]/pt;
  p[3] = dir[2]/pt;
  //insert everything back into track
  const Double_t* pcovtmp = tr->GetCovariance();
  p[4] = tr->GetSigned1Pt(); //copy the momentum
  memcpy(pcov,pcovtmp,15*sizeof(Double_t));
  tr->Set(x,alpha,p,pcov);
  return kTRUE;

  ////put params back into track and propagate to ref
  //Double_t p[5];
  //Double_t pcov[15];
  //p[0] = point[1];
  //p[1] = point[2];
  //Double_t xnew = point[0];
  //Double_t pt = TMath::Sqrt(dir[0]*dir[0]+dir[1]*dir[1]);
  //if (pt==0.0) return kFALSE;
  //p[2] = dir[1]/pt;
  //p[3] = dir[2]/pt;
  //p[4] = tr->GetSigned1Pt(); //copy the momentum
  //const Double_t* pcovtmp = tr->GetCovariance();
  //memcpy(pcov,pcovtmp,15*sizeof(Double_t));
  //printf("x before: %.5f, after: %.5f\n",x, xnew);
  //printf("before: %.4f %.4f %.4f %.4f %.4f \n",tr->GetParameter()[0],tr->GetParameter()[1],tr->GetParameter()[2],tr->GetParameter()[3],tr->GetParameter()[4]);
  //printf("after:  %.4f %.4f %.4f %.4f %.4f \n",p[0],p[1],p[2],p[3],p[4]);
  //tr->Set(xnew,alpha,p,pcov);
  //return tr->PropagateTo(x,fMagField);
}

//______________________________________________________________________________
void AliRelAlignerKalman::Reset()
{
  //full reset to defaults
  fPX->Zero();
  (*fPX)(6)=1.;
  ResetCovariance();

  //initialize the differentials per parameter
  for (Int_t i=0;i<4;i++) 
  {
    delete [] (fResArrSigma2Median[i]);
  }
  fRejectOutliersSigma2Median=kFALSE;

  fNMatchedCosmics=0;
  fNMatchedTPCtracklets=0;
  fNUpdates=0;
  fNOutliers=0;
  fNTracks=0;
  fNProcessedEvents=0;
  fRunNumber=0;
  fTimeStamp=0;
}

//______________________________________________________________________________
void AliRelAlignerKalman::ResetCovariance( const Double_t number )
{
  //Resets the covariance to the default if arg=0 or resets the off diagonals
  //to zero and releases the diagonals by factor arg.
  if (number!=0.)
  {
    for (Int_t z=0;z<6;z++)
    {
      for (Int_t zz=0;zz<6;zz++)
      {
        if (zz==z) continue; //don't touch diagonals
        (*fPXcov)(zz,z) = 0.;
        (*fPXcov)(z,zz) = 0.;
      }
      (*fPXcov)(z,z) = (*fPXcov)(z,z) * number;
    }
  }
  else
  {
    //Resets the covariance of the fit to a default value
    fPXcov->Zero();
    (*fPXcov)(0,0) = .08*.08; //psi (rad)
    (*fPXcov)(1,1) = .08*.08; //theta (rad
    (*fPXcov)(2,2) = .08*.08; //phi (rad)
    (*fPXcov)(3,3) = .3*.3; //x (cm)
    (*fPXcov)(4,4) = .3*.3; //y (cm)
    (*fPXcov)(5,5) = .3*.3; //z (cm)
  }
  ResetTPCparamsCovariance(number); 
}

//______________________________________________________________________________
void AliRelAlignerKalman::ResetTPCparamsCovariance( const Double_t number )
{
  //Resets the covariance to the default if arg=0 or resets the off diagonals
  //to zero and releases the diagonals by factor arg.
  
  //release diagonals
  if (number==0.)
  {
    if (fgkNSystemParams>6) (*fPXcov)(6,6) = .1*.1;
    if (fgkNSystemParams>7) (*fPXcov)(7,7) = 1.*1.;
    if (fgkNSystemParams>8) (*fPXcov)(8,8) = .1*.1;
  }
  else
  {
    if (fgkNSystemParams>6) (*fPXcov)(6,6) = number * (*fPXcov)(6,6);
    if (fgkNSystemParams>7) (*fPXcov)(7,7) = number * (*fPXcov)(7,7);
    if (fgkNSystemParams>8) (*fPXcov)(8,8) = number * (*fPXcov)(8,8);
  }
  
  //set crossterms to zero
  for (Int_t i=0;i<fgkNSystemParams;i++)
  {
    for (Int_t j=6;j<fgkNSystemParams;j++) //TPC params
    {
      if (i==j) continue; //don't touch diagonals
      (*fPXcov)(i,j) = 0.;
      (*fPXcov)(j,i) = 0.;
    }
  }
}

//______________________________________________________________________________
Bool_t AliRelAlignerKalman::Merge( const AliRelAlignerKalman* al )
{
  //Merge two aligners
  
  if (!al) return kFALSE;
  if (al==this) return kTRUE;
  if (al->fNUpdates == 0) return kTRUE; //no point in merging with an empty one
  
  //store the pointers to current stuff
  TVectorD* pmes = fPMeasurement;
  TMatrixDSym* pmescov = fPMeasurementCov;
  TVectorD* pmespred = fPMeasurementPrediction;
  TMatrixD* ph = fPH;

  //make a unity system matrix
  TMatrixD tmp(fgkNSystemParams,fgkNSystemParams);
  fPH = new TMatrixD(TMatrixD::kUnit, tmp);

  //mesurement is the state of the new aligner
  fPMeasurement = al->fPX;
  fPMeasurementCov = al->fPXcov;

  //the mesurement prediction is the state
  fPMeasurementPrediction = fPX; //this is safe as fPX doesn't change until end
  
  //do the merging
  Bool_t success = UpdateEstimateKalman();
  
  //restore pointers to old stuff
  fPMeasurement = pmes;
  fPMeasurementCov = pmescov;
  fPMeasurementPrediction = pmespred;
  delete fPH;
  fPH = ph;

  //merge stats
  if (!success)    
  {
    fNMergesFailed++;
    //printf("AliRelAlignerKalman::Merge failed\n");
    return kFALSE; //no point in merging stats if merge not succesful
  }
  fNProcessedEvents += al->fNProcessedEvents;
  fNUpdates += al->fNUpdates;
  fNOutliers += al->fNOutliers;
  fNOutliersSigma2Median += al->fNOutliersSigma2Median;
  fNTracks += al->fNTracks;
  fNMatchedTPCtracklets += al->fNMatchedTPCtracklets;
  fNMatchedCosmics += al->fNMatchedCosmics;
  if (fNMerges==0 || al->fNMerges==0) fNMerges++;
  else fNMerges += al->fNMerges;
  if (fTimeStamp < al->fTimeStamp) fTimeStamp = al->fTimeStamp; //take the newer one

  return success;
}

//______________________________________________________________________________
Long64_t AliRelAlignerKalman::Merge( TCollection* list )
{
  //merge all aligners in the collection
  Long64_t numberOfMerges=0;
  AliRelAlignerKalman* alignerFromList;
  if (!list) return 0;
  TIter next(list);
  while ( (alignerFromList = dynamic_cast<AliRelAlignerKalman*>(next())) )
  {
    if (alignerFromList == this) continue;
    if (Merge(alignerFromList)) numberOfMerges++;
  }
  return numberOfMerges;
}

//______________________________________________________________________________
Int_t AliRelAlignerKalman::Compare(const TObject *obj) const
{
  if (this == obj) return 0;
  const AliRelAlignerKalman* aobj = dynamic_cast<const AliRelAlignerKalman*>(obj);
  if (!aobj) return 0;
  if (fTimeStamp < aobj->fTimeStamp) return -1;
  else if (fTimeStamp > aobj->fTimeStamp) return 1;
  else return 0;
}

//________________________________________________________________________
Int_t AliRelAlignerKalman::FindMatchingTracks(TObjArray& arrITS, TObjArray& arrTPC, AliESDEvent* pESD)
{
  //find matching tracks and return tobjarrays with the params
  
  Int_t ntracks = pESD->GetNumberOfTracks();
  Int_t* matchedArr = new Int_t[ntracks]; //storage for index of ITS track for which a match was found
  for (Int_t i=0;i<ntracks;i++)
  {
    matchedArr[i]=0;
  }

  Int_t iMatched=-1;
  for (Int_t i=0; i<ntracks; i++)
  {
    //get track1 and friend
    AliESDtrack* track1 = pESD->GetTrack(i);
    if (!track1) continue;

    if (track1->GetNcls(0) < fMinPointsVol1) continue;

    if (!( ( track1->IsOn(AliESDtrack::kITSrefit)) &&
           (!track1->IsOn(AliESDtrack::kTPCin)) )) continue;

    const AliESDfriendTrack* constfriendtrack1 = track1->GetFriendTrack();
    if (!constfriendtrack1) continue;
    AliESDfriendTrack friendtrack1(*constfriendtrack1);
    
    if (!friendtrack1.GetITSOut()) continue;
    AliExternalTrackParam params1(*(friendtrack1.GetITSOut()));

    Double_t bestd = 1000.; //best distance
    Bool_t newi = kTRUE; //whether we start with a new i
    for (Int_t j=0; j<ntracks; j++)
    {
      if (matchedArr[j]>0 && matchedArr[j]!=i) continue; //already matched, everything tried 
      //get track2 and friend
      AliESDtrack* track2 = pESD->GetTrack(j);
      if (!track2) continue;
      if (track1==track2) continue;
      //if ( ( ( track2->IsOn(AliESDtrack::kITSout)) &&
      //       (!track2->IsOn(AliESDtrack::kTPCin)) )) continue; //all but ITS standalone

      if (track2->GetNcls(0) != track1->GetNcls(0)) continue;
      if (track2->GetITSClusterMap() != track1->GetITSClusterMap()) continue;
      if (track2->GetNcls(1) < fMinPointsVol2) continue; //min 80 clusters in TPC
      if (track2->GetTgl() > 1.) continue; //acceptance
      //cut crossing tracks
      if (track2->GetOuterParam()->GetZ()*track2->GetInnerParam()->GetZ()<0) continue;
      if (track2->GetInnerParam()->GetX()>90) continue;
      if (TMath::Abs(track2->GetInnerParam()->GetZ())<10.) continue; //too close to membrane?


      if (!track2->GetInnerParam()) continue;
      AliExternalTrackParam params2(*(track2->GetInnerParam()));

      //bring to same reference plane
      if (!params2.Rotate(params1.GetAlpha())) continue;
      if (!params2.PropagateTo(params1.GetX(), pESD->GetMagneticField())) continue;

      //pt cut
      if (params2.Pt() < fMinPt) continue;

      const Double32_t*	p1 = params1.GetParameter();
      const Double32_t*	p2 = params2.GetParameter();

      //hard cuts
      Double_t dy = TMath::Abs(p2[0]-p1[0]);
      Double_t dz = TMath::Abs(p2[1]-p1[1]);
      Double_t dphi = TMath::Abs(p2[2]-p1[2]);
      Double_t dlam = TMath::Abs(p2[3]-p1[3]);
      if (dy > 2.0) continue;
      if (dz > 10.0) continue;
      if (dphi > 0.1 ) continue;
      if (dlam > 0.1 ) continue;

      //best match only
      Double_t d = TMath::Sqrt(dy*dy+dz*dz+dphi*dphi+dlam*dlam);
      if ( d >= bestd) continue;
      bestd = d;
      matchedArr[j]=i; //j-th track matches i-th (ITS) track
      if (newi) iMatched++; newi=kFALSE; //increment at most once per i
      if (arrITS[iMatched] && arrTPC[iMatched])
      {
        *(arrITS[iMatched]) = params1;
        *(arrTPC[iMatched]) = params2;
      }
      else
      {
        arrITS[iMatched] = new AliExternalTrackParam(params1);
        arrTPC[iMatched] = new AliExternalTrackParam(params2);
      }//else
    }//for j
  }//for i
  delete [] matchedArr;
  return iMatched;
}

//________________________________________________________________________
void AliRelAlignerKalman::SetRejectOutliersSigma2Median(const Bool_t set )
{
  //Sets up or destroys the memory hungry array to hold the statistics
  //for data rejection with median
  if (set)
  {
    for (Int_t i=0;i<4;i++) 
    {
      if (!fResArrSigma2Median[i]) fResArrSigma2Median[i] = 
                                   new Double_t[fgkNtracksSigma2Median];
    }
    fRejectOutliersSigma2Median = kTRUE;
  }//else
  else
  {
    // it probably doesn't make sense to delete the arrays, they are not streamed
    //if (fRejectOutliersSigma2Median)
    //for (Int_t i=0;i<4;i++) 
    //{
    //  delete [] (fResArrSigma2Median[i]);
    //}
    fRejectOutliersSigma2Median = kFALSE;
  }//if
}
