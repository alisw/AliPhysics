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
#include <TH1D.h>
#include <TF1.h>

#include "AliESDtrack.h"
#include "AliTrackPointArray.h"
#include "AliGeomManager.h"
#include "AliTrackFitterKalman.h"
#include "AliTrackFitterRieman.h"
#include "AliESDfriendTrack.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliExternalTrackParam.h"

#include "AliRelAlignerKalman.h"

ClassImp(AliRelAlignerKalman)

//______________________________________________________________________________
AliRelAlignerKalman::AliRelAlignerKalman():
    TObject(),
    fAlpha(0.),
    fLocalX(80.),
    fPTrackParamArr1(new AliExternalTrackParam[fgkNTracksPerMeasurement]),
    fPTrackParamArr2(new AliExternalTrackParam[fgkNTracksPerMeasurement]),
    fMagField(0.),
    fPX(new TVectorD( fgkNSystemParams )),
    fPXcov(new TMatrixDSym( fgkNSystemParams )),
    fPH(new TMatrixD( fgkNMeasurementParams, fgkNSystemParams )),
    fQ(1.e-15),
    fPMeasurement(new TVectorD( fgkNMeasurementParams )),
    fPMeasurementCov(new TMatrixDSym( fgkNMeasurementParams )),
    fOutRejSigmas(1.),
    fDelta(new Double_t[fgkNSystemParams]),
    fNumericalParanoia(kTRUE),
    fRejectOutliers(kTRUE),
    fRequireMatchInTPC(kFALSE),
    fCuts(kFALSE),
    fMinPointsVol1(3),
    fMinPointsVol2(50),
    fMinMom(0.),
    fMaxMom(1.e100),
    fMaxMatchingAngle(0.1),
    fMaxMatchingDistance(10.),  //in cm
    fCorrectionMode(kFALSE),
    fNTracks(0),
    fNUpdates(0),
    fNOutliers(0),
    fNMatchedCosmics(0),
    fNMatchedTPCtracklets(0),
    fNProcessedEvents(0),
    fTrackInBuffer(0),
    fTPCvd(2.64),
    fTPCZLengthA(2.4972500e02),
    fTPCZLengthC(2.4969799e02)
{
  //Default constructor

  //default seed: zero, reset errors to large default
  Reset();
}

//______________________________________________________________________________
AliRelAlignerKalman::AliRelAlignerKalman(const AliRelAlignerKalman& a):
    TObject(static_cast<TObject>(a)),
    fAlpha(a.fAlpha),
    fLocalX(a.fLocalX),
    fPTrackParamArr1(new AliExternalTrackParam[fgkNTracksPerMeasurement]),
    fPTrackParamArr2(new AliExternalTrackParam[fgkNTracksPerMeasurement]),
    fMagField(a.fMagField),
    fPX(new TVectorD( *a.fPX )),
    fPXcov(new TMatrixDSym( *a.fPXcov )),
    fPH(new TMatrixD( *a.fPH )),
    fQ(a.fQ),
    fPMeasurement(new TVectorD( *a.fPMeasurement )),
    fPMeasurementCov(new TMatrixDSym( *a.fPMeasurementCov )),
    fOutRejSigmas(a.fOutRejSigmas),
    fDelta(new Double_t[fgkNSystemParams]),
    fNumericalParanoia(a.fNumericalParanoia),
    fRejectOutliers(a.fRejectOutliers),
    fRequireMatchInTPC(a.fRequireMatchInTPC),
    fCuts(a.fCuts),
    fMinPointsVol1(a.fMinPointsVol1),
    fMinPointsVol2(a.fMinPointsVol2),
    fMinMom(a.fMinMom),
    fMaxMom(a.fMaxMom),
    fMaxMatchingAngle(a.fMaxMatchingAngle),
    fMaxMatchingDistance(a.fMaxMatchingDistance),  //in cm
    fCorrectionMode(a.fCorrectionMode),
    fNTracks(a.fNTracks),
    fNUpdates(a.fNUpdates),
    fNOutliers(a.fNOutliers),
    fNMatchedCosmics(a.fNMatchedCosmics),
    fNMatchedTPCtracklets(a.fNMatchedTPCtracklets),
    fNProcessedEvents(a.fNProcessedEvents),
    fTrackInBuffer(a.fTrackInBuffer),
    fTPCvd(a.fTPCvd),
    fTPCZLengthA(a.fTPCZLengthA),
    fTPCZLengthC(a.fTPCZLengthC)
    //fApplyCovarianceCorrection(a.fApplyCovarianceCorrection),
    //fCalibrationMode(a.fCalibrationMode),
    //fFillHistograms(a.fFillHistograms),
    //fNHistogramBins(a.fNHistogramBins),
    //fPMes0Hist(new TH1D(*a.fPMes0Hist)),
    //fPMes1Hist(new TH1D(*a.fPMes1Hist)),
    //fPMes2Hist(new TH1D(*a.fPMes2Hist)),
    //fPMes3Hist(new TH1D(*a.fPMes3Hist)),
    //fPMesErr0Hist(new TH1D(*a.fPMesErr0Hist)),
    //fPMesErr1Hist(new TH1D(*a.fPMesErr1Hist)),
    //fPMesErr2Hist(new TH1D(*a.fPMesErr2Hist)),
    //fPMesErr3Hist(new TH1D(*a.fPMesErr3Hist)),
    //fPMeasurementCovCorr(new TMatrixDSym(*a.fPMeasurementCovCorr)),
{
  //copy constructor
  memcpy(fDelta,a.fDelta,fgkNSystemParams*sizeof(Double_t));
}

//______________________________________________________________________________
AliRelAlignerKalman& AliRelAlignerKalman::operator=(const AliRelAlignerKalman& a)
{
  //assignment operator
  //fAlpha=a.fAlpha;
  //fLocalX=a.fLocalX;
  //memcpy(fPTrackParamArr1,a.fPTrackParamArr1,fgkNTracksPerMeasurement*sizeof(AliExternalTrackParam));
  //memcpy(fPTrackParamArr2,a.fPTrackParamArr2,fgkNTracksPerMeasurement*sizeof(AliExternalTrackParam));
  fMagField=a.fMagField,
  *fPX = *a.fPX;
  *fPXcov = *a.fPXcov;
  //*fPH = *a.fPH;
  fQ=a.fQ;
  //*fPMeasurement=*a.fPMeasurement;
  //*fPMeasurementCov=*a.fPMeasurementCov;
  fOutRejSigmas=a.fOutRejSigmas;
  memcpy(fDelta,a.fDelta,fgkNSystemParams*sizeof(Double_t));
  fNumericalParanoia=a.fNumericalParanoia;
  fRejectOutliers=a.fRejectOutliers;
  fRequireMatchInTPC=a.fRequireMatchInTPC;
  fCuts=a.fCuts;
  fMinPointsVol1=a.fMinPointsVol1;
  fMinPointsVol2=a.fMinPointsVol2;
  fMinMom=a.fMinMom;
  fMaxMom=a.fMaxMom;
  fMaxMatchingAngle=a.fMaxMatchingAngle;
  fMaxMatchingDistance=a.fMaxMatchingDistance;  //in c;
  fCorrectionMode=a.fCorrectionMode;
  fNTracks=a.fNTracks;
  fNUpdates=a.fNUpdates;
  fNOutliers=a.fNOutliers;
  fNMatchedCosmics=a.fNMatchedCosmics;
  fNMatchedTPCtracklets=a.fNMatchedTPCtracklets;
  fNProcessedEvents=a.fNProcessedEvents;
  fTrackInBuffer=a.fTrackInBuffer;
  //fApplyCovarianceCorrection=a.fApplyCovarianceCorrection;
  //fCalibrationMode=a.fCalibrationMode;
  //fFillHistograms=a.fFillHistograms;
  //fNHistogramBins=a.fNHistogramBins;
  //*fPMes0Hist=*(a.fPMes0Hist);
  //*fPMes1Hist=*(a.fPMes1Hist);
  //*fPMes2Hist=*(a.fPMes2Hist);
  //*fPMes3Hist=*(a.fPMes3Hist);
  //*fPMesErr0Hist=*(a.fPMesErr0Hist);
  //*fPMesErr1Hist=*(a.fPMesErr1Hist);
  //*fPMesErr2Hist=*(a.fPMesErr2Hist);
  //*fPMesErr3Hist=*(a.fPMesErr3Hist);
  //*fPMeasurementCovCorr=*(a.fPMeasurementCovCorr);
  fTPCvd=a.fTPCvd;
  fTPCZLengthA=a.fTPCZLengthA;
  fTPCZLengthC=a.fTPCZLengthC;
  return *this;
}

//______________________________________________________________________________
AliRelAlignerKalman::~AliRelAlignerKalman()
{
  //destructor
  if (fPTrackParamArr1) delete [] fPTrackParamArr1;
  if (fPTrackParamArr2) delete [] fPTrackParamArr2;
  if (fPX) delete fPX;
  if (fPXcov) delete fPXcov;
  if (fPH) delete fPH;
  if (fPMeasurement) delete fPMeasurement;
  if (fPMeasurementCov) delete fPMeasurementCov;
  if (fDelta) delete [] fDelta;
  //if (fPMes0Hist) delete fPMes0Hist;
  //if (fPMes1Hist) delete fPMes1Hist;
  //if (fPMes2Hist) delete fPMes2Hist;
  //if (fPMes3Hist) delete fPMes3Hist;
  //if (fPMesErr0Hist) delete fPMesErr0Hist;
  //if (fPMesErr1Hist) delete fPMesErr1Hist;
  //if (fPMesErr2Hist) delete fPMesErr2Hist;
  //if (fPMesErr3Hist) delete fPMesErr3Hist;
  //if (fPMeasurementCovCorr) delete fPMeasurementCovCorr;
}

//______________________________________________________________________________
Bool_t AliRelAlignerKalman::AddESDevent( const AliESDEvent* pEvent )
{
  //Add all tracks in an ESD event

  fNProcessedEvents++; //update the counter
  
  Bool_t success=kFALSE;
  SetMagField( pEvent->GetMagneticField() );
  AliESDtrack* track;
  
  for (Int_t i=0; i<pEvent->GetNumberOfTracks(); i++)
  {
    track = pEvent->GetTrack(i);
    if (!track) continue;
    if ( ((track->GetStatus()&AliESDtrack::kTPCin)>0)&&
         ((track->GetStatus()&AliESDtrack::kITSout)>0)&&
         (track->GetNcls(0)>=fMinPointsVol1)&&
         (track->GetNcls(1)>=fMinPointsVol2) )
    { 
      success = ( AddESDtrack( track ) || success );
    }
  }
  return success;
}

//______________________________________________________________________________
Bool_t AliRelAlignerKalman::AddESDtrack( const AliESDtrack* pTrack )
{
  ////Adds a full track, returns true if results in a new estimate

  //const AliExternalTrackParam* pconstparamsITS = pTrack->GetOuterParam();
  //if (!pconstparamsITS) return kFALSE;
  //const AliExternalTrackParam* pconstparamsTPC = pTrack->GetInnerParam();
  //if (!pconstparamsTPC) return kFALSE;
  //
  ////TPC part
  //AliExternalTrackParam paramsTPC = (*pconstparamsTPC);
  //paramsTPC.Rotate(pconstparamsITS->GetAlpha());
  //paramsTPC.PropagateTo(pconstparamsITS->GetX(), pconstparamsITS->GetAlpha());

  //if (!SetTrackParams( pconstparamsITS, &paramsTPC )) return kFALSE;
  //
  ////do some accounting and update
  //return (Update());

  const AliESDtrack* p;
  p=pTrack; //shuts up the compiler
  
  return kFALSE;
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
  AliESDtrack* ptrack;
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
    params2.PropagateTo( params1.GetX(), params1.GetAlpha() );

    if (!SetTrackParams( &params1, &params2 )) continue;
    
    //do some accounting and update
    if (Update())
      success = kTRUE;
    else
      continue;
  }
  return success;
}

//______________________________________________________________________________
void AliRelAlignerKalman::Print(Option_t*) const
{
  //Print some useful info
  Double_t rad2deg = 180./TMath::Pi();
  printf("\nAliRelAlignerKalman\n");
  if (fCorrectionMode) printf("(Correction mode)\n");
  printf("  %i pairs, %i updates, %i outliers,\n", fNTracks, fNUpdates, fNOutliers );
  printf("  %i TPC matches, %i ITS-TPC matches in %i events\n", fNMatchedTPCtracklets, fNMatchedCosmics, fNProcessedEvents );
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
  for ( Int_t i=0; i<fgkNMeasurementParams; i++ )
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
  
  fPTrackParamArr1[fTrackInBuffer] = *exparam1;
  fPTrackParamArr2[fTrackInBuffer] = *exparam2;
  
  fTrackInBuffer++;

  if (fTrackInBuffer == fgkNTracksPerMeasurement)
  {
    fTrackInBuffer = 0;
    return kTRUE;
  }
  return kFALSE;
}

//______________________________________________________________________________
void AliRelAlignerKalman::SetRefSurface( const Double_t radius, const Double_t alpha )
{
  //sets the reference surface by setting the radius (localx)
  //and rotation angle wrt the global frame of reference
  //locally the reference surface becomes a plane with x=r;
  fLocalX = radius;
  fAlpha = alpha;
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
  
  for (Int_t i=0;i<fgkNTracksPerMeasurement;i++)
  {
    const Double_t* pararr1 = fPTrackParamArr1[i].GetParameter();
    const Double_t* pararr2 = fPTrackParamArr2[i].GetParameter();

    //Take the track parameters and calculate the input to the Kalman filter
    Int_t x = i*4;
    (*fPMeasurement)(x+0) = pararr2[0]-pararr1[0];
    (*fPMeasurement)(x+1) = pararr2[1]-pararr1[1];
    (*fPMeasurement)(x+2) = pararr2[2]-pararr1[2];
    (*fPMeasurement)(x+3) = pararr2[3]-pararr1[3];

    //the covariance
    const Double_t* parcovarr1 = fPTrackParamArr1[i].GetCovariance();
    const Double_t* parcovarr2 = fPTrackParamArr2[i].GetCovariance();
    (*fPMeasurementCov)(x+0,x+0)=parcovarr1[0];
    (*fPMeasurementCov)(x+0,x+1)=parcovarr1[1];
    (*fPMeasurementCov)(x+0,x+2)=parcovarr1[3];
    (*fPMeasurementCov)(x+0,x+3)=parcovarr1[6];
    (*fPMeasurementCov)(x+1,x+0)=parcovarr1[1];
    (*fPMeasurementCov)(x+1,x+1)=parcovarr1[2];
    (*fPMeasurementCov)(x+1,x+2)=parcovarr1[4];
    (*fPMeasurementCov)(x+1,x+3)=parcovarr1[7];
    (*fPMeasurementCov)(x+2,x+0)=parcovarr1[3];
    (*fPMeasurementCov)(x+2,x+1)=parcovarr1[4];
    (*fPMeasurementCov)(x+2,x+2)=parcovarr1[5];
    (*fPMeasurementCov)(x+2,x+3)=parcovarr1[8];
    (*fPMeasurementCov)(x+3,x+0)=parcovarr1[6];
    (*fPMeasurementCov)(x+3,x+1)=parcovarr1[7];
    (*fPMeasurementCov)(x+3,x+2)=parcovarr1[8];
    (*fPMeasurementCov)(x+3,x+3)=parcovarr1[9];
    (*fPMeasurementCov)(x+0,x+0)+=parcovarr2[0];
    (*fPMeasurementCov)(x+0,x+1)+=parcovarr2[1];
    (*fPMeasurementCov)(x+0,x+2)+=parcovarr2[3];
    (*fPMeasurementCov)(x+0,x+3)+=parcovarr2[6];
    (*fPMeasurementCov)(x+1,x+0)+=parcovarr2[1];
    (*fPMeasurementCov)(x+1,x+1)+=parcovarr2[2];
    (*fPMeasurementCov)(x+1,x+2)+=parcovarr2[4];
    (*fPMeasurementCov)(x+1,x+3)+=parcovarr2[7];
    (*fPMeasurementCov)(x+2,x+0)+=parcovarr2[3];
    (*fPMeasurementCov)(x+2,x+1)+=parcovarr2[4];
    (*fPMeasurementCov)(x+2,x+2)+=parcovarr2[5];
    (*fPMeasurementCov)(x+2,x+3)+=parcovarr2[8];
    (*fPMeasurementCov)(x+3,x+0)+=parcovarr2[6];
    (*fPMeasurementCov)(x+3,x+1)+=parcovarr2[7];
    (*fPMeasurementCov)(x+3,x+2)+=parcovarr2[8];
    (*fPMeasurementCov)(x+3,x+3)+=parcovarr2[9];
    
    fNTracks++; //count added track sets
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

  TVectorD z1( fgkNMeasurementParams );
  TVectorD z2( fgkNMeasurementParams );
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
    for (Int_t j=0; j<fgkNMeasurementParams; j++ )
    {
      (*fPH)(j,i) = ( z2(j)-z1(j) ) / fDelta[i];
    }
  }
  return kTRUE;
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
    for (Int_t i=0;i<fgkNTracksPerMeasurement;i++)
    {
      AliExternalTrackParam track(fPTrackParamArr2[i]); //make a copy track
      if (!CorrectTrack( &track, state )) return kFALSE; //predict what the ideal track would be by applying correction
      
      const Double_t* oldparam = fPTrackParamArr2[i].GetParameter();
      const Double_t* newparam = track.GetParameter();

      Int_t x = 4*i;
      //calculate the predicted residual
      pred(x+0) = oldparam[0] - newparam[0];
      pred(x+1) = oldparam[1] - newparam[1];
      pred(x+2) = oldparam[2] - newparam[2];
      pred(x+3) = oldparam[3] - newparam[3];
      return kTRUE;
    }
  }
  else
  {
    for (Int_t i=0;i<fgkNTracksPerMeasurement;i++)
    {
      AliExternalTrackParam track(fPTrackParamArr1[i]); //make a copy track
      if (!MisalignTrack( &track, state )) return kFALSE; //predict what the measured track would be by applying misalignment

      const Double_t* oldparam = fPTrackParamArr1[i].GetParameter();
      const Double_t* newparam = track.GetParameter();

      Int_t x = 4*i;
      //calculate the predicted residual
      pred(x+0) = newparam[0] - oldparam[0];
      pred(x+1) = newparam[1] - oldparam[1];
      pred(x+2) = newparam[2] - oldparam[2];
      pred(x+3) = newparam[3] - oldparam[3];
      return kTRUE;
    }
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
    Double_t det;
    hpht.InvertFast(&det); //since the matrix is small...
    if (det < 2e-55) return kFALSE; //we need some sort of protection even in this case....
//  }
  //printf("KalmanUpdate: det(hpht): %.4g\n",det);

  TMatrixD k(pht, TMatrixD::kMult, hpht ); //compute K (hpht is already inverted)

  // update the state and its covariance matrix
  TVectorD xupdate(fgkNSystemParams);
  TVectorD hx(fgkNMeasurementParams);
  PredictMeasurement( hx, (*fPX) );
  xupdate = k*((*fPMeasurement)-hx);

  //SIMPLE OUTLIER REJECTION
  if ( IsOutlier( xupdate, (*fPXcov) ) && fRejectOutliers )
  {
    fNOutliers++;
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
Bool_t AliRelAlignerKalman::IsPositiveDefinite( const TMatrixD& mat ) const
{
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
      if (pTrack->GetP()<fMinMom || pTrack->GetP()>fMaxMom) continue;
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
  for (Int_t i=0;i<fgkNSystemParams;i++) fDelta[i] = 1.e-6;

  fNMatchedCosmics=0;
  fNMatchedTPCtracklets=0;
  fNUpdates=0;
  fNOutliers=0;
  fNTracks=0;
  fNProcessedEvents=0;
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
    if (fgkNSystemParams>8) (*fPXcov)(8,8) = 1.*1.;
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
//void AliRelAlignerKalman::PrintCovarianceCorrection()
//{
//  //Print the measurement covariance correction matrix
//  printf("Covariance correction matrix:\n");
//  for ( Int_t i=0; i<fgkNMeasurementParams; i++ )
//  {
//    for ( Int_t j=0; j<i+1; j++ )
//    {
//      printf("% -2.2f  ", (*fPMeasurementCovCorr)(i,j) );
//    }//for i
//    printf("\n");
//  }//for j
//  printf("\n");
//  return;
//}

//_______________________________________________________________________________
//Bool_t AliRelAlignerKalman::UpdateCalibration()
//{
//  //Update the calibration with new data (in calibration mode)
//
//  fPMes0Hist->Fill( (*fPMeasurement)(0) );
//  fPMes1Hist->Fill( (*fPMeasurement)(1) );
//  fPMes2Hist->Fill( (*fPMeasurement)(2) );
//  fPMes3Hist->Fill( (*fPMeasurement)(3) );
//  fPMesErr0Hist->Fill( TMath::Sqrt((*fPMeasurementCov)(0,0)) );
//  fPMesErr1Hist->Fill( TMath::Sqrt((*fPMeasurementCov)(1,1)) );
//  fPMesErr2Hist->Fill( TMath::Sqrt((*fPMeasurementCov)(2,2)) );
//  fPMesErr3Hist->Fill( TMath::Sqrt((*fPMeasurementCov)(3,3)) );
//  return kTRUE;
//}

//______________________________________________________________________________
//Bool_t AliRelAlignerKalman::SetCalibrationMode( const Bool_t cp )
//{
//  //sets the calibration mode
//  if (cp)
//  {
//    fCalibrationMode=kTRUE;
//    return kTRUE;
//  }//if (cp)
//  else
//  {
//    if (fCalibrationMode) // do it only after the calibration pass
//    {
//      CalculateCovarianceCorrection();
//      SetApplyCovarianceCorrection();
//      fCalibrationMode=kFALSE;
//      return kTRUE;
//    }//if (fCalibrationMode)
//  }//else (cp)
//  return kFALSE;
//}

//______________________________________________________________________________
//Bool_t AliRelAlignerKalman::CalculateCovarianceCorrection()
//{
//  //Calculates the correction to the measurement covariance
//  //using the calibration histograms
//
//  fPMeasurementCovCorr->Zero(); //reset the correction
//
//  Double_t s,m,c;  //sigma,meansigma,correction
//
//  //TF1* fitformula;
//  //fPMes0Hist->Fit("gaus");
//  //fitformula = fPMes0Hist->GetFunction("gaus");
//  //s = fitformula->GetParameter(2);   //spread of the measurement
//  //fPMesErr0Hist->Fit("gaus");
//  //fitformula = fPMesErr0Hist->GetFunction("gaus"); //average error from cov matrices
//  //m = fitformula->GetParameter(1);
//  s = fPMes0Hist->GetRMS();
//  m = fPMesErr0Hist->GetMean();
//  c = s-m; //the difference between the average error and real spread of the data
//  if (c>0) //only correct is spread bigger than average error
//    (*fPMeasurementCovCorr)(0,0) = c*c;
//
//  //fPMes1Hist->Fit("gaus");
//  //fitformula = fPMes1Hist->GetFunction("gaus");
//  //s = fitformula->GetParameter(2);
//  //fPMesErr1Hist->Fit("gaus");
//  //fitformula = fPMesErr1Hist->GetFunction("gaus");
//  //m = fitformula->GetParameter(1);
//  s = fPMes1Hist->GetRMS();
//  m = fPMesErr1Hist->GetMean();
//  c = s-m;
//  if (c>0) //only correct is spread bigger than average error
//    (*fPMeasurementCovCorr)(1,1) = c*c;
//
//  //fPMes2Hist->Fit("gaus");
//  //fitformula = fPMes2Hist->GetFunction("gaus");
//  //s = fitformula->GetParameter(2);
//  //fPMesErr2Hist->Fit("gaus");
//  //fitformula = fPMesErr2Hist->GetFunction("gaus");
//  //m = fitformula->GetParameter(1);
//  s = fPMes2Hist->GetRMS();
//  m = fPMesErr2Hist->GetMean();
//  c = s-m;
//  if (c>0) //only correct is spread bigger than average error
//    (*fPMeasurementCovCorr)(2,2) = c*c;
//
//  //fPMes3Hist->Fit("gaus");
//  //fitformula = fPMes3Hist->GetFunction("gaus");
//  //s = fitformula->GetParameter(2);
//  //fPMesErr3Hist->Fit("gaus");
//  //fitformula = fPMesErr3Hist->GetFunction("gaus");
//  //m = fitformula->GetParameter(1);
//  s = fPMes3Hist->GetRMS();
//  m = fPMesErr3Hist->GetMean();
//  c = s-m;
//  if (c>0) //only correct is spread bigger than average error
//    (*fPMeasurementCovCorr)(3,3) = c*c;
//
//  return kTRUE;
//}

