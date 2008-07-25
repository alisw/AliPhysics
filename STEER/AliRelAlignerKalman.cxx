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
//    Determines the transformation of the second volume (TPC) with respect to
//    the first (ITS) by measuring the residual between the 2 tracks.
//    Additionally calculates some callibration parameters for TPC
//    Fit parameters are:
//    - 3 shifts, x,y,z
//    - 3 Cardan angles, psi, theta, phi (see definition in alignment docs),
//    - TPC drift velocity slope correction (vel(y) = vel(y0) + (1+correction)*y),
//    - TPC offset correction.
//
//    Basic usage:
//    When aligning two volumes at any given time a single instance of
//    the class should be active. The fit of the parameters is updated
//    by adding new data using one of the Add.... methods.
//
//    User methods:
//    In collision events add an ESD track to update the fit,
//    
//        Bool_t AddESDTrack( AliESDtrack* pTrack );
//    
//    For cosmic data, the assumption is that the tracking is done twice:
//    once global and once only ITS and the tracklets are saved inside 
//    one AliESDEvent. The method
//    
//        Bool_t AddCosmicEventSeparateTracking( AliESDEvent* pEvent );
//    
//    then searches the event for matching tracklets and upon succes it updates.
//    One cosmic ideally triggers two updates: for the upper and lower half of
//    the cosmic (upper ITS tracklet+upper TPC tracklet, idem dito for lower)
//    
//    _________________________________________________________________________
//    Expert options:
//    look at AddCosmicEventSeparateTracking(); to get the idea of how the
//    aligner works.
//    
//    The following is dangerous!! Cripples the outlier rejection!
//    In the calibration mode set by
//
//      void SetCalibrationMode( const Bool_t cal=kTRUE );
//
//    a correction for the covariance matrix of the measurement can be calculated
//    in case the errors estimated by the track fit do not correspond to the
//    actual spread of the residuals.
//    In calibration mode the aligner fills histograms of the residuals and of 
//    the errors of the residuals.
//    Setting the calibration mode to false:
//      void SetCalibrationMode( const Bool_t cal=kFALSE );
//    automatically enables the correction.
//
//    Pointers to the histograms are available with apropriate getters for
//    plotting and analysis.
//    
//
//    Origin: Mikolaj Krzewicki, Nikhef, Mikolaj.Krzewicki@cern.ch
//
//////////////////////////////////////////////////////////////////////////////

#include "AliRelAlignerKalman.h"

ClassImp(AliRelAlignerKalman)

//______________________________________________________________________________
AliRelAlignerKalman::AliRelAlignerKalman():
    fAlpha(0.),
    fLocalX(80.),
    fPTrackParam1(NULL),
    fPTrackParam2(NULL),
    fPX(new TVectorD( fgkNSystemParams )),
    fPXcov(new TMatrixDSym( fgkNSystemParams )),
    fPH(new TMatrixD( fgkNMeasurementParams, fgkNSystemParams )),
    fQ(1e-10),
    fPMeasurement(new TVectorD( fgkNMeasurementParams )),
    fPMeasurementCov(new TMatrixDSym( fgkNMeasurementParams )),
    fOutRejSigmas(1.),
    fRejectOutliers(kTRUE),
    fCalibrationMode(kFALSE),
    fFillHistograms(kTRUE),
    fRequireDoubleTPCtrack(kFALSE),
    fApplyCovarianceCorrection(kFALSE),
    fCuts(kFALSE),
    fMinPointsVol1(2),
    fMinPointsVol2(50),
    fMinMom(0.),
    fMaxMom(1e100),
    fMinAbsSinPhi(0.),
    fMaxAbsSinPhi(1.),
    fMinSinTheta(-1.),
    fMaxSinTheta(1.),
    fMaxMatchingAngle(0.2),
    fMaxMatchingDistance(2.),  //in cm
    fNTracks(0),
    fNUpdates(0),
    fNOutliers(0),
    fNMatchedCosmics(0),
    fNProcessedEvents(0),
    fPMes0Hist(new TH1D("y","y", 50, 0, 0)),
    fPMes1Hist(new TH1D("z","z", 50, 0, 0)),
    fPMes2Hist(new TH1D("sinphi","sinphi", 50, 0, 0)),
    fPMes3Hist(new TH1D("tanlambda","tanlambda", 50, 0, 0)),
    fPMesErr0Hist(new TH1D("mescov11","mescov11", 50, 0, 0)),
    fPMesErr1Hist(new TH1D("mescov22","mescov22", 50, 0, 0)),
    fPMesErr2Hist(new TH1D("mescov33","mescov33", 50, 0, 0)),
    fPMesErr3Hist(new TH1D("mescov44","mescov44", 50, 0, 0)),
    fPMeasurementCovCorr(new TMatrixDSym(fgkNMeasurementParams))
{
    //Default constructor
    
    //default seed: zero, reset errors to large default
    ResetCovariance();
    //initialize the differentials per parameter
    for (Int_t i=0;i<fgkNSystemParams;i++) fDelta[i] = 1.e-6;
    //fDelta[0] = 3e-8;
    //fDelta[1] = 3e-8;
    //fDelta[2] = 3e-8;
    //fDelta[3] = 3e-6;
    //fDelta[4] = 3e-6;
    //fDelta[5] = 3e-6;
    //fDelta[6] = 3e-12;
    //fDelta[7] = 3e-8;
    (*fPX)(6)=1.;
}

//______________________________________________________________________________
AliRelAlignerKalman::~AliRelAlignerKalman()
{
    //destructor
    delete fPX;
    delete fPXcov;
    delete fPH;
    delete fPMeasurement;
    delete fPMeasurementCov;
    delete fPMes0Hist;
    delete fPMes1Hist;
    delete fPMes2Hist;
    delete fPMes3Hist;
    delete fPMesErr0Hist;
    delete fPMesErr1Hist;
    delete fPMesErr2Hist;
    delete fPMesErr3Hist;
    delete fDelta;
}

//______________________________________________________________________________
Bool_t AliRelAlignerKalman::AddESDTrack( AliESDtrack* pTrack )
{
    //Adds a full track, to be implemented when data format clear
    if (pTrack) return kFALSE;
    return kFALSE;
}

//______________________________________________________________________________
Bool_t AliRelAlignerKalman::AddCosmicEventSeparateTracking( AliESDEvent* pEvent )
{
    //Add an cosmic with separately tracked ITS and TPC parts, do trackmatching

    fNProcessedEvents++; //update the counter

    Int_t iits1,iits2,itpc1,itpc2;
    if (!FindCosmicTrackletNumbersInEvent( iits1, itpc1, iits2, itpc2, pEvent )) return kFALSE;
   // printf("Found tracks: %d, %d,    %d, %d\n",iits1,itpc1,iits2,itpc2);
    Double_t field = pEvent->GetMagneticField();
    AliESDtrack* ptrack;
    const AliExternalTrackParam* constpparams;
    AliExternalTrackParam* pparams;
    
    ////////////////////////////////
    //first pair:
    if (iits1>=0 || itpc2>=0)
    {
        //ITS track
        ptrack = pEvent->GetTrack(iits1);
        constpparams = ptrack->GetOuterParam();
        if (!constpparams) return kFALSE;
        pparams = const_cast<AliExternalTrackParam*>(constpparams);
        SetRefSurface( pparams->GetX(), pparams->GetAlpha() );
        //pparams->PropagateTo(fLocalX, field);
        SetTrackParams1(pparams);
        //TPC track
        ptrack = pEvent->GetTrack(itpc1);
        constpparams = ptrack->GetInnerParam();
        if (!constpparams) return kFALSE;
        pparams = const_cast<AliExternalTrackParam*>(constpparams);
        pparams->Rotate(fAlpha);
        pparams->PropagateTo(fLocalX, field);
        SetTrackParams2(pparams);
        //do some accounting and update
        if (PrepareUpdate()) Update();
    }
    ////////////////////////////////
    //second pair:
    if (iits2>=0 || itpc2>=0)
    {
        //ITS track
        ptrack = pEvent->GetTrack(iits2);
        constpparams = ptrack->GetOuterParam();
        if (!constpparams) return kFALSE;
        pparams = const_cast<AliExternalTrackParam*>(constpparams);
        SetRefSurface( pparams->GetX(), pparams->GetAlpha() );
        //pparams->PropagateTo(fLocalX, field);
        SetTrackParams1(pparams);
        //TPC track
        ptrack = pEvent->GetTrack(itpc2);
        constpparams = ptrack->GetInnerParam();
        if (!constpparams) return kFALSE;
        pparams = const_cast<AliExternalTrackParam*>(constpparams);
        pparams->Rotate(fAlpha);
        pparams->PropagateTo(fLocalX, field);
        SetTrackParams2(pparams);
        //do some accounting and update
        if (PrepareUpdate()) Update();
    } 
    return kTRUE;
}

//______________________________________________________________________________
void AliRelAlignerKalman::Print(Option_t*) const
{
    //Print some useful info
    printf("\nAliRelAlignerKalman:\n");
    printf("  %i tracks, %i updates, %i outliers,", fNTracks, fNUpdates, fNOutliers );
    printf(" %i found cosmics in %i events\n", fNMatchedCosmics, fNProcessedEvents );
    printf("  psi(x):           % .3f ± (%.2f) mrad\n", 1e3*(*fPX)(0),1e3*TMath::Sqrt((*fPXcov)(0,0)));
    printf("  theta(y):         % .3f ± (%.2f) mrad\n", 1e3*(*fPX)(1),1e3*TMath::Sqrt((*fPXcov)(1,1)));
    printf("  phi(z):           % .3f ± (%.2f) mrad\n", 1e3*(*fPX)(2),1e3*TMath::Sqrt((*fPXcov)(2,2)));
    printf("  x:                % .3f ± (%.2f) micron\n", 1e4*(*fPX)(3),1e4*TMath::Sqrt((*fPXcov)(3,3)));
    printf("  y:                % .3f ± (%.2f) micron\n", 1e4*(*fPX)(4),1e4*TMath::Sqrt((*fPXcov)(4,4)));
    printf("  z:                % .3f ± (%.2f) micron\n", 1e4*(*fPX)(5),1e4*TMath::Sqrt((*fPXcov)(5,5)));
    printf("  TPC dftcorr       % .3g ± (%.2g) factor\n", (*fPX)(6),TMath::Sqrt((*fPXcov)(6,6)));
    printf("  TPC T0 offset     % .3f ± (%.2f) micron\n\n", 1e4*(*fPX)(7),1e4*TMath::Sqrt((*fPXcov)(7,7)));
    return;
}

//______________________________________________________________________________
void AliRelAlignerKalman::PrintCovarianceCorrection()
{
    //Print the measurement covariance correction matrix
    printf("Covariance correction matrix:\n");
    for ( Int_t i=0; i<fgkNMeasurementParams; i++ )
    {
        for ( Int_t j=0; j<i+1; j++ )
        {
            printf("% -2.2f  ", (*fPMeasurementCovCorr)(i,j) );
        }//for i
        printf("\n");
    }//for j
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
void AliRelAlignerKalman::SetTrackParams1( const AliExternalTrackParam* exparam )
{
    //Set the parameters for track in first volume
    fPTrackParam1 = exparam;
}

//______________________________________________________________________________
void AliRelAlignerKalman::SetTrackParams2( const AliExternalTrackParam* exparam )
{
    //Set the parameters for track in second volume
    fPTrackParam2 = exparam;
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
Bool_t AliRelAlignerKalman::PrepareUpdate()
{
    //Cast the extrapolated data (points and directions) into
    //the internal Kalman filter data representation.
    //takes the 3d coordinates of the points of intersection with
    //the reference surface and projects them onto a 2D plane.
    //does the same for angles, combines the result in one vector

    if (!PrepareMeasurement()) return kFALSE;
    if (!PrepareSystemMatrix()) return kFALSE;
    return kTRUE;
}

//______________________________________________________________________________
Bool_t AliRelAlignerKalman::Update()
{
    //perform the update - either kalman or calibration
    if (fCalibrationMode) return UpdateCalibration();
    if (fFillHistograms)
    {
        if (!UpdateCalibration()) return kFALSE;
        return UpdateEstimateKalman();
    }
    else return UpdateEstimateKalman();
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
    return;
}

//______________________________________________________________________________
Bool_t AliRelAlignerKalman::PrepareMeasurement()
{
    //Calculate the residuals and their covariance matrix
    if (!fPTrackParam1) return kFALSE;
    if (!fPTrackParam2) return kFALSE;
    const Double_t* pararr1 = fPTrackParam1->GetParameter();
    const Double_t* pararr2 = fPTrackParam2->GetParameter();

    //Take the track parameters and calculate the input to the Kalman filter
    (*fPMeasurement)(0) = pararr2[0]-pararr1[0];
    (*fPMeasurement)(1) = pararr2[1]-pararr1[1];
    (*fPMeasurement)(2) = pararr2[2]-pararr1[2];
    (*fPMeasurement)(3) = pararr2[3]-pararr1[3];
    fNTracks++; //count added track sets

    //the covariance
    const Double_t* parcovarr1 = fPTrackParam1->GetCovariance();
    const Double_t* parcovarr2 = fPTrackParam2->GetCovariance();
    (*fPMeasurementCov)(0,0)=parcovarr1[0];(*fPMeasurementCov)(0,1)=parcovarr1[1];(*fPMeasurementCov)(0,2)=parcovarr1[3];(*fPMeasurementCov)(0,3)=parcovarr1[6];
    (*fPMeasurementCov)(1,0)=parcovarr1[1];(*fPMeasurementCov)(1,1)=parcovarr1[2];(*fPMeasurementCov)(1,2)=parcovarr1[4];(*fPMeasurementCov)(1,3)=parcovarr1[7];
    (*fPMeasurementCov)(2,0)=parcovarr1[3];(*fPMeasurementCov)(2,1)=parcovarr1[4];(*fPMeasurementCov)(2,2)=parcovarr1[5];(*fPMeasurementCov)(2,3)=parcovarr1[8];
    (*fPMeasurementCov)(3,0)=parcovarr1[6];(*fPMeasurementCov)(3,1)=parcovarr1[7];(*fPMeasurementCov)(3,2)=parcovarr1[8];(*fPMeasurementCov)(3,3)=parcovarr1[9];
    (*fPMeasurementCov)(0,0)+=parcovarr2[0];(*fPMeasurementCov)(0,1)+=parcovarr2[1];(*fPMeasurementCov)(0,2)+=parcovarr2[3];(*fPMeasurementCov)(0,3)+=parcovarr2[6];
    (*fPMeasurementCov)(1,0)+=parcovarr2[1];(*fPMeasurementCov)(1,1)+=parcovarr2[2];(*fPMeasurementCov)(1,2)+=parcovarr2[4];(*fPMeasurementCov)(1,3)+=parcovarr2[7];
    (*fPMeasurementCov)(2,0)+=parcovarr2[3];(*fPMeasurementCov)(2,1)+=parcovarr2[4];(*fPMeasurementCov)(2,2)+=parcovarr2[5];(*fPMeasurementCov)(2,3)+=parcovarr2[8];
    (*fPMeasurementCov)(3,0)+=parcovarr2[6];(*fPMeasurementCov)(3,1)+=parcovarr2[7];(*fPMeasurementCov)(3,2)+=parcovarr2[8];(*fPMeasurementCov)(3,3)+=parcovarr2[9];
    if (fApplyCovarianceCorrection)
        *fPMeasurementCov += *fPMeasurementCovCorr;
    return kTRUE;
}

//______________________________________________________________________________
Bool_t AliRelAlignerKalman::PrepareSystemMatrix()
{
    //Calculate the system matrix for the Kalman filter
    //approximate the system using as reference the track in the first volume

    TVectorD z1( fgkNMeasurementParams );
    TVectorD z2( fgkNMeasurementParams );
    TVectorD x1( *fPX );
    TVectorD x2( *fPX );
    TMatrixD D( fgkNMeasurementParams, 1 );
    for ( Int_t i=0; i<fgkNSystemParams; i++ )
    {
        x1 = *fPX;
        x2 = *fPX;
        x1(i) -= fDelta[i]/2.;
        x2(i) += fDelta[i]/2.;
        if (!PredictMeasurement( z1, x1 )) return kFALSE;
        if (!PredictMeasurement( z2, x2 )) return kFALSE;
        for (Int_t j=0; j<fgkNMeasurementParams; j++ )
            D.GetMatrixArray()[j] = (z2.GetMatrixArray()[j]-z1.GetMatrixArray()[j])/fDelta[i];
        fPH->SetSub( 0, i, D );
    }
    return kTRUE;
}

//______________________________________________________________________________
Bool_t AliRelAlignerKalman::PredictMeasurement( TVectorD& pred, const TVectorD& state )
{
    // Implements a system model for the Kalman fit
    // pred is [dy,dz,dsinphi,dtanlambda]
    // state is [psi,theta,phi,x,y,z,driftTPC,offsetTPC]
    // note: the measurement is in a local frame, so the prediction also has to be
    // note: state is the misalignment in global reference system

    AliExternalTrackParam track(*fPTrackParam1); //make a copy of first track
    if (!MisalignTrack( &track, state )) return kFALSE;              //apply misalignments to get a prediction

    const Double_t* oldparam = fPTrackParam1->GetParameter();
    const Double_t* newparam = track.GetParameter();

    pred(0) = newparam[0] - oldparam[0];
    pred(1) = newparam[1] - oldparam[1];
    pred(2) = newparam[2] - oldparam[2];
    pred(3) = newparam[3] - oldparam[3];
    return kTRUE;
}

//______________________________________________________________________________
Bool_t AliRelAlignerKalman::UpdateEstimateKalman()
{
    //Kalman estimation of noisy constants: in the model A=1
    //The arguments are (following the usual convention):
    //  x - the state vector (parameters)
    //  P - the state covariance matrix (parameter errors)
    //  z - measurement vector
    //  R - measurement covariance matrix
    //  H - measurement model matrix ( z = Hx + v ) v being measurement noise (error fR)
    TVectorD *x = fPX;
    TMatrixDSym *P = fPXcov;
    TVectorD *z = fPMeasurement;
    TMatrixDSym *R = fPMeasurementCov;
    TMatrixD *H = fPH;
    
    TMatrixDSym I(TMatrixDSym::kUnit, *P);            //unit matrix
        
    //predict the state
    *P = *P + fQ*I;  //add some process noise (diagonal)
    
    // update prediction with measurement
        // calculate Kalman gain
        // K = PH/(HPH+R)
    TMatrixD PHT( *P, TMatrixD::kMultTranspose, *H );  //common factor (used twice)
    TMatrixD HPHT( *H, TMatrixD::kMult, PHT );
    HPHT += *R;
    TMatrixD K(PHT, TMatrixD::kMult, HPHT.Invert());                 //compute K
  
        // update the state and its covariance matrix
    TVectorD xupdate(*x);
    TVectorD Hx(*z);
    PredictMeasurement( Hx, *x );
    xupdate = K*((*z)-Hx);
    
    //SIMPLE OUTLIER REJECTION
    if ( IsOutlier( xupdate, *P ) && fRejectOutliers )
    {
        fNOutliers++;
        return kFALSE;
    }
    
    *x += xupdate;
    TMatrixD KH( K, TMatrixD::kMult, *H );
    TMatrixD IKH(I);
    IKH = I - KH;
    TMatrixD IKHP( IKH, TMatrixD::kMult, *P ); // (I-KH)P
    TMatrixDSym_from_TMatrixD( *P, IKHP );
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
        is = (is) || (TMath::Abs(update(i)) > fOutRejSigmas*TMath::Sqrt((covmatrix)(i,i)));
    }
    return is;
}

//______________________________________________________________________________
void AliRelAlignerKalman::TMatrixDSym_from_TMatrixD( TMatrixDSym& matsym, const TMatrixD& mat )
{
    //Produce a valid symmetric matrix out of an almost symmetric TMatrixD

    //not very efficient, diagonals are computer twice
    for (Int_t i=0; i<mat.GetNcols(); i++)
    {
        for (Int_t j=i; j<mat.GetNcols(); j++)
        {
            Double_t average = (mat(i,j)+mat(j,i))/2.;
            matsym(i,j)=average;
            matsym(j,i)=average;
        }
    }
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
            printf("% -1.2f  ", (*fPXcov)(i,j)/TMath::Sqrt( (*fPXcov)(i,i) * (*fPXcov)(j,j) ) );
        }//for j
        printf("\n");
    }//for i
    printf("\n");
    return;
}

//______________________________________________________________________________
Bool_t AliRelAlignerKalman::FindCosmicTrackletNumbersInEvent( Int_t& ITSgood1, Int_t& TPCgood1, Int_t& ITSgood2, Int_t& TPCgood2, const AliESDEvent* pEvent )
{
    //Find track point arrays belonging to one cosmic in a separately tracked ESD
    //and put them in the apropriate data members

    //Sanity cuts on tracks + check which tracks are ITS which are TPC
    Int_t ntracks = pEvent->GetNumberOfTracks(); //printf("number of tracks in event: %i\n", ntracks);
    if(ntracks<2) return kFALSE;
    Float_t* phiArr = new Float_t[ntracks];
    Float_t* thetaArr = new Float_t[ntracks];
    Double_t* distanceFromVertexArr = new Double_t[ntracks];
    Int_t* goodtracksArr = new Int_t[ntracks];
    Int_t* ITStracksArr = new Int_t[ntracks];
    Int_t* TPCtracksArr = new Int_t[ntracks];
    Int_t* nPointsArr = new Int_t[ntracks];
    Int_t nITStracks = 0;
    Int_t nTPCtracks = 0;
    Int_t nGoodTracks = 0;
    AliESDtrack* pTrack = NULL;
    
    const AliESDVertex* pVertex = pEvent->GetVertex();
    Double_t vertexposition[3];
    pVertex->GetXYZ(vertexposition);
    

    //select sane tracks
    for (Int_t itrack=0; itrack < ntracks; itrack++)
    {
        pTrack = pEvent->GetTrack(itrack);
        if (!pTrack) {std::cout<<"no track!"<<std::endl;continue;}
        if(pTrack->GetNcls(0)+pTrack->GetNcls(1) < fMinPointsVol1) continue;
        Float_t phi = pTrack->GetAlpha()+TMath::ASin(pTrack->GetSnp());
        Float_t theta = 0.5*TMath::Pi()-TMath::ATan(pTrack->GetTgl());
        //printf("phi: %4.2f theta: %4.2f\n", phi, theta);
        if(fCuts)
        {
            if(pTrack->GetP()<fMinMom || pTrack->GetP()>fMaxMom) continue;
            Float_t abssinphi = TMath::Abs(TMath::Sin(phi));
            if(abssinphi<fMinAbsSinPhi || abssinphi>fMaxAbsSinPhi) continue;
            Float_t sintheta = TMath::Sin(theta);
            if(sintheta<fMinSinTheta || sintheta>fMaxSinTheta) continue;
        }
        goodtracksArr[nGoodTracks]=itrack;
        phiArr[nGoodTracks]=phi;
        thetaArr[nGoodTracks]=theta;

        //Double_t magfield = pEvent->GetMagneticField();
        //pTrack->RelateToVertex( pVertex, magfield, 10000. );
        //Double_t trackposition[3];
        //pTrack->GetXYZ( trackposition );
        //distanceFromVertexArr[nGoodTracks] = 
        //      TMath::Sqrt((trackposition[0]-vertexposition[0])*(trackposition[0]-vertexposition[0])
        //               + (trackposition[1]-vertexposition[1])*(trackposition[1]-vertexposition[1])
        //               + (trackposition[2]-vertexposition[2])*(trackposition[2]-vertexposition[2]));

        //check if track is ITS or TPC
        Int_t nClsITS = pTrack->GetNcls(0);
        Int_t nClsTPC = pTrack->GetNcls(1);
        if ( ((pTrack->GetStatus()&AliESDtrack::kITSout)>0)&&
             !((pTrack->GetStatus()&AliESDtrack::kTPCin)>0)&&
             !(nClsITS<fMinPointsVol1) )  //enough points
        {
            ITStracksArr[nITStracks] = nGoodTracks;
            nITStracks++;
        }

        if ( ((pTrack->GetStatus()&AliESDtrack::kTPCin)>0)&&
             !(nClsTPC<fMinPointsVol2) )  //enough points
        {
            TPCtracksArr[nTPCtracks] = nGoodTracks;
            nTPCtracks++;
        }
        //if(nClsITS>=2 && nClsTPC==0)
        //{ // ITS SA
        //    ITStracksArr[nITStracks] = nGoodTracks;
        //    nITStracks++;
        //}
        //if(nClsTPC>=50)
        //{ // TPC
        //    TPCtracksArr[nTPCtracks] = nGoodTracks;
        //    nTPCtracks++;
        //}

        nGoodTracks++;
    }//for itrack   -   sanity cuts

    //printf("TrackFinder: there are good tracks, %d in ITS and %d TPC.\n", nITStracks, nTPCtracks);
    
    if( nITStracks < 2 || nTPCtracks < 2 )
    {
        delete [] goodtracksArr; goodtracksArr=0;
        delete [] ITStracksArr; ITStracksArr=0;
        delete [] TPCtracksArr; TPCtracksArr=0;
        delete [] nPointsArr; nPointsArr=0;
        delete [] phiArr; phiArr=0;
        delete [] thetaArr; thetaArr=0;
        delete [] distanceFromVertexArr; distanceFromVertexArr=0;
        //printf("TrackFinder: not enough tracks in ITS or TPC\n");
        return kFALSE;
    }

    //find matching in TPC
    Float_t min = 10000000.;
    TPCgood1 = -1;
    TPCgood2 = -1;
    for(Int_t itr1=0; itr1<nTPCtracks; itr1++)
    {
        for(Int_t itr2=itr1+1; itr2<nTPCtracks; itr2++)
        {
            Float_t deltatheta = TMath::Abs(TMath::Pi()-thetaArr[TPCtracksArr[itr1]]-thetaArr[TPCtracksArr[itr2]]);
            if(deltatheta > fMaxMatchingAngle) continue;
            Float_t deltaphi = TMath::Abs(TMath::Abs(phiArr[TPCtracksArr[itr1]]-phiArr[TPCtracksArr[itr2]])-TMath::Pi());
            if(deltaphi > fMaxMatchingAngle) continue;
            //printf("ITS: %f  %f     %f  %f\n",deltaphi,deltatheta,thetaArr[ITStracksArr[itr1]],thetaArr[ITStracksArr[itr2]]);
            if(deltatheta+deltaphi<min) //only the best matching pair
            {
                min=deltatheta+deltaphi;
                TPCgood1 = TPCtracksArr[itr1];  //store the index of track in goodtracksArr[]
                TPCgood2 = TPCtracksArr[itr2];
            }
        }
    }
    if (TPCgood1 < 0) //no dubble cosmic track
    {
        delete [] goodtracksArr; goodtracksArr=0;
        delete [] ITStracksArr; ITStracksArr=0;
        delete [] TPCtracksArr; TPCtracksArr=0;
        delete [] nPointsArr; nPointsArr=0;
        delete [] phiArr; phiArr=0;
        delete [] thetaArr; thetaArr=0;
        delete [] distanceFromVertexArr; distanceFromVertexArr=0;
        //printf("TrackFinder: no cosmic pair inside ITS\n");
        return kFALSE;
    }

    //find for the first TPC track the matching ITS track
    ITSgood1 = -1;
    min = 10000000.;
    for(Int_t i=0; i<nITStracks; i++)
    {
        Float_t deltatheta = TMath::Abs(thetaArr[TPCgood1]-thetaArr[ITStracksArr[i]]);
        if(deltatheta > fMaxMatchingAngle) continue;
        Float_t deltaphi = TMath::Abs(phiArr[TPCgood1]-phiArr[ITStracksArr[i]]);
        if(deltaphi > fMaxMatchingAngle) continue;
        //printf("ITS: %f  %f     %f  %f\n",deltaphi,deltatheta,thetaArr[ITStracksArr[itr1]],thetaArr[ITStracksArr[itr2]]);
        if(deltatheta+deltaphi<min) //only the best matching pair
        {
            min=deltatheta+deltaphi;
            ITSgood1 = ITStracksArr[i];  //store the index of track in goodtracksArr[]
        }
    }

    //find for the second TPC track the matching ITS track
    ITSgood2 = -1;
    min = 10000000.;
    for(Int_t i=0; i<nITStracks; i++)
    {
        Float_t deltatheta = TMath::Abs(thetaArr[TPCgood2]-thetaArr[ITStracksArr[i]]);
        if(deltatheta > fMaxMatchingAngle) continue;
        Float_t deltaphi = TMath::Abs(phiArr[TPCgood2]-phiArr[ITStracksArr[i]]);
        if(deltaphi > fMaxMatchingAngle) continue;
        //printf("ITS: %f  %f     %f  %f\n",deltaphi,deltatheta,thetaArr[ITStracksArr[itr1]],thetaArr[ITStracksArr[itr2]]);
        if(deltatheta+deltaphi<min) //only the best matching pair
        {
            min=deltatheta+deltaphi;
            ITSgood2 = ITStracksArr[i];  //store the index of track in goodtracksArr[]
        }
    }
    
    if ((ITSgood1 < 0) && (ITSgood2 < 0))
    {
        delete [] goodtracksArr; goodtracksArr=0;
        delete [] ITStracksArr; ITStracksArr=0;
        delete [] TPCtracksArr; TPCtracksArr=0;
        delete [] nPointsArr; nPointsArr=0;
        delete [] phiArr; phiArr=0;
        delete [] thetaArr; thetaArr=0;
        delete [] distanceFromVertexArr; distanceFromVertexArr=0;
        return kFALSE;
    }

    //we found a cosmic
    fNMatchedCosmics++;
    
    ///////////////////////////////////////////////////////////////////////////
    // convert indexes from local goodtrackarrays to global track index
    TPCgood1 = goodtracksArr[TPCgood1];
    TPCgood2 = goodtracksArr[TPCgood2];
    ITSgood1 = (ITSgood1==-1) ? -1 : goodtracksArr[ITSgood1];
    ITSgood2 = (ITSgood2==-1) ? -1 : goodtracksArr[ITSgood2];
    ///////////////////////////////////////////////////////////////////////////

    delete [] goodtracksArr;
    delete [] ITStracksArr;
    delete [] TPCtracksArr;
    delete [] nPointsArr;
    delete [] phiArr;
    delete [] thetaArr;
    delete [] distanceFromVertexArr;
    return kTRUE;
}
//_______________________________________________________________________________
Bool_t AliRelAlignerKalman::UpdateCalibration()
{
    //Update the calibration with new data (in calibration mode)

    fPMes0Hist->Fill( (*fPMeasurement)(0) );
    fPMes1Hist->Fill( (*fPMeasurement)(1) );
    fPMes2Hist->Fill( (*fPMeasurement)(2) );
    fPMes3Hist->Fill( (*fPMeasurement)(3) );
    fPMesErr0Hist->Fill( TMath::Sqrt((*fPMeasurementCov)(0,0)) );
    fPMesErr1Hist->Fill( TMath::Sqrt((*fPMeasurementCov)(1,1)) );
    fPMesErr2Hist->Fill( TMath::Sqrt((*fPMeasurementCov)(2,2)) );
    fPMesErr3Hist->Fill( TMath::Sqrt((*fPMeasurementCov)(3,3)) );
    return kTRUE;
}

//______________________________________________________________________________
Bool_t AliRelAlignerKalman::SetCalibrationMode( Bool_t cp )
{
    //sets the calibration mode
    if (cp)
    {
        fCalibrationMode=kTRUE;
        return kTRUE;
    }//if (cp)
    else
    {
        if (fCalibrationMode) // do it only after the calibration pass
        {
            CalculateCovarianceCorrection();
            SetApplyCovarianceCorrection();
            fCalibrationMode=kFALSE;
            return kTRUE;
        }//if (fCalibrationMode)
    }//else (cp)
    return kFALSE;
}

//______________________________________________________________________________
Bool_t AliRelAlignerKalman::CalculateCovarianceCorrection()
{
    //Calculates the correction to the measurement covariance
    //using the calibration histograms

    fPMeasurementCovCorr->Zero(); //reset the correction

    Double_t s,m,c;  //sigma,meansigma,correction

    //TF1* fitformula;
    //fPMes0Hist->Fit("gaus");
    //fitformula = fPMes0Hist->GetFunction("gaus");
    //s = fitformula->GetParameter(2);   //spread of the measurement
    //fPMesErr0Hist->Fit("gaus");
    //fitformula = fPMesErr0Hist->GetFunction("gaus"); //average error from cov matrices
    //m = fitformula->GetParameter(1);
    s = fPMes0Hist->GetRMS();
    m = fPMesErr0Hist->GetMean();
    c = s-m; //the difference between the average error and real spread of the data
    if (c>0) //only correct is spread bigger than average error
        (*fPMeasurementCovCorr)(0,0) = c*c;

    //fPMes1Hist->Fit("gaus");
    //fitformula = fPMes1Hist->GetFunction("gaus");
    //s = fitformula->GetParameter(2);
    //fPMesErr1Hist->Fit("gaus");
    //fitformula = fPMesErr1Hist->GetFunction("gaus");
    //m = fitformula->GetParameter(1);
    s = fPMes1Hist->GetRMS();
    m = fPMesErr1Hist->GetMean();
    c = s-m;
    if (c>0) //only correct is spread bigger than average error
        (*fPMeasurementCovCorr)(1,1) = c*c;

    //fPMes2Hist->Fit("gaus");
    //fitformula = fPMes2Hist->GetFunction("gaus");
    //s = fitformula->GetParameter(2);
    //fPMesErr2Hist->Fit("gaus");
    //fitformula = fPMesErr2Hist->GetFunction("gaus");
    //m = fitformula->GetParameter(1);
    s = fPMes2Hist->GetRMS();
    m = fPMesErr2Hist->GetMean();
    c = s-m;
    if (c>0) //only correct is spread bigger than average error
        (*fPMeasurementCovCorr)(2,2) = c*c;

    //fPMes3Hist->Fit("gaus");
    //fitformula = fPMes3Hist->GetFunction("gaus");
    //s = fitformula->GetParameter(2);
    //fPMesErr3Hist->Fit("gaus");
    //fitformula = fPMesErr3Hist->GetFunction("gaus");
    //m = fitformula->GetParameter(1);
    s = fPMes3Hist->GetRMS();
    m = fPMesErr3Hist->GetMean();
    c = s-m;
    if (c>0) //only correct is spread bigger than average error
        (*fPMeasurementCovCorr)(3,3) = c*c;

    return kTRUE;
}

//______________________________________________________________________________
void AliRelAlignerKalman::PrintDebugInfo()
{
    //prints some debug info
    Print();
    std::cout<<"AliRelAlignerKalman debug info"<<std::endl;
    printf("TrackParams1:");
    fPTrackParam1->Print();
    printf("TrackParams2:");
    fPTrackParam2->Print();
    printf("Measurement:");
    fPMeasurement->Print();
    printf("Measurement covariance:");
    fPMeasurementCov->Print();
}

//______________________________________________________________________________
AliRelAlignerKalman::AliRelAlignerKalman(const AliRelAlignerKalman& a):
    TObject(a),
    fAlpha(a.fAlpha),
    fLocalX(a.fLocalX),
    fPTrackParam1(a.fPTrackParam1),
    fPTrackParam2(a.fPTrackParam2),
    fPX(new TVectorD( *a.fPX )),
    fPXcov(new TMatrixDSym( *a.fPXcov )),
    fPH(new TMatrixD( *a.fPH )),
    fQ(a.fQ),
    fPMeasurement(new TVectorD( *a.fPMeasurement )),
    fPMeasurementCov(new TMatrixDSym( *a.fPMeasurementCov )),
    fOutRejSigmas(a.fOutRejSigmas),
    fRejectOutliers(a.fRejectOutliers),
    fCalibrationMode(a.fCalibrationMode),
    fFillHistograms(a.fFillHistograms),
    fRequireDoubleTPCtrack(a.fRequireDoubleTPCtrack),
    fApplyCovarianceCorrection(a.fApplyCovarianceCorrection),
    fCuts(a.fCuts),
    fMinPointsVol1(a.fMinPointsVol1),
    fMinPointsVol2(a.fMinPointsVol2),
    fMinMom(a.fMinMom),
    fMaxMom(a.fMaxMom),
    fMinAbsSinPhi(a.fMinAbsSinPhi),
    fMaxAbsSinPhi(a.fMaxAbsSinPhi),
    fMinSinTheta(a.fMinSinTheta),
    fMaxSinTheta(a.fMaxSinTheta),
    fMaxMatchingAngle(a.fMaxMatchingAngle),
    fMaxMatchingDistance(a.fMaxMatchingDistance),  //in cm
    fNTracks(a.fNTracks),
    fNUpdates(a.fNUpdates),
    fNOutliers(a.fNOutliers),
    fNMatchedCosmics(a.fNMatchedCosmics),
    fNProcessedEvents(a.fNProcessedEvents),
    fPMes0Hist(new TH1D(*a.fPMes0Hist)),
    fPMes1Hist(new TH1D(*a.fPMes1Hist)),
    fPMes2Hist(new TH1D(*a.fPMes2Hist)),
    fPMes3Hist(new TH1D(*a.fPMes3Hist)),
    fPMesErr0Hist(new TH1D(*a.fPMesErr0Hist)),
    fPMesErr1Hist(new TH1D(*a.fPMesErr1Hist)),
    fPMesErr2Hist(new TH1D(*a.fPMesErr2Hist)),
    fPMesErr3Hist(new TH1D(*a.fPMesErr3Hist)),
    fPMeasurementCovCorr(new TMatrixDSym(*a.fPMeasurementCovCorr))
{
    //copy constructor
    memcpy(fDelta,a.fDelta,fgkNSystemParams*sizeof(Double_t));
}

//______________________________________________________________________________
AliRelAlignerKalman& AliRelAlignerKalman::operator=(const AliRelAlignerKalman& a)
{
    //assignment operator
    fAlpha=a.fAlpha;
    fLocalX=a.fLocalX;
    fPTrackParam1=a.fPTrackParam1;
    fPTrackParam2=a.fPTrackParam2;
    *fPX = *a.fPX;
    *fPXcov = *a.fPXcov;
    *fPH = *a.fPH;
    fQ=a.fQ;
    *fPMeasurement=*a.fPMeasurement;
    *fPMeasurementCov=*a.fPMeasurementCov;
    fOutRejSigmas=a.fOutRejSigmas;
    memcpy(fDelta,a.fDelta,fgkNSystemParams*sizeof(Double_t));
    fRejectOutliers=a.fRejectOutliers;
    fCalibrationMode=a.fCalibrationMode;
    fFillHistograms=a.fFillHistograms;
    fRequireDoubleTPCtrack=a.fRequireDoubleTPCtrack;
    fApplyCovarianceCorrection=a.fApplyCovarianceCorrection;
    fCuts=a.fCuts;
    fMinPointsVol1=a.fMinPointsVol1;
    fMinPointsVol2=a.fMinPointsVol2;
    fMinMom=a.fMinMom;
    fMaxMom=a.fMaxMom;
    fMinAbsSinPhi=a.fMinAbsSinPhi;
    fMaxAbsSinPhi=a.fMaxAbsSinPhi;
    fMinSinTheta=a.fMinSinTheta;
    fMaxSinTheta=a.fMaxSinTheta;
    fMaxMatchingAngle=a.fMaxMatchingAngle;
    fMaxMatchingDistance=a.fMaxMatchingDistance,  //in c;
    fNTracks=a.fNTracks;
    fNUpdates=a.fNUpdates;
    fNOutliers=a.fNOutliers;
    fNMatchedCosmics=a.fNMatchedCosmics;
    fNProcessedEvents=a.fNProcessedEvents;
    *fPMes0Hist=*a.fPMes0Hist;
    *fPMes1Hist=*a.fPMes1Hist;
    *fPMes2Hist=*a.fPMes2Hist;
    *fPMes3Hist=*a.fPMes3Hist;
    *fPMesErr0Hist=*a.fPMesErr0Hist;
    *fPMesErr1Hist=*a.fPMesErr1Hist;
    *fPMesErr2Hist=*a.fPMesErr2Hist;
    *fPMesErr3Hist=*a.fPMesErr3Hist;
    *fPMeasurementCovCorr=*a.fPMeasurementCovCorr;
    return *this;
}

//______________________________________________________________________________
Bool_t AliRelAlignerKalman::MisalignTrack( AliExternalTrackParam* tr, const TVectorD& misal )
{
    //Misalign the track
    
    Double_t x = tr->GetX();
    Double_t alpha = tr->GetAlpha();
    Double_t point[3],dir[3];
    tr->GetXYZ(point);
    tr->GetDirection(dir);
    TVector3 Point(point);
    TVector3 Dir(dir);
    
    //Misalign track
    //TPC drift correction
    Point(2) *= misal(6);
    Dir(2) *= misal(6);
    Dir=Dir.Unit(); //to be safe
    //TPC offset
    if (Point(2)>0) Point(2) += misal(7);
    else Point(2) -= misal(7);
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
    if (dir[0]==0) return kFALSE;
    Double_t s=(x-point[0])/dir[0];
    p[0] = point[1]+s*dir[1];
    p[1] = point[2]+s*dir[2];
    Double_t pt = TMath::Sqrt(dir[0]*dir[0]+dir[1]*dir[1]);
    if (pt==0) return kFALSE;
    p[2] = dir[1]/pt;
    p[3] = dir[2]/pt;

    const Double_t* pcovtmp = tr->GetCovariance();
    memcpy(pcov,pcovtmp,15*sizeof(Double_t));

    tr->Set(x,alpha,p,pcov);
    return kTRUE;
}

//______________________________________________________________________________
void AliRelAlignerKalman::ResetCovariance()
{
    //Resets the covariance of the fit to a default value
    fPXcov->UnitMatrix();
    (*fPXcov)(0,0) = .1*.1; //psi (rad)
    (*fPXcov)(1,1) = .1*.1; //theta (rad
    (*fPXcov)(2,2) = .1*.1; //phi (rad)
    (*fPXcov)(3,3) = 1.*1.; //x (cm)
    (*fPXcov)(4,4) = 1.*1.; //y (cm)
    (*fPXcov)(5,5) = 1.*1.; //z (cm)
    (*fPXcov)(6,6) = .1*.1;//drift slope correction (fraction: 1. is 100%)
    (*fPXcov)(7,7) = 1.*1.; //offset (cm)
}
