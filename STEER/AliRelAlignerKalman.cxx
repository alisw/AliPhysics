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
//    Aligns two tracking volumes (by default ITS and TPC)
//    Determines the transformation of the second volume (TPC) wrt the first (ITS)
//    additionally calculates some callibration parameters for TPC
//    Fit parameters are:
//    - 3 shifts, x,y,z
//    - 3 Cardan angles, psi, theta, phi (see definition in alignment docs),
//    - TPC drift velocity correction (vel_true = vel*(1+correction)),
//    - TPC offset correction.
//
//    Basic usage:
//    When aligning two volumes at any given time a single instance of
//    the class should be active. The fit of the parameters is updated
//    by adding new data using one of the Add.... methods.
//
//    User methods:
//    In collision events add an ESD track to update the fit,
//    track will be automatically split in ITS and TPC tracklets
//    and the fit of the alignment constants will be updated:
//    
//        Bool_t AddESDTrack( AliESDtrack* pTrack );
//    
//    You may also give it the AliTrackPointArray,
//    the fit will be updated if there are enough ITS and TPC
//    points in it. Again, all is done automatically:
//    
//        Bool_t AddTrackPointArray( AliTrackPointArray* pTrackPointArr );
//    
//    For cosmic data, the assumption is that the tracking is done separately
//    for ITS and TPC and the tracklets are saved inside one AliESDEvent. The
//    method
//    
//        Bool_t AddCosmicEventSeparateTracking( AliESDEvent* pEvent );
//    
//    then searches the event for matching tracklets and upon succes it updates
//    the fit. For cosmics the fitter switches to a 3 tracks mode, the 3 tracks
//    being: 1 ITS and 2 TPC tracklets (upper and lower). 
//    
//    Experts only:
//    
//
//    Origin: Mikolaj Krzewicki, Nikhef, Mikolaj.Krzewicki@cern.ch
//
//////////////////////////////////////////////////////////////////////////////

#include "AliRelAlignerKalman.h"

ClassImp(AliRelAlignerKalman)

//########################################################
//
//    Control section
//
//########################################################

AliRelAlignerKalman::AliRelAlignerKalman():
    fPTrackPointArray1(new AliTrackPointArray(1)),
    fPTrackPointArray2(new AliTrackPointArray(1)),
    fPTrackPointArray2b(new AliTrackPointArray(1)),
    fPVolIDArray1(new TArrayI(1)),
    fPVolIDArray2(new TArrayI(1)),
    fPVolIDArray2b(new TArrayI(1)),
    fQ(1e-10),
    fNMeasurementParams(fgkNMeasurementParams2TrackMode),
    fNSystemParams(fgkNSystemParams),
    fOutRejSigmas(2.),
    f3TracksMode(kFALSE),
    fSortTrackPointsWithY(kTRUE),
    fFixedMeasurementCovariance(kTRUE),
    fCalibrationPass(kFALSE),
    fCalibrationUpdated(kFALSE),
    fFitted(kFALSE),
    fCuts(kFALSE),
    fNTracks(0),
    fNUpdates(0),
    fNOutliers(0),
    fNUnfitted(0),
    fNMatchedCosmics(0),
    fMinPointsVol1(2),
    fMinPointsVol2(10),
    fMinMom(0.),
    fMaxMom(1e100),
    fMinAbsSinPhi(0.),
    fMaxAbsSinPhi(1.),
    fMinSinTheta(-1.),
    fMaxSinTheta(1.),
    fMaxMatchingAngle(0.1),
    fMaxMatchingDistance(1.),  //in cm
    fFirstLayerVol1(1),
    fLastLayerVol1(6),
    fFirstLayerVol2(7),
    fLastLayerVol2(8)
{
    //Default constructor
    
    fPDir1 = new TVector3; //Direction
    fPDir2 = new TVector3;
    fPDir2b = new TVector3;
    fPDir1Cov = new TMatrixDSym(3);
    fPDir2Cov = new TMatrixDSym(3);
    fPDir2bCov = new TMatrixDSym(3);
    fPDir1ThetaPhiCov = new TMatrixDSym(2); //Direction vectors are unit vectors, therefore 2 independent components!
    fPDir2ThetaPhiCov = new TMatrixDSym(2);
    fPDir2bThetaPhiCov = new TMatrixDSym(2);
    fPPoint1 = new TVector3;
    fPPoint2 = new TVector3;
    fPPoint2b = new TVector3;
    fPPoint1Cov = new TMatrixDSym(3);
    fPPoint2Cov = new TMatrixDSym(3);
    fPPoint2bCov = new TMatrixDSym(3);
    fPRefPoint = new TVector3(0.,0.,0.);
    fPRefPointDir = new TVector3(0.,1.,0.);
    Float_t refpointcov[6] = { 1., 0., 0., 
                                   0., 0.,
                                       1.  };
    fPRefAliTrackPoint = new AliTrackPoint( 0.,0.,0.,refpointcov,0 );
    fPMeasurement = new TVectorD( fNMeasurementParams );
    fPMeasurementCov = new TMatrixDSym( fNMeasurementParams );
    fPX = new TVectorD( fNSystemParams );
    fPXcov = new TMatrixDSym( fNSystemParams );
    fPH = new TMatrixD( fNMeasurementParams, fNSystemParams );
    fPRotMat = new TMatrixD(3,3);
    fPVec010 = new TVector3(0,1,0);

    //default seed: zero, unit(large) errors
    fPXcov->UnitMatrix();
    (*fPXcov) *= 1.;

    fPPhiMesHist = new TH1D("phi","phi", 50, 0, 0);
    fPThetaMesHist = new TH1D("theta","theta", 50, 0, 0);
    fPXMesHist = new TH1D("x","x", 50, 0, 0);
    fPZMesHist = new TH1D("z","z", 50, 0, 0);
    fPPhiMesHist2 = new TH1D("phi_","phi_", 50, 0, 0);
    fPThetaMesHist2 = new TH1D("theta_","theta_", 50, 0, 0);
    fPXMesHist2 = new TH1D("x_","x_", 50, 0, 0);
    fPZMesHist2 = new TH1D("z_","z_", 50, 0, 0);
}

Bool_t AliRelAlignerKalman::AddESDTrack( AliESDtrack* pTrack )
{
    //Add an ESD track from a central event
    //from the ESDtrack extract the pointarray

    SetCosmicEvent(kFALSE);
    
    const AliTrackPointArray *pTrackPointArrayConst = pTrack->GetTrackPointArray();
    AliTrackPointArray* pTrackPointArray = const_cast < AliTrackPointArray* > (pTrackPointArrayConst);    
    if (pTrackPointArray == 0) return kFALSE;

    if (!AddTrackPointArray( pTrackPointArray )) return kFALSE;
    return kTRUE;
}

Bool_t AliRelAlignerKalman::AddTrackPointArray( AliTrackPointArray* pTrackPointArray )
{
    //Add a trackpointarray from a central event
    
    SetCosmicEvent(kFALSE);
    fNTracks++; 

    if (!ExtractSubTrackPointArray( fPTrackPointArray1, fPVolIDArray1, pTrackPointArray, 1,6 ))
        return kFALSE;
    if (!ExtractSubTrackPointArray( fPTrackPointArray2, fPVolIDArray2, pTrackPointArray, 7,8 ))
        return kFALSE;
   
    if (!ProcessTrackPointArrays()) return kFALSE; 
    if (!PrepareUpdate()) return kFALSE;
    if (!Update()) return kFALSE;
    
    return kTRUE;
}

Bool_t AliRelAlignerKalman::AddCosmicEventSeparateTracking( AliESDEvent* pEvent )
{
    //Add an cosmic with separately tracked ITS and TPC parts, do trackmatching

    SetCosmicEvent(kTRUE);  //set the cosmic flag
    fNTracks++; 

    if (!FindCosmicInEvent( pEvent )) return kFALSE;
    if (!ProcessTrackPointArrays()) return kFALSE; 
    //PrintDebugInfo();
    if (!PrepareUpdate()) return kFALSE;
    if (!Update()) return kFALSE;
    printf("fpMeasurementCov:\n");
    fPMeasurementCov->Print();

    return kTRUE;
}

void AliRelAlignerKalman::Print()
{
    //Print some useful info
    printf("\nAliRelAlignerKalman:\n");
    printf("  %i tracks, %i updates, %i outliers, ", fNTracks, fNUpdates, fNOutliers );
    printf("%i not fitted, %i matched cosmic tracklets\n", fNUnfitted, fNMatchedCosmics );
    printf("  psi(x):         %.1f ± (%.1f) mrad\n", 1e3*(*fPX)(0),1e3*TMath::Sqrt((*fPXcov)(0,0)));
    printf("  theta(y):       %.1f ± (%.1f) mrad\n", 1e3*(*fPX)(1),1e3*TMath::Sqrt((*fPXcov)(1,1)));
    printf("  phi(z):         %.1f ± (%.1f) mrad\n", 1e3*(*fPX)(2),1e3*TMath::Sqrt((*fPXcov)(2,2)));
    printf("  x:              %.1f ± (%.1f) micron\n", 1e4*(*fPX)(3),1e4*TMath::Sqrt((*fPXcov)(3,3)));
    printf("  y:              %.1f ± (%.1f) micron\n", 1e4*(*fPX)(4),1e4*TMath::Sqrt((*fPXcov)(4,4)));
    printf("  z:              %.1f ± (%.1f) micron\n", 1e4*(*fPX)(5),1e4*TMath::Sqrt((*fPXcov)(5,5)));
    printf("  TPC(drift corr) %.1f ± (%.1f) percent\n", 1e2*(*fPX)(6),1e2*TMath::Sqrt((*fPXcov)(6,6)));
    printf("  TPC(offset)     %.1f ± (%.1f) micron\n", 1e4*(*fPX)(7),1e4*TMath::Sqrt((*fPXcov)(7,7)));
    return;
}

void AliRelAlignerKalman::SetTrackParams1( const TVector3& point, const TMatrixDSym& pointcov, const TVector3& dir, const TMatrixDSym& dircov )
{
    //Set the trackparameters for track in the first volume at reference
    *fPPoint1 = point;
    *fPPoint1Cov = pointcov;
    *fPDir1 = dir;
    *fPDir1Cov = dircov;
}

void AliRelAlignerKalman::SetTrackParams2( const TVector3& point, const TMatrixDSym& pointcov, const TVector3& dir, const TMatrixDSym& dircov )
{
    //Set the trackparameters for track in the second volume at reference
    *fPPoint2 = point;
    *fPPoint2Cov = pointcov;
    *fPDir2 = dir;
    *fPDir2Cov = dircov;
}

void AliRelAlignerKalman::SetTrackParams2b( const TVector3& point, const TMatrixDSym& pointcov, const TVector3& dir, const TMatrixDSym& dircov )
{
    //Set the trackparameters for second track in the second volume at reference
    *fPPoint2b = point;
    *fPPoint2bCov = pointcov;
    *fPDir2b = dir;
    *fPDir2bCov = dircov;
}   

void AliRelAlignerKalman::SetRefSurface( const Double_t yref )
{
    //set the reference surface
    fPRefPoint->SetXYZ( 0., 1.*yref, 0. );
    fPRefPointDir->SetXYZ( 0., 1., 0. );
    Float_t refpointcov[6] = { 1., 0., 0., 
                                   0., 0.,
                                       1.  };
    fPRefAliTrackPoint->SetXYZ( 0., 1.*yref, 0., refpointcov );
}

void AliRelAlignerKalman::SetRefSurface( const TVector3& point, const TVector3& dir )
{
    //set the reference surface
    *fPRefPoint = point;
    *fPRefPointDir = dir;
}

void AliRelAlignerKalman::SetRefSurface( const AliTrackPoint& point, const Bool_t horizontal )
{
    //set the reference surface
    (*fPRefAliTrackPoint) = point;

    (*fPRefPoint)(0) = point.GetX();
    (*fPRefPoint)(1) = point.GetY();
    (*fPRefPoint)(2) = point.GetZ();
    
    if (horizontal)
    {
        (*fPRefPointDir)(0) = 0.;
        (*fPRefPointDir)(1) = 1.;
        (*fPRefPointDir)(2) = 0.;
        
        Float_t refpointcov[6] = { 1., 0., 0., 
                                       0., 0.,
                                           1.  };
        fPRefAliTrackPoint->SetXYZ( point.GetX(), point.GetY() , point.GetZ(), refpointcov );

    }
    else
    {
        Double_t angle = point.GetAngle();
        (*fPRefPointDir)(0) = TMath::Cos(angle);
        (*fPRefPointDir)(1) = TMath::Sin(angle);
        (*fPRefPointDir)(2) = 0.;
        *fPRefPointDir = fPRefPointDir->Unit();
    }
}

void AliRelAlignerKalman::Set3TracksMode( Bool_t mode )
{
    //set the mode where 2 tracklets are allowed in volume 2 (default:TPC)
    //used for alignment with cosmics

    f3TracksMode = mode;
    if (mode)
    {
        if (fNMeasurementParams!=fgkNMeasurementParams3TrackMode) //do something only if necessary
        {
            fNMeasurementParams = fgkNMeasurementParams3TrackMode;
            delete fPMeasurement;
            delete fPMeasurementCov;
            delete fPH;
            fPMeasurement = new TVectorD( fNMeasurementParams );
            fPMeasurementCov = new TMatrixDSym( fNMeasurementParams );
            fPH = new TMatrixD( fNMeasurementParams, fNSystemParams );
        }
    }
    else
    {
        if (fNMeasurementParams!=fgkNMeasurementParams2TrackMode)
        {
            fNMeasurementParams = fgkNMeasurementParams2TrackMode;
            delete fPMeasurement;
            delete fPMeasurementCov;
            delete fPH;
            fPMeasurement = new TVectorD( fNMeasurementParams );
            fPMeasurementCov = new TMatrixDSym( fNMeasurementParams );
            fPH = new TMatrixD( fNMeasurementParams, fNSystemParams );
        }
    }
}
    
Bool_t AliRelAlignerKalman::PrepareUpdate()
{
    //Cast the extrapolated data (points and directions) into
    //the internal Kalman filter data representation.
    //takes the 3d coordinates of the points of intersection with
    //the reference surface and projects them onto a 2D plane.
    //does the same for angles, combines the result in one vector

    if (!FillMeasurement()) return kFALSE;
    if (!FillMeasurementMatrix()) return kFALSE;
    return kTRUE;
}

Bool_t AliRelAlignerKalman::Update()
{
    //perform the update - either kalman or calibration
    if (fCalibrationPass) return UpdateCalibration();
    else return UpdateEstimateKalman();
}

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

Bool_t AliRelAlignerKalman::FillMeasurement()
{
    //Take the track parameters and calculate the input to the Kalman filter

    //Direction vectors 
    //Measured is the difference between the polar angles - convert to that
    (*fPMeasurement)(0) = fPDir2->Theta() - fPDir1->Theta();
    (*fPMeasurement)(1) = fPDir2->Phi() - fPDir1->Phi();
    
    //Points
    //Measured is the distance in XZ plane at Y of the reference point
    (*fPMeasurement)(2) = fPPoint2->X() - fPPoint1->X();
    (*fPMeasurement)(3) = fPPoint2->Z() - fPPoint1->Z();

    //if 3 track mode set, set also the stuff for 2nd TPC track
    if (f3TracksMode)
    {
        (*fPMeasurement)(4) = fPDir2b->Theta() - fPDir1->Theta();
        (*fPMeasurement)(5) = fPDir2b->Phi() - fPDir1->Phi();
        (*fPMeasurement)(6) = fPPoint2b->X() - fPPoint1->X();
        (*fPMeasurement)(7) = fPPoint2b->Z() - fPPoint1->Z();
    }

    //Fill the covariance
    if (fFixedMeasurementCovariance) return kTRUE;
    //the dirs
    GetThetaPhiCov( *fPDir1ThetaPhiCov, *fPDir1, *fPDir1Cov );
    GetThetaPhiCov( *fPDir2ThetaPhiCov, *fPDir2, *fPDir2Cov );
    fPMeasurementCov->SetSub( 0, (*fPDir1ThetaPhiCov) + (*fPDir2ThetaPhiCov) );
    
    //the points
    (*fPMeasurementCov)(2,2) = (*fPPoint1Cov)(0,0)+(*fPPoint2Cov)(0,0);
    (*fPMeasurementCov)(2,3) = (*fPPoint1Cov)(0,2)+(*fPPoint2Cov)(0,2);
    (*fPMeasurementCov)(3,2) = (*fPPoint1Cov)(2,0)+(*fPPoint2Cov)(2,0);
    (*fPMeasurementCov)(3,3) = (*fPPoint1Cov)(2,2)+(*fPPoint2Cov)(2,2);

    //if 3 track mode set, set also the stuff for 2nd TPC track
    if (f3TracksMode)
    {
        GetThetaPhiCov( *fPDir2bThetaPhiCov, *fPDir2b, *fPDir2bCov );
        fPMeasurementCov->SetSub( 4, (*fPDir1ThetaPhiCov) + (*fPDir2bThetaPhiCov) );
        (*fPMeasurementCov)(6,6) = (*fPPoint1Cov)(0,0)+(*fPPoint2bCov)(0,0);
        (*fPMeasurementCov)(6,7) = (*fPPoint1Cov)(0,2)+(*fPPoint2bCov)(0,2);
        (*fPMeasurementCov)(7,6) = (*fPPoint1Cov)(2,0)+(*fPPoint2bCov)(2,0);
        (*fPMeasurementCov)(7,7) = (*fPPoint1Cov)(2,2)+(*fPPoint2bCov)(2,2);
    }
    return kTRUE;
}

Bool_t AliRelAlignerKalman::FillMeasurementMatrix()
{
    //Calculate the system matrix for the Kalman filter
    //approximate the system using as reference the track in the first volume

    Double_t delta = 1.e-8;
    TVectorD z1( fNMeasurementParams );
    TVectorD z2( fNMeasurementParams );
    TVectorD x1( *fPX );
    TVectorD x2( *fPX );
    TMatrixD D( fNMeasurementParams, 1 );
    for ( Int_t i=0; i<fNSystemParams; i++ )
    {
        x1 = *fPX;
        x2 = *fPX;
        x1(i) -= delta;
        x2(i) += delta;
        if (!PredictMeasurement( z1, x1 )) return kFALSE;
        if (!PredictMeasurement( z2, x2 )) return kFALSE;
        for (Int_t j=0; j<fNMeasurementParams; j++ )
            D.GetMatrixArray()[j] = (z2.GetMatrixArray()[j]-z1.GetMatrixArray()[j])/(2.*delta);
        fPH->SetSub( 0, i, D );
    }
    return kTRUE;
}

Bool_t AliRelAlignerKalman::PredictMeasurement( TVectorD& pred, const TVectorD& state )
{
    // Implements a system model for the Kalman fit
    // pred is [dtheta, dphi, dx, dz, [dtheta', dphi', dx', dz'] ] - second part is for 2nd TPC track in cosmics
    // state is [psi,theta,phi,x,y,z,driftTPC,offsetTPC]

    TVector3 newPoint = *fPPoint1;
    TVector3 newDir = *fPDir1;
    TVector3 shift;
    //TPC drift correction
    newDir(2) *= (1.+state(6));
    newPoint(2) *= (1.+state(6));
    //TPC offset
    if (newPoint(2)>0) newPoint(2) += state(7);
    else newPoint(2) -= state(7);
    //rotate
    TMatrixD rotmat(3,3);
    RotMat( rotmat, state );
    newDir = rotmat * newDir;
    newPoint = rotmat * newPoint;
    //shift
    shift(0) = state(3);
    shift(1) = state(4);
    shift(2) = state(5);
    newPoint += shift;
    //Predicted measurement
    if (!IntersectionLinePlane( newPoint, newPoint, newDir, *fPRefPoint, *fPVec010 )) return kFALSE;
    pred(0) = newDir.Theta() - fPDir1->Theta();
    pred(1) = newDir.Phi() - fPDir1->Phi();
    pred(2) = newPoint(0) - fPPoint1->X();
    pred(3) = newPoint(2) - fPPoint1->Z();

    if (f3TracksMode)
    {
        //Do the same for second TPC track
        //should be the same as the first TPC track
        pred(4) = pred(0);
        pred(5) = pred(1);
        pred(6) = pred(2);
        pred(7) = pred(3);
    }
    return kTRUE;
}

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
    if (
            (TMath::Abs(xupdate(0)) > fOutRejSigmas*TMath::Sqrt((*P)(0,0))) ||
            (TMath::Abs(xupdate(1)) > fOutRejSigmas*TMath::Sqrt((*P)(1,1))) ||
            (TMath::Abs(xupdate(2)) > fOutRejSigmas*TMath::Sqrt((*P)(2,2))) ||
            (TMath::Abs(xupdate(3)) > fOutRejSigmas*TMath::Sqrt((*P)(3,3))) ||
            (TMath::Abs(xupdate(4)) > fOutRejSigmas*TMath::Sqrt((*P)(4,4))) ||
            (TMath::Abs(xupdate(5)) > fOutRejSigmas*TMath::Sqrt((*P)(5,5))) ||
            (TMath::Abs(xupdate(6)) > fOutRejSigmas*TMath::Sqrt((*P)(6,6))) ||
            (TMath::Abs(xupdate(7)) > fOutRejSigmas*TMath::Sqrt((*P)(7,7))) 
        )
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

void AliRelAlignerKalman::TMatrixDSym_from_TMatrixD( TMatrixDSym& matsym, const TMatrixD& mat )
{
    //Produce a valid symmetric matrix out of a TMatrixD

    //not very efficient, diagonals are computer twice, but who cares
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

void AliRelAlignerKalman::GetThetaPhiCov( TMatrixDSym& cov, const TVector3& vec, const TMatrixDSym& vecCov )
{
    //Given the covariance matrix of a vector in euclid.space
    //give cov matrix of the 2 spherical angles
    const Double_t delta = 1e-8;
    TVector3 vec1;
    TVector3 vec2;
    TMatrixD T(2,3);
    for (Int_t i=0; i<3; i++)
    {
        vec1 = vec;
        vec2 = vec;
        vec1(i) -= delta/2.;
        vec2(i) += delta/2.;
        T(0,i) = (vec2.Theta() - vec1.Theta())/delta;
        T(1,i) = (vec2.Phi() - vec1.Phi())/delta;
    }
    TMatrixD tmp( T, TMatrixD::kMult, vecCov );
    TMatrixD nscov( tmp, TMatrixD::kMultTranspose, T );
    cov(0,0) = nscov(0,0); cov(0,1) = nscov(0,1);
    cov(1,0) = nscov(0,1); cov(1,1) = nscov(1,1);
    return;
}

Bool_t AliRelAlignerKalman::IntersectionLinePlane( TVector3& intersection, const TVector3& linebase, const TVector3& linedir, const TVector3& planepoint, const TVector3& planenormal )
{
    //calculate the intersection point between a straight line and a plane
    //plane is defined with a point in it and a normal, line has a point and direction
    if (planenormal==TVector3(0,1,0))
    {
        if (linedir(1)==0) return kFALSE;
        Double_t t = ( planepoint(1) - linebase(1) ) / linedir(1);
        intersection = linebase + t*linedir;
        return kTRUE;
    }
    else
    {
        Double_t dirdotn = linedir.Dot(planenormal);
        if (dirdotn==0) return kFALSE;
        Double_t t = ( planepoint.Dot(planenormal) - linebase.Dot(planenormal) ) / dirdotn;
        intersection = linebase + t*linedir;
        return kTRUE;
    }
}

void AliRelAlignerKalman::Angles( TVectorD &angles, const TMatrixD &rotmat )
{
//Calculate the Cardan angles (psi,theta,phi) from rotation matrix
//b = R*a 
//
        angles(0) = TMath::ATan2( -rotmat(1,2), rotmat(2,2) );
        angles(1) = TMath::ASin( rotmat(0,2) );
        angles(2) = TMath::ATan2( -rotmat(0,1), rotmat(0,0) );
        return;
}

void AliRelAlignerKalman::PrintCorrelationMatrix()
{
    //Print the correlation matrix for the fitted parameters

    for ( Int_t j=0; j<fNSystemParams; j++ )
    {
        for ( Int_t i=0; i<fNSystemParams; i++ )
        {
            printf("%1.2f  ", (*fPXcov)(i,j)/TMath::Sqrt( (*fPXcov)(i,i) * (*fPXcov)(j,j) ) );
        }//for i
        printf("\n");
    }//for j
    printf("\n");
    return;
}

Bool_t AliRelAlignerKalman::FindCosmicInEvent( AliESDEvent* pEvent )
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
    
    Double_t magfield = pEvent->GetMagneticField();

    //select sane tracks
    for (Int_t itrack=0; itrack < ntracks; itrack++)
    {
        pTrack = pEvent->GetTrack(itrack);
        if (!pTrack) {cout<<"no track!"<<endl;continue;}
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

        pTrack->RelateToVertex( pVertex, magfield, 10000. );
        Double_t trackposition[3];
        pTrack->GetXYZ( trackposition );
        distanceFromVertexArr[nGoodTracks] = 
              TMath::Sqrt((trackposition[0]-vertexposition[0])*(trackposition[0]-vertexposition[0])
                       + (trackposition[1]-vertexposition[1])*(trackposition[1]-vertexposition[1])
                       + (trackposition[2]-vertexposition[2])*(trackposition[2]-vertexposition[2]));

        //check if track is ITS or TPC
        if ( ((pTrack->GetStatus()&AliESDtrack::kITSin)>0)&&
             !((pTrack->GetStatus()&AliESDtrack::kTPCin)>0) )
        {
            ITStracksArr[nITStracks] = nGoodTracks;
            nITStracks++;
        }

        if ( ((pTrack->GetStatus()&AliESDtrack::kTPCin)>0)&&
             !((pTrack->GetStatus()&AliESDtrack::kITSin)>0) )
        {
            TPCtracksArr[nTPCtracks] = nGoodTracks;
            nTPCtracks++;
        }

        nGoodTracks++;
    }//for itrack   -   sanity cuts

    //printf("there are good tracks, %d in ITS and %d TPC.\n", nITStracks, nTPCtracks);
    
    if( nITStracks < 2 || nTPCtracks < 2 )
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

    //find matching in TPC
    Float_t min = 10000000.;
    Int_t TPCgood1 = -1, TPCgood2 = -1;
    for(Int_t itr1=0; itr1<nTPCtracks; itr1++)
    {
        for(Int_t itr2=itr1+1; itr2<nTPCtracks; itr2++)
        {
            Float_t deltatheta = TMath::Abs(TMath::Pi()-thetaArr[TPCtracksArr[itr1]]-thetaArr[TPCtracksArr[itr2]]);
            if(deltatheta > fMaxMatchingAngle) continue;
            Float_t deltaphi = TMath::Abs(TMath::Abs(phiArr[TPCtracksArr[itr1]]-phiArr[TPCtracksArr[itr2]])-TMath::Pi());
            if(deltaphi > fMaxMatchingAngle) continue;
            //printf("TPC: %f  %f     %f  %f\n",deltaphi,deltatheta,thetaArr[TPCtracksArr[itr1]],thetaArr[TPCtracksArr[itr2]]);
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
        return kFALSE;
    }

    //in ITS find 2 tracks that best match the found TPC cosmic
    Double_t min1 = 10000000.;
    Double_t min2 = 10000000.;
    Int_t ITSgood1 = -1, ITSgood2 = -1;
    for(Int_t itr1=0; itr1<nITStracks; itr1++)
    {
        if(distanceFromVertexArr[ITStracksArr[itr1]]>fMaxMatchingDistance) continue;
        Float_t deltatheta = TMath::Abs(TMath::Pi()-thetaArr[ITStracksArr[itr1]]-thetaArr[TPCgood1]);
        if(deltatheta > fMaxMatchingAngle) deltatheta = TMath::Abs( deltatheta - TMath::Pi() );
        if(deltatheta > fMaxMatchingAngle) continue;
        Float_t deltaphi = TMath::Abs(TMath::Abs(phiArr[ITStracksArr[itr1]]-phiArr[TPCgood1])-TMath::Pi());
        if(deltaphi > fMaxMatchingAngle) deltaphi = TMath::Abs( deltaphi - TMath::Pi() );
        if(deltaphi > fMaxMatchingAngle) continue;
        //printf("ITS %i dtheta, dphi, vrtx: %.4f, %.4f, %.4f\n",ITStracksArr[itr1], deltatheta, deltaphi,distanceFromVertexArr[ITStracksArr[itr1]]);
        if(deltatheta+deltaphi<min1) //only the best one
        {
            min1=deltatheta+deltaphi;
            //ITSgood2 = ITSgood1;
            ITSgood1 = ITStracksArr[itr1];
            continue;
        }
        if(deltatheta+deltaphi<min2) //the second best
        {
            min2=deltatheta+deltaphi;
            ITSgood2 = ITStracksArr[itr1];
            continue;
        }
    }
    if (ITSgood2 < 0) //no ITS cosmic track
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
    
    //////////////////////////////////////////////////////////////////////////////////////
    // convert indexes from local goodtrackarrays to global track index
    TPCgood1 = goodtracksArr[TPCgood1];
    TPCgood2 = goodtracksArr[TPCgood2];
    ITSgood1 = goodtracksArr[ITSgood1];
    ITSgood2 = goodtracksArr[ITSgood2];
    /////////////////////////////////////////////////////////////////////////////////////

    AliESDtrack * pTPCtrack1 = pEvent->GetTrack(TPCgood1);
    AliESDtrack * pTPCtrack2 = pEvent->GetTrack(TPCgood2);
    AliESDtrack * pITStrack1 = pEvent->GetTrack(ITSgood1);
    AliESDtrack * pITStrack2 = pEvent->GetTrack(ITSgood2);
    const AliTrackPointArray* pTPCArray1 = pTPCtrack1->GetTrackPointArray();
    const AliTrackPointArray* pTPCArray2 = pTPCtrack2->GetTrackPointArray();
    const AliTrackPointArray* pITSArray1 = pITStrack1->GetTrackPointArray();
    const AliTrackPointArray* pITSArray2 = pITStrack2->GetTrackPointArray();

    AliTrackPointArray* pfullITStrackArray = new AliTrackPointArray(1);
    JoinTrackArrays( pfullITStrackArray, pITSArray1, pITSArray2 );

    //cout<<"selected tracks: TPC: "<<TPCgood1<<" "<<TPCgood2<<" ITS: "<<ITSgood1<<" "<<ITSgood2<<endl;

    //Fill the final trackpointarrays
    if (!ExtractSubTrackPointArray( fPTrackPointArray1,  fPVolIDArray1,  pfullITStrackArray, 1, 6 )) return kFALSE;
    if (!ExtractSubTrackPointArray( fPTrackPointArray2,  fPVolIDArray2,  pTPCArray1,         7, 8 )) return kFALSE;
    if (!ExtractSubTrackPointArray( fPTrackPointArray2b, fPVolIDArray2b, pTPCArray2,         7, 8 )) return kFALSE;

    delete pfullITStrackArray; pfullITStrackArray=NULL;
    delete [] goodtracksArr; goodtracksArr=NULL;
    delete [] ITStracksArr; ITStracksArr=NULL;
    delete [] TPCtracksArr; TPCtracksArr=NULL;
    delete [] nPointsArr; nPointsArr=NULL;
    delete [] phiArr; phiArr=NULL;
    delete [] thetaArr; thetaArr=NULL;
    delete [] distanceFromVertexArr; distanceFromVertexArr=NULL;
    //printf("FindCosmicInEvent END\n");
    return kTRUE;
}

Bool_t AliRelAlignerKalman::ExtractSubTrackPointArray( AliTrackPointArray* pTrackOut, TArrayI* pVolIDArrayOut, const AliTrackPointArray* pTrackIn, const Int_t firstlayer, const Int_t lastlayer )
{
    //printf("ExtractSubTrackPointArray\n");
    //From pTrackIn select points between firstlayer and lastlayer and put the result in pTrackOut.
    //VolID's of selected points go in pVolIDArrayOut
    //only good points selected.
    //points can be ordered with z by setting fSortTrackPointsWithY

    AliTrackPoint point;
    Int_t nPointsIn = pTrackIn->GetNPoints();
    if (nPointsIn<TMath::Min(fMinPointsVol1,fMinPointsVol2)) return kFALSE;
    Int_t iPointArr[1000];
    Int_t VolIDArr[1000];
    Float_t yArr[1000];
    Int_t iOuter=-1;
    Int_t OuterLayerID=0;
    Int_t nGood=0;
    
    //loop over trackpoints and record the good ones
    for (Int_t i=0; i<nPointsIn; i++)
    {
        UShort_t VolIDshort = pTrackIn->GetVolumeID()[i];
        Int_t layerID = AliGeomManager::VolUIDToLayer(VolIDshort);
        //some points are very dirty: have nan's in coordinates or are otherwise sick
        {
        Float_t x = pTrackIn->GetX()[i];
        Float_t y = pTrackIn->GetY()[i];
        Float_t z = pTrackIn->GetZ()[i];
        if ( ((TMath::Abs(x)<1.)||(TMath::Abs(x)>1000.))&&
             ((TMath::Abs(y)<1.)||(TMath::Abs(y)>1000.))||
             ((TMath::Abs(z)>1000.))
           )
            continue;
        }
        if ( (layerID >= firstlayer) && (layerID <= lastlayer) )
        {
            
            iPointArr[nGood] = i;
            VolIDArr[nGood] = VolIDshort;
            yArr[nGood] = pTrackIn->GetY()[i];
            if (layerID>OuterLayerID)
            {
                OuterLayerID=layerID;
                iOuter = nGood;
            }
            nGood++;
        }
        else continue;
    }//for i=0..nPointsIn
    if (nGood<TMath::Min(fMinPointsVol1,fMinPointsVol2)) return kFALSE;
    //Fill the output array of VolID's
    delete pVolIDArrayOut;
    pVolIDArrayOut = new TArrayI( nGood );

    //fill the output trackarray
    Int_t* iSortedYArr = new Int_t[nGood];
    if (fSortTrackPointsWithY)
    {
        TMath::Sort( nGood, yArr, iSortedYArr, kFALSE );
    }
    delete pTrackOut;
    pTrackOut = new AliTrackPointArray( nGood );

    for ( Int_t i=0; i<nGood; i++)
    {
        pTrackIn->GetPoint( point, iPointArr[(fSortTrackPointsWithY)?iSortedYArr[i]:i] );
        pTrackOut->AddPoint( i, &point );
        pVolIDArrayOut->AddAt( VolIDArr[(fSortTrackPointsWithY)?iSortedYArr[i]:i], i );
    }
    //printf("ExtractSubTrackPointArray END\n");
    return kTRUE;
}

void AliRelAlignerKalman::SortTrackPointArrayWithY( AliTrackPointArray* pout, const AliTrackPointArray* pinn, const Bool_t descending )
{
    //sorts a given trackpointarray by y coordinate

    Int_t npoints = pinn->GetNPoints();
    const AliTrackPointArray* pin;
    if (pinn==pout)
        pin = new AliTrackPointArray(*pinn);
    else
    {
        pin = pinn;
        if (pout!=NULL) delete pout;
        pout = new AliTrackPointArray(npoints);
    }
    
    Int_t* iarr = new Int_t[npoints];

    TMath::Sort( npoints, pin->GetY(), iarr, descending );

    AliTrackPoint p;
    for (Int_t i=0; i<npoints; i++)
    {
        pin->GetPoint( p, iarr[i] );
        pout->AddPoint( i, &p );
    }
    delete [] iarr;
}

Bool_t AliRelAlignerKalman::FitTrackPointArrayRieman( TVector3* pPoint, TMatrixDSym* pPointCov, TVector3* pDir, TMatrixDSym* pDirCov, AliTrackPointArray* pTrackPointArray, TArrayI* pVolIDs, AliTrackPoint* pRefPoint )
{
    //printf("FitTrackPointArrayRieman\n");
    //Use AliTrackFitterRieman to fit a track through given point array
    
    AliTrackFitterRieman fitter;
    AliTrackPoint trackPoint;
    const Float_t* pointcov;
    ////here's an ugly hack://///////
    //AliTrackPointArray* pTrackPointArray = const_cast < AliTrackPointArray* > (pTrackPointArrayIn);
    /////////////////////////////////
    fitter.SetTrackPointArray( pTrackPointArray, kFALSE );
    if ( !fitter.Fit( pVolIDs ) ) return kFALSE;
    if ( !fitter.GetPCA( *fPRefAliTrackPoint, trackPoint ) ) return kFALSE;
    Double_t xPoint = trackPoint.GetX();
    pPoint->SetXYZ( xPoint, trackPoint.GetY(), trackPoint.GetZ() );
    pDir->SetXYZ( 1., fitter.GetDYat(xPoint), fitter.GetDZat(xPoint) );
    *pDir = pDir->Unit();
    //TODO cov of dir!
    
    pointcov = trackPoint.GetCov();
    (*pPointCov)(0,0) = pointcov[0];(*pPointCov)(0,1) = pointcov[1];(*pPointCov)(0,2) = pointcov[2];
    (*pPointCov)(1,0) = pointcov[1];(*pPointCov)(1,1) = pointcov[3];(*pPointCov)(1,2) = pointcov[4];
    (*pPointCov)(2,0) = pointcov[2];(*pPointCov)(2,1) = pointcov[4];(*pPointCov)(2,2) = pointcov[5];

    //printf("FitTrackPointArrayRieman END\n");
    return kTRUE;
}

Bool_t AliRelAlignerKalman::FitTrackPointArrayKalman( TVector3* pPoint, TMatrixDSym* pPointCov, TVector3* pDir, TMatrixDSym* pDirCov, AliTrackPointArray* pTrackPointArray, TArrayI* pVolIDs, AliTrackPoint* pRefPoint )
{
    //printf("FitTrackPointArrayKalman\n");
    //Use AliTrackFitterKalman to fit a track through given point array
    //still needs work
    
    AliTrackFitterKalman fitter;
    AliTrackPoint trackPoint, trackPoint2;

    pTrackPointArray->GetPoint( trackPoint, 0 );
    pTrackPointArray->GetPoint( trackPoint2, 1 );
    //Make sure all points count
    fitter.SetMaxChi2(100000000.);
    if (!fitter.MakeSeed( &trackPoint, &trackPoint2 )) return kFALSE;
    Int_t npoints = pTrackPointArray->GetNPoints();
    for (Int_t i=2; i<npoints; i++)
    {
        pTrackPointArray->GetPoint( trackPoint, i );
        if (!fitter.AddPoint( &trackPoint )) continue;
    }//for i
    TMatrixDSym fullcov = fitter.GetCovariance();
    printf("FitTrackPointArrayKalman: the covariance matrix of the fit");
    fullcov.Print();
    
    //now propagate to ref plane - its a hack that depends on the implementation of AddPoint()!
    // - first set the maxchi2 to 0 so addpoint returns false, but
    //actually the trackparameters have already been been propagated to the new planei, it's safe
    //because AliTrackFitterKalman::Propagate() is always kTRUE.
    fitter.SetMaxChi2(0.);
    fitter.AddPoint( fPRefAliTrackPoint );
        
    const Double_t* pparams = fitter.GetParam();
    fullcov = fitter.GetCovariance();

    pDir->SetXYZ(  pparams[3], 1., pparams[4] );
    pPoint->SetXYZ( pparams[0], pparams[1], pparams[2] );
    
    (*pPointCov)(0,0) = fullcov(0,0);(*pPointCov)(0,1) = fullcov(0,1);(*pPointCov)(0,2) = fullcov(0,2);
    (*pPointCov)(1,0) = fullcov(1,0);(*pPointCov)(1,1) = fullcov(1,1);(*pPointCov)(1,2) = fullcov(1,2);
    (*pPointCov)(2,0) = fullcov(2,0);(*pPointCov)(2,1) = fullcov(2,1);(*pPointCov)(2,2) = fullcov(2,2);

    (*pDirCov)(0,0)=fullcov(3,3); (*pDirCov)(0,1)=0.; (*pDirCov)(0,2)=fullcov(3,4);
    (*pDirCov)(1,0)=0.;           (*pDirCov)(1,1)=0.; (*pDirCov)(1,2)=0.;
    (*pDirCov)(2,0)=fullcov(4,3); (*pDirCov)(2,1)=0.; (*pDirCov)(2,2)=fullcov(4,4);
    return kTRUE;
}

void AliRelAlignerKalman::SanitizeExtendedCovMatrix( TMatrixDSym& covout, const TMatrixDSym& pointcov, const TVector3& dir )
{
    //reverse the effect of extending of the covariance matrix along the track
    //direction as done by the AliTrackFitters

    //idea is to rotate to a system where track direction is y axis,
    //in that refsystem setting the y related components of covariance
    //then rotate back.
    TVector3 yaxis(0,1,0);
    Double_t angle = dir.Angle(yaxis);
    TVector3 rotaxisvec = dir.Cross(yaxis);
    TVector3 rotaxis = rotaxisvec.Unit();
    TMatrixD R(3,3);

    //Fill the rotation matrix.
    Double_t ex=rotaxis(0);
    Double_t ey=rotaxis(1);
    Double_t ez=rotaxis(2);
    Double_t sin = TMath::Sin(angle);
    Double_t cos = TMath::Cos(angle);
    Double_t icos = 1 - cos;
    R(0,0) = cos + (ex*ex)*icos;
    R(0,1) = icos*ex*ey - ez*sin;
    R(0,2) = icos*ex*ez + ey*sin;
    R(1,0) = icos*ey*ex + ez*sin;
    R(1,1) = cos + (ey*ey)*icos;
    R(1,2) = icos*ey*ez - ex*sin;
    R(2,0) = icos*ez*ex - ey*sin;
    R(2,1) = icos*ez*ey + ex*sin;
    R(2,2) = cos + (ez*ez)*icos;

    //Transform covariance:
    TMatrixD covR( pointcov, TMatrixD::kMultTranspose, R);
    TMatrixD RcovR( R, TMatrixD::kMult, covR );
    TMatrixD newcov(3,3);
    newcov(0,0)=RcovR(0,0);
    newcov(0,2)=RcovR(0,2);
    newcov(2,0)=RcovR(2,0);
    newcov(2,2)=RcovR(2,2);
    newcov(1,1)=(newcov(0,0)+newcov(2,2))/2.;  //this is a bit insane!!
    covR = newcov *R;
    RcovR = R.T() * covR;
    TMatrixDSym_from_TMatrixD( covout, RcovR );
}

void AliRelAlignerKalman::JoinTrackArrays( AliTrackPointArray* pout, const AliTrackPointArray* pin1, const AliTrackPointArray* pin2 )
{
    //Join two AliTrackPointArrays

    Int_t np1 = pin1->GetNPoints();
    Int_t np2 = pin2->GetNPoints();
    if (pout!=NULL) delete pout;
    pout = new AliTrackPointArray(np1+np2);
    AliTrackPoint p;
    for (Int_t i=0;i<np1;i++)
    {
        pin1->GetPoint( p, i );
        pout->AddPoint( i, &p );
    }
    for (Int_t i=0;i<np2;i++)
    {
        pin2->GetPoint( p, i );
        pout->AddPoint( i+np1, &p );
    }
}

Bool_t AliRelAlignerKalman::ProcessTrackPointArrays()
{
    //Fit the track point arrays and update some household info
    
    fFitted=kFALSE; //not fitted yet
    if ( !FitTrackPointArrayKalman( fPPoint2, fPPoint2Cov, fPDir2, fPDir2Cov,
                                     fPTrackPointArray2, fPVolIDArray2, fPRefAliTrackPoint ) )
    {
        fNUnfitted++;
        return kFALSE;
    }

    if ( f3TracksMode )
    {
        if ( !FitTrackPointArrayKalman( fPPoint2b, fPPoint2bCov, fPDir2b, fPDir2bCov,
                                         fPTrackPointArray2b, fPVolIDArray2b, fPRefAliTrackPoint ) )
        {
            fNUnfitted++;
            return kFALSE;
        }
    }

    if ( !FitTrackPointArrayKalman( fPPoint1, fPPoint1Cov, fPDir1, fPDir1Cov,
                                     fPTrackPointArray1, fPVolIDArray1, fPRefAliTrackPoint ) )
    {
        fNUnfitted++;
        return kFALSE;
    }

    fFitted=kTRUE;
    //printf("ProcessTrackPointArrays: fitted!\n");
    return kTRUE;
}

Bool_t AliRelAlignerKalman::UpdateCalibration()
{
    //Update the calibration with new data, used in calibration pass

    fPThetaMesHist->Fill( (*fPMeasurement)(0) );
    fPPhiMesHist->Fill( (*fPMeasurement)(1) );
    fPXMesHist->Fill( (*fPMeasurement)(2) );
    fPZMesHist->Fill( (*fPMeasurement)(3) );
    if (f3TracksMode)
    {
        fPThetaMesHist2->Fill( (*fPMeasurement)(4) );
        fPPhiMesHist2->Fill( (*fPMeasurement)(5) );
        fPXMesHist2->Fill( (*fPMeasurement)(6) );
        fPZMesHist2->Fill( (*fPMeasurement)(7) );
    }
    fCalibrationUpdated=kTRUE;
    return kTRUE;
}

Bool_t AliRelAlignerKalman::SetCalibrationPass( Bool_t cp )
{
    //sets the calibration mode
    if (cp)
    {
        fCalibrationPass=kTRUE;
        return kTRUE;
    }//if (cp)
    else
    {
        if (!fCalibrationUpdated) return kFALSE;
        if (fCalibrationPass) // do it only after the calibration pass
        {
            Double_t tmp;
            fPMeasurementCov->Zero(); //reset the covariance

            fPThetaMesHist->Fit("gaus");
            TF1* fitformula = fPThetaMesHist->GetFunction("gaus");
            tmp = fitformula->GetParameter(2);
            (*fPMeasurementCov)(0,0) = tmp*tmp;

            fPPhiMesHist->Fit("gaus");
            fitformula = fPPhiMesHist->GetFunction("gaus");
            tmp = fitformula->GetParameter(2);
            (*fPMeasurementCov)(1,1) = tmp*tmp;

            fPXMesHist->Fit("gaus");
            fitformula = fPXMesHist->GetFunction("gaus");
            tmp = fitformula->GetParameter(2);
            (*fPMeasurementCov)(2,2) = tmp*tmp;

            fPZMesHist->Fit("gaus");
            fitformula = fPZMesHist->GetFunction("gaus");
            tmp = fitformula->GetParameter(2);
            (*fPMeasurementCov)(3,3) = tmp*tmp;

            if (f3TracksMode)
            {
                fPThetaMesHist2->Fit("gaus");
                fitformula = fPThetaMesHist2->GetFunction("gaus");
                tmp = fitformula->GetParameter(2);
                (*fPMeasurementCov)(4,4) = tmp*tmp;

                fPPhiMesHist2->Fit("gaus");
                fitformula = fPPhiMesHist2->GetFunction("gaus");
                tmp = fitformula->GetParameter(2);
                (*fPMeasurementCov)(5,5) = tmp*tmp;

                fPXMesHist2->Fit("gaus");
                fitformula = fPXMesHist2->GetFunction("gaus");
                tmp = fitformula->GetParameter(2);
                (*fPMeasurementCov)(6,6) = tmp*tmp;

                fPZMesHist2->Fit("gaus");
                fitformula = fPZMesHist2->GetFunction("gaus");
                tmp = fitformula->GetParameter(2);
                (*fPMeasurementCov)(7,7) = tmp*tmp;
            }

            fCalibrationPass=kFALSE;
            fFixedMeasurementCovariance=kTRUE;
            fPMeasurementCov->Print();
            return kTRUE;
        }//if (fCalibrationPass)
    return kFALSE;
    }//else (cp)
    return kFALSE;
}

void AliRelAlignerKalman::PrintDebugInfo()
{
    //prints debug info
    cout<<"AliRelAlignerKalman debug info"<<endl;
    printf("fPTrackPointArray1 y: ");
    for (Int_t i=0; i<fPTrackPointArray1->GetNPoints();i++)
    {
        printf("%.2f ",fPTrackPointArray1->GetY()[i]);
    }
    printf("\n");
    printf("fPTrackPointArray2 y: ");
    for (Int_t i=0; i<fPTrackPointArray2->GetNPoints();i++)
    {
        printf("%.2f ",fPTrackPointArray2->GetY()[i]);
    }
    printf("\n");
    printf("fPTrackPointArray2b y: ");
    for (Int_t i=0; i<fPTrackPointArray2b->GetNPoints();i++)
    {
        printf("%.2f ",fPTrackPointArray2b->GetY()[i]);
    }
    printf("\n");

    printf("Measurement covariance");
    fPMeasurementCov->Print();
}

