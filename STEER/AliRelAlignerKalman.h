#ifndef ALIRELALIGNERKALMAN_H
#define ALIRELALIGNERKALMAN_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
//
//     Relative alignment of two tracking volumes (default ITS and TPC)
//     (see AliRelAlignerKalman.cxx for details)
//
//     Origin: Mikolaj Krzewicki, Nikhef, Mikolaj.Krzewicki@cern.ch
//
//////////////////////////////////////////////////////////////////////////////

#include <iostream>
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

class AliRelAlignerKalman {

public:
    AliRelAlignerKalman();
    virtual ~AliRelAlignerKalman() {}

    //User methods:
    Bool_t AddESDTrack( AliESDtrack* pTrack );
    Bool_t AddTrackPointArray( AliTrackPointArray* pTrackPointArr );
    Bool_t AddCosmicEventSeparateTracking( AliESDEvent* pEvent );
    
    void Print();
    void PrintCorrelationMatrix();

    void GetMeasurementCov( TMatrixDSym& cov ) { cov = *fPMeasurementCov; }
    void GetMeasurement( TVectorD& mes ) { mes = *fPMeasurement; }
    void GetState( TVectorD& state ) { state = *fPX; }
    void GetStateCov ( TMatrixDSym& cov ) { cov = *fPXcov; }
    void GetSeed( TVectorD& seed, TMatrixDSym& seedCov ) { seed = *fPX; seedCov = *fPXcov; }
    Double_t* GetStateArr() { return fPX->GetMatrixArray(); }
    Double_t* GetStateCovArr() { return fPXcov->GetMatrixArray(); }
    Double_t* GetMeasurementArr() { return fPMeasurement->GetMatrixArray(); }
    Double_t* GetMeasurementCovArr() { return fPMeasurementCov->GetMatrixArray(); }
    
    Bool_t SetCalibrationPass( Bool_t cp=kTRUE );

    //Experts only:
    void SetOutRejSigma( const Double_t a=2. ) { fOutRejSigmas = a; }
    void SetMeasurementCov( const TMatrixDSym& cov ) {*fPMeasurementCov = cov;}
    void SetMeasurement( const TVectorD& mes ) {*fPMeasurement = mes;}
    void SetState( const TVectorD& param ) {*fPX = param;}
    void SetStateCov (const TMatrixDSym& cov ) {*fPXcov = cov;}
    void SetSeed( const TVectorD& seed, const TMatrixDSym& seedCov ) {*fPX = seed; *fPXcov = seedCov; }
    void SetTrackParams1( const TVector3& point, const TMatrixDSym& pointcov, const TVector3& dir, const TMatrixDSym& dircov );
    void SetTrackParams2( const TVector3& point, const TMatrixDSym& pointcov, const TVector3& dir, const TMatrixDSym& dircov );
    void SetTrackParams2b( const TVector3& point, const TMatrixDSym& pointcov, const TVector3& dir, const TMatrixDSym& dircov );
    void SetRefSurface( const TVector3& point, const TVector3& dir );
    void SetRefSurface( const AliTrackPoint& point, const Bool_t horizontal=kTRUE );
    void SetRefSurface( const Double_t y );
    Bool_t PrepareUpdate();
    Bool_t Update();
    Bool_t IsReady() {return fFitted;}
    void SetQ( const Double_t Q = 1e-10 ) { fQ = Q; }
    void Set3TracksMode( Bool_t mode=kTRUE );
    void SetCosmicEvent( Bool_t cosmic=kTRUE ) {Set3TracksMode(cosmic);}
    void PrintDebugInfo();
    
protected:
    Bool_t UpdateEstimateKalman();
    Bool_t FillMeasurement();
    Bool_t FillMeasurementMatrix();
    Bool_t PredictMeasurement( TVectorD& z, const TVectorD& x );
    void GetThetaPhiCov( TMatrixDSym& cov, const TVector3& vec, const TMatrixDSym& vecCov );
    void Angles( TVectorD &angles, const TMatrixD &rotmat );
    void RotMat( TMatrixD& R, const TVectorD& angles );
    Bool_t IntersectionLinePlane( TVector3& intersection, const TVector3& base, const TVector3& dir, const TVector3& p, const TVector3& n );
    void TMatrixDSym_from_TMatrixD( TMatrixDSym& matsym, const TMatrixD& mat );
    void SanitizeExtendedCovMatrix( TMatrixDSym& covout, const TMatrixDSym& pointcov, const TVector3& dir );
    Bool_t ProcessTrackPointArrays();

    Bool_t ExtractSubTrackPointArray( AliTrackPointArray* pTrackOut, TArrayI* pVolIDArrayOut, const AliTrackPointArray* pTrackIn, const Int_t firstlayer, const Int_t lastlayer );
    Bool_t FitTrackPointArrayRieman( TVector3* pPoint, TMatrixDSym* pPointCov, TVector3* pDir, TMatrixDSym* pDirCov, AliTrackPointArray* pTrackPointArray, TArrayI* pVolIDs, AliTrackPoint* pRefPoint );
    Bool_t FitTrackPointArrayKalman( TVector3* pPoint, TMatrixDSym* pPointCov, TVector3* pDir, TMatrixDSym* pDirCov, AliTrackPointArray* pTrackPointArray, TArrayI* pVolIDs, AliTrackPoint* pRefPoint );
    Bool_t FindCosmicInEvent( AliESDEvent* pEvent );
    void SortTrackPointArrayWithY( AliTrackPointArray* pout, const AliTrackPointArray* pin, const Bool_t descending=kFALSE );
    void JoinTrackArrays( AliTrackPointArray* pout, const AliTrackPointArray* pin1, const AliTrackPointArray* pin2 );
    Bool_t UpdateCalibration();

private:
    AliTrackPointArray* fPTrackPointArray1;  //track point array in first volume (ITS)
    AliTrackPointArray* fPTrackPointArray2;  //track point array in second volume (TPC)
    AliTrackPointArray* fPTrackPointArray2b; //optionally second trackpoint array in second volume (TPC)
    TArrayI* fPVolIDArray1;                 //array with VolID's of the points in the trackpointarray 
    TArrayI* fPVolIDArray2;                 //array with VolID's of the points in the trackpointarray
    TArrayI* fPVolIDArray2b;                //array with VolID's of the points in the trackpointarray
    TVector3* fPDir1;  //track 1 direction
    TVector3* fPDir2;  //track 2 direction
    TVector3* fPDir2b; //track 2b direction (for cosmics)
    TMatrixDSym* fPDir1Cov;     //Covariance matrix of dir vec track 1
    TMatrixDSym* fPDir2Cov;     //Covariance matrix of dir vec track 2
    TMatrixDSym* fPDir2bCov;    //Covariance matrix of dir vec track 2b
    TMatrixDSym* fPDir1ThetaPhiCov; //Covariance matrix of the dir vector in spherical coord (angles only)
    TMatrixDSym* fPDir2ThetaPhiCov; //Covariance matrix of the dir vector in spherical coord (angles only)
    TMatrixDSym* fPDir2bThetaPhiCov; //Covariance matrix of the dir vector in spherical coord (angles only)
    TVector3* fPPoint1;  //track 1 extrapolated point at reference surface
    TVector3* fPPoint2;  //track 2 extrapolated point at reference surface
    TVector3* fPPoint2b; //track 2b extrapol. (for cosmics)
    TMatrixDSym* fPPoint1Cov; //covariance matrix of point 1
    TMatrixDSym* fPPoint2Cov; //covariance matrix of point 2
    TMatrixDSym* fPPoint2bCov; //covariance matrix of point 2b
    TVector3* fPRefPoint;    //a point on the reference surface
    TVector3* fPRefPointDir; //reference surfaces' direction vec
    AliTrackPoint* fPRefAliTrackPoint; //reference point as an AliTrackPoint (orientation stored in cov matrix)
    TMatrixD* fPRotMat; // system rotation matrix for estimated alignment angles
    
    //
    //Kalman filter related stuff
    //
    TVectorD* fPX; //System (fit) parameters (phi, theta, psi, x, y, z, driftcorr, driftoffset )
    TMatrixDSym* fPXcov; //covariance matrix of system parameters

    TMatrixD* fPH;      //System measurement matrix
    Double_t fQ;        //measure for system noise
    
    TVectorD* fPMeasurement; //the measurement vec for Kalman filter (theta,phi,x,z)
    TMatrixDSym* fPMeasurementCov; //measurement vec cvariance

    Int_t fNMeasurementParams; //how many numbers do we measure?
    Int_t fNSystemParams;      //how many parameters do we fit?

    Double_t fOutRejSigmas; //number of sigmas for outlier rejection
    //
    //
    //

    Bool_t f3TracksMode; //are we using 3 tracklets?
    Bool_t fSortTrackPointsWithY; //whether to sort the points after processing
    Bool_t fFixedMeasurementCovariance; //don't fiddle with measurement cov - supply it externally
    Bool_t fCalibrationPass;
    Bool_t fCalibrationUpdated;
    
    Bool_t fFitted; //did fitting finish without error?
    Bool_t fCuts;    //track cuts?

    Int_t fNTracks; //number of processed tracks
    Int_t fNUpdates; //number of successful Kalman update
    Int_t fNOutliers; //number of outliers
    Int_t fNUnfitted; //number of trackpointarrays that were not fitted for some reason
    Int_t fNMatchedCosmics; //number of matched tracklets for cosmics

    Int_t fMinPointsVol1;   //mininum number of points in volume 1
    Int_t fMinPointsVol2;   //mininum number of points in volume 2
    Double_t fMinMom;       //min momentum of track for track cuts
    Double_t fMaxMom;       //max momentum of track for track cuts
    Double_t fMinAbsSinPhi; //more cuts
    Double_t fMaxAbsSinPhi; //more cuts
    Double_t fMinSinTheta;  //cuts
    Double_t fMaxSinTheta;  //cuts
    Double_t fMaxMatchingAngle; //cuts
    Double_t fMaxMatchingDistance; //cuts
    Int_t fFirstLayerVol1; //first layer of volume 1
    Int_t fLastLayerVol1;  //last layer of volume 1
    Int_t fFirstLayerVol2; //first layer of volume 2
    Int_t fLastLayerVol2;  //last layer of volume 2

    TVector3 * fPVec010; //vector pointing up

    //Control and calibration histograms
    TH1D* fPXMesHist;  //histo of x measurement
    TH1D* fPZMesHist;  //histo of y measurement
    TH1D* fPPhiMesHist; //histo of phi measurement
    TH1D* fPThetaMesHist; //histo of theta measurement
    TH1D* fPXMesHist2;  //histo of x measurement (3tracks mode)
    TH1D* fPZMesHist2;  //histo of y measurement (3tracks mode)
    TH1D* fPPhiMesHist2; //histo of phi measurement (3tracks mode)
    TH1D* fPThetaMesHist2; //histo of theta measurement (3tracks mode)
    TH1D* fPMesCov11Hist;  //histogram of the covariance of a fit parameter, used in calibration
    TH1D* fPMesCov22Hist;  //histogram of the covariance of a fit parameter, used in calibration
    TH1D* fPMesCov33Hist;  //histogram of the covariance of a fit parameter, used in calibration
    TH1D* fPMesCov44Hist;  //histogram of the covariance of a fit parameter, used in calibration
    TH1D* fPMesCov55Hist;  //histogram of the covariance of a fit parameter, used in calibration
    TH1D* fPMesCov66Hist;  //histogram of the covariance of a fit parameter, used in calibration
    TH1D* fPMesCov77Hist;  //histogram of the covariance of a fit parameter, used in calibration
    TH1D* fPMesCov88Hist;  //histogram of the covariance of a fit parameter, used in calibration
    TMatrixDSym* fPMeasurementCovCorr; //correction to be applied to the measurement covariance
    
    static const Int_t fgkNMeasurementParams2TrackMode = 4; //how many measurables in 2 track mode
    static const Int_t fgkNMeasurementParams3TrackMode = 8; //how many measurables in 3 track mode
    static const Int_t fgkNSystemParams = 8;                //how many fit parameters
    
    AliRelAlignerKalman& operator= (const AliRelAlignerKalman& aligner );
    AliRelAlignerKalman(const AliRelAlignerKalman&);

    ClassDef(AliRelAlignerKalman,1)     //AliRelAlignerKalman class
};

#endif
