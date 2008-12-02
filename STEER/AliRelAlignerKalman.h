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

#include <TMath.h>
#include <TObject.h>
#include <TVectorD.h>
#include <TMatrix.h>

class AliExternalTrackParam;
class AliESDEvent;
class AliESDtrack;

class AliRelAlignerKalman : public TObject {

public:
    AliRelAlignerKalman();
    virtual ~AliRelAlignerKalman();
    AliRelAlignerKalman& operator= (const AliRelAlignerKalman& a );
    AliRelAlignerKalman(const AliRelAlignerKalman& a);

    //User methods:
    Bool_t AddESDTrack(  const AliESDtrack* pTrack );
    Bool_t AddCosmicEvent( const AliESDEvent* pEvent );
    
    void Print(Option_t* option="") const;

    Double_t GetPsi() const {return (*fPX)(0);}
    Double_t GetTheta() const {return (*fPX)(1);}
    Double_t GetPhi() const {return (*fPX)(2);}
    Double_t GetX() const {return (*fPX)(3);}
    Double_t GetY() const {return (*fPX)(4);}
    Double_t GetZ() const {return (*fPX)(5);}
    Double_t GetTPCvdCorr() const {return (*fPX)(6);}
    Double_t GetTPCvdY() const {return (*fPX)(7);}
    Double_t GetTPCt0() const {return (*fPX)(8);}
    Double_t GetPsiErr() const {return TMath::Sqrt((*fPXcov)(0,0));}
    Double_t GetThetaErr() const {return TMath::Sqrt((*fPXcov)(1,1));}
    Double_t GetPhiErr() const {return TMath::Sqrt((*fPXcov)(2,2));}
    Double_t GetXErr() const {return TMath::Sqrt((*fPXcov)(3,3));}
    Double_t GetYErr() const {return TMath::Sqrt((*fPXcov)(4,4));}
    Double_t GetZErr() const {return TMath::Sqrt((*fPXcov)(5,5));}
    Double_t GetTPCvdCorrErr() const {return TMath::Sqrt((*fPXcov)(6,6));}
    Double_t GetTPCvdYErr() const {return TMath::Sqrt((*fPXcov)(7,7));}
    Double_t GetTPCt0Err() const {return TMath::Sqrt((*fPXcov)(8,8));}
    void GetMeasurement( TVectorD& mes ) const { mes = *fPMeasurement; }
    void GetMeasurementCov( TMatrixDSym& cov ) const { cov = *fPMeasurementCov; }
    void GetState( TVectorD& state ) const { state = *fPX; }
    void GetStateCov ( TMatrixDSym& cov ) const { cov = *fPXcov; }
    void GetSeed( TVectorD& seed, TMatrixDSym& seedCov ) const { seed = *fPX; seedCov = *fPXcov; }
    void SetMeasurement( const TVectorD& mes ) {*fPMeasurement = mes;}
    void SetMeasurementCov( const TMatrixDSym& cov ) {*fPMeasurementCov = cov;}
    void SetState( const TVectorD& param ) {*fPX = param;}
    void SetStateCov (const TMatrixDSym& cov ) {*fPXcov = cov;}
    void SetSeed( const TVectorD& seed, const TMatrixDSym& seedCov ) {*fPX = seed; *fPXcov = seedCov; }

    //Expert methods:
    Bool_t FindCosmicTrackletNumbersInEvent( TArrayI& outITStracksTArr, TArrayI& outTPCtracksTArr, const AliESDEvent* pEvent );
    Bool_t PrepareUpdate();
    Bool_t Update();
    void SetRefSurface( const Double_t x, const Double_t alpha );
    void PrintDebugInfo();
    void PrintCorrelationMatrix();
    void PrintCovarianceCorrection();
    void PrintSystemMatrix();
    void Reset();
    void ResetCovariance( const Double_t number=0. );
    void ResetTPCparamsCovariance( const Double_t number=0. );
    Double_t* GetStateArr() const { return fPX->GetMatrixArray(); }
    Double_t* GetStateCovArr() const { return fPXcov->GetMatrixArray(); }
    Double_t* GetMeasurementArr() const { return fPMeasurement->GetMatrixArray(); }
    Double_t* GetMeasurementCovArr() const { return fPMeasurementCov->GetMatrixArray(); }
    const Double_t* GetDeltaArr() const {return fDelta;}
    TH1D* GetMes0Hist() const {return fPMes0Hist;}
    TH1D* GetMes1Hist() const {return fPMes1Hist;} 
    TH1D* GetMes2Hist() const {return fPMes2Hist;}
    TH1D* GetMes3Hist() const {return fPMes3Hist;}
    TH1D* GetMesErr0Hist() const {return fPMesErr0Hist;}
    TH1D* GetMesErr1Hist() const {return fPMesErr1Hist;}
    TH1D* GetMesErr2Hist() const {return fPMesErr2Hist;}
    TH1D* GetMesErr3Hist() const {return fPMesErr3Hist;}
    Bool_t SetCalibrationMode( const Bool_t cp=kTRUE );
    void SetCorrectionMode( const Bool_t mode=kTRUE ){ fCorrectionMode=mode; }
    void SetApplyCovarianceCorrection( const Bool_t s=kTRUE ) {fApplyCovarianceCorrection = s;}
    void SetOutRejSigma( const Double_t a=2. ) { fOutRejSigmas = a; }
    void SetRejectOutliers( const Bool_t r=kTRUE ) {fRejectOutliers = r;}
    void SetCovarianceCorrection( const TMatrixDSym& c ) {*fPMeasurementCovCorr = c;}
    void GetCovarianceCorrection( TMatrixDSym& c ) {c=*fPMeasurementCovCorr;}
    void SetTrackParams1( const AliExternalTrackParam* exparam );
    void SetTrackParams2( const AliExternalTrackParam* exparam );
    void SetMinPointsVol1( const Int_t min ) {fMinPointsVol1=min;}
    void SetMinPointsVol2( const Int_t min ) {fMinPointsVol2=min;}
    void SetRequireMatchInTPC( const Bool_t s=kTRUE ) {fRequireMatchInTPC = s;}
    void SetNHistogramBins( const Int_t n ) {fNHistogramBins = n;}
    void SetQ( const Double_t Q = 1e-10 ) { fQ = Q; }
    void SetMaxMatchingDistance( const Double_t m ) {fMaxMatchingDistance=m;}
    void SetMaxMatchingAngle( const Double_t m ) {fMaxMatchingAngle=m;}
    void SetTPCvd( const Float_t v ) {fTPCvd=v;}
    void SetTPCZLengthA( const Double_t l ) {fTPCZLengthA=l;}
    void SetTPCZLengthC( const Double_t l ) {fTPCZLengthC=l;}
    Bool_t CorrectTrack( AliExternalTrackParam* tr, const TVectorD& misalignment );
    Bool_t MisalignTrack( AliExternalTrackParam* tr, const TVectorD& misalignment );
    static void Angles( TVectorD &angles, const TMatrixD &rotmat );
    static void RotMat( TMatrixD& R, const TVectorD& angles );
    static void TMatrixDSymFromTMatrixD( TMatrixDSym& matsym, const TMatrixD& mat );
    
protected:
    Bool_t UpdateCalibration();
    Bool_t UpdateEstimateKalman();
    Bool_t PrepareMeasurement();
    Bool_t PrepareSystemMatrix();
    Bool_t PredictMeasurement( TVectorD& z, const TVectorD& x );
    Bool_t IsOutlier( const TVectorD& update, const TMatrixDSym& covmatrix );
    Bool_t CalculateCovarianceCorrection();

private:
    static const Int_t fgkNMeasurementParams = 4;           //how many measurables
    static const Int_t fgkNSystemParams = 9;                //how many fit parameters
    
    //Track parameters
    Double_t fAlpha;       //rotation angle between the local and global coordinate system like in AliExternalTrackParam
    Double_t fLocalX;      //local x coordinate of reference plane = r(global)
    const AliExternalTrackParam* fkPTrackParam1;   //local track parameters (theta,phi,y,z)
    const AliExternalTrackParam* fkPTrackParam2;   //local track parameters

    //Kalman filter related stuff
    TVectorD* fPX; //System (fit) parameters (phi, theta, psi, x, y, z, driftcorr, driftoffset )
    TMatrixDSym* fPXcov; //covariance matrix of system parameters
    TMatrixD* fPH;      //System measurement matrix
    Double_t fQ;        //measure for system noise
    TVectorD* fPMeasurement; //the measurement vec for Kalman filter (theta,phi,x,z)
    TMatrixDSym* fPMeasurementCov; //measurement vec cvariance
    Double_t fOutRejSigmas; //number of sigmas for outlier rejection
    Double_t fDelta[fgkNSystemParams]; //array with differentials for calculating derivatives for every parameter(see PrepareSystemMatrix())

    //Control
    Bool_t fRejectOutliers; //whether to do outlier rejection in the Kalman filter
    Bool_t fCalibrationMode;            //are we running in calibration mode?
    Bool_t fFillHistograms;     //whether to fill the histograms with residuals
    Bool_t fRequireMatchInTPC;  //when looking for a cosmic in event, require that TPC has 2 matching segments
    Bool_t fApplyCovarianceCorrection;         //add the correction to the covariance of measurement
    Bool_t fCuts;    //track cuts?
    Int_t fMinPointsVol1;   //mininum number of points in volume 1
    Int_t fMinPointsVol2;   //mininum number of points in volume 2
    Double_t fMinMom;       //min momentum of track for track cuts
    Double_t fMaxMom;       //max momentum of track for track cuts
    Double_t fMaxMatchingAngle; //cuts
    Double_t fMaxMatchingDistance; //cuts
    Bool_t fCorrectionMode; //calculate corrective transform for TPC (or monitor actual TPC misal params)
    
    //Counters
    Int_t fNTracks; //number of processed tracks
    Int_t fNUpdates; //number of successful Kalman updates
    Int_t fNOutliers; //number of outliers
    Int_t fNMatchedCosmics; //number of cosmic events with matching tracklets (good cosmics)
    Int_t fNMatchedTPCtracklets;//number of cosmic events with 2 matching TPC tracklets
    Int_t fNProcessedEvents; //number of processed events

    //Calibration histograms
    Int_t fNHistogramBins;   //how many bins in control histograms
    TH1D* fPMes0Hist;  //histo of x measurement
    TH1D* fPMes1Hist;  //histo of y measurement
    TH1D* fPMes2Hist; //histo of phi measurement
    TH1D* fPMes3Hist; //histo of theta measurement
    TH1D* fPMesErr0Hist;  //histogram of the covariance of a fit parameter, used in calibration
    TH1D* fPMesErr1Hist;  //histogram of the covariance of a fit parameter, used in calibration
    TH1D* fPMesErr2Hist;  //histogram of the covariance of a fit parameter, used in calibration
    TH1D* fPMesErr3Hist;  //histogram of the covariance of a fit parameter, used in calibration
    TMatrixDSym* fPMeasurementCovCorr; //correction to be added to the measurement covariance

    //TPC stuff
    Double_t fTPCvd; //TPC drift velocity
    Double_t fTPCZLengthA; //TPC length side A
    Double_t fTPCZLengthC; //TPC length side C
    
    ClassDef(AliRelAlignerKalman,2)     //AliRelAlignerKalman class
};

#endif

