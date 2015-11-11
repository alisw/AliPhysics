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
#include <TVectorD.h>
#include <TMatrix.h>

class TObject;
class AliExternalTrackParam;
class AliESDEvent;
class AliESDtrack;
class TArrayI;
class TObjArray;

class AliRelAlignerKalman : public TObject {

public:
    AliRelAlignerKalman();
    AliRelAlignerKalman(const char* name);
    virtual ~AliRelAlignerKalman();
    AliRelAlignerKalman& operator= (const AliRelAlignerKalman& a );
    AliRelAlignerKalman(const AliRelAlignerKalman& a);

    //User methods:
    Bool_t AddCosmicEvent( const AliESDEvent* pEvent );
    Bool_t AddTrackParams( const AliExternalTrackParam* p1, const AliExternalTrackParam* p2 );
    
    void Print(Option_t* option="") const;

    Double_t GetPsi() const {return (*fPX)(0);}
    Double_t GetTheta() const {return (*fPX)(1);}
    Double_t GetPhi() const {return (*fPX)(2);}
    Double_t GetX() const {return (*fPX)(3);}
    Double_t GetY() const {return (*fPX)(4);}
    Double_t GetZ() const {return (*fPX)(5);}
    Double_t GetTPCvdCorr() const {return (*fPX)(6);}
    Double_t GetTPCt0() const {return (*fPX)(7);}
    Double_t GetTPCvdY() const {if (fgkNSystemParams>8) return (*fPX)(8); else return 0.0;}
    Double_t GetPsiErr() const {return TMath::Sqrt((*fPXcov)(0,0));}
    Double_t GetThetaErr() const {return TMath::Sqrt((*fPXcov)(1,1));}
    Double_t GetPhiErr() const {return TMath::Sqrt((*fPXcov)(2,2));}
    Double_t GetXErr() const {return TMath::Sqrt((*fPXcov)(3,3));}
    Double_t GetYErr() const {return TMath::Sqrt((*fPXcov)(4,4));}
    Double_t GetZErr() const {return TMath::Sqrt((*fPXcov)(5,5));}
    Double_t GetTPCvdCorrErr() const {return TMath::Sqrt((*fPXcov)(6,6));}
    Double_t GetTPCt0Err() const {return TMath::Sqrt((*fPXcov)(7,7));}
    Double_t GetTPCvdYErr() const {if (fgkNSystemParams>8) return TMath::Sqrt((*fPXcov)(8,8)); else return 0.0;}
    void GetMeasurement( TVectorD& mes ) const { mes = *fPMeasurement; }
    TVectorD* GetMeasurement() { return fPMeasurement; }
    void GetMeasurementCov( TMatrixDSym& cov ) const { cov = *fPMeasurementCov; }
    void SetMeasurement( const TVectorD& mes ) {*fPMeasurement = mes;}
    void SetMeasurementCov( const TMatrixDSym& cov ) {*fPMeasurementCov = cov;}
    TMatrixDSym* GetMeasurementCov() const { return fPMeasurementCov; }
    void GetState( TVectorD& state ) const { state = *fPX; }
    TVectorD* GetState() const { return fPX; }
    void GetStateCov ( TMatrixDSym& cov ) const { cov = *fPXcov; }
    void SetState( const TVectorD& param ) {*fPX = param;}
    void SetStateCov (const TMatrixDSym& cov ) {*fPXcov = cov;}
    TMatrixDSym* GetStateCov() const { return fPXcov; }
    void GetSeed( TVectorD& seed, TMatrixDSym& seedCov ) const { seed = *fPX; seedCov = *fPXcov; }
    void SetSeed( const TVectorD& seed, const TMatrixDSym& seedCov ) {*fPX = seed; *fPXcov = seedCov; }
    Bool_t Merge( const AliRelAlignerKalman* al );
    Long64_t Merge( TCollection* list );

    //Expert methods:
    Bool_t AddESDevent( const AliESDEvent* pEvent );
    Bool_t AddESDtrack( const AliESDtrack* pTrack );
    void SetMagField( const Double_t f ) { fMagField=f; }
    Double_t GetMagField() const { return fMagField; }
    Bool_t FindCosmicTrackletNumbersInEvent( TArrayI& outITStracksTArr, TArrayI& outTPCtracksTArr, const AliESDEvent* pEvent );
    Int_t FindMatchingTracks(TObjArray& arrITS, TObjArray& arrTPC, AliESDEvent* pESD);
    Bool_t Update();
    void PrintCorrelationMatrix();
    //void PrintCovarianceCorrection();
    void PrintSystemMatrix();
  Int_t CheckCovariance(); // check covariance matrix
    void Reset();
    void ResetCovariance( const Double_t number=0. );
    void ResetTPCparamsCovariance( const Double_t number=0. );
    Double_t* GetStateArr() const { return fPX->GetMatrixArray(); }
    Double_t* GetStateCovArr() const { return fPXcov->GetMatrixArray(); }
    Double_t* GetMeasurementArr() const { return fPMeasurement->GetMatrixArray(); }
    Double_t* GetMeasurementCovArr() const { return fPMeasurementCov->GetMatrixArray(); }
    TMatrixD* GetH() const { return fPH; }
    TVectorD* GetMeasurementPrediction() const {return fPMeasurementPrediction;}
    const Double_t* GetDeltaArr() const {return fDelta;}
    void SetNumericalParanoia( const Bool_t mode=kFALSE ) { fNumericalParanoia=mode; }
    void SetCorrectionMode( const Bool_t mode=kTRUE ) { fCorrectionMode=mode; }
    void SetOutRejSigma( const Double_t a=2. ) { fOutRejSigmas = a; }
    void SetRejectOutliers( const Bool_t r=kTRUE ) {fRejectOutliers = r;}
    void SetRejectOutliersSigma2Median( const Bool_t b=kTRUE );
    void SetOutRejSigma2Median( const Double_t s ) {fOutRejSigma2Median = s;}
    Bool_t SetTrackParams( const AliExternalTrackParam* exparam1, const AliExternalTrackParam* exparam2 );
    const AliExternalTrackParam* GetTrackParams1() const {return fPTrackParam1;}
    const AliExternalTrackParam* GetTrackParams2() const {return fPTrackParam2;}
    void SetMinPointsVol1( const Int_t min ) {fMinPointsVol1=min;}
    void SetMinPointsVol2( const Int_t min ) {fMinPointsVol2=min;}
    void SetRequireMatchInTPC( const Bool_t s=kTRUE ) {fRequireMatchInTPC = s;}
    void SetQ( const Double_t Q = 1e-10 ) { fQ = Q; }
    void SetMaxMatchingDistance( const Double_t m ) {fMaxMatchingDistance=m;}
    void SetMaxMatchingAngle( const Double_t m ) {fMaxMatchingAngle=m;}
    void SetTPCvd( const Float_t v ) {fTPCvd=v;}
    void SetTPCZLengthA( const Double_t l ) {fTPCZLengthA=l;}
    void SetTPCZLengthC( const Double_t l ) {fTPCZLengthC=l;}
    Bool_t CorrectTrack( AliExternalTrackParam* tr, const TVectorD& misalignment ) const;
    Bool_t MisalignTrack( AliExternalTrackParam* tr, const TVectorD& misalignment ) const;
    static void Angles( TVectorD &angles, const TMatrixD &rotmat );
    static void RotMat( TMatrixD& R, const TVectorD& angles );
    static void TMatrixDSymFromTMatrixD( TMatrixDSym& matsym, const TMatrixD& mat );
    Bool_t IsPositiveDefinite( const TMatrixD& mat ) const;
    void SetTimeStamp( const UInt_t ts ) { fTimeStamp = ts; }
    UInt_t GetTimeStamp() const {return fTimeStamp;}
    void SetRunNumber( const Int_t rn ) { fRunNumber = rn; }
    Int_t GetRunNumber() const {return fRunNumber;}
    Int_t Compare(const TObject *obj) const;
    Bool_t IsSortable() const { return kTRUE; }
    Int_t GetNTracks() const {return fNTracks;}
    Int_t GetNUpdates() const {return fNUpdates;}
    Int_t GetNOutliers() const {return fNOutliers;}
    Int_t GetNOutliersSigma2Median() const {return fNOutliersSigma2Median;}
    Int_t GetNMerges() const {return fNMerges;}
    Int_t GetNMergesFailed() const {return fNMergesFailed;}
    void SetPoint2Track( Bool_t o );
    
protected:
    Bool_t UpdateEstimateKalman();
    Bool_t PrepareMeasurement();
    Bool_t PrepareSystemMatrix();
    Bool_t PreparePrediction();
    Bool_t PredictMeasurement( TVectorD& z, const TVectorD& x );
    Bool_t IsOutlier( const TVectorD& update, const TMatrixDSym& covmatrix );
    Bool_t IsOutlierSigma2Median( const AliExternalTrackParam* pITS, const AliExternalTrackParam* pTPC );

private:
    static const Int_t fgkNSystemParams=9;                //how many fit parameters
    static const Int_t fgkNtracksSigma2Median=500;        //how many sets for median and rms
    
    //Track parameters
    AliExternalTrackParam* fPTrackParam1;   //!local track parameters
    AliExternalTrackParam* fPTrackParam2;   //!local track parameters
    Double_t fMagField; //magnetic field

    //Kalman filter related stuff
    Int_t fNMeasurementParams;           //how many measurables
    TVectorD* fPX; //System (fit) parameters (phi, theta, psi, x, y, z, driftcorr, driftoffset )
    TMatrixDSym* fPXcov; //covariance matrix of system parameters
    TMatrixD* fPH;      //!System measurement matrix
    Double_t fQ;        //!measure for system noise
    TVectorD* fPMeasurement; //!the measurement vec for Kalman filter (theta,phi,x,z)
    TMatrixDSym* fPMeasurementCov; //!measurement vec cvariance
    TVectorD* fPMeasurementPrediction; //!prediction of the measurement
    Double_t fOutRejSigmas; //number of sigmas for outlier rejection
    Double_t fOutRejSigma2Median; //nsigmas to median of input residual distribution
    Double_t fDelta[fgkNSystemParams]; //array with differentials for calculating derivatives for every parameter(see PrepareSystemMatrix())
    Double_t* fResArrSigma2Median[4]; //!holds residuals for median based outlier removal

    //Control
    Bool_t fYZOnly;   //whether to consider only yz without directions.
    Bool_t fNumericalParanoia; //whether to perform additional checks for numerical stability
    Bool_t fRejectOutliers; //whether to do outlier rejection in the Kalman filter
    Bool_t fRejectOutliersSigma2Median; //whether to reject input based on distance to median
    Bool_t fRequireMatchInTPC;  //when looking for a cosmic in event, require that TPC has 2 matching segments
    Bool_t fCuts;    //track cuts?
    Int_t fMinPointsVol1;   //mininum number of points in volume 1
    Int_t fMinPointsVol2;   //mininum number of points in volume 2
    Double_t fMinPt;       //min momentum of track for track cuts
    Double_t fMaxPt;       //max momentum of track for track cuts
    Double_t fMaxMatchingAngle; //cuts
    Double_t fMaxMatchingDistance; //cuts
    Bool_t fCorrectionMode; //calculate corrective transform for TPC (or monitor actual TPC misal params)
    
    //Counters
    Int_t fNTracks; //number of processed tracks
    Int_t fNUpdates; //number of successful Kalman updates
    Int_t fNOutliers; //number of outliers
    Int_t fNOutliersSigma2Median; //number of rejected inputs
    Int_t fNMatchedCosmics; //number of cosmic events with matching tracklets (good cosmics)
    Int_t fNMatchedTPCtracklets;//number of cosmic events with 2 matching TPC tracklets
    Int_t fNProcessedEvents; //number of processed events
    UInt_t fTimeStamp;    //time stamp
    Int_t fRunNumber;    //run number
    Int_t fNMerges;      //how many succesful merges
    Int_t fNMergesFailed; //how many merges failed

    //TPC stuff
    Double_t fTPCvd; //TPC drift velocity
    Double_t fTPCZLengthA; //TPC length side A
    Double_t fTPCZLengthC; //TPC length side C
    
    ClassDef(AliRelAlignerKalman,4)     //AliRelAlignerKalman class
};

#endif


