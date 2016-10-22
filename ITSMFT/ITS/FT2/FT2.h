#ifndef FT2_H
#define FT2_H
#include <TObject.h>
#include <TBits.h>
#include "AliExternalTrackParam.h"
#include "AliITSMFTAux.h"
#include "AliESDpid.h"
#include "AliPIDResponse.h"
#include "AliITSURecoLayer.h"
#include "AliNDLocalRegression.h"
#include "AliTPCParam.h"
#include "AliTPCRecoParam.h"
#include <TGenPhaseSpace.h>

//#define DEBUG 11
//#define DEBUG 1

class TTreeSRedirector;
class AliITSUReconstructor;
class AliITSURecoDet;
class AliITSURecoSens;
class TParticle;
class AliVertex;
class TF3;

const float kRidiculous = 999999;
const float kMaxZ = 250.;
const double kTrackToler = 1.e-6; // tracking tolerance

class FTProbe : public AliExternalTrackParam
{
public:
	FTProbe();
	virtual ~FTProbe(){};
	virtual Double_t GetTPCmomentum() const {return fTPCmomentum;}
	virtual Double_t GetTPCsignal() const {return fTPCSignal;}
	virtual UShort_t GetTPCsignalN() const {return fTPCSignalN;}
	virtual Int_t GetTPCTrackingPID() const {return fAbsPdgCodeForTracking;}
	virtual Double_t GetDecayInfo() const {return fIsDecayed;}
	virtual Double_t GetAbsorbtionInfo() const {return fIsAbsorbed;}
	virtual void GetInnerTrackParam(Double_t iTP[7]) const {
		for(Int_t i=0;i<7;i++){
			iTP[i]=fInnerTrackParameters[i];
		}
	}
	//
	Double_t GetElossNew() {return ft2ProbeElossNew;}
	Double_t GetElossOld() {return ft2ProbeElossOld;}
protected:
	Double_t fProbeMass;				// true mass
	Double_t fTPCmomentum;				// momentum after TPC reconstruction
	Double_t fTPCSignal;				// TPC signal
	UShort_t fTPCSignalN;
	Int_t fAbsPdgCode;				// |pdg code| of particle
	Int_t fPdgCode;					// pdg code of particle
	Int_t fAbsPdgCodeForTracking;		// |pdg code| used for tracking of particle
	Double_t fTrueMass;					// true mass of the particle
	Double_t fInnerTrackParameters[7];	// (fAlpha,fX,fP[5]) at inner TPC radius
	
	Double_t fProbeITSClusterIsCls[8];		// (ITScluster?,x,y,z)
	Double_t fProbeITSClusterX[8];				// (ITScluster?,x,y,z)
	Double_t fProbeITSClusterY[8];				// (ITScluster?,x,y,z)
	Double_t fProbeITSClusterZ[8];				// (ITScluster?,x,y,z)
	Double_t fProbeTPCClusterIsCls[160];	// (TPCcluster?,x,y,z)
	Double_t fProbeTPCClusterX[160];			// (TPCcluster?,x,y,z)
	Double_t fProbeTPCClusterY[160];			// (TPCcluster?,x,y,z)
	Double_t fProbeTPCClusterZ[160];			// (TPCcluster?,x,y,z)
	
	Double_t ft2ProbeElossOld;
	Double_t ft2ProbeElossNew;
	Int_t		fProbeNClTPC;			// N used TPC clusters
	Int_t		fProbeNClITS;			// N used ITS clusters
	Int_t		fProbeNClITSFakes;		// N used ITS Fake clusters
	Int_t		fProbeITSPatternFake;	// fakes pattern for ITS
	Int_t		fProbeITSPattern;		// pattern for ITS clusters
	Double_t	fProbeChi2TPC;			// total chi2 in TPC
	Double_t	fProbeChi2ITS;						// total chi2 in ITS
	Bool_t fIsDecayed;									// is particle decayed?
	Bool_t fIsAbsorbed;									// is particle absorbed?
	Double_t fDecayRadius;								// radius when particle decayed
	Double_t fAbsorbtionRadius;						// radius when particle was absorbed
	Double_t fAbsorbtionX;						// X when particle was absorbed
	Double_t fAbsorbtionY;						// Y when particle was absorbed
	Double_t fDecayX;						// X when particle decayed
	Double_t fDecayY;						// Y when particle decayed
	Bool_t fLostInItsTpcMatching;					// was track lost due to ITS-TPC matching efficiency?
	Int_t fTrackToClusterChi2CutITS;		// is particle rejected by ITS track to cluster chi2?
	Double_t chiwITS[7];								// chi2 for each layer of the ITS
	Double_t fProbeZAtCutOffCheck;	// z coordinate when the probe is check for z
	Int_t fProbeReconstructionError; //
	ClassDef(FTProbe,1)
};


class FT2 : public TObject
{
public:
	enum {kMaxITSLr=7, kMaxHitPerLr=2};
	struct FT2TPCLayer {
  FT2TPCLayer(int id=-1,float xr=0,float x2x=0, float pitchL =0, float pitchA =0) :
		rowId(id),x(xr),x2x0(x2x),isEdge(0),pitch(pitchL),pitchAng(pitchA),
		hitY(0),hitZ(0),hitSect(-1) {}
		Int_t    rowId;
		Float_t  x;
		Float_t  x2x0;
		Bool_t   isEdge;
		Float_t  pitch;
		Float_t  pitchAng;
		//
		Float_t  hitY;
		Float_t  hitZ;
		Int_t  hitSect;
		//
		void GetXYZLab(Float_t xyz[3]) const
		{
			float phi = (hitSect*20+10)*TMath::DegToRad(), cs=TMath::Cos(phi), sn=TMath::Sin(phi);
			xyz[0] = x*cs - hitY*sn; xyz[1] = x*sn + hitY*cs; xyz[2] = hitZ;
		}
	};
	typedef struct FT2TPCLayer FT2TPCLayer_t;
	
	
	FT2();
	virtual ~FT2();
	
	void InitEnvLocal();
	void InitTPCParaFile(const char *TPCparaFile);
	void InitXSectionFile(const char *XSectionFile);
	void InitTPCPIDResponse();
	void InitDetector(Bool_t addTPC=kTRUE,Float_t scEdge=2.0); //used to be 2.6; 2.0 determined by Marian and Johannes from 2013 data
	//
	void PrintLayout();
	//
	Bool_t ProcessTrack(TParticle* trc, AliVertex* vtx);
	void   SetSimMat(Bool_t v=kTRUE) {fSimMat = v;}
	void   SetTuneOnDataOrMC(Bool_t t=kFALSE) {fTuneOnDataOrMC = t;}
	void   SetAllowDecay(Bool_t d=kTRUE) {fAllowDecay = d;}
	void   SetAllowAbsorbtion(Bool_t f=kFALSE) {fAllowAbsorbtion = f;}
	void   SetUseConversionExtension(Bool_t z=kFALSE) {fUseConverisons = z;}
	void   SetMinTPCHits(int v=0)  {fMinTPCHits=v;} // 60
	void   SetMinITSLrHit(int v=0)  {fMinITSLrHit=v;}
	void   SetUsePIDForTracking(Bool_t usePID) {fUsePIDForTracking=usePID;}
	void	 SetRunNumber( Int_t runnumber) {fRunNumber = runnumber;}
	void   SetStandaloneTune( Int_t stdTune) {fStandaloneTune = stdTune;}
	void   SetStreamLevel(Int_t level) {fStreamLevel = level;}
	//Int_t  ProbeDecay(double* posIni);
	Int_t		ProbeDecay(double* posIni,double *params);
	Int_t		ProbeAbsorbtion(double* posIni,double *params);
	Bool_t	Pion2Muon();
	//
	Int_t    GetNClITSFakes()  const {return fNClITSFakes;}
	Int_t    GetNClITS()  const {return fNClITS;}
	Int_t    GetNClTPC()  const {return fNClTPC;}
	Double_t GetChi2ITS() const {return fChi2ITS;}
	Double_t GetChi2TPC() const {return fChi2TPC;}
	const TBits&   GetTPCHitMap() const {return fTPCMap;}
	FTProbe& GetProbe() const {return (FTProbe&)fProbe;}
	AliExternalTrackParam& GetKalmanOut(int i) {return (AliExternalTrackParam&)fKalmanOutward[i];}
	void   SetUseKalmanOut(Bool_t v=kTRUE)  {fUseKalmanOut = v;}
	Bool_t GetUseKalmanOut()   const {return fUseKalmanOut;}
	//
	const Double_t* GetDCA()    const {return &fDCA[0];}
	const Double_t* GetDCACov() const {return &fDCACov[0];}
	static float GetMaxStepTGeo() {return fgMaxStepTGeo;}
	static void  SetMaxStepTGeo(float stp=1.) {fgMaxStepTGeo = stp;}
	//
	void     SetAllocCorrelatedITSFakes(Bool_t v=kTRUE) {fAllocCorrelatedITSFakes=v;}
	Bool_t   GetAllocCorrelatedITSFakes()               {return fAllocCorrelatedITSFakes;}
	Double_t GetMagneticField()													{return fBz;}
	//
	Double_t GetNCorrelITSFakes(Int_t lr) const         {return fNCorrelITSFakes[lr];}
	Double_t GetCorrelITSFakesSigY(Int_t lr) const      {return fCorrelITSFakesSigY[lr];}
	Double_t GetCorrelITSFakesSigZ(Int_t lr) const      {return fCorrelITSFakesSigZ[lr];}
	//
	void     SetNCorrelITSFakes(Int_t lr, double v)     {fNCorrelITSFakes[lr] = v;}
	void     SetCorrelITSFakesSigY(Int_t lr, double v)  {fCorrelITSFakesSigY[lr] = v;}
	void     SetCorrelITSFakesSigZ(Int_t lr, double v)  {fCorrelITSFakesSigZ[lr] = v;}
	//
	void		SetMCTrueTrackMultiplicity(Int_t mult)			{fTrueMCtrackMult = mult;}
	//
	void     SetdNdY(double v=-1) {fdNdY = v;}
	Double_t GetdNdY() const {return fdNdY;}
	Double_t HitDensity(double r2, double tgl) const;
	Bool_t   BiasAsFake(double yz[2], const double* extyz, const double *cov) const;
	Bool_t DiagonalizeErrors(const double *cov, double &sy2d, double &sz2d);
	Int_t  GetITSPattern() const {return fITSPattern;}
	Int_t  GetITSPatternFakes() const {return fITSPatternFake;}
	//
	Double_t GetITSRMin() const;
	Double_t ParticleDecayProbability(Double_t step);
	Double_t ParticleAbsorptionProbability(Double_t length,Double_t rho, Double_t A, Double_t Z);
	//
	Int_t CutOnTrackToClusterChi2ITS();
	void AddCovariance();
	void AddCovarianceAdd();
	Int_t FindNextLayer(Double_t radius, Bool_t insideITS); // function to check radial origin of conversions
protected:
	void AddTPC(Float_t scEdge=2.6);
	void AddTPCLayer(Int_t rowID, Float_t x, Float_t x2x0, Float_t pitch, Float_t pitchAng);
	Bool_t InitProbe(TParticle* trc);
	Bool_t MakeITSKalmanOut();
	Bool_t PrepareProbe();
	Bool_t ApplyMSEloss(double x2X0, double xrho, double A, double Z, int dEcheck);
	Bool_t PropagateToX(double xTgt, int dir,Bool_t propErr,Bool_t simMat,Bool_t useTGeo,Int_t dECheckptx, Bool_t decAbsAllowed);
	Bool_t PropagateToR(double xTgt, int dir,Bool_t propErr,Bool_t simMat,Bool_t useTGeo, Bool_t decAbsAllowed);
	Bool_t IsZero(double val, double tol=1e-9) const {return TMath::Abs(val)<tol;}
	Bool_t PassActiveITSLayer(AliITSURecoLayer* lr);
	Bool_t GetRoadWidth(AliITSURecoLayer* lrA,double *pr,Int_t nstd = 3);
	void   ResetCovMat(AliExternalTrackParam* trc);
	Double_t UpdateKalman(AliExternalTrackParam* trc, double y,double z,double sigY,double sigZ,
						  Bool_t randomize=kTRUE,double scly=1.,double sclz=1.);
	const AliExternalTrackParam* GetSmoothedEstimate(int ilr,const AliExternalTrackParam* trcInw, double* trPos,double* trCov);
	Int_t  ReconstructOnITSLayer(int ilr, double chi2Cut=70.);
	
	Bool_t ReconstructProbe(TParticle *part);
	AliPIDResponse::EDetPidStatus GetComputeTPCProbability (const AliVTrack *track, Int_t nSpecies, Double_t p[]) const;
	//
protected:
	AliITSUReconstructor* fITSRec; // interface for reconstructor
	AliITSURecoDet* fITS;          // interface for ITS
	Bool_t fUsePIDForTracking;
	Int_t fRunNumber;
	Bool_t fStandaloneTune;
	Bool_t fIsTPC;                 // TPC added
	AliPIDResponse* fPIDResponse;
	TFile *fTPCParaFile;
	TFile *fXSectionFile;
	TH1F* fXSectionHp[6];
	Float_t fTPCSectorEdge;        // cut in cm on sector edge
	Double_t fMaxSnpTPC;           // stop particle if snp>fMaxSnpTPC
	Int_t fTPCInnerRows;
	Int_t fTPCMiddleRows;
	Int_t fTPCOuterRows;
	std::vector<FT2TPCLayer_t> fTPCLayers;
	std::vector<int>           fTPCHitLr;
	Int_t                      fNTPCHits; //!
	Int_t                      fMinTPCHits; // require min amount of TPC hits
	Int_t                      fMinITSLrHit; // require min amount of ITS lr hit
	Int_t											 fTrueMCtrackMult; // mc track multiplicity
	Double_t fDCA[2],fDCACov[3];   //! dca to vertex and its covariance
	//
	FTProbe fProbe;  // track
	AliExternalTrackParam fProbeIni;				//! initial probe kinematics
	AliExternalTrackParam* fKalmanOutward;	//! parameters of outward kalman
	Bool_t  fUseKalmanOut;									//! use KalmanOut estimate for fakes
	//Double_t              fProbeMass;			// probe mass
	Double_t              fBz;							// bz
	Bool_t								fTuneOnDataOrMC;	// used to tune some parameters on data or MC input; can only be used in second iteration
	Bool_t                fSimMat;					// simulate material effects in probe preparation
	Bool_t								fAllowDecay;			// necessary for standlone FT2 mode
	Bool_t								fAllowAbsorbtion; // necessary for standlone FT2 mode
	Bool_t								fUseConverisons;  // activates some extra checks needed for handling conversions 
	Int_t                 fCurrITSLr;				//! current ITS layer under tracking
	Int_t                 fNClTPC;					//! N used TPC clusters
	Int_t                 fNClITS;					//! N used ITS clusters
	Int_t                 fNClITSFakes;			//! N used ITS Fake clusters
  Int_t                 fITSPatternFake;	//! fakes pattern for ITS
	Int_t                 fITSPattern;			//! pattern for ITS clusters
	Double_t              fChi2TPC;					//! total chi2 in TPC
	Double_t              fChi2ITS;					//! total chi2 in ITS
	TBits                 fTPCMap;					//! tpc hit map
	//
	// hit info in the ITS
	Double_t fSigYITS,fSigZITS;       // nominal ITS layer resolution
	Double_t fITSerrSclY,fITSerrSclZ; // scaling factor for assigned errors
	Int_t fNITSHits[kMaxITSLr]; //! n hits per ITS layer
	Int_t fNITSSensCand[kMaxITSLr]; //! n sensor candidates per ITS layer
	Int_t fNITSLrHit;           //! n of ITS layers whith hit
	AliITSURecoSens* fITSSensHit[kMaxITSLr][2]; //! hit sensors
	AliITSURecoSens* fITSSensCand[kMaxITSLr][AliITSURecoLayer::kMaxSensMatching]; //! hit sensor candidates
	Double_t fITSHitYZ[kMaxITSLr][kMaxHitPerLr][2]; //! tracking Y,Z of each hit
	//
	Double_t fdNdY;             // if positive, use it for fakes simulation
	AliNDLocalRegression *fTPCClsLossProbIROC;					// parameterization of the TPC cluster pick up probability in the IROC
	AliNDLocalRegression *fTPCClsLossProbOROCmedium;		// parameterization of the TPC cluster pick up probability in the OROC (medium)
	AliNDLocalRegression *fTPCClsLossProbOROClong;			// parameterization of the TPC cluster pick up probability in the OROC (long)
	AliNDLocalRegression *fTpcPidSignal[5];							// parameterization of the TPC signal (e,mu,pi,K,p)
	AliNDLocalRegression *fTPCDistortionRPhi;						// parameterization of the TPC field distortions in rphi
	AliNDLocalRegression *fTPCDistortionR;							// parameterization of the TPC field distortions in r
	AliNDLocalRegression *fTPCfMCChi2;									// parameterization of the TPC standardized chi2
	AliNDLocalRegression *fTPCft2Chi2;									// parameterization of the FT2 standardized chi2
	AliNDLocalRegression *fItsTpcMatchingfMC;						// parameterization of the ITS+TPC matching efficiency in full MC
	AliNDLocalRegression *fItsTpcMatchingft2;						// parameterization of the ITS+TPC matching efficiency in FT2
	TGenPhaseSpace* fPionDecayGen;

	TParticle *fParticle;
	
  TTreeSRedirector *fDebugStreamer;     //!debug streamer
	Int_t fStreamLevel;
	TF1 *fTpcClusterAcc;
	TF1 *fTPCClusterErrorParamY; // parameterization of the TPC cluster error in y and z
	TF1 *fTPCClusterErrorParamZ; // parameterization of the TPC cluster error in y and z
	const Double_t *fTPCClusterErrInner; // systematic error of the cluster - used to downscale the information
	const Double_t *fTPCSystematicErr; // systematic errors in the track parameters - to be added to TPC covariance matrix
	const Double_t *fTPCClusterErr; // systematic error of the cluster - used e.g in OpenGG run to provide better cluster to track association efficiency
	Bool_t fTPCUseSystematicCorrelation; // switch to use the correlation for the sys
	AliTPCParam *fTPCParam; // TPC param
	AliTPCRecoParam *fTPCRecoParam; // TPC reco param
	Bool_t   fAllocCorrelatedITSFakes; // simulate noise hits accompanying the probe
	Double_t fNCorrelITSFakes[kMaxITSLr]; // av.number of accompanying hits
	Double_t fCorrelITSFakesSigY[kMaxITSLr]; // their width in Y
	Double_t fCorrelITSFakesSigZ[kMaxITSLr]; // their width in Z
	//
	Float_t fC0tr2clChi2[kMaxITSLr];	// cut on cluster to track chi2
	Float_t fC0gloChi2[kMaxITSLr];		// cut on seed global norm chi2
	Float_t fC0missPen[kMaxITSLr];		// missing cluster penalty
	//
	static float fgMaxStepTGeo; // max step for tracking accounting for TGeo materials
	
	ClassDef(FT2,1)
};

#endif
