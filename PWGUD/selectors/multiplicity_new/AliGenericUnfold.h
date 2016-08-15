#ifndef ALIGENERICUNFOLD_H
#define ALIGENERICUNFOLD_H

#include <TH1.h>
#include <TH2.h>
#include <TVectorD.h>
#include <TMatrixD.h>
#include <Rtypes.h>
#include <TVectorT.h>
#include <TMatrixT.h>
#include <TGraph.h>
#include <TObject.h>
#include <TCanvas.h>

const char* optbayes = "Bayes";
const char* optchi2 = "Chi2";
const char* optsd = "TK";
const char* optsdw = "SDW";
const char* optsd1w = "SD1W";
const char* optsd2w = "SD2W";
const char* opttd = "TD";
const char* optsdp = "SDP";

const Double_t biasTreshold = 0.2;

class AliGenericUnfold : public TObject {

public:
	AliGenericUnfold ( Option_t* options = "SD Chi2" );
	virtual ~AliGenericUnfold();

	static AliGenericUnfold* GetInstance();

	//get-set
	void SetOptions ( Option_t* options );
	Option_t* GetOptions( );

	void SetWeight ( Double_t weight, Option_t* type );
	void SetWeight ( TVectorD* weights );
	Double_t GetWeight ( Option_t* type );

	void SetBayesLimit ( unsigned int limit );
	unsigned int GetBayesLimit();

	void SetData ( TVectorD* m, TVectorD* me, TMatrixD* r, Int_t stop = -1 );

	void SetDampMeasured ( Double_t threshold );

	void SetRegStart ( Int_t start );
	void SetmDampStart ( Int_t start );

	void ResetUnfolded ();

	void SetStrategy ( Int_t s );
	void SetPrint ( Int_t p );

	TH2D* GetC();
	TH2D* GetCR();

	//utility
	Double_t ChiSquared ( Double_t* params );
	TVectorD* Unfold();
	TVectorD* GetUnfoldErrorMinuit();
	TVectorD* Fold();
	TVectorD* GetResiduals();

	TMatrixD* CalculateCMatrix();
	TMatrixD* CalculateCMatrixDirect();
	TVectorD* GetUnfoldingBias();
	TMatrixD* GetUnfoldingVariance();
	TMatrixD* GetUnfoldingBiasVariance();
	TMatrixD* GetFullVariance(Bool_t minus = kTRUE);
	
	TMatrixD* GetBayesianVariance();
	
	Double_t GetDeltaChi2();
	Double_t GetDeltaBeta();
	
	Int_t GetLastBin();
	
	void SetUseMinuit2();

	void ScanWeight ( Double_t wlow, Double_t whi, Int_t nsteps, char* filename, Bool_t draw, Bool_t blog = kFALSE );
	void ScanWeightBayes ( Int_t nsteps, char* filename, Bool_t draw );


protected:
	TMatrixD* CalculateA();
	TMatrixD* CalculateB();
	TMatrixD* CalculateG();
	
	Double_t SecondDerivativeWeightedRegularization();
	
private:
	void LogMessage ( const char* message, Bool_t bcerr = kFALSE );

	static AliGenericUnfold* fInstance;

	//config
	TString fOptions;
	Bool_t fuseMinuit2;

	unsigned int fBayesIterationLimit;

	Double_t fBayesWeight;
	Double_t fChi2Weight;
	TVectorD* fChi2Weights;

	Int_t nX;
	Int_t nXstop;
	Int_t nY;

	Bool_t bUseMeasuredDamp;
	Double_t fDampThreshold;
	Int_t regStart;
	Int_t mDampStart;

	Int_t mStrategy;
	Int_t mPrint;
	
	//reg
	TMatrixD *G;
	TMatrixD *W;

	//data
	TVectorD* fMeasured;
	TVectorD* fMeasuredCopy;
	TVectorD* fMeasuredErr;
	TMatrixD* fMeasuredCov;
	TMatrixD* fMeasuredCovInv;
	TMatrixD* fResponse;
	TVectorD* fUnfolded;
	TVectorD* fUnfoldedErr;
	
	Double_t fMeasuredIntegral;

	TMatrixD* C;
	TMatrixD* A;
	TMatrixD* B;
// 	TMatrixD* G;
	TVectorD* bias;

	TMatrixD* Delta;

	Double_t deltaChi2;
	Double_t totalChi2;
	Int_t lastbin;
	
	TMatrixD* M; //Bayes migration matrix

	ClassDef ( AliGenericUnfold, 7 );

};

#endif // ALIGENERICUNFOLD_H
