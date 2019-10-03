//-------------------------------------------------------------------------
//    Description: 
//    This class is included into LRC library for Long-Range Correlation analysis
//    it is base class for NN, PtN, PtPt
//    implements base methods for thees classes
//    Origin: Petr Naumenko, SPbSU-CERN, Petr.Naoumenko@cern.ch,
//    Andrey Ivanov (SPbSU-CERN), Igor Altsebeev (SPbSU-CERN) 
//-------------------------------------------------------------------------

#ifndef ALILRCANALYSIS_H
#define ALILRCANALYSIS_H

/*  See cxx source for full Copyright notice */

#include "TObject.h"



class TFile;
class TProfile;
class TH1D;
class TH2D;


class AliLRCAnalysis :public  TObject{
    public:
	void DrawAbs();
	void DrawAbs(const int * const mas );
	void DrawAbsPure( const int * const mas, bool drawPaveLabel );
	void DrawRel();
	void DrawRel( const int * const mas );
	void DrawRelPure( const int * const mas, bool drawPaveLabel );
	void DrawHist( const int * const mDrawArray, bool drawPaveLabel, double aCoef, double bCoef
		, double aCoefError, double bCoefError, TH1D* profToDraw, int histType );

	void SetXmin(double xMin);
	void SetXmax(double xMax);
	void SetNsigma(double nSigma);
	void SetBinsRange(int binMin, int binMax);
	double GetArel() const;
	double GetBrel() const;
	double GetArelError() const;
	double GetBrelError() const;
	double GetXi2rel() const;
	double GetAabs() const;
	double GetBabs() const;
	double GetAabsError() const;
	double GetBabsError() const;
	double GetXi2abs() const;
	double GetFitXmin() const;
	double GetFitXmax() const;
	void Calculate();
	bool SetFitRange(double xMin, double xMax);
	bool SetFitRangeMin(double xMin);
	bool SetFitRangeMax(double xMax);
	void SetFullFitRange();
	void SetFitMethod(int id);

	double GetRoundWithError( double value, double error ) const ;
	double GetRoundWithError( double value, double error, int pres ) const;
	double GetRoundWithPrecision( double value, int pres ) const;
	
	double GetRoundValueErrorPrecision( double value, double error, int pres ) const;
	TH1D* GetAbsHisto() const;
	TH1D* GetRelHisto() const;
	TH1D* GetForwardValueDist() const;
	TH1D* GetBackwardValueDist() const;
	
	AliLRCAnalysis();
	AliLRCAnalysis(const AliLRCAnalysis& a);
	AliLRCAnalysis& operator= (const AliLRCAnalysis& a);
	virtual ~AliLRCAnalysis();

protected:
	TH1D* fPrAbs; //Work 1d histogramm in absolute var
	TH1D* fPrRel; //Work 1d histogramm in rellation var
	TH1D* fPrf; //Forward distribution
	TH1D* fPrb; //Backward distribution
	TFile* fileHist; // File with histrogramms
	double fdptb; //Work var for error calculation
	int fEntries; //Number of bins
	double Integral(TH2D* source, int nbin) const;
	//Creating profile from histogramm
	void CreateHist(char *name, char *nameAbs, char *nameRel, char *atitleF, char *atitleB,char *rtitleF, char *rtitleB,TH2D* sourceHist);
	void SetGraphics() const;
	void SetErrors(TH2D* source, const char *name);
	void SetErrors(TH2D* source, const char *name, double ptd, TH2D* nb);
	void SetErrors(TH2D* source, const char *name, double ptd, const TProfile* nb);

    private:
	static const int fgkPlotFlags = 10;  // Size of flags array used to chouse vat variables should be shown on a plot
	char*  fSx; 		// Title of x axis
	char*  fSy; 		// Title of y axis
	double fxFitMin; 	// FitMin minimum of fit baundary
	double fxFitMax; 	// FitMax maximum of fit baundary
	double fNsigma; 	// N sigma for fit range
	double farel; 		// ax = b the a relative
	double fbrel; 		// ax = b the b reletive
	double farelError; 	// a relative error
	double fbrelError; 	// b reletive error
	double fXi2rel; 	// chi square reletive
	double faabs; 		// ax = b the a absolut
	double fbabs; 		// ax = b the b absolut
	double faabsError; 	// a absolut error
	double fbabsError; 	// b absolut error
	double fXi2abs; 	// chi square absolut  
	int fFitMethod; 	// 0 - 1st variant, 1 - 2nd variant

	ClassDef(AliLRCAnalysis,0)                 // macro for rootcint
};

#endif

