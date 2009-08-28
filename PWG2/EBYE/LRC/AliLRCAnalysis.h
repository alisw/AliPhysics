//-------------------------------------------------------------------------
//    Description: 
//    This class include into LRC library for Long-Range Correlation analysis
//    it is base class for NN, PtN, PtPt
//    implements base methods for thees classes
//    Origin: Petr Naumenko, SPbSU-CERN, Petr.Naoumenko@cern.ch
//-------------------------------------------------------------------------

#ifndef ALILRCANALYSIS_H
#define ALILRCANALYSIS_H

/*  See cxx source for full Copyright notice */


/* $Id$ */

//-------------------------------------------------------------------------
//         LRC library for Long-Range Correlation analysis
//
//    Origin: Petr Naumenko, SPbSU-CERN, Petr.Naoumenko@cern.ch
//-------------------------------------------------------------------------

#include "TFile.h"
#include "AliLRCFit.h"
#include "TProfile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TPaveText.h"
#include "TF1.h"
#include "math.h"
#include "TStyle.h"


class TH1D;
class TH2D;
class TFile;

class AliLRCAnalysis{
    public:
	void DrawAbs();
	void DrawRel();
	void SetXmin(double xMin);
	void SetXmax(double xMax);
	void SetBinsRange(int binMin, int binMax);
	double GetArel();
	double GetBrel();
	double GetXi2rel();
	double GetAabs();
	double GetBabs();
	double GetXi2abs();
	void Calculate();
	bool SetFitRange(double xMin, double xMax);
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
	double HI2(TH1D *h, double a, double b, double xmin, double xmax) const;
	double HI2(TH1D *h, double a, double b) const;
	double Integral(TH2D* source, int nbin) const;
	//Creating profile from histogramm
	void CreateHist(char *name, char *nameAbs, char *nameRel, char *atitleF, char *atitleB,char *rtitleF, char *rtitleB,TH2D* sourceHist);
	void SetGraphics() const;
	void SetErrors(TH2D* source, const char *name);
	void SetErrors(TH2D* source, const char *name, double ptd, TH2D* nb);
	void SetErrors(TH2D* source, const char *name, double ptd, TProfile* nb);
	

    private:
	char*  fSx; //Title of x axis
	char*  fSy; //Title of y axis
	double fxFitMin;
	double fxFitMax;
	double fa_rel;
	double fb_rel;
	double fXi2_rel;
	double fa_abs;
	double fb_abs;
	double fXi2_abs;


	ClassDef(AliLRCAnalysis,0)                 // macro for rootcint
};

#endif

