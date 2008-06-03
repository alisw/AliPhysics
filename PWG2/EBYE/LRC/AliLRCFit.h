//-------------------------------------------------------------------------
//    Description: 
//    This class include into LRC library for Long-Range Correlation analysis
//    it makes fit of the 1d histogramm
//    calculates ax+b coefficients with error and hi square
//    Origin: Petr Naumenko, SPbSU-CERN, Petr.Naoumenko@cern.ch
//-------------------------------------------------------------------------

#ifndef ALILRCFIT_H
#define ALILRCFIT_H

/*  See cxx source for full Copyright notice */


/* $Id$ */

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "math.h"
#include "TProfile.h"
#include "TPaveText.h"
#include "TF1.h"
#include "TStyle.h"


class TH1D;
class TH2D;
class TFile;


class AliLRCFit{
    public:
    	AliLRCFit();
    	AliLRCFit(TH1D* h);
    	AliLRCFit(TH1D* h, double xmin, double xmax);
   	virtual ~AliLRCFit();
	double Geta() const;
	double Getb() const;
	double Getda() const;
	double Getdb() const;
	double Getda1() const;
	double Getdb1() const;
	double Gethi2() const;
	double Getf() const;
	double Getxmin() const;
	double Getxmax() const;
	int GetN() const;

    private:
    	int fN; //Number of bins
	int fTrueN; //Number of bins between xmin and xmax
    	int fNmin; //xmin bin
	int fNmax; //xmax bin 
    	double fS1; // work wariable
    	double fSz; // work wariable
    	double fSfz; // work wariable
    	double fSf; // work wariable
    	double fSf2; // work wariable
	double fhi2; // hi square
    	double fw; // work wariable
	double fz; // work wariable
    	double f; // work wariable
	double fnum; // work wariable
	double fdf; // work wariable
	double fdelta; // work wariable
	double fa; // coefficient a of ax+b
	double fb; // coefficient b of ax+b
	double fda; //error of coefficient a of ax+b
	double fdb; //error of coefficient b of ax+b
	double fda1;  // work wariable
	double fdb1; // work wariable
	double fxmin; // xmin
	double fxmax; // xmax
	ClassDef(AliLRCFit,0)                 // macro for rootcint

};

#endif

