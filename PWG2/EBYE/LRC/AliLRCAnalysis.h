#ifndef ALILRCANALYSIS_H
#define ALILRCANALYSIS_H

//-------------------------------------------------------------------------
//   Description: 
//   This class is included into LRC library for Long-Range Correlation analysis
//   it is base class for NN, PtN, PtPt
//   implements base methods for thees classes
//   Origin: Petr Naumenko, SPbSU-CERN, Petr.Naoumenko@cern.ch
//-------------------------------------------------------------------------


/*  See cxx source for full Copyright notice */


/* $Id$ */

//-------------------------------------------------------------------------
//         LRC library for Long-Range Correlation analysis
//
//    Origin: Petr Naumenko, SPbSU-CERN, Petr.Naoumenko@cern.ch
//-------------------------------------------------------------------------

#define N_PL_FLAGS 10

#include "TObject.h"

class TFile;
class AliLRCFit;
class TProfile;
class TH1D;
class TH2D;
class TPaveText;
class TF1;
class math;
class TStyle;

class AliLRCAnalysis{
 public:
  void DrawAbs();
  void DrawAbs( int * mas );
  void DrawAbsPure( const int * const mas, bool drawPaveLabel );
  void DrawRel();
  void DrawRel( int * mas );
  void DrawRelPure( const int * const mas, bool drawPaveLabel );
  
  void SetXmin(double xMin);
  void SetXmax(double xMax);
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
  void Calculate();
  bool SetFitRange(double xMin, double xMax);
  void SetFullFitRange();
  
  double GetRoundWithError( double value, double error );
  double GetRoundWithError( double value, double error, int pres );
  double GetRoundWithPrecision( double value, int pres );
  
  double GetRoundValueErrorPrecision( double value, double error, int pres ) const;
  
  
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
  double HI2(TH1D * const h, double a, double b, double xmin, double xmax) const;
  double HI2(TH1D * const h, double a, double b) const;
  double Integral(TH2D* source, int nbin) const;
  //Creating profile from histogramm
  void CreateHist(char *name, char *nameAbs, char *nameRel, char *atitleF, char *atitleB,char *rtitleF, char *rtitleB,TH2D* sourceHist);
  void SetGraphics() const;
  void SetErrors(TH2D* source, const char *name);
  void SetErrors(TH2D* source, const char *name, double ptd, TH2D* nb);
  void SetErrors(TH2D* source, const char *name, double ptd, TProfile* nb);
  
  
 private:
  char*  fSx; 		// Title of x axis
  char*  fSy; 		// Title of y axis
  double fxFitMin; 	// FitMin minimum of fit baundary
  double fxFitMax; 	// FitMax maximum of fit baundary
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
  
  ClassDef(AliLRCAnalysis,0)                 // macro for rootcint
};

#endif

