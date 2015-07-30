#ifndef ALINDLOCALREGRESSION_H
#define ALINDLOCALREGRESSION_H
/* Copyright(c) 2006-07, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "TGraph.h"
#include "TGraphSmooth.h"
#include "TRandom.h"
#include "TSpline.h"
#include "TLinearFitter.h"
#include "TDecompSVD.h"
#include "TDecompSparse.h"
#include "TMatrixDSparse.h"
#include "TF1.h"
#include "TH1F.h"
#include "TNamed.h"
#include "TClonesArray.h"
#include "TMath.h"


#include "TTreeStream.h"
class TTreeSRedirector;
class THn;
class TObjString; 

class AliNDLocalRegression : public TNamed {
 public:
  AliNDLocalRegression();
  ~AliNDLocalRegression();

  Bool_t MakeFit(TTree * tree , const char *formulaVal, const char * formulaVar, const char*selection, const char * formulaKernel,  const char * dimensionFormula, Double_t weightCut=0.00001, Int_t entries=1000000000);

  Double_t Eval(Double_t *point);
  Double_t EvalError(Double_t *point);
  const THn *GetHistogram() {return fHistPoints;}
  void SetCuts(Double_t nSigma=6, Double_t robustFraction=0.95, Int_t estimator=1);
  void SetHistogram(THn* histo );
  void SetTree(TTree * tree) {fInputTree = tree;}
  TTreeSRedirector *GetStreamer(){return fStreamer;}
  void SetStreamer( TTreeSRedirector *streamer){ fStreamer=streamer;}
  //
  // function to access the Local Regression from the TFormula
  static void AddVisualCorrection(AliNDLocalRegression* corr, Int_t position=0);
  static AliNDLocalRegression*  GetVisualCorrection(Int_t position);
  static AliNDLocalRegression*  GetVisualCorrection(const char *corName){return (fgVisualCorrection==NULL) ? 0: ( AliNDLocalRegression*) fgVisualCorrection->FindObject(corName);}
  static Int_t GetVisualCorrectionIndex(const char *corName);
  Int_t GetVisualCorrectionIndex(){return GetVisualCorrectionIndex(GetName());}
  static TObjArray*  GetVisualCorrections() { return fgVisualCorrection;}
  //
  static Double_t GetCorrND(Double_t index, Double_t par0);
  static Double_t GetCorrNDError(Double_t index, Double_t par0);
  static Double_t GetCorrND(Double_t index, Double_t par0,Double_t par1);
  static Double_t GetCorrNDError(Double_t index, Double_t par0,Double_t par1);
  static Double_t GetCorrND(Double_t index, Double_t par0,Double_t par1, Double_t par2);
  static Double_t GetCorrNDError(Double_t index, Double_t par0,Double_t par1, Double_t par2);
  static Double_t GetCorrND(Double_t index, Double_t par0,Double_t par1, Double_t par2, Double_t par3);
  static Double_t GetCorrNDError(Double_t index, Double_t par0,Double_t par1, Double_t par2, Double_t par3);

 protected:
  Bool_t MakeRobustStatistic(TVectorD &values,TVectorD &errors,  TObjArray &pointArray,  TObjArray &kernelArray, Double_t weightCut, Double_t robustFraction);

  THn *fHistPoints;                     //  histogram local point distoribution
  Double_t fRobustFractionLTS;          //  fraction of data used for the robust mean and robust rms estimator (LTS https://en.wikipedia.org/wiki/Least_trimmed_squares)
  Double_t fRobustRMSLTSCut;            //  cut on the robust RMS  |value-localmean|<fRobustRMSLTSCut*localRMS
  Int_t    fCutType;                    //  type of the cut 0- no cut 1-cut localmean=median, 2-cut localmen=rosbut mean 
  TTree * fInputTree;                 // ! input tree - object not owner  
  TTreeSRedirector * fStreamer;       // ! streamer to keep - test intermediate data
  //
  TObjString *fFormulaVal;            // value:error definition formula
  TObjString *fSelection;             // point selector formula
  //
  TObjString *fFormulaVar;            // :separated variable definition formula
  TObjString *fKernelWidthFormula;    //: separated - kernel width for the regression
  TObjString *fPolDimensionFormula;   //: separated  - polynom for the regression
  //
  Int_t fNVar;                        // number of variables
  Int_t fNParameters;                 // number of local parameters to fit
  //
  TObjArray *fLocalFitParam;          // local fit parameters + RMS + chi2
  TObjArray *fLocalFitQuality;        // local fit npoints chi2
  TObjArray *fLocalFitCovar;          // local fit covariance matrix  
  //
  TMatrixD  *fLocalRobustStat;        // local robust statistic
private:
  static TObjArray *fgVisualCorrection; ///< array of orrection for visualization
  Int_t    *fBinIndex;                  //[fNParameters] working arrays current bin index
  Double_t *fBinCenter;                 //[fNParameters] working current local variables - bin center
  Double_t *fBinDelta;                  //[fNParameters] working current local variables - bin delta
  AliNDLocalRegression& operator=(const AliNDLocalRegression&);
  AliNDLocalRegression(const AliNDLocalRegression&);
  ClassDef(AliNDLocalRegression, 1);
};
#endif
