
#ifndef ALIPERFANALYZEINVPT_H
#define ALIPERFANALYZEINVPT_H

//----------------------------------------------------------------------------------------------------
// Class to analyse output of AliPerformancePtCalib.cxx and AliPerformancePtCalibMC.cxx:
// Projection of charge/pt vs theta and vs phi resp. histoprams will be fitted with either
// polynomial or gaussian fit function to extract minimum position of charge/pt.
// Fit options and theta, phi bins can be set by user.
// Attention: use the Set* functions of AliPerformancePtCalib.h and AliPerformancePtCalibMC.h
// when running AliPerformancePtCalib*::Analyse().
//
// Author: S. Schuchmann 11/13/2009 
//----------------------------------------------------------------------------------------------------

class TNamed;
class TString;
class TList;
class TF1;
class TH1F; 
class TH2F;
class TH1D;
class TGraphErrors;
//class TMath;
class TCanvas;
class TObjArray;


class AliPerfAnalyzeInvPt : public TNamed {

public:
   AliPerfAnalyzeInvPt();
   AliPerfAnalyzeInvPt(Char_t* name, Char_t* title);
   virtual ~AliPerfAnalyzeInvPt(){;}
   void InitGraphs(Double_t *binsXTheta,Double_t *fitParamTheta,Double_t *errFitParamTheta,Double_t *binsXPhi,Double_t *fitParamPhi,Double_t *errFitParamPhi);
   void InitFitFcn();
   void StartAnalysis(const TH2F *histThetaInvPt, const TH2F *histPhiInvPt, TObjArray* aFolderObj);
   void SetProjBinsPhi(const Double_t *pBins,const Int_t sizep);
   void SetProjBinsTheta(const Double_t *tBins, const Int_t sizet);
   void SetMakeFitOption(const Bool_t setGausFit, const Double_t exclusionR,const Double_t fitR );
   void SetDoRebin(const Int_t rebin){if(rebin) {fDoRebin = kTRUE; fRebin = rebin;}}
  
protected:
   Double_t fThetaBins[100];// array of theta bins for projection and fit
   Double_t fPhiBins[100];// array of phi bins for projection and fit
  
   Int_t fNThetaBins;// number of theta bins for projection and fit
   Int_t fNPhiBins; // number of phi bins for projection and fit
   Double_t fRange;// fit range of 1/pt spectrum
   Double_t fExclRange ;// range of exlusion of points around 0 when fitting 1/pt
   Bool_t fFitGaus;// set this flag for usage of gaussian fit function instead of polynomial (default)
   Bool_t  fDoRebin;
   Int_t  fRebin;
   // projection histos
   TH1D *fHistFitTheta[100];// projection histos for analysis in theta bins
   TH1D *fHistFitPhi[100];// projection histos for analysis in phi bins

   //histos for copy of input histos
   TH2F *fHistH2InvPtTheta; // for copy of input histos 
   TH2F *fHistH2InvPtPhi; // for copy of input histos
   
   //graphs for displaying the fit parameters (minimum position of 1/pt and error)
   TGraphErrors *fGrMinPosTheta;
   TGraphErrors *fGrMinPosPhi;

 

private:
   //fit functions
   TF1 *fFitMinPos;//fit function for first fit with polynomial
   TF1 *fFitMinPosRejP;//fit function for second fit with polynomial to improve result
   TF1 *fFitInvGauss;//fit functionfor first fit with gaussian
   TF1 *fFitInvGaussRejP;//fit function for second fit with gaussian to improve result

   // infput for fit functions
   static Double_t Polynomial(Double_t *x, const Double_t *par);
   static Double_t PolynomialRejP(Double_t *x, const Double_t *par);
   static Double_t InvGauss(Double_t *x, const Double_t *par);
   static Double_t InvGaussRejP(Double_t *x, const Double_t *par);

   //functions for fitting procedure
   void MakeFit(TH1D *dmproy, TF1 * fitpb, Double_t &mean, Double_t &ErrMean,  Double_t &excl,Double_t &range);
   void MakeFitBetter(TH1D *dmproy, TF1 * fitpb2, Double_t &mean, Double_t &ErrMean, Double_t &f, Double_t &excl, Double_t &range);
   void MakeFitInvGauss(TH1D *dmproy, TF1 * fitpb2, Double_t &mean, Double_t &ErrMean,Double_t &excl , Double_t &range);
   void MakeFitInvGaussBetter(TH1D *dmproy, TF1 * fitpb2, Double_t &mean, Double_t &ErrMean, Double_t &f,Double_t &excl,  Double_t &range);
 

  

   AliPerfAnalyzeInvPt(const AliPerfAnalyzeInvPt&);            // not implemented 
   AliPerfAnalyzeInvPt& operator=(const AliPerfAnalyzeInvPt&); // not implemented 

   ClassDef(AliPerfAnalyzeInvPt, 1); 
};

#endif
