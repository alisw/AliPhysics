#ifndef ALIBLASTWAVEFITTER_H
#define ALIBLASTWAVEFITTER_H


/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliBlastwaveFitter.h 49869 2012-05-17 04:49:51Z fnoferin $ */

/////////////////////////////////////////////////
//                                             //
//          Blastwave Fitter Class             //
//           noferini@bo.infn.it               //
/////////////////////////////////////////////////

#include "TF1.h"
#include "TF2.h"
#include "TGraph.h"
#include "TMinuit.h"

#include "AliBlastwaveFit.h"

class AliBlastwaveFitter : public TNamed
{
public:
    AliBlastwaveFitter(const char *name);
    AliBlastwaveFitter();
    ~AliBlastwaveFitter();
    
    Int_t AddFitFunction(AliBlastwaveFit *fitf){if(fNfunction < fgNmaxFunction){fFunc[fNfunction]=fitf;fNfunction++;return 0;} else return 1;};

    void ResetFunctions(){fNfunction=0; for(Int_t i=0;i<fgNmaxFunction;i++) fFunc[i] = NULL;};

    Int_t CheckAvailability();
    Int_t PrepareToFit();
    Int_t Fit();

    static void FCN(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
    static Double_t GetChi2(){return fgChi2;};
    static Int_t GetNDGF(){return fgNDGF;};

    TMinuit *GetMinuit(){return fMinuit;};

    Int_t GetNFunctions(){return fNfunction;};
    AliBlastwaveFit *GetFunction(Int_t i){if(i>= 0 && i < fNfunction) return fFunc[i]; else return NULL;};

    void SetMinos(Bool_t flag=kTRUE){fMinos=flag;};

    TGraph* DoContour(Int_t np,Int_t ip1,Int_t ip2,Float_t nsigma=1);

    TGraph* DoContourBetaT(Int_t np,Int_t iBoostOrBeta,Int_t iT,Float_t nsigma=1);
    TGraph* ConvertContourFromBoostToBeta(TGraph *g,Int_t iBoostOrBeta,Int_t iT);

private:
    AliBlastwaveFitter(const AliBlastwaveFitter & old); 
    AliBlastwaveFitter& operator=(const AliBlastwaveFitter & source); // ass. op. 

    static const Int_t fgNmaxFunction = 20;   // max number of function allowed
    Int_t fNfunction;                         // number of functions to be fitted
    AliBlastwaveFit *fFunc[fgNmaxFunction];   // functions to be fitted

    Int_t fSpectraFlag[fgNmaxFunction];       // availability for spectra fit (1=TH1, 2=TGraphErrors)
    Int_t fV2Flag[fgNmaxFunction];       // availability for v2 fit (1=TH1, 2=TGraphErrors)

    TMinuit *fMinuit;                    // minuit object

    // information copied in the static variable to perform the fit
    static Int_t fgNparReal;
    static Int_t fgNDGF;
    static Int_t fgNfunctionCurrent;
    static AliBlastwaveFit *fgFuncC[50];   // functions to be fitted
    static Int_t fgSpectraFlagC[50];       // availability for spectra fit (1=TH1, 2=TGraphErrors)
    static Int_t fgV2FlagC[50];       // availability for v2 fit (1=TH1, 2=TGraphErrors)
    static Double_t fgChi2;            // Chi2

    Bool_t fMinos;                   // switch to try also minos

  ClassDef(AliBlastwaveFitter,1)  // fitter bases on minuit for blastwave
};

#endif


