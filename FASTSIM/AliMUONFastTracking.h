#ifndef ALIMUONFASTTRACKING
#define ALIMUONFASTTRACKING
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
class TF1;
class TSpline3;
class TFile;
class AliMUONFastTrackingEntry;


#include <TObject.h>

class AliMUONFastTracking :  public TObject {
 public:
    static  AliMUONFastTracking* Instance();
    ~AliMUONFastTracking(){;}
    void Init(Float_t bkg);
    void ReadLUT(TFile *file);
    void GetBinning(Int_t &nbinp, Float_t &pmin, Float_t &pmax,
		    Int_t &nbintheta, Float_t &thetamin, Float_t &thetamax,
		    Int_t &nbinphi, Float_t &phimin, Float_t &phimax);
    void GetIpIthetaIphi(Float_t p, Float_t theta, Float_t phi, Int_t charge,
			 Int_t &ip, Int_t &itheta, Int_t &iphi);
    void GetSplit(Int_t ip, Int_t itheta, Int_t &nSplitP, Int_t &nSplitTheta);
    Float_t Efficiency(Float_t p, Float_t theta, Float_t phi, Int_t charge);
    Float_t Acceptance(Float_t p, Float_t theta, Float_t phi, Int_t charge);
    Float_t MeanP(Float_t p, Float_t theta, Float_t phi, Int_t charge);
    Float_t SigmaP(Float_t p, Float_t theta, Float_t phi, Int_t charge);
    Float_t Sigma1P(Float_t p, Float_t theta, Float_t phi, Int_t charge);
    Float_t NormG2(Float_t p, Float_t theta, Float_t phi, Int_t charge);
    Float_t MeanG2(Float_t p, Float_t theta, Float_t phi, Int_t charge);
    Float_t SigmaG2(Float_t p, Float_t theta, Float_t phi, Int_t charge);
    Float_t MeanTheta(Float_t p, Float_t theta, Float_t phi, Int_t charge);
    Float_t SigmaTheta(Float_t p, Float_t theta, Float_t phi, Int_t charge);  
    Float_t MeanPhi(Float_t p, Float_t theta, Float_t phi, Int_t charge);
    Float_t SigmaPhi(Float_t p, Float_t theta, Float_t phi, Int_t charge);

    void SetSpline();
    Float_t GetBackground() {return fBkg;}
    void SetBackground(Float_t bkg);
    void UseSpline (Int_t splineSwitch=1) {fSpline = splineSwitch;}
    void SmearMuon(Float_t pgen, Float_t thetagen, Float_t phigen, Int_t charge,
		   Float_t &psmear, Float_t &thetasmear, Float_t &phismear,
		   Float_t &eff, Float_t &acc);
    TF1* GetFitP() {return fFitp;}
 private:
    AliMUONFastTracking();
    AliMUONFastTracking(Float_t bkg){;}
 protected:
    Int_t   fNentries;
    Int_t   fNbinp; 
    Float_t fPmin;
    Float_t fPmax;
    Float_t fDeltaP;
    Int_t   fNbintheta;
    Float_t fThetamin;
    Float_t fThetamax;
    Float_t fDeltaTheta;
    Int_t   fNbinphi;
    Float_t fPhimin;
    Float_t fPhimax;
    Float_t fDeltaPhi;
    Int_t   fPrintLevel;
    Float_t fBkg;
    TF1 *fFitp;                                   // func for psmear-pgen distr
    AliMUONFastTrackingEntry *fEntry[20][20][20][4]; // array of LUT parameters
    AliMUONFastTrackingEntry *fCurrentEntry[20][20][20]; // array of LUT parameters
 public:
    TSpline3 *fSplineEff[200][3];                 // spline funcs for efficiency
    TSpline3 *fSplineAcc[200][3];                 // spline funcs for acceptance
    TSpline3 *fSplineSigmap[200][3];              // 
    TSpline3 *fSplineSigma1p[200][3];             //!
    TSpline3 *fSplineSigmatheta[200][3];          //!
    TSpline3 *fSplineSigmaphi[200][3];            //!
 protected: 
    Int_t fSpline;
    static AliMUONFastTracking*    fgMUONFastTracking; //!Pointer to single instance
    ClassDef(AliMUONFastTracking,1)                    // Fast MUON Tracking Data Handler
};

#endif

