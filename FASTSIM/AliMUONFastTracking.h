#ifndef ALIMUONFASTTRACKING_H
#define ALIMUONFASTTRACKING_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//        Class AliMUONFastTracking 
//
//  Manager for the fast simulation of tracking in the muon spectrometer
//  This class reads the lookup tables containing the parameterization 
//  of the deltap, deltatheta, deltaphi for different background levels
//  and provides the related smeared parameters   
//-------------------------------------------------------------------------

class TF1;
class TSpline3;
class TFile;
class AliMUONFastTrackingEntry;


#include <TObject.h>

enum LUTClusterType {kOld, kNew};

class AliMUONFastTracking :  public TObject {
 public:
    static  AliMUONFastTracking* Instance();
    ~AliMUONFastTracking(){;}
    void Init(Float_t bkg);
    void ReadLUT(TFile *file);
    void GetBinning(Int_t &nbinp, Float_t &pmin, Float_t &pmax,
		    Int_t &nbintheta, Float_t &thetamin, Float_t &thetamax,
		    Int_t &nbinphi, Float_t &phimin, Float_t &phimax) const;
    void GetIpIthetaIphi(Float_t p, Float_t theta, Float_t phi, Int_t charge,
			 Int_t &ip, Int_t &itheta, Int_t &iphi) const;
    void GetSplit(Int_t ip, Int_t itheta, Int_t &nSplitP, Int_t &nSplitTheta) const;
    Float_t Efficiency(Float_t p, Float_t theta, Float_t phi, Int_t charge);
    Float_t Acceptance(Float_t p, Float_t theta, Float_t phi, Int_t charge); 
    Float_t MeanP(Float_t p, Float_t theta, Float_t phi, Int_t charge)     const;
    Float_t SigmaP(Float_t p, Float_t theta, Float_t phi, Int_t charge)    const;
    Float_t Sigma1P(Float_t p, Float_t theta, Float_t phi, Int_t charge)   const;
    Float_t NormG2(Float_t p, Float_t theta, Float_t phi, Int_t charge)    const;
    Float_t MeanG2(Float_t p, Float_t theta, Float_t phi, Int_t charge)    const;
    Float_t SigmaG2(Float_t p, Float_t theta, Float_t phi, Int_t charge)   const;
    Float_t MeanTheta(Float_t p, Float_t theta, Float_t phi, Int_t charge) const;
    Float_t SigmaTheta(Float_t p, Float_t theta, Float_t phi, Int_t charge)const;  
    Float_t MeanPhi(Float_t p, Float_t theta, Float_t phi, Int_t charge)   const;
    Float_t SigmaPhi(Float_t p, Float_t theta, Float_t phi, Int_t charge);

    void SetSpline();
    Float_t GetBackground() const {return fBkg;}
    void SetLUTClusterFinder(LUTClusterType clusterFinder) { fClusterFinder = clusterFinder;}
    void SetBackground(Float_t bkg);
    void UseSpline (Int_t splineSwitch=1) {fSpline = splineSwitch;}
    TF1* GetFitP(Int_t ip, Int_t itheta, Int_t iphi); 
 protected:
    Int_t   fNbinp;         // n. of momentum bins in the lookup table 
    Float_t fPmin;          // min. value of momentum parameterized in LUT
    Float_t fPmax;          // max. value of momentum parameterized in LUT
    Float_t fDeltaP;        // momentum bin width 
    Int_t   fNbintheta;     // n. of theta bins in the lookup table 
    Float_t fThetamin;      // min. value of theta parameterized in LUT
    Float_t fThetamax;      // max. value of theta parameterized in LUT
    Float_t fDeltaTheta;    // theta bin width
    Int_t   fNbinphi;       // n. of phi bins in the lookup table 
    Float_t fPhimin;        // min. value of phi parameterized in LUT
    Float_t fPhimax;        // min. value of phi parameterized in LUT
    Float_t fDeltaPhi;      // phi bin width
    Int_t   fPrintLevel;    // level of information printed for debugging
    Float_t fBkg;           // soft background level  
    TF1 *fFitp[20][20][20];                    // func for psmear-pgen distr
    AliMUONFastTrackingEntry *fEntry[20][20][20][4]; // array of LUT parameters
    AliMUONFastTrackingEntry *fCurrentEntry[20][20][20]; // array of LUT parameters
    TSpline3 *fSplineEff[200][3];        // spline funcs for efficiency
    TSpline3 *fSplineAcc[200][3];        // spline funcs for acceptance
    TSpline3 *fSplineSigmap[200][3];     // spl.funcs for dp distribution width
    TSpline3 *fSplineSigma1p[200][3];    // spl.funcs for dp distr. width correction (see function FitP)
    TSpline3 *fSplineSigmatheta[200][3]; // spl.funcs for dtheta distr. width
    TSpline3 *fSplineSigmaphi[200][3];   // spl.funcs for dphi distr. width
    Int_t fSpline;                       // switches on/off the use of spline
    LUTClusterType fClusterFinder;       // type of cluster finder (old/new)
    static AliMUONFastTracking*    fgMUONFastTracking; //!Pointer to single instance
    ClassDef(AliMUONFastTracking,1)      // Fast MUON Tracking Data Handler
 private:
    AliMUONFastTracking();
    AliMUONFastTracking(Float_t /*bkg*/){;}
    AliMUONFastTracking(const AliMUONFastTracking &ft);
    void Copy(TObject &) const;
    AliMUONFastTracking& operator=(const AliMUONFastTracking & rhs);
};

#endif

