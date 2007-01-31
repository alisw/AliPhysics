#ifndef ALIGENSLOWNUCLEONS_H
#define ALIGENSLOWNUCLEONS_H
/* Copyright(c) 198-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
//
//  Generator for slow nucleons in pA interactions. 
//  Source is modelled by a relativistic Maxwell distributions.
//  Original code by  Ferenc Sikler  <sikler@rmki.kfki.hu>
//  This class: andreas.morsch@cern.ch
//
#include "AliGenerator.h"
class AliSlowNucleonModel;
class TH2F;
class TH1F;
class TF1;

class AliGenSlowNucleons : public AliGenerator
{
public:
    AliGenSlowNucleons();
    AliGenSlowNucleons(Int_t npart);
    virtual ~AliGenSlowNucleons();
    virtual void Init();
    virtual void FinishRun();
    virtual void Generate();
    virtual void SetPmax(Float_t pmax = 10.) {fPmax = pmax;}
    virtual void SetNominalCmsEnergy(Float_t energy = 14000.) {fCMS = energy;}
    virtual void SetTarget(Float_t a=208, Float_t z=82) {fATarget = a; fZTarget = z;}
    virtual void SetProtonDirection(Float_t dir = 1.);
    virtual void SetCharge(Int_t c = 1) {fCharge = c;}
    virtual void SetTemperature(Double_t t1 = 0.04, Double_t t2 = 0.004)
	{fTemperatureG = t1; fTemperatureB = t2;}
    virtual void SetBetaSource(Double_t b1 = 0.05, Double_t b2 = 0.)
	{fBetaSourceG = b1; fBetaSourceB = b2;}
    //
    virtual void SetSlowNucleonModel(AliSlowNucleonModel* model) 
	{fSlowNucleonModel = model;}
    virtual Bool_t NeedsCollisionGeometry() const {return kTRUE;}
    virtual void   SetCollisionGeometry(AliCollisionGeometry* geom)
	{fCollisionGeometry = geom;}
    virtual void   SetDebug(Int_t flag = 0) {fDebug = flag;}
    virtual void   SetNumbersOfSlowNucleons(Int_t ngp, Int_t ngn, Int_t nbp, Int_t nbn)
	{fNgp = ngp; fNgn = ngn; fNbp = nbp; fNbn = nbn;}
    //
    // Added by Chiara to take into account angular distribution 4 gray tracks
    virtual void   SetThetaDist(Int_t flag=0) {fThetaDistribution = flag;}
    
 protected:
    void     GenerateSlow(Int_t charge, Double_t T, Double_t beta, Float_t* q, Float_t &theta);
    Double_t Maxwell(Double_t m, Double_t p, Double_t t);
    void     Lorentz(Double_t m, Double_t beta, Float_t* q);
 protected:
    Float_t  fCMS;             // Center of mass energy
    Float_t  fMomentum;        // Target nucleus momentum
    Float_t  fBeta;            // Target nucleus beta
    Float_t  fPmax;            // Maximum slow nucleon momentum
    Float_t  fATarget;         // Target nucleus mass number
    Float_t  fZTarget;         // Target nucleus charge number
    Int_t    fCharge;          // Slow nucleon charge
    Float_t  fProtonDirection; // Direction of the proton
    Float_t  fTemperatureG;    // Source Temperature for gray nucleons
    Float_t  fBetaSourceG;     // Source beta for gray nucleons
    Float_t  fTemperatureB;    // Source Temperature for black nucleons
    Float_t  fBetaSourceB;     // Source beta for black nucleons
    Int_t    fNgp;             // Number of gray  protons
    Int_t    fNgn;             // Number of gray  neutrons
    Int_t    fNbp;             // Number of black protons
    Int_t    fNbn;             // Number of black neutrons
    Int_t    fDebug;           // Debug flag
    TH2F*    fDebugHist1;      // Histogram for debugging
    TH2F*    fDebugHist2;      // Histogram for debugging
    // Added by Chiara to take into account angular distribution 4 gray tracks
    Int_t    fThetaDistribution;// 0 -> flat dist., 1 -> fwd. peaked distribution
    TH1F*    fCosThetaGrayHist; // Histogram for debugging
    TF1*     fCosTheta;         // Function for non-uniform cos(theta) distribution
    
    //
    AliSlowNucleonModel* fSlowNucleonModel; // The slow nucleon model

 private:
    AliGenSlowNucleons(const AliGenSlowNucleons &sn);
    AliGenSlowNucleons & operator=(const AliGenSlowNucleons & rhs);

    ClassDef(AliGenSlowNucleons,2) // Slow Nucleon Generator
};
#endif






