#ifndef ALIFASTGLAUBER_H
#define ALIFASTGLAUBER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
//
// Utility class to make simple Glauber type calculations for collision geometries:
// Impact parameter, production points, reaction plane dependence
//
// Author: andreas.morsch@cern.ch

#include <TObject.h>
#include <TF2.h>
class TF1;
class TF2;

class AliFastGlauber : public TObject {
 public:
    AliFastGlauber();
    virtual ~AliFastGlauber(){;}
    void SetWoodSaxonParameters(Double_t r0, Double_t d, Double_t w, Double_t n)
	{fWSr0 = r0; fWSd = d; fWSw = w; fWSn = n;}
    void SetMaxImpact(Float_t bmax = 20.) {fgBMax = bmax;};
    void SetHardCrossSection(Float_t xs = 6.6) {fSigmaHard = xs;}

    static Double_t WSb            (Double_t *xx, Double_t *par);
    static Double_t WSbz           (Double_t *xx, Double_t *par);
    static Double_t WSz            (Double_t *xx, Double_t *par);
    static Double_t WSta           (Double_t *xx, Double_t *par);
    static Double_t WStarfi        (Double_t *xx, Double_t *par);
    static Double_t WStaa          (Double_t *xx, Double_t *par);
    static Double_t WSgeo          (Double_t *xx, Double_t *par);
    static Double_t WSbinary       (Double_t *xx, Double_t *par);
    static Double_t WSN            (Double_t *xx, Double_t *par);
    static Double_t WAlmond        (Double_t *xx, Double_t *par);
    static Double_t WPathLength0   (Double_t *xx, Double_t *par);
    static Double_t WPathLength    (Double_t *xx, Double_t *par);
    static Double_t WIntRadius     (Double_t *xx, Double_t *par);
    static Double_t WEnergyDensity (Double_t *xx, Double_t *par);
    
    void Init(Int_t mode = 0);
    void DrawWSb();
    void DrawThickness();
    void DrawOverlap();
    void DrawGeo();
    void DrawBinary();
    void DrawN();    
    void DrawKernel(Double_t b = 0.);
    void DrawAlmond(Double_t b = 0.);
    void DrawPathLength0(Double_t b = 0., Int_t iopt = 0);
    void DrawPathLength(Double_t b, Int_t ni = 1000, Int_t iopt = 0);
    void DrawIntRadius(Double_t b = 0.);
    void DrawEnergyDensity();
    
    Double_t CrossSection(Double_t b1, Double_t b2);
    Double_t FractionOfHardCrossSection(Double_t b1, Double_t b2);
    Double_t Binaries(Double_t b);
    TF2* Kernel()  {return fgWStarfi;}
    TF1* Overlap() {return fgWStaa;}
    void SimulateTrigger(Int_t n);
    void GetRandom(Float_t& b, Float_t& p, Float_t& mult);
    void GetRandom(Int_t& bin, Bool_t& hard);
    Float_t GetRandomImpactParameter(Float_t bmin, Float_t bmax);


    void SetLengthDefinition(Int_t def=1) { fEllDef=def; }
    void SetCentralityClass(Double_t xsecFrLow=0.0,Double_t xsecFrUp=0.1);    
    void StoreAlmonds();
    void GetRandomBHard(Double_t& b);
    void GetRandomXY(Double_t& x,Double_t& y);
    void GetRandomPhi(Double_t& phi);
    Double_t CalculateLength(Double_t b=0.,Double_t x0=0.,Double_t y0=0.,
			     Double_t phi0=0.);
    void GetLength(Double_t& ell,Double_t b=-1.);
    void GetLengthsBackToBack(Double_t& ell1,Double_t& ell2,Double_t b=-1.);
    void GetLengthsForPythia(Int_t n,Double_t* phi,Double_t* ell,
			     Double_t b=-1.);
    void PlotBDistr(Int_t n=1000);
    void PlotLengthDistr(Int_t n=1000,Bool_t save=kFALSE,
			 Char_t *fname="length.root");
    void PlotLengthB2BDistr(Int_t n=1000,Bool_t save=kFALSE,
			    Char_t *fname="lengthB2B.root");
    void PlotAlmonds();

 protected:
    static TF1*    fgWSb;            // Wood-Saxon Function (b)
    static TF2*    fgWSbz;           // Wood-Saxon Function (b, z)
    static TF1*    fgWSz;            // Wood-Saxon Function (b = b0, z)
    static TF1*    fgWSta;           // Thickness Function
    static TF2*    fgWStarfi;        // Kernel for Overlap Function
    static TF1*    fgWStaa;          // Overlap Function
    static TF2*    fgWAlmond;        // Interaction Almond
    static TF1*    fgWPathLength0;   // Path Length as a function of phi
    static TF1*    fgWPathLength;    // Path Length as a function of phi
    static TF1*    fgWIntRadius;     // Interaction Radius
    static TF1*    fgWSgeo;          // dSigma/db geometric
    static TF1*    fgWSbinary;       // dSigma/db binary
    static TF1*    fgWSN;            // dN/db binary
    static TF1*    fgWEnergyDensity; // Energy density as a function of impact parameter
    TF2  fWAlmondFixedB[40]; // Interaction Almonds read from file
    TF2*    fWAlmondCurrent; // Interaction Almond used for length
    
    Float_t fWSr0;           // Wood-Saxon Parameter r0
    Float_t fWSd;            // Wood-Saxon Parameter d
    Float_t fWSw;            // Wood-Saxon Parameter w
    Float_t fWSn;            // Wood-Saxon Parameter n
    Float_t fSigmaHard;      // Hard Cross Section
    static Float_t fgBMax;   // Maximum Impact Parameter
    
    Int_t fEllDef;           // definition of length (see CalculateLength())

    ClassDef(AliFastGlauber,1) // Event geometry simulation in the Glauber Model
};

#endif 



