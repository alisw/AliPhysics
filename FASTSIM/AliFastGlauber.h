#ifndef ALIFASTGLAUBER_H
#define ALIFASTGLAUBER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TObject.h>
class TF1;
class TF2;

class AliFastGlauber : public TObject {
 public:
    AliFastGlauber();
    virtual ~AliFastGlauber(){;}
    void SetWoodSaxonParameters(Double_t r0, Double_t d, Double_t w, Double_t n)
	{fWSr0 = r0; fWSd = d; fWSw = w; fWSn = n;}
    void SetMaxImpact(Float_t bmax = 20.) {fbMax = bmax;};
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
    TF2* Kernel()  {return fWStarfi;}
    TF1* Overlap() {return fWStaa;}
    void SimulateTrigger(Int_t n);
    void GetRandom(Float_t& b, Float_t& p, Float_t& mult);
    void GetRandom(Int_t& bin, Bool_t& hard);
    Float_t GetRandomImpactParameter(Float_t bmin, Float_t bmax);
 protected:
    static TF1*    fWSb;            // Wood-Saxon Function (b)
    static TF2*    fWSbz;           // Wood-Saxon Function (b, z)
    static TF1*    fWSz;            // Wood-Saxon Function (b = b0, z)
    static TF1*    fWSta;           // Thickness Function
    static TF2*    fWStarfi;        // Kernel for Overlap Function
    static TF1*    fWStaa;          // Overlap Function
    static TF2*    fWAlmond;        // Interaction Almond
    static TF1*    fWPathLength0;   // Path Length as a function of phi
    static TF1*    fWPathLength;    // Path Length as a function of phi
    static TF1*    fWIntRadius;     // Interaction Radius
    static TF1*    fWSgeo;          // dSigma/db geometric
    static TF1*    fWSbinary;       // dSigma/db binary
    static TF1*    fWSN;            // dN/db binary
    static TF1*    fWEnergyDensity; // Energy density as a function of impact parameter
    
    Float_t fWSr0;           // Wood-Saxon Parameter r0
    Float_t fWSd;            // Wood-Saxon Parameter d
    Float_t fWSw;            // Wood-Saxon Parameter w
    Float_t fWSn;            // Wood-Saxon Parameter n
    Float_t fSigmaHard;      // Hard Cross Section
    static Float_t fbMax;    // Maximum Impact Parameter
    
    ClassDef(AliFastGlauber,1) // Event geometry simulation in the Glauber Model
};

#endif 



