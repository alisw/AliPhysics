#ifndef ALIFASTGLAUBER_H
#define ALIFASTGLAUBER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

// Utility class to make simple Glauber type calculations for collision geometries:
// Impact parameter, production points, reaction plane dependence
// Author: Andreas Morsch
// andreas.morsch@cern.ch

#include <TObject.h>
#include <TString.h>
class TF1;
class TF2;


class AliFastGlauber : public TObject {
 public:
    static AliFastGlauber* Instance();
    virtual ~AliFastGlauber();
    void Init(Int_t mode = 0);

    void SetWoodSaxonParameters(Double_t r0, Double_t d, Double_t w, Double_t n)
	{fWSr0 = r0; fWSd = d; fWSw = w; fWSn = n;}
    void SetWoodSaxonParametersAu()
	{fWSr0 = 6.38; fWSd = 0.535; fWSw = 0.; fWSn = 8.59e-4;}
    void SetWoodSaxonParametersPb()
	{fWSr0 = 6.78; fWSd = 0.54; fWSw = 0.; fWSn = 7.14e-4;}
    void SetMaxImpact(Float_t bmax = 20.) {fgBMax = bmax;};
    void SetHardCrossSection(Float_t xs = 1.0) {fSigmaHard = xs;}
    void SetNNCrossSection  (Float_t xs = 55.6) {fSigmaNN = xs;}
    void SetNucleus(Int_t n=208) {fA=n;}
    void SetAuAuRhic();
    void SetPbPbLHC();
    void SetFileName(const TString &fn){fName=fn;}
    void SetFileName(const char *fn="$(ALICE_ROOT)/FASTSIM/data/glauberPbPb.root"){fName=fn;}

    const TF1* GetWSB()            const {return fgWSb;}
    const TF1* GetRWSB()           const {return fgRWSb;}
    const TF2* GetWSbz()           const {return fgWSbz;}
    const TF1* GetWSz()            const {return fgWSz;} 
    const TF1* GetWSta()           const {return fgWSta;}
    const TF2* Kernel()            const {return fgWStarfi;}
    const TF2* GetWStarfi()        const {return fgWStarfi;}
    const TF2* GetWKParticipants() const {return fgWKParticipants;}
    const TF1* GetWParticipants()  const {return fgWParticipants;} 
    const TF1* Overlap()           const {return fgWStaa;}
    const TF1* GetWStaa()          const {return fgWStaa;} 
    const TF2* GetWAlmond()        const {return fgWAlmond;}
    const TF1* GetWPathLength0()   const {return fgWPathLength0;} 
    const TF1* GetWPathLength()    const {return fgWPathLength;}
    const TF1* GetWIntRadius()     const {return fgWIntRadius;}
    const TF1* GetWSgeo()          const {return fgWSgeo;}
    const TF1* GetWSbinary()       const {return fgWSbinary;}
    const TF1* GetWSN()            const {return fgWSN;}     
    const TF1* GetWEnergyDensity() const {return fgWEnergyDensity;} 
    const TF2* GetWAlmondFixedB(Int_t i) const {return fgWAlmondFixedB[i];}
    
    Float_t GetWr0() const {return fWSr0;}
    Float_t GetWSd() const {return fWSd;}
    Float_t GetWSw() const {return fWSw;}
    Float_t GetWSn() const {return fWSn;}
    Float_t GetSigmaHard()       const {return fSigmaHard;}
    Float_t GetSigmaNN()         const {return fSigmaNN;}
    Int_t GetA()                 const {return fA;}
    const TString* GetFileName() const {return &fName;}
    Float_t GetBmin() const {return fBmin;}
    Float_t GetBmax() const {return fBmax;}

    void DrawWSb()          const;
    void DrawThickness()    const;
    void DrawOverlap()      const;
    void DrawParticipants() const;
    void DrawGeo()          const;
    void DrawBinary()       const;
    void DrawN()            const;    
    void DrawKernel(Double_t b = 0.) const;
    void DrawAlmond(Double_t b = 0.) const;
    void DrawPathLength0(Double_t b = 0., Int_t iopt = 0)            const;
    void DrawPathLength(Double_t b, Int_t ni = 1000, Int_t iopt = 0) const;
    void DrawIntRadius(Double_t b = 0.) const;
    void DrawEnergyDensity()            const;
    
    Double_t CrossSection(Double_t b1, Double_t b2)               const;
    Double_t HardCrossSection(Double_t b1, Double_t b2)           const;
    Double_t NHard(Double_t b1, Double_t b2)                      const;
    Double_t FractionOfHardCrossSection(Double_t b1, Double_t b2) const;
    Double_t Binaries(Double_t b)                 const;
    Double_t GetNumberOfBinaries(Double_t b)      const;
    Double_t Participants(Double_t b)             const;
    Double_t GetNumberOfParticipants(Double_t  b) const;
    Double_t GetNumberOfCollisions(Double_t  b)   const;
    Double_t GetNumberOfCollisionsPerEvent(Double_t  b) const;
    Double_t MeanOverlap(Double_t b1, Double_t b2);
    Double_t MeanNumberOfCollisionsPerEvent(Double_t b1, Double_t b2);
    void SimulateTrigger(Int_t n);
    void GetRandom(Float_t& b, Float_t& p, Float_t& mult);
    void GetRandom(Int_t& bin, Bool_t& hard);
    Double_t GetRandomImpactParameter(Double_t bmin, Double_t bmax);

    void StoreFunctions() const;
    void StoreAlmonds()   const;

    void SetLengthDefinition(Int_t def=1) {fEllDef=def;}
    Int_t GetLengthDef() const {return fEllDef;}
    void SetCentralityClass(Double_t xsecFrLow=0.0,Double_t xsecFrUp=0.1);    
    void GetRandomBHard(Double_t& b);
    void GetRandomXY(Double_t& x,Double_t& y);
    void GetSavedXY(Double_t xy[2]) const {xy[0] = fXY[0]; xy[1] = fXY[1];} 
    void GetSavedI0I1(Double_t i0i1[2]) const {i0i1[0] = fI0I1[0]; i0i1[1] = fI0I1[1];}
    void SaveXY(Double_t x, Double_t y) {fXY[0] = x; fXY[1] = y;}
    void SaveI0I1(Double_t i0, Double_t i1) {fI0I1[0] = i0; fI0I1[1] = i1;}

    void GetRandomPhi(Double_t& phi);
    Double_t CalculateLength(Double_t b=0.,Double_t x0=0.,Double_t y0=0.,
                             Double_t phi0=0.);
    void GetLengthAndPhi(Double_t& ell,Double_t &phi,Double_t b=-1.);
    void GetLength(Double_t& ell,Double_t b=-1.);
    void GetLengthsBackToBackAndPhi(Double_t& ell1,Double_t& ell2,
				    Double_t &phi,
				    Double_t b=-1.);
    void GetLengthsBackToBack(Double_t& ell1,Double_t& ell2,
			      Double_t b=-1.);
    void GetLengthsForPythia(Int_t n,Double_t* const phi,Double_t* ell,
			     Double_t b=-1.);
    void PlotBDistr(Int_t n=1000);
    void PlotLengthDistr(Int_t n=1000,Bool_t save=kFALSE,
			 const char *fname="length.root");
    void PlotLengthB2BDistr(Int_t n=1000,Bool_t save=kFALSE,
			    const char *fname="lengthB2B.root");
    void CalculateI0I1(Double_t& integral0,Double_t& integral1,
		       Double_t b=0.,
		       Double_t x0=0.,Double_t y0=0.,Double_t phi0=0.,
		       Double_t ellCut=20.) const;
    void GetI0I1AndPhi(Double_t& integral0,Double_t& integral1,Double_t &phi,
		 Double_t ellCut=20.,Double_t b=-1.);
    void GetI0I1(Double_t& integral0,Double_t& integral1,
		 Double_t ellCut=20.,Double_t b=-1.);
    void GetI0I1BackToBackAndPhi(Double_t& integral01,Double_t& integral11,
				 Double_t& integral02,Double_t& integral12,
				 Double_t& phi,
				 Double_t ellCut=20.,Double_t b=-1.);
    void GetI0I1BackToBackAndPhiAndXY(Double_t& integral01,Double_t& integral11,
				      Double_t& integral02,Double_t& integral12,
				      Double_t& phi,Double_t& x,Double_t&y,
				      Double_t ellCut=20.,Double_t b=-1.);
    void GetI0I1BackToBack(Double_t& integral01,Double_t& integral11,
			   Double_t& integral02,Double_t& integral12,
			   Double_t ellCut=20.,Double_t b=-1.);
    void GetI0I1ForPythia(Int_t n,Double_t* phi,
			  Double_t* integral0,Double_t* integral1,
			  Double_t ellCut=20.,Double_t b=-1.);
    void GetI0I1ForPythiaAndXY(Int_t n,Double_t* phi,
			       Double_t* integral0,Double_t* integral1,
			       Double_t&x, Double_t &y,
			       Double_t ellCut=20.,Double_t b=-1.);
    void PlotI0I1Distr(Int_t n=1000,Double_t ellCut=20.,Bool_t save=kFALSE,
		       const char *fname="i0i1.root");
    void PlotI0I1B2BDistr(Int_t n=1000,Double_t ellCut=20.,Bool_t save=kFALSE,
			  const char *fname="i0i1B2B.root");
    void PlotAlmonds() const;
    // Copy
    AliFastGlauber& operator=(const AliFastGlauber & rhs);
    void Copy(TObject&) const;
 protected:
    static Double_t RWSb           (const Double_t *xx, const Double_t *par);
    static Double_t WSb            (const Double_t *xx, const Double_t *par);
    static Double_t WSbz           (const Double_t *xx, const Double_t *par);
    static Double_t WSz            (const Double_t *xx, const Double_t *par);
    static Double_t WSta           (const Double_t *xx, const Double_t *par);
    static Double_t WStarfi        (const Double_t *xx, const Double_t *par);
    static Double_t WStaa          (const Double_t *xx, const Double_t *par);
    static Double_t WKParticipants (const Double_t *xx, const Double_t *par);
    static Double_t WParticipants  (const Double_t *xx, const Double_t *par);    
    static Double_t WSgeo          (const Double_t *xx, const Double_t *par);
    static Double_t WSbinary       (const Double_t *xx, const Double_t *par);
    static Double_t WSN            (const Double_t *xx, const Double_t *par);
    static Double_t WAlmond        (const Double_t *xx, const Double_t *par);
    static Double_t WPathLength0   (const Double_t *xx, const Double_t *par);
    static Double_t WPathLength    (const Double_t *xx, const Double_t *par);
    static Double_t WIntRadius     (const Double_t *xx, const Double_t *par);
    static Double_t WEnergyDensity (const Double_t *xx, const Double_t *par);

    void Reset() const;
 private:
    AliFastGlauber();
    AliFastGlauber(const AliFastGlauber& glauber);

    static Float_t         fgBMax;       // Maximum Impact Parameter
    static const Int_t     fgkMCInts;    // Number of MC integrations
    static AliFastGlauber* fgGlauber;    // Singleton instance
     
    
    static TF1*    fgWSb;            // Wood-Saxon Function (b)
    static TF1*    fgRWSb;           // Wood-Saxon Function (b) with phase space factor
    static TF2*    fgWSbz;           // Wood-Saxon Function (b, z)
    static TF1*    fgWSz;            // Wood-Saxon Function (b = b0, z)
    static TF1*    fgWSta;           // Thickness Function
    static TF2*    fgWStarfi;        // Kernel for Overlap Function
    static TF2*    fgWKParticipants; // Kernel for number of participants
    static TF1*    fgWParticipants;  // Number of participants
    static TF1*    fgWStaa;          // Overlap Function
    static TF2*    fgWAlmond;        // Interaction Almond
    static TF1*    fgWPathLength0;   // Path Length as a function of phi
    static TF1*    fgWPathLength;    // Path Length as a function of phi
    static TF1*    fgWIntRadius;     // Interaction Radius
    static TF1*    fgWSgeo;          // dSigma/db geometric
    static TF1*    fgWSbinary;       // dSigma/db binary
    static TF1*    fgWSN;            // dN/db binary
    static TF1*    fgWEnergyDensity; // Energy density as a function of impact parameter
    static TF2*    fgWAlmondFixedB[40]; // Interaction Almonds read from file
    static TF2*    fgWAlmondCurrent;    // Interaction Almond used for length
    
    Float_t fWSr0;      // Wood-Saxon Parameter r0
    Float_t fWSd;       // Wood-Saxon Parameter d
    Float_t fWSw;       // Wood-Saxon Parameter w
    Float_t fWSn;       // Wood-Saxon Parameter n
                        // (chosen such that integral is one)
    Float_t fSigmaHard; // Hard Cross Section [mbarn]
    Float_t fSigmaNN;   // NN Cross Section [mbarn]   
    Int_t fA;           // Nucleon number of nucleus A

    Float_t fBmin;      // Minimum b (set through centrality selection)
    Float_t fBmax;      // Coresponding maximum b
    Double_t fXY[2];    // Current generated production point
    Double_t fI0I1[2];  // Current integrals I0 and I1
    Int_t fEllDef;      // definition of length (see CalculateLength())
    TString fName;      // filename of stored distributions
    ClassDef(AliFastGlauber,2) // Event geometry simulation in the Glauber Model
};

#endif 
