#ifndef ALIITSBASEGEOMETRY_H
#define ALIITSBASEGEOMETRY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
  $Id$
 */

/////////////////////////////////////////////////////////////////////////
//  A basic geometry class for the ITS simulation geometry stucture
/////////////////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TArrayI.h>
#include "AliModule.h"
class TString;

class AliITSBaseGeometry : public TObject {
 public:
    AliITSBaseGeometry(); // Default constructor
    AliITSBaseGeometry(AliModule *its,Int_t iflag); // Standard Constructor
    virtual ~AliITSBaseGeometry(); // Destructor
    virtual void BuildDisplayGeometry(){}; // Calls ROOT geometry interface
                                      // to AliRoot display
    virtual void CreateG3Geometry(){}; // Calls Geant3 interface geometry routines
    virtual void CreateG3Materials(){}; // Calls Geant3 interface for materials
    virtual Int_t IsVersion() const{return 11;}// return version of geometry.

    Int_t ITSG3VnameToIndex(const char name[3])const; // Get Index for Geant3 v name
    char* ITSIndexToITSG3name(const Int_t i); // Get Geant3 volume name
    Int_t AddVolName(const TString name); // Add volumen name to list
    TString GetVolName(const Int_t i)const; // Return volume name at index
    Int_t GetVolumeIndex(const TString &a);
    void SetScalecm(){fScale = 1.0;}// Sets scale factor for centemeters
    void SetScalemm(){fScale = 0.10;}// Sets scale factor for milimeters
    void SetScalemicrons(){fScale = 1.0E-04;}// Sets scale factor for microns
    void SetScale(Double_t s=1.0){fScale = s;}// Sets scale factor
    Double_t GetScale()const{return fScale;}// Returns the scale factor
    Bool_t IsScalecm()const{// Returens kTRUE if scale factor is set of [cm]
        if(fScale==1.0) return kTRUE; return kFALSE;}
    // Create a Box
    void Box(const char gnam[3],const TString &dis,
             Double_t dx,Double_t dy,Double_t dz,Int_t med);
    // Greate A Trapizoid with the x dimension varing along z.
    void Trapezoid1(const char gnam[3],const TString &dis,Double_t dxn,
                    Double_t dxp,Double_t dy,Double_t dz,Int_t med);
    // Greate A Trapizoid with the x and y dimension varing along z.
    void Trapezoid2(const char gnam[3],const TString &dis,Double_t dxn,
                    Double_t dxp,Double_t dyn,Double_t dyp,Double_t dz,
                    Int_t med);
    // General trapazoid.
    void Trapezoid(const char gnam[3],const TString &dis,Double_t dz,
                   Double_t thet,Double_t phi,Double_t h1,Double_t bl1,
                   Double_t tl1,Double_t alp1,Double_t h2,Double_t bl2,
                   Double_t tl2,Double_t alp2,Int_t med);
    // Simple Tube.
    void Tube(const char gnam[3],const TString &dis,Double_t rmin,
              Double_t rmax,Double_t dz,Int_t med);
    // Tube segment.
    void TubeSegment(const char gnam[3],const TString &dis,Double_t rmin,
                     Double_t rmax,Double_t dz,Double_t phi1,Double_t phi2,
                     Int_t med);
    // Simple Cone.
    void Cone(const char gnam[3],const TString &dis,Double_t dz,Double_t rmin1,
              Double_t rmax1,Double_t rmin2,Double_t rmax2,Int_t med);
    // Segment of a Cone.
    void ConeSegment(const char gnam[3],const TString &dis,Double_t dz,
                     Double_t rmin1,Double_t rmax1,Double_t rmin2,
                     Double_t rmax2,Double_t phi1,Double_t phi2,Int_t med);
    // Spherical shell segment.
    void Sphere(const char gnam[3],const TString &dis,Double_t rmin,
                Double_t rmax,Double_t the1,Double_t the2,Double_t phi1,
                Double_t phi2,Int_t med);
    // Parallelepiped.
    void Parallelepiped(const char gnam[3],const TString &dis,Double_t dx,
                        Double_t dy,Double_t dz,Double_t alph,Double_t thet,
                        Double_t phi,Int_t med);
    // Polygon.
    void Polygon(const char gnam[3],const TString &dis,Double_t phi1,
                 Double_t dphi,Int_t npdv,Int_t nz,Double_t *z,Double_t *rmin,
                 Double_t *rmax,Int_t med);
    //Poly-Cone
    void PolyCone(const char gnam[3],const TString &dis,Double_t phi1,
                  Double_t dphi,Int_t nz,Double_t *z,Double_t *rmin,
                  Double_t *rmax,Int_t med);
    // Ellliptical cross-sectino tube
    void TubeElliptical(const char gnam[3],const TString &dis,Double_t p1,
                        Double_t p2,Double_t dz,Int_t med);
    // Hyperbolic tube
    void HyperbolicTube(const char gnam[3],const TString &dis,Double_t rmin,
                        Double_t rmax,Double_t dz,Double_t thet,Int_t med);
    // Twisted genral trapezoid.
    void TwistedTrapezoid(const char gnam[3],const TString &dis,Double_t dz,
                          Double_t thet,Double_t phi,Double_t twist,
                          Double_t h1,Double_t bl1,Double_t tl1,
                          Double_t apl1,Double_t h2,Double_t bl2,
                          Double_t tl2,Double_t apl2,Int_t med);
    // Cut tube.
    void CutTube(const char gnam[3],const TString &dis,Double_t rmin,
                 Double_t rmax,Double_t dz,Double_t phi1,Double_t phi2,
                 Double_t lx,Double_t ly,Double_t lz,Double_t hx,Double_t hy,
                 Double_t hz,Int_t med);
    // Position one volume inside another
    void Pos(const char vol[3],Int_t cn,const char moth[3],Double_t x,
             Double_t y,Double_t z,Int_t irot);
    void SetMedArray(){// Sets up the array of media
        fidmed = ((fits->GetIdtmed())->GetArray())-199;}// Define rotation matrix
    void Matrix(Int_t irot,Double_t thet1,Double_t phi1,Double_t thet2,
                Double_t phi2,Double_t thet3,Double_t phi3);
    // Defube ritatuib matrix
    void Matrix(Int_t irot,Double_t rot[3][3]);
    // Rotation matrix about axis i (i=0=>x, i=1=>y, i=2=>z).
    void Matrix(Int_t irot,Int_t axis,Double_t thet);
    // Rotation matrix about x axis
    void XMatrix(Int_t irot,Double_t thet){Matrix(irot,0,thet);}
    // Rotation matrix about y axis
    void YMatrix(Int_t irot,Double_t thet){Matrix(irot,1,thet);}
    // Rotation matrix about z axis
    void ZMatrix(Int_t irot,Double_t thet){Matrix(irot,2,thet);}
    // Define Element material and medium
    void Element(Int_t imat,const char *name,Int_t z,Double_t dens,Int_t istd);
    // Define Material by constituant weights
    void MixtureByWeight(Int_t imat,const char *name,Int_t *z,Double_t *w,
                         Double_t dens,Int_t nelments,Int_t istd);
    // Define Material by constituant relative number
    void MixtureByNumber(Int_t imat,const char *name,Int_t *z,Int_t *i,
                         Double_t dens,Int_t nelments,Int_t istd);
    // Returns standard radiation lenghts of elements.
    Float_t GetRadLength(Int_t z){return RadLength(z,(Double_t)GetA(z));}
    // Returns natrual abundance atomic mass numbers for a given element
    Float_t GetA(Int_t z);
    // Returns ITS standard Theata Max transport cut values
    Float_t GetStandardThetaMax(Int_t istd);
    // Returns ITS standard Max step size transport cut values
    Float_t GetStandardMaxStepSize(Int_t istd);
    // Returns ITS standard frational energy transport cut values
    Float_t GetStandardEfraction(Int_t istd);
    // Returns ITS standard epsilon transport cut values
    Float_t GetStandardEpsilon(Int_t istd);
    // Degree Versions of TMath functions (as needed)
    Double_t Sind(Double_t t){return TMath::Sin(TMath::Pi()*t/180.);}
    Double_t Cosd(Double_t t){return TMath::Cos(TMath::Pi()*t/180.);}
    Double_t Tand(Double_t t){return TMath::Tan(TMath::Pi()*t/180.);}
    Double_t ASind(Double_t t){return 180.0*TMath::ASin(t)/TMath::Pi();}
    Double_t ACosd(Double_t t){return 180.0*TMath::ACos(t)/TMath::Pi();}
    Double_t ATand(Double_t t){return 180.0*TMath::ATan(t)/TMath::Pi();}
    Double_t ATand2(Double_t y,Double_t x){return 180.0*TMath::ATan2(y,x)/
					       TMath::Pi();}
    Double_t RadLength(Int_t iz,Double_t a); // Computes radiation length
                                             // for an element
 private:
    static Int_t fNCreates; //! Counts the number of time this class has
    // been created.
    static const Double_t fAlpha = 7.297352533e-3; //! find structure constant
    static const Double_t fRe = 2.81794028e-13;//![cm]classical electron radius
    static const Double_t fNa = 6.02214199e+23; //! [#/mole] Avogadro's number
    static Int_t *fidrot;
    static Int_t fidrotsize;
    static Int_t fidrotlast;
    static TString *fVolName; // Array of ITS Volumen names.
    static Int_t fVolNameSize; // Size of Array fVolName
    static Int_t fVolNameLast; // Last filled element of fVolName
    Double_t fScale; // Scale factor (=1=>[cm]).
    Int_t *fidmed; // pointer to array of medium numbers
    AliModule *fits; // local pointer to ITS module needed for AliMixture...

    ClassDef(AliITSBaseGeometry,1) // Basic ITS Geometry class
};

#endif
