#ifndef ALIITSBASEVOLPARAMS_H
#define ALIITSBASEVOLPARAMS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
  $Id$
 */

/////////////////////////////////////////////////////////////////////////
//  Geant 3 base class for volume parameters
/////////////////////////////////////////////////////////////////////////
#include <TObject.h>
class TString;

class AliITSBaseVolParams : public TObject {
 public:
    AliITSBaseVolParams(){fVol=-1;} // Default constructor
    virtual ~AliITSBaseVolParams(){} // Destructor
    virtual void SetVid(Int_t i){fVol = i;}// Sets the volume id from the next
	                                   // available one
    virtual Int_t GetVid(){return fVol;} // Returns volume id
    virtual void SetName(const char *c){fName=c;}//Sets name of this volume
    virtual TString* GetName(){return &fName;} // Returns volume name
    virtual void Print(ostream *os){} // Prints output content of this class
    virtual void Read(istream *is){} // Reads output created by Print above.
 private:
    Int_t   fVol;  // Volume index number
    TString fName; // Volume name

    ClassDef(AliITSBaseVolParams,1) // Basic ITS volume parameters class
};

#endif


#ifndef ALIITSTUBEDATA_H
#define ALIITSTUBEDATA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
  $Id$
 */

/////////////////////////////////////////////////////////////////////////
//  Geant 3 Poly-Cone data structure.
/////////////////////////////////////////////////////////////////////////
#include <math.h> // for the definision of NAN.
class TString;

class AliITSTubeData : public AliITSBaseVolParams {
 public:
    AliITSTubeData() : AliITSBaseVolParams() // Default constructor
	{fDz=0;fRmin=0;fRmax=0;}
    virtual ~AliITSTubeData(){} // Standard destructor
    // Getters
    Double_t DzAt() // Return fDz coordiante
	{return fDz;}
    Double_t Rmin() {return fRmin;} // Returm arrau of rmin values
    Double_t Rmax() {return fRmax;} // Returm arrau of rmax values
    Double_t& Z()// Returns address of fDz
	{return fDz;}
    Double_t& Rn()// Returns address of fRmin
	{return fRmin;}
    Double_t& Rx()// Returns address of fRmax
	{return fRmax;}
    void Print(ostream *os); // Prints output content of this class
    void Read(istream *is); // Reads output created by Print above.
 private:
    Double_t fRmin; // Inner Radius
    Double_t fRmax; // Outer Radius
    Double_t fDz;   // Z half length.

    ClassDef(AliITSTubeData,1) // Poly Cone data class
};
// Input and output function for standard C++ input/output.
ostream &operator<<(ostream &os,AliITSTubeData &source);
istream &operator>>(istream &os,AliITSTubeData &source);


#endif
   
#ifndef ALIITSPCONEDATA_H
#define ALIITSPCONEDATA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
  $Id$
 */

/////////////////////////////////////////////////////////////////////////
//  Geant 3 Poly-Cone data structure.
/////////////////////////////////////////////////////////////////////////
#include <math.h> // for the definision of NAN.
class TString;

class AliITSPConeData : public AliITSBaseVolParams {
 public:
    AliITSPConeData() : AliITSBaseVolParams() // Default constructor
	{fNz=0;fPhi0=0.0;fDphi=0.0;fZ=0;fRmin=0;fRmax=0;}
    AliITSPConeData(Int_t n) : AliITSBaseVolParams() // Standard constructor
	{fNz=n;fPhi0=0.0;fDphi=360.0;fZ=new Double_t[n];
	fRmin=new Double_t[n];fRmax=new Double_t[n];}
    AliITSPConeData(AliITSPConeData &source) // Copy constructor
	{ *this = source;}
    virtual ~AliITSPConeData() // Standard destructor
	{delete[] fZ;delete[] fRmin;delete[] fRmax;fNz=0;}
    AliITSPConeData& operator=(AliITSPConeData &source) // Equals operator
	{this->SetVid(source.GetVid());
	this->SetName((source.GetName())->Data());
	this->fNz = source.fNz;this->fPhi0=source.fPhi0;
	this->fDphi=source.fDphi;
	if(this->fZ!=0) delete[] this->fZ;
	if(this->fRmin!=0) delete[] this->fRmin;
	if(this->fRmax!=0) delete[] this->fRmax;
	this->fZ=0;this->fRmin=0;this->fRmax=0;if(this->fNz<=0) return *this;
	this->fZ=new Double_t[this->fNz];this->fRmin=new Double_t[this->fNz];
	this->fRmax=new Double_t[this->fNz];for(Int_t i=0;i<this->fNz;i++){
	    this->fZ[i]=source.fZ[i];this->fRmin[i]=source.fRmin[i];
	    fRmax[i]=source.fRmax[i];}return *this;}
    void Size(Int_t n)// Sets the number of Z,Rmin,Rmax parameters
	{if(fZ!=0) delete fZ;if(fRmin!=0) delete fRmin;
	if(fRmax!=0) delete fRmax; fNz=n;fPhi0=0.0;fDphi=360.0;
	fZ=new Double_t[n];fRmin=new Double_t[n];fRmax=new Double_t[n];}
    void Size(Int_t n,const char *c){//Sets the number of Z,Rmin, and Rmax
	// parameters as well as the volume name.
	Size(n);SetName(c);}
    // Getters
    Int_t Nz() {return fNz;} // Returns fNz
    Double_t Phi0(){return fPhi0;}// Return starting phi value
    Double_t DPhi(){return fDphi;}// Return delta phi value
    Double_t ZAt(Int_t i) // Return Z coordiante
	{if(i<0||i>=fNz) return NAN;else return fZ[i];}
    Double_t *Z(){return fZ;} // Returns array of z values
    Double_t *Rmin() {return fRmin;} // Returm arrau of rmin values
    Double_t *Rmax() {return fRmax;} // Returm arrau of rmax values
    Double_t Rmin(Int_t i) // Return Inner radius value
	{if(i<0||i>=fNz) return NAN;else return fRmin[i];}
    Double_t Rmax(Int_t i) // Return Outer radius value
	{if(i<0||i>=fNz) return NAN;else return fRmax[i];}
    // Setters
    Double_t& P0() // Returns the address of fPhi0
	{return fPhi0;}
    Double_t& dP() // Returns the address of fDphi
	{return fDphi;}
    Double_t& Z(Int_t i)// Returns address of fZ
	{/*if(i<0||i>=fNz) return 0;*/return fZ[i];}
    Double_t& Rn(Int_t i)// Returns address of fRmin
	{/*if(i<0||i>=fNz) return 0;*/return fRmin[i];}
    Double_t& Rx(Int_t i)// Returns address of fRmax
	{/*if(i<0||i>=fNz) return 0;*/return fRmax[i];}
    void Print(ostream *os); // Prints output content of this class
    void Read(istream *is); // Reads output created by Print above.
 private:
    Int_t fNz; // Number of z sections
    Double_t fPhi0,fDphi; // Starting phi angle and delta phi [degrees]
    Double_t *fZ;    //[n] Z coordiantes
    Double_t *fRmin; //[n] Inner radius
    Double_t *fRmax; //[n] Outer radius

    ClassDef(AliITSPConeData,1) // Poly Cone data class
};
// Input and output function for standard C++ input/output.
ostream &operator<<(ostream &os,AliITSPConeData &source);
istream &operator>>(istream &os,AliITSPConeData &source);


#endif
   

#ifndef ALIITSPGONDATA_H
#define ALIITSPGONDATA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
  $Id$
 */

/////////////////////////////////////////////////////////////////////////
//  Geant 3 Poly-Gon data structure.
/////////////////////////////////////////////////////////////////////////
#include <math.h> // for the definision of NAN.
class TString;

class AliITSPGonData : public AliITSBaseVolParams {
 public:
    AliITSPGonData() : AliITSBaseVolParams() // Default constructor
	{fNz=0;fNphi=0;fPhi0=0.0;fDphi=0.0;fZ=0;fRmin=0;fRmax=0;}
    AliITSPGonData(Int_t n) : AliITSBaseVolParams() // Standard constructor
	{fNz=n;fNphi=0;fPhi0=0.0;fDphi=360.0;fZ=new Double_t[n];
	fRmin=new Double_t[n];fRmax=new Double_t[n];}
    AliITSPGonData(AliITSPGonData &source) // Copy constructor
	{ *this = source;}
    virtual ~AliITSPGonData() // Standard destructor
	{delete[] fZ;delete[] fRmin;delete[] fRmax;fNz=0;}
    AliITSPGonData& operator=(AliITSPGonData &source) // Equals operator
	{this->SetVid(source.GetVid());
	this->SetName((source.GetName())->Data());
	this->fNz = source.fNz;this->fNphi=source.fNphi;
	this->fPhi0=source.fPhi0;this->fDphi=source.fDphi;
	if(this->fZ!=0) delete[] this->fZ;
	if(this->fRmin!=0) delete[] this->fRmin;
	if(this->fRmax!=0) delete[] this->fRmax;
	this->fZ=0;this->fRmin=0;this->fRmax=0;if(this->fNz<=0) return *this;
	this->fZ=new Double_t[this->fNz];this->fRmin=new Double_t[this->fNz];
	this->fRmax=new Double_t[this->fNz];for(Int_t i=0;i<this->fNz;i++){
	    this->fZ[i]=source.fZ[i];this->fRmin[i]=source.fRmin[i];
	    fRmax[i]=source.fRmax[i];}return *this;}
    void Size(Int_t n)// Sets the number of Z,Rmin,Rmax parameters
	{if(fZ!=0) delete fZ;if(fRmin!=0) delete fRmin;
	if(fRmax!=0) delete fRmax; fNz=n;fPhi0=0.0;fDphi=360.0;
	fZ=new Double_t[n];fRmin=new Double_t[n];fRmax=new Double_t[n];}
    void Size(Int_t n,const char *c){//Sets the number of Z,Rmin, and Rmax
	// parameters as well as the volume name.
	Size(n);SetName(c);}
    // Getters
    Int_t NPhi() {return fNz;} // Returns fNphi
    Int_t Nz() {return fNz;} // Returns fNz
    Double_t Phi0(){return fPhi0;}// Return starting phi value
    Double_t DPhi(){return fDphi;}// Return delta phi value
    Double_t *Z(){return fZ;} // Returns array of z values
    Double_t *Rmin() {return fRmin;} // Returm arrau of rmin values
    Double_t *Rmax() {return fRmax;} // Returm arrau of rmax values
    Double_t ZAt(Int_t i) // Return Z coordiante
	{if(i<0||i>=fNz) return NAN;else return fZ[i];}
    Double_t Rmin(Int_t i) // Return Inner radius value
	{if(i<0||i>=fNz) return NAN;else return fRmin[i];}
    Double_t Rmax(Int_t i) // Return Outer radius value
	{if(i<0||i>=fNz) return NAN;else return fRmax[i];}
    // Setters
    void Nphi(Int_t i) {fNphi = i;} // Sets fNphi
    Double_t& P0() // Returns the address of fPhi0
	{return fPhi0;}
    Double_t& dP() // Returns the address of fDphi
	{return fDphi;}
    Double_t& Z(Int_t i)// Returns address of fZ
	{/*if(i<0||i>=fNz) return 0;*/return fZ[i];}
    Double_t& Rn(Int_t i)// Returns address of fRmin
	{/*if(i<0||i>=fNz) return 0;*/return fRmin[i];}
    Double_t& Rx(Int_t i)// Returns address of fRmax
	{/*if(i<0||i>=fNz) return 0;*/return fRmax[i];}
    void Print(ostream *os); // Prints output content of this class
    void Read(istream *is); // Reads output created by Print above.
 private:
    Int_t fNphi;  // Number of sections in phi.
    Int_t fNz;    // Number of Z sections
    Double_t fPhi0,fDphi; // Starting phi angle and delta phi [degrees]
    Double_t *fZ;    //[n] Z coordiantes
    Double_t *fRmin; //[n] Inner radius
    Double_t *fRmax; //[n] Outer radius

    ClassDef(AliITSPGonData,1) // Poly Gon Data Class
};
// Input and output function for standard C++ input/output.
ostream &operator<<(ostream &os,AliITSPGonData &source);
istream &operator>>(istream &os,AliITSPGonData &source);

#endif
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
#include "AliModule.h"
class TString;
class TVector3;

class AliITSBaseGeometry : public TObject {
 public:
    AliITSBaseGeometry(); // Default constructor
    AliITSBaseGeometry(AliModule *its,Int_t iflag); // Standard Constructor
    virtual ~AliITSBaseGeometry(); // Destructor
    virtual void BuildDisplayGeometry(){}; // Calls ROOT geometry interface
                                      // to AliRoot display
     // Calls Geant3 interface geometry routines
    virtual void CreateG3Geometry(){};
    virtual void CreateG3Materials(){};//Calls Geant3 interface for materials
    virtual Int_t IsVersion() const{return 11;}// return version of geometry.
    // Get Index for Geant3 v name
    Int_t ITSG3VnameToIndex(const char name[3])const; 
    Int_t ITSIndexToITSG3name(const Int_t i); // Get Geant3 volume name
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
    virtual void Box(const char *gnam,const TString &dis,
             Double_t dx,Double_t dy,Double_t dz,Int_t med);
    // Greate A Trapizoid with the x dimension varing along z.
    virtual void Trapezoid1(const char *gnam,const TString &dis,Double_t dxn,
                    Double_t dxp,Double_t dy,Double_t dz,Int_t med);
    // Greate A Trapizoid with the x and y dimension varing along z.
    virtual void Trapezoid2(const char *gnam,const TString &dis,Double_t dxn,
                    Double_t dxp,Double_t dyn,Double_t dyp,Double_t dz,
                    Int_t med);
    // General trapazoid.
    virtual void Trapezoid(const char *gnam,const TString &dis,Double_t dz,
                   Double_t thet,Double_t phi,Double_t h1,Double_t bl1,
                   Double_t tl1,Double_t alp1,Double_t h2,Double_t bl2,
                   Double_t tl2,Double_t alp2,Int_t med);
    // Simple Tube.
    virtual void Tube(const char *gnam,const TString &dis,Double_t rmin,
              Double_t rmax,Double_t dz,Int_t med);
    virtual void Tube(AliITSTubeData &d,Int_t med);
    // Tube segment.
    virtual void TubeSegment(const char *gnam,const TString &dis,Double_t rmin,
                     Double_t rmax,Double_t dz,Double_t phi1,Double_t phi2,
                     Int_t med);
    // Simple Cone.
    virtual void Cone(const char *gnam,const TString &dis,Double_t dz,
	      Double_t rmin1,Double_t rmax1,Double_t rmin2,Double_t rmax2,
	      Int_t med);
    // Segment of a Cone.
    virtual void ConeSegment(const char *gnam,const TString &dis,Double_t dz,
                     Double_t rmin1,Double_t rmax1,Double_t rmin2,
                     Double_t rmax2,Double_t phi1,Double_t phi2,Int_t med);
    // Spherical shell segment.
    virtual void Sphere(const char *gnam,const TString &dis,Double_t rmin,
                Double_t rmax,Double_t the1,Double_t the2,Double_t phi1,
                Double_t phi2,Int_t med);
    // Parallelepiped.
    virtual void Parallelepiped(const char *gnam,const TString &dis,
				Double_t dx,Double_t dy,Double_t dz,
				Double_t alph,Double_t thet,
				Double_t phi,Int_t med);
    // Polygon.
    virtual void PolyGon(const char *gnam,const TString &dis,Double_t phi1,
                 Double_t dphi,Int_t npdv,Int_t nz,Double_t *z,Double_t *rmin,
                 Double_t *rmax,Int_t med);
    virtual void PolyGon(const char *gnam,const TString &dis,AliITSPGonData &d,
		  Int_t med){PolyGon(gnam,dis,d.Phi0(),d.DPhi(),d.NPhi(),
				      d.Nz(),d.Z(),d.Rmin(),d.Rmax(),med);}
    virtual void PolyGon(AliITSPGonData &d,Int_t med);
    //Poly-Cone
    virtual void PolyCone(const char *gnam,const TString &dis,Double_t phi1,
                  Double_t dphi,Int_t nz,Double_t *z,Double_t *rmin,
                  Double_t *rmax,Int_t med);
    virtual void PolyCone(const char *gnam,const TString &dis,
			  AliITSPConeData &d,
			  Int_t med){PolyCone(gnam,dis,d.Phi0(),
					      d.DPhi(),d.Nz(),d.Z(),
					      d.Rmin(),d.Rmax(),med);}
    virtual void PolyCone(AliITSPConeData &d,Int_t med);
    // Ellliptical cross-sectino tube
    virtual void TubeElliptical(const char *gnam,const TString &dis,
				Double_t p1,Double_t p2,Double_t dz,Int_t med);
    // Hyperbolic tube
    virtual void HyperbolicTube(const char *gnam,const TString &dis,
				Double_t rmin,Double_t rmax,Double_t dz,
				Double_t thet,Int_t med);
    // Twisted genral trapezoid.
    virtual void TwistedTrapezoid(const char *gnam,const TString &dis,
				  Double_t dz,
				  Double_t thet,Double_t phi,Double_t twist,
				  Double_t h1,Double_t bl1,Double_t tl1,
				  Double_t apl1,Double_t h2,Double_t bl2,
				  Double_t tl2,Double_t apl2,Int_t med);
    // Cut tube.
    virtual void CutTube(const char *gnam,const TString &dis,Double_t rmin,
                 Double_t rmax,Double_t dz,Double_t phi1,Double_t phi2,
                 Double_t lx,Double_t ly,Double_t lz,Double_t hx,Double_t hy,
                 Double_t hz,Int_t med);
    // Position one volume inside another
    virtual void Pos(const char *vol,Int_t cn,const char *moth,Double_t x,
             Double_t y,Double_t z,Int_t irot);
    // Position one volume inside another
    virtual void Pos(AliITSBaseVolParams &v,Int_t cn,
	     AliITSBaseVolParams &m,TVector3 &t,Int_t irot);
    Int_t GetMed(Int_t med){return (fits->GetIdtmed())->At(med);}
    // Define rotation matrix
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
    void Element(Int_t imat,const char *name,Int_t z,Double_t den,Int_t istd);
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
    void G3name(const char *gname,char *name)//Add's I to name and ending null
	{for(Int_t i=0;i<3;i++) name[i+1] = gname[i];
	name[0]='I';name[4]='\0';}
    //
 protected:
    static Int_t fNCreates; //! Counts the number of time this class has
    // been created.
    static const Double_t fAlpha = 7.297352533e-3; //! find structure constant
    static const Double_t fRe = 2.81794028e-13; //![cm]classical elect. radius
    static const Double_t fNa = 6.02214199e+23; //! [#/mole] Avogadro's number
    static Int_t *fidrot;
    static Int_t fidrotsize;
    static Int_t fidrotlast;
    static TString *fVolName; // Array of ITS Volumen names.
    static Int_t fVolNameSize; // Size of Array fVolName
    static Int_t fVolNameLast; // Last filled element of fVolName
    Double_t fScale; // Scale factor (=1=>[cm]).
    AliModule *fits; // local pointer to ITS module needed for AliMixture...

    ClassDef(AliITSBaseGeometry,1) // Basic ITS Geometry class
};

#endif
     
