//
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
    AliITSBaseVolParams(){fVol=-1;fCpn=-1;} // Default constructor
    virtual ~AliITSBaseVolParams(){} // Destructor
    virtual void SetVid(Int_t i){fVol = i;}// Sets the volume id from the next
	                                   // available one
    virtual Int_t GetVid(){return fVol;} // Returns volume id
    virtual void SetName(const char *c){fName=c;}//Sets name of this volume
    virtual TString* GetName(){return &fName;} // Returns volume name
    virtual void Print(ostream *os); // Prints output content of this class
    virtual void Read(istream *is); // Reads output created by Print above.
    virtual Int_t GetCopyNumber(){return fCpn;}//Returns existing max copy no.
    // incoment and return copy number.
    virtual Int_t GetNextcpn(){fCpn++;return fCpn;}
    // Same as above, but add's 1 for FORTRAN/Geant3 compliance.
    virtual Int_t GetG3cpn(){return GetNextcpn()+1;}
    /*
    Int_t ITSIndexToITSG3name(const Int_t i)
	{return AliITSBaseGeometry::ITSIndexToITSG3name(i);};
    Int_t ITSG3VnameToIndex(const char *name)const
	{return AliITSBaseGeometry::ITSG3VnameToIndex(name);};
    */
 private:
    Int_t   fVol;  // Volume index number
    Int_t   fCpn;  // max Copy number
    TString fName; // Volume name

    ClassDef(AliITSBaseVolParams,1) // Basic ITS volume parameters class
};
// Input and output function for standard C++ input/output.
ostream &operator<<(ostream &os,AliITSBaseVolParams &source);
istream &operator>>(istream &os,AliITSBaseVolParams &source);
#endif

#ifndef ALIITSBOXDATA_H
#define ALIITSBOXDATA_H
/////////////////////////////////////////////////////////////////////////
//  Geant 3 Box data structure.
/////////////////////////////////////////////////////////////////////////

class AliITSBoxData : public AliITSBaseVolParams {
 public:
    AliITSBoxData() : AliITSBaseVolParams() // Default constructor
	{fDx=0.0;fDy=0.0;fDz=0;}
    virtual ~AliITSBoxData(){} // Standard destructor
    // Getters
    Double_t DxAt() {return fDx;} // Returm arrau of rmin values
    Double_t DyAt() {return fDy;} // Returm arrau of rmax values
    Double_t DzAt() {return fDz;} // Return fDz coordiante
    Double_t& Dx(){return fDx;}// Returns address of fRmax
    Double_t& Dy(){return fDy;}// Returns address of fRmin
    Double_t& Dz(){return fDz;}// Returns address of fDz
    void Print(ostream *os); // Prints output content of this class
    void Read(istream *is); // Reads output created by Print above.
 private:
    Double_t fDx; // X half length
    Double_t fDy; // Y half length
    Double_t fDz; // Z half length.

    ClassDef(AliITSBoxData,1) // Box data class
};
// Input and output function for standard C++ input/output.
ostream &operator<<(ostream &os,AliITSBoxData &source);
istream &operator>>(istream &os,AliITSBoxData &source);
#endif

#ifndef ALIITSTRAPEZOID1DATA_H
#define ALIITSTRAPEZOID1DATA_H
/////////////////////////////////////////////////////////////////////////
//  Geant 3 Trapezoid 1 data structure.
/////////////////////////////////////////////////////////////////////////

class AliITSTrapezoid1Data : public AliITSBaseVolParams {
 public:
    AliITSTrapezoid1Data() : AliITSBaseVolParams() // Default constructor
	{fDx[0]=0.0;fDx[1]=0.0;fDy=0.0;fDz=0;}
    virtual ~AliITSTrapezoid1Data(){} // Standard destructor
    // Getters
    Double_t DxAt(Int_t i) {return fDx[i];} // Returm arrau of rmin values
    Double_t DyAt() {return fDy;} // Returm arrau of rmax values
    Double_t DzAt() {return fDz;} // Return fDz coordiante
    Double_t& Dx(Int_t i){return fDx[i];}// Returns address of fRmax
    Double_t& Dy(){return fDy;}// Returns address of fRmin
    Double_t& Dz(){return fDz;}// Returns address of fDz
    void Print(ostream *os); // Prints output content of this class
    void Read(istream *is); // Reads output created by Print above.
 private:
    Double_t fDx[2]; // X half length
    Double_t fDy; // Y half length
    Double_t fDz; // Z half length.

    ClassDef(AliITSTrapezoid1Data,1) // Trapezoid 1 data class
};
// Input and output function for standard C++ input/output.
ostream &operator<<(ostream &os,AliITSTrapezoid1Data &source);
istream &operator>>(istream &os,AliITSTrapezoid1Data &source);
#endif

#ifndef ALIITSTRAPEZOID2DATA_H
#define ALIITSTRAPEZOID2DATA_H
/////////////////////////////////////////////////////////////////////////
//  Geant 3 Trapezoid 2 data structure.
/////////////////////////////////////////////////////////////////////////

class AliITSTrapezoid2Data : public AliITSBaseVolParams {
 public:
    AliITSTrapezoid2Data() : AliITSBaseVolParams() // Default constructor
	{fDx[0]=0.0;fDx[1]=0.0;fDy[0]=0.0;fDy[1]=0.0;fDz=0.0;}
    virtual ~AliITSTrapezoid2Data(){} // Standard destructor
    // Getters
    Double_t DxAt(Int_t i) {return fDx[i];} // Returm arrau of rmin values
    Double_t DyAt(Int_t i) {return fDy[i];} // Returm arrau of rmax values
    Double_t DzAt() {return fDz;} // Return fDz coordiante
    Double_t& Dx(Int_t i){return fDx[i];}// Returns address of fRmax
    Double_t& Dy(Int_t i){return fDy[i];}// Returns address of fRmin
    Double_t& Dz(){return fDz;}// Returns address of fDz
    void Print(ostream *os); // Prints output content of this class
    void Read(istream *is); // Reads output created by Print above.
 private:
    Double_t fDx[2]; // X half length
    Double_t fDy[2]; // Y half length
    Double_t fDz; // Z half length.

    ClassDef(AliITSTrapezoid2Data,1) // Trapezoid 2 data class
};
// Input and output function for standard C++ input/output.
ostream &operator<<(ostream &os,AliITSTrapezoid2Data &source);
istream &operator>>(istream &os,AliITSTrapezoid2Data &source);
#endif

#ifndef ALIITSTRAPEZOIDDATA_H
#define ALIITSTRAPEZOIDDATA_H
/////////////////////////////////////////////////////////////////////////
//  Geant 3 Trapezoid General data structure.
/////////////////////////////////////////////////////////////////////////;

class AliITSTrapezoidData : public AliITSBaseVolParams {
 public:
    AliITSTrapezoidData() : AliITSBaseVolParams() // Default constructor
	{fDz=0.0;fTheta=0.0;fPhi=0.0;fH[0]=0.0;fH[1]=0.0;fBl[0]=0.0;
	fBl[1]=0.0;fTl[0]=0.0;fTl[1]=0.0;fAlp[0]=0.0;fAlp[1]=0.0;}
    virtual ~AliITSTrapezoidData(){} // Standard destructor
    // Getters
    Double_t DzAt() {return fDz;} // Returm arrau of dZ values
    Double_t Theta() {return fTheta;} // Returm arrau of Theta values
    Double_t Phi() {return fPhi;} // Return Phi coordiante
    Double_t HAt(Int_t i){return fH[i];} // Return fH[]
    Double_t Bl(Int_t i){return fBl[i];} // Return fBl[]
    Double_t Tl(Int_t i){return fTl[i];} // Return fTl[]
    Double_t Alpha(Int_t i){return fAlp[i];} // Return fAlpha[]
    Double_t& Dz(){return fDz;}// Returns address of fDz
    Double_t& Th(){return fTheta;}// Returns address of fTheta
    Double_t& Ph(){return fPhi;}// Returns address of fPhi
    Double_t& H(Int_t i){return fH[i];}// Returns address of fH
    Double_t& B(Int_t i){return fBl[i];}// Returns address of fBl
    Double_t& T(Int_t i){return fTl[i];}// Returns address of fTl
    Double_t& A(Int_t i){return fAlp[i];}// Returns address of fAlp
    void Print(ostream *os); // Prints output content of this class
    void Read(istream *is); // Reads output created by Print above.
 private:
    Double_t fTheta; // Polar angle of the line joing the center of the 
                     // face at -dz to the center of the one at dz [degree].
    Double_t fPhi;   // aximuthal angle of the line joing the center of the
                     //  face at -dz to the center of the one at +dz [degree].
    Double_t fDz;    // Half-length along the z-asix
    Double_t fH[2];  // half-length along y of the face at -dz & _dz.
    Double_t fBl[2]; // half-length along x of the side at -h1 in y of 
                     // the face at -dz & +dz in z.
    Double_t fTl[2]; // half-length along x of teh side at +h1 in y of 
                     // the face at -dz & +dz in z.
    Double_t fAlp[2];// angle with respect to the y axis from the center of
                     // the side at -fH[0] & +fH[1] in y to the cetner of the 
                     // side at +h1 in y of the face at -dz in z [degree].

    ClassDef(AliITSTrapezoidData,1) // General Trapezoid data class
};
// Input and output function for standard C++ input/output.
ostream &operator<<(ostream &os,AliITSTrapezoidData &source);
istream &operator>>(istream &os,AliITSTrapezoidData &source);
#endif

#ifndef ALIITSTRAPEZOIDTWISTEDDATA_H
#define ALIITSTRAPEZOIDTWISTEDDATA_H
/////////////////////////////////////////////////////////////////////////
//  Geant 3 Trapezoid General data structure.
/////////////////////////////////////////////////////////////////////////;

class AliITSTrapezoidTwistedData : public AliITSBaseVolParams {
 public:
    AliITSTrapezoidTwistedData() : AliITSBaseVolParams() //Default constructor
	{fTwist=0.0;fDz=0.0;fTheta=0.0;fPhi=0.0;fH[0]=0.0;fH[1]=0.0;
	fBl[0]=0.0;fBl[1]=0.0;fTl[0]=0.0;fTl[1]=0.0;fAlp[0]=0.0;fAlp[1]=0.0;}
    virtual ~AliITSTrapezoidTwistedData(){} // Standard destructor
    // Getters
    Double_t DzAt() {return fDz;} // Returm arrau of dZ values
    Double_t Theta() {return fTheta;} // Returm arrau of Theta values
    Double_t Phi() {return fPhi;} // Return Phi angle
    Double_t Twist(){return fTwist;} // Returns Twist angle
    Double_t HAt(Int_t i){return fH[i];} // Return fH[]
    Double_t Bl(Int_t i){return fBl[i];} // Return fBl[]
    Double_t Tl(Int_t i){return fTl[i];} // Return fTl[]
    Double_t Alpha(Int_t i){return fAlp[i];} // Return fAlpha[]
    Double_t& Dz(){return fDz;}// Returns address of fDz
    Double_t& Th(){return fTheta;}// Returns address of fTheta
    Double_t& Ph(){return fPhi;}// Returns address of fPhi
    Double_t& Tw(){return fTwist;}// Returns address of fTwist
    Double_t& H(Int_t i){return fH[i];}// Returns address of fH
    Double_t& B(Int_t i){return fBl[i];}// Returns address of fBl
    Double_t& T(Int_t i){return fTl[i];}// Returns address of fTl
    Double_t& A(Int_t i){return fAlp[i];}// Returns address of fAlp
    void Print(ostream *os); // Prints output content of this class
    void Read(istream *is); // Reads output created by Print above.
 private:
    Double_t fTwist; // Twist angle of the faces parallel to the x-y 
                     // plane at z=+-dz around an axis parallel to z
    Double_t fTheta; // Polar angle of the line joing the center of the 
                     // face at -dz to the center of the one at dz [degree].
    Double_t fPhi;   // aximuthal angle of the line joing the center of the
                     //  face at -dz to the center of the one at +dz [degree].
    Double_t fDz;    // Half-length along the z-asix
    Double_t fH[2];  // half-length along y of the face at -dz & _dz.
    Double_t fBl[2]; // half-length along x of the side at -h1 in y of 
                     // the face at -dz & +dz in z.
    Double_t fTl[2]; // half-length along x of teh side at +h1 in y of 
                     // the face at -dz & +dz in z.
    Double_t fAlp[2];// angle with respect to the y axis from the center of
                     // the side at -fH[0] & +fH[1] in y to the cetner of the 
                     // side at +h1 in y of the face at -dz in z [degree].

    ClassDef(AliITSTrapezoidTwistedData,1) // Twisted Trapezoid data class
};
// Input and output function for standard C++ input/output.
ostream &operator<<(ostream &os,AliITSTrapezoidTwistedData &source);
istream &operator>>(istream &os,AliITSTrapezoidTwistedData &source);
#endif

#ifndef ALIITSTUBEDATA_H
#define ALIITSTUBEDATA_H
/////////////////////////////////////////////////////////////////////////
//  Geant 3 Tube data structure.
/////////////////////////////////////////////////////////////////////////

class AliITSTubeData : public AliITSBaseVolParams {
 public:
    AliITSTubeData() : AliITSBaseVolParams() // Default constructor
	{fDz=0;fRmin=0;fRmax=0;}
    virtual ~AliITSTubeData(){} // Standard destructor
    // Getters
    Double_t DzAt(){return fDz;} // Return fDz coordiante
    Double_t Rmin() {return fRmin;} // Returm arrau of rmin values
    Double_t Rmax() {return fRmax;} // Returm arrau of rmax values
    Double_t& Z(){return fDz;}// Returns address of fDz
    Double_t& Rn(){return fRmin;}// Returns address of fRmin
    Double_t& Rx(){return fRmax;}// Returns address of fRmax
    void Print(ostream *os); // Prints output content of this class
    void Read(istream *is); // Reads output created by Print above.
 private:
    Double_t fRmin; // Inner Radius
    Double_t fRmax; // Outer Radius
    Double_t fDz;   // Z half length.

    ClassDef(AliITSTubeData,1) // Simple Tube data class
};
// Input and output function for standard C++ input/output.
ostream &operator<<(ostream &os,AliITSTubeData &source);
istream &operator>>(istream &os,AliITSTubeData &source);
#endif

#ifndef ALIITSTUBESEGDATA_H
#define ALIITSTUBESEGDATA_H
/////////////////////////////////////////////////////////////////////////
//  Geant 3 Tube Segment data structure.
/////////////////////////////////////////////////////////////////////////

class AliITSTubeSegData : public AliITSBaseVolParams {
 public:
    AliITSTubeSegData() : AliITSBaseVolParams() // Default constructor
	{fDz=0;fRmin=0;fRmax=0;}
    virtual ~AliITSTubeSegData(){} // Standard destructor
    // Getters
    Double_t Phi0() {return fPhi0;} // Returns starting Phi value
    Double_t Phi1() {return fPhi1;} // Returns endinging Phi value
    Double_t DzAt(){return fDz;} // Return fDz coordiante
    Double_t Rmin() {return fRmin;} // Returm arrau of rmin values
    Double_t Rmax() {return fRmax;} // Returm arrau of rmax values
    Double_t& P0(){return fPhi0;}// Returns addres of fPhi0
    Double_t& P1(){return fPhi1;}// Returns addres of fPhi1
    Double_t& Z(){return fDz;}// Returns address of fDz
    Double_t& Rn(){return fRmin;}// Returns address of fRmin
    Double_t& Rx(){return fRmax;}// Returns address of fRmax
    void Print(ostream *os); // Prints output content of this class
    void Read(istream *is); // Reads output created by Print above.
 private:
    Double_t fPhi0,fPhi1;   // Starting and ending phi angles [degrees]
    Double_t fRmin; // Inner Radius
    Double_t fRmax; // Outer Radius
    Double_t fDz;   // Z half length.

    ClassDef(AliITSTubeSegData,1) // Segment of a Tube data class
};
// Input and output function for standard C++ input/output.
ostream &operator<<(ostream &os,AliITSTubeSegData &source);
istream &operator>>(istream &os,AliITSTubeSegData &source);
#endif

#ifndef ALIITSTUBECUTDATA_H
#define ALIITSTUBECUTDATA_H

#include <TVector3.h>

/////////////////////////////////////////////////////////////////////////
//  Geant 3 Tube Cut data structure.
/////////////////////////////////////////////////////////////////////////

class AliITSTubeCutData : public AliITSBaseVolParams {
 public:
    AliITSTubeCutData() : AliITSBaseVolParams() // Default constructor
	{fDz=0;fRmin=0;fRmax=0;fPhi0=0.0,fPhi1=0.0;}
    virtual ~AliITSTubeCutData(){} // Standard destructor
    // Getters
    Double_t Phi0() {return fPhi0;} // Returns starting Phi value
    Double_t Phi1() {return fPhi1;} // Returns endinging Phi value
    Double_t DzAt(){return fDz;} // Return fDz coordiante
    Double_t Rmin() {return fRmin;} // Returm arrau of rmin values
    Double_t Rmax() {return fRmax;} // Returm arrau of rmax values
    TVector3 Normal(Int_t i){return (fNorm[i]).Unit();}// Returns unit normal
                                                     //at z serface -dz & +dz.
    // Returns jth unit normal at z serface -dz & +dz.
    Double_t Normal(Int_t i,Int_t j){return ((fNorm[i]).Unit())[j];}
    Double_t& P0(){return fPhi0;}// Returns addres of fPhi0
    Double_t& P1(){return fPhi1;}// Returns addres of fPhi1
    Double_t& Z(){return fDz;}// Returns address of fDz
    Double_t& Rn(){return fRmin;}// Returns address of fRmin
    Double_t& Rx(){return fRmax;}// Returns address of fRmax
    // Returns address of normal at z serface -dz & +dz.
    TVector3& Nv(Int_t i){return fNorm[i];}
    Double_t& Nx(Int_t i){return (fNorm[i])[0];}// return address of x
                                                //component of the nomal vect.
    Double_t& Ny(Int_t i){return (fNorm[i])[1];}// return address of y
                                                //component of the nomal vect.
    Double_t& Nz(Int_t i){return (fNorm[i])[2];}// return address of z
                                                //component of the nomal vect.
    void Print(ostream *os); // Prints output content of this class
    void Read(istream *is); // Reads output created by Print above.
 private:
    Double_t fPhi0,fPhi1;   // Starting and ending phi angles [degrees]
    Double_t fRmin; // Inner Radius
    Double_t fRmax; // Outer Radius
    Double_t fDz;   // Z half length.
    TVector3 fNorm[2]; // unit vector normal to -dz and +dz surfaces.

    ClassDef(AliITSTubeCutData,1) // A tube segment with cut ends data class
};
// Input and output function for standard C++ input/output.
ostream &operator<<(ostream &os,AliITSTubeCutData &source);
istream &operator>>(istream &os,AliITSTubeCutData &source);
#endif

#ifndef ALIITSTUBEELLIPTICALDATA_H
#define ALIITSTUBEELLIPTICALDATA_H

/////////////////////////////////////////////////////////////////////////
//  Geant 3 Tube Elliptical data structure.
/////////////////////////////////////////////////////////////////////////

class AliITSTubeEllipticalData : public AliITSBaseVolParams {
 public:
    AliITSTubeEllipticalData() : AliITSBaseVolParams() // Default constructor
	{fDz=0;fP0=0;fP1=0;}
    virtual ~AliITSTubeEllipticalData(){} // Standard destructor
    // Getters
    Double_t DzAt(){return fDz;} // Return fDz coordiante
    Double_t P0At() {return fP0;} // Returm fP0 values
    Double_t P1At() {return fP1;} // Returm fP1 values
    Double_t& Z(){return fDz;}// Returns address of fDz
    Double_t& P0(){return fP0;}// Returns address of fP0
    Double_t& P1(){return fP1;}// Returns address of fP1
    void Print(ostream *os); // Prints output content of this class
    void Read(istream *is); // Reads output created by Print above.
 private:
    Double_t fP0; // semi-axis of the elipse along x
    Double_t fP1; // semi-asis of the elipse along y
    Double_t fDz;   // Z half length.

    ClassDef(AliITSTubeEllipticalData,1) // Tube with an elliptical cross
	                                 // section data class
};
// Input and output function for standard C++ input/output.
ostream &operator<<(ostream &os,AliITSTubeEllipticalData &source);
istream &operator>>(istream &os,AliITSTubeEllipticalData &source);
#endif

#ifndef ALIITSTUBEHYPERBOLICDATA_H
#define ALIITSTUBEHYPERBOLICDATA_H

/////////////////////////////////////////////////////////////////////////
//  Geant 3 Tube Hyperbolic data structure.
/////////////////////////////////////////////////////////////////////////

class AliITSTubeHyperbolicData : public AliITSBaseVolParams {
 public:
    AliITSTubeHyperbolicData() : AliITSBaseVolParams() // Default constructor
	{fTheta=0.0;fDz=0;fRmin=0;fRmax=0;}
    virtual ~AliITSTubeHyperbolicData(){} // Standard destructor
    // Getters
    Double_t ThetaAt() {return fTheta;} // Returns starting Theta value
    Double_t DzAt(){return fDz;} // Return fDz coordiante
    Double_t Rmin() {return fRmin;} // Returm arrau of rmin values
    Double_t Rmax() {return fRmax;} // Returm arrau of rmax values
    Double_t& Theta(){return fTheta;}// Returns addres of fTheta
    Double_t& Z(){return fDz;}// Returns address of fDz
    Double_t& Rn(){return fRmin;}// Returns address of fRmin
    Double_t& Rx(){return fRmax;}// Returns address of fRmax
    void Print(ostream *os); // Prints output content of this class
    void Read(istream *is); // Reads output created by Print above.
 private:
    Double_t fTheta;   // stero angle of rotation of the tetwo faces [degrees]
    Double_t fRmin; // Inner Radius
    Double_t fRmax; // Outer Radius
    Double_t fDz;   // Z half length.

    ClassDef(AliITSTubeHyperbolicData,1) // Hyperbolic Tube data class
};
// Input and output function for standard C++ input/output.
ostream &operator<<(ostream &os,AliITSTubeHyperbolicData &source);
istream &operator>>(istream &os,AliITSTubeHyperbolicData &source);
#endif

#ifndef ALIITSCONEDATA_H
#define ALIITSCONEDATA_H
/////////////////////////////////////////////////////////////////////////
//  Geant 3 Cone data structure.
/////////////////////////////////////////////////////////////////////////

class AliITSConeData : public AliITSBaseVolParams {
 public:
    AliITSConeData() : AliITSBaseVolParams() // Default constructor
	{fDz=0;fRmin0=0;fRmax0=0;fRmin1=0;fRmax1=0;}
    virtual ~AliITSConeData(){} // Standard destructor
    // Getters
    Double_t DzAt(){return fDz;} // Return fDz coordiante
    Double_t Rmin0() {return fRmin0;} // Returm arrau of rmin0 values
    Double_t Rmax0() {return fRmax0;} // Returm arrau of rmax0 values
    Double_t Rmin1() {return fRmin1;} // Returm arrau of rmin1 values
    Double_t Rmax1() {return fRmax1;} // Returm arrau of rmax1 values
    Double_t& Z(){return fDz;}// Returns address of fDz
    Double_t& Rn0(){return fRmin0;}// Returns address of fRmin0
    Double_t& Rx0(){return fRmax0;}// Returns address of fRmax0
    Double_t& Rn1(){return fRmin1;}// Returns address of fRmin1
    Double_t& Rx1(){return fRmax1;}// Returns address of fRmax1
    void Print(ostream *os); // Prints output content of this class
    void Read(istream *is); // Reads output created by Print above.
 private:
    Double_t fRmin0,fRmin1; // Inner Radius
    Double_t fRmax0,fRmax1; // Outer Radius
    Double_t fDz;   // Z half length.

    ClassDef(AliITSConeData,1) // Simple Cone data class
};
// Input and output function for standard C++ input/output.
ostream &operator<<(ostream &os,AliITSConeData &source);
istream &operator>>(istream &os,AliITSConeData &source);
#endif
  
#ifndef ALIITSCONESEGDATA_H
#define ALIITSCONESEGDATA_H
/////////////////////////////////////////////////////////////////////////
//  Geant 3 Cone Segment data structure.
/////////////////////////////////////////////////////////////////////////

class AliITSConeSegData : public AliITSBaseVolParams {
 public:
    AliITSConeSegData() : AliITSBaseVolParams() // Default constructor
	{fPhi0=0.0;fPhi1=360.0;fDz=0;fRmin0=0;fRmax0=0;fRmin1=0;fRmax1=0;}
    virtual ~AliITSConeSegData(){} // Standard destructor
    // Getters
    Double_t Phi0() {return fPhi0;} // Returns starting Phi value
    Double_t Phi1() {return fPhi1;} // Returns endinging Phi value
    Double_t DzAt() {return fDz;} // Return fDz coordiante
    Double_t Rmin0() {return fRmin0;} // Returm arrau of rmin0 values
    Double_t Rmax0() {return fRmax0;} // Returm arrau of rmax0 values
    Double_t Rmin1() {return fRmin1;} // Returm arrau of rmin1 values
    Double_t Rmax1() {return fRmax1;} // Returm arrau of rmax1 values
    Double_t& P0(){return fPhi0;}// Returns addres of fPhi0
    Double_t& P1(){return fPhi1;}// Returns addres of fPhi1
    Double_t& Z(){return fDz;}// Returns address of fDz
    Double_t& Rn0(){return fRmin0;}// Returns address of fRmin0
    Double_t& Rx0(){return fRmax0;}// Returns address of fRmax0
    Double_t& Rn1(){return fRmin1;}// Returns address of fRmin1
    Double_t& Rx1(){return fRmax1;}// Returns address of fRmax1
    void Print(ostream *os); // Prints output content of this class
    void Read(istream *is); // Reads output created by Print above.
 private:
    Double_t fPhi0,fPhi1;   // Starting and ending phi angles [degrees]
    Double_t fRmin0,fRmin1; // Inner Radius
    Double_t fRmax0,fRmax1; // Outer Radius
    Double_t fDz;   // Z half length.

    ClassDef(AliITSConeSegData,1) // Cone segment data class
};
// Input and output function for standard C++ input/output.
ostream &operator<<(ostream &os,AliITSConeSegData &source);
istream &operator>>(istream &os,AliITSConeSegData &source);
#endif

#ifndef ALIITSPCONEDATA_H
#define ALIITSPCONEDATA_H
/////////////////////////////////////////////////////////////////////////
//  Geant 3 Poly-Cone data structure.
/////////////////////////////////////////////////////////////////////////
#include <math.h> // for the definision of NAN.

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
	{/*if(i<0||i>=fNz) return NAN;else*/ return fZ[i];}
    Double_t *Z(){return fZ;} // Returns array of z values
    Double_t *Rmin() {return fRmin;} // Returm arrau of rmin values
    Double_t *Rmax() {return fRmax;} // Returm arrau of rmax values
    Double_t Rmin(Int_t i) // Return Inner radius value
	{/*if(i<0||i>=fNz) return NAN;else*/ return fRmin[i];}
    Double_t Rmax(Int_t i) // Return Outer radius value
	{/*if(i<0||i>=fNz) return NAN;else*/ return fRmax[i];}
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

#ifndef ALIITSSPHEREDATA_H
#define ALIITSSPHEREDATA_H
/////////////////////////////////////////////////////////////////////////
//  Geant 3 Tube Segment data structure.
/////////////////////////////////////////////////////////////////////////

class AliITSSphereData : public AliITSBaseVolParams {
 public:
    AliITSSphereData() : AliITSBaseVolParams() // Default constructor
	{fTheta[0]=0.0,fTheta[1]=0.0;fPhi[0]=0.0,fPhi[1]=0.0;fRmin=0;fRmax=0;}
    virtual ~AliITSSphereData(){} // Standard destructor
    // Getters
    Double_t Theta0() {return fTheta[0];} // Returns starting Phi value
    Double_t Theta1() {return fTheta[1];} // Returns endinging Phi value
    Double_t Phi0() {return fPhi[0];} // Returns starting Phi value
    Double_t Phi1() {return fPhi[1];} // Returns endinging Phi value
    Double_t Rmin() {return fRmin;} // Returm arrau of rmin values
    Double_t Rmax() {return fRmax;} // Returm arrau of rmax values
    Double_t& T0(){return fTheta[0];}// Returns addres of fTheta[0]
    Double_t& T1(){return fTheta[1];}// Returns addres of fTheta[1]
    Double_t& P0(){return fPhi[0];}// Returns addres of fPhi[0]
    Double_t& P1(){return fPhi[1];}// Returns addres of fPhi[1]
    Double_t& Rn(){return fRmin;}// Returns address of fRmin
    Double_t& Rx(){return fRmax;}// Returns address of fRmax
    void Print(ostream *os); // Prints output content of this class
    void Read(istream *is); // Reads output created by Print above.
 private:
    Double_t fTheta[2]; // Starting and ending theta angles [degree]
    Double_t fPhi[2];   // Starting and ending phi angles [degrees]
    Double_t fRmin; // Inner Radius
    Double_t fRmax; // Outer Radius

    ClassDef(AliITSSphereData,1) // sphere data class
};
// Input and output function for standard C++ input/output.
ostream &operator<<(ostream &os,AliITSSphereData &source);
istream &operator>>(istream &os,AliITSSphereData &source);
#endif

#ifndef ALIITSPARALLELEPIPEDDATA_H
#define ALIITSPARALLELEPIPEDDATA_H
/////////////////////////////////////////////////////////////////////////
//  Geant 3 Box data structure.
/////////////////////////////////////////////////////////////////////////

class AliITSParallelpipedData : public AliITSBaseVolParams {
 public:
    AliITSParallelpipedData() : AliITSBaseVolParams() // Default constructor
	{fDx=0.0;fDy=0.0;fDz=0;}
    virtual ~AliITSParallelpipedData(){} // Standard destructor
    // Getters
    Double_t DxAt() {return fDx;} // Returm arrau of rmin values
    Double_t DyAt() {return fDy;} // Returm arrau of rmax values
    Double_t DzAt() {return fDz;} // Return fDz coordiante
    Double_t Theta() {return fTheta;} // Returm arrau of Theta values
    Double_t Phi() {return fPhi;} // Return Phi coordiante
    Double_t Alpha(){return fAlpha;} // Return fAlph
    Double_t& Dx(){return fDx;}// Returns address of fRmax
    Double_t& Dy(){return fDy;}// Returns address of fRmin
    Double_t& Dz(){return fDz;}// Returns address of fDz
    Double_t& Th(){return fTheta;}// Returns address of fTheta
    Double_t& Ph(){return fPhi;}// Returns address of fPhi
    Double_t& A(){return fAlpha;}// Returns address of fAlpha
    void Print(ostream *os); // Prints output content of this class
    void Read(istream *is); // Reads output created by Print above.
 private:
    Double_t fDx;    // X half length
    Double_t fDy;    // Y half length
    Double_t fDz;    // Z half length.
    Double_t fAlpha; // angle formed by the y axis and by the plane 
                     // joining the center of teh faces parallel to the 
                     // z-x plane at -dY and +dy [degree].
    Double_t fTheta; //polar angle of the line joining the centers of 
                     // the faces at -dz and +dz in z [degree].
    Double_t fPhi;   // azimuthal angle of teh line joing the centers 
                     // of the faaces at -dz and +dz in z [degree].

    ClassDef(AliITSParallelpipedData,1) // Parallel piped data class
};
// Input and output function for standard C++ input/output.
ostream &operator<<(ostream &os,AliITSParallelpipedData &source);
istream &operator>>(istream &os,AliITSParallelpipedData &source);
#endif   

#ifndef ALIITSPGONDATA_H
#define ALIITSPGONDATA_H
/////////////////////////////////////////////////////////////////////////
//  Geant 3 Poly-Gon data structure.
/////////////////////////////////////////////////////////////////////////
#include <math.h> // for the definision of NAN.

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
	{/*if(i<0||i>=fNz) return NAN;else*/ return fZ[i];}
    Double_t Rmin(Int_t i) // Return Inner radius value
	{/*if(i<0||i>=fNz) return NAN;else*/ return fRmin[i];}
    Double_t Rmax(Int_t i) // Return Outer radius value
	{/*if(i<0||i>=fNz) return NAN;else*/ return fRmax[i];}
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
/////////////////////////////////////////////////////////////////////////
//  A basic geometry class for the ITS simulation geometry stucture
/////////////////////////////////////////////////////////////////////////
#include <TObject.h>

class AliModule;
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
    static Int_t ITSG3VnameToIndex(const char name[3]);
    static Int_t ITSIndexToITSG3name(const Int_t i); // Get Geant3 volume name
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
    virtual void Box(AliITSBoxData &d,Int_t med);
    // Greate A Trapizoid with the x dimension varing along z.
    virtual void Trapezoid1(const char *gnam,const TString &dis,Double_t dxn,
                    Double_t dxp,Double_t dy,Double_t dz,Int_t med);
    virtual void Trapezoid1(AliITSTrapezoid1Data &d,Int_t med);
    // Greate A Trapizoid with the x and y dimension varing along z.
    virtual void Trapezoid2(const char *gnam,const TString &dis,Double_t dxn,
                    Double_t dxp,Double_t dyn,Double_t dyp,Double_t dz,
                    Int_t med);
    virtual void Trapezoid2(AliITSTrapezoid2Data &d,Int_t med);
    // General trapazoid.
    virtual void Trapezoid(const char *gnam,const TString &dis,Double_t dz,
                   Double_t thet,Double_t phi,Double_t h1,Double_t bl1,
                   Double_t tl1,Double_t alp1,Double_t h2,Double_t bl2,
                   Double_t tl2,Double_t alp2,Int_t med);
    virtual void Trapezoid(AliITSTrapezoidData &d,Int_t med);
    // Twisted genral trapezoid.
    virtual void TwistedTrapezoid(const char *gnam,const TString &dis,
				  Double_t dz,
				  Double_t thet,Double_t phi,Double_t twist,
				  Double_t h1,Double_t bl1,Double_t tl1,
				  Double_t apl1,Double_t h2,Double_t bl2,
				  Double_t tl2,Double_t apl2,Int_t med);
    virtual void TwistedTrapezoid(AliITSTrapezoidTwistedData &d,Int_t med);
    // Simple Tube.
    virtual void Tube(const char *gnam,const TString &dis,Double_t rmin,
              Double_t rmax,Double_t dz,Int_t med);
    virtual void Tube(AliITSTubeData &d,Int_t med);
    // Tube segment.
    virtual void TubeSegment(const char *gnam,const TString &dis,
			     Double_t rmin,Double_t rmax,Double_t dz,
			     Double_t phi1,Double_t phi2,Int_t med);
    virtual void TubeSegment(AliITSTubeSegData &v,Int_t med);
    // Cut tube.
    virtual void CutTube(const char *gnam,const TString &dis,Double_t rmin,
                 Double_t rmax,Double_t dz,Double_t phi1,Double_t phi2,
                 Double_t lx,Double_t ly,Double_t lz,Double_t hx,Double_t hy,
                 Double_t hz,Int_t med);
    virtual void CutTube(AliITSTubeCutData &d,Int_t med);
    // Ellliptical cross-sectino tube
    virtual void TubeElliptical(const char *gnam,const TString &dis,
				Double_t p1,Double_t p2,Double_t dz,
				Int_t med);
    virtual void TubeElliptical(AliITSTubeEllipticalData &v,Int_t med);
    // Hyperbolic tube
    virtual void HyperbolicTube(const char *gnam,const TString &dis,
				Double_t rmin,Double_t rmax,Double_t dz,
				Double_t thet,Int_t med);
    virtual void HyperbolicTube(AliITSTubeHyperbolicData &d,Int_t med);
    // Simple Cone.
    virtual void Cone(const char *gnam,const TString &dis,Double_t dz,
	      Double_t rmin1,Double_t rmax1,Double_t rmin2,Double_t rmax2,
	      Int_t med);
    virtual void Cone(AliITSConeData &d,Int_t med);
    // Segment of a Cone.
    virtual void ConeSegment(const char *gnam,const TString &dis,Double_t dz,
                     Double_t rmin1,Double_t rmax1,Double_t rmin2,
                     Double_t rmax2,Double_t phi1,Double_t phi2,Int_t med);
    virtual void ConeSegment(AliITSConeSegData &d,Int_t med);
    //Poly-Cone
    virtual void PolyCone(const char *gnam,const TString &dis,Double_t phi1,
                  Double_t dphi,Int_t nz,Double_t *z,Double_t *rmin,
                  Double_t *rmax,Int_t med);
    virtual void PolyCone(AliITSPConeData &d,Int_t med);
    // Spherical shell segment.
    virtual void Sphere(const char *gnam,const TString &dis,Double_t rmin,
                Double_t rmax,Double_t the1,Double_t the2,Double_t phi1,
                Double_t phi2,Int_t med);
    virtual void Sphere(AliITSSphereData &d,Int_t med);
    // Parallelepiped.
    virtual void Parallelepiped(const char *gnam,const TString &dis,
				Double_t dx,Double_t dy,Double_t dz,
				Double_t alph,Double_t thet,
				Double_t phi,Int_t med);
    virtual void Parallelepiped(AliITSParallelpipedData &d,Int_t med);
    // Polygon.
    virtual void PolyGon(const char *gnam,const TString &dis,Double_t phi1,
                 Double_t dphi,Int_t npdv,Int_t nz,Double_t *z,Double_t *rmin,
                 Double_t *rmax,Int_t med);
    virtual void PolyGon(AliITSPGonData &d,Int_t med);
    // Position one volume inside another
    virtual void Pos(const char *vol,Int_t cn,const char *moth,Double_t x,
             Double_t y,Double_t z,Int_t irot);
    // Position one volume inside another
    virtual void Pos(AliITSBaseVolParams &v,Int_t cn,
	     AliITSBaseVolParams &m,TVector3 &t,Int_t irot);
    // Position one volume inside another
    virtual void Pos(AliITSBaseVolParams &v,AliITSBaseVolParams &m,
		     TVector3 &t,Int_t irot){Pos(v,m.GetG3cpn(),m,t,irot);};
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
    // gives angles in degree between 0.<=t<360.
    Double_t Mod360(Double_t t){if(t>=360.) return Mod360(t-360.);
        if(t<0.0) return Mod360(t+360.);return t;}
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
    static const Double_t fAlpha; //! find structure constant
    static const Double_t fRe; //![cm]classical elect. radius
    static const Double_t fNa; //! [#/mole] Avogadro's number
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
     
