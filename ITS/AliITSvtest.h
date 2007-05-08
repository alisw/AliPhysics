#ifndef ALIITSVTEST_H
#define ALIITSVTEST_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////
//     Manager and hits classes for  ITS version 5
////////////////////////////////////////////////////////////////////////

#include "AliITSInitGeometry.h"
#include "AliITS.h"

class TBRIK;
class AliITSvtest : public AliITS {
 public:
    AliITSvtest(); // Default constructor
    AliITSvtest(const Char_t *title,Int_t version=0); // Standard Constructor
    // Constructor for Euclid geometries
    AliITSvtest(const char *fileeuc,const char *filetme,
                const char *name, const char *title);
    virtual ~AliITSvtest(); // Destructor
    void CreateGeometry(); // Create test geometry 
    void CreateMaterials(); // Create test geometry materials
    void Init(); // Initlizie test geometry for transport
    Int_t IsVersion() const {// returns the ITS version number 
        return kvtest;}
    void StepManager(); // Transport step manager
    void SetWriteDet(Bool_t det=kTRUE){ // set .det write
        fGeomDetOut = det;}
    void SetWriteDet(const char *f){ // set write file
        fWrite=f;fGeomDetOut = kTRUE;}
    void SetReadDet(Bool_t det=kTRUE){ //set .det read
        fGeomDetIn = det;}
    void SetReadDet(const char *f){ // set read file
        fRead=f;fGeomDetIn = kTRUE;}
    void SetEUCLIDFileName(const char *f){ // set write file
        fEuclidGeometry=f;fEuclidOut = kTRUE;}
    void SetMinorVersion(Int_t v){ // Choose between existing minor versions
	fMinorVersion = 1;if(v==1) fMinorVersion = 1;
	else if(v==2) fMinorVersion = 2;
	else Warning("SetMinorVersion","Undefined Minor Version setting =1");}
    Bool_t GetWriteDet() const { // returns value GeomDetOut flag.
        return fGeomDetOut;}
    Bool_t GetReadDet() const { // returns value GeomDetIn flag.
        return fGeomDetIn;}
    Char_t *GetEULIIDFileName() const{ // return .euc file name
        return (Char_t*)(fEuclidGeometry.Data());}
    Char_t *GetReadDetFileName()const{ // return .det read file name
        if(!fRead.IsNull()) return (Char_t*)(fRead.Data());
        else return GetEULIIDFileName();}
    Char_t *GetWriteDetFileName()const{ // return .det write file name
        if(!fWrite.IsNull()) return (Char_t*)(fWrite.Data());
        else return GetEULIIDFileName();}
    Int_t GetMajorVersion() const {// return Major Version Number
	return fMajorVersion;}
    Int_t GetMinorVersion() const {// return Major Version Number
	return fMinorVersion;}
    // Print class in ascii form to stream
    void PrintAscii(ostream *os)const;
    // Read in class in ascii form from stream
    void ReadAscii(istream *is);

  protected:
    // Units, Convert from k?? to cm,degree,GeV,seconds,
    static const Double_t fgkmicron; // Convert micron to TGeom's cm.
    static const Double_t fgkmm; // Convert mm to TGeom's cm.
    static const Double_t fgkcm; // Convert cm to TGeom's cm.
    static const Double_t fgkDegree; //Convert degrees to TGeom's degrees
    static const Double_t fgkRadian; //To Radians
    static const Double_t fgkgcm3;   // Density in g/cm^3
    static const Double_t fgkCelsius; // Temperature in degrees Celcius
    static const Double_t fgkPascal;  // Preasure in Pascal
    static const Double_t fgkKPascal;  // Preasure in KPascal
    static const Double_t fgkeV;  // Energy in eV
    static const Double_t fgkKeV;  // Energy in KeV
    static const Double_t fgkMeV;  // Energy in MeV
    static const Double_t fgkGeV;  // Energy in GeV


 private:
    AliITSvtest(const AliITSvtest &source); // copy constructor
    AliITSvtest& operator=(const AliITSvtest &source); // assignment operator
    void CreateMaterialsEuclid(); // Create test geometry materials from Euclid
    void CreateGeometryEuclid(); // Create test geometry Geometry from Euclid
    void InitEuclid(); // Create test geometry Init for Euclid
    void CreateMaterials2(); // Create test geometry materials from geometry2
    void CreateGeometry2(); // Create test geometry Geometry from geometry2
    void Init2(); // Create test geometry Init for geometry2
    // TString fEuclidGeomtery,fEuclidMaterial defined in AliModule.
    Bool_t  fGeomDetOut;       // Flag to write .det file out
    Bool_t  fGeomDetIn;        // Flag to read .det file or directly from Geat.
    Int_t   fMajorVersion;     // Major version number == IsVersion
    Int_t   fMinorVersion;     // Minor version number
    TString fEuclidGeomDet;// file where detector transormation are define.
    TString fRead;         //! file name to read .det file
    TString fWrite;        //! file name to write .det file
    AliITSInitGeometry fIgm;//! Geometry initilization object

    ClassDef(AliITSvtest,2)  //Hits manager for ITS test version, Private ITS class for different test geometries
};
// Input and output function for standard C++ input/output.
ostream &operator<<(ostream &os,const AliITSvtest &s);
istream &operator>>(istream &is,AliITSvtest &s);
 
#endif
