#ifndef ALIITSV11_H
#define ALIITSV11_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
  $Id$
 */

/////////////////////////////////////////////////////////////////////////
//  Manager and hits classes for set: ITS version 11, 2003 geometry    //
/////////////////////////////////////////////////////////////////////////
 
#include "AliITS.h"
 
class AliITSv11 : public AliITS {

 public:
    AliITSv11();
    AliITSv11(const char *title);
    AliITSv11(const AliITSv11 &source); // copy constructor
    AliITSv11& operator=(const AliITSv11 &source); // assignment operator
    virtual       ~AliITSv11() ;
    virtual void   BuildGeometry();
    virtual void   CreateGeometry();
    virtual void   CreateMaterials();
    virtual Int_t  IsVersion() const {// returns the ITS version number 
                                      return 11;} 
    virtual void   Init(); 
    virtual void   SetDefaults();
    virtual void   DrawModule();
    virtual void   StepManager();
    virtual void   SetWriteDet(Bool_t det=kTRUE){ // set .det write
	                                         fGeomDetOut = det;}
    virtual void   SetWriteDet(const char *f){ // set write file
	                             strncpy(fWrite,f,60);fGeomDetOut = kTRUE;}
    virtual void   SetReadDet(Bool_t det=kTRUE){ //set .det read
	                                        fGeomDetIn = det;}
    virtual void   SetReadDet(const char *f){ // set read file
	                               strncpy(fRead,f,60);fGeomDetIn = kTRUE;}
    virtual void   SetEUCLID(Bool_t euclid=kTRUE){ // set write Euclid file
	                                          fEuclidOut = euclid;}
    virtual void   SetEUCLIDFileName(const char *f){ // set write file
	                     fEuclidGeometry=f;fEuclidOut = kTRUE;}
    virtual void   SetMinorVersion(Int_t v=22){ // Choose between existing minor versions
	fMinorVersion = v;}
    virtual void   SetThicknessDet1(Float_t v=200.){ 
	 // Set detector thickness in layer 1
	 fDet1 = v;}
    virtual void   SetThicknessDet2(Float_t v=200.){ 
	 // Set detector thickness in layer 2
	 fDet2 = v;}
    virtual void   SetThicknessChip1(Float_t v=300.){ 
	 // Set chip thickness in layer 1
	 fChip1 = v;}	 	 
    virtual void   SetThicknessChip2(Float_t v=200.){ 
	 // Set chip thickness in layer 2
	 fChip2 = v;}
    virtual void   SetRails(Int_t v=1){ 
	 // Set flag for rails
	 fRails = v;}	 
    virtual void   SetCoolingFluid(Int_t v=1){ 
	 // Set flag for cooling fluid
	 fFluid = v;}	 	 
    virtual Bool_t GetEUCLID(){return fEuclidOut;}// returns value Euclid flag.
    virtual const char  *GetEULIIDFileName() const{ // return .euc file name
	                               return fEuclidGeometry.Data();}
    virtual Bool_t GetWriteDet() { // returns value GeomDetOut flag.
	                          return fGeomDetOut;}
    virtual Bool_t GetReadDet() { // returns value GeomDetIn flag.
	                         return fGeomDetIn;}
    virtual char  *GetReadDetFileName(){ // return .det read file name
	          if(fRead[0]!='\0') return fRead; else return fEuclidGeomDet;}
    virtual char  *GetWriteDetFileName(){ // return .det write file name
	        if(fWrite[0]!='\0') return fWrite; else return fEuclidGeomDet;}
    virtual Int_t GetMajorVersion(){// return Major Version Number
	return fMajorVersion;}
    virtual Int_t GetMinorVersion(){// return Major Version Number
	return fMinorVersion;}
    virtual Float_t GetThicknessDet1(){ 
	 // Get detector thickness in layer 1
	 return fDet1;}
    virtual Float_t GetThicknessDet2(){ 
	 // Get detector thickness in layer 2
	 return fDet2;}
    virtual Float_t GetThicknessChip1(){ 
	 // Get chip thickness in layer 1
	 return fChip1;}	 	 
    virtual Float_t GetThicknessChip2(){ 
	 // Get chip thickness in layer 2
	 return fChip2;}
    virtual Int_t GetRails(){ 
	 // Get flag for rails
	 return fRails;}	 
    virtual Int_t GetCoolingFluid(){ 
	 // Get flag for cooling fluid
	 return fFluid;}	 	 	 
	 	 
 private:
    void InitAliITSgeom();
    void SetScalecm(){// Sets scale factor for centemeters
    fScale = 1.0;}
    void SetScalemm(){// Sets scale factor for milimeters
    fScale = 0.10;}
    void SetScalemicrons(){// Sets scale factor for micronsmeters
    fScale = 1.0E-04;}
    void SetScale(Double_t s=1.0){// Sets scale factor
    fScale = s;}
    Double_t GetScale(){// Returns the scale factor
    return fScale;}
    Bool_t IsScalecm(){// Returens kTRUE if scale factor is set of [cm]
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
    void TubeElliptical(const char gnam[3],const TString &dis,Double_t p1,
			Double_t p2,Double_t dz,Int_t med);
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
	fidmed = fIdtmed->GetArray()-199;}
    // Define rotation matrix
    void Matrix(Int_t irot,Double_t thet1,Double_t phi1,Double_t thet2,
		Double_t phi2,Double_t thet3,Double_t phi3);
    // Defube ritatuib matrix
    void Matrix(Int_t irot,Double_t rot[3][3]);
    // Rotation matrix about axis i (i=0=>x, i=1=>y, i=2=>z).
    void Matrix(Int_t irot,Int_t axis,Double_t thet);
    // Rotation matrix about x axis
    void XMatrix(Int_t irot,Double_t thet){
	Matrix(irot,0,thet);}
    // Rotation matrix about y axis
    void YMatrix(Int_t irot,Double_t thet){
	Matrix(irot,1,thet);}
    // Rotation matrix about z axis
    void ZMatrix(Int_t irot,Double_t thet){
	Matrix(irot,2,thet);}
    // Define Element material and medium
    void Element(Int_t imat,const char *name,Int_t z,Double_t dens,Int_t istd);
    // Returns standard radiation lenghts of elements.
    Float_t GetRadLength(Int_t z);
    // Returns natrual abundance atomic mass numbers for a given element
    Float_t GetA(Int_t z);
    // Returns ITS standard Theata Max transport cut values
    Float_t GetStandardThetaMax(Int_t istd);
    // Returns ITS standard Theata Max transport cut values
    Float_t GetStandardMaxStepSize(Int_t istd);
    // Returns ITS standard Theata Max transport cut values
    Float_t GetStandardEfraction(Int_t istd);
    // Returns ITS standard Theata Max transport cut values
    Float_t GetStandardEpsilon(Int_t istd);
    // Degree Versions of TMath functions (as needed)
    Double_t Sind(Double_t t){return TMath::Sin(TMath::Pi()*t/180.);}
    Double_t Cosd(Double_t t){return TMath::Cos(TMath::Pi()*t/180.);}
    Double_t Tand(Double_t t){return TMath::Tan(TMath::Pi()*t/180.);}
    Double_t ASind(Double_t t){return 180.0*TMath::ASin(t)/TMath::Pi();}
    Double_t ACosd(Double_t t){return 180.0*TMath::ACos(t)/TMath::Pi();}
    Double_t ATand(Double_t t){return 180.0*TMath::ATan(t)/TMath::Pi();}
    Double_t ATand2(Double_t y,Double_t x){return 180.0*TMath::ATan2(y,x)/TMath::Pi();}

    // TString fEuclidGeomtery,fEuclidMaterial defined in AliModule.
    Bool_t fEuclidOut;        // Flag to write geometry in euclid format
    Bool_t fGeomDetOut;       // Flag to write .det file out
    Bool_t fGeomDetIn;        // Flag to read .det file or directly from Geat.
    Int_t  fMajorVersion;     // Major version number == IsVersion
    Int_t  fMinorVersion;     // Minor version number
    char   fEuclidGeomDet[60];// file where detector transormation are define.
    char   fRead[60];         //! file name to read .det file
    char   fWrite[60];        //! file name to write .det file
    Float_t  fDet1;	      // thickness of detector in SPD layer 1
    Float_t  fDet2;	      // thickness of detector in SPD layer 2
    Float_t  fChip1;	      // thickness of chip in SPD layer 1   
    Float_t  fChip2;	      // thickness of chip in SPD layer 2   
    Int_t    fRails;          // flag to switch rails on (=1) and off (=0)
    Int_t    fFluid;          // flag to switch between water (=1) and freon (=0)
    Int_t fIDMother;          //! ITS Mother Volume id.
    //
    Int_t *fidmed;            //! array of media indexes.
    Int_t *fidrot;            //! array of rotation matrixies indexes.
    Double_t fScale;          //! scale factor (=1=>[cm])

    ClassDef(AliITSv11,1)  //Hits manager for set:ITS version 11 
                                 // PPR detailed Geometry asymmetric
};
 
#endif
