#ifndef TRD_H
#define TRD_H
////////////////////////////////////////////////
//  Manager and hits classes for set:TRD     //
////////////////////////////////////////////////
 
#include "AliDetector.h"
#include "AliHit.h" 
 
class AliTRD : public AliDetector {
 
public:
  AliTRD();
  AliTRD(const char *name, const char *title);
  virtual      ~AliTRD() {}
  virtual void  AddHit(Int_t, Int_t*, Float_t*);
  virtual void  BuildGeometry();
  virtual void  CreateGeometry() {}
  virtual void  CreateMaterials();
  Int_t         DistancetoPrimitive(Int_t px, Int_t py);
  virtual void  Init();
  virtual Int_t IsVersion() const =0;
  virtual void  StepManager()=0; 
  virtual void  DrawModule() {}
  
  ClassDef(AliTRD,1)       // Transition Radiation Detector base class
};

const Int_t   nsect   = 18;      //Number of sectors in the full detector
const Int_t   nmodul  = 6;       //Number of modules in each sector
const Float_t rmin    = 281;     //r-Coordinates of the TRD-frame
const Float_t rmax    = 350.282;
const Float_t zmax1   = 351.2;   //z-Coordinates of the TRD-frame
const Float_t zmax2   = 292.35;
const Float_t alframe = 1.0;     //Thickness of the aluminium of the support frame
const Float_t alfram1 = 1.0;
const Float_t alfram2 = 0.5;
const Float_t inframe = 3.0;     //Thickness of the interior of the support frame
const Float_t ccframe = 1.0;     //Thickness of the carbon chamber frame
const Float_t pethick = 0.15;    //Thickness of the PE-layer in the radiator
const Float_t pezpos  =  0.;     //z-position of the PE-layer in the radiator
const Float_t rathick = 6.23;    //Thickness of the radiator
const Float_t razpos  = -2.6585; //z-position of the radiator
const Float_t mythick = 0.005;   //Thickness of the mylar-layer
const Float_t myzpos  =  0.459;  //z-position of the mylar-layer
const Float_t xethick = 3.6;     //Thickness of the Xe/C02-layer
const Float_t xezpos  =  2.2615; //z-position of the Xe/C02-layer
const Float_t cuthick = 0.002;   //Thickness of the Cu-layer (Pads)
const Float_t cuzpos  =  4.0625; //z-position of the Cu-layer (Pads)
const Float_t kathick = 0.01;    //Thickness of the kapton-layer
const Float_t kazpos  =  4.0695; //z-position of the kapton-layer
const Float_t nothick = 0.05;    //Thickness of the NOMEX-layer
const Float_t nozpos  =  4.8235; //z-position of the NOMEX-layer
const Float_t rothick = 0.018;   //Thickness of the readout-layer
const Float_t rozpos  =  5.2;    //z-position of the readout-layer
 
 
//_____________________________________________________________________________
class AliTRDhit : public AliHit {

public:
  Int_t        fSector;     // TRD sector number
  Int_t        fChamber;    // TRD chamber number
  Int_t        fPlane;      // TRD plane number 
  Float_t      fQ;          // Charge created by a hit (geometry 2)
 
public:
  AliTRDhit() {}
  AliTRDhit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits);
  virtual ~AliTRDhit() {}
 
  ClassDef(AliTRDhit,1)     // Hits for Transition Radiation Detector
};

#endif
