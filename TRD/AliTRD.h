#ifndef TRD_H
#define TRD_H
////////////////////////////////////////////////
//  Manager and hits classes for set: TRD     //
////////////////////////////////////////////////
 
#include "AliDetector.h"
#include "AliHit.h" 
 
class AliTRD : public AliDetector {
 
protected:

  Int_t         fGasMix;         // Gas mixture. 0: Xe/Isobutane 1: Xe/CO2
  Int_t         fSensSelect;     // Switch to select only parts of the detector
  Int_t         fSensPlane;      // Sensitive detector plane
  Int_t         fSensChamber;    // Sensitive detector chamber
  Int_t         fSensSector;     // Sensitive detector sector

public:
  AliTRD();
  AliTRD(const char *name, const char *title);
  virtual      ~AliTRD() {}
  virtual void  AddHit(Int_t, Int_t*, Float_t*);
  virtual void  BuildGeometry();
  virtual void  CreateMaterials();
  Int_t         DistancetoPrimitive(Int_t px, Int_t py);
  virtual void  Init();
  virtual Int_t IsVersion() const =0;
  virtual void  StepManager()=0; 
  virtual void  SetGasMix(Int_t imix = 0);
  virtual void  SetSensPlane(Int_t iplane = 0);
  virtual void  SetSensChamber(Int_t ichamber = 0);
  virtual void  SetSensSector(Int_t isector = 0);

  ClassDef(AliTRD,1)             // Transition Radiation Detector base class
};

////////////////////////////////////////////////
//  Geometry parameter
////////////////////////////////////////////////

const Int_t   nsect   = 18;      // Number of sectors in the full detector
const Int_t   nmodul  = 6;       // Number of modules in each sector
const Int_t   ncham   = 6;       // Number of different chambers
const Int_t   narmsec = 5;       // Number of sectors in one arm (geometry 1)

const Float_t rmin    = 294.0;   // r-dimensions of the TRD-frame
const Float_t rmax    = 368.0;

const Float_t zmax1   = 378.35;  // z-dimensions of the TRD-frame
const Float_t zmax2   = 302.0;

const Float_t zleni   = 110.0;   // z-dimension of the inner chamber
const Float_t zlenn   = 156.0;   // z-dimension of the neighbouring chambers
const Float_t zleno   = 156.0;   // z-dimension of the outer chambers

const Float_t widpl1  = 99.6;    // rphi-dimension of the innermost plane

const Float_t alframe = 1.0;     // Thickness of the aluminium of the support frame
 
const Float_t ccframe = 1.0;     // Thickness of the carbon chamber frame

// Thicknesses of the the material layers
const Float_t sethick = 0.02;    // Radiator seal
const Float_t rathick = 4.2;     // Radiator
const Float_t pethick = 0.20;    // PE-layer in the radiator
const Float_t mythick = 0.005;   // Mylar-layer
const Float_t xethick = 3.5;     // Gas mixture
const Float_t cuthick = 0.001;   // Pad plane
const Float_t suthick = 0.06;    // HEXCEL+G10 support structure (= 0.31% X0)
const Float_t fethick = 0.0044;  // FEE + signal lines (= 0.31% X0)
const Float_t cothick = 0.02;    // PE of the cooling device
const Float_t wathick = 0.01;    // Cooling water

// Position of the material layers
const Float_t sezpos  = -5.657;  // Radiator seal
const Float_t razpos  = -3.557;  // Radiator
const Float_t pezpos  =  0.0;    // PE-layer in the radiator
const Float_t myzpos  = -1.4545; // Mylar-layer
const Float_t xezpos  =  0.298;  // Gas mixture
const Float_t cuzpos  =  2.047;  // Pad plane
const Float_t suzpos  =  3.046;  // HEXCEL+G10 support structure
const Float_t fezpos  =  4.0482; // FEE + signal lines
const Float_t cozpos  =  4.1504; // PE of the cooling devices
const Float_t wazpos  =  4.3004; // Cooling water

////////////////////////////////////////////////
// Parameter for the energy loss calculation
////////////////////////////////////////////////

const Float_t kPoti = 12.1;      // First ionization potential (eV)
const Float_t kEend = 50000.0;   // Maximum energy (50 keV);

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
