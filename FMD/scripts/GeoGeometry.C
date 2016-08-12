//____________________________________________________________________
//
//
// $Id$
//
// Script I used for rapid prototyping of the FMD3 geometry - in
// particular the support cone 
//
/** @defgroup FMD_geo_geom Simple geometry
    @ingroup FMD_script
*/
#ifndef __CINT___
# ifndef ROOT_TGeoVolume
#  include <TGeoVolume.h>
# endif
# ifndef ROOT_TGeoMaterial
#  include <TGeoMaterial.h>
# endif
# ifndef ROOT_TGeoMedium
#  include <TGeoMedium.h>
# endif
# ifndef ROOT_TGeoMatrix
#  include <TGeoMatrix.h>
# endif
# ifndef ROOT_TGeoXtru
#  include <TGeoXtru.h>
# endif
# ifndef ROOT_TGeoPcon
#  include <TGeoPcon.h>
# endif
# ifndef ROOT_TGeoTube
#  include <TGeoTube.h>
# endif
# ifndef ROOT_TGeoManager
#  include <TGeoManager.h>
# endif
# ifndef ROOT_TGeoPhysicalNode
#  include <TGeoPhysicalNode.h>
# endif
# ifndef ROOT_TVector3
#  include <TVector3.h>
# endif
# ifndef ROOT_TString
#  include <TString.h>
# endif
# ifndef ROOT_TRandom
#  include <TRandom.h>
# endif
# ifndef ROOT_TBrowser
#  include <TBrowser.h>
# endif
# ifndef ROOT_TCanvas
#  include <TCanvas.h>
# endif
# ifndef __CSTDARG__
#  include <cstdarg>
# endif
# ifndef __IOSTREAM__
#  include <iostream>
# endif
#endif

//____________________________________________________________________
/** @brief Simple geometry 
    @ingroup FMD_geo_geom
    @code 
    gSystem->Load("libPhysics");
    gSystem->Load("libGeom");
    gROOT->LoadMacro("GeoGeometry.C+");
    Geometry g;
    g.Exec();
    @endcode 
 */
class Geometry 
{
public:
  /** Constructor */
  Geometry();
  /** Destructor */
  virtual ~Geometry() {}
  /** Initialize  */
  virtual void Initialize();
  /** Register  */
  virtual void Register();
  /** @e Do-It member function  */
  virtual void Exec();
  /** Try to align */
  virtual void Align();
  /** Convert detector coordinates to spatial coordinates 
      @param sector Sector number 
      @param strip  Strip number
      @param xyz    Spatial coordinates
  */
  virtual void Detector2XYZ(UInt_t sector, UInt_t strip, TVector3& xyz);
  /** Debug level */
  static Int_t fgDebug;
  /** Print debugging messages 
      @param lvl Acceptance level 
      @param where Where it happened
      @param format Message format. */
  static void Debug(Int_t lvl, const char* where, const char* format, ...) 
  {
    if (lvl > fgDebug) return;
    static char buf[2048];
    va_list ap;
    va_start(ap, format);
    vsnprintf(buf, 2048, format, ap);
    std::cout << "D: " << where << ": " << buf << std::endl;
    va_end(ap);
  }
protected:
  /** List of translations */
  TObjArray*    fMatricies;
  /** Shape parameter */
  TVector2    	fA;
  /** Shape parameter */
  TVector2    	fB;
  /** Shape parameter */
  TVector2    	fC;
  /** Spacing between sensor and hybrid */
  Double_t      fSpacer;
  /** Thickness of aluminum plates in honey comb */
  Double_t    	fAlThickness;
  /** Width of bonding pads */
  Double_t      fBondingWidth;
  /** chip layer thickenss */
  Double_t      fChipThickness;
  /** Copper layer thickness */
  Double_t      fCopperThickness;
  /** Upper radious */
  Double_t 	fHighR;
  /** Thickness of honey comb */
  Double_t	fHoneycombThickness;
  /** Inner honey comb inner radius  */
  Double_t	fInnerHoneyHighR;
  /** Inner honey comb outer radius  */
  Double_t	fInnerHoneyLowR;
  /** Z coordinate of inner ring  */
  Double_t	fInnerZ;
  /** Length of support legs */
  Double_t	fLegLength;
  /** Offset from edge of legs */
  Double_t	fLegOffset;
  /** Radius of support legs */
  Double_t	fLegRadius;
  /** Inner radius */
  Double_t	fLowR;
  /** Spacing between modules */
  Double_t	fModuleSpacing;
  /** Number of strps */
  Int_t		fNStrips;
  /** Outer honey comb inner radius */
  Double_t	fOuterHoneyHighR;
  /** Outer honey comb outer radius */
  Double_t	fOuterHoneyLowR;
  /** Z coordinate of outer */
  Double_t	fOuterZ;
  /** Thickness of print board */
  Double_t	fPrintboardThickness;
  /** Cache of ring depth */
  Double_t	fRingDepth;
  /** Thickness of silicon sensor */
  Double_t	fSiThickness;
  /** ??? */
  Double_t	fSpacerHeight;
  /** Opening angle of sensor */
  Double_t	fTheta;
  /** Radius of wafer the sensor is cut out of */
  Double_t	fWaferRadius;
  /** Air tracking medium */
  TGeoMedium* 	fAir;
  /** Aluminum tracking medium */
  TGeoMedium* 	fAl;
  /** Carbon tracking medium */
  TGeoMedium*   fCarbon;
  /** Chip tracking medium */
  TGeoMedium*   fChip;
  /** Copper tracking medium */
  TGeoMedium*   fCopper;
  /** Kapton tracking medium */
  TGeoMedium*   fKapton;
  /** PCB tracking medium */
  TGeoMedium*	fPCB;
  /** Plastic tracking medium */
  TGeoMedium*	fPlastic;
  /** Active silicon tracking medium */
  TGeoMedium*	fSi;
  /** Vacuum  tracking medium */
  TGeoMedium*	fVacuum;
  /** Use assemblies */
  Bool_t        fDoubleAssembly;
  /** Make a detailed geometry */
  void MakeDetailed(TGeoVolume* mother);
};

//____________________________________________________________________
Int_t Geometry::fgDebug = 10;

//____________________________________________________________________
Geometry::Geometry()
  : fMatricies(0),
    fA( 4.3000,   1.3972),
    fB(17.2000,   2.1452),
    fC(15.3327,   4.9819)
{
  fBondingWidth           =   0.5000;
  fWaferRadius            =   6.7000;
  fSiThickness            =   0.0300;
  fLowR                   =   4.3000;
  fHighR                  =  17.2000;
  fTheta                  =  18.0000;
  fNStrips                =      512;
  fRingDepth              =   2.1600;
  fLegRadius              =   0.5000;
  fLegLength              =   1.0000;
  fLegOffset              =   2.0000;
  fModuleSpacing          =   1.0000;
  fPrintboardThickness    =   0.1000;
  fSpacerHeight           =   0.0300;

  fInnerZ                 =  83.4000;
  fOuterZ                 =  75.2000;
  fHoneycombThickness     =   1.0000;
  fAlThickness            =   0.1000;
  fInnerHoneyLowR         =   5.3000;
  fInnerHoneyHighR        =  18.2000;
  fOuterHoneyLowR         =  16.6000;
  fOuterHoneyHighR        =  29.0000;

  fCopperThickness        =   0.01;
  fChipThickness          =   0.01;
  fAlThickness            =   0.10;
  fHoneycombThickness     =   1.00;
  fSpacer                 =   0.10;
  fLegLength              =   1.00;
}

//____________________________________________________________________
void 
Geometry::Exec()
{
  Initialize();
  Register();
  // gGeoManager->DefaultColors();
  gGeoManager->ViewLeaves(kFALSE);
  TCanvas* c = new TCanvas;
  c->SetFillColor(0);
  gGeoManager->GetTopVolume()->Draw();
  new TBrowser("b", "Browser");
}

//____________________________________________________________________
void 
Geometry::Initialize()
{
  Double_t mMax  = 10;
  Double_t mType = 2;
  
  // Air 
  Double_t pAir[]     = { 0., mType, mMax, 1.,  .001, 1., .001, .001 };
  TGeoMixture* fmdAir = new TGeoMixture("FMD air", 4, .00120479);
  fmdAir->DefineElement(0, 12.0107,  6., 0.000124);
  fmdAir->DefineElement(1, 14.0067,  7., 0.755267);
  fmdAir->DefineElement(2, 15.9994,  8., 0.231781);
  fmdAir->DefineElement(3, 39.948,  18., 0.012827);
  fmdAir->SetTransparency('0');
  fAir = new TGeoMedium("FMD Air", 1, fmdAir,  pAir);

  // Silicon 
  Double_t pSi[]      = { 1., mType, mMax, 1.,  .001, 1., .001, .001 };
  TGeoMaterial* fmdSi = new TGeoMaterial("FMD Si", 28.0855, 14, 2.33);
  fmdSi->SetFillColor(2);
  fSi = new TGeoMedium("FMD Si", 1, fmdSi,  pSi);

  // Vacumm 
  Double_t pVacuum[]  = { 0., mType, mMax, 10,  .01, .1, .003, .003 };
  TGeoMaterial* fmdVacuum = new TGeoMaterial("FMD Vacuum",1e-16,1e-16,1e-16);
  fmdVacuum->SetTransparency('0');
  fVacuum = new TGeoMedium("FMD Vacuum", 1, fmdVacuum,pVacuum);


  // PCB 
  Double_t pPCB[]     = { 0., mType, mMax, 1.,  .001, 1., .001, .001 };
  TGeoMixture* fmdPCB = new TGeoMixture("FMD PCB", 14, 1.8);
  fmdPCB->DefineElement(0,  28.0855,   14, 0.15144894);
  fmdPCB->DefineElement(1,  40.078,    20, 0.08147477);
  fmdPCB->DefineElement(2,  26.981538, 13, 0.04128158);
  fmdPCB->DefineElement(3,  24.305,    12, 0.00904554);
  fmdPCB->DefineElement(4,  10.811,     5, 0.01397570);
  fmdPCB->DefineElement(5,  47.867,    22, 0.00287685);
  fmdPCB->DefineElement(6,  22.98977,  11, 0.00445114);
  fmdPCB->DefineElement(7,  39.0983,   19, 0.00498089);
  fmdPCB->DefineElement(8,  55.845,    26, 0.00209828);
  fmdPCB->DefineElement(9,  18.9984,    9, 0.00420000);
  fmdPCB->DefineElement(10, 15.9994,    8, 0.36043788);
  fmdPCB->DefineElement(11, 12.0107,    6, 0.27529425);
  fmdPCB->DefineElement(12, 14.0067,    7, 0.01415852);
  fmdPCB->DefineElement(13,  1.00794,   1, 0.03427566);
  fmdPCB->SetFillColor(7);
  fPCB = new TGeoMedium("FMD PCB", 1, fmdPCB,  pPCB);

  // Chip 
  Double_t pChip[]  = { 0., mType, mMax, 10., .01,  1., .003, .003 };
  TGeoMixture* fmdChip = new TGeoMixture("FMD Chip", 6, 2.36436);
  fmdChip->DefineElement(0,  12.0107,   6, 0.039730642);
  fmdChip->DefineElement(1,  14.0067,   7, 0.001396798);
  fmdChip->DefineElement(2,  15.9994,   8, 0.01169634);
  fmdChip->DefineElement(3,   1.00794,  1, 0.004367771);
  fmdChip->DefineElement(4,  28.0855,  14, 0.844665);
  fmdChip->DefineElement(5, 107.8682,  47, 0.0981434490);
  fmdChip->SetFillColor(4);
  fChip = new TGeoMedium("FMD Chip", 1, fmdChip, pChip);

  // Carbon 
  Double_t pC[]       = { 0., mType, mMax, 10., .01,  1., .003, .003 };
  TGeoMaterial* fmdC = new TGeoMaterial("FMD C", 12.011, 6, 2.265);
  fmdC->SetFillColor(6);
  fCarbon = new TGeoMedium("FMD C", 1, fmdC, pC);
  
  // Kapton (inside of Honeycombs)
  Double_t pKapton[]  = { 0., mType, mMax, 1.,  .001, 1., .001, .001 };
  TGeoMaterial* fmdKapton = new TGeoMaterial("FMD Kapton", 12.011,   6.,0.01);
  fmdKapton->SetFillColor(3);
  fKapton                 = new TGeoMedium("FMD Kapton",1,fmdKapton,pKapton);

  // Plastic 
  Double_t pPlastic[] = { 0., mType, mMax, 10., .01,  1., .003, .003 };
  TGeoMixture* fmdPlastic = new TGeoMixture("FMD Plastic", 2, 1.03);
  fmdPlastic->DefineElement(0,  1.01,   1, .5);
  fmdPlastic->DefineElement(1,  12.01,  6, .5);
  fmdPlastic->SetFillColor(4);
  fPlastic = new TGeoMedium("FMD Plastic", 1, fmdPlastic,  pPlastic);

  // Aluminium 
  Double_t pAl[]  = { 0., mType, mMax, 10., .001,  -1., .003, .003 };
  TGeoMaterial* fmdAl = new TGeoMaterial("FMD Al", 26.981539, 13, 2.7);
  fmdAl->SetFillColor(3);
  fAl = new TGeoMedium("FMD Al", 1, fmdAl, pAl);

  // Copper 
  Double_t pCopper[]  = { 0., mType, mMax, 10., .01,  1., .003, .003 };
  TGeoMaterial* fmdCopper = new TGeoMaterial("FMD Copper", 63.546,  29.,8.96);
  fmdCopper->SetFillColor(3);
  fCopper = new TGeoMedium("FMD Copper",  1, fmdCopper, pCopper);
}

#define DEGRAD TMath::Pi()/180.
//____________________________________________________________________
void 
Geometry::Detector2XYZ(UInt_t sector, UInt_t strip, TVector3& xyz)
{
  UInt_t   mod      = sector / 2;
  if (!fMatricies) {
    Warning("Detector2XYZ", "No matricies");
    return;
  }
  TGeoMatrix* m = static_cast<TGeoMatrix*>(fMatricies->At(mod));
  if (!m) {
    Warning("Detector2XYZ", "No matrix found for module %d", mod);
    return;
  }
  Debug(10, "Detector2XYZ", "Transforming (%2d,%3d)", sector, strip);
  Double_t rmax     = fB.Mod();
  Double_t stripoff = fA.Mod();
  Double_t dstrip   = (rmax - stripoff) / fNStrips;
  Double_t r        = (strip + .5) * dstrip + stripoff; // fLowR
  Double_t theta    = ((sector % 2) - .5) * fTheta;
  Double_t modThick = (fSiThickness 
		       + fPrintboardThickness 
		       + fCopperThickness
		       + fChipThickness 
		       + fSpacer);
  Debug(10,"Detector2XYZ", "Radius %7.3f, angle %7.3f (%f, %f)", r, theta, 
       fLowR, stripoff);
  Double_t local[] = {
    r * TMath::Cos(theta * DEGRAD), 
    r * TMath::Sin(theta * DEGRAD), 
    -modThick + fSiThickness / 2
  };
  Double_t master[3];
  Debug(10, "Detector2XYZ", "Local (%7.3f,%7.3f,%7.3f)", 
       local[0], local[1], local[2]);
  m->LocalToMaster(local, master);
  Debug(10, "Detector2XYZ", "Master (%7.3f,%7.3f,%7.3f)", 
       master[0], master[1], master[2]);
  xyz.SetXYZ(master[0], master[1], master[2]);
}
//____________________________________________________________________
void
Geometry::MakeDetailed(TGeoVolume* caveVolume)
{
  Info("MakeSimple", "Using a detailed geometry");
  Double_t xv[6] = { fA.X(), fC.X(), fB.X(),  fB.X(),  fC.X(),  fA.X() };
  Double_t yv[6] = { fA.Y(), fC.Y(), fB.Y(), -fB.Y(), -fC.Y(), -fA.Y() };
  Double_t rmax     = fB.Mod();
  Double_t stripoff = fA.Mod();
  Double_t dstrip   = (rmax - stripoff) / fNStrips;

  // Double_t hybridThick = (fPrintboardThickness + fCopperThickness 
  //			     + fChipThickness);
  // Double_t modThick    = fSiThickness + fSpacer + hybridThick;
  // Double_t fmdWidth    = (modThick + fHoneycombThickness + fModuleSpacing
  // 			     + fLegLength + fSpacer);
  
  // Top
  // TGeoTube*   fmdShape     = new TGeoTube(fLowR-.1,fHighR+.1,fmdWidth/2+.1);
  // TGeoPcon*   fmdShape     = new TGeoPcon(0, 360, 2);
  // fmdShape->DefineSection(0, -fmdWidth / 2 - .1, fLowR-.1, fHighR+.1);
  // fmdShape->DefineSection(1,  fmdWidth / 2 + .1, fLowR-.1, fHighR+.1);
  // TGeoVolume* fmdVolume    = new TGeoVolume("FMD1", fmdShape, fVacuum);
  
  
  // Sensor 
  TGeoXtru* sensorShape = new TGeoXtru(2);
  sensorShape->DefinePolygon(6, xv, yv);
  sensorShape->DefineSection(0, - fSiThickness/2);
  sensorShape->DefineSection(1, fSiThickness/2);
  TGeoVolume* sensorVolume= new TGeoVolume("FISE",sensorShape,fSi);
  sensorVolume->SetLineColor(2);
  // sensorVolume->VisibleDaughters(kFALSE);
  // sensorVolume->SetVisibility(kTRUE);

  // Virtual volume shape to divide 
  TGeoTubeSeg* activeShape   = new TGeoTubeSeg(fLowR, rmax, fSiThickness/2, 
					       - fTheta, fTheta);
  TGeoVolume*  activeVolume  = new TGeoVolume("FIAC", activeShape,fSi);
  activeVolume->SetLineColor(3);
  sensorVolume->AddNodeOverlap(activeVolume, 0);
  // activeVolume->SetTransparency(0x3f);
  TGeoVolume* sectorVolume   = activeVolume->Divide("FISC",2,2,-fTheta,
						    0,0,"N");
  TGeoVolume* stripVolume    = sectorVolume->Divide("FIST",1,fNStrips,
						    stripoff,dstrip,0,"SX");
  (void)stripVolume;
  
  // Position
  Double_t x, y, z;
  
  // Make PCB volume 
  for (Int_t i = 0; i < 3; i++) yv[i] -= fBondingWidth;
  for (Int_t i = 3; i < 6; i++) yv[i] += fBondingWidth;
  Double_t off = (TMath::Tan(TMath::Pi() * fTheta / 180) * fBondingWidth);

  // PCB layer 
  TGeoXtru* pcbShape      = new TGeoXtru(2);
  pcbShape->DefinePolygon(6, xv, yv);
  pcbShape->DefineSection(0, - fPrintboardThickness/2);
  pcbShape->DefineSection(1, fPrintboardThickness/2);
  TGeoVolume* pcbVolume   = new TGeoVolume("FPCB",pcbShape,fPCB);
  pcbVolume->SetLineColor(4);

  // Copper layer
  TGeoXtru* cuShape       = new TGeoXtru(2);
  cuShape->DefinePolygon(6, xv, yv);
  cuShape->DefineSection(0, - fCopperThickness/2);
  cuShape->DefineSection(1, fCopperThickness/2);
  TGeoVolume* cuVolume    = new TGeoVolume("FCOP",cuShape,fCopper);
  cuVolume->SetLineColor(4);

  // Chip layer
  TGeoXtru*   chipShape   = new TGeoXtru(2);
  chipShape->DefinePolygon(6, xv, yv);
  chipShape->DefineSection(0, - fChipThickness/2);
  chipShape->DefineSection(1, fChipThickness/2);
  TGeoVolume* chipVolume = new TGeoVolume("FCHI",chipShape,fChip);
  chipVolume->SetLineColor(4);

  // Short leg shape 
  TGeoTube*   shortLegShape  = new TGeoTube(0, fLegRadius, fLegLength / 2);
  TGeoVolume* shortLegVolume = new TGeoVolume("FISL", shortLegShape, fPlastic);
  shortLegVolume->SetLineColor(2);
  // Long leg shape
  TGeoTube*   longLegShape   = new TGeoTube(0, fLegRadius, 
					    (fLegLength + fModuleSpacing) / 2);
  TGeoVolume* longLegVolume  = new TGeoVolume("FILL", longLegShape, fPlastic);
  longLegVolume->SetLineColor(2);
  
  // Make a front volume 
  TGeoVolume* frontVolume    =  new TGeoVolumeAssembly("FIFV");
  z                          =  fPrintboardThickness / 2;
  frontVolume->AddNode(pcbVolume, 0, new TGeoTranslation(0, 0, z));
  z                          += (fPrintboardThickness + fCopperThickness) / 2;
  frontVolume->AddNode(cuVolume, 0, new TGeoTranslation(0, 0, z));
  z                          += (fCopperThickness + fChipThickness) / 2;
  frontVolume->AddNode(chipVolume, 0, new TGeoTranslation(0, 0, z));
  z                          += (fChipThickness + fModuleSpacing 
				 + fLegLength) / 2;
  x                          =  fA.X() + fLegOffset + fLegRadius;
  y                          =  0;
  frontVolume->AddNode(longLegVolume, 0, new TGeoTranslation(x, y, z));
  x                          =  fC.X();
  y                          =  fC.Y() - fLegOffset - fLegRadius - off;
  frontVolume->AddNode(longLegVolume, 1, new TGeoTranslation(x,y,z));
  y                          =  -y;
  frontVolume->AddNode(longLegVolume, 2, new TGeoTranslation(x,y,z));

  // Make a back volume 
  TGeoVolume* backVolume     =  new TGeoVolumeAssembly("FIBV");
  z                          =  fPrintboardThickness / 2;
  backVolume->AddNode(pcbVolume, 0, new TGeoTranslation(0, 0, z));
  z                          += (fPrintboardThickness + fCopperThickness) / 2;
  backVolume->AddNode(cuVolume, 0, new TGeoTranslation(0, 0, z));
  z                          += (fCopperThickness + fChipThickness) / 2;
  backVolume->AddNode(chipVolume, 0, new TGeoTranslation(0, 0, z));
  z                          += (fChipThickness + fLegLength) / 2;
  x                          =  fA.X() + fLegOffset + fLegRadius;
  y                          =  0;
  backVolume->AddNode(shortLegVolume, 0, new TGeoTranslation(x, y, z));
  x                          =  fC.X();
  y                          =  fC.Y() - fLegOffset - fLegRadius - off;
  backVolume->AddNode(shortLegVolume, 1, new TGeoTranslation(x,y,z));
  y                          =  -y;
  backVolume->AddNode(shortLegVolume, 2, new TGeoTranslation(x,y,z));

#if 1
  TGeoVolume* ringTopVolume  = new TGeoVolumeAssembly("FITV");
  TGeoVolume* ringBotVolume  = new TGeoVolumeAssembly("FIBV");
  Double_t dz = 0;
#else
  Double_t dz = (fSiThickness + fSpacer + fModuleSpacing + 
		 fPrintboardThickness + fCopperThickness + 
		 fChipThickness + fLegLength);
  TGeoTubeSeg* topRingShape = new TGeoTubeSeg(fA.X(), rmax, dz/2, 0, 180);
  TGeoTubeSeg* botRingShape = new TGeoTubeSeg(fA.X(), rmax, dz/2, 180, 360);
  TGeoVolume* ringTopVolume = new TGeoVolume("FITV", topRingShape, fAir);
  TGeoVolume* ringBotVolume = new TGeoVolume("FIBV", botRingShape, fAir);
#endif

  TGeoVolume* half = ringTopVolume;
  // Place in mother 
  for (Int_t i = 0; i < 10; i++) {
    if (i > 4)  half  = ringBotVolume;
    Bool_t      front = ((i % 2) == 0);
    z                 =  -dz/2 + fSiThickness / 2 + (i % 2) * fModuleSpacing;
    TGeoMatrix* m1    = new TGeoCombiTrans(0,0,z,0);
    m1->RotateZ((2 * i + 1) * fTheta);
    half->AddNode(sensorVolume, i, m1);
    TGeoVolume* vol   = (front ? frontVolume : backVolume);
    z                 += fSpacer + fSiThickness / 2;
    TGeoMatrix* m2    = new TGeoCombiTrans(0,0,z,0);
    m2->RotateZ((2 * i + 1) * fTheta);
    half->AddNode(vol, i, m2);    
  }

  TGeoVolume* fmdTopVolume  = new TGeoVolumeAssembly("FMT1");
  TGeoVolume* fmdBotVolume  = new TGeoVolumeAssembly("FMB1");
  fmdTopVolume->AddNode(ringTopVolume, 0, new TGeoTranslation(0,0,dz/2));
  fmdBotVolume->AddNode(ringBotVolume, 0, new TGeoTranslation(0,0,dz/2));
  
  // Top of Honeycomb
  TGeoTubeSeg* hcShape       = new TGeoTubeSeg(fLowR, fHighR,
					       fHoneycombThickness / 2,
					       0, 180);
  TGeoVolume*  hcVolume      = new TGeoVolume("F1II", hcShape, fAl);
  hcVolume->VisibleDaughters(kFALSE);
  hcVolume->SetLineColor(5);
  // Air in top of honeycomb
  TGeoTubeSeg* ihcShape      = new TGeoTubeSeg(fLowR + fAlThickness, 
					       fHighR - fAlThickness, 
					       (fHoneycombThickness 
						- fAlThickness) / 2,
					       0, 180);
  TGeoVolume*  ihcVolume     = new TGeoVolume("F1IK", ihcShape, fKapton);
  hcVolume->AddNode(ihcVolume, 0);

  // Add honey comb to mothers. 
  z = (fSiThickness + fSpacer + fPrintboardThickness + fCopperThickness
       + fChipThickness + fModuleSpacing + fLegLength + fHoneycombThickness/2);
  fmdTopVolume->AddNode(hcVolume, 0, new TGeoTranslation(0,0,z));
  TGeoMatrix* r = new TGeoCombiTrans(0,0,z, 0); r->RotateZ(180);
  fmdBotVolume->AddNode(hcVolume, 1, r);

  z = 0;
  caveVolume->AddNode(fmdTopVolume, 1, new TGeoTranslation(0,0,z));
  caveVolume->AddNode(fmdBotVolume, 1, new TGeoTranslation(0,0,z));
}

//____________________________________________________________________
void
Geometry::Register()
{
  // Framework::Task::Register(option);
  TGeoShape* caveShape = new TGeoBBox(fHighR+10, fHighR+10, fHighR+10);
  TGeoVolume* caveVolume = new TGeoVolume("Cave", caveShape, fVacuum);
  gGeoManager->SetTopVolume(caveVolume);

  if (!gGeoManager->GetVolume("FM1T")) MakeDetailed(caveVolume);
  gGeoManager->CloseGeometry();
}

//____________________________________________________________________
void
Geometry::Align()
{
  if (!gGeoManager)                 return;
  if (!gGeoManager->GetTopVolume()) return;
  if (!gGeoManager->IsClosed())     return;
  
  if (!fMatricies) fMatricies = new TObjArray;
  TGeoIterator next(gGeoManager->GetTopVolume());
  TGeoNode* node = 0;

  while ((node = static_cast<TGeoNode*>(next()))) {
    // node->Print();
    if (node->GetName()[0] == 'F' && node->GetName()[2] == 'S' && 
	node->GetName()[3] == 'E') {
      // We go an FMD module 
      TString path("/");
      path.Append(gGeoManager->GetNode(0)->GetName());
      Int_t nlevel = next.GetLevel();
      Debug(10, "Exec", "Level is %d", nlevel);
      for (int lvl = 0; lvl <= nlevel; lvl++) {
	TGeoNode* p = next.GetNode(lvl);
	if (!p) continue;
	Debug(10, "Exec", "Adding '%s' to path '%s'",
	      p->GetName(), path.Data());
	if (!path.IsNull()) path.Append("/");
	path.Append(p->GetName());
      }
      TGeoPhysicalNode* pnode = gGeoManager->MakePhysicalNode(path.Data());
      TGeoHMatrix* matrix = new TGeoHMatrix(*node->GetMatrix());
      TGeoRotation pertub;
      Double_t angles[] = { gRandom->Uniform(0, 3), 
			    gRandom->Uniform(0, 3), 
			    gRandom->Uniform(0, 3) };
      pertub.RotateX(angles[0]);
      pertub.RotateY(angles[1]);
      pertub.RotateZ(angles[1]);
      *matrix *= pertub;
      Debug(5, "Exec", "Aliging %s (%f,%f,%f)",
	    pnode->GetName(), angles[0], angles[1], angles[2]);
      pnode->Align(matrix);
    }
  }
}

//____________________________________________________________________
//
// EOF
//

  
