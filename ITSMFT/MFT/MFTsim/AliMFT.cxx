// **************************************************************************
// * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
// *                                                                        *
// * Author: The ALICE Off-line Project.                                    *
// * Contributors are mentioned in the code where appropriate.              *
// *                                                                        *
// * Permission to use, copy, modify and distribute this software and its   *
// * documentation strictly for non-commercial purposes is hereby granted   *
// * without fee, provided that the above copyright notice appears in all   *
// * copies and that both the copyright notice and this permission notice   *
// * appear in the supporting documentation. The authors make no claims     *
// * about the suitability of this software for any purpose. It is          *
// * provided "as is" without express or implied warranty.                  *
// **************************************************************************

//====================================================================================================================================================
//
//      Main class of the ALICE Muon Forward Tracker
//
//      Contact author: antonio.uras@cern.ch
//
//====================================================================================================================================================

#include "AliLog.h"
#include "AliCodeTimer.h"
#include "TFile.h"
#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoMatrix.h"
#include "TVirtualMC.h"
#include "TClonesArray.h"
#include "TGeoGlobalMagField.h"
#include "AliRun.h"
#include "AliLoader.h"
#include "AliDetector.h"
#include "AliMC.h"
#include "AliMagF.h"
#include "AliMFT.h"
#include "AliMFTHit.h"
#include "AliMFTDigit.h"
#include "AliMFTCluster.h"
#include "AliTrackReference.h"
#include "AliMFTSegmentation.h"
#include "AliMFTDigitizer.h"
#include "AliMFTPlane.h"
#include "AliMFTConstants.h"
#include "TString.h"
#include "TObjArray.h"
#include "AliMFTGeometry.h"

ClassImp(AliMFT)

//====================================================================================================================================================

AliMFT::AliMFT():
AliDetector(),
fVersion(1),
fNPlanes(0),
fSDigitsPerPlane(0),
fDigitsPerPlane(0),
fRecPointsPerPlane(0),
fSideDigits(0),
fSegmentation(0),
fNameGeomFile(0),
fChargeDispersion(25.e-4),
fSingleStepForChargeDispersion(0),
fNStepForChargeDispersion(4),
fDensitySupportOverSi(0.036),
fFileNameForUnderyingEvent(0),
fFileNameForPileUpEvents(0),
fNPileUpEvents(0),
fUnderlyingEventID(-1),
fGeomTGeo(0)
{
  
  // default constructor
  
  for (Int_t iPileUp=0; iPileUp<AliMFTConstants::fNMaxPileUpEvents; iPileUp++) fPileUpEventsIDs[iPileUp] = -1;
  
  
}

//====================================================================================================================================================

AliMFT::AliMFT(const Char_t *name, const Char_t *title):
AliDetector(name, title),
fVersion(1),
fNPlanes(0),
fSDigitsPerPlane(0),
fDigitsPerPlane(0),
fRecPointsPerPlane(0),
fSideDigits(0),
fSegmentation(0),
fNameGeomFile(0),
fChargeDispersion(25.e-4),
fSingleStepForChargeDispersion(0),
fNStepForChargeDispersion(4),
fDensitySupportOverSi(0.036),
fFileNameForUnderyingEvent(0),
fFileNameForPileUpEvents(0),
fNPileUpEvents(0),
fUnderlyingEventID(-1),
fGeomTGeo(0)
{
  
  for (Int_t iPileUp=0; iPileUp<AliMFTConstants::fNMaxPileUpEvents; iPileUp++) fPileUpEventsIDs[iPileUp] = -1;
  
  fNameGeomFile = "AliMFTGeometry.root";
  
  SetGeometry();
  
}

//====================================================================================================================================================

AliMFT::AliMFT(const Char_t *name, const Char_t *title, Char_t *nameGeomFile):
AliDetector(name, title),
fVersion(1),
fNPlanes(0),
fSDigitsPerPlane(0),
fDigitsPerPlane(0),
fRecPointsPerPlane(0),
fSideDigits(0),
fSegmentation(0),
fNameGeomFile(0),
fChargeDispersion(25.e-4),
fSingleStepForChargeDispersion(0),
fNStepForChargeDispersion(4),
fDensitySupportOverSi(0.036),
fFileNameForUnderyingEvent(0),
fFileNameForPileUpEvents(0),
fNPileUpEvents(0),
fUnderlyingEventID(-1)
{
  
  for (Int_t iPileUp=0; iPileUp<AliMFTConstants::fNMaxPileUpEvents; iPileUp++) fPileUpEventsIDs[iPileUp] = -1;
  
  fNameGeomFile = nameGeomFile;
  
  SetGeometry();
  
}

//====================================================================================================================================================

AliMFT::~AliMFT() {
  
  if (fSDigitsPerPlane)   { fSDigitsPerPlane->Delete();    delete fSDigitsPerPlane;   }
  if (fDigitsPerPlane)    { fDigitsPerPlane->Delete();     delete fDigitsPerPlane;    }
  if (fRecPointsPerPlane) { fRecPointsPerPlane->Delete();  delete fRecPointsPerPlane; }

  delete fGeomTGeo;  

}

//====================================================================================================================================================

void AliMFT::Init() {

  fGeomTGeo = new AliMFTGeomTGeo();

}
  
//====================================================================================================================================================

void AliMFT::CreateMaterials() {
  
  /// \todo Check all materials Remove the one unneeded
  
  
  // Definition of MFT materials  - to be updated to the most recent values
  
  AliDebug(1,"Start MFT materials");
  
    
  // data from PDG booklet 2002                 density [gr/cm^3]     rad len [cm]           abs len [cm]
  Float_t   aSi = 28.085 ,    zSi   = 14. ,     dSi      =  2.329 ,   radSi   =  21.82/dSi , absSi   = 108.4/dSi  ;    // Silicon
  Float_t   aCarb = 12.01 ,   zCarb =  6. ,     dCarb    =  2.265 ,   radCarb =  18.8 ,      absCarb = 49.9       ;    // Carbon
  Float_t   aAlu = 26.98 ,    zAlu  = 13. ,     dAlu     =  2.70  ,   radAlu  =  8.897 ,     absAlu  = 39.70      ;    // Aluminum
  Float_t   aBe = 9.012182 ,  zBe   = 4. ,      dBe      =  1.85 ,    radBe   =  65.19/dBe , absBe   = 77.8/dBe  ;     // Beryllium
  Float_t   aCu = 63.546 ,    zCu  = 29.  ,     dCu      =  8.96  ,   radCu   =  1.436 ,     absCu   = 15.32      ;    // Copper

  // Air mixture
  const Int_t nAir = 4;
  Float_t   aAir[nAir] = {12, 14, 16, 36} ,  zAir[nAir] = {6, 7, 8, 18} ,   wAir[nAir]={0.000124, 0.755267, 0.231781, 0.012827} , dAir=0.00120479, dAirVacuum=0.00120479e-4;

  // Water mixture
  const Int_t nWater = 2;
  Float_t   aWater[nWater] = {1.00794, 15.9994} ,  zWater[nWater] = {1, 8} ,   wWater[nWater] = {0.111894, 0.888106} , dWater=1.;
  
  // SiO2 mixture
  const Int_t nSiO2 = 2;
  Float_t   aSiO2[nSiO2] = {15.9994, 28.0855} ,   zSiO2[nSiO2] = {8., 14.} ,   wSiO2[nSiO2] = {0.532565, 0.467435} , dSiO2 = 2.20;
  
  // Inox mixture
  const Int_t nInox = 9;
  Float_t   aInox[nInox] = {12.0107, 54.9380, 28.0855, 30.9738, 32.0660, 58.6928, 51.9961, 95.9400, 55.8450} ;
  Float_t   zInox[nInox] = { 6,      25,      14,      15,      16,      28,      24,      42,      26     } ;
  Float_t   wInox[nInox] = {0.0003,  0.02,    0.01,    0.00045, 0.0003,  0.12,    0.17,    0.025,   0.65395} ;
  Float_t   dInox = 8.03;
  
  // Kapton polyimide film (from SPD AliITSv11.cxx)  and http://physics.nist.gov/cgi-bin/Star/compos.pl?matno=179
  Float_t aKapton[4]={1.00794,12.0107, 14.010,15.9994};
  Float_t zKapton[4]={1.,6.,7.,8.};
  Float_t wKapton[4]={0.026362,0.69113,0.07327,0.209235};
  Float_t dKapton   = 1.42;
  
  //--- EPOXY  --- C18 H19 O3 from ITS AliITSv11.cxx
  Float_t aEpoxy[3] = {15.9994, 1.00794, 12.0107} ;
  Float_t zEpoxy[3] = {     8.,      1.,      6.} ;
  Float_t wEpoxy[3] = {     3.,     19.,     18.} ;
  Float_t dEpoxy = 1.23; //  1.8 very high value from ITS! ou 1.23 from eccobond 45 lv datasheet

  //--- Silicone SE4445 Dow Corning  
  // Si, Al, O, C, H
  Float_t aSE4445[5] = {28.0855, 26.981538, 15.9994, 12.0107, 1.00794} ;
  Float_t zSE4445[5] = {    14.,       13.,      8.,      6.,      1.} ;
  Float_t wSE4445[5] = {  5.531,    45.222,  43.351,   4.717,   1.172} ;
  Float_t dSE4445 = 2.36; //from LBNL file, priv. comm.
  
  //--- CARBON FIBER CM55J --- from ITS AliITSv11.cxx
  Float_t aCM55J[4]={12.0107,14.0067,15.9994,1.00794};
  Float_t zCM55J[4]={6.,7.,8.,1.};
  Float_t wCM55J[4]={0.908508078,0.010387573,0.055957585,0.025146765};
  Float_t dCM55J = 1.33; // new value for MFT, from J.M. Buhour infos

  // Rohacell mixture
  const Int_t nRohacell = 3;
  Float_t aRohacell[nRohacell] = {1.00794, 12.0107, 15.9994};
  Float_t zRohacell[nRohacell] = {1., 6., 8.};
  Float_t wRohacell[nRohacell] = {0.0858, 0.5964, 0.3178};
  Float_t dRohacell = 0.032;  // 0.032 g/cm3 rohacell 31, 0.075 g/cm3 rohacell 71;
  
  // Polyimide pipe mixture
  const Int_t nPolyimide = 4;
  Float_t aPolyimide[nPolyimide] = {1.00794, 12.0107, 14.0067, 15.9994};
  Float_t zPolyimide[nPolyimide] = {1, 6, 7, 8};
  Float_t wPolyimide[nPolyimide] = {0.00942, 0.56089, 0.13082, 0.29887};
  Float_t dPolyimide = 1.4;   

  // PEEK mixture (Polyether Ether Ketone)
  const Int_t nPEEK = 3;
  Float_t   aPEEK[nPEEK] = {1.00794, 12.0107, 15.9994} ;
  Float_t   zPEEK[nPEEK] = {1,       6,        8} ;
  Float_t   wPEEK[nPEEK] = {0.06713, 0.40001,  0.53285} ;
  Float_t   dPEEK = 1.32;
  
  // (Printed Circuit Board), material type FR4
  const Int_t nFR4 = 5;
  Float_t   aFR4[nFR4] = {1.00794,    12.0107, 15.9994, 28.0855,   79.904} ;
  Float_t   zFR4[nFR4] = {1,          6,       8,       14,   35} ;
  Float_t   wFR4[nFR4] = {0.0684428,  0.278042,0.405633, 0.180774,    0.0671091} ;
  Float_t   dFR4 = 1.7; //Density FR4= 1.7 Cu=8.96


  //======================== From ITS code ===================================
  //X7R capacitors - updated from F.Tosello's web page - M.S. 18 Oct 10
  // 58.6928 --> innner electrodes (mainly Ni)
  // 63.5460 --> terminaisons (Cu) 
  // 118.710 --> terminaisons (Sn)
  // 137.327 Ba, 47.867 Ti, 15.9994 O  (mainly BaTiO3)
  Float_t aX7R[6]={137.327,47.867,15.9994,58.6928,63.5460,118.710};
  Float_t zX7R[6]={56.,22.,8.,28.,29.,50.};
  Float_t wX7R[6]={0.524732,0.176736,0.179282,0.079750,0.019750,0.019750};
  Float_t dX7R = 6.07914;
  
  //X7R weld, i.e. Sn 60% Pb 40% (from F.Tosello's web page - M.S. 15 Oct 10)
  
  Float_t aX7Rweld[2]={118.71 , 207.20};
  Float_t zX7Rweld[2]={ 50.   ,  82.  };
  Float_t wX7Rweld[2]={  0.60 ,   0.40};
  Float_t dX7Rweld   = 8.52358;
  //==========================================================================
  
  Int_t   matId  = 0;                        // tmp material id number
  Int_t   unsens = 0, sens=1;                // sensitive or unsensitive medium
  Int_t   itgfld = 3;			     // type of field intergration 0 no field -1 user in guswim 1 Runge Kutta 2 helix 3 const field along z
  Float_t maxfld = 5.; 		             // max field value
  
  Float_t tmaxfd = -10.0;                    // max deflection angle due to magnetic field in one step
  Float_t stemax =  0.001;                   // max step allowed [cm]
  Float_t deemax = -0.2;                     // maximum fractional energy loss in one step 0<deemax<=1
  Float_t epsil  =  0.001;                   // tracking precision [cm]
  Float_t stmin  = -0.001;                   // minimum step due to continuous processes [cm] (negative value: choose it automatically)
  
  Float_t tmaxfdSi =  0.1;                   // max deflection angle due to magnetic field in one step
  Float_t stemaxSi =  5.0e-4;                // maximum step allowed [cm]
  Float_t deemaxSi =  0.1;                   // maximum fractional energy loss in one step 0<deemax<=1
  Float_t epsilSi  =  0.5e-4;                // tracking precision [cm]
  Float_t stminSi  = -0.001;                 // minimum step due to continuous processes [cm] (negative value: choose it automatically)
  
  Int_t    fieldType        = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->Integ();     // Field type
  Double_t maxField         = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->Max();     // Field max.
  
  AliMixture(kAir,"Air$", aAir, zAir, dAir, nAir, wAir);
  AliMedium(kAir,    "Air$", kAir, unsens, fieldType, maxField, tmaxfd, stemax, deemax, epsil, stmin);
  //
  //     Vacuum
  AliMixture(kVacuum, "Vacuum$", aAir, zAir, dAirVacuum, nAir, wAir);
  AliMedium(kVacuum,  "Vacuum$", kVacuum, unsens, itgfld, maxfld, tmaxfd, stemax, deemax, epsil, stmin);

  AliMaterial(++matId, "Si$", aSi, zSi, dSi, radSi, absSi);
  AliMedium(kSi,       "Si$", matId, sens, fieldType, maxField, tmaxfdSi, stemaxSi, deemaxSi, epsilSi, stminSi);
  
  AliMaterial(++matId, "Readout$", aSi, zSi, dSi, radSi, absSi);
  AliMedium(kReadout,  "Readout$", matId, unsens, fieldType, maxField, tmaxfdSi, stemaxSi, deemaxSi, epsilSi, stminSi);
  
  AliMaterial(++matId, "Support$", aSi, zSi, dSi*fDensitySupportOverSi, radSi/fDensitySupportOverSi, absSi/fDensitySupportOverSi);
  AliMedium(kSupport,  "Support$", matId, unsens, fieldType, maxField, tmaxfdSi, stemaxSi, deemaxSi, epsilSi, stminSi);
  
  Double_t maxBending       = 0;     // Max Angle
  Double_t maxStepSize      = 0.001; // Max step size
  Double_t maxEnergyLoss    = 1;     // Max Delta E
  Double_t precision        = 0.001; // Precision
  Double_t minStepSize      = 0.001; // Minimum step size
 // Carbon
  aCarb                = 12.011;
  zCarb                = 6.;
  dCarb          = 2.265;
  radCarb  = 18.8;
  absCarb = 999;
  maxBending       = 10;
  maxStepSize      = .01;
  precision        = .003;
  minStepSize      = .003;
  AliMaterial(matId, "Carbon$", aCarb, zCarb, dCarb, radCarb, absCarb);
  AliMedium(kCarbon, "Carbon$", matId,0,fieldType,maxField,maxBending,
            maxStepSize,maxEnergyLoss,precision,minStepSize);

  //  AliMaterial(++matId, "Carbon$", aCarb, zCarb, dCarb, radCarb, absCarb );
//  AliMedium(kCarbon,   "Carbon$", matId, unsens, fieldType,  maxField, tmaxfd, stemax, deemax, epsil, stmin);
  
  AliMaterial(++matId, "Be$", aBe, zBe, dBe, radBe, absBe );
  AliMedium(kBe,   "Be$", matId, unsens, fieldType,  maxField, tmaxfd, stemax, deemax, epsil, stmin);
  
  AliMaterial(++matId, "Alu$", aAlu, zAlu, dAlu, radAlu, absAlu);
  AliMedium(kAlu,      "Alu$", matId, unsens, fieldType,  maxField, tmaxfd, stemax, deemax, epsil, stmin);
  
  
  AliMixture(++matId, "Water$", aWater, zWater, dWater, nWater, wWater);
  AliMedium(kWater,   "Water$", matId, unsens, itgfld, maxfld, tmaxfd, stemax, deemax, epsil, stmin);
  
  AliMixture(++matId, "SiO2$", aSiO2, zSiO2, dSiO2, nSiO2, wSiO2);
  AliMedium(kSiO2,    "SiO2$", matId, unsens, itgfld, maxfld, tmaxfd, stemax, deemax, epsil, stmin);
  
  AliMixture(++matId, "Inox$", aInox, zInox, dInox, nInox, wInox);
  AliMedium(kInox,    "Inox$", matId, unsens, itgfld, maxfld, tmaxfd, stemax, deemax, epsil, stmin);
  
  AliMixture(++matId, "Kapton$", aKapton, zKapton, dKapton, 4, wKapton);
  AliMedium(kKapton,"Kapton$", matId, unsens, itgfld, maxfld, tmaxfd, stemax, deemax, epsil, stmin);
  
  AliMixture(++matId, "Epoxy$", aEpoxy, zEpoxy, dEpoxy, -3, wEpoxy);
  AliMedium(kEpoxy,"Epoxy$", matId, unsens, itgfld, maxfld, tmaxfd, stemax, deemax, epsil, stmin);
  
  AliMixture(++matId, "SE4445$", aSE4445, zSE4445, dSE4445, -5, wSE4445);
  AliMedium(kSE4445,"SE4445$", matId, unsens, itgfld, maxfld, tmaxfd, stemax, deemax, epsil, stmin);
  
  AliMixture(++matId,"CarbonFiber$",aCM55J,zCM55J,dCM55J,4,wCM55J);
  AliMedium(kCarbonEpoxy,"CarbonFiber$", matId, unsens, itgfld, maxfld, tmaxfd, stemax, deemax, epsil, stmin);
  
  AliMixture(++matId,  "Rohacell", aRohacell, zRohacell, dRohacell, nRohacell, wRohacell);
  AliMedium(kRohacell, "Rohacell", matId, unsens, itgfld, maxfld, tmaxfd, stemax, deemax, epsil, stmin);
  
  AliMixture(++matId,  "Polyimide", aPolyimide, zPolyimide, dPolyimide, nPolyimide, wPolyimide);
  AliMedium(kPolyimide, "Polyimide", matId, unsens, itgfld, maxfld, tmaxfd, stemax, deemax, epsil, stmin);
	
  AliMixture(++matId, "PEEK$", aPEEK, zPEEK, dPEEK, nPEEK, wPEEK);
  AliMedium(kPEEK,    "PEEK$", matId, unsens, itgfld, maxfld, tmaxfd, stemax, deemax, epsil, stmin);
  
  AliMixture(++matId, "FR4$", aFR4, zFR4, dFR4, nFR4, wFR4);
  AliMedium(kFR4,    "FR4$", matId, unsens, itgfld, maxfld, tmaxfd, stemax, deemax, epsil, stmin);
  
  AliMaterial(++matId, "Cu$", aCu, zCu, dCu, radCu, absCu);
  AliMedium(kCu,       "Cu$", matId, unsens, itgfld, maxfld, tmaxfd, stemax, deemax, epsil, stmin);
 
  AliMixture(++matId, "X7Rcapacitors$",aX7R,zX7R,dX7R,6,wX7R);
  AliMedium(kX7R,     "X7Rcapacitors$",matId, unsens, itgfld, maxfld, tmaxfd, stemax, deemax, epsil, stmin);

  AliMixture(++matId, "X7Rweld$",aX7Rweld,zX7Rweld,dX7Rweld,2,wX7Rweld);
  AliMedium(kX7Rw,    "X7Rweld$",matId, unsens, itgfld, maxfld, tmaxfd, stemax, deemax, epsil, stmin);

 // Carbon fleece from AliITSSUv2.cxx
  AliMaterial(++matId,"CarbonFleece$",12.0107,6,0.4,radCarb,absCarb);          // 999,999);  why 999???
  AliMedium(kCarbonFleece,  "CarbonFleece$",matId, unsens, itgfld, maxfld, tmaxfd, stemax, deemax, epsil, stmin);


  AliDebug(1,"End MFT materials");
  
}

//====================================================================================================================================================

void AliMFT::CreateGeometry() {
  
  // Creates detailed geometry simulation (currently GEANT volumes tree)
  
  if(!TVirtualMC::GetMC()->IsRootGeometrySupported()) return;
  
  AliMFTGeometry* mftGeom =   AliMFTGeometry::Instance();
  
  mftGeom->Build();
  
  if (fNStepForChargeDispersion) fSingleStepForChargeDispersion = fChargeDispersion/Double_t(fNStepForChargeDispersion);
  
  Init();

}

//====================================================================================================================================================

void AliMFT::AddAlignableVolumes() {
  
  // Create entries for alignable volumes associating the symbolic volume
  // name with the corresponding volume path. Needs to be syncronized with
  // eventual changes in the geometry.
  
  TString sysName = "MFT";
  TString volPath = "/ALIC_1/MFT_0";
  
  if (!gGeoManager->SetAlignableEntry(sysName.Data(),volPath.Data())) {
    AliFatal(Form("Alignable entry %s not created. Volume path %s not valid", sysName.Data(), volPath.Data()));
  }
  
}

//====================================================================================================================================================

void AliMFT::StepManager() {
  
  // If MFT is not active do nothing
  if (!(this->IsActive())) return;
  
  AliMFTGeometry *mftGeo = AliMFTGeometry::Instance();
  AliMFTSegmentation * seg = mftGeo->GetSegmentation();
  if (!seg) AliFatal("No segmentation available");

  TVirtualMC* mc = fMC;
  
  Double_t absQ = TMath::Abs(mc->TrackCharge());
  if (absQ <= 0) return;

  Int_t copy;
  // Check if hit is into a MFT sensor volume
  if(mc->CurrentVolID(copy) != mftGeo->GetSensorVolumeID() ) return;
  // Get The Sensor Unique ID
  int chipId=-1,ladderId=-1,diskId=-1,halfId=-1,level=0;
  mc->CurrentVolOffID(++level,chipId);
  mc->CurrentVolOffID(++level,ladderId);
  mc->CurrentVolOffID(++level,diskId);
  mc->CurrentVolOffID(++level,halfId);
  Int_t detElemID = mftGeo->GetObjectID(AliMFTGeometry::kSensorType,halfId,diskId,ladderId,chipId);
  AliDebug(1,Form("Found hit into half = %d; disk = %d; ladder = %d; chip = %d",halfId,diskId,ladderId,chipId));

  if (mc->IsTrackExiting()) {
    AddTrackReference(gAlice->GetMCApp()->GetCurrentTrackNumber(), AliTrackReference::kMFT);
  }
  
  static TLorentzVector position, momentum;
  static AliMFTHit hit;
  
  Int_t  status = 0;
  
  // Track status
  if (mc->IsTrackInside())      status += 0x1<<0;
  if (mc->IsTrackEntering())    status += 0x1<<1;
  if (mc->IsTrackExiting())     status += 0x1<<2;
  if (mc->IsTrackOut())         status += 0x1<<3;
  if (mc->IsTrackDisappeared()) status += 0x1<<4;
  if (mc->IsTrackStop())        status += 0x1<<5;
  if (mc->IsTrackAlive())       status += 0x1<<6;
  
  // ---------- Fill hit structure
  
  hit.SetDetElemID(detElemID);
  hit.SetPlane(diskId);
  hit.SetTrack(gAlice->GetMCApp()->GetCurrentTrackNumber());
  
  mc->TrackPosition(position);
  mc->TrackMomentum(momentum);
  
  AliDebug(1, Form(" %s Hit #%06d (x=%f, y=%f, z=%f) belongs to track %02d",
                   mc->CurrentVolName(), fNhits, position.X(), position.Y(), position.Z(), gAlice->GetMCApp()->GetCurrentTrackNumber()));
  
  hit.SetPosition(position);
  hit.SetTOF(mc->TrackTime());
  hit.SetMomentum(momentum);
  hit.SetStatus(status);
  hit.SetEloss(mc->Edep());
  //  hit.SetShunt(GetIshunt());
  //   if (mc->IsTrackEntering()) {
  //     hit.SetStartPosition(position);
  //     hit.SetStartTime(mc->TrackTime());
  //     hit.SetStartStatus(status);
  //     return; // don't save entering hit.
  //   }
  
  // Fill hit structure with this new hit.
  new ((*fHits)[fNhits++]) AliMFTHit(hit);
  
  // Save old position... for next hit.
  //   hit.SetStartPosition(position);
  //   hit.SetStartTime(mc->TrackTime());
  //   hit.SetStartStatus(status);
  
  return;
  
}


//====================================================================================================================================================

void AliMFT::Hits2SDigits(){
  
  // Interface method invoked from AliSimulation to create a list of sdigits corresponding to list of hits. Every hit generates one sdigit.
  AliCodeTimerAuto("",0);

  AliDebug(1, "Start Hits2SDigits.");
  
  AliMFTGeometry *mftGeo = AliMFTGeometry::Instance();
  AliMFTSegmentation * seg = mftGeo->GetSegmentation();
  if (!seg) AliFatal("No segmentation available");
  
  if (!fLoader->TreeH()) fLoader->LoadHits();
  
  if (!fLoader->TreeS()) {
    
    for (Int_t iEvt=0;iEvt<fLoader->GetRunLoader()->GetNumberOfEvents(); iEvt++) {
      
      fLoader->GetRunLoader()->GetEvent(iEvt);
      fLoader->MakeTree("S");
      MakeBranch("S");
      SetTreeAddress();
      
      AliDebug(1, Form("Event %03d: fLoader->TreeH()->GetEntries() = %2d", iEvt, Int_t(fLoader->TreeH()->GetEntries())));

      for (Int_t iTrack=0; iTrack<fLoader->TreeH()->GetEntries(); iTrack++) {
        fLoader->TreeH()->GetEntry(iTrack);
        Hits2SDigitsLocal(Hits(), GetSDigitsList(), iTrack);    // convert these hits to a list of sdigits
      }
      
      fLoader->TreeS()->Fill();
      fLoader->WriteSDigits("OVERWRITE");
      ResetSDigits();
      
    }
  }
  
  fLoader->UnloadHits();
  fLoader->UnloadSDigits();
  
  AliDebug(1,"Stop Hits2SDigits.");
  
}

//====================================================================================================================================================

void AliMFT::Hits2SDigitsLocal(TClonesArray *hits, const TObjArray *pSDig, Int_t track) {
  
  //  Add sdigits of these hits to the list
  
  AliDebug(1, "Entering Hits2SDigitsLocal");
  
  TClonesArray *pSDigList[AliMFTConstants::kNDisks];
  for (Int_t iPlane=0; iPlane<AliMFTConstants::kNDisks; iPlane++) pSDigList[iPlane] = NULL;
  for (Int_t iPlane=0; iPlane<AliMFTConstants::kNDisks;    iPlane++) {
    pSDigList[iPlane] = (TClonesArray*) (*pSDig)[iPlane];
    AliDebug(1,Form("Entries of pSDigList %3d; plane: %02d,",pSDigList[iPlane]->GetEntries(),iPlane));
    if (!track && pSDigList[iPlane]->GetEntries()!=0) AliErrorClass("Some of sdigits lists is not empty");
  }
  
  for (Int_t iHit=0; iHit<hits->GetEntries(); iHit++) {

    AliMFTHit *hit = (AliMFTHit*) hits->At(iHit);

//    AliDebug(1,Form("\n--- New hit x,y,z %f %f %f  ",hit->X(), hit->Y(), hit->Z()));

    // Creating "main digit"
    
    AliMFTDigit *mainSDigit = new AliMFTDigit();
    mainSDigit->SetEloss(hit->GetEloss());
    mainSDigit->SetDetElemID(hit->GetDetElemID());
    mainSDigit->SetPlane(hit->GetPlane());
    mainSDigit->AddMCLabel(hit->GetTrack());
    Int_t xPixel = -1;
    Int_t yPixel = -1;
    AliMFTGeometry * mftGeom = AliMFTGeometry::Instance();

    AliDebug(2,Form("Hit at x,y,z = %f %f %f ",hit->X(), hit->Y(), hit->Z()));
    if(mftGeom->Hit2PixelID(hit->X(), hit->Y(), hit->Z(), mainSDigit->GetDetElemID(), xPixel, yPixel)){
      mainSDigit->SetPixID(xPixel, yPixel, 0);
      mainSDigit->SetPixWidth(AliMFTConstants::kXPixelPitch, AliMFTConstants::kYPixelPitch,AliMFTConstants::kSensorThickness);
      Double_t xCenter,  yCenter, zCenter;
      mftGeom->GetPixelCenter( xPixel,  yPixel, mainSDigit->GetDetElemID(),  xCenter,  yCenter, zCenter );
      mainSDigit->SetPixCenter(xCenter,  yCenter, zCenter );
      new ((*fSideDigits)[fSideDigits->GetEntries()]) AliMFTDigit(*mainSDigit);
      AliDebug(2, Form("Created new sdigit (%f, %f, %f) from hit (%f, %f, %f)",
                       mainSDigit->GetPixelCenterX(), mainSDigit->GetPixelCenterY(), mainSDigit->GetPixelCenterZ(), hit->X(), hit->Y(), hit->Z()));
    } else
      AliDebug(1,Form("Hit outside active area : hit x,y,z = %f ; %f ; %f --> Pixel %d ; %d ",hit->X(), hit->Y(), hit->Z(),xPixel,yPixel));
    
    
    

    // creating "side digits" to simulate the effect of charge dispersion
    
    Double_t pi4 = TMath::Pi()/4.;
    for (Int_t iStep=0; iStep<fNStepForChargeDispersion; iStep++) {
      Double_t shift = (iStep+1) * fSingleStepForChargeDispersion;
      for (Int_t iAngle=0; iAngle<8; iAngle++) {
        Double_t shiftX = shift*TMath::Cos(iAngle*pi4);
        Double_t shiftY = shift*TMath::Sin(iAngle*pi4);
        if (mftGeom->Hit2PixelID(hit->X()+shiftX, hit->Y()+shiftY, hit->Z(), mainSDigit->GetDetElemID(), xPixel, yPixel) ){
          Bool_t digitExists = kFALSE;
          for (Int_t iSideDigit=0; iSideDigit<fSideDigits->GetEntries(); iSideDigit++) {
            if (xPixel==((AliMFTDigit*) fSideDigits->At(iSideDigit))->GetPixelX() &&
                yPixel==((AliMFTDigit*) fSideDigits->At(iSideDigit))->GetPixelY() &&
                mainSDigit->GetDetElemID() == ((AliMFTDigit*) fSideDigits->At(iSideDigit))->GetDetElemID()) {
              digitExists = kTRUE;
              break;
            }
          }
          if (!digitExists) {
            AliMFTDigit *sideSDigit = new AliMFTDigit();
            sideSDigit->SetEloss(0.);
            sideSDigit->SetDetElemID(hit->GetDetElemID());
            sideSDigit->SetPlane(hit->GetPlane());
            sideSDigit->AddMCLabel(hit->GetTrack());
            sideSDigit->SetPixID(xPixel, yPixel, 0);
            mainSDigit->SetPixWidth(AliMFTConstants::kXPixelPitch, AliMFTConstants::kYPixelPitch,AliMFTConstants::kSensorThickness);
            Double_t xCenter,  yCenter, zCenter;
            mftGeom->GetPixelCenter( xPixel,  yPixel, mainSDigit->GetDetElemID(),  xCenter,  yCenter, zCenter );
            mainSDigit->SetPixCenter(xCenter,  yCenter, zCenter );
            new ((*fSideDigits)[fSideDigits->GetEntries()]) AliMFTDigit(*sideSDigit);
          }
        }
      }
    }


    // -------- checking which pixels switched on have their diode actually within the charge dispersion radius
    
    for (Int_t iSDigit=0; iSDigit<fSideDigits->GetEntries(); iSDigit++) {
      AliMFTDigit *mySDig = (AliMFTDigit*) (fSideDigits->At(iSDigit));
      Double_t distance = TMath::Sqrt(TMath::Power(mySDig->GetPixelCenterX()-hit->X(),2) + TMath::Power(mySDig->GetPixelCenterY()-hit->Y(),2));
        if (distance<fChargeDispersion) {
          AliDebug(1,Form("Created new side sdigit (%f, %f, %f) from hit (%f, %f, %f)",
                           mySDig->GetPixelCenterX(), mySDig->GetPixelCenterY(), mySDig->GetPixelCenterZ(), hit->X(), hit->Y(), hit->Z()));
          new ((*pSDigList[mySDig->GetPlane()])[pSDigList[mySDig->GetPlane()]->GetEntries()]) AliMFTDigit(*mySDig);
        }
    }
    
    fSideDigits->Delete();
    
  }
  
  AliDebug(1,"Exiting Hits2SDigitsLocal");
  
}

//====================================================================================================================================================

void AliMFT::MakeBranch(Option_t *option) {
  
  // Create Tree branches
  AliDebug(1, Form("Start with option= %s.",option));
  
  const Int_t kBufSize = 4000;
  
  const Char_t *cH = strstr(option,"H");
  const Char_t *cD = strstr(option,"D");
  const Char_t *cS = strstr(option,"S");
  
  if (cH && fLoader->TreeH()) {
    CreateHits();
    MakeBranchInTree(fLoader->TreeH(), "MFT", &fHits, kBufSize, 0);
  }
  
  if (cS && fLoader->TreeS()) {
    CreateSDigits();
    for(Int_t iPlane=0; iPlane<AliMFTConstants::kNDisks; iPlane++) MakeBranchInTree(fLoader->TreeS(),
                                                                    Form("Plane_%02d",iPlane),
                                                                    &((*fSDigitsPerPlane)[iPlane]),
                                                                    kBufSize, 0);
  }
  
  if (cD && fLoader->TreeD()) {
    CreateDigits();
    for(Int_t iPlane=0; iPlane<AliMFTConstants::kNDisks; iPlane++) MakeBranchInTree(fLoader->TreeD(),
                                                                    Form("Plane_%02d",iPlane),
                                                                    &((*fDigitsPerPlane)[iPlane]),
                                                                    kBufSize, 0);
  }
  
  AliDebug(1,"Stop.");
  
}

//====================================================================================================================================================

void AliMFT::SetTreeAddress() {
  
  
  //Set branch address for the Hits and Digits Tree.
  AliDebug(1, "Start.");
  
  AliDebug(1, Form("AliMFT::SetTreeAddress Hits  fLoader->TreeH() = %p\n", fLoader->TreeH()));
  if (fLoader->TreeH() && fLoader->TreeH()->GetBranch("MFT")) {
    CreateHits();
    fLoader->TreeH()->SetBranchAddress("MFT", &fHits);
  }
  
  AliDebug(1, Form("AliMFT::SetTreeAddress SDigits  fLoader->TreeS() = %p\n", fLoader->TreeS()));
  if (fLoader->TreeS() && fLoader->TreeS()->GetBranch("Plane_00")) {
    CreateSDigits();
    for(Int_t iPlane=0; iPlane<AliMFTConstants::kNDisks; iPlane++) {
      fLoader->TreeS()->SetBranchAddress(Form("Plane_%02d",iPlane), &((*fSDigitsPerPlane)[iPlane]));
    }
  }
  
  AliDebug(1, Form("AliMFT::SetTreeAddress Digits  fLoader->TreeD() = %p\n", fLoader->TreeD()));
  if (fLoader->TreeD() && fLoader->TreeD()->GetBranch("Plane_00")) {
    CreateDigits();
    for(Int_t iPlane=0; iPlane<AliMFTConstants::kNDisks; iPlane++) {
      fLoader->TreeD()->SetBranchAddress(Form("Plane_%02d",iPlane), &((*fDigitsPerPlane)[iPlane]));
    }
  }
  
  AliDebug(1, Form("AliMFT::SetTreeAddress RecPoints  fLoader->TreeR() = %p\n", fLoader->TreeR()));
  if (fLoader->TreeR() && fLoader->TreeR()->GetBranch("Plane_00")) {
    CreateRecPoints();
    for(Int_t iPlane=0; iPlane<AliMFTConstants::kNDisks; iPlane++) {
      fLoader->TreeR()->SetBranchAddress(Form("Plane_%02d",iPlane), &((*fRecPointsPerPlane)[iPlane]));
    }
  }
  
  AliDebug(1,"Stop.");
  
}

//====================================================================================================================================================

void AliMFT::SetGeometry() {
  
  //  AliInfo("AliMFT::SetGeometry\n");
  //
  //  fSegmentation = new AliMFTSegmentation(fNameGeomFile.Data());
  //
  //  fNPlanes = fSegmentation->GetNPlanes();
  
  fNPlanes = AliMFTConstants::kNDisks;

}

//====================================================================================================================================================

void AliMFT::CreateHits() { 
  
  // create array of hits
  
  AliDebug(1, "AliMFT::CreateHits()");
  
  if (fHits) return;    
  fHits = new TClonesArray("AliMFTHit");  
  
}

//====================================================================================================================================================

void AliMFT::CreateSDigits() { 
  
  // create sdigits list
  
  AliDebug(1, "AliMFT::CreateSDigits()");
  
  if (fSDigitsPerPlane) return; 
  fSDigitsPerPlane = new TObjArray(AliMFTConstants::kNDisks);
  for (Int_t iPlane=0; iPlane<AliMFTConstants::kNDisks; iPlane++) fSDigitsPerPlane->AddAt(new TClonesArray("AliMFTDigit"), iPlane);
  
  fSideDigits = new TClonesArray("AliMFTDigit");
  
}

//====================================================================================================================================================

void AliMFT::CreateDigits() {
  
  // create digits list
  
  AliDebug(1, "AliMFT::CreateDigits()");
  
  if (fDigitsPerPlane) return; 
  fDigitsPerPlane = new TObjArray(AliMFTConstants::kNDisks);
  for(Int_t iPlane=0; iPlane<AliMFTConstants::kNDisks; iPlane++) fDigitsPerPlane->AddAt(new TClonesArray("AliMFTDigit"), iPlane);
  
}

//====================================================================================================================================================

void AliMFT::CreateRecPoints() {
  
  // create recPoints list
  
  AliDebug(1, "AliMFT::CreateRecPoints()");
  
  if (fRecPointsPerPlane) return; 
  fRecPointsPerPlane = new TObjArray(AliMFTConstants::kNDisks);
  for(Int_t iPlane=0; iPlane<AliMFTConstants::kNDisks; iPlane++) fRecPointsPerPlane->AddAt(new TClonesArray("AliMFTCluster"), iPlane);
  
}

//====================================================================================================================================================
