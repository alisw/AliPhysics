/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id: AliEMCALEMCGeometry.cxx 29514 2008-10-26 10:24:38Z hristov $*/

//_________________________________________________________________________
// Geometry class  for EMCAL : singleton  
// EMCAL consists of layers of scintillator and lead
// with scintillator fiber arranged as "shish-kebab" skewers 
// Places the the Barrel Geometry of The EMCAL at Midrapidity
// between 80 and 180(or 190) degrees of Phi and
// -0.7 to 0.7 in eta 
//
//     EMCAL geometry tree:
//     EMCAL -> superModule -> module -> tower(cell)
//     Indexes
//     absId -> nSupMod     -> nModule -> (nIphi,nIeta)
//
//   Name choices: 
//   EMCAL_PDC06 (geometry used for PDC06 simulations, kept for backward compatibility)
//      = equivalent to SHISH_77_TRD1_2X2_FINAL_110DEG in old notation
//   EMCAL_COMPLETE (geometry for expected complete detector)
//      = equivalent to SHISH_77_TRD1_2X2_FINAL_110DEG scTh=0.176 pbTh=0.144
//          in old notation
//
//   EMCAL_FIRSTYEAR - geometry for December 2009 to December 2010 run period
//                     with four Super Modules
//  
//   Adding V1 (EMCAL_FIRSTYEARV1, EMCAL_COMPLETEV1) - geometry from December 2009 ; 
//                1. Fixed bug for positions of modules inside SM
//                   (first module has tilt 0.75 degree);
//                2. Added Al front plate (width 1 cm) and 2 paper sheets per sampling
//                   layer (additional 0.2 mm) 
//                   The sizes have updated with last information from production
//                   drawing (end of October 2010). 
//                3. COMPLETEV1 contains now only 10 SM for runs from 2011
//
//   EMCAL_WSUC (Wayne State test stand)
//      = no definite equivalent in old notation, was only used by
//          Aleksei, but kept for testing purposes
//
//   etc.
//
//
//
//*-- Author: Sahal Yacoob (LBL / UCT)
//     and  : Yves Schutz (SUBATECH)
//     and  : Jennifer Klay (LBL)
//     and  : Aleksei Pavlinov (WSU) 
//     and  : Magali Estienne (SUBATECH)

// --- Root header files ---
#include <TObjArray.h>
#include <TObjString.h>
#include <TRegexp.h>

// -- ALICE Headers.
#include "AliLog.h"

// --- EMCAL headers
#include "AliEMCALEMCGeometry.h"
#include <cassert>

ClassImp(AliEMCALEMCGeometry)

// these initialisations are needed for a singleton
Bool_t    AliEMCALEMCGeometry::fgInit      = kFALSE;
const Char_t*   AliEMCALEMCGeometry::fgkDefaultGeometryName = "EMCAL_COMPLETEV1";


AliEMCALEMCGeometry::AliEMCALEMCGeometry() 
  : TNamed(),
    fGeoName(0),fArrayOpts(0),fNAdditionalOpts(0),fECPbRadThickness(0.),fECScintThick(0.),
    fNECLayers(0),fArm1PhiMin(0.),fArm1PhiMax(0.),fArm1EtaMin(0.),fArm1EtaMax(0.),fIPDistance(0.),
    fShellThickness(0.),fZLength(0.),fNZ(0),fNPhi(0),fSampling(0.),fNumberOfSuperModules(0),
    fFrontSteelStrip(0.),fLateralSteelStrip(0.),fPassiveScintThick(0.),fPhiModuleSize(0.),
    fEtaModuleSize(0.),fPhiTileSize(0.),fEtaTileSize(0.),fLongModuleSize(0.),fNPhiSuperModule(0),
    fNPHIdiv(0),fNETAdiv(0), fNCells(0),fNCellsInSupMod(0),fNCellsInModule(0),
    // Trigger staff
    fNTRUEta(0), fNTRUPhi(0), fNModulesInTRUEta(0), fNModulesInTRUPhi(0), fNEtaSubOfTRU(0),
    // 
    fTrd1Angle(0.),f2Trd1Dx2(0.),
    fPhiGapForSM(0.),fKey110DEG(0),fPhiBoundariesOfSM(0), fPhiCentersOfSM(0),fEtaMaxOfTRD1(0),
    fTrd1AlFrontThick(0.0), fTrd1BondPaperThick(0.),
    fCentersOfCellsEtaDir(0), fCentersOfCellsXDir(0),fCentersOfCellsPhiDir(0),
    fEtaCentersOfCells(0),fPhiCentersOfCells(0),fShishKebabTrd1Modules(0),
    fParSM(), fILOSS(-1), fIHADR(-1),
    //obsolete member data
     fGap2Active(0.), fSteelFrontThick(0.), fTrd2AngleY(0.),
    f2Trd2Dy2(0.), fEmptySpace(0.), fTubsR(0.), fTubsTurnAngle(0.)
{ 
  // default ctor only for internal usage (singleton)
  // must be kept public for root persistency purposes, 
  // but should never be called by the outside world    
  fParSM[0]=0; fParSM[1]=0; fParSM[2]=0;
  fEnvelop[0] = 0; fEnvelop[1] = 0; fEnvelop[2] = 0;
  for(Int_t i = 0; i < 6; i++) fkAdditionalOpts[i] = "";
  
  AliDebug(2, "AliEMCALEMCGeometry : default ctor ");
}
//______________________________________________________________________
AliEMCALEMCGeometry::AliEMCALEMCGeometry(const Text_t* name, const Text_t* title,
                                         const Text_t* mcname, const Text_t* mctitle ) :
  TNamed(name,title),
    fGeoName(0),fArrayOpts(0),fNAdditionalOpts(0),fECPbRadThickness(0.),fECScintThick(0.),
    fNECLayers(0),fArm1PhiMin(0.),fArm1PhiMax(0.),fArm1EtaMin(0.),fArm1EtaMax(0.),fIPDistance(0.),
    fShellThickness(0.),fZLength(0.),fNZ(0),fNPhi(0),fSampling(0.),fNumberOfSuperModules(0),
    fFrontSteelStrip(0.),fLateralSteelStrip(0.),fPassiveScintThick(0.),fPhiModuleSize(0.),
    fEtaModuleSize(0.),fPhiTileSize(0.),fEtaTileSize(0.),fLongModuleSize(0.),fNPhiSuperModule(0),
    fNPHIdiv(0),fNETAdiv(0), fNCells(0),fNCellsInSupMod(0),fNCellsInModule(0),
    // Trigger staff
    fNTRUEta(0), fNTRUPhi(0), fNModulesInTRUEta(0), fNModulesInTRUPhi(0), fNEtaSubOfTRU(0),
    // 
    fTrd1Angle(0.),f2Trd1Dx2(0.),
    fPhiGapForSM(0.),fKey110DEG(0),fPhiBoundariesOfSM(0), fPhiCentersOfSM(0), fEtaMaxOfTRD1(0),
    fTrd1AlFrontThick(0.0), fTrd1BondPaperThick(0.),
    fCentersOfCellsEtaDir(0),fCentersOfCellsXDir(0),fCentersOfCellsPhiDir(0),
    fEtaCentersOfCells(0),fPhiCentersOfCells(0),fShishKebabTrd1Modules(0),
    fParSM(),fILOSS(-1), fIHADR(-1), 
    //obsolete member data
    fGap2Active(0.), fSteelFrontThick(0.), fTrd2AngleY(0.),
    f2Trd2Dy2(0.), fEmptySpace(0.), fTubsR(0.), fTubsTurnAngle(0.)
{
  // ctor only for internal usage (singleton)
  AliDebug(2, Form("AliEMCALEMCGeometry(%s,%s,%s,%s) ", name,title,mcname,mctitle));

  Init(mcname,mctitle);

  //  CreateListOfTrd1Modules();

  if (AliDebugLevel()>=2) {
    PrintGeometry();
  }

}
//______________________________________________________________________
AliEMCALEMCGeometry::AliEMCALEMCGeometry(const AliEMCALEMCGeometry& geom)
  : TNamed(geom),
    fGeoName(geom.fGeoName),
    fArrayOpts(geom.fArrayOpts),
    fNAdditionalOpts(geom.fNAdditionalOpts),
    fECPbRadThickness(geom.fECPbRadThickness),
    fECScintThick(geom.fECScintThick),
    fNECLayers(geom.fNECLayers),
    fArm1PhiMin(geom.fArm1PhiMin),
    fArm1PhiMax(geom.fArm1PhiMax),
    fArm1EtaMin(geom.fArm1EtaMin),
    fArm1EtaMax(geom.fArm1EtaMax),
    fIPDistance(geom.fIPDistance),
    fShellThickness(geom.fShellThickness),
    fZLength(geom.fZLength),
    fNZ(geom.fNZ),
    fNPhi(geom.fNPhi),
    fSampling(geom.fSampling),
    fNumberOfSuperModules(geom.fNumberOfSuperModules),
    fFrontSteelStrip(geom.fFrontSteelStrip),
    fLateralSteelStrip(geom.fLateralSteelStrip),
    fPassiveScintThick(geom.fPassiveScintThick),
    fPhiModuleSize(geom.fPhiModuleSize),
    fEtaModuleSize(geom.fEtaModuleSize),
    fPhiTileSize(geom.fPhiTileSize),
    fEtaTileSize(geom.fEtaTileSize),
    fLongModuleSize(geom.fLongModuleSize),
    fNPhiSuperModule(geom.fNPhiSuperModule),
    fNPHIdiv(geom.fNPHIdiv),
    fNETAdiv(geom.fNETAdiv),
    fNCells(geom.fNCells),
    fNCellsInSupMod(geom.fNCellsInSupMod),
    fNCellsInModule(geom.fNCellsInModule),
    // Trigger staff
    fNTRUEta(geom.fNTRUEta),
    fNTRUPhi(geom.fNTRUPhi),
    fNModulesInTRUEta(geom.fNModulesInTRUEta),
    fNModulesInTRUPhi(geom.fNModulesInTRUPhi),
    fNEtaSubOfTRU(geom.fNEtaSubOfTRU),
    //
    fTrd1Angle(geom.fTrd1Angle),
    f2Trd1Dx2(geom.f2Trd1Dx2),
    fPhiGapForSM(geom.fPhiGapForSM),
    fKey110DEG(geom.fKey110DEG),
    fPhiBoundariesOfSM(geom.fPhiBoundariesOfSM),
    fPhiCentersOfSM(geom.fPhiCentersOfSM),
    fEtaMaxOfTRD1(geom.fEtaMaxOfTRD1),
    fTrd1AlFrontThick(geom.fTrd1AlFrontThick),
    fTrd1BondPaperThick(geom.fTrd1BondPaperThick),
    fCentersOfCellsEtaDir(geom.fCentersOfCellsEtaDir),
    fCentersOfCellsXDir(geom.fCentersOfCellsXDir),
    fCentersOfCellsPhiDir(geom.fCentersOfCellsPhiDir),
    fEtaCentersOfCells(geom.fEtaCentersOfCells),
    fPhiCentersOfCells(geom.fPhiCentersOfCells),
    fShishKebabTrd1Modules(geom.fShishKebabTrd1Modules),
    fILOSS(geom.fILOSS), fIHADR(geom.fIHADR),
    //obsolete member data
    fGap2Active(geom.fGap2Active),
    fSteelFrontThick(geom.fSteelFrontThick),
    fTrd2AngleY(geom.fTrd2AngleY),
    f2Trd2Dy2(geom.f2Trd2Dy2),
    fEmptySpace(geom.fEmptySpace),
    fTubsR(geom.fTubsR),
    fTubsTurnAngle(geom.fTubsTurnAngle)
{
  //copy ctor
  fParSM[0]=geom.fParSM[0]; 
  fParSM[1]=geom.fParSM[1]; 
  fParSM[2]=geom.fParSM[2];
  fEnvelop[0] = geom.fEnvelop[0]; 
  fEnvelop[1] = geom.fEnvelop[1]; 
  fEnvelop[2] = geom.fEnvelop[2];
  for(Int_t i = 0; i < 6; i++) fkAdditionalOpts[i] = geom.fkAdditionalOpts[i];

}

//______________________________________________________________________
AliEMCALEMCGeometry::~AliEMCALEMCGeometry(void){
    // dtor
}

//______________________________________________________________________
void AliEMCALEMCGeometry::Init(const Text_t* mcname, const Text_t* mctitle){
  //
  // Initializes the EMCAL parameters based on the name
  // Only Shashlyk geometry is available, but various combinations of
  // layers and number of supermodules can be selected with additional
  // options or geometry name
  //

  fkAdditionalOpts[0] = "nl=";       // number of sampling layers (fNECLayers)
  fkAdditionalOpts[1] = "pbTh=";     // cm, Thickness of the Pb   (fECPbRadThick)
  fkAdditionalOpts[2] = "scTh=";     // cm, Thickness of the Sc    (fECScintThick)
  fkAdditionalOpts[3] = "latSS=";    // cm, Thickness of lateral steel strip (fLateralSteelStrip)
  fkAdditionalOpts[4] = "allILOSS="; // = 0,1,2,3,4 (4 - energy loss without fluctuation)
  fkAdditionalOpts[5] = "allIHADR="; // = 0,1,2 (0 - no hadronic interaction)

  fNAdditionalOpts = sizeof(fkAdditionalOpts) / sizeof(char*);

  // geometry
  fgInit = kFALSE; // Assume failed until proven otherwise.
  fGeoName   = GetName();
  fGeoName.ToUpper();

  //Convert old geometry names to new ones
  if(fGeoName.Contains("SHISH_77_TRD1_2X2_FINAL_110DEG")) {
    if(fGeoName.Contains("PBTH=0.144") && fGeoName.Contains("SCTH=0.176")) {
      fGeoName = "EMCAL_COMPLETE";
    } else {
      fGeoName = "EMCAL_PDC06";
    }
  }
  if(fGeoName.Contains("WSUC")) fGeoName = "EMCAL_WSUC";

  //check that we have a valid geometry name
  if(!(fGeoName.Contains("EMCAL_PDC06") || fGeoName.Contains("EMCAL_COMPLETE") || fGeoName.Contains("EMCAL_WSUC") 
       || fGeoName.Contains("EMCAL_FIRSTYEAR") || fGeoName.Contains("EMCAL_FIRSTYEARV1") || fGeoName.Contains("EMCAL_COMPLETEV1"))) {
    Fatal("Init", "%s is an undefined geometry!", fGeoName.Data()) ; 
  }

  // Option to know whether we have the "half" supermodule(s) or not
  fKey110DEG = 0;
  if(fGeoName.Contains("COMPLETE") || fGeoName.Contains("PDC06")) fKey110DEG = 1; // for GetAbsCellId
  if(fGeoName.Contains("COMPLETEV1"))  fKey110DEG = 0; 
  fShishKebabTrd1Modules = 0;

  // JLK 13-Apr-2008
  //default parameters are those of EMCAL_COMPLETE geometry
  //all others render variations from these at the end of
  //geometry-name specific options

  fNumberOfSuperModules = 12;       // 12 = 6 * 2 (6 in phi, 2 in Z)
  fNPhi                 = 12;	    // module granularity in phi within smod (azimuth)
  fNZ                   = 24;       // module granularity along Z within smod (eta)
  fNPHIdiv = fNETAdiv   = 2;        // tower granularity within module
  fArm1PhiMin           = 80.0;	    // degrees, Starting EMCAL Phi position
  fArm1PhiMax           = 200.0;    // degrees, Ending EMCAL Phi position
  fArm1EtaMin           = -0.7;	    // pseudorapidity, Starting EMCAL Eta position
  fArm1EtaMax           = +0.7;	    // pseudorapidity, Ending EMCAL Eta position
  fIPDistance           = 428.0;    // cm, radial distance to front face from nominal vertex point
  fPhiGapForSM          = 2.;       // cm, only for final TRD1 geometry
  fFrontSteelStrip      = 0.025;    // 0.025cm = 0.25mm  (13-may-05 from V.Petrov)
  fPassiveScintThick    = 0.8;      // 0.8cm   = 8mm     (13-may-05 from V.Petrov)
  fLateralSteelStrip    = 0.01;     // 0.01cm  = 0.1mm   (13-may-05 from V.Petrov) - was 0.025
  fTrd1Angle            = 1.5;      // in degrees	

  fSampling             = 1.;       // should be calculated with call to DefineSamplingFraction()
  fNECLayers            = 77;       // (13-may-05 from V.Petrov) - can be changed with additional options
  fECScintThick         = 0.176;    // scintillator layer thickness
  fECPbRadThickness     = 0.144;    // lead layer thickness

  fPhiModuleSize = 12.26 - fPhiGapForSM / Float_t(fNPhi); // first assumption
  fEtaModuleSize = fPhiModuleSize;

  fZLength              = 700.;     // Z coverage (cm)


  //needs to be called for each geometry and before setting geometry
  //parameters which can depend on the outcome
  CheckAdditionalOptions();

  //modifications to the above for PDC06 geometry
  if(fGeoName.Contains("PDC06")){ // 18-may-05 - about common structure
    fECScintThick  = fECPbRadThickness = 0.16;// (13-may-05 from V.Petrov)    
    CheckAdditionalOptions();
  }

  //modifications to the above for WSUC geometry
  if(fGeoName.Contains("WSUC")){ // 18-may-05 - about common structure
    fNumberOfSuperModules = 2;  // 27-may-05; Nov 24,2010 for TB
    fNPhi = fNZ = 4; 
    fTrd1AlFrontThick   = 1.0;  // one cm
    // Bond paper - two sheets around Sc tile
    fTrd1BondPaperThick = 0.01; // 0.01cm = 0.1 mm
    
    fPhiModuleSize = 12.0;
    fEtaModuleSize = fPhiModuleSize;
    fLateralSteelStrip = 0.015; // 0.015cm  = 0.15mm

    CheckAdditionalOptions();
  }

  //In 2009-2010 data taking runs only 4 SM, in the upper position.
  if(fGeoName.Contains("FIRSTYEAR")){	
    fNumberOfSuperModules = 4;	
    fArm1PhiMax           = 120.0; 
    CheckAdditionalOptions();	
  }	
  
  if(fGeoName.Contains("FIRSTYEARV1") || fGeoName.Contains("COMPLETEV1") ){
    // Oct 26,2010 : First module has tilt = 0.75 degree : 
    // look to AliEMCALShishKebabTrd1Module::DefineFirstModule(key)
    // New sizes from production drawing, added Al front plate.
    // The thickness of sampling is change due to existing two sheets of paper.
    
    // Will replace fFrontSteelStrip
    fTrd1AlFrontThick   = 1.0;  // one cm
    // Bond paper - two sheets around Sc tile
    fTrd1BondPaperThick = 0.01; // 0.01cm = 0.1 mm
    
    fPhiModuleSize = 12.0;
    fEtaModuleSize = fPhiModuleSize;
    fLateralSteelStrip = 0.015; // 0.015cm  = 0.15mm
    
    if(fGeoName.Contains("COMPLETEV1")){
      fNumberOfSuperModules = 10;	
      fArm1PhiMax           = 180.0; 
    }
    
    CheckAdditionalOptions();	
  }

  // constant for transition absid <--> indexes
  fNCellsInModule = fNPHIdiv*fNETAdiv;
  fNCellsInSupMod = fNCellsInModule*fNPhi*fNZ;
  fNCells         = fNCellsInSupMod*fNumberOfSuperModules;
  if(GetKey110DEG()) fNCells -= fNCellsInSupMod;

  fNPhiSuperModule = fNumberOfSuperModules/2;
  if(fNPhiSuperModule < 1) fNPhiSuperModule = 1;
    
  fPhiTileSize = fPhiModuleSize/double(fNPHIdiv) - fLateralSteelStrip; // 13-may-05 
  fEtaTileSize = fEtaModuleSize/double(fNETAdiv) - fLateralSteelStrip; // 13-may-05 

  fLongModuleSize = fNECLayers*(fECScintThick + fECPbRadThickness);  
  if(fGeoName.Contains("V1")){
    Double_t ws = fECScintThick + fECPbRadThickness + 2.*fTrd1BondPaperThick; // sampling width
    // Number of Pb tiles = Number of Sc tiles - 1
    fLongModuleSize = fTrd1AlFrontThick + (ws*fNECLayers - fECPbRadThickness);
  }
  f2Trd1Dx2 = fEtaModuleSize + 2.*fLongModuleSize*TMath::Tan(fTrd1Angle*TMath::DegToRad()/2.);

  if(!fGeoName.Contains("WSUC")) fShellThickness  = TMath::Sqrt(fLongModuleSize*fLongModuleSize + f2Trd1Dx2*f2Trd1Dx2);

  //These parameters are used to create the mother volume to hold the supermodules
  //2cm padding added to allow for misalignments - JLK 30-May-2008
  fEnvelop[0]     = fIPDistance - 1.; // mother volume inner radius
  fEnvelop[1]     = fIPDistance + fShellThickness + 1.; // mother volume outer r.
  fEnvelop[2]     = fZLength + 2.; //mother volume length 

  // Local coordinates
  fParSM[0] = GetShellThickness()/2.;        
  fParSM[1] = GetPhiModuleSize() * GetNPhi()/2.;
  fParSM[2] = fZLength/4.;  //divide by 4 to get half-length of SM

  // SM phi boundaries - (0,1),(2,3) .. (10,11) - has the same boundaries; Nov 7, 2006 
  fPhiBoundariesOfSM.Set(fNumberOfSuperModules);
  fPhiCentersOfSM.Set(fNumberOfSuperModules/2);
  fPhiBoundariesOfSM[0] = TMath::PiOver2() - TMath::ATan2(fParSM[1] , fIPDistance); // 1th and 2th modules)
  fPhiCentersOfSM[0]     = TMath::PiOver2();
  if(fNumberOfSuperModules > 1) 
    fPhiBoundariesOfSM[1] = TMath::PiOver2() + TMath::ATan2(fParSM[1] , fIPDistance);
  if(fNumberOfSuperModules > 2) {
	Int_t maxPhiBlock =fNumberOfSuperModules/2-1;
	if(fNumberOfSuperModules > 10) maxPhiBlock = 4;
    for(int i=1; i<=maxPhiBlock; i++) { // from 2th ro 9th
      fPhiBoundariesOfSM[2*i]   = fPhiBoundariesOfSM[0] + 20.*TMath::DegToRad()*i;
      fPhiBoundariesOfSM[2*i+1] = fPhiBoundariesOfSM[1] + 20.*TMath::DegToRad()*i;
      fPhiCentersOfSM[i]        = fPhiCentersOfSM[0]     + 20.*TMath::DegToRad()*i;
    }
  }
  if(fNumberOfSuperModules > 10) {
    fPhiBoundariesOfSM[11] = 190.*TMath::DegToRad();
    fPhiBoundariesOfSM[10] = fPhiBoundariesOfSM[11] - TMath::ATan2((fParSM[1]) , fIPDistance);
    fPhiCentersOfSM[5]     = (fPhiBoundariesOfSM[10]+fPhiBoundariesOfSM[11])/2.; 
  }

  //called after setting of scintillator and lead layer parameters
  DefineSamplingFraction(mcname,mctitle);

  
  // TRU parameters - Apr 29,08 by PAI. 
  // These parameters values was updated at Nov 05, 2007
  // As is on Olivier  BOURRION (LPSC) ppt preasentation 
  // at ALICE trigger meeting at 13th-14th March
  fNTRUEta = 1;           // was 3
  fNTRUPhi = 3;           // was 1
  fNModulesInTRUEta = 24; // was 8
  fNModulesInTRUPhi = 4;  // was 12
  // Jet trigger 
  // 3*6*10 + 2*6*2 = 204 -> matrix (nphi(17), neta(12))
  fNEtaSubOfTRU     = 6;  

  fgInit = kTRUE; 
}

//___________________________________________________________________
void AliEMCALEMCGeometry::PrintGeometry()
{
  // Separate routine is callable from broswer; Nov 7,2006
  printf("\nInit: geometry of EMCAL named %s :\n", fGeoName.Data());
  if(fArrayOpts) {
    for(Int_t i=0; i<fArrayOpts->GetEntries(); i++){
      TObjString *o = (TObjString*)fArrayOpts->At(i);
      printf(" %i : %s \n", i, o->String().Data());
    }
  }
  printf("Granularity: %d in eta and %d in phi\n", GetNZ(), GetNPhi()) ;
  printf("Layout: phi = (%7.1f, %7.1f), eta = (%5.2f, %5.2f), IP = %7.2f -> for EMCAL envelope only\n",  
	   GetArm1PhiMin(), GetArm1PhiMax(),GetArm1EtaMin(), GetArm1EtaMax(), GetIPDistance() );

  printf( "               ECAL      : %d x (%f cm Pb, %f cm Sc) \n", 
  GetNECLayers(), GetECPbRadThick(), GetECScintThick() ) ; 
  printf("                fSampling %5.2f \n",  fSampling );
  printf(" fIPDistance       %6.3f cm \n", fIPDistance);
  printf(" fNPhi %i   |  fNZ %i \n", fNPhi, fNZ);
  printf(" fNCellsInModule %i : fNCellsInSupMod %i : fNCells %i\n",fNCellsInModule, fNCellsInSupMod, fNCells);
  printf(" X:Y module size     %6.3f , %6.3f cm \n", fPhiModuleSize, fEtaModuleSize);
  printf(" X:Y   tile size     %6.3f , %6.3f cm \n", fPhiTileSize, fEtaTileSize);
  printf(" #of sampling layers %i(fNECLayers) \n", fNECLayers);
  printf(" fLongModuleSize     %6.3f cm \n", fLongModuleSize);
  printf(" #supermodule in phi direction %i \n", fNPhiSuperModule );
  printf(" fILOSS %i : fIHADR %i \n", fILOSS, fIHADR);
  printf(" fTrd1Angle %7.4f\n", fTrd1Angle);
  printf(" f2Trd1Dx2  %7.4f\n",  f2Trd1Dx2);
  printf(" fTrd1AlFrontThick   %7.4f \n", fTrd1AlFrontThick);
  printf(" fTrd1BondPaperThick %5.4f \n", fTrd1BondPaperThick);
  printf("SM dimensions(TRD1) : dx %7.2f dy %7.2f dz %7.2f (SMOD, BOX)\n", 
	 fParSM[0],fParSM[1],fParSM[2]);
  printf(" fPhiGapForSM  %7.4f cm (%7.4f <- phi size in degree)\n",  
	 fPhiGapForSM, TMath::ATan2(fPhiGapForSM,fIPDistance)*TMath::RadToDeg());
  if(GetKey110DEG()) printf(" Last two modules have size 10 degree in  phi (180<phi<190)\n");
  printf(" phi SM boundaries \n"); 
  for(int i=0; i<fPhiBoundariesOfSM.GetSize()/2.; i++) {
    printf(" %i : %7.5f(%7.2f) -> %7.5f(%7.2f) : center %7.5f(%7.2f) \n", i, 
	   fPhiBoundariesOfSM[2*i], fPhiBoundariesOfSM[2*i]*TMath::RadToDeg(),
	   fPhiBoundariesOfSM[2*i+1], fPhiBoundariesOfSM[2*i+1]*TMath::RadToDeg(),
	   fPhiCentersOfSM[i], fPhiCentersOfSM[i]*TMath::RadToDeg());
  }

}

//______________________________________________________________________
void AliEMCALEMCGeometry::CheckAdditionalOptions()
{
  // Feb 06,2006
  // Additional options that
  // can be used to select
  // the specific geometry of 
  // EMCAL to run
  // Dec 27,2006
  // adeed allILOSS= and allIHADR= for MIP investigation
  fArrayOpts = new TObjArray;
  Int_t nopt = ParseString(fGeoName, *fArrayOpts);
  if(nopt==1) { // no aditional option(s)
    fArrayOpts->Delete();
    delete fArrayOpts;
    fArrayOpts = 0; 
    return;
  }  		 
  for(Int_t i=1; i<nopt; i++){
    TObjString *o = (TObjString*)fArrayOpts->At(i); 

    TString addOpt = o->String();
    Int_t indj=-1;
    for(Int_t j=0; j<fNAdditionalOpts; j++) {
      TString opt = fkAdditionalOpts[j];
      if(addOpt.Contains(opt,TString::kIgnoreCase)) {
	  indj = j;
        break;
      }
    }
    if(indj<0) {
      AliDebug(2,Form("<E> option |%s| unavailable : ** look to the file AliEMCALGeometry.h **\n", 
		      addOpt.Data()));
      assert(0);
    } else {
      AliDebug(2,Form("<I> option |%s| is valid : number %i : |%s|\n", 
		      addOpt.Data(), indj, fkAdditionalOpts[indj]));
      if       (addOpt.Contains("NL=",TString::kIgnoreCase))   {// number of sampling layers
        sscanf(addOpt.Data(),"NL=%i", &fNECLayers);
        AliDebug(2,Form(" fNECLayers %i (new) \n", fNECLayers));
      } else if(addOpt.Contains("PBTH=",TString::kIgnoreCase)) {//Thickness of the Pb(fECPbRadThicknes)
        sscanf(addOpt.Data(),"PBTH=%f", &fECPbRadThickness);
      } else if(addOpt.Contains("SCTH=",TString::kIgnoreCase)) {//Thickness of the Sc(fECScintThick)
        sscanf(addOpt.Data(),"SCTH=%f", &fECScintThick);
      } else if(addOpt.Contains("LATSS=",TString::kIgnoreCase)) {// Thickness of lateral steel strip (fLateralSteelStrip)
        sscanf(addOpt.Data(),"LATSS=%f", &fLateralSteelStrip);
        AliDebug(2,Form(" fLateralSteelStrip %f (new) \n", fLateralSteelStrip));
      } else if(addOpt.Contains("ILOSS=",TString::kIgnoreCase)) {// As in Geant
        sscanf(addOpt.Data(),"ALLILOSS=%i", &fILOSS);
        AliDebug(2,Form(" fILOSS %i \n", fILOSS));
      } else if(addOpt.Contains("IHADR=",TString::kIgnoreCase)) {// As in Geant
        sscanf(addOpt.Data(),"ALLIHADR=%i", &fIHADR);
        AliDebug(2,Form(" fIHADR %i \n", fIHADR));
      }
    }
  }
}

//__________________________________________________________________
void AliEMCALEMCGeometry::DefineSamplingFraction(const Text_t* mcname, const Text_t* mctitle)
{
  // Jun 05,2006
  // Look http://rhic.physics.wayne.edu/~pavlinov/ALICE/SHISHKEBAB/RES/linearityAndResolutionForTRD1.html
  // Keep for compatibilty
  //
  
  // Sampling factor for G3
  fSampling = 10.87; // Default value - Nov 25,2010
  if(fNECLayers == 69) {        // 10% layer reduction
    fSampling = 12.55;
  } else if(fNECLayers == 61) { // 20% layer reduction
    fSampling = 12.80;
  } else if(fNECLayers == 77) {
    if(fGeoName.Contains("V1")){
      fSampling = 10.87; //Adding paper sheets and cover plate; Nov 25,2010
    } else if   (fECScintThick>0.159 && fECScintThick<0.161) { // original sampling fraction, equal layers
      fSampling = 12.327; // fECScintThick = fECPbRadThickness = 0.160;
    } else if (fECScintThick>0.175 && fECScintThick<0.177) { // 10% Pb thicknes reduction
      fSampling = 10.5; // fECScintThick = 0.176, fECPbRadThickness=0.144;
    } else if(fECScintThick>0.191 && fECScintThick<0.193) { // 20% Pb thicknes reduction
      fSampling = 8.93; // fECScintThick = 0.192, fECPbRadThickness=0.128;
    }
  }

  // Default sampling factor for G3, modify it for other transport model
  TString mcName  = mcname;
  TString mcTitle = mctitle;

  Float_t samplingFactorTranportModel = 1. ;
  if     (mcName.Contains("Geant3")) samplingFactorTranportModel = 1.;//0.988 // Do nothing
  else if(mcName.Contains("Fluka") ) samplingFactorTranportModel = 1.; // To be set
  else if(mcName.Contains("Geant4")){
    if(mcTitle.Contains("EMV"))      samplingFactorTranportModel = 1.096; // 0.906, 0.896 (OPT)
    else                             samplingFactorTranportModel = 0.86; // 1.15 (CHIPS), 1.149 (BERT), 1.147 (BERT_CHIPS) 
  }      
  
  AliDebug(2,Form("MC modeler <%s>, Title <%s>: Sampling %f, model fraction with respect to G3 %f, final sampling %f \n",
               mcName.Data(),mcTitle.Data(),fSampling,samplingFactorTranportModel,fSampling*samplingFactorTranportModel));

  
  fSampling*=samplingFactorTranportModel;
  
}

//________________________________________________________________________________________________
Double_t AliEMCALEMCGeometry::GetPhiCenterOfSM(Int_t nsupmod) const
{
  //returns center of supermodule in phi
  int i = nsupmod/2;
  return fPhiCentersOfSM[i];

}

//________________________________________________________________________________________________
Bool_t AliEMCALEMCGeometry::GetPhiBoundariesOfSM(Int_t nSupMod, Double_t &phiMin, Double_t &phiMax) const
{
  // 0<= nSupMod <=11; phi in rad
    static int i;
  if(nSupMod<0 || nSupMod >11) return kFALSE;
  i = nSupMod/2;
  phiMin = (Double_t)fPhiBoundariesOfSM[2*i];
  phiMax = (Double_t)fPhiBoundariesOfSM[2*i+1];
  return kTRUE;
}

//________________________________________________________________________________________________
Bool_t AliEMCALEMCGeometry::GetPhiBoundariesOfSMGap(Int_t nPhiSec, Double_t &phiMin, Double_t &phiMax) const
{
  // 0<= nPhiSec <=4; phi in rad
  // 0;  gap boundaries between  0th&2th  | 1th&3th SM
  // 1;  gap boundaries between  2th&4th  | 3th&5th SM
  // 2;  gap boundaries between  4th&6th  | 5th&7th SM
  // 3;  gap boundaries between  6th&8th  | 7th&9th SM
  // 4;  gap boundaries between  8th&10th | 9th&11th SM
  if(nPhiSec<0 || nPhiSec >4) return kFALSE;
  phiMin = fPhiBoundariesOfSM[2*nPhiSec+1];
  phiMax = fPhiBoundariesOfSM[2*nPhiSec+2];
  return kTRUE;
}

//________________________________________________________________________________________________
int AliEMCALEMCGeometry::ParseString(const TString &topt, TObjArray &Opt)
{ 
	//Parse string, does what? GCB 08/09
	Ssiz_t begin, index, end, end2;
	begin = index = end = end2 = 0;
	TRegexp separator("[^ ;,\\t\\s/]+");
	while ( (begin < topt.Length()) && (index != kNPOS) ) {
		// loop over given options
		index = topt.Index(separator,&end,begin);
		if (index >= 0 && end >= 1) {
			TString substring(topt(index,end));
			Opt.Add(new TObjString(substring.Data()));
		}
		begin += end+1;
	}
	return Opt.GetEntries();
}

