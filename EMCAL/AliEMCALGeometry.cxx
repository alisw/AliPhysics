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

/* $Id$*/

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
//

#include <cassert>

// --- Root header files ---
#include <Riostream.h>
#include <TBrowser.h>
#include <TClonesArray.h>
#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include <TGeoNode.h>
#include <TList.h>
#include <TMatrixD.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TVector2.h>
#include <TVector3.h>
#include <TParticle.h>
// -- ALICE Headers.
#include "AliLog.h"

// --- EMCAL headers
#include "AliEMCALGeometry.h"
#include "AliEMCALShishKebabTrd1Module.h"
#include "AliEMCALRecPoint.h"
#include "AliEMCALDigit.h"
#include "AliEMCALHistoUtilities.h"

ClassImp(AliEMCALGeometry)

// these initialisations are needed for a singleton
AliEMCALGeometry  *AliEMCALGeometry::fgGeom      = 0;
Bool_t             AliEMCALGeometry::fgInit      = kFALSE;
const Char_t*            AliEMCALGeometry::fgDefaultGeometryName = "EMCAL_COMPLETE";
//
// Usage: 
//        You can create the AliEMCALGeometry object independently from anything.
//        You have to use just the correct name of geometry. If name is empty string the
//        default name of geometry will be used.
//         
//  AliEMCALGeometry* g = AliEMCALGeometry::GetInstance(name,title); // first time
//  ..
//  g = AliEMCALGeometry::GetInstance();                             // after first time
//
//  MC:   If you work with MC data you have to get geometry the next way: 
//  ==                                      =============================
//  AliRunLoader    *rl   = AliRunLoader::Instance();
//  AliEMCALGeometry *geom = dynamic_cast<AliEMCAL*>(rl->GetAliRun()->GetDetector("EMCAL"))->GetGeometry();
//  TGeoManager::Import("geometry.root");

AliEMCALGeometry::AliEMCALGeometry() 
  : AliGeometry(),
    fGeoName(0),fArrayOpts(0),fNAdditionalOpts(0),fECPbRadThickness(0.),fECScintThick(0.),
    fNECLayers(0),fArm1PhiMin(0.),fArm1PhiMax(0.),fArm1EtaMin(0.),fArm1EtaMax(0.),fIPDistance(0.),
    fShellThickness(0.),fZLength(0.),fNZ(0),fNPhi(0),fSampling(0.),fNumberOfSuperModules(0),
    fFrontSteelStrip(0.),fLateralSteelStrip(0.),fPassiveScintThick(0.),fPhiModuleSize(0.),
    fEtaModuleSize(0.),fPhiTileSize(0.),fEtaTileSize(0.),fLongModuleSize(0.),fNPhiSuperModule(0),
    fNPHIdiv(0),fNETAdiv(0), fNCells(0),fNCellsInSupMod(0),fNCellsInModule(0),
    // Trigger staff
    fNTRUEta(0),
    fNTRUPhi(0), 
    fNModulesInTRUEta(0), 
    fNModulesInTRUPhi(0), 
    fNEtaSubOfTRU(0),
    // 
    fTrd1Angle(0.),f2Trd1Dx2(0.),
    fPhiGapForSM(0.),fKey110DEG(0),fPhiBoundariesOfSM(0), fPhiCentersOfSM(0),fEtaMaxOfTRD1(0),
    fCentersOfCellsEtaDir(0), fCentersOfCellsXDir(0),fCentersOfCellsPhiDir(0),
    fEtaCentersOfCells(0),fPhiCentersOfCells(0),fShishKebabTrd1Modules(0),
    fILOSS(-1), fIHADR(-1),
    //obsolete member data
    fAlFrontThick(0.), fGap2Active(0.), fSteelFrontThick(0.), fTrd2AngleY(0.),
    f2Trd2Dy2(0.), fEmptySpace(0.), fTubsR(0.), fTubsTurnAngle(0.)
{ 
  // default ctor only for internal usage (singleton)
  // must be kept public for root persistency purposes, 
  // but should never be called by the outside world    

  AliDebug(2, "AliEMCALGeometry : default ctor ");
}
//______________________________________________________________________
AliEMCALGeometry::AliEMCALGeometry(const Text_t* name, const Text_t* title) 
  : AliGeometry(name, title),
    fGeoName(0),fArrayOpts(0),fNAdditionalOpts(0),fECPbRadThickness(0.),fECScintThick(0.),
    fNECLayers(0),fArm1PhiMin(0.),fArm1PhiMax(0.),fArm1EtaMin(0.),fArm1EtaMax(0.),fIPDistance(0.),
    fShellThickness(0.),fZLength(0.),fNZ(0),fNPhi(0),fSampling(0.),fNumberOfSuperModules(0),
    fFrontSteelStrip(0.),fLateralSteelStrip(0.),fPassiveScintThick(0.),fPhiModuleSize(0.),
    fEtaModuleSize(0.),fPhiTileSize(0.),fEtaTileSize(0.),fLongModuleSize(0.),fNPhiSuperModule(0),
    fNPHIdiv(0),fNETAdiv(0), fNCells(0),fNCellsInSupMod(0),fNCellsInModule(0),
    // Trigger staff
    fNTRUEta(0),
    fNTRUPhi(0), 
    fNModulesInTRUEta(0), 
    fNModulesInTRUPhi(0), 
    fNEtaSubOfTRU(0),
    // 
    fTrd1Angle(0.),f2Trd1Dx2(0.),
    fPhiGapForSM(0.),fKey110DEG(0),fPhiBoundariesOfSM(0), fPhiCentersOfSM(0), fEtaMaxOfTRD1(0),
    fCentersOfCellsEtaDir(0),fCentersOfCellsXDir(0),fCentersOfCellsPhiDir(0),
    fEtaCentersOfCells(0),fPhiCentersOfCells(0),fShishKebabTrd1Modules(0),
    fILOSS(-1), fIHADR(-1), 
    //obsolete member data
    fAlFrontThick(0.), fGap2Active(0.), fSteelFrontThick(0.), fTrd2AngleY(0.),
    f2Trd2Dy2(0.), fEmptySpace(0.), fTubsR(0.), fTubsTurnAngle(0.)
{
  // ctor only for internal usage (singleton)
  AliDebug(2, Form("AliEMCALGeometry(%s,%s) ", name,title));

  Init();

  CreateListOfTrd1Modules();

  if (AliDebugLevel()>=2) {
    PrintGeometry();
  }

}
//______________________________________________________________________
AliEMCALGeometry::AliEMCALGeometry(const AliEMCALGeometry& geom)
  : AliGeometry(geom),
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
    fCentersOfCellsEtaDir(geom.fCentersOfCellsEtaDir),
    fCentersOfCellsXDir(geom.fCentersOfCellsXDir),
    fCentersOfCellsPhiDir(geom.fCentersOfCellsPhiDir),
    fEtaCentersOfCells(geom.fEtaCentersOfCells),
    fPhiCentersOfCells(geom.fPhiCentersOfCells),
    fShishKebabTrd1Modules(geom.fShishKebabTrd1Modules),
    fILOSS(geom.fILOSS), fIHADR(geom.fIHADR),
    //obsolete member data
    fAlFrontThick(geom.fAlFrontThick),
    fGap2Active(geom.fGap2Active),
    fSteelFrontThick(geom.fSteelFrontThick),
    fTrd2AngleY(geom.fTrd2AngleY),
    f2Trd2Dy2(geom.f2Trd2Dy2),
    fEmptySpace(geom.fEmptySpace),
    fTubsR(geom.fTubsR),
    fTubsTurnAngle(geom.fTubsTurnAngle)
{
  //copy ctor
}

//______________________________________________________________________
AliEMCALGeometry::~AliEMCALGeometry(void){
    // dtor
}

//______________________________________________________________________
void AliEMCALGeometry::Init(void){
  //
  // Initializes the EMCAL parameters based on the name
  // Only Shashlyk geometry is available, but various combinations of
  // layers and number of supermodules can be selected with additional
  // options or geometry name
  //

  fAdditionalOpts[0] = "nl=";       // number of sampling layers (fNECLayers)
  fAdditionalOpts[1] = "pbTh=";     // cm, Thickness of the Pb   (fECPbRadThick)
  fAdditionalOpts[2] = "scTh=";     // cm, Thickness of the Sc    (fECScintThick)
  fAdditionalOpts[3] = "latSS=";    // cm, Thickness of lateral steel strip (fLateralSteelStrip)
  fAdditionalOpts[4] = "allILOSS="; // = 0,1,2,3,4 (4 - energy loss without fluctuation)
  fAdditionalOpts[5] = "allIHADR="; // = 0,1,2 (0 - no hadronic interaction)

  fNAdditionalOpts = sizeof(fAdditionalOpts) / sizeof(char*);

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
  if(!(fGeoName.Contains("EMCAL_PDC06") || fGeoName.Contains("EMCAL_COMPLETE") || fGeoName.Contains("EMCAL_WSUC") || fGeoName.Contains("EMCAL_1stYear"))) {
    Fatal("Init", "%s is an undefined geometry!", fGeoName.Data()) ; 
  }

  // Option to know whether we have the "half" supermodule(s) or not
  fKey110DEG = 0;
  if(fGeoName.Contains("COMPLETE") || fGeoName.Contains("PDC06")) fKey110DEG = 1; // for GetAbsCellId
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
    fPhiModuleSize = 12.5;     // 20-may-05 - rectangular shape
    fEtaModuleSize = 11.9;
    fECScintThick  = fECPbRadThickness = 0.16;// (13-may-05 from V.Petrov)
    fNumberOfSuperModules = 1; // 27-may-05
    fShellThickness = 30.;       // should be change 
    fNPhi = fNZ = 4; 
    CheckAdditionalOptions();
  }

  if(fGeoName.Contains("1stYear")){	
	fNumberOfSuperModules = 2;	
	 
	if(fGeoName.Contains("LowerEta")) {
		fNPhiSuperModule = 1;		
	}
	else if(fGeoName.Contains("LowerPhi_SideA")){
	fNPhiSuperModule = 2;	
	fArm1EtaMax=0;		
	}
	else if(fGeoName.Contains("LowerPhi_SideC")){
	fNPhiSuperModule = 2;		
	fArm1EtaMin=0;	
	}
		
      CheckAdditionalOptions();	
  }	

  // constant for transition absid <--> indexes
  fNCellsInModule  = fNPHIdiv*fNETAdiv;
  fNCellsInSupMod = fNCellsInModule*fNPhi*fNZ;
  fNCells         = fNCellsInSupMod*fNumberOfSuperModules;
  if(GetKey110DEG()) fNCells -= fNCellsInSupMod;

  fNPhiSuperModule = fNumberOfSuperModules/2;
  if(fNPhiSuperModule < 1) fNPhiSuperModule = 1;
    
  fPhiTileSize = fPhiModuleSize/double(fNPHIdiv) - fLateralSteelStrip; // 13-may-05 
  fEtaTileSize = fEtaModuleSize/double(fNETAdiv) - fLateralSteelStrip; // 13-may-05 

  fLongModuleSize = fNECLayers*(fECScintThick + fECPbRadThickness);  
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
    for(int i=1; i<=4; i++) { // from 2th ro 9th
      fPhiBoundariesOfSM[2*i]   = fPhiBoundariesOfSM[0] + 20.*TMath::DegToRad()*i;
      fPhiBoundariesOfSM[2*i+1] = fPhiBoundariesOfSM[1] + 20.*TMath::DegToRad()*i;
      fPhiCentersOfSM[i]         = fPhiCentersOfSM[0]     + 20.*TMath::DegToRad()*i;
    }
  }
  if(fNumberOfSuperModules > 10) {
    fPhiBoundariesOfSM[11] = 190.*TMath::DegToRad();
    fPhiBoundariesOfSM[10] = fPhiBoundariesOfSM[11] - TMath::ATan2((fParSM[1]) , fIPDistance);
    fPhiCentersOfSM[5]      = (fPhiBoundariesOfSM[10]+fPhiBoundariesOfSM[11])/2.; 
  }

  //called after setting of scintillator and lead layer parameters
  DefineSamplingFraction();

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
void AliEMCALGeometry::PrintGeometry()
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
  printf(" fShishKebabTrd1Modules has %i modules : max eta %5.4f \n", 
	 fShishKebabTrd1Modules->GetSize(),fEtaMaxOfTRD1);
  
  printf("\n Cells grid in eta directions : size %i\n", fCentersOfCellsEtaDir.GetSize());
  for(Int_t i=0; i<fCentersOfCellsEtaDir.GetSize(); i++) {
    printf(" ind %2.2i : z %8.3f : x %8.3f \n", i, 
	   fCentersOfCellsEtaDir.At(i),fCentersOfCellsXDir.At(i));
    int ind=0; // Nov 21,2006
    for(Int_t iphi=0; iphi<fCentersOfCellsPhiDir.GetSize(); iphi++) {
      ind = iphi*fCentersOfCellsEtaDir.GetSize() + i;
      printf("%6.4f ", fEtaCentersOfCells[ind]);
      if((iphi+1)%12 == 0) printf("\n");
    }
    printf("\n");
    
  }

  printf("\n Cells grid in phi directions : size %i\n", fCentersOfCellsPhiDir.GetSize());
  for(Int_t i=0; i<fCentersOfCellsPhiDir.GetSize(); i++) {
    double phi=fPhiCentersOfCells.At(i);
    printf(" ind %2.2i : y %8.3f : phi %7.5f(%6.2f) \n", i, fCentersOfCellsPhiDir.At(i), 
	   phi, phi*TMath::RadToDeg());
  }

}

//______________________________________________________________________
void AliEMCALGeometry::PrintCellIndexes(Int_t absId, int pri, const char *tit)
{
  // Service methods
  Int_t nSupMod, nModule, nIphi, nIeta;
  Int_t iphi, ieta;
  TVector3 vg;

  GetCellIndex(absId,  nSupMod, nModule, nIphi, nIeta);
  printf(" %s | absId : %i -> nSupMod %i nModule %i nIphi %i nIeta %i \n", tit, absId,  nSupMod, nModule, nIphi, nIeta);
  if(pri>0) {
    GetCellPhiEtaIndexInSModule(nSupMod,nModule,nIphi,nIeta, iphi,ieta);
    printf(" local SM index : iphi %i : ieta %i \n", iphi,ieta);
    GetGlobal(absId, vg);
    printf(" vglob : mag %7.2f : perp %7.2f : z %7.2f : eta %6.4f : phi %6.4f(%6.2f) \n", 
	   vg.Mag(), vg.Perp(), vg.Z(), vg.Eta(), vg.Phi(), vg.Phi()*TMath::RadToDeg());
  }
}

//______________________________________________________________________
void AliEMCALGeometry::CheckAdditionalOptions()
{
  // Feb 06,2006
  // Additional options that
  // can be used to select
  // the specific geometry of 
  // EMCAL to run
  // Dec 27,2006
  // adeed allILOSS= and allIHADR= for MIP investigation
  fArrayOpts = new TObjArray;
  Int_t nopt = AliEMCALHistoUtilities::ParseString(fGeoName, *fArrayOpts);
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
      TString opt = fAdditionalOpts[j];
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
		      addOpt.Data(), indj, fAdditionalOpts[indj]));
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
void AliEMCALGeometry::DefineSamplingFraction()
{
  // Jun 05,2006
  // Look http://rhic.physics.wayne.edu/~pavlinov/ALICE/SHISHKEBAB/RES/linearityAndResolutionForTRD1.html
  // Keep for compatibilty
  //
  if(fNECLayers == 69) {        // 10% layer reduction
    fSampling = 12.55;
  } else if(fNECLayers == 61) { // 20% layer reduction
    fSampling = 12.80;
  } else if(fNECLayers == 77) {
    if       (fECScintThick>0.159 && fECScintThick<0.161) { // original sampling fraction, equal layers
      fSampling = 12.327; // fECScintThick = fECPbRadThickness = 0.160;
    } else if (fECScintThick>0.175 && fECScintThick<0.177) { // 10% Pb thicknes reduction
      fSampling = 10.5; // fECScintThick = 0.176, fECPbRadThickness=0.144;
    } else if(fECScintThick>0.191 && fECScintThick<0.193) { // 20% Pb thicknes reduction
      fSampling = 8.93; // fECScintThick = 0.192, fECPbRadThickness=0.128;
    }

  }
}

//______________________________________________________________________
void AliEMCALGeometry::GetModulePhiEtaIndexInSModuleFromTRUIndex(Int_t itru, Int_t iphitru, Int_t ietatru, Int_t &iphiSM, Int_t &ietaSM) const 
{
  
  // This method transforms the (eta,phi) index of module in a 
  // TRU matrix into Super Module (eta,phi) index.
  
  // Calculate in which row and column where the TRU are 
  // ordered in the SM

  Int_t col = itru/ fNTRUPhi ; // indexes of TRU in SM
  Int_t row = itru - col*fNTRUPhi ;
   
  iphiSM = fNModulesInTRUPhi*row + iphitru  ;
  ietaSM = fNModulesInTRUEta*col + ietatru  ; 
  //printf(" GetModulePhiEtaIndexInSModuleFromTRUIndex : itru %2i iphitru %2i ietatru %2i iphiSM %2i ietaSM %2i \n", 
  // itru, iphitru, ietatru, iphiSM, ietaSM);
}

//______________________________________________________________________
AliEMCALGeometry *  AliEMCALGeometry::GetInstance(){ 
  // Returns the pointer of the unique instance
  
  AliEMCALGeometry * rv = static_cast<AliEMCALGeometry *>( fgGeom );
  return rv; 
}

//______________________________________________________________________
AliEMCALGeometry* AliEMCALGeometry::GetInstance(const Text_t* name,
						const Text_t* title){
    // Returns the pointer of the unique instance

    AliEMCALGeometry * rv = 0; 
    if ( fgGeom == 0 ) {
      if ( strcmp(name,"") == 0 ) { // get default geometry
	 fgGeom = new AliEMCALGeometry(fgDefaultGeometryName, title);
      } else {
	 fgGeom = new AliEMCALGeometry(name, title);
      }  // end if strcmp(name,"")
      if ( fgInit ) rv = (AliEMCALGeometry * ) fgGeom;
      else {
	 rv = 0; 
	 delete fgGeom; 
	 fgGeom = 0; 
      } // end if fgInit
    }else{
	if ( strcmp(fgGeom->GetName(), name) != 0) {
	  printf("\ncurrent geometry is %s : ", fgGeom->GetName());
	  printf(" you cannot call %s ",name);  
	}else{
	  rv = (AliEMCALGeometry *) fgGeom; 
	} // end 
    }  // end if fgGeom
    return rv; 
}

//_____________________________________________________________________________
Bool_t AliEMCALGeometry::IsInEMCAL(Double_t x, Double_t y, Double_t z) const {
  // Checks whether point is inside the EMCal volume, used in AliEMCALv*.cxx
  //
  // Code uses cylindrical approximation made of inner radius (for speed)
  //
  // Points behind EMCAl, i.e. R > outer radius, but eta, phi in acceptance 
  // are considered to inside

  Double_t r=sqrt(x*x+y*y);

  if ( r > fEnvelop[0] ) {
     Double_t theta;
     theta  =    TMath::ATan2(r,z);
     Double_t eta;
     if(theta == 0) 
       eta = 9999;
     else 
       eta    =   -TMath::Log(TMath::Tan(theta/2.));
     if (eta < fArm1EtaMin || eta > fArm1EtaMax)
       return 0;
 
     Double_t phi = TMath::ATan2(y,x) * 180./TMath::Pi();
     if (phi < 0) phi += 360;  // phi should go from 0 to 360 in this case
     if (phi > fArm1PhiMin && phi < fArm1PhiMax)
       return 1;
  }
  return 0;
}

//
// == Shish-kebab cases ==
//
//________________________________________________________________________________________________
Int_t AliEMCALGeometry::GetAbsCellId(Int_t nSupMod, Int_t nModule, Int_t nIphi, Int_t nIeta) const
{ 
  // 27-aug-04; 
  // corr. 21-sep-04; 
  //       13-oct-05; 110 degree case
  // May 31, 2006; ALICE numbering scheme:
  // 0 <= nSupMod < fNumberOfSuperModules
  // 0 <= nModule  < fNPHI * fNZ ( fNPHI * fNZ/2 for fKey110DEG=1)
  // 0 <= nIphi   < fNPHIdiv
  // 0 <= nIeta   < fNETAdiv
  // 0 <= absid   < fNCells
  static Int_t id=0; // have to change from 0 to fNCells-1
  if(fKey110DEG == 1 && nSupMod >= 10) { // 110 degree case; last two supermodules
    id  = fNCellsInSupMod*10 + (fNCellsInSupMod/2)*(nSupMod-10);
  } else {
    id  = fNCellsInSupMod*nSupMod;
  }
  id += fNCellsInModule *nModule;
  id += fNPHIdiv *nIphi;
  id += nIeta;
  if(id<0 || id >= fNCells) {
//     printf(" wrong numerations !!\n");
//     printf("    id      %6i(will be force to -1)\n", id);
//     printf("    fNCells %6i\n", fNCells);
//     printf("    nSupMod %6i\n", nSupMod);
//     printf("    nModule  %6i\n", nModule);
//     printf("    nIphi   %6i\n", nIphi);
//     printf("    nIeta   %6i\n", nIeta);
    id = -TMath::Abs(id); // if negative something wrong
  }
  return id;
}

//________________________________________________________________________________________________
Bool_t  AliEMCALGeometry::CheckAbsCellId(Int_t absId) const
{ 
  // May 31, 2006; only trd1 now
  if(absId<0 || absId >= fNCells) return kFALSE;
  else                            return kTRUE;
}

//________________________________________________________________________________________________
Bool_t AliEMCALGeometry::GetCellIndex(Int_t absId,Int_t &nSupMod,Int_t &nModule,Int_t &nIphi,Int_t &nIeta) const
{ 
  // 21-sep-04; 19-oct-05;
  // May 31, 2006; ALICE numbering scheme:
  // 
  // In:
  // absId   - cell is as in Geant,     0<= absId   < fNCells;
  // Out:
  // nSupMod - super module(SM) number, 0<= nSupMod < fNumberOfSuperModules;
  // nModule  - module number in SM,     0<= nModule  < fNCellsInSupMod/fNCellsInSupMod or(/2) for tow last SM (10th and 11th);
  // nIphi   - cell number in phi driection inside module; 0<= nIphi < fNPHIdiv; 
  // nIeta   - cell number in eta driection inside module; 0<= nIeta < fNETAdiv; 
  // 
  static Int_t tmp=0, sm10=0;
  if(!CheckAbsCellId(absId)) return kFALSE;

  sm10 = fNCellsInSupMod*10;
  if(fKey110DEG == 1 && absId >= sm10) { // 110 degree case; last two supermodules  
    nSupMod = (absId-sm10) / (fNCellsInSupMod/2) + 10;
    tmp     = (absId-sm10) % (fNCellsInSupMod/2);
  } else {
    nSupMod = absId / fNCellsInSupMod;
    tmp     = absId % fNCellsInSupMod;
  }

  nModule  = tmp / fNCellsInModule;
  tmp     = tmp % fNCellsInModule;
  nIphi   = tmp / fNPHIdiv;
  nIeta   = tmp % fNPHIdiv;

  return kTRUE;
}

//________________________________________________________________________________________________
void AliEMCALGeometry::GetModulePhiEtaIndexInSModule(Int_t nSupMod, Int_t nModule,  int &iphim, int &ietam) const
{ 
  // added nSupMod; - 19-oct-05 !
  // Alice numbering scheme        - Jun 01,2006 
  // ietam, iphi - indexes of module in two dimensional grid of SM
  // ietam - have to change from 0 to fNZ-1
  // iphim - have to change from 0 to nphi-1 (fNPhi-1 or fNPhi/2-1)
  static Int_t nphi;

  if(fKey110DEG == 1 && nSupMod>=10) nphi = fNPhi/2;
  else                               nphi = fNPhi;

  ietam = nModule/nphi;
  iphim = nModule%nphi;
}

//________________________________________________________________________________________________
void AliEMCALGeometry::GetCellPhiEtaIndexInSModule(Int_t nSupMod, Int_t nModule, Int_t nIphi, Int_t nIeta, 
int &iphi, int &ieta) const
{ 
  // 
  // Added nSupMod; Nov 25, 05
  // Alice numbering scheme  - Jun 01,2006 
  // IN:
  // nSupMod - super module(SM) number, 0<= nSupMod < fNumberOfSuperModules;
  // nModule  - module number in SM,     0<= nModule  < fNCellsInSupMod/fNCellsInSupMod or(/2) for tow last SM (10th and 11th);
  // nIphi   - cell number in phi driection inside module; 0<= nIphi < fNPHIdiv; 
  // nIeta   - cell number in eta driection inside module; 0<= nIeta < fNETAdiv; 
  // 
 // OUT:
  // ieta, iphi - indexes of cell(tower) in two dimensional grid of SM
  // ieta - have to change from 0 to (fNZ*fNETAdiv-1)
  // iphi - have to change from 0 to (fNPhi*fNPHIdiv-1 or fNPhi*fNPHIdiv/2-1)
  //
  static Int_t iphim, ietam;

  GetModulePhiEtaIndexInSModule(nSupMod,nModule, iphim, ietam); 
  //  ieta  = ietam*fNETAdiv + (1-nIeta); // x(module) = -z(SM) 
  ieta  = ietam*fNETAdiv + (fNETAdiv - 1 - nIeta); // x(module) = -z(SM) 
  iphi  = iphim*fNPHIdiv + nIphi;     // y(module) =  y(SM) 

  if(iphi<0 || ieta<0)
  AliDebug(1,Form(" nSupMod %i nModule %i nIphi %i nIeta %i => ieta %i iphi %i\n", 
  nSupMod, nModule, nIphi, nIeta, ieta, iphi));
}

//________________________________________________________________________________________________
Int_t  AliEMCALGeometry::GetSuperModuleNumber(Int_t absId)  const
{
  // Return the number of the  supermodule given the absolute
  // ALICE numbering id

  static Int_t nSupMod, nModule, nIphi, nIeta;
  GetCellIndex(absId, nSupMod, nModule, nIphi, nIeta);
  return nSupMod;
} 

//________________________________________________________________________________________________
void  AliEMCALGeometry::GetModuleIndexesFromCellIndexesInSModule(Int_t nSupMod, Int_t iphi, Int_t ieta, 
			Int_t &iphim, Int_t &ietam, Int_t &nModule) const
{
  // Transition from cell indexes (ieta,iphi) to module indexes (ietam,iphim, nModule)
  static Int_t nphi;
  nphi  = GetNumberOfModuleInPhiDirection(nSupMod);  

  ietam  = ieta/fNETAdiv;
  iphim  = iphi/fNPHIdiv;
  nModule = ietam * nphi + iphim; 
}

//________________________________________________________________________________________________
Int_t  AliEMCALGeometry::GetAbsCellIdFromCellIndexes(Int_t nSupMod, Int_t iphi, Int_t ieta) const
{
  // Transition from super module number(nSupMod) and cell indexes (ieta,iphi) to absId
  static Int_t ietam, iphim, nModule;
  static Int_t nIeta, nIphi; // cell indexes in module

  GetModuleIndexesFromCellIndexesInSModule(nSupMod, iphi, ieta, ietam, iphim, nModule);

  nIeta = ieta%fNETAdiv;
  nIeta = fNETAdiv - 1 - nIeta;
  nIphi = iphi%fNPHIdiv;

  return GetAbsCellId(nSupMod, nModule, nIphi, nIeta);
}


// Methods for AliEMCALRecPoint - Feb 19, 2006
//________________________________________________________________________________________________
Bool_t AliEMCALGeometry::RelPosCellInSModule(Int_t absId, Double_t &xr, Double_t &yr, Double_t &zr) const
{
  // Look to see what the relative
  // position inside a given cell is
  // for a recpoint.
  // Alice numbering scheme - Jun 08, 2006
  // In:
  // absId   - cell is as in Geant,     0<= absId   < fNCells;
  // OUT:
  // xr,yr,zr - x,y,z coordinates of cell with absId inside SM 

  // Shift index taking into account the difference between standard SM 
  // and SM of half size in phi direction
  const Int_t kphiIndexShift = fCentersOfCellsPhiDir.GetSize()/4; // Nov 22, 2006; was 6 for cas 2X2
  static Int_t nSupMod, nModule, nIphi, nIeta, iphi, ieta;
  if(!CheckAbsCellId(absId)) return kFALSE;

  GetCellIndex(absId, nSupMod, nModule, nIphi, nIeta);
  GetCellPhiEtaIndexInSModule(nSupMod,nModule,nIphi,nIeta, iphi, ieta); 
 
  xr = fCentersOfCellsXDir.At(ieta);
  zr = fCentersOfCellsEtaDir.At(ieta);

  if(nSupMod<10) {
    yr = fCentersOfCellsPhiDir.At(iphi);
  } else {
    yr = fCentersOfCellsPhiDir.At(iphi + kphiIndexShift);
  }
  AliDebug(1,Form("absId %i nSupMod %i iphi %i ieta %i xr %f yr %f zr %f ",absId,nSupMod,iphi,ieta,xr,yr,zr));

  return kTRUE;
}

//________________________________________________________________________________________________
Bool_t AliEMCALGeometry::RelPosCellInSModule(Int_t absId, Double_t loc[3]) const
{
  // Alice numbering scheme - Jun 03, 2006
  loc[0] = loc[1] = loc[2]=0.0;
  if(RelPosCellInSModule(absId, loc[0],loc[1],loc[2])) {
    return kTRUE;
  }
  return kFALSE;
}

//________________________________________________________________________________________________
Bool_t AliEMCALGeometry::RelPosCellInSModule(Int_t absId, TVector3 &vloc) const
{
  static Double_t loc[3];
  if(RelPosCellInSModule(absId,loc)) {
    vloc.SetXYZ(loc[0], loc[1], loc[2]);
    return kTRUE;
  } else {
    vloc.SetXYZ(0,0,0);
    return kFALSE;
  }
  // Alice numbering scheme - Jun 03, 2006
}

//________________________________________________________________________________________________
Bool_t AliEMCALGeometry::RelPosCellInSModule(Int_t absId, Double_t distEff, Double_t &xr, Double_t &yr, Double_t &zr) const
{
  // Jul 30, 2007 - taking into account position of shower max
  // Look to see what the relative
  // position inside a given cell is
  // for a recpoint.
  // In:
  // absId   - cell is as in Geant,     0<= absId   < fNCells;
  // e       - cluster energy
  // OUT:
  // xr,yr,zr - x,y,z coordinates of cell with absId inside SM 

  // Shift index taking into account the difference between standard SM 
  // and SM of half size in phi direction
  const  Int_t kphiIndexShift = fCentersOfCellsPhiDir.GetSize()/4; // Nov 22, 2006; was 6 for cas 2X2
  static Int_t nSupMod, nModule, nIphi, nIeta, iphi, ieta;
  static Int_t iphim, ietam;
  static AliEMCALShishKebabTrd1Module *mod = 0;
  static TVector2 v;
  if(!CheckAbsCellId(absId)) return kFALSE;

  GetCellIndex(absId, nSupMod, nModule, nIphi, nIeta);
  GetModulePhiEtaIndexInSModule(nSupMod, nModule, iphim, ietam);
  GetCellPhiEtaIndexInSModule(nSupMod,nModule,nIphi,nIeta, iphi, ieta); 
 
  mod = GetShishKebabModule(ietam);
  mod->GetPositionAtCenterCellLine(nIeta, distEff, v); 
  xr = v.Y() - fParSM[0];
  zr = v.X() - fParSM[2];

  if(nSupMod<10) {
    yr = fCentersOfCellsPhiDir.At(iphi);
  } else {
    yr = fCentersOfCellsPhiDir.At(iphi + kphiIndexShift);
  }
  AliDebug(1,Form("absId %i nSupMod %i iphi %i ieta %i xr %f yr %f zr %f ",absId,nSupMod,iphi,ieta,xr,yr,zr));

  return kTRUE;
}

//________________________________________________________________________________________________
Bool_t AliEMCALGeometry::RelPosCellInSModule(Int_t absId, Int_t maxAbsId, Double_t distEff, Double_t &xr, Double_t &yr, Double_t &zr) const
{
  // Jul 31, 2007 - taking into account position of shower max and apply coor2.
  // Look to see what the relative
  // position inside a given cell is
  // for a recpoint.
  // In:
  // absId     - cell is as in Geant,     0<= absId   < fNCells;
  // maxAbsId  - abs id of cell with highest energy
  // e         - cluster energy
  // OUT:
  // xr,yr,zr - x,y,z coordinates of cell with absId inside SM 

  // Shift index taking into account the difference between standard SM 
  // and SM of half size in phi direction
  const  Int_t kphiIndexShift = fCentersOfCellsPhiDir.GetSize()/4; // Nov 22, 2006; was 6 for cas 2X2
  static Int_t nSupMod, nModule, nIphi, nIeta, iphi, ieta;
  static Int_t iphim, ietam;
  static AliEMCALShishKebabTrd1Module *mod = 0;
  static TVector2 v;

  static Int_t nSupModM, nModuleM, nIphiM, nIetaM, iphiM, ietaM;
  static Int_t iphimM, ietamM, maxAbsIdCopy=-1;
  static AliEMCALShishKebabTrd1Module *modM = 0;
  static Double_t distCorr;

  if(!CheckAbsCellId(absId)) return kFALSE;

  GetCellIndex(absId, nSupMod, nModule, nIphi, nIeta);
  GetModulePhiEtaIndexInSModule(nSupMod, nModule, iphim, ietam);
  GetCellPhiEtaIndexInSModule(nSupMod,nModule,nIphi,nIeta, iphi, ieta); 
  mod = GetShishKebabModule(ietam);

  if(absId != maxAbsId) {
    distCorr = 0.;
    if(maxAbsIdCopy != maxAbsId) {
      GetCellIndex(maxAbsId, nSupModM, nModuleM, nIphiM, nIetaM);
      GetModulePhiEtaIndexInSModule(nSupModM, nModuleM, iphimM, ietamM);
      GetCellPhiEtaIndexInSModule(nSupModM,nModuleM,nIphiM,nIetaM, iphiM, ietaM); 
      modM = GetShishKebabModule(ietamM); // do I need this ?
      maxAbsIdCopy = maxAbsId;
    }
    if(ietamM !=0) {
      distCorr = GetEtaModuleSize()*(ietam-ietamM)/TMath::Tan(modM->GetTheta()); // Stay here
      //printf(" distCorr %f | dist %f | ietam %i -> etamM %i\n", distCorr, dist, ietam, ietamM);  
    }
    // distEff += distCorr;
  }
  // Bad resolution in this case, strong bias vs phi
  // distEff = 0.0; 
  mod->GetPositionAtCenterCellLine(nIeta, distEff, v); // Stay here
  xr = v.Y() - fParSM[0];
  zr = v.X() - fParSM[2];

  if(nSupMod<10) {
    yr = fCentersOfCellsPhiDir.At(iphi);
  } else {
    yr = fCentersOfCellsPhiDir.At(iphi + kphiIndexShift);
  }
  AliDebug(1,Form("absId %i nSupMod %i iphi %i ieta %i xr %f yr %f zr %f ",absId,nSupMod,iphi,ieta,xr,yr,zr));

  return kTRUE;
}

//________________________________________________________________________________________________
void AliEMCALGeometry::CreateListOfTrd1Modules()
{
  // Generate the list of Trd1 modules
  // which will make up the EMCAL
  // geometry

  AliDebug(2,Form(" AliEMCALGeometry::CreateListOfTrd1Modules() started "));

  AliEMCALShishKebabTrd1Module *mod=0, *mTmp=0; // current module
  if(fShishKebabTrd1Modules == 0) {
    fShishKebabTrd1Modules = new TList;
    fShishKebabTrd1Modules->SetName("ListOfTRD1");
    for(int iz=0; iz< GetNZ(); iz++) { 
      if(iz==0) { 
        mod  = new AliEMCALShishKebabTrd1Module(TMath::Pi()/2.,this);
      } else {
        mTmp  = new AliEMCALShishKebabTrd1Module(*mod);
        mod   = mTmp;
      }
      fShishKebabTrd1Modules->Add(mod);
    }
  } else {
    AliDebug(2,Form(" Already exits : "));
  }
  mod = (AliEMCALShishKebabTrd1Module*)fShishKebabTrd1Modules->At(fShishKebabTrd1Modules->GetSize()-1);
  fEtaMaxOfTRD1 = mod->GetMaxEtaOfModule(0);

  AliDebug(2,Form(" fShishKebabTrd1Modules has %i modules : max eta %5.4f \n", 
		  fShishKebabTrd1Modules->GetSize(),fEtaMaxOfTRD1));
  // Feb 20,2006;
  // Jun 01, 2006 - ALICE numbering scheme
  // define grid for cells in eta(z) and x directions in local coordinates system of SM
  // Works just for 2x2 case only -- ?? start here
  // 
  //
  // Define grid for cells in phi(y) direction in local coordinates system of SM
  // as for 2X2 as for 3X3 - Nov 8,2006
  // 
  AliDebug(2,Form(" Cells grid in phi directions : size %i\n", fCentersOfCellsPhiDir.GetSize()));
  Int_t ind=0; // this is phi index
  Int_t ieta=0, nModule=0, iphiTemp;
  Double_t xr=0., zr=0., theta=0., phi=0., eta=0., r=0., x=0.,y=0.;
  TVector3 vglob;
  Double_t ytCenterModule=0.0, ytCenterCell=0.0;

  fCentersOfCellsPhiDir.Set(fNPhi*fNPHIdiv);
  fPhiCentersOfCells.Set(fNPhi*fNPHIdiv);

  Double_t r0 = GetIPDistance() + GetLongModuleSize()/2.;
  for(Int_t it=0; it<fNPhi; it++) { // cycle on modules
    ytCenterModule = -fParSM[1] + fPhiModuleSize*(2*it+1)/2;  // center of module
    for(Int_t ic=0; ic<fNPHIdiv; ic++) { // cycle on cells in module
      if(fNPHIdiv==2) {
        ytCenterCell = ytCenterModule + fPhiTileSize *(2*ic-1)/2.;
      } else if(fNPHIdiv==3){
        ytCenterCell = ytCenterModule + fPhiTileSize *(ic-1);
      } else if(fNPHIdiv==1){
        ytCenterCell = ytCenterModule;
      }
      fCentersOfCellsPhiDir.AddAt(ytCenterCell,ind);
      // Define grid on phi direction
      // Grid is not the same for different eta bin;
      // Effect is small but is still here
      phi = TMath::ATan2(ytCenterCell, r0);
      fPhiCentersOfCells.AddAt(phi, ind);

      AliDebug(2,Form(" ind %2.2i : y %8.3f ", ind, fCentersOfCellsPhiDir.At(ind))); 
      ind++;
    }
  }

  fCentersOfCellsEtaDir.Set(fNZ *fNETAdiv);
  fCentersOfCellsXDir.Set(fNZ *fNETAdiv);
  fEtaCentersOfCells.Set(fNZ *fNETAdiv * fNPhi*fNPHIdiv);
  AliDebug(2,Form(" Cells grid in eta directions : size %i\n", fCentersOfCellsEtaDir.GetSize()));
  for(Int_t it=0; it<fNZ; it++) {
    AliEMCALShishKebabTrd1Module *trd1 = GetShishKebabModule(it);
    nModule = fNPhi*it;
    for(Int_t ic=0; ic<fNETAdiv; ic++) {
      if(fNPHIdiv==2) {
        trd1->GetCenterOfCellInLocalCoordinateofSM(ic, xr, zr);      // case of 2X2
        GetCellPhiEtaIndexInSModule(0, nModule, 0, ic, iphiTemp, ieta); 
      } if(fNPHIdiv==3) {
        trd1->GetCenterOfCellInLocalCoordinateofSM_3X3(ic, xr, zr);  // case of 3X3
        GetCellPhiEtaIndexInSModule(0, nModule, 0, ic, iphiTemp, ieta); 
      } if(fNPHIdiv==1) {
        trd1->GetCenterOfCellInLocalCoordinateofSM_1X1(xr, zr);      // case of 1X1
        GetCellPhiEtaIndexInSModule(0, nModule, 0, ic, iphiTemp, ieta); 
      }
      fCentersOfCellsXDir.AddAt(float(xr) - fParSM[0],ieta);
      fCentersOfCellsEtaDir.AddAt(float(zr) - fParSM[2],ieta);
      // Define grid on eta direction for each bin in phi
      for(int iphi=0; iphi<fCentersOfCellsPhiDir.GetSize(); iphi++) {
        x = xr + trd1->GetRadius();
        y = fCentersOfCellsPhiDir[iphi];
        r = TMath::Sqrt(x*x + y*y + zr*zr);
        theta = TMath::ACos(zr/r);
        eta   = AliEMCALShishKebabTrd1Module::ThetaToEta(theta);
	//        ind   = ieta*fCentersOfCellsPhiDir.GetSize() + iphi;
        ind   = iphi*fCentersOfCellsEtaDir.GetSize() + ieta;
        fEtaCentersOfCells.AddAt(eta, ind);
      }
      //printf(" ieta %i : xr + trd1->GetRadius() %f : zr %f : eta %f \n", ieta, xr + trd1->GetRadius(), zr, eta);
    }
  }
  for(Int_t i=0; i<fCentersOfCellsEtaDir.GetSize(); i++) {
    AliDebug(2,Form(" ind %2.2i : z %8.3f : x %8.3f", i+1, 
                    fCentersOfCellsEtaDir.At(i),fCentersOfCellsXDir.At(i)));
  }

}

//________________________________________________________________________________________________
void AliEMCALGeometry::GetGlobal(const Double_t *loc, Double_t *glob, int ind) const
{
  // Figure out the global numbering
  // of a given supermodule from the
  // local numbering and the transformation
  // matrix stored by the geometry manager (allows for misaligned
  // geometry)

  if(ind>=0 && ind < GetNumberOfSuperModules()) {
    TString volpath = "ALIC_1/XEN1_1/SMOD_";
    volpath += ind+1;

    if(GetKey110DEG() && ind>=10) {
      volpath = "ALIC_1/XEN1_1/SM10_";
      volpath += ind-10+1;
    }

    if(!gGeoManager->cd(volpath.Data()))
      AliFatal(Form("AliEMCALGeometry::GeoManager cannot find path %s!",volpath.Data()));

    TGeoHMatrix* m = gGeoManager->GetCurrentMatrix();
    if(m) {
      m->LocalToMaster(loc, glob);
    } else {
      AliFatal("Geo matrixes are not loaded \n") ;
    }
  }
}

//________________________________________________________________________________________________
void AliEMCALGeometry::GetGlobal(const TVector3 &vloc, TVector3 &vglob, int ind) const
{
  //Figure out the global numbering
  //of a given supermodule from the
  //local numbering given a 3-vector location

  static Double_t tglob[3], tloc[3];
  vloc.GetXYZ(tloc);
  GetGlobal(tloc, tglob, ind);
  vglob.SetXYZ(tglob[0], tglob[1], tglob[2]);
}

//________________________________________________________________________________________________
void AliEMCALGeometry::GetGlobal(Int_t absId , double glob[3]) const
{ 
  // Alice numbering scheme - Jun 03, 2006
  static Int_t nSupMod, nModule, nIphi, nIeta;
  static double loc[3];

  if (!gGeoManager || !gGeoManager->IsClosed()) {
    AliError("Can't get the global coordinates! gGeoManager doesn't exist or it is still open!");
    return;
  }

  glob[0]=glob[1]=glob[2]=0.0; // bad case
  if(RelPosCellInSModule(absId, loc)) {
    GetCellIndex(absId, nSupMod, nModule, nIphi, nIeta);

    TString volpath = "ALIC_1/XEN1_1/SMOD_";
    volpath += (nSupMod+1);

    if(GetKey110DEG() && nSupMod>=10) {
      volpath = "ALIC_1/XEN1_1/SM10_";
      volpath += (nSupMod-10+1);
    }
    if(!gGeoManager->cd(volpath.Data()))
      AliFatal(Form("GeoManager cannot find path %s!",volpath.Data()));

    TGeoHMatrix* m = gGeoManager->GetCurrentMatrix();
    if(m) {
      m->LocalToMaster(loc, glob);
    } else {
      AliFatal("Geo matrixes are not loaded \n") ;
    }
  }
}

//___________________________________________________________________
void AliEMCALGeometry::GetGlobal(Int_t absId , TVector3 &vglob) const
{ 
  // Alice numbering scheme - Jun 03, 2006
  static Double_t glob[3];

  GetGlobal(absId, glob);
  vglob.SetXYZ(glob[0], glob[1], glob[2]);

}

//____________________________________________________________________________
void AliEMCALGeometry::GetGlobal(const AliRecPoint* /*rp*/, TVector3& /* vglob */) const
{
  AliFatal(Form("Please use GetGlobalEMCAL(recPoint,gpos) instead of GetGlobal!"));
}

//_________________________________________________________________________________
void AliEMCALGeometry::GetGlobalEMCAL(const AliEMCALRecPoint *rp, TVector3 &vglob) const
{
  // Figure out the global numbering
  // of a given supermodule from the
  // local numbering for RecPoints

  static TVector3 vloc;
  static Int_t nSupMod, nModule, nIphi, nIeta;

  const AliEMCALRecPoint *rpTmp = rp;
  const AliEMCALRecPoint *rpEmc = rpTmp;

  GetCellIndex(rpEmc->GetAbsId(0), nSupMod, nModule, nIphi, nIeta);
  rpTmp->GetLocalPosition(vloc);
  GetGlobal(vloc, vglob, nSupMod);
}

//________________________________________________________________________________________________
void AliEMCALGeometry::EtaPhiFromIndex(Int_t absId,Double_t &eta,Double_t &phi) const
{
  // Nov 16, 2006- float to double
  // version for TRD1 only
  static TVector3 vglob;
  GetGlobal(absId, vglob);
  eta = vglob.Eta();
  phi = vglob.Phi();
}

//________________________________________________________________________________________________
void AliEMCALGeometry::EtaPhiFromIndex(Int_t absId,Float_t &eta,Float_t &phi) const
{
  // Nov 16,2006 - should be discard in future
  static TVector3 vglob;
  GetGlobal(absId, vglob);
  eta = float(vglob.Eta());
  phi = float(vglob.Phi());
}

//________________________________________________________________________________________________
Bool_t AliEMCALGeometry::GetPhiBoundariesOfSM(Int_t nSupMod, Double_t &phiMin, Double_t &phiMax) const
{
  // 0<= nSupMod <=11; phi in rad
  static int i;
  if(nSupMod<0 || nSupMod >11) return kFALSE; 
  i = nSupMod/2;
  phiMin = fPhiBoundariesOfSM[2*i];
  phiMax = fPhiBoundariesOfSM[2*i+1];
  return kTRUE; 
}

//________________________________________________________________________________________________
Bool_t AliEMCALGeometry::GetPhiBoundariesOfSMGap(Int_t nPhiSec, Double_t &phiMin, Double_t &phiMax) const
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
Bool_t AliEMCALGeometry::SuperModuleNumberFromEtaPhi(Double_t eta, Double_t phi, Int_t &nSupMod) const
{ 
  // Return false if phi belongs a phi cracks between SM
 
  static Int_t i;

  if(TMath::Abs(eta) > fEtaMaxOfTRD1) return kFALSE;

  phi = TVector2::Phi_0_2pi(phi); // move phi to (0,2pi) boundaries
  for(i=0; i<6; i++) {
    if(phi>=fPhiBoundariesOfSM[2*i] && phi<=fPhiBoundariesOfSM[2*i+1]) {
      nSupMod = 2*i;
      if(eta < 0.0) nSupMod++;
      AliDebug(1,Form("eta %f phi %f(%5.2f) : nSupMod %i : #bound %i", eta,phi,phi*TMath::RadToDeg(), nSupMod,i));
      return kTRUE;
    }
  }
  return kFALSE;
}

//________________________________________________________________________________________________
Bool_t AliEMCALGeometry::GetAbsCellIdFromEtaPhi(Double_t eta, Double_t phi, Int_t &absId) const
{
  // Nov 17,2006
  // stay here - phi problem as usual 
  static Int_t nSupMod, i, ieta, iphi, etaShift, nphi;
  static Double_t absEta=0.0, d=0.0, dmin=0.0, phiLoc;
  absId = nSupMod = - 1;
  if(SuperModuleNumberFromEtaPhi(eta, phi, nSupMod)) {
    // phi index first
    phi    = TVector2::Phi_0_2pi(phi);
    phiLoc = phi - fPhiCentersOfSM[nSupMod/2];
    nphi   = fPhiCentersOfCells.GetSize();
    if(nSupMod>=10) {
      phiLoc = phi - 190.*TMath::DegToRad();
      nphi  /= 2;
    }

    dmin   = TMath::Abs(fPhiCentersOfCells[0]-phiLoc);
    iphi   = 0;
    for(i=1; i<nphi; i++) {
      d = TMath::Abs(fPhiCentersOfCells[i] - phiLoc);
      if(d < dmin) {
        dmin = d;
        iphi = i;
      }
      //      printf(" i %i : d %f : dmin %f : fPhiCentersOfCells[i] %f \n", i, d, dmin, fPhiCentersOfCells[i]);
    }
    // odd SM are turned with respect of even SM - reverse indexes
    AliDebug(2,Form(" iphi %i : dmin %f (phi %f, phiLoc %f ) ", iphi, dmin, phi, phiLoc));
    // eta index
    absEta   = TMath::Abs(eta);
    etaShift = iphi*fCentersOfCellsEtaDir.GetSize();
    dmin     = TMath::Abs(fEtaCentersOfCells[etaShift]-absEta);
    ieta     = 0;
    for(i=1; i<fCentersOfCellsEtaDir.GetSize(); i++) {
      d = TMath::Abs(fEtaCentersOfCells[i+etaShift] - absEta);
      if(d < dmin) {
        dmin = d;
        ieta = i;
      }
    }
    AliDebug(2,Form(" ieta %i : dmin %f (eta=%f) : nSupMod %i ", ieta, dmin, eta, nSupMod));

    if(eta<0) iphi = (nphi-1) - iphi;
    absId = GetAbsCellIdFromCellIndexes(nSupMod, iphi, ieta);

    return kTRUE;
  }
  return kFALSE;
}

//________________________________________________________________________________________________
AliEMCALShishKebabTrd1Module* AliEMCALGeometry::GetShishKebabModule(Int_t neta) const
{
  //This method was too long to be
  //included in the header file - the
  //rule checker complained about it's
  //length, so we move it here.  It returns the
  //shishkebabmodule at a given eta index point.

  static AliEMCALShishKebabTrd1Module* trd1=0;
  if(fShishKebabTrd1Modules && neta>=0 && neta<fShishKebabTrd1Modules->GetSize()) {
    trd1 = (AliEMCALShishKebabTrd1Module*)fShishKebabTrd1Modules->At(neta);
  } else trd1 = 0;
  return trd1;
}

//________________________________________________________________________________________________
Int_t AliEMCALGeometry::GetAbsTRUNumberFromNumberInSm(const Int_t row, const Int_t col, const Int_t sm)
{ // Nov 6, 2007
  Int_t itru = row + col*GetNModulesInTRUPhi() + sm*GetNTRU();
  // printf("  GetAbsTRUNumberFromNumberInSm : row %2i col %2i sm %2i -> itru %2i\n", row, col, sm, itru); 
  return itru;
}

//________________________________________________________________________________________________
void AliEMCALGeometry::Browse(TBrowser* b)
{
  //Browse the modules
  if(fShishKebabTrd1Modules) b->Add(fShishKebabTrd1Modules);
}

//________________________________________________________________________________________________
Bool_t AliEMCALGeometry::IsFolder() const
{
  //Check if fShishKebabTrd1Modules is in folder
  if(fShishKebabTrd1Modules) return kTRUE;
  else                       return kFALSE;
}

//________________________________________________________________________________________________
Double_t AliEMCALGeometry::GetPhiCenterOfSM(Int_t nsupmod) const
{
  //returns center of supermodule in phi 
  int i = nsupmod/2;
  return fPhiCentersOfSM[i];

}
//____________________________________________________________________________
Bool_t  AliEMCALGeometry::Impact(const TParticle * particle) const 
{
  // Tells if a particle enters EMCAL
  Bool_t in=kFALSE;
  Int_t AbsID=0;
  TVector3 vtx(particle->Vx(),particle->Vy(),particle->Vz());
  TVector3 vimpact(0,0,0);
  ImpactOnEmcal(vtx,particle->Theta(),particle->Phi(),AbsID,vimpact);
  if(AbsID!=0) 
    in=kTRUE;
  return in;
}
//____________________________________________________________________________
void AliEMCALGeometry::ImpactOnEmcal(TVector3 vtx, Double_t theta, Double_t phi, 
				     Int_t & absId, TVector3 & vimpact) const
{
  // calculates the impact coordinates on EMCAL (centre of a tower/not on EMCAL surface) 
  // of a neutral particle  
  // emitted in the vertex vtx[3] with direction theta and phi in the ALICE global coordinate system

  TVector3 p(TMath::Sin(theta)*TMath::Cos(phi),TMath::Sin(theta)*TMath::Sin(phi),TMath::Cos(theta)) ;

  vimpact.SetXYZ(0,0,0);
  absId=-1;
  if(phi==0 || theta==0) return;

   TVector3 direction;
   Double_t factor = (GetIPDistance()-vtx[1])/p[1];
  direction = vtx + factor*p;

  if (!gGeoManager){
    AliFatal("Geo manager not initialized\n");
  }
  //from particle direction -> tower hitted
  GetAbsCellIdFromEtaPhi(direction.Eta(),direction.Phi(),absId);
  
  //tower absID hitted -> tower/module plane (evaluated at the center of the tower)
  Int_t nSupMod, nModule, nIphi, nIeta;
  Double_t loc[3],loc2[3],loc3[3];
  Double_t glob[3]={},glob2[3]={},glob3[3]={};
  
  if(!RelPosCellInSModule(absId,loc)) return;
  
  //loc is cell center of tower
  GetCellIndex(absId, nSupMod, nModule, nIphi, nIeta);

  //look at 2 neighbours-s cell using nIphi={0,1} and nIeta={0,1}
  Int_t nIphi2,nIeta2,absId2,absId3;
  if(nIeta==0) nIeta2=1;
  else nIeta2=0;
  absId2=GetAbsCellId(nSupMod,nModule,nIphi,nIeta2);  
  if(nIphi==0) nIphi2=1;
  else nIphi2=0;
  absId3=GetAbsCellId(nSupMod,nModule,nIphi2,nIeta);

  //2nd point on emcal cell plane
  if(!RelPosCellInSModule(absId2,loc2)) return;
    
  //3rd point on emcal cell plane
  if(!RelPosCellInSModule(absId3,loc3)) return;
    
  TString volpath = "ALIC_1/XEN1_1/SMOD_";
  volpath += (nSupMod+1);
  
  if(GetKey110DEG() && nSupMod>=10) {
    volpath = "ALIC_1/XEN1_1/SM10_";
    volpath += (nSupMod-10+1);
  }
  if(!gGeoManager->cd(volpath.Data())){
    AliFatal(Form("GeoManager cannot find path %s!",volpath.Data()))
    return;
  }
  TGeoHMatrix* m = gGeoManager->GetCurrentMatrix();
  if(m) {
    m->LocalToMaster(loc, glob);
    m->LocalToMaster(loc2, glob2);
    m->LocalToMaster(loc3, glob3);
  } else {
    AliFatal("Geo matrixes are not loaded \n") ;
  }

  //Equation of Plane from glob,glob2,glob3 (Ax+By+Cz+D=0)
   Double_t A = glob[1]*(glob2[2]-glob3[2]) + glob2[1]*(glob3[2]-glob[2]) + glob3[1]*(glob[2]-glob2[2]);
   Double_t B = glob[2]*(glob2[0]-glob3[0]) + glob2[2]*(glob3[0]-glob[0]) + glob3[2]*(glob[0]-glob2[0]);
   Double_t C = glob[0]*(glob2[1]-glob3[1]) + glob2[0]*(glob3[1]-glob[1]) + glob3[0]*(glob[1]-glob2[1]);
   Double_t D = glob[0]*(glob2[1]*glob3[2]-glob3[1]*glob2[2]) + glob2[0]*(glob3[1]*glob[2]-glob[1]*glob3[2]) + glob3[0]*(glob[1]*glob2[2]-glob2[1]*glob[2]);
  D=-D;
  
  //shift equation of plane from tower/module center to surface along vector (A,B,C) normal to tower/module plane
  Double_t dist = GetLongModuleSize()/2.;
  Double_t norm = TMath::Sqrt(A*A+B*B+C*C);
  Double_t glob4[3]={};
  TVector3 dir(A,B,C);
  TVector3 point(glob[0],glob[1],glob[2]); 
  if(point.Dot(dir)<0) dist*=-1;
  glob4[0]=glob[0]-dist*A/norm;
  glob4[1]=glob[1]-dist*B/norm;
  glob4[2]=glob[2]-dist*C/norm;
  D = glob4[0]*A +  glob4[1]*B +  glob4[2]*C ;
  D = -D;

  //Line determination (2 points for equation of line : vtx and direction)
  //impact between line (particle) and plane (module/tower plane)
  Double_t den = A*(vtx(0)-direction(0)) + B*(vtx(1)-direction(1)) + C*(vtx(2)-direction(2));
  if(den==0){
    printf("ImpactOnEmcal() No solution :\n");
    return;
  }
  
  Double_t length = A*vtx(0)+B*vtx(1)+C*vtx(2)+D;
  length /=den;
  
  vimpact.SetXYZ(vtx(0)+length*(direction(0)-vtx(0)),vtx(1)+length*(direction(1)-vtx(1)),vtx(2)+length*(direction(2)-vtx(2)));
  
  //shift vimpact from tower/module surface to center along vector (A,B,C) normal to tower/module plane
  vimpact.SetXYZ(vimpact(0)+dist*A/norm,vimpact(1)+dist*B/norm,vimpact(2)+dist*C/norm);
  
  return;
}
