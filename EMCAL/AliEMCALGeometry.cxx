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
// Places the the Barrel Geometry of The EMCAL at Midrapidity
// between 80 and 180(or 190) degrees of Phi and
// -0.7 to 0.7 in eta 
// Number of Modules and Layers may be controlled by 
// the name of the instance defined               
//*-- Author: Sahal Yacoob (LBL / UCT)
//     and  : Yves Schutz (SUBATECH)
//     and  : Jennifer Klay (LBL)
//     SHASHLYK : Aleksei Pavlinov (WSU)
//     SuperModules -> module(or tower) -> cell

// --- AliRoot header files ---
#include <TMath.h>
#include <TVector3.h>
#include <TRegexp.h>
#include <TObjArray.h>
#include <TObjString.h>

// -- ALICE Headers.
//#include "AliConst.h"

// --- EMCAL headers
#include "AliEMCALGeometry.h"

ClassImp(AliEMCALGeometry)

AliEMCALGeometry *AliEMCALGeometry::fgGeom = 0;
Bool_t            AliEMCALGeometry::fgInit = kFALSE;
TString name; // contains name of geometry

char *additionalOpts[]={"nl=",   // number of sampling layers
		       "pbTh=", // cm, Thickness of the Pb
		       "scTh="  // cm, Thickness of the Sc
};
int  nAdditionalOpts = sizeof(additionalOpts) / sizeof(char*);

//______________________________________________________________________
AliEMCALGeometry::~AliEMCALGeometry(void){
    // dtor
}

//______________________________________________________________________
Bool_t AliEMCALGeometry::AreInSameTower(Int_t id1, Int_t id2) const {
  // Find out whether two hits are in the same tower - have to be change
  Int_t idmax = TMath::Max(id1, id2) ; 
  Int_t idmin = TMath::Min(id1, id2) ;
  if ( ((idmax - GetNZ() * GetNPhi()) == idmin ) || 
       ((idmax - 2 * GetNZ() * GetNPhi()) == idmin ) )
    return kTRUE ; 
  else 
    return kFALSE ; 
}

//______________________________________________________________________
void AliEMCALGeometry::Init(void){
  // Initializes the EMCAL parameters
  // naming convention : GUV_WX_N_ gives the composition of a tower
  // WX inform about the composition of the EM calorimeter section: 
  //   thickness in mm of Pb radiator (W) and of scintillator (X), and number of scintillator layers (N)
  // New geometry: EMCAL_55_25
  // 24-aug-04 for shish-kebab
  // SHISH_25 or SHISH_62
  // 11-oct-05   - correction for pre final design
  // Feb 06,2006 - decrease the weight of EMCAL
  fgInit = kFALSE; // Assume failed until proven otherwise.
  name   = GetName();
  name.ToUpper();
  fKey110DEG = 0;
  if(name.Contains("110DEG")) fKey110DEG = 1; // for GetAbsCellId

  fNZ             = 114;	// granularity along Z (eta) 
  fNPhi           = 168;	// granularity in phi (azimuth)
  fArm1PhiMin     = 60.0;	// degrees, Starting EMCAL Phi position
  fArm1PhiMax     = 180.0;	// degrees, Ending EMCAL Phi position
  fArm1EtaMin     = -0.7;	// pseudorapidity, Starting EMCAL Eta position
  fArm1EtaMax     = +0.7;	// pseudorapidity, Ending EMCAL Eta position
  fIPDistance     = 454.0;      // cm, Radial distance to inner surface of EMCAL
  fPhiGapForSM    = 0.;         // cm, only for final TRD1 geometry

  // geometry
  if(name.Contains("SHISH")){ // Only shahslyk now
    // 7-sep-05; integration issue
    fArm1PhiMin     = 80.0;	// 60  -> 80
    fArm1PhiMax     = 180.0;	// 180 -> 190

    fNumberOfSuperModules = 10; // 12 = 6 * 2 (6 in phi, 2 in Z);
    fSteelFrontThick = 2.54;    //  9-sep-04
    fIPDistance      = 460.0;
    fFrontSteelStrip = fPassiveScintThick = 0.0; // 13-may-05
    fLateralSteelStrip = 0.025; // before MAY 2005 
    fPhiModuleSize   = fEtaModuleSize   = 11.4;
    fPhiTileSize = fEtaTileSize      = 5.52; // (11.4-5.52*2)/2. = 0.18 cm (wall thickness)
    fNPhi            = 14;
    fNZ              = 30;
    fAlFrontThick    = fGap2Active = 0;
    fNPHIdiv = fNETAdiv = 2;

    fNECLayers       = 62;
    fECScintThick    = fECPbRadThickness = 0.2;
    fSampling        = 1.;  // 30-aug-04 - should be calculated
    if(name.Contains("TWIST")) { // all about EMCAL module
      fNZ             = 27;  // 16-sep-04
    } else if(name.Contains("TRD")) {
      fIPDistance      = 428.0;  //  11-may-05
      fSteelFrontThick = 0.0;    // 3.17 -> 0.0; 28-mar-05 : no stell plate
      fNPhi            = 12;
      fSampling       = 12.327;
      fPhiModuleSize = fEtaModuleSize = 12.26;
      fNZ            = 26;     // 11-oct-04
      fTrd1Angle     = 1.3;    // in degree
// 18-nov-04; 1./0.08112=12.327
// http://pdsfweb01.nersc.gov/~pavlinov/ALICE/SHISHKEBAB/RES/linearityAndResolutionForTRD1.html
      if(name.Contains("TRD1")) {       // 30-jan-05
	// for final design
        fPhiGapForSM    = 2.;         // cm, only for final TRD1 geometry
        if(name.Contains("MAY05") || name.Contains("WSUC") || name.Contains("FINAL")){
          fNumberOfSuperModules = 12; // 20-may-05
          if(name.Contains("WSUC")) fNumberOfSuperModules = 1; // 27-may-05
          fNECLayers     = 77;       // (13-may-05 from V.Petrov)
          fPhiModuleSize = 12.5;     // 20-may-05 - rectangular shape
          fEtaModuleSize = 11.9;
          fECScintThick  = fECPbRadThickness = 0.16;// (13-may-05 from V.Petrov)
          fFrontSteelStrip   = 0.025;// 0.025cm = 0.25mm  (13-may-05 from V.Petrov)
          fLateralSteelStrip = 0.01; // 0.01cm  = 0.1mm   (13-may-05 from V.Petrov) - was 0.025
          fPassiveScintThick = 0.8;  // 0.8cm   = 8mm     (13-may-05 from V.Petrov)
          fNZ                = 24;
          fTrd1Angle         = 1.5;  // 1.3 or 1.5

          if(name.Contains("FINAL")) { // 9-sep-05
            fNumberOfSuperModules = 10;
            if(name.Contains("110DEG")) {
              fNumberOfSuperModules = 12;// last two modules have size 10 degree in phi (180<phi<190)
              fArm1PhiMax = 200.0; // for XEN1 and turn angle of super modules
	    }
            fPhiModuleSize = 12.26 - fPhiGapForSM / Float_t(fNPhi); // first assumption
            fEtaModuleSize = fPhiModuleSize;
            if(name.Contains("HUGE")) fNECLayers *= 3; // 28-oct-05 for analysing leakage    
          }
	}
      } else if(name.Contains("TRD2")) {       // 30-jan-05
        fSteelFrontThick = 0.0;         // 11-mar-05
        fIPDistance+= fSteelFrontThick; // 1-feb-05 - compensate absence of steel plate
        fTrd1Angle  = 1.64;             // 1.3->1.64
        fTrd2AngleY = fTrd1Angle;       //  symmetric case now
        fEmptySpace    = 0.2; // 2 mm
        fTubsR         = fIPDistance; // 31-jan-05 - as for Fred case

        fPhiModuleSize  = fTubsR*2.*TMath::Tan(fTrd2AngleY*TMath::DegToRad()/2.);
        fPhiModuleSize -= fEmptySpace/2.; // 11-mar-05  
        fEtaModuleSize  = fPhiModuleSize; // 20-may-05 
        fTubsTurnAngle  = 3.;
      }
      fNPHIdiv = fNETAdiv  = 2;   // 13-oct-04 - division again
      if(name.Contains("3X3")) {   // 23-nov-04
        fNPHIdiv = fNETAdiv  = 3;
      } else if(name.Contains("4X4")) {
        fNPHIdiv = fNETAdiv  = 4;
      }
    }
    fPhiTileSize = fPhiModuleSize/2. - fLateralSteelStrip; // 13-may-05 
    fEtaTileSize = fEtaModuleSize/2. - fLateralSteelStrip; // 13-may-05 

    if(name.Contains("25")){
      fNECLayers     = 25;
      fECScintThick  = fECPbRadThickness = 0.5;
    }
    if(name.Contains("WSUC")){ // 18-may-05 - about common structure
      fShellThickness = 30.; // should be change 
      fNPhi = fNZ = 4; 
    }

    CheckAditionalOptions();

    // constant for transition absid <--> indexes
    fNCellsInTower  = fNPHIdiv*fNETAdiv;
    fNCellsInSupMod = fNCellsInTower*fNPhi*fNZ;
    fNCells         = fNCellsInSupMod*fNumberOfSuperModules;
    if(name.Contains("110DEG")) fNCells -= fNCellsInSupMod;

    fLongModuleSize = fNECLayers*(fECScintThick + fECPbRadThickness);
    if(name.Contains("MAY05")) fLongModuleSize += (fFrontSteelStrip + fPassiveScintThick);

    // 30-sep-04
    if(name.Contains("TRD")) {
      f2Trd1Dx2 = fEtaModuleSize + 2.*fLongModuleSize*TMath::Tan(fTrd1Angle*TMath::DegToRad()/2.);
      if(name.Contains("TRD2")) {  // 27-jan-05
        f2Trd2Dy2 = fPhiModuleSize + 2.*fLongModuleSize*TMath::Tan(fTrd2AngleY*TMath::DegToRad()/2.);
      }
    }
  } else Fatal("Init", "%s is an undefined geometry!", name.Data()) ; 

  fNPhiSuperModule = fNumberOfSuperModules/2;
  if(fNPhiSuperModule<1) fNPhiSuperModule = 1;
  //There is always one more scintillator than radiator layer because of the first block of aluminium
  fShellThickness = fAlFrontThick + fGap2Active + fNECLayers*GetECScintThick()+(fNECLayers-1)*GetECPbRadThick();
  if(name.Contains("SHISH")) {
    fShellThickness = fSteelFrontThick + fLongModuleSize;
    if(name.Contains("TWIST")) { // 13-sep-04
      fShellThickness  = TMath::Sqrt(fLongModuleSize*fLongModuleSize + fPhiModuleSize*fEtaModuleSize);
      fShellThickness += fSteelFrontThick;
    } else if(name.Contains("TRD")) { // 1-oct-04
      fShellThickness  = TMath::Sqrt(fLongModuleSize*fLongModuleSize + f2Trd1Dx2*f2Trd1Dx2);
      fShellThickness += fSteelFrontThick;
    }
  }

  fZLength        = 2.*ZFromEtaR(fIPDistance+fShellThickness,fArm1EtaMax); // Z coverage
  fEnvelop[0]     = fIPDistance; // mother volume inner radius
  fEnvelop[1]     = fIPDistance + fShellThickness; // mother volume outer r.
  fEnvelop[2]     = 1.00001*fZLength; // add some padding for mother volume. 
  
  fgInit = kTRUE; 
  
  if (kTRUE) {
    printf("Init: geometry of EMCAL named %s is as follows:\n", name.Data());
    printf( "               ECAL      : %d x (%f cm Pb, %f cm Sc) \n", GetNECLayers(), GetECPbRadThick(), GetECScintThick() ) ; 
    if(name.Contains("SHISH")){
      printf(" fIPDistance       %6.3f cm \n", fIPDistance);
      if(fSteelFrontThick>0.) 
      printf(" fSteelFrontThick  %6.3f cm \n", fSteelFrontThick);
      printf(" fNPhi %i   |  fNZ %i \n", fNPhi, fNZ);
      printf(" fNCellsInTower %i : fNCellsInSupMod %i : fNCells %i\n",fNCellsInTower, fNCellsInSupMod, fNCells);
      if(name.Contains("MAY05")){
	printf(" fFrontSteelStrip         %6.4f cm (thickness of front steel strip)\n", 
        fFrontSteelStrip);
	printf(" fLateralSteelStrip       %6.4f cm (thickness of lateral steel strip)\n", 
        fLateralSteelStrip);
 	printf(" fPassiveScintThick  %6.4f cm (thickness of front passive Sc tile)\n",
        fPassiveScintThick);
      }
      printf(" X:Y module size     %6.3f , %6.3f cm \n", fPhiModuleSize, fEtaModuleSize);
      printf(" X:Y   tile size     %6.3f , %6.3f cm \n", fPhiTileSize, fEtaTileSize);
      printf(" #of sampling layers %i(fNECLayers) \n", fNECLayers);
      printf(" fLongModuleSize     %6.3f cm \n", fLongModuleSize);
      printf(" #supermodule in phi direction %i \n", fNPhiSuperModule );
    }
    if(name.Contains("TRD")) {
      printf(" fTrd1Angle %7.4f\n", fTrd1Angle);
      printf(" f2Trd1Dx2  %7.4f\n",  f2Trd1Dx2);
      if(name.Contains("TRD2")) {
        printf(" fTrd2AngleY     %7.4f\n", fTrd2AngleY);
        printf(" f2Trd2Dy2       %7.4f\n", f2Trd2Dy2);
        printf(" fTubsR          %7.2f cm\n", fTubsR);
        printf(" fTubsTurnAngle  %7.4f\n", fTubsTurnAngle);
        printf(" fEmptySpace     %7.4f cm\n", fEmptySpace);
      } else if(name.Contains("TRD1") && name.Contains("FINAL")){
        printf(" fPhiGapForSM  %7.4f cm \n",  fPhiGapForSM);
        if(name.Contains("110DEG"))printf(" Last two modules have size 10 degree in  phi (180<phi<190)\n");
      }
    }
    printf("Granularity: %d in eta and %d in phi\n", GetNZ(), GetNPhi()) ;
    printf("Layout: phi = (%7.1f, %7.1f), eta = (%5.2f, %5.2f), IP = %7.2f\n",  
	   GetArm1PhiMin(), GetArm1PhiMax(),GetArm1EtaMin(), GetArm1EtaMax(), GetIPDistance() );
  }
}

//______________________________________________________________________

void AliEMCALGeometry::CheckAditionalOptions()
{ // Feb 06,2006
  fArrayOpts = new TObjArray;
  Int_t nopt = ParseString(name, *fArrayOpts);
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
    for(Int_t j=0; j<nAdditionalOpts; j++) {
      TString opt = additionalOpts[j];
      if(addOpt.Contains(opt,TString::kIgnoreCase)) {
	  indj = j;
        break;
      }
    }
    if(indj<0) {
      printf("<E> option |%s| unavailable : ** look to the file AliEMCALGeometry.h **\n", 
      addOpt.Data());
      assert(0);
    } else {
      printf("<I> option |%s| is valid : number %i : |%s|\n", 
	     addOpt.Data(), indj, additionalOpts[indj]);
      if       (addOpt.Contains("NL=",TString::kIgnoreCase))   {// number of sampling layers
        sscanf(addOpt.Data(),"NL=%i", &fNECLayers);
        printf(" fNECLayers %i (new) \n", fNECLayers);
      } else if(addOpt.Contains("PBTH=",TString::kIgnoreCase)) {//Thickness of the Pb
        sscanf(addOpt.Data(),"PBTH=%f", &fECPbRadThickness);
      } else if(addOpt.Contains("SCTH=",TString::kIgnoreCase)) {//Thickness of the Sc
        sscanf(addOpt.Data(),"SCTH=%f", &fECScintThick);
      }
    }
  }
}

//______________________________________________________________________
AliEMCALGeometry *  AliEMCALGeometry::GetInstance(){ 
  // Returns the pointer of the unique instance
  
  return static_cast<AliEMCALGeometry *>( fgGeom ) ; 
}

//______________________________________________________________________
AliEMCALGeometry* AliEMCALGeometry::GetInstance(const Text_t* name,
						const Text_t* title){
    // Returns the pointer of the unique instance

    AliEMCALGeometry * rv = 0; 
    if ( fgGeom == 0 ) {
	if ( strcmp(name,"") == 0 ) rv = 0;
	else {    
	    fgGeom = new AliEMCALGeometry(name, title);
	    if ( fgInit ) rv = (AliEMCALGeometry * ) fgGeom;
	    else {
		rv = 0; 
		delete fgGeom; 
		fgGeom = 0; 
	    } // end if fgInit
	} // end if strcmp(name,"")
    }else{
	if ( strcmp(fgGeom->GetName(), name) != 0 ) {
	  printf("\ncurrent geometry is ") ;  
	  printf(fgGeom->GetName());
	  printf("\n                      you cannot call     "); 
	  printf(name);  
	}else{
	  rv = (AliEMCALGeometry *) fgGeom; 
	} // end if
    }  // end if fgGeom
    return rv; 
}

// These methods are obsolete but use in AliEMCALRecPoint - keep it now
//______________________________________________________________________
Int_t AliEMCALGeometry::TowerIndex(Int_t ieta,Int_t iphi) const {
  // Returns the tower index number from the based on the Z and Phi
  // index numbers.
  // Inputs:
  //   Int_t ieta    // index along z axis [1-fNZ]
  //   Int_t iphi  // index along phi axis [1-fNPhi]
  // Outputs:
  //   none.
  // Returned
  //   Int_t index // Tower index number 
  
  if ( (ieta <= 0 || ieta>GetNEta()) || 
       (iphi <= 0 || iphi>GetNPhi())) {
    Error("TowerIndex", "Unexpected parameters eta = %d phi = %d!", ieta, iphi) ; 
    return -1;
  }
  return ( (iphi - 1)*GetNEta() + ieta ); 
}

//______________________________________________________________________
void AliEMCALGeometry::TowerIndexes(Int_t index,Int_t &ieta,Int_t &iphi) const {
  // Inputs:
  //   Int_t index // Tower index number [1-fNZ*fNPhi]
  // Outputs:
  //   Int_t ieta    // index allong z axis [1-fNZ]
  //   Int_t iphi  // index allong phi axis [1-fNPhi]
  // Returned
  //   none.

  Int_t nindex = 0;

  if ( IsInECA(index) ) { // ECAL index
    nindex = index ;
  }
  else {
    Error("TowerIndexes", "Unexpected Id number!") ;
    ieta = -1;
    iphi = -1;
    return;
  }   

  if (nindex%GetNZ()) 
    iphi = nindex / GetNZ() + 1 ; 
  else 
    iphi = nindex / GetNZ() ; 
  ieta = nindex - (iphi - 1) * GetNZ() ; 

  if (gDebug==2)
    printf("TowerIndexes: index=%d,%d, ieta=%d, iphi = %d", index, nindex,ieta, iphi) ; 
  return;
  
}

//______________________________________________________________________
void AliEMCALGeometry::EtaPhiFromIndex(Int_t index,Float_t &eta,Float_t &phi) const {
    // given the tower index number it returns the based on the eta and phi
    // of the tower.
    // Inputs:
    //   Int_t index // Tower index number [1-fNZ*fNPhi]
    // Outputs:
    //   Float_t eta  // eta of center of tower in pseudorapidity
    //   Float_t phi  // phi of center of tower in degrees
    // Returned
    //   none.
    Int_t ieta, iphi;
    Float_t deta, dphi ;

    TowerIndexes(index,ieta,iphi);
    
    if (gDebug == 2) 
      printf("EtaPhiFromIndex: index = %d, ieta = %d, iphi = %d", index, ieta, iphi) ;

    deta = (GetArm1EtaMax()-GetArm1EtaMin())/(static_cast<Float_t>(GetNEta()));
    eta  = GetArm1EtaMin() + ((static_cast<Float_t>(ieta) - 0.5 ))*deta;

    dphi = (GetArm1PhiMax() - GetArm1PhiMin())/(static_cast<Float_t>(GetNPhi()));  // in degrees.
    phi  = GetArm1PhiMin() + dphi*(static_cast<Float_t>(iphi) - 0.5);//iphi range [1-fNphi].
}

//______________________________________________________________________
Int_t AliEMCALGeometry::TowerIndexFromEtaPhi(Float_t eta,Float_t phi) const {
    // returns the tower index number based on the eta and phi of the tower.
    // Inputs:
    //   Float_t eta  // eta of center of tower in pseudorapidity
    //   Float_t phi  // phi of center of tower in degrees
    // Outputs:
    //   none.
    // Returned
    //   Int_t index // Tower index number [1-fNZ*fNPhi]

    Int_t ieta,iphi;

    ieta = static_cast<Int_t> ( 1 + (static_cast<Float_t>(GetNEta()) * (eta - GetArm1EtaMin()) / (GetArm1EtaMax() - GetArm1EtaMin())) ) ;

    if( ieta <= 0 || ieta > GetNEta() ) { 
      Error("TowerIndexFromEtaPhi", "Unexpected (eta, phi) = (%f, %f) value, outside of EMCAL!", eta, phi) ; 
      return -1 ; 
    }

    iphi = static_cast<Int_t> ( 1 + (static_cast<Float_t>(GetNPhi()) * (phi - GetArm1PhiMin()) / (GetArm1PhiMax() - GetArm1PhiMin())) ) ;

    if( iphi <= 0 || iphi > GetNPhi() ) { 
      Error("TowerIndexFromEtaPhi", "Unexpected (eta, phi) = (%f, %f) value, outside of EMCAL!", eta, phi) ; 
      return -1 ; 
    }

    return TowerIndex(ieta,iphi);
}

//______________________________________________________________________
Bool_t AliEMCALGeometry::AbsToRelNumbering(Int_t AbsId, Int_t *relid) const {
    // Converts the absolute numbering into the following array/
    //  relid[0] = Row number inside EMCAL
    //  relid[1] = Column number inside EMCAL
    // Input:
    //   Int_t AbsId // Tower index number [1-2*fNZ*fNPhi]
    // Outputs:
    //   Int_t *relid // array of 2. Described above.
    Bool_t rv  = kTRUE ;
    Int_t ieta=0,iphi=0,index=AbsId;

    TowerIndexes(index,ieta,iphi);
    relid[0] = ieta;
    relid[1] = iphi;

    return rv;
}

//______________________________________________________________________
void AliEMCALGeometry::PosInAlice(const Int_t *relid, Float_t &theta, Float_t &phi) const 
{
  // Converts the relative numbering into the local EMCAL-module (x, z)
  // coordinates
  Int_t ieta = relid[0]; // offset along x axis
  Int_t iphi = relid[1]; // offset along z axis
  Int_t index;
  Float_t eta;
  
  index = TowerIndex(ieta,iphi);
  EtaPhiFromIndex(index,eta,phi);
  //theta = 180.*(2.0*TMath::ATan(TMath::Exp(-eta)))/TMath::Pi();
  theta = 2.0*TMath::ATan(TMath::Exp(-eta));

  // correct for distance to IP
  Float_t d = GetIP2ECASection() - GetIPDistance() ;  

  Float_t correction = 1 + d/GetIPDistance() ; 
  Float_t tantheta = TMath::Tan(theta) * correction ; 
  theta = TMath::ATan(tantheta) * TMath::RadToDeg() ; 
  if (theta < 0 ) 
    theta += 180. ; 
  
  return;
}

//______________________________________________________________________
void AliEMCALGeometry::PosInAlice(Int_t absid, Float_t &theta, Float_t &phi) const 
{
  // Converts the relative numbering into the local EMCAL-module (x, z)
  // coordinates
  Int_t relid[2] ; 
  AbsToRelNumbering(absid, relid) ;
  Int_t ieta = relid[0]; // offset along x axis
  Int_t iphi = relid[1]; // offset along z axis
  Int_t index;
  Float_t eta;
  
  index = TowerIndex(ieta,iphi);
  EtaPhiFromIndex(index,eta,phi);
  theta = 2.0*TMath::ATan(TMath::Exp(-eta)) ;
  
  // correct for distance to IP
  Float_t d = 0. ; 
  if (IsInECA(absid))
    d = GetIP2ECASection() - GetIPDistance() ; 
  else {
    Error("PosInAlice", "Unexpected id # %d!", absid) ; 
    return;
  }

  Float_t correction = 1 + d/GetIPDistance() ; 
  Float_t tantheta = TMath::Tan(theta) * correction ; 
  theta = TMath::ATan(tantheta) * TMath::RadToDeg() ; 
  if (theta < 0 ) 
    theta += 180. ; 
  
  return;
}

//______________________________________________________________________
void AliEMCALGeometry::XYZFromIndex(const Int_t *relid,Float_t &x,Float_t &y, Float_t &z) const {
    // given the tower relative number it returns the X, Y and Z
    // of the tower.
    
    // Outputs:
    //   Float_t x  // x of center of tower in cm
    //   Float_t y  // y of center of tower in cm
    //   Float_t z  // z of centre of tower in cm
    // Returned
    //   none.
    
    Float_t eta,theta, phi,cylradius=0. ;
    
    Int_t ieta = relid[0]; // offset along x axis
    Int_t iphi = relid[1]; // offset along z axis.
    Int_t index;
    
    index = TowerIndex(ieta,iphi);
    EtaPhiFromIndex(index,eta,phi);
    theta = 180.*(2.0*TMath::ATan(TMath::Exp(-eta)))/TMath::Pi();
    
    cylradius = GetIP2ECASection() ;  

    Double_t  kDeg2Rad = TMath::DegToRad() ; 
    x =  cylradius * TMath::Cos(phi * kDeg2Rad ) ;
    y =  cylradius * TMath::Sin(phi * kDeg2Rad ) ; 
    z =  cylradius / TMath::Tan(theta * kDeg2Rad ) ; 
 
 return;
} 

//______________________________________________________________________
void AliEMCALGeometry::XYZFromIndex(Int_t absid,  TVector3 &v) const {
    // given the tower relative number it returns the X, Y and Z
    // of the tower.
    
    // Outputs:
    //   Float_t x  // x of center of tower in cm
    //   Float_t y  // y of center of tower in cm
    //   Float_t z  // z of centre of tower in cm
    // Returned
    //   none.
    
    Float_t theta, phi,cylradius=0. ;
        
    PosInAlice(absid, theta, phi) ; 
    
    if ( IsInECA(absid) ) 
      cylradius = GetIP2ECASection() ;
    else {
      Error("XYZFromIndex", "Unexpected Tower section") ;
      return;
    }

    Double_t  kDeg2Rad = TMath::DegToRad() ; 
    v.SetX(cylradius * TMath::Cos(phi * kDeg2Rad ) );
    v.SetY(cylradius * TMath::Sin(phi * kDeg2Rad ) ); 
    v.SetZ(cylradius / TMath::Tan(theta * kDeg2Rad ) ) ; 
 
 return;
} 

Bool_t AliEMCALGeometry::IsInEMCAL(Double_t x, Double_t y, Double_t z) const {
  // Checks whether point is inside the EMCal volume
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
     if (phi > fArm1PhiMin && phi < fArm1PhiMax)
       return 1;
  }
  return 0;
}
// ==

//
// == Shish-kebab cases ==
//
Int_t AliEMCALGeometry::GetAbsCellId(Int_t nSupMod, Int_t nTower, Int_t nIphi, Int_t nIeta)
{ // 27-aug-04; 
  // corr. 21-sep-04; 
  //       13-oct-05; 110 degree case
  // 1 <= nSupMod <= fNumberOfSuperModules
  // 1 <= nTower  <= fNPHI * fNZ ( fNPHI * fNZ/2 for fKey110DEG=1)
  // 1 <= nIphi   <= fNPHIdiv
  // 1 <= nIeta   <= fNETAdiv
  // 1 <= absid   <= fNCells
  static Int_t id=0; // have to change from 1 to fNCells
  if(fKey110DEG == 1 && nSupMod > 10) { // 110 degree case; last two supermodules
    id  = fNCellsInSupMod*10 + (fNCellsInSupMod/2)*(nSupMod-11);
  } else {
    id  = fNCellsInSupMod*(nSupMod-1);
  }
  id += fNCellsInTower *(nTower-1);
  id += fNPHIdiv *(nIphi-1);
  id += nIeta;
  if(id<=0 || id > fNCells) {
//     printf(" wrong numerations !!\n");
//     printf("    id      %6i(will be force to -1)\n", id);
//     printf("    fNCells %6i\n", fNCells);
//     printf("    nSupMod %6i\n", nSupMod);
//     printf("    nTower  %6i\n", nTower);
//     printf("    nIphi   %6i\n", nIphi);
//     printf("    nIeta   %6i\n", nIeta);
    id = -TMath::Abs(id);
  }
  return id;
}

Bool_t  AliEMCALGeometry::CheckAbsCellId(Int_t ind)
{ // 17-niv-04 - analog of IsInECA
   if(name.Contains("TRD")) {
     if(ind<=0 || ind > fNCells) return kFALSE;
     else                        return kTRUE;
   } else return IsInECA(ind);
}

Bool_t AliEMCALGeometry::GetCellIndex(Int_t absId,Int_t &nSupMod,Int_t &nTower,Int_t &nIphi,Int_t &nIeta)
{ // 21-sep-04
  // 19-oct-05;
  static Int_t tmp=0, sm10=0;
  if(absId<=0 || absId>fNCells) {
//     Info("GetCellIndex"," wrong abs Id %i !! \n", absId); 
    return kFALSE;
  }
  sm10 = fNCellsInSupMod*10;
  if(fKey110DEG == 1 && absId > sm10) { // 110 degree case; last two supermodules  
    nSupMod = (absId-1-sm10) / (fNCellsInSupMod/2) + 11;
    tmp     = (absId-1-sm10) % (fNCellsInSupMod/2);
  } else {
    nSupMod = (absId-1) / fNCellsInSupMod + 1;
    tmp     = (absId-1) % fNCellsInSupMod;
  }

  nTower  = tmp / fNCellsInTower + 1;
  tmp     = tmp % fNCellsInTower;
  nIphi   = tmp / fNPHIdiv + 1;
  nIeta   = tmp % fNPHIdiv + 1;

  return kTRUE;
}

void AliEMCALGeometry::GetTowerPhiEtaIndexInSModule(Int_t nSupMod, Int_t nTower,  int &iphit, int &ietat)
{ // added nSupMod; have to check  - 19-oct-05 ! 
  static Int_t nphi;

  if(fKey110DEG == 1 && nSupMod>=11) nphi = fNPhi/2;
  else                               nphi = fNPhi;

  ietat = (nTower-1)/nphi + 1; // have to change from 1 to fNZ
  iphit = (nTower-1)%nphi + 1; // have to change from 1 to fNPhi
}

void AliEMCALGeometry::GetCellPhiEtaIndexInSModule(Int_t nSupMod, Int_t nTower, Int_t nIphi, Int_t nIeta, 
int &iphi, int &ieta)
{ // added nSupMod; Nov 25, 05
  static Int_t iphit, ietat;

  GetTowerPhiEtaIndexInSModule(nSupMod,nTower, iphit, ietat); 
  // have to change from 1 to fNZ*fNETAdiv
  ieta  = (ietat-1)*fNETAdiv + (3-nIeta); // x(module) = -z(SM) 
  // iphi - have to change from 1 to fNPhi*fNPHIdiv
  iphi  = (iphit-1)*fNPHIdiv + nIphi;     // y(module) =  y(SM) 
}
// Service routine 
int  AliEMCALGeometry::ParseString(const TString &topt, TObjArray &Opt)
{ // Feb 06, 2006
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
