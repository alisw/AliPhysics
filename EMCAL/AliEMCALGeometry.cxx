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
// between 0 and 120 degrees of Phi and
// -0.7 to 0.7 in eta 
// Number of Modules and Layers may be controlled by 
// the name of the instance defined               
// EMCALArch2x has more modules along both phi and eta
// EMCALArchxa has less Layers in the Radial Direction
//*-- Author: Sahal Yacoob (LBL / UCT)
//     and  : Yves Schutz (SUBATECH)
//     and  : Jennifer Klay (LBL)

// --- ROOT system ---

// --- Standard library ---
//#include <stdlib.h> 

// --- AliRoot header files ---
//#include <TError.h>
#include <TMath.h>
#include <TVector3.h>

// -- ALICE Headers.
//#include "AliConst.h"

// --- EMCAL headers
#include "AliEMCALGeometry.h"

ClassImp(AliEMCALGeometry);

AliEMCALGeometry *AliEMCALGeometry::fgGeom = 0;
Bool_t            AliEMCALGeometry::fgInit = kFALSE;

//______________________________________________________________________
AliEMCALGeometry::~AliEMCALGeometry(void){
    // dtor
}

//______________________________________________________________________
const Bool_t AliEMCALGeometry::AreInSameTower(Int_t id1, Int_t id2) const {
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
  // naming convention : GUV_L_WX_N_YZ_M gives the composition of a tower
  // UV inform about the compsition of the pre-shower section: 
  //   thickness in mm of Pb radiator (U) and of scintillator (V), and number of scintillator layers (L)
  // WX inform about the composition of the EM calorimeter section: 
  //   thickness in mm of Pb radiator (W) and of scintillator (X), and number of scintillator layers (N) 
  // YZ inform about the composition of the hadron calorimeter section: 
  //   thickness in mm of Cu radiator (Y) and of scintillator (Z), and number of scintillator layers (M) 
  // Valid geometries are G56_2_55_19_104_14
  //                      G56_2_55_19 or EMCAL_5655_21
  //                      G65_2_64_19 or EMCAL_6564_21 

  fgInit = kFALSE; // Assume failer untill proven otherwise.
  TString name(GetName()) ; 

  if ( name == "G56_2_55_19_104_14" ) {
    fPRPbRadThickness  = 0.5;  // cm, Thickness of the Pb radiators for the preshower section 
    fPRScintThick      = 0.6;  // cm, Thickness of the sintilator for the preshower section of the tower
    fNPRLayers         = 2;    // number of scintillator layers in the preshower section 
    
    fECPbRadThickness  = 0.5;  // cm, Thickness of the Pb radiators for the EM calorimeter  section 
    fECScintThick      = 0.5;  // cm, Thickness of the sintilator for the EM alorimeter section of the tower  
    fNECLayers         = 19;   // number of scintillator layers in the EM calorimeter section 
    
    fHCCuRadThickness  = 1.0;  // cm, Thickness of the Cu radiators.
    fHCScintThick      = 0.4;  // cm, Thickness of the sintilator for the hadronic alorimeter section of the tower  
    fNHCLayers         = 14;   // number of scintillator layers in the hadronic calorimeter section
    
    fSampling          = 11.3 ; 
    fSummationFraction = 0.8 ;

    fAlFrontThick      = 3.0;  // cm, Thickness of front Al layer
    fGap2Active        = 1.0;  // cm, Gap between Al and 1st Scintillator
  }
  else if ( name == "G56_2_55_19" || name == "EMCAL_5655_21" ) {
    fPRPbRadThickness  = 0.5;  // cm, Thickness of the Pb radiators for the preshower section 
    fPRScintThick      = 0.6;  // cm, Thickness of the sintilator for the preshower section of the tower
    fNPRLayers         = 2;    // number of scintillator layers in the preshower section 
    
    fECPbRadThickness  = 0.5;  // cm, Thickness of the Pb radiators for the EM calorimeter  section 
    fECScintThick      = 0.5;  // cm, Thickness of the sintilator for the EM alorimeter section of the tower  
    fNECLayers         = 19;   // number of scintillator layers in the EM calorimeter section 
    
    fHCCuRadThickness  = 0.0;  // cm, Thickness of the Cu radiators.
    fHCScintThick      = 0.0;  // cm, Thickness of the sintilator for the hadronic alorimeter section of the tower  
    fNHCLayers         = 0;    // number of scintillator layers in the hadronic calorimeter section
    
    fSampling          = 11.3 ; 
    fSummationFraction = 0.8 ;
 
    fAlFrontThick      = 3.0;  // cm, Thickness of front Al layer
    fGap2Active        = 1.0;  // cm, Gap between Al and 1st Scintillator
  }
  else if ( name == "G65_2_64_19" || name == "EMCAL_6564_21" ) {
    fPRPbRadThickness  = 0.6;  // cm, Thickness of the Pb radiators for the preshower section 
    fPRScintThick      = 0.5;  // cm, Thickness of the sintilator for the preshower section of the tower
    fNPRLayers         = 2;    // number of scintillator layers in the preshower section 
    
    fECPbRadThickness  = 0.6;  // cm, Thickness of the Pb radiators for the EM calorimeter  section 
    fECScintThick      = 0.4;  // cm, Thickness of the sintilator for the EM alorimeter section of the tower  
    fNECLayers         = 19;   // number of scintillator layers in the EM calorimeter section 
    
    fHCCuRadThickness  = 0.0;  // cm, Thickness of the Cu radiators.
    fHCScintThick      = 0.0;  // cm, Thickness of the sintilator for the hadronic alorimeter section of the tower  
    fNHCLayers         = 0;    // number of scintillator layers in the hadronic calorimeter section
    
    fSampling          = 16. ; 
    fSummationFraction = 0.8 ;
 
    fAlFrontThick      = 3.0;  // cm, Thickness of front Al layer
    fGap2Active        = 1.0;  // cm, Gap between Al and 1st Scintillator
  }
  else
    Fatal("Init", "%s is an undefined geometry!", name.Data()) ; 
		 
   //  if( name != "EMCALArch1a" &&
// 	name != "EMCALArch1b" && 
// 	name != "EMCALArch2a" && 
// 	name != "EMCALArch2b" && 
// 	name != "EMCALArch1aN" ){
//       Fatal("Init", "%s is not a known geometry (choose among EMCALArch1a, EMCALArch1b, EMCALArch2a and EMCALArch2b, EMCALArch1aN)",  name.Data()) ;  
//     } // end if
//     //
//     if ( name == "EMCALArch1a"  ||
// 	 name == "EMCALArch1b"  || 
// 	 name == "EMCALArch1aN") {
//       fNZ         = 96;
//       fNPhi       = 144;
//     } // end if
//     if ( name == "EMCALArch2a"  ||
// 	 name == "EMCALArch2b" ) {
// 	fNZ         = 112;
// 	fNPhi       = 168;
//     } // end if
//     if ( name == "EMCALArch1a"  ||
// 	 name == "EMCALArch2a" ) {
//       fNPRLayers  = 2;
//       fNECLayers  = 19;
//       fNHCLayers  = 0;
//     } // end if
//     if ( name == "EMCALArch1b"  ||
// 	 name == "EMCALArch2b" ) {
// 	fNPRLayers  = 2;
// 	fNECLayers  = 23;
// 	fNHCLayers  = 0;
//     } // end if
//     if ( name == "EMCALArch1aN") { 
//       fNPRLayers   = 2;
//       fNECLayers  = 19;
//       fNHCLayers  = 14;
//     }

  // geometry
  fNZ             = 96;    // granularity along Z (eta) 
  fNPhi           = 144;   // granularity in phi (azimuth)
  fArm1PhiMin     =  60.0; // degrees, Starting EMCAL Phi position
  fArm1PhiMax     = 180.0; // degrees, Ending EMCAL Phi position
  fArm1EtaMin     = -0.7;  // pseudorapidity, Starting EMCAL Eta position
  fArm1EtaMax     = +0.7;  // pseudorapidity, Ending EMCAL Eta position
  
  fIPDistance     = 454.0; // cm, Radial distance to inner surface of EMCAL
  fShellThickness = fAlFrontThick + fGap2Active + 2.*(GetPRScintThick() + GetPRPbRadThick()) + // pre shower 
    (fNECLayers-1)*(GetECScintThick()+ GetECPbRadThick()) + // E cal -1 because the last element is a scintillator
    fNHCLayers*(GetHCScintThick()+ GetHCCuRadThick()) + // H cal
    GetHCScintThick() ; // last scintillator
  fZLength        = 2.*ZFromEtaR(fIPDistance+fShellThickness,fArm1EtaMax); // Z coverage
  fEnvelop[0]     = fIPDistance; // mother volume inner radius
  fEnvelop[1]     = fIPDistance + fShellThickness; // mother volume outer r.
  fEnvelop[2]     = 1.00001*fZLength; // add some padding for mother volume. 
  
  fgInit = kTRUE; 
  
  if (gDebug) {
    Info("Init", "geometry of EMCAL named %s is as follows:", name.Data());
    printf( "Tower geometry pre-shower: %d x (%f mm Pb, %f mm Sc) \n", GetNPRLayers(), GetPRPbRadThick(), GetPRScintThick() ) ; 
    printf( "               ECAL      : %d x (%f mm Pb, %f mm Sc) \n", GetNECLayers(), GetECPbRadThick(), GetECScintThick() ) ; 
    if ( GetNHCLayers() > 0 )
      printf( "               HCAL      : %d x (%f mm Pb, %f mm Sc) \n", GetNHCLayers(), GetHCCuRadThick(), GetHCScintThick() ) ; 
    printf("Granularity: %d in eta and %d in phi\n", GetNZ(), GetNPhi()) ;
    printf("Layout: phi = (%f, %f), eta = (%f, %f), y = %f\n",  
	   GetArm1PhiMin(), GetArm1PhiMax(),GetArm1EtaMin(), GetArm1EtaMax(), GetIPDistance() ) ;    
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
	  TString message("\n") ; 
	  message += "current geometry is " ;  
	  message += fgGeom->GetName() ;
	  message += "\n                      you cannot call     " ; 
	  message += name ;  
	  ::Info("GetGeometry", message.Data() ) ; 
	}else{
	  rv = (AliEMCALGeometry *) fgGeom; 
	} // end if
    }  // end if fgGeom
    return rv; 
}

//______________________________________________________________________
Int_t AliEMCALGeometry::TowerIndex(Int_t ieta,Int_t iphi) const {
  // Returns the tower index number from the based on the Z and Phi
  // index numbers. There are 2 times the number of towers to separate
  // out the full towers from the pre-showers.
  // Inputs:
  //   Int_t ieta    // index allong z axis [1-fNZ]
  //   Int_t iphi  // index allong phi axis [1-fNPhi]
  //   Int_t where // 1 = PRE section, 0 = EC section, 2 = HC section
  // Outputs:
  //   none.
  // Returned
  //   Int_t index // Tower index number 
  
  if ( (ieta <= 0 || ieta>GetNEta()) || 
       (iphi <= 0 || iphi>GetNPhi())) 
    Fatal("TowerIndex", "Unexpected parameters eta = %d phi = %d!", ieta, iphi) ; 
  
  return ( (iphi - 1)*GetNEta() + ieta ); 
}

//______________________________________________________________________
void AliEMCALGeometry::TowerIndexes(Int_t index,Int_t &ieta,Int_t &iphi,
				    Int_t &ipre) const {
  // Inputs:
  //   Int_t index // Tower index number [1-i*fNZ*fNPhi] PRE(i=1)/ECAL(i=2)/HCAL(i=3)
  // Outputs:
  //   Int_t ieta    // index allong z axis [1-fNZ]
  //   Int_t iphi  // index allong phi axis [1-fNPhi]
  //   Int_t ipre  // 0 = ECAL section, 1 = Pre-shower section, 2 = HCAL section
  // Returned
  //   none.
  

  Int_t nindex = 0, itowers = GetNEta() * GetNPhi();

  if ( IsInPRE(index) ) {       // PRE index
    nindex = index - itowers;
    ipre = 1 ; 
  }
  else  if ( IsInECA(index) ) { // ECAL index
    nindex = index ;
    ipre = 0 ; 
  }
  else  if ( IsInHCA(index) ) { // HCAL index
    nindex = index - 2*itowers;
    ipre = 2 ; 
  }
  else 
    Fatal("TowerIndexes", "Unexpected Id number!") ;
   
  if (nindex%GetNZ()) 
    iphi = nindex / GetNZ() + 1 ; 
  else 
    iphi = nindex / GetNZ() ; 
  ieta = nindex - (iphi - 1) * GetNZ() ; 

  if (gDebug==2)
    Info("TowerIndexes", "index=%d,%d, ieta=%d, iphi = %d", index, nindex,ieta, iphi) ; 
  return;
  
}

//______________________________________________________________________
void AliEMCALGeometry::EtaPhiFromIndex(Int_t index,Float_t &eta,Float_t &phi) const {
    // given the tower index number it returns the based on the eta and phi
    // of the tower.
    // Inputs:
    //   Int_t index // Tower index number [1-i*fNZ*fNPhi] PRE(i=1)/ECAL(i=2)/HCAL(i=3)
    // Outputs:
    //   Float_t eta  // eta of center of tower in pseudorapidity
    //   Float_t phi  // phi of center of tower in degrees
    // Returned
    //   none.
    Int_t ieta, iphi, ipre ;
    Float_t deta, dphi ;

    TowerIndexes(index,ieta,iphi,ipre);
    
    if (gDebug == 2) 
      Info("EtaPhiFromIndex","index = %d, ieta = %d, iphi = %d", index, ieta, iphi) ;

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
Int_t AliEMCALGeometry::PreTowerIndexFromEtaPhi(Float_t eta,Float_t phi) const {
    // returns the pretower index number based on the eta and phi of the tower.
    // Inputs:
    //   Float_t eta  // eta of center of tower in pseudorapidity
    //   Float_t phi  // phi of center of tower in degrees
    // Outputs:
    //   none.
    // Returned
    //   Int_t index // PreTower index number [fNZ*fNPhi-2*fNZ*fNPhi]

    return GetNEta()*GetNPhi()+TowerIndexFromEtaPhi(eta,phi);
}

//______________________________________________________________________
Bool_t AliEMCALGeometry::AbsToRelNumbering(Int_t AbsId, Int_t *relid) const {
    // Converts the absolute numbering into the following array/
    //  relid[0] = EMCAL Arm number 1:1 
    //  relid[1] = 0  ECAL section ; = 1  PRE section; = 2 HCA section
    //  relid[2] = Row number inside EMCAL
    //  relid[3] = Column number inside EMCAL
    // Input:
    //   Int_t AbsId // Tower index number [1-2*fNZ*fNPhi]
    // Outputs:
    //   Int_t *relid // array of 5. Discribed above.
    Bool_t rv  = kTRUE ;
    Int_t ieta=0,iphi=0,ipre=0,index=AbsId;

    TowerIndexes(index,ieta,iphi,ipre);
    relid[0] = 1;
    relid[1] = ipre; 
    relid[2] = ieta;
    relid[3] = iphi;

    return rv;
}

//______________________________________________________________________
void AliEMCALGeometry::PosInAlice(const Int_t *relid, Float_t &theta, Float_t &phi) const 
{
  // Converts the relative numbering into the local EMCAL-module (x, z)
  // coordinates
  Int_t sect = relid[1]; // PRE/ECAL/HCAL section 1/0/2
  Int_t ieta = relid[2]; // offset along x axis
  Int_t iphi = relid[3]; // offset along z axis
  Int_t index;
  Float_t eta;
  
  index = TowerIndex(ieta,iphi);
  EtaPhiFromIndex(index,eta,phi);
  theta = 180.*(2.0*TMath::ATan(TMath::Exp(-eta)))/TMath::Pi();

    // correct for distance to IP different in PRE/ECAL/HCAL
  Float_t d = 0. ; 
  if (sect == 1)
    d = GetIP2PRESection() -  GetIPDistance() ; 
  else if (sect == 0)
    d = GetIP2ECASection() - GetIPDistance() ; 
  else if (sect == 2) 
    d = GetIP2HCASection() - GetIPDistance() ;
  else 
    Fatal("PosInAlice", "Unexpected tower section!") ; 

  Float_t correction = 1 + d/GetIPDistance() ; 
  Float_t tantheta = TMath::Tan(theta) * correction ; 
  theta = TMath::ATan(tantheta) * TMath::RadToDeg() ; 
  if (theta < 0 ) 
    theta += 180. ; 
  
  return;
}

//______________________________________________________________________
void AliEMCALGeometry::PosInAlice(const Int_t absid, Float_t &theta, Float_t &phi) const 
{
  // Converts the relative numbering into the local EMCAL-module (x, z)
  // coordinates
  
  Int_t relid[4] ; 
  AbsToRelNumbering(absid, relid) ;
  Int_t ieta = relid[2]; // offset along x axis
  Int_t iphi = relid[3]; // offset along z axis
  Int_t index;
  Float_t eta;
  
  index = TowerIndex(ieta,iphi);
  EtaPhiFromIndex(index,eta,phi);
  theta = 2.0*TMath::ATan(TMath::Exp(-eta)) ;
  
  // correct for distance to IP different in PRE/ECAL/HCAL
  Float_t d = 0. ; 
  if (IsInPRE(absid))
    d = GetIP2PRESection() -  GetIPDistance() ; 
  else if (IsInECA(absid))
    d = GetIP2ECASection() - GetIPDistance() ; 
  else if (IsInHCA(absid)) 
    d = GetIP2HCASection() - GetIPDistance() ;
  else 
    Fatal("PosInAlice", "Unexpected id # %d!", absid) ; 

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
    
    Float_t eta,theta, phi,cyl_radius=0. ;
    
    Int_t ieta   = relid[2]; // offset along x axis
    Int_t iphi = relid[3]; // offset along z axis
    Int_t ipre = relid[1]; // indicates 0 ECAL section, 1 PRE section, 2 HCAL section.
    Int_t index;
    
    index = TowerIndex(ieta,iphi);
    EtaPhiFromIndex(index,eta,phi);
    theta = 180.*(2.0*TMath::ATan(TMath::Exp(-eta)))/TMath::Pi();
    
    if ( ipre == 0 ) 
      cyl_radius = GetIP2ECASection() ;
    else if ( ipre == 1 ) 
      cyl_radius = GetIP2PRESection() ;
    else if ( ipre == 2 ) 
      cyl_radius = GetIP2HCASection() ;
    else 
      Fatal("XYZFromIndex", "Unexpected Tower section # %d", ipre) ;  

    Double_t  kDeg2Rad = TMath::DegToRad() ; 
    x =  cyl_radius * TMath::Cos(phi * kDeg2Rad ) ;
    y =  cyl_radius * TMath::Sin(phi * kDeg2Rad ) ; 
    z =  cyl_radius / TMath::Tan(theta * kDeg2Rad ) ; 
 
 return;
} 

//______________________________________________________________________
void AliEMCALGeometry::XYZFromIndex(const Int_t absid,  TVector3 &v) const {
    // given the tower relative number it returns the X, Y and Z
    // of the tower.
    
    // Outputs:
    //   Float_t x  // x of center of tower in cm
    //   Float_t y  // y of center of tower in cm
    //   Float_t z  // z of centre of tower in cm
    // Returned
    //   none.
    
    Float_t theta, phi,cyl_radius=0. ;
        
    PosInAlice(absid, theta, phi) ; 
    
    if ( IsInECA(absid) ) 
      cyl_radius = GetIP2ECASection() ;
    else if ( IsInPRE(absid) ) 
      cyl_radius = GetIP2PRESection() ;
    else if ( IsInHCA(absid) ) 
      cyl_radius = GetIP2HCASection() ;
    else 
      Fatal("XYZFromIndex", "Unexpected Tower section") ;  

    Double_t  kDeg2Rad = TMath::DegToRad() ; 
    v.SetX(cyl_radius * TMath::Cos(phi * kDeg2Rad ) );
    v.SetY(cyl_radius * TMath::Sin(phi * kDeg2Rad ) ); 
    v.SetZ(cyl_radius / TMath::Tan(theta * kDeg2Rad ) ) ; 
 
 return;
} 

//______________________________________________________________________
/*
Boot_t AliEMCALGeometry::AreNeighbours(Int_t index1,Int_t index2) const {
    // Returns kTRUE if the two towers are neighbours or not, including
    // diagonals. Both indexes are required to be either towers or preshower.
    // Inputs:
    //   Int_t index1  // index of tower 1
    //   Int_t index2  // index of tower 2
    // Outputs:
    //   none.
    // Returned
    //   Boot_t kTRUE if the towers are neighbours otherwise false.
    Boot_t anb = kFALSE;
    Int_t ieta1 = 0, ieta2 = 0, iphi1 = 0, iphi2 = 0, ipre1 = 0, ipre2 = 0;

    TowerIndexes(index1,ieta1,iphi1,ipre1);
    TowerIndexes(index2,ieta2,iphi2,ipre2);
    if(ipre1!=ipre2) return anb;
    if((ieta1>=ieta2-1 && ieta1<=ieta2+1) && (iphi1>=iphi2-1 &&iphi1<=iphi2+1))
                                                                 anb = kTRUE;
    return anb;
}
 */
