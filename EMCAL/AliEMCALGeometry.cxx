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
//*-- Author: Sahal Yacoob (LBL / UCT)
//     and  : Yves Schutz (SUBATECH)
//     and  : Jennifer Klay (LBL)

// --- AliRoot header files ---
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
Bool_t AliEMCALGeometry::AreInSameTower(Int_t id1, Int_t id2) const {
  // Find out whether two hits are in the same tower
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

  fgInit = kFALSE; // Assume failed until proven otherwise.
  TString name(GetName()) ; 
  if (name == "EMCAL_55_25") {
    fECPbRadThickness  = 0.5;  // cm, Thickness of the Pb radiators
    fECScintThick      = 0.5;  // cm, Thickness of the scintillator
    fNECLayers         = 25;   // number of scintillator layers
    
    fSampling          = 11.8; 
 
    fAlFrontThick      = 3.5;  // cm, Thickness of front Al layer
    fGap2Active        = 1.0;  // cm, Gap between Al and 1st Scintillator
  }
  else if( name == "G56_2_55_19" || name == "EMCAL_5655_21" || name == "G56_2_55_19_104_14"|| name == "G65_2_64_19" || name == "EMCAL_6564_21"){
    Fatal("Init", "%s is an old geometry! Please update your Config file", name.Data()) ;
  }
  else
    Fatal("Init", "%s is an undefined geometry!", name.Data()) ; 
		 
  // geometry
  fNZ             = 114;	// granularity along Z (eta) 
  fNPhi           = 168;	// granularity in phi (azimuth)
  fArm1PhiMin     = 60.0;	// degrees, Starting EMCAL Phi position
  fArm1PhiMax     = 180.0;	// degrees, Ending EMCAL Phi position
  fArm1EtaMin     = -0.7;	// pseudorapidity, Starting EMCAL Eta position
  fArm1EtaMax     = +0.7;	// pseudorapidity, Ending EMCAL Eta position
  
  fIPDistance     = 454.0; // cm, Radial distance to inner surface of EMCAL

  //There is always one more scintillator than radiator layer because of the first block of aluminium
  fShellThickness = fAlFrontThick + fGap2Active + fNECLayers*GetECScintThick()+(fNECLayers-1)*GetECPbRadThick();

  fZLength        = 2.*ZFromEtaR(fIPDistance+fShellThickness,fArm1EtaMax); // Z coverage
  fEnvelop[0]     = fIPDistance; // mother volume inner radius
  fEnvelop[1]     = fIPDistance + fShellThickness; // mother volume outer r.
  fEnvelop[2]     = 1.00001*fZLength; // add some padding for mother volume. 
  
  fgInit = kTRUE; 
  
  if (gDebug) {
    printf("Init: geometry of EMCAL named %s is as follows:", name.Data());
    printf( "               ECAL      : %d x (%f mm Pb, %f mm Sc) \n", GetNECLayers(), GetECPbRadThick(), GetECScintThick() ) ; 
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
