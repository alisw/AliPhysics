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
#include <stdlib.h> 

// --- AliRoot header files ---
#include <TMath.h>

// -- ALICE Headers.
#include "AliConst.h"

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
void AliEMCALGeometry::Init(void){
    // Initializes the EMCAL parameters

    fgInit = kFALSE; // Assume failer untill proven otherwise.

    TString name(GetName()) ; 
		 
    if( name != "EMCALArch1a" &&
	name != "EMCALArch1b" && 
	name != "EMCALArch2a" && 
	name != "EMCALArch2b" && 
	name != "EMCALArch1aN" ){
      Fatal("Init", "%s is not a known geometry (choose among EMCALArch1a, EMCALArch1b, EMCALArch2a and EMCALArch2b, EMCALArch1aN)",  name.Data()) ;  
    } // end if
    //
    if ( name == "EMCALArch1a"  ||
	 name == "EMCALArch1b"  || 
	 name == "EMCALArch1aN") {
      fNZ         = 96;
      fNPhi       = 144;
    } // end if
    if ( name == "EMCALArch2a"  ||
	 name == "EMCALArch2b" ) {
	fNZ         = 112;
	fNPhi       = 168;
    } // end if
    if ( name == "EMCALArch1a"  ||
	 name == "EMCALArch2a" ) {
      fNPRLayers  = 2;
      fNECLayers  = 19;
      fNHCLayers  = 0;
    } // end if
    if ( name == "EMCALArch1b"  ||
	 name == "EMCALArch2b" ) {
	fNPRLayers  = 2;
	fNECLayers  = 23;
	fNHCLayers  = 0;
    } // end if
    if ( name == "EMCALArch1aN") { 
      fNPRLayers   = 2;
      fNECLayers  = 19;
      fNHCLayers  = 14;
    }

    // geometry
    fArm1PhiMin     =  60.0; // degrees, Starting EMCAL Phi position
    fArm1PhiMax     = 180.0; // degrees, Ending EMCAL Phi position
    fArm1EtaMin     = -0.7; // pseudorapidity, Starting EMCAL Eta position
    fArm1EtaMax     = +0.7; // pseudorapidity, Ending EMCAL Eta position

    fAlFrontThick        = 3.18; // cm, Thickness of front Al layer
    fGap2Active          = 1.0;  // cm, Gap between Al and 1st Scintillator
    fPbRadThickness      = 0.5;  // cm, Thickness of the Pb radiators.
    fPreShowerSintThick  = 0.6;  // cm, Thickness of the sintilator for the preshower part of the calorimeter
    fFullShowerSintThick = 0.5;  // cm, Thickness of the sintilator for the dull shower part of the calorimeter
    fCuRadThickness      = 0.0;  // cm, Thickness of the Cu radiators.

    if (name ==  "EMCALArch1aN") {
      fAlFrontThick        = 3.0;  // cm, Thickness of front Al layer
      fGap2Active          = 1.0;  // cm, Gap between Al and 1st Scintillator
      fPbRadThickness      = 0.6;  // cm, Thickness of the Pb radiators.
      fPreShowerSintThick  = 0.5;  // cm, Thickness of the sintilator for the preshower part of the calorimeter
      fFullShowerSintThick = 0.4;  // cm, Thickness of the sintilator for the full shower part of the calorimeter
      fCuRadThickness      = 1.0;  // cm, Thickness of the Cu radiators.
   }

    fIPDistance     = 454.0; // cm, Radial distance to inner surface of EMCAL
    fShellThickness = fAlFrontThick + fGap2Active + 2.*(GetPreSintThick() + GetPbRadThick()) + // pre shower 
      (fNECLayers-1)*(GetFullSintThick()+ GetPbRadThick()) + // E cal -1 because the last element is a scintillator
      fNHCLayers*(GetFullSintThick()+ GetCuRadThick()) + // H cal
      GetFullSintThick() ; // last scintillator
    fZLength        = 2.*ZFromEtaR(fIPDistance+fShellThickness,fArm1EtaMax); // Z coverage
    fEnvelop[0]     = fIPDistance; // mother volume inner radius
    fEnvelop[1]     = fIPDistance + fShellThickness; // mother volume outer r.
    fEnvelop[2]     = 1.00001*fZLength; // add some padding for mother volume. 
    fgInit = kTRUE; 
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
Int_t AliEMCALGeometry::TowerIndex(Int_t ieta,Int_t iphi,Int_t ipre) const {
    // Returns the tower index number from the based on the Z and Phi
    // index numbers. There are 2 times the number of towers to separate
    // out the full towers from the pre-showers.
    // Inputs:
    //   Int_t ieta    // index allong z axis [1-fNZ]
    //   Int_t iphi  // index allong phi axis [1-fNPhi]
    //   Int_t ipre  // 0 = Full tower, 1 = Pre-shower tower only. [0,1]
    // Outputs:
    //   none.
    // Returned
    // Int_t the absoulute tower index. [1-2*fNZ*fNPhi]
    Int_t index;

    if((ieta<=0 || ieta>GetNEta()) || (iphi<=0 || iphi>GetNPhi()) ||
       (ipre<0 || ipre>1) ){
      TString message ("\n") ; 
      message += "inputs out of range ieta= " ; 
      message += ieta ; 
      message += " [1-" ; 
      message += GetNEta() ;
      message += "] iphi= " ; 
      message += iphi ; 
      message += " [1-" ; 
      message += GetNPhi() ; 
      message += "] ipre= " ;
      message += ipre ; 
      message += "[0,1]. returning -1" ; 
      Warning("TowerIndex", message.Data() ) ; 
      return -1;
    } // end if
    index = iphi + GetNPhi()*(ieta-1) + ipre*(GetNPhi()*GetNEta());
    return index;
}

//______________________________________________________________________
void AliEMCALGeometry::TowerIndexes(Int_t index,Int_t &ieta,Int_t &iphi,
				    Int_t &ipre) const {
    // given the tower index number it returns the based on the Z and Phi
    // index numbers and if it is for the full tower or the pre-tower number.
    // There are 2 times the number of towers to separate
    // out the full towsers from the pre-towsers.
    // Inputs:
    //   Int_t index // Tower index number [1-2*fNZ*fNPhi]
    // Outputs:
    //   Int_t ieta    // index allong z axis [1-fNZ]
    //   Int_t iphi  // index allong phi axis [1-fNPhi]
    //   Int_t ipre  // 0 = Full tower, 1 = Pre-shower tower only. [0,1]
    // Returned
    //   none.
    Int_t itowers;

    itowers = GetNEta()*GetNPhi();
    if(index<1 || index>2*itowers){
      TString message("\n") ; 
      message += "index= " ; 
      message += index ; 
      message += " is out of range [1-" ;
      message += 2*itowers ; 
      message += "], returning -1 for all." ;
      Warning("TowerIndex", message.Data() ) ; 
      ieta = -1; iphi = -1; ipre = -1;
      return ;
    } // end if
    ipre = 0;
    if(index>itowers){ // pre shower indexs
	ipre = 1;
	index = index - itowers;
    } // end if
    ieta = 1+ (Int_t)((index-1)/GetNPhi());
    iphi = index - GetNPhi()*(ieta-1);
    return;
}

//______________________________________________________________________
void AliEMCALGeometry::EtaPhiFromIndex(Int_t index,Float_t &eta,Float_t &phi) const {
    // given the tower index number it returns the based on the eta and phi
    // of the tower.
    // Inputs:
    //   Int_t index // Tower index number [1-2*fNZ*fNPhi]
    // Outputs:
    //   Float_t eta  // eta of center of tower in pseudorapidity
    //   Float_t phi  // phi of center of tower in degrees
    // Returned
    //   none.
    Int_t ieta,iphi,ipre;
    Double_t deta,dphi,phid;

    TowerIndexes(index,ieta,iphi,ipre);
    deta = (GetArm1EtaMax()-GetArm1EtaMin())/((Float_t)GetNEta());
    eta  = GetArm1EtaMin() + (((Float_t)ieta)-0.5)*deta;
    dphi = (GetArm1PhiMax() - GetArm1PhiMin())/((Float_t)GetNPhi());  // in degrees.
    phid = GetArm1PhiMin() + dphi*((Float_t)iphi -0.5);//iphi range [1-fNphi].
    phi  = phid;
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

    ieta = 1 + (Int_t)(((Float_t)GetNEta())*(eta-GetArm1EtaMin())/
		  (GetArm1EtaMax() - GetArm1EtaMin()));
    if(ieta<=0 || ieta>GetNEta()){
      TString message("\n") ; 
      message += "ieta = " ; 
      message += ieta ; 
      message += " eta=" ; 
      message += eta ; 
      message += " is outside of EMCAL. etamin=" ;
      message += GetArm1EtaMin() ;
      message += " to etamax=" ; 
      message += GetArm1EtaMax();
      message += " returning -1";
      Warning("TowerIndexFromEtaPhi", message.Data() ) ; 
      return -1;
    } // end if
    iphi = 1 + (Int_t)(((Float_t)GetNPhi())*(phi-GetArm1PhiMin())/
		  ((Float_t)(GetArm1PhiMax() - GetArm1PhiMin())));
    if(iphi<=0 || iphi>GetNPhi()){
      TString message("\n") ; 
      message += "iphi=" ; 
      message += iphi ;  
      message += "phi= " ; 
      message += phi ; 
      message += " is outside of EMCAL." ;
      message += " Phimin=" ; 
      message += GetArm1PhiMin() ; 
      message += " PhiMax=" ; 
      message += GetArm1PhiMax() ;
      message += " returning -1" ;
      Warning("TowerIndexFromEtaPhi", message.Data() ) ; 
      return -1;
    } // end if
    return TowerIndex(ieta,iphi,0);
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
    //  relid[1] = 0  Not in Pre Shower layers
    //           = -1 In Pre Shower
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
    relid[1] = 0;
    if(ipre==1) 
      relid[1] = -1;
    relid[2] = ieta;
    relid[3] = iphi;

    return rv;
}

//______________________________________________________________________
void AliEMCALGeometry::PosInAlice(const Int_t *relid,Float_t &theta,
				     Float_t &phi) const {
    // Converts the relative numbering into the local EMCAL-module (x, z)
    // coordinates
    Int_t ieta   = relid[2]; // offset along x axis
    Int_t iphi = relid[3]; // offset along z axis
    Int_t ipre = relid[1]; // indicates -1 preshower, or 0 full tower.
    Int_t index;
    Float_t eta;

    if(ipre==-1) ipre = 1;
    index = TowerIndex(ieta,iphi,ipre);
    EtaPhiFromIndex(index,eta,phi);
    theta = 180.*(2.0*TMath::ATan(TMath::Exp(-eta)))/TMath::Pi();

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
    
    Float_t eta,theta, phi,cyl_radius,kDeg2Rad;
    
    Int_t ieta   = relid[2]; // offset along x axis
    Int_t iphi = relid[3]; // offset along z axis
    Int_t ipre = relid[1]; // indicates -1 preshower, or 0 full tower.
    Int_t index;
    

    if(ipre==-1) ipre = 1;
    index = TowerIndex(ieta,iphi,ipre);
    EtaPhiFromIndex(index,eta,phi);
    theta = 180.*(2.0*TMath::ATan(TMath::Exp(-eta)))/TMath::Pi();
    
    kDeg2Rad = TMath::Pi() / static_cast<Double_t>(180) ; 
    if ( ipre == -1 ) 
      cyl_radius = GetIP2PreShower() ;
    else 
      cyl_radius = GetIP2Tower() ;

    x =  cyl_radius * TMath::Cos(phi * kDeg2Rad ) ;
    y =  cyl_radius * TMath::Sin(phi * kDeg2Rad ) ; 
    z =  cyl_radius / TMath::Tan(theta * kDeg2Rad ) ; 
 
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
