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

#include <iostream.h>

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

    if(!(  (strcmp( fName, "EMCALArch1a" ) == 0) |
	   (strcmp( fName, "EMCALArch1b" ) == 0) | 
	   (strcmp( fName, "EMCALArch2a" ) == 0) | 
	   (strcmp( fName, "EMCALArch2b" ) == 0) )){
	fgInit = kFALSE;
	cout <<"Instance " << fName << " undefined" << endl;
    } // end if
    fgInit = kTRUE; 

    // geometry 
    fAirGap     = 5.0; 
    fArm1PhiMin = 0.0; 
    fArm1PhiMax = 120.0; 

    fIPDistance = 454.0; 
    fZLength = 817.0; 
    fEnvelop[0] = fIPDistance; 
    fEnvelop[2] = fZLength; 
    fGap2Active = 1.0; 
    fShellThickness = 3.18 + 1.2 + (double)((2*fNLayers -3)/2);   
    fEnvelop[1] = fIPDistance + fShellThickness;

    if (((strcmp( fName, "EMCALArch1a" ))    == 0) |
	((strcmp( fName, "EMCALArch1b" ))    == 0)){
	fNZ         = 96;
	fNPhi       = 144;
    } // end if
    if (((strcmp( fName, "EMCALArch2a" ))    == 0) |
	((strcmp( fName, "EMCALArch2b" ))    == 0)){
	fNZ         = 112;
	fNPhi       = 168;
    } // end if
    if (((strcmp( fName, "EMCALArch1a" ))    == 0) |
	((strcmp( fName, "EMCALArch2a" ))    == 0)){
	fNLayers    = 21;
    } // end if
    if (((strcmp( fName, "EMCALArch1b" ))    == 0) |
	((strcmp( fName, "EMCALArch2b" ))    == 0)){
	fNLayers    = 25;
    } // end if
}
//______________________________________________________________________
AliEMCALGeometry *  AliEMCALGeometry::GetInstance(){ 
    // Returns the pointer of the unique instance

    return (AliEMCALGeometry *) fgGeom; 
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
	    cout << "AliEMCALGeometry <E> : current geometry is " 
		 << fgGeom->GetName() << endl
		 << "                      you cannot call     " << name 
		 << endl; 
	}else{
	    rv = (AliEMCALGeometry *) fgGeom; 
	} // end if
    }  // end if fgGeom
    return rv; 
}
//______________________________________________________________________
Int_t AliEMCALGeometry::TowerIndex(Int_t iz,Int_t iphi,Int_t ipre){
    // Returns the tower index number from the based on the Z and Phi
    // index numbers. There are 2 times the number of towers to separate
    // out the full towsers from the pre-towsers.
    // Inputs:
    //   Int_t iz    // index allong z axis [1-fNZ]
    //   Int_t iphi  // index allong phi axis [1-fNPhi]
    //   Int_t ipre  // 0 = Full tower, 1 = Pre-shower tower only. [0,1]
    // Outputs:
    //   none.
    // Returned
    // Int_t the absoulute tower index. [1-2*fNZ*fNPhi]
    Int_t index;

    if((iz<=0 || iz>GetNZ()) || (iphi<=0 || iphi>GetNPhi()) ||
       (ipre<0 || ipre>1) ){
	cout << "inputs out of range iz=" << iz << "[1-" << GetNZ();
	cout << "] iPhi=" << iphi << "[1-" << GetNPhi() << "] ipre=";
	cout << ipre << "[0,1]. returning -1" << endl;
	return -1;
    } // end if
    index = iphi + GetNPhi()*(iz-1) + ipre*(GetNPhi()*GetNZ());
    return index;
}
//______________________________________________________________________
void AliEMCALGeometry::TowerIndexes(Int_t index,Int_t &iz,Int_t &iphi,
				    Int_t &ipre){
    // given the tower index number it returns the based on the Z and Phi
    // index numbers and if it is for the full tower or the pre-tower number.
    // There are 2 times the number of towers to separate
    // out the full towsers from the pre-towsers.
    // Inputs:
    //   Int_t index // Tower index number [1-2*fNZ*fNPhi]
    // Outputs:
    //   Int_t iz    // index allong z axis [1-fNZ]
    //   Int_t iphi  // index allong phi axis [1-fNPhi]
    //   Int_t ipre  // 0 = Full tower, 1 = Pre-shower tower only. [0,1]
    // Returned
    //   none.
    Int_t itowers;

    itowers = GetNZ()*GetNPhi();
    if(index<1 || index>2*itowers){
	cout << "index=" << index <<" is out of range [1-";
	cout << 2*itowers << "], returning -1 for all." << endl;
	iz = -1; iphi = -1; ipre = -1;
	return ;
    } // end if
    ipre = 0;
    if(index>itowers){ // pre shower indexs
	ipre = 1;
	index = index - itowers;
    } // end if
    iz = 1+ (Int_t)(index/GetNPhi());
    iphi = index - GetNPhi()*(iz-1);
    return;
}
//______________________________________________________________________
void AliEMCALGeometry::EtaPhiFromIndex(Int_t index,Float_t &eta,Float_t &phi){
    // given the tower index number it returns the based on the eta and phi
    // of the tower.
    // Inputs:
    //   Int_t index // Tower index number [1-2*fNZ*fNPhi]
    // Outputs:
    //   Float_t eta  // eta of center of tower in pseudorapidity
    //   Float_t phi  // phi of center of tower in degrees
    // Returned
    //   none.
    Int_t iz,iphi,ipre;
    Double_t dz,dphi,zmax,z,phid,r;

    TowerIndexes(index,iz,iphi,ipre);
    zmax = (Double_t) GetZLength();
    dz   = zmax/((Double_t)GetNZ());
    r    = GetIPDistance();
    z    = dz*((Double_t)iz - 0.5);  // iz range [1-fNZ].
    eta  = -TMath::Log(TMath::Tan(0.5*TMath::ATan2(r,z)));
    dphi = GetArm1PhiMax() - GetArm1PhiMin();  // in degrees.
    phid = GetArm1PhiMin() + dphi*((Double_t)iphi -0.5);//iphi range [1-fNphi].
    phi  = phid;
}
//______________________________________________________________________
Int_t AliEMCALGeometry::TowerIndexFromEtaPhi(Float_t eta,Float_t phi){
    // returns the tower index number based on the eta and phi of the tower.
    // Inputs:
    //   Float_t eta  // eta of center of tower in pseudorapidity
    //   Float_t phi  // phi of center of tower in degrees
    // Outputs:
    //   none.
    // Returned
    //   Int_t index // Tower index number [1-fNZ*fNPhi]
    Int_t iz,iphi;
    Double_t z,zp,r,zl;

    r = GetIPDistance();
    z = 2.0*TMath::ATan(TMath::Exp(-eta));
    z = TMath::Tan(z);
    if(z!=0.0) z = r/z;
    else z = 0.0;
    zp = z;
    zl = GetZLength();
    z = 0.5*zl+z;
    iz = (Int_t)(((Double_t)GetNZ())*z/zl);
    if(iz<=0 || iz>GetNZ()){
	cout << "z=" << zp << " is outside of EMCAL. r=" << r << " eta =";
	cout << eta << " z length =" << zl ;
	cout << " returning -1" << endl;
	return -1;
    } // end if
    iphi =(Int_t)(((Double_t)GetNPhi())*((Double_t)phi)/
		  ((Double_t)(GetArm1PhiMax() - GetArm1PhiMin())));
    if(iphi<=0 || iphi>GetNPhi()){
	cout << "phi=" << phi << " is outside of EMCAL. r=" << r;
	cout << " Phimin=" << GetArm1PhiMin() << " PhiMax=" << GetArm1PhiMax();
	cout << " returning -1" << endl;
	return -1;
    } // end if
    return TowerIndex(iz,iphi,0);
}
//______________________________________________________________________
Int_t AliEMCALGeometry::PreTowerIndexFromEtaPhi(Float_t eta,Float_t phi){
    // returns the pretower index number based on the eta and phi of the tower.
    // Inputs:
    //   Float_t eta  // eta of center of tower in pseudorapidity
    //   Float_t phi  // phi of center of tower in degrees
    // Outputs:
    //   none.
    // Returned
    //   Int_t index // PreTower index number [fNZ*fNPhi-2*fNZ*fNPhi]

    return GetNZ()*GetNPhi()+TowerIndexFromEtaPhi(eta,phi);
}
//______________________________________________________________________
Bool_t AliEMCALGeometry::AbsToRelNumbering(Int_t AbsId, Int_t *relid){
    // Converts the absolute numbering into the following array/
    //  relid[0] = EMCAL Module number 1:1 (EMCAL arm number)
    //  relid[1] = 0  Not in Pre Shower layers
    //           = -1 In Pre Shower
    //  relid[2] = Row number inside EMCAL
    //  relid[3] = Column number inside EMCAL
    // Input:
    //   Int_t AbsId // Tower index number [1-2*fNZ*fNPhi]
    // Outputs:
    //   Int_t *relid // array of 5. Discribed above.
    Bool_t rv  = kTRUE ;
    Int_t iz=0,iphi=0,ipre=0,index=AbsId;

    TowerIndexes(index,iz,iphi,ipre);
    relid[0] = 1;
    relid[1] = 0;
    if(ipre==1) relid[1] = -1;
    relid[2] = iz;
    relid[3] = iphi;

    return rv;
}
//______________________________________________________________________
void AliEMCALGeometry::RelPosInModule(const Int_t *relid,Float_t &theta,
				     Float_t &phi){
    // Converts the relative numbering into the local PHOS-module (x, z)
    // coordinates
    Int_t iz   = relid[2]; // offset along x axis
    Int_t iphi = relid[3]; // offset along z axis
    Int_t ipre = relid[1]; // indecates -1 preshower, or 0 full tower.
    Int_t index;
    Float_t eta;

    if(ipre==-1) ipre = 1;
    index = TowerIndex(iz,iphi,ipre);
    EtaPhiFromIndex(index,eta,phi);
    theta = 180.*(2.0*TMath::ATan(TMath::Exp(-eta)))/TMath::Pi();

    return;
}
//______________________________________________________________________
/*
Boot_t AliEMCALGeometry::AreNeighbours(Int_t index1,Int_t index2){
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
    Int_t iz1 = 0, iz2 = 0, iphi1 = 0, iphi2 = 0, ipre1 = 0, ipre2 = 0;

    TowerIndexes(index1,iz1,iphi1,ipre1);
    TowerIndexes(index2,iz2,iphi2,ipre2);
    if(ipre1!=ipre2) return anb;
    if((iz1>=iz2-1 && iz1<=iz2+1) && (iphi1>=iphi2-1 && iphi1<=iphi2+1))
                                                                 anb = kTRUE;
    return anb;
}
 */
