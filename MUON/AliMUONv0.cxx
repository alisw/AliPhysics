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

/* $Id$ */

/////////////////////////////////////////////////////////
//  Manager and hits classes for set:MUON version 0    //
/////////////////////////////////////////////////////////

#include <TLorentzVector.h> 
#include <TVirtualMC.h>

#include "AliConst.h" 
#include "AliMUONChamber.h"
#include "AliMUONGeometryModule.h"
#include "AliMUONConstants.h"
#include "AliMUONFactory.h"
#include "AliMUONHit.h"
#include "AliMUONv0.h"
#include "AliMagF.h"
#include "AliRun.h"
#include "AliMC.h"

ClassImp(AliMUONv0)
 
//___________________________________________
AliMUONv0::AliMUONv0() 
  : AliMUON()
{
// Constructor
    fChambers = 0;
}
 
//___________________________________________
AliMUONv0::AliMUONv0(const char *name, const char *title)
  : AliMUON(name,title)
{
// Constructor
    AliMUONFactory factory;
    factory.Build(this, title);
}

//___________________________________________
void AliMUONv0::CreateGeometry()
{
// Creates coarse geometry for hit density simulations
    Int_t *idtmed = fIdtmed->GetArray()-1099;
//
     Float_t zpos, dAlu, tpar[3];
     Int_t idAir=idtmed[1100];
     Int_t idAlu=idtmed[1103];     

     AliMUONChamber *iChamber;
     // Loop over all chambers (tracking and trigger)
     for (Int_t ch = 0; ch < AliMUONConstants::NCh(); ch++) {
	 char alu[8];
	 char gas[8];
     
	 iChamber=(AliMUONChamber*) (*fChambers)[ch];
	 // Z of the chamber
	 zpos=iChamber->Z(); 
	 dAlu=iChamber->DAlu();
	 if (ch < AliMUONConstants::NTrackingCh()) {
	   // tracking chambers
	     sprintf(alu,"SA0%1d",ch);
	     sprintf(gas,"SG0%1d",ch);	 
	 } else {
	   // trigger chambers
	     sprintf(alu,"SA%2d",ch);
	     sprintf(gas,"SG%2d",ch);	 
	 }
//
	 tpar[0] = iChamber->RInner(); 
	 tpar[1] = iChamber->ROuter();
	 tpar[2] = (dAlu+0.2)/2.;
	 if (ch !=4 && ch !=5) {
	     gMC->Gsvolu(alu, "TUBE", idAlu, tpar, 3);
	     tpar[2] = 0.1;
	     gMC->Gsvolu(gas, "TUBE", idAir, tpar, 3);
	 } else {
	     gMC->Gsvolu(alu, "TUBE", idAlu, tpar, 3);
	     tpar[2] = 0.1;
	     gMC->Gsvolu(gas, "TUBE", idAir, tpar, 3);
	 }
	 gMC->Gspos(gas, 1, alu,  0., 0., 0., 0, "ONLY");
	 if (ch == 4 || ch ==5) {
	     if (gMC->VolId("DDIP")) {
		 gMC->Gspos(alu, 1, "DDIP", 0., 0., zpos, 0, "ONLY");
	     } else {
		 gMC->Gspos(alu, 1, "ALIC", 0., 0., zpos, 0, "ONLY");
	     }
	 } else {
	     gMC->Gspos(alu, 1, "ALIC", 0., 0., zpos, 0, "ONLY");
	 }
     }
}

//___________________________________________
void AliMUONv0::CreateMaterials()
{
// Creates materials for coarse geometry
// Air
  Float_t aAir[4]={12.0107,14.0067,15.9994,39.948};
  Float_t zAir[4]={6.,7.,8.,18.};
  Float_t wAir[4]={0.000124,0.755267,0.231781,0.012827};
  Float_t dAir = 1.20479E-3;
  AliMixture(15, "AIR$      ", aAir,  zAir, dAir,4, wAir);
  //  AliMaterial(15, "AIR$      ", 14.61,  7.3, .001205, 30423.24, 67500);
  AliMaterial( 9, "ALUMINIUM$", 26.98, 13. , 2.7, 8.9, 37.2);

  Float_t epsil  = .001; // Tracking precision, 
  Float_t stemax = -1.;  // Maximum displacement for multiple scat 
  Float_t tmaxfd = -20.; // Maximum angle due to field deflection 
  Float_t deemax = -.3;  // Maximum fractional energy loss, DLS 
  Float_t stmin  = -.8;
  Int_t isxfld   = gAlice->Field()->Integ();
  Float_t sxmgmx = gAlice->Field()->Max();
  
  //
  //    Air 
  AliMedium(1, "AIR_CH_US         ", 15, 1, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(4, "ALU_CH_US         ",  9, 0, isxfld, sxmgmx, tmaxfd, fMaxStepAlu, 
	    fMaxDestepAlu, epsil, stmin);
  
}

void AliMUONv0::Init()
{
   // 
   // Initialize Tracking Chambers
   //
    char vName[8];
    printf("\n\n\n Start Init for version 0 - CPC chamber type\n\n\n");
    for (Int_t i=0; i<AliMUONConstants::NCh(); i++) {
// Initialise chamber
	((AliMUONChamber*) (*fChambers)[i])->Init();
// Set sensitive volume Id
	if (i < AliMUONConstants::NTrackingCh()) {
	    // tracking chambers
	    sprintf(vName,"SG0%1d",i);	 
	} else {
	    // trigger chambers
	    sprintf(vName,"SG%2d",i);	 
	}
	//((AliMUONChamber*) (*fChambers)[i])->SetGid(gMC->VolId(vName));
	((AliMUONChamber*) (*fChambers)[i])
	   ->GetGeometry()->SetSensitiveVolume(gMC->VolId(vName));
    }
}

void AliMUONv0::StepManager()
{
//
// Step manager for hit density simulations
  Int_t          copy, id;
  static Int_t   idvol;
  static Int_t   vol[2];
  Int_t          ipart;
  TLorentzVector pos;
  TLorentzVector mom;
  Float_t        theta,phi;
  
  //  modifs perso
  static Float_t hits[15];

  TClonesArray &lhits = *fHits;
  //
  // Only gas gap inside chamber
  // Tag chambers and record hits when track enters 
  idvol=-1;
  id=gMC->CurrentVolID(copy);
  
    for (Int_t i=1; i<=AliMUONConstants::NCh(); i++) {
      //if(id==((AliMUONChamber*)(*fChambers)[i-1])->GetGid()){ 
      if ( ((AliMUONChamber*)(*fChambers)[i-1])->IsSensId(id) ) {
	  vol[0]=i; 
	  idvol=i-1;
      }
    }
    if (idvol == -1) return;
  //
  // Get current particle id (ipart), track position (pos)  and momentum (mom) 
  gMC->TrackPosition(pos);
  gMC->TrackMomentum(mom);

  ipart  = gMC->TrackPid();
  //
  // record hits when track enters ...
  if( !(gMC->TrackCharge()) ) return; 
  if( gMC->IsTrackEntering()) {
//      printf("\n Particle entering %f %f %f", pos[0], pos[1], pos[2] );
      
      Double_t tc = mom[0]*mom[0]+mom[1]*mom[1];
      Double_t rt = TMath::Sqrt(tc);
      theta   = Float_t(TMath::ATan2(rt,Double_t(mom[2])))*kRaddeg;
      phi     = Float_t(TMath::ATan2(Double_t(mom[1]),Double_t(mom[0])))*kRaddeg;
      hits[0] = Float_t(ipart);             // Geant3 particle type
      hits[1] = pos[0];                     // X-position for hit
      hits[2] = pos[1];                     // Y-position for hit
      hits[3] = pos[2];                     // Z-position for hit
      hits[4] = theta;                      // theta angle of incidence
      hits[5] = phi;                        // phi angle of incidence 
      hits[8] = -1;                         // first padhit
      hits[9] = -1;                         // last pad hit
      hits[10] = mom[3];                    // hit Energy
      hits[11] = mom[0];                    // Px
      hits[12] = mom[1];                    // Py
      hits[13] = mom[2];                    // Pz
      hits[14] = gMC->TrackTime();          // time of flight
      new(lhits[fNhits++]) 
	  AliMUONHit(fIshunt,gAlice->GetMCApp()->GetCurrentTrackNumber(),vol,hits);

  }
//  if( gMC->IsTrackExiting()) gMC->StopTrack(); 
}









