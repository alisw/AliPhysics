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

/*
$Log$
Revision 1.25  2003/01/14 10:50:19  alibrary
Cleanup of STEER coding conventions

Revision 1.24  2002/11/21 17:01:56  alibrary
Removing AliMCProcess and AliMC

Revision 1.23  2002/10/23 07:24:57  alibrary
Introducing Riostream.h

Revision 1.22  2002/10/14 14:57:29  hristov
Merging the VirtualMC branch to the main development branch (HEAD)

Revision 1.20.6.2  2002/07/24 10:07:21  alibrary
Updating VirtualMC

Revision 1.21  2002/07/23 10:02:46  morsch
All volume names start with "S".

Revision 1.20  2001/10/30 12:18:45  morsch
Place station 3 into DDIP only if DDIP is present.

Revision 1.19  2001/07/17 09:51:38  morsch
Place station 3 inside Dipole.

Revision 1.18  2001/04/06 11:24:43  morsch
Dependency on implementations of AliSegmentation and AliMUONResponse moved to AliMUONFactory class.
Static method Build() builds the MUON system out of chambers, segmentation and response.

Revision 1.17  2001/03/17 10:07:20  morsch
Correct inconsistent variable name / method name / comments.

Revision 1.16  2001/01/27 08:50:50  morsch
Call non default constructors of segmentation classes.

Revision 1.15  2001/01/17 20:57:45  hristov
Unused variable removed

Revision 1.14  2000/12/21 22:42:55  morsch
Constructor contains default set-up for segmentation.
Record charged particles only.

Revision 1.13  2000/10/06 10:03:38  morsch
Call to gMC->VolId() moved to Init()

Revision 1.12  2000/10/02 21:28:09  fca
Removal of useless dependecies via forward declarations

Revision 1.11  2000/06/27 07:31:07  morsch
fChambers = 0; deleted from constructor.

Revision 1.10  2000/06/26 14:02:38  morsch
Add class AliMUONConstants with MUON specific constants using static memeber data and access methods.

Revision 1.9  2000/06/15 07:58:49  morsch
Code from MUON-dev joined

Revision 1.8.4.9  2000/06/12 19:20:49  morsch
Constructor sets default geometry, segmentation and response parameters.

Revision 1.8.4.8  2000/06/09 21:55:28  morsch
Most coding rule violations corrected.

Revision 1.8.4.7  2000/05/02 13:15:18  morsch
Coding rule violations RS3, RN13 corected

Revision 1.8.4.6  2000/05/02 10:24:26  morsch
Public access to fdGas and fdAlu of AliMUONChamber replaced by getters.

Revision 1.8.4.5  2000/04/26 19:58:47  morsch
Obsolete reference to trig_ removed.

Revision 1.8.4.4  2000/04/19 19:42:47  morsch
change NCH to kNCH

Revision 1.8.4.3  2000/02/17 08:17:43  morsch
Gammas and neutrons are also scored in the stepmanager
*/

/////////////////////////////////////////////////////////
//  Manager and hits classes for set:MUON version 0    //
/////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TClonesArray.h>
#include <TLorentzVector.h> 
#include <TNode.h> 
#include <TRandom.h> 
#include <TTUBE.h>

#include "AliMUONv0.h"
#include "AliMUONChamber.h"
#include "AliRun.h"
#include "AliMagF.h"
#include "AliMUONHit.h"
#include "AliMUONPadHit.h"
#include "AliCallf77.h"
#include "AliConst.h" 
#include "AliMUONConstants.h"
#include "AliMUONFactory.h"
ClassImp(AliMUONv0)
 
//___________________________________________
AliMUONv0::AliMUONv0() : AliMUON()
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
    AliMaterial(15, "AIR$      ", 14.61,  7.3, .001205, 30423.24, 67500);
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
	((AliMUONChamber*) (*fChambers)[i])->SetGid(gMC->VolId(vName));
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
      if(id==((AliMUONChamber*)(*fChambers)[i-1])->GetGid()){ 
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
	  AliMUONHit(fIshunt,gAlice->CurrentTrack(),vol,hits);

  }
//  if( gMC->IsTrackExiting()) gMC->StopTrack(); 
}









