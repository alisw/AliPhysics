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

//_________________________________________________________________________
// Implementation version v0 of PHOS Manager class 
// Layout EMC + PPSD has name GPS2  
// The main goal of this version of AliPHOS is to calculte the 
//  induced charged in the PIN diode, taking into account light
//  tracking in the PbWO4 crystal, induced signal in the 
//  PIN due to MIPS particle and electronic noise.
// This is done in the StepManager 
//                  
//*-- Author:  Odd Harald Oddland & Gines Martinez (SUBATECH)


// --- ROOT system ---
#include "TRandom.h"

// --- Standard library ---

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <strstream>

// --- AliRoot header files ---

#include "AliPHOSv1.h"
#include "AliPHOSHit.h"
#include "AliPHOSDigit.h"
#include "AliRun.h"
#include "AliConst.h"

ClassImp(AliPHOSv1)

//____________________________________________________________________________
AliPHOSv1::AliPHOSv1(const char *name, const char *title):
  AliPHOSv0(name,title)
{
  // ctor 

  // Number of electrons created in the PIN due to light collected in the PbWo4 crystal is calculated using 
  // following formula
  // NumberOfElectrons = EnergyLost * LightYield * PINEfficiency * 
  //                     exp (-LightYieldAttenuation * DistanceToPINdiodeFromTheHit) *
  //                     RecalibrationFactor ;
  // LightYield is obtained as a Poissonian distribution with a mean at 700000 photons per GeV fromValery Antonenko
  // PINEfficiency is 0.1875 from Odd Harald Odland work
  // k_0 is 0.0045 from Valery Antonenko 

  fLightYieldMean = 700000. ;
  fIntrinsicPINEfficiency = 0.1875 ;
  fLightYieldAttenuation = 0.0045 ;
  fRecalibrationFactor = 6.2 / fLightYieldMean ;
  fElectronsPerGeV = 2.77e+8 ; 
}

//____________________________________________________________________________
AliPHOSv1::AliPHOSv1(AliPHOSReconstructioner * Reconstructioner, const char *name, const char *title):
  AliPHOSv0(Reconstructioner,name,title)
{
  // ctor 

  // Number of electrons created in the PIN due to light collected in the PbWo4 crystal is calculated using 
  // following formula
  // NumberOfElectrons = EnergyLost * LightYield * PINEfficiency * 
  //                     exp (-LightYieldAttenuation * DistanceToPINdiodeFromTheHit) *
  //                     RecalibrationFactor ;
  // LightYield is obtained as a Poissonian distribution with a mean at 700000 photons per GeV fromValery Antonenko
  // PINEfficiency is 0.1875 from Odd Harald Odland work
  // k_0 is 0.0045 from Valery Antonenko 

  fLightYieldMean = 700000.;
  fIntrinsicPINEfficiency = 0.1875 ;
  fLightYieldAttenuation = 0.0045 ;
  fRecalibrationFactor = 6.2 / fLightYieldMean ;
  fElectronsPerGeV = 2.77e+8 ;
}

//____________________________________________________________________________
void AliPHOSv1::StepManager(void)
{
  // Accumulates hits as long as the track stays in a single crystal or PPSD gas Cell
  // Adds the energy deposited in the PIN diode

  Int_t          relid[4] ;      // (box, layer, row, column) indices
  Float_t        xyze[4] ;       // position wrt MRS and energy deposited
  TLorentzVector pos ;
  Int_t copy;
  Float_t        lightyield ;   // Light Yield per GeV
  Float_t        nElectrons ;   // Number of electrons in the PIN diode
  TString name = fGeom->GetName() ; 
  Float_t        global[3] ;
  Float_t        local[3] ;
  Float_t        lostenergy ;

  Int_t primary =  gAlice->GetPrimary( gAlice->CurrentTrack() ); 

  if ( name == "GPS2" ) { // the CPV is a PPSD
    if( gMC->CurrentVolID(copy) == gMC->VolId("GCEL") ) // We are inside a gas cell 
    {
      gMC->TrackPosition(pos) ;
      xyze[0] = pos[0] ;
      xyze[1] = pos[1] ;
      xyze[2] = pos[2] ;
      xyze[3] = gMC->Edep() ;


      if ( xyze[3] != 0 ) { // there is deposited energy 
       	gMC->CurrentVolOffID(5, relid[0]) ;  // get the PHOS Module number
       	gMC->CurrentVolOffID(3, relid[1]) ;  // get the Micromegas Module number 
      // 1-> Geom->GetNumberOfModulesPhi() *  fGeom->GetNumberOfModulesZ() upper                         
      //  >  fGeom->GetNumberOfModulesPhi()  *  fGeom->GetNumberOfModulesZ() lower
       	gMC->CurrentVolOffID(1, relid[2]) ;  // get the row number of the cell
        gMC->CurrentVolID(relid[3]) ;        // get the column number 

	// get the absolute Id number

	Int_t absid ; 
       	fGeom->RelToAbsNumbering(relid,absid) ; 
	

	AddHit(primary, absid, xyze );

      } // there is deposited energy 
     } // We are inside the gas of the CPV  
   } // GPS2 configuration
  
  if(gMC->CurrentVolID(copy) == gMC->VolId("PXTL") )//  We are inside a PBWO crystal 
     {
       gMC->TrackPosition(pos) ;
       xyze[0] = pos[0] ;
       xyze[1] = pos[1] ;
       xyze[2] = pos[2] ;
       lostenergy = gMC->Edep() ; 
       xyze[3] = gMC->Edep() ;

       global[0] = pos[0] ;
       global[1] = pos[1] ;
       global[2] = pos[2] ;

       if ( xyze[3] != 0 ) {
          gMC->CurrentVolOffID(10, relid[0]) ; // get the PHOS module number ;
          relid[1] = 0   ;                    // means PW04
          gMC->CurrentVolOffID(4, relid[2]) ; // get the row number inside the module
          gMC->CurrentVolOffID(3, relid[3]) ; // get the cell number inside the module

      // get the absolute Id number

          Int_t absid ; 
          fGeom->RelToAbsNumbering(relid,absid) ; 
	  gMC->Gmtod(global, local, 1) ;
	  
	  // calculating number of electrons in the PIN diode asociated to this hit
	  lightyield = gRandom->Poisson(fLightYieldMean) ;
	  nElectrons = lostenergy * lightyield * fIntrinsicPINEfficiency *
	    exp(-fLightYieldAttenuation * (local[1]+fGeom->GetCrystalSize(1)/2.0 ) ) ;

	  xyze[3] = nElectrons * fRecalibrationFactor ;
	  // add current hit to the hit list
          AddHit(primary, absid, xyze);
    
       } // there is deposited energy
    } // we are inside a PHOS Xtal

   if(gMC->CurrentVolID(copy) == gMC->VolId("PPIN") ) // We are inside de PIN diode 
     {
       gMC->TrackPosition(pos) ;
       xyze[0] = pos[0] ;
       xyze[1] = pos[1] ;
       xyze[2] = pos[2] ;
       lostenergy = gMC->Edep() ;
       xyze[3] = gMC->Edep() ;

       if ( xyze[3] != 0 ) {
          gMC->CurrentVolOffID(11, relid[0]) ; // get the PHOS module number ;
          relid[1] = 0   ;                    // means PW04 and PIN
          gMC->CurrentVolOffID(5, relid[2]) ; // get the row number inside the module
          gMC->CurrentVolOffID(4, relid[3]) ; // get the cell number inside the module

      // get the absolute Id number

          Int_t absid ; 
          fGeom->RelToAbsNumbering(relid,absid) ;
	  
	  // calculating number of electrons in the PIN diode asociated to this hit
	  nElectrons = lostenergy * fElectronsPerGeV ;
	  xyze[3] = nElectrons * fRecalibrationFactor ;

	  // add current hit to the hit list
          AddHit(primary, absid, xyze);
	  //printf("PIN volume is  %d, %d, %d, %d \n",relid[0],relid[1],relid[2],relid[3]);
	  //printf("Lost energy in the PIN is %f \n",lostenergy) ;
       } // there is deposited energy
    } // we are inside a PHOS XtalPHOS PIN diode
}

