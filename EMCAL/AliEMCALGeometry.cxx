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
// EMCAL consists of a shell of Pb 
//                  
//*-- Author: Yves Schutz (SUBATECH)

// --- ROOT system ---

// --- Standard library ---

#include <iostream.h>

// --- AliRoot header files ---

#include "AliEMCALGeometry.h"
#include "AliConst.h"

ClassImp(AliEMCALGeometry) ;

AliEMCALGeometry * AliEMCALGeometry::fgGeom = 0 ;
Bool_t            AliEMCALGeometry::fgInit = kFALSE ;

//____________________________________________________________________________
AliEMCALGeometry::~AliEMCALGeometry(void)
{
  // dtor

}

//____________________________________________________________________________

void AliEMCALGeometry::Init(void)
{
  // Initializes the EMCAL parameters

  fgInit = kTRUE ; 

  // geometry 
  fAirGap     = 5.0 ; 
  fArm1PhiMin = 130.0 ; 
  fArm1PhiMax = 210.0 ; 
  fArm2PhiMin = 330.0 ; 
  fArm2PhiMax = 410.0 ; 
  fIPDistance = 423.0 ; 
  fShellThickness = 50.0 ; 
  fZLength = 817.0 ; 
  fEnvelop[0] = fIPDistance ; 
  fEnvelop[1] = fIPDistance + fShellThickness ;
  fEnvelop[2] = fZLength ; 

  // material
  fAmat = 207.2;  
  fZmat = 82.;   
  fDmat = 5.798167 ;  
  fRmat = 1.261061;  
  fEmat = 23 ;  
 
}

//____________________________________________________________________________
AliEMCALGeometry *  AliEMCALGeometry::GetInstance() 
{ 
  // Returns the pointer of the unique instance
  return (AliEMCALGeometry *) fgGeom ; 
}

//____________________________________________________________________________
AliEMCALGeometry *  AliEMCALGeometry::GetInstance(const Text_t* name, const Text_t* title) 
{
  // Returns the pointer of the unique instance
  AliEMCALGeometry * rv = 0  ; 
  if ( fgGeom == 0 ) {
    if ( strcmp(name,"") == 0 ) 
      rv = 0 ;
    else {    
      fgGeom = new AliEMCALGeometry(name, title) ;
      if ( fgInit )
	rv = (AliEMCALGeometry * ) fgGeom ;
      else {
	rv = 0 ; 
	delete fgGeom ; 
	fgGeom = 0 ; 
      }
    }
  }
  else {
    if ( strcmp(fgGeom->GetName(), name) != 0 ) {
      cout << "AliEMCALGeometry <E> : current geometry is " << fgGeom->GetName() << endl
	   << "                      you cannot call     " << name << endl ; 
    }
    else
      rv = (AliEMCALGeometry *) fgGeom ; 
  } 
  return rv ; 
}

