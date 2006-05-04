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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Transition Radiation Detector version 2 -- slow simulator with           //
//  detailed geometry                                                        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <stdlib.h> 

#include <TMath.h>
#include <TVirtualMC.h>

#include "AliConst.h"
#include "AliRun.h"
#include "AliTRDgeometry.h"
#include "AliTRDv2.h"

ClassImp(AliTRDv2)
 
//_____________________________________________________________________________
AliTRDv2::AliTRDv2():AliTRDv1()
{
  //
  // Default constructor
  //

}

//_____________________________________________________________________________
AliTRDv2::AliTRDv2(const char *name, const char *title) 
         :AliTRDv1(name, title) 
{
  //
  // Standard constructor for Transition Radiation Detector version 2
  //

  // Check that FRAME is there otherwise we have no place where to
  // put TRD
  AliModule* frame = gAlice->GetModule("FRAME");
  if (!frame) {
    Error("Ctor","TRD needs FRAME to be present\n");
    exit(1);
  } 

  if (frame->IsVersion() == 1) {
    // Detailed geometry without hole
    if (fGeometry) delete fGeometry;
    fGeometry = new AliTRDgeometry();
  }
  else {
    Error("Ctor","Could not find valid FRAME version 1\n");
    exit(1);
  }

}

//_____________________________________________________________________________
AliTRDv2::AliTRDv2(const AliTRDv2 &trd):AliTRDv1(trd)
{
  //
  // Copy constructor
  //

  ((AliTRDv2 &) trd).Copy(*this);

}

//_____________________________________________________________________________
AliTRDv2::~AliTRDv2()
{
  //
  // AliTRDv2 destructor
  //

}
 
//_____________________________________________________________________________
AliTRDv2 &AliTRDv2::operator=(const AliTRDv2 &trd)
{
  //
  // Assignment operator
  //

  if (this != &trd) ((AliTRDv2 &) trd).Copy(*this);
  return *this;

}
 
//_____________________________________________________________________________
void AliTRDv2::Copy(TObject &trd) const
{
  //
  // Copy function
  //

  AliTRDv1::Copy(trd); 

}

//_____________________________________________________________________________
void AliTRDv2::CreateGeometry()
{
  //
  // Create the geometry for the Transition Radiation Detector version 2
  //

  // Check that FRAME is there otherwise we have no place where to put the TRD
  AliModule* frame = gAlice->GetModule("FRAME");
  if (!frame) return;

  // Define the chambers
  AliTRD::CreateGeometry();

}

//_____________________________________________________________________________
void AliTRDv2::CreateMaterials()
{
  //
  // Create materials for the Transition Radiation Detector version 2
  //

  AliTRD::CreateMaterials();

}

