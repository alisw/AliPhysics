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
//  Simplified TRD slow simulator                                            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
 
#include <stdlib.h>
 
#include "AliTRDsimple.h"
#include "AliTRDsimpleGen.h" 
#include "AliTRDsimpleMC.h"

ClassImp(AliTRDsimple)
 
//_____________________________________________________________________________
AliTRDsimple::AliTRDsimple():TObject()
{                       
  //
  // AliTRDsimple default constructor
  //

  fGenerator = NULL;
                                                         
}                                                                               
 
//_____________________________________________________________________________
AliTRDsimple::AliTRDsimple(const AliTRDsimple &s):TObject(s)
{
  //
  // AliTRDsimple copy constructor
  //
 
  ((AliTRDsimple &) s).Copy(*this);
 
}
 
//_____________________________________________________________________________
AliTRDsimple::~AliTRDsimple()
{
  //
  // AliTRDsimple destructor
  //
 
  if (fGenerator) {
    delete fGenerator;
  }

}                                                                               
 
//_____________________________________________________________________________
AliTRDsimple &AliTRDsimple::operator=(const AliTRDsimple &s)
{
  //
  // Assignment operator
  //
 
  if (this != &s) ((AliTRDsimple &) s).Copy(*this);
  return *this;
 
}
 
//_____________________________________________________________________________
void AliTRDsimple::Init()
{
  //
  // Initialization
  //

  fGenerator = new AliTRDsimpleGen();

  // Create the MC object
  new AliTRDsimpleMC("simple","Simplified Monte Carlo");
                                                         
}
 
//_____________________________________________________________________________
void AliTRDsimple::Copy(TObject &s) const
{
  //
  // Copy function
  //                             
                  
  fGenerator->Copy(*((AliTRDsimple &) s).fGenerator);  
                              
}
                                                                                
//_____________________________________________________________________________
void AliTRDsimple::ProcessEvent(Int_t ievent)
{
  //
  // Runs a single event
  //

  Int_t copy = 0;

  // Generate a new particle
  fGenerator->NewParticle(ievent);

  // Track the event
  do {
    gMC->ProcessEvent();
  } 
  while (gMC->CurrentVolID(copy) != -1);

}
                                                                                
