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
// Class for a QA checkable that is an Int  
// To be used with AliPHOSChecker
// or any derived class
// ..
//*-- Author :  Yves Schutz (SUBATECH) 
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---

// --- Standard library ---

// --- AliRoot header files ---

#include "AliPHOSQAIntCheckable.h"

ClassImp(AliPHOSQAIntCheckable)


//____________________________________________________________________________ 
  AliPHOSQAIntCheckable::AliPHOSQAIntCheckable(const char * name) : AliPHOSQAVirtualCheckable(name) 
{
  //ctor initial value is zero
  fType  = "I" ; 
  fValue = 0 ; 
}


//____________________________________________________________________________ 
  AliPHOSQAIntCheckable::~AliPHOSQAIntCheckable()
{
  // dtor
}

//____________________________________________________________________________ 
void AliPHOSQAIntCheckable::Print() const 
{
  // Print the chekable name and its value
  Info("Print", "Checkable-> %s : value = %d", GetName(), fValue) ; 
}

//____________________________________________________________________________ 
void AliPHOSQAIntCheckable::Set(Int_t value)
{
  // Changes the value of the checkable and modifies its status
  fValue = value ; 
  fChange = kTRUE ; 
}

//____________________________________________________________________________ 
void AliPHOSQAIntCheckable::Update(Int_t value) 
{
  // Increments the value of the checkable and modifies its status
  fValue += value ; 
  fChange = kTRUE ; 
}
