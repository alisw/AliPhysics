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
// Class for a QA checkable that is a Float  
//
//*-- Author :  Yves Schutz (SUBATECH) 
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---

#include "TClass.h"

// --- Standard library ---

// --- AliRoot header files ---

#include "AliPHOSQAFloatCheckable.h"

ClassImp(AliPHOSQAFloatCheckable)


//____________________________________________________________________________ 
  AliPHOSQAFloatCheckable::AliPHOSQAFloatCheckable(const char * name) : AliPHOSQAVirtualCheckable(name) 
{
  //ctor initial value is zero
  fType = "F" ;
  fValue = 0. ; 
}


//____________________________________________________________________________ 
  AliPHOSQAFloatCheckable::~AliPHOSQAFloatCheckable()
{
 // dtor
}


//____________________________________________________________________________ 
void AliPHOSQAFloatCheckable::Print() const
{
  // Print the chekable name and its value
  Info("Print", "Checkable-> %s : value = %f", GetName(), fValue) ; 
}

//____________________________________________________________________________ 
void AliPHOSQAFloatCheckable::Set(Float_t value)
{
  // Changes the value of the checkable and modifies its status
  fValue = value ; 
  fChange = kTRUE ; 
}

//____________________________________________________________________________ 
void AliPHOSQAFloatCheckable::Update(Float_t value)
{
  // Increments the value of the checkable and modifies its status
  fValue += value ; 
  fChange = kTRUE ; 
}
