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
//
//*-- Author :  Yves Schutz (SUBATECH) 
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---

#include "TClass.h"

// --- Standard library ---

// --- AliRoot header files ---

#include "AliPHOSQAObjectCheckable.h"

ClassImp(AliPHOSQAObjectCheckable)


//____________________________________________________________________________ 
  AliPHOSQAObjectCheckable::AliPHOSQAObjectCheckable(const char * name) : AliPHOSQAVirtualCheckable(name) 
{
  fType   = "O" ; 
  fObject = 0 ; 
}

//____________________________________________________________________________ 
  AliPHOSQAObjectCheckable::~AliPHOSQAObjectCheckable()
{

}

//____________________________________________________________________________ 
void AliPHOSQAObjectCheckable::Print() const
{
  Info("Print", "Checkable-> %s : value = ", GetName()) ; 
  if( fObject ) 
    fObject->Print() ;
  else
    Info("Print", "no object specified yet" ) ; 
}

