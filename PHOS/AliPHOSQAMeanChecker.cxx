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
// QA checker that compares a number with an average value plus or minus
// a width 
//*-- Author :  Yves Schutz (SUBATECH) 
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---

#include "TDatime.h" 
#include "TFolder.h" 

// --- Standard library ---

// --- AliRoot header files ---

#include "AliPHOSQAMeanChecker.h" 
#include "AliPHOSQAAlarm.h"

ClassImp(AliPHOSQAMeanChecker)


//____________________________________________________________________________ 
  AliPHOSQAMeanChecker::AliPHOSQAMeanChecker(const char * name) : AliPHOSQAChecker(name,"") 
{
  // ctor
  SetTitle("checks against average value +/- width") ; 
}

//____________________________________________________________________________ 
  AliPHOSQAMeanChecker::AliPHOSQAMeanChecker(const char * name, Float_t mean, Float_t rms) : AliPHOSQAChecker(name,"") 
{
  // ctor
  SetTitle("checks against average value +/- width") ; 
  fMean = mean ; 
  fRms  = rms ; 
}

//____________________________________________________________________________ 
  AliPHOSQAMeanChecker::~AliPHOSQAMeanChecker()
{
  // dtor
}

//____________________________________________________________________________ 
TString AliPHOSQAMeanChecker::CheckingOperation()
{
  // The user defined checking operation
  // Return a non empty string in case the check was not satisfied

  TString rv ; 

  Float_t checked = 0. ;  
  if ( (fCheckable->HasA() == "I") &&  (fCheckable->HasA() == "F") ) {
    Error("CheckingOperation", "checker %s says you got the wrong checkable %s or the checkable has no value !", GetName(), fCheckable->GetName()) ;  
  } else {
    checked = fCheckable->GetValue(); 
    if (checked < fMean-fRms || checked > fMean+fRms) {
      char * tempo = new char[110] ;
      sprintf(tempo, "-->Checkable : %s :: Checker : %s :: Message : %f outside bond %f +/- %f\n", 
	      fCheckable->GetName(), GetName(), checked, fMean, fRms) ; 
      rv = tempo ;
      delete [] tempo ; 
    } 
  }  
  return rv ; 
} 

//____________________________________________________________________________ 
  void AliPHOSQAMeanChecker::Print()
{
  // print the name 
  
  Info("Print", "Checker : %s : %s : Mean = %f Rms = %f", GetName(), GetTitle(), fMean, fRms) ;  
}
