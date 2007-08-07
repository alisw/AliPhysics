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

/*
  Checks the quality assurance. 
  By comparing with reference data
  Y. Schutz CERN July 2007
*/

// --- ROOT system ---
#include <TClass.h>
#include <TH1F.h> 
#include <TH1I.h> 
#include <TIterator.h> 
#include <TKey.h> 
#include <TFile.h> 

// --- Standard library ---

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliQualAss.h"
#include "AliQualAssChecker.h"
#include "AliPHOSQualAssChecker.h"

ClassImp(AliPHOSQualAssChecker)

//__________________________________________________________________
AliPHOSQualAssChecker& AliPHOSQualAssChecker::operator = (const AliPHOSQualAssChecker& qadm )
{
  // Equal operator.
  this->~AliPHOSQualAssChecker();
  new(this) AliPHOSQualAssChecker(qadm);
  return *this;
}

//____________________________________________________________________________
const Double_t AliPHOSQualAssChecker::Check(const Option_t * opt) 
{
  // Performs the checking

  TDirectory * wRefDir = fDetectorDir->GetDirectory(opt) ; 
  TDirectory * wInDir  = fDetectorDir ->GetDirectory(opt) ; 
  Double_t test = 0.0  ;

  if (!wRefDir || !wInDir) 
    test = -1. ;
  else {
    TList * keyList = wRefDir->GetListOfKeys() ; 
    TIter next(keyList) ; 
    TKey * key ;
    Int_t count = 0 ; 
    while ( (key = static_cast<TKey *>(next())) ) {
      TObject * oref = wRefDir->Get(key->GetName()) ; 
      if ( oref->IsA()->InheritsFrom("TH1") ) {
	TH1 * href = static_cast<TH1F*>(oref) ; 
	TH1 * hin  = static_cast<TH1F*>(wInDir->Get(key->GetName())) ; 
	test += DiffK(href, hin) ;
	AliInfo(Form("test = %f", test)) ; 
	count++ ; 
	  } else
	  AliError(Form("%s is a class name that cannot be processed", key->GetClassName())) ;
    }
    if (count != 0) 
      test /= count ;
  }
    
  return test ; 
}  
