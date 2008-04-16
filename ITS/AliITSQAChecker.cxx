/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
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

// *****************************************
//  Checks the quality assurance 
//  by comparing with reference data
//  W. Ferrarese Oct 2007
//  INFN Torino

// --- ROOT system ---
#include "TH1.h"

// --- AliRoot header files ---
#include "AliITSQAChecker.h"
#include "AliITSQASPDChecker.h"
#include "AliITSQASDDChecker.h"
#include "AliITSQASSDChecker.h"

ClassImp(AliITSQAChecker)

//____________________________________________________________________________
AliITSQAChecker::AliITSQAChecker(Bool_t kMode, Short_t subDet, Short_t ldc) :
AliQACheckerBase("ITS","SDD Quality Assurance Checker")
{
  fkOnline = kMode; fDet = subDet; fLDC = ldc;
  if(fDet == 0 || fDet == 1) {
    AliDebug(1,"AliITSQAChecker::Create SPD Checker\n");
  }
  if(fDet == 0 || fDet == 2) {
    AliDebug(1,"AliITSQAChecker::Create SDD Checker\n");
  }
  if(fDet == 0 || fDet == 3) {
    AliDebug(1,"AliITSQAChecker::Create SSD Checker\n");
  }

}

//__________________________________________________________________
AliITSQAChecker& AliITSQAChecker::operator = (const AliITSQAChecker& qac )
{
  // Equal operator.
  this->~AliITSQAChecker();
  new(this) AliITSQAChecker(qac);
  return *this;
}

//____________________________________________________________________________
const Double_t AliITSQAChecker::Check(TObjArray * list)
{

  // Super-basic check on the QA histograms on the input list:
  // look whether they are empty!
  if(fDet == 0 || fDet == 1) {
    AliDebug(1,"AliITSQAChecker::Create SPD Checker\n");
	if(!fSPDChecker) fSPDChecker = new AliITSQASPDChecker();
	Double_t SPDcheck = fSPDChecker->Check();
  }
  if(fDet == 0 || fDet == 2) {
    AliDebug(1,"AliITSQAChecker::Create SDD Checker\n");
	if(!fSDDChecker) fSDDChecker = new AliITSQASDDChecker();
	Double_t SDDcheck = fSDDChecker->Check();
  }
  if(fDet == 0 || fDet == 3) {
    AliDebug(1,"AliITSQAChecker::Create SSD Checker\n");
	if(!fSSDChecker) fSSDChecker = new AliITSQASSDChecker();
	Double_t SSDcheck = fSSDChecker->Check();
  }
  // here merging part for common ITS QA result
  return 0;
}


