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
AliQACheckerBase("ITS","SDD Quality Assurance Checker"),
fkOnline(0),
fDet(0),  
fLDC(0),
fSPDOffset(0), 
fSDDOffset(0), 
fSSDOffset(0),
fSPDChecker(0),  // SPD Checker
fSDDChecker(0),  // SDD Checker
fSSDChecker(0)  // SSD Checker
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
const Double_t AliITSQAChecker::Check(AliQA::ALITASK_t index, TObjArray * list)
{

  // Super-basic check on the QA histograms on the input list:
  // look whether they are empty!
  Double_t spdCheck, sddCheck, ssdCheck;
  Double_t retval = 0.;
  if(fDet == 0 || fDet == 1) {
    AliDebug(1,"AliITSQAChecker::Create SPD Checker\n");
    if(!fSPDChecker) fSPDChecker = new AliITSQASPDChecker();
    spdCheck = fSPDChecker->Check(index, list, fSPDOffset);
    if(spdCheck>retval)retval = spdCheck;
  }
  if(fDet == 0 || fDet == 2) {
    AliDebug(1,"AliITSQAChecker::Create SDD Checker\n");
    if(!fSDDChecker) fSDDChecker = new AliITSQASDDChecker();
    sddCheck = fSDDChecker->Check(index, list, fSDDOffset);
    if(sddCheck>retval)retval = sddCheck;
  }
  if(fDet == 0 || fDet == 3) {
    AliDebug(1,"AliITSQAChecker::Create SSD Checker\n");
    if(!fSSDChecker) fSSDChecker = new AliITSQASSDChecker();
    //AliInfo(Form("Number of monitored objects SSD: %d", list->GetEntries()));	     
    ssdCheck = fSSDChecker->Check(index, list, fSSDOffset);
    if(ssdCheck>retval)retval = ssdCheck;
  }
  // here merging part for common ITS QA result
  // 
  return retval;
}


//____________________________________________________________________________
void AliITSQAChecker::SetTaskOffset(Int_t SPDOffset, Int_t SDDOffset, Int_t SSDOffset)
{
  //Setting the 3 offsets for each task called
  fSPDOffset = SPDOffset;
  fSDDOffset = SDDOffset;
  fSSDOffset = SSDOffset;
}
