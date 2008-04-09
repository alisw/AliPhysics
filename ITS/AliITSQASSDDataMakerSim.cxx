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

/* $Id$   */

//  *************************************************************
//  Checks the quality assurance 
//  by comparing with reference data
//  contained in a DB
//  -------------------------------------------------------------
//  W. Ferrarese + P. Cerello Feb 2008
//  INFN Torino

// --- ROOT system ---
#include <TTree.h>
// --- Standard library ---

// --- AliRoot header files ---
#include "AliITSQADataMakerSim.h"
#include "AliITSQASSDDataMakerSim.h"
#include "AliLog.h"
#include "AliQA.h"
#include "AliQAChecker.h"
#include "AliRawReader.h"

ClassImp(AliITSQASSDDataMakerSim)

//____________________________________________________________________________ 
AliITSQASSDDataMakerSim::AliITSQASSDDataMakerSim(AliITSQADataMakerSim *aliITSQADataMakerSim) :
TObject(),
fAliITSQADataMakerSim(aliITSQADataMakerSim),
fSSDhDigits(0),
fSSDhSDigits(0),
fSSDhHits(0),
fDigitsOffset(0),
fSDigitsOffset(0),
fHitsOffset(0)
{
  //ctor used to discriminate OnLine-Offline analysis   
}

//____________________________________________________________________________ 
AliITSQASSDDataMakerSim::AliITSQASSDDataMakerSim(const AliITSQASSDDataMakerSim& qadm) :
TObject(),
fAliITSQADataMakerSim(qadm.fAliITSQADataMakerSim),
fSSDhDigits(qadm.fSSDhDigits),
fSSDhSDigits(qadm.fSSDhSDigits),
fSSDhHits(qadm.fSSDhHits),
fDigitsOffset(qadm.fDigitsOffset),
fSDigitsOffset(qadm.fSDigitsOffset),
fHitsOffset(qadm.fHitsOffset)
{
  //copy ctor 
  fAliITSQADataMakerSim->SetName((const char*)qadm.fAliITSQADataMakerSim->GetName()) ; 
  fAliITSQADataMakerSim->SetTitle((const char*)qadm.fAliITSQADataMakerSim->GetTitle());
  }

//__________________________________________________________________
AliITSQASSDDataMakerSim& AliITSQASSDDataMakerSim::operator = (const AliITSQASSDDataMakerSim& qac )
{
  // Equal operator.
  this->~AliITSQASSDDataMakerSim();
  new(this) AliITSQASSDDataMakerSim(qac);
  return *this;
}

//____________________________________________________________________________ 
void AliITSQASSDDataMakerSim::StartOfDetectorCycle()
{
  //Detector specific actions at start of cycle
  AliDebug(1,"AliITSQADM::Start of SSD Cycle\n");
}

//____________________________________________________________________________ 
void AliITSQASSDDataMakerSim::EndOfDetectorCycle(AliQA::TASKINDEX_t /*task*/, TObjArray* /*list*/)
{
  // launch the QA checking
  AliDebug(1,"AliITSDM instantiates checker with Run(AliQA::kITS, task, list)\n"); 
  
  //AliQAChecker::Instance()->Run( AliQA::kITS , task, list);
}

//____________________________________________________________________________ 
void AliITSQASSDDataMakerSim::InitDigits()
{ 
  // Initialization for DIGIT data - SSD -
  fDigitsOffset = (fAliITSQADataMakerSim->fDigitsQAList)->GetEntries();

  // custom code here

  //fSSDhDigits must be incremented by one unit every time a histogram is ADDED to the QA List

  AliDebug(1,Form("%d SSD Digits histograms booked\n",fSSDhDigits));

}


//____________________________________________________________________________
void AliITSQASSDDataMakerSim::MakeDigits(TTree * /*digits*/)
{ 
  // Fill QA for DIGIT - SSD -
}




//____________________________________________________________________________ 
void AliITSQASSDDataMakerSim::InitSDigits()
{ 
  // Initialization for SDIGIT data - SSD -
  fSDigitsOffset = (fAliITSQADataMakerSim->fSDigitsQAList)->GetEntries();

  // custom code here

  //fSSDhSDigits must be incremented by one unit every time a histogram is ADDED to the QA List

  AliDebug(1,Form("%d SSD SDigits histograms booked\n",fSSDhSDigits));

}


//____________________________________________________________________________
void AliITSQASSDDataMakerSim::MakeSDigits(TTree * /*sdigits*/)
{ 
  // Fill QA for SDIGIT - SSD -
}





//____________________________________________________________________________ 
void AliITSQASSDDataMakerSim::InitHits()
{ 
  // Initialization for HITS data - SSD -
  fHitsOffset = (fAliITSQADataMakerSim->fHitsQAList)->GetEntries();

  // custom code here

  //fSSDhHits must be incremented by one unit every time a histogram is ADDED to the QA List

  AliDebug(1,Form("%d SSD Hits histograms booked\n",fSSDhHits));

}


//____________________________________________________________________________
void AliITSQASSDDataMakerSim::MakeHits(TTree * /*hits*/)
{ 
  // Fill QA for HITS - SSD -
}
