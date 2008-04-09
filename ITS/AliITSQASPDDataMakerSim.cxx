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

/* $Id$  */
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
#include "AliITSQASPDDataMakerSim.h"
#include "AliLog.h"
#include "AliQA.h"
#include "AliQAChecker.h"
#include "AliRawReader.h"

ClassImp(AliITSQASPDDataMakerSim)

//____________________________________________________________________________ 
AliITSQASPDDataMakerSim::AliITSQASPDDataMakerSim(AliITSQADataMakerSim *aliITSQADataMakerSim) :
TObject(),
fAliITSQADataMakerSim(aliITSQADataMakerSim),
fSPDhDigits(0),
fSPDhSDigits(0),
fSPDhHits(0),
fDigitsOffset(0),
fSDigitsOffset(0),
fHitsOffset(0)
{
  //ctor used to discriminate OnLine-Offline analysis   
}

//____________________________________________________________________________ 
AliITSQASPDDataMakerSim::AliITSQASPDDataMakerSim(const AliITSQASPDDataMakerSim& qadm) :
TObject(),
fAliITSQADataMakerSim(qadm.fAliITSQADataMakerSim),
fSPDhDigits(qadm.fSPDhDigits),
fSPDhSDigits(qadm.fSPDhSDigits),
fSPDhHits(qadm.fSPDhHits),
fDigitsOffset(qadm.fDigitsOffset),
fSDigitsOffset(qadm.fSDigitsOffset),
fHitsOffset(qadm.fHitsOffset)
{
  //copy ctor 
  fAliITSQADataMakerSim->SetName((const char*)qadm.fAliITSQADataMakerSim->GetName()) ; 
  fAliITSQADataMakerSim->SetTitle((const char*)qadm.fAliITSQADataMakerSim->GetTitle());
  }

//__________________________________________________________________
AliITSQASPDDataMakerSim& AliITSQASPDDataMakerSim::operator = (const AliITSQASPDDataMakerSim& qac )
{
  // Equal operator.
  this->~AliITSQASPDDataMakerSim();
  new(this) AliITSQASPDDataMakerSim(qac);
  return *this;
}

//____________________________________________________________________________ 
void AliITSQASPDDataMakerSim::StartOfDetectorCycle()
{
  //Detector specific actions at start of cycle
  AliDebug(1,"AliITSQADM::Start of SPD Cycle\n");
}

//____________________________________________________________________________ 
void AliITSQASPDDataMakerSim::EndOfDetectorCycle(AliQA::TASKINDEX_t /*task*/, TObjArray* /*list*/)
{
  // launch the QA checking
  AliDebug(1,"AliITSDM instantiates checker with Run(AliQA::kITS, task, list)\n"); 
  
  //AliQAChecker::Instance()->Run( AliQA::kITS , task, list);
}

//____________________________________________________________________________ 
void AliITSQASPDDataMakerSim::InitDigits()
{ 
  // Initialization for DIGIT data - SPD -
  fDigitsOffset = (fAliITSQADataMakerSim->fDigitsQAList)->GetEntries();

  // custom code here

  //fSPDhDigits must be incremented by one unit every time a histogram is ADDED to the QA List

  AliDebug(1,Form("%d SPD Digits histograms booked\n",fSPDhDigits));

}


//____________________________________________________________________________
void AliITSQASPDDataMakerSim::MakeDigits(TTree * /*digits*/)
{ 
  // Fill QA for DIGIT - SPD -
}




//____________________________________________________________________________ 
void AliITSQASPDDataMakerSim::InitSDigits()
{ 
  // Initialization for SDIGIT data - SPD -
  fSDigitsOffset = (fAliITSQADataMakerSim->fSDigitsQAList)->GetEntries();

  // custom code here

  //fSPDhSDigits must be incremented by one unit every time a histogram is ADDED to the QA List

  AliDebug(1,Form("%d SPD SDigits histograms booked\n",fSPDhSDigits));

}


//____________________________________________________________________________
void AliITSQASPDDataMakerSim::MakeSDigits(TTree * /*sdigits*/)
{ 
  // Fill QA for SDIGIT - SPD -
}





//____________________________________________________________________________ 
void AliITSQASPDDataMakerSim::InitHits()
{ 
  // Initialization for HITS data - SPD -
  fHitsOffset = (fAliITSQADataMakerSim->fHitsQAList)->GetEntries();

  // custom code here

  //fSPDhHits must be incremented by one unit every time a histogram is ADDED to the QA List

  AliDebug(1,Form("%d SPD Hits histograms booked\n",fSPDhHits));

}


//____________________________________________________________________________
void AliITSQASPDDataMakerSim::MakeHits(TTree * /*hits*/)
{ 
  // Fill QA for HITS - SPD -
}
