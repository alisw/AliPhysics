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

//  *************************************************************
//  Checks the quality assurance 
//  by comparing with reference data
//  contained in a DB
//  -------------------------------------------------------------
//  W. Ferrarese Nov 2007
//  INFN Torino

// --- ROOT system ---
#include <TH2D.h>
#include <TTree.h>

// --- Standard library ---

// --- AliRoot header files ---
#include "AliITSQADataMakerSim.h"
#include "AliLog.h"
#include "AliQAChecker.h"


ClassImp(AliITSQADataMakerSim)

//____________________________________________________________________________ 
AliITSQADataMakerSim::AliITSQADataMakerSim() : 
  AliQADataMakerSim(AliQA::GetDetName(AliQA::kITS), "SDD Quality Assurance Data Maker")
{ 
  // ctor 
}

//____________________________________________________________________________ 
AliITSQADataMakerSim::AliITSQADataMakerSim(Int_t /*ldc */, Bool_t /*kMode */) :
  AliQADataMakerSim(AliQA::GetDetName(AliQA::kITS), "SDD Quality Assurance Data Maker")
{
  //ctor used to discriminate OnLine-Offline analysis
}

//____________________________________________________________________________ 
AliITSQADataMakerSim::AliITSQADataMakerSim(const AliITSQADataMakerSim& qadm) :
  AliQADataMakerSim()
{
  //copy ctor 
  SetName((const char*)qadm.GetName()) ; 
  SetTitle((const char*)qadm.GetTitle()); 
}

//__________________________________________________________________
AliITSQADataMakerSim& AliITSQADataMakerSim::operator = (const AliITSQADataMakerSim& qac )
{
  // Equal operator.
  this->~AliITSQADataMakerSim();
  new(this) AliITSQADataMakerSim(qac);
  return *this;
}

//____________________________________________________________________________ 
void AliITSQADataMakerSim::StartOfDetectorCycle() const
{
  //Detector specific actions at start of cycle
  AliDebug(1,"AliITSQADM::Start of ITS Cycle\n");
}

//____________________________________________________________________________ 
void AliITSQADataMakerSim::EndOfDetectorCycle(AliQA::TASKINDEX task, TObjArray *list)
{
  // launch the QA checking
  AliDebug(1,"AliITSDM instantiates checker with Run(AliQA::kITS, task, list)\n"); 
  
  AliQAChecker::Instance()->Run( AliQA::kITS , task, list);
}

