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
#include "AliITSQASPDDataMakerRec.h"
#include "AliLog.h"
#include "AliQA.h"
#include "AliQAChecker.h"
#include "AliRawReader.h"

ClassImp(AliITSQASPDDataMakerRec)

//____________________________________________________________________________ 
AliITSQASPDDataMakerRec::AliITSQASPDDataMakerRec(AliITSQADataMakerRec *aliITSQADataMakerRec, Bool_t kMode, Short_t ldc) :
TObject(),
fAliITSQADataMakerRec(aliITSQADataMakerRec),
fkOnline(kMode),
fLDC(ldc),
fSPDhRaws(0),
fSPDhRecs(0),
fRawsOffset(0),
fRecsOffset(0)
{
  //ctor used to discriminate OnLine-Offline analysis   
}

//____________________________________________________________________________ 
AliITSQASPDDataMakerRec::AliITSQASPDDataMakerRec(const AliITSQASPDDataMakerRec& qadm) :
TObject(),
fAliITSQADataMakerRec(qadm.fAliITSQADataMakerRec),
fkOnline(qadm.fkOnline),
fLDC(qadm.fLDC),
fSPDhRaws(qadm.fSPDhRaws),
fSPDhRecs(qadm.fSPDhRecs),
fRawsOffset(qadm.fRawsOffset),
fRecsOffset(qadm.fRecsOffset)
{
  //copy ctor 
  fAliITSQADataMakerRec->SetName((const char*)qadm.fAliITSQADataMakerRec->GetName()) ; 
  fAliITSQADataMakerRec->SetTitle((const char*)qadm.fAliITSQADataMakerRec->GetTitle());
  }

//__________________________________________________________________
AliITSQASPDDataMakerRec& AliITSQASPDDataMakerRec::operator = (const AliITSQASPDDataMakerRec& qac )
{
  // Equal operator.
  this->~AliITSQASPDDataMakerRec();
  new(this) AliITSQASPDDataMakerRec(qac);
  return *this;
}

//____________________________________________________________________________ 
void AliITSQASPDDataMakerRec::StartOfDetectorCycle()
{
  //Detector specific actions at start of cycle
  AliDebug(1,"AliITSQADM::Start of SPD Cycle\n");
}

//____________________________________________________________________________ 
void AliITSQASPDDataMakerRec::EndOfDetectorCycle(AliQA::TASKINDEX /*task*/, TObjArray* /*list*/)
{
  // launch the QA checking
  AliDebug(1,"AliITSDM instantiates checker with Run(AliQA::kITS, task, list)\n"); 
  
  //AliQAChecker::Instance()->Run( AliQA::kITS , task, list);
}

//____________________________________________________________________________ 
void AliITSQASPDDataMakerRec::InitRaws()
{ 
  // Initialization for RAW data - SPD -
  fRawsOffset = (fAliITSQADataMakerRec->fRawsQAList)->GetEntries();

  // custom code here

  //fSPDhRaws must be incremented by one unit every time a histogram is ADDED to the QA List

  AliDebug(1,Form("%d SPD Raws histograms booked\n",fSPDhRaws));

}


//____________________________________________________________________________
void AliITSQASPDDataMakerRec::MakeRaws(AliRawReader* /*rawReader*/)
{ 
  // Fill QA for RAW - SPD -
}

//____________________________________________________________________________ 
void AliITSQASPDDataMakerRec::InitRecPoints()
{
  // Initialization for RECPOINTS - SPD -
  fRecsOffset = (fAliITSQADataMakerRec->fRecPointsQAList)->GetEntries();

  // custom code here

  //fSPDhRecs must be incremented by one unit every time a histogram is ADDED to the QA List
  AliDebug(1,Form("%d SPD Recs histograms booked\n",fSPDhRecs));

}

//____________________________________________________________________________ 
void AliITSQASPDDataMakerRec::MakeRecPoints(TTree * /*clustersTree*/)
{
  // Fill QA for RecPoints - SPD -
}

