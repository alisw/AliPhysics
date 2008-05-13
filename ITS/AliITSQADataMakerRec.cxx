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
//  W. Ferrarese + P. Cerello Feb 2008
//  INFN Torino

// --- ROOT system ---
#include <TTree.h>
// --- Standard library ---

// --- AliRoot header files ---
#include "AliITSQADataMakerRec.h"
#include "AliITSQASPDDataMakerRec.h"
#include "AliITSQASDDDataMakerRec.h"
#include "AliITSQASSDDataMakerRec.h"
#include "AliLog.h"
#include "AliQA.h"
#include "AliQAChecker.h"
#include "AliITSQAChecker.h"
#include "AliRawReader.h"

ClassImp(AliITSQADataMakerRec)

//____________________________________________________________________________ 
AliITSQADataMakerRec::AliITSQADataMakerRec(Bool_t kMode, Short_t subDet, Short_t ldc) :
AliQADataMakerRec(AliQA::GetDetName(AliQA::kITS), "ITS Quality Assurance Data Maker"),
fkOnline(kMode),
fSubDetector(subDet),
fLDC(ldc),
fSPDDataMaker(NULL),
fSDDDataMaker(NULL),
fSSDDataMaker(NULL)
{
  //ctor used to discriminate OnLine-Offline analysis
  if(fSubDetector < 0 || fSubDetector > 3) {
	AliError("Error: fSubDetector number out of range; return\n");
  }

  // Initialization for RAW data 
  if(fSubDetector == 0 || fSubDetector == 1) {
    AliDebug(1,"AliITSQADM::Create SPD DataMakerRec\n");
	fSPDDataMaker = new AliITSQASPDDataMakerRec(this,fkOnline);
  }
  if(fSubDetector == 0 || fSubDetector == 2) {
	AliDebug(1,"AliITSQADM::Create SDD DataMakerRec\n");
	fSDDDataMaker = new AliITSQASDDDataMakerRec(this,fkOnline);
  }
  if(fSubDetector == 0 || fSubDetector == 3) {
	AliDebug(1,"AliITSQADM::Create SSD DataMakerRec\n");
	fSSDDataMaker = new AliITSQASSDDataMakerRec(this,fkOnline);
  }
}

//____________________________________________________________________________ 
AliITSQADataMakerRec::~AliITSQADataMakerRec(){
  // destructor
  if(fSPDDataMaker)delete fSPDDataMaker;
  if(fSDDDataMaker)delete fSDDDataMaker;
  if(fSSDDataMaker)delete fSSDDataMaker;
}

//____________________________________________________________________________ 
AliITSQADataMakerRec::AliITSQADataMakerRec(const AliITSQADataMakerRec& qadm) :
AliQADataMakerRec(),
fkOnline(qadm.fkOnline),
fSubDetector(qadm.fSubDetector),
fLDC(qadm.fLDC),
fSPDDataMaker(NULL),
fSDDDataMaker(NULL),
fSSDDataMaker(NULL)
{
  //copy ctor 
  SetName((const char*)qadm.GetName()) ; 
  SetTitle((const char*)qadm.GetTitle());
}

//__________________________________________________________________
AliITSQADataMakerRec& AliITSQADataMakerRec::operator = (const AliITSQADataMakerRec& qac )
{
  // Equal operator.
  this->~AliITSQADataMakerRec();
  new(this) AliITSQADataMakerRec(qac);
  return *this;
}

//____________________________________________________________________________ 
void AliITSQADataMakerRec::StartOfDetectorCycle()
{
  //Detector specific actions at start of cycle
  AliDebug(1,"AliITSQADM::Start of ITS Cycle\n");
  if(fSubDetector == 0 || fSubDetector == 1) fSPDDataMaker->StartOfDetectorCycle();
  if(fSubDetector == 0 || fSubDetector == 2) fSDDDataMaker->StartOfDetectorCycle();
  if(fSubDetector == 0 || fSubDetector == 3) fSSDDataMaker->StartOfDetectorCycle();
}

//____________________________________________________________________________ 
void AliITSQADataMakerRec::EndOfDetectorCycle(AliQA::TASKINDEX_t task, TObjArray* list)
{
  // launch the QA checking
  AliDebug(1,"AliITSDM instantiates checker with Run(AliQA::kITS, task, list)\n"); 
  if(fSubDetector == 0 || fSubDetector == 1) fSPDDataMaker->EndOfDetectorCycle(task, list);
  if(fSubDetector == 0 || fSubDetector == 2) fSDDDataMaker->EndOfDetectorCycle(task, list);
  if(fSubDetector == 0 || fSubDetector == 3) fSSDDataMaker->EndOfDetectorCycle(task, list);
  
  AliQAChecker *qac = AliQAChecker::Instance();
  AliITSQAChecker *qacb = (AliITSQAChecker *) qac->GetDetQAChecker(0);
  qacb->SetTaskOffset(fSPDDataMaker->GetOffset(), fSDDDataMaker->GetOffset(), fSSDDataMaker->GetOffset()); //Setting the offset for the QAChecker list
  qac->Run( AliQA::kITS , task, list);  //temporary skipping the checking
}

//____________________________________________________________________________ 
void AliITSQADataMakerRec::EndOfDetectorCycle(const char * /*fgDataName*/)
{
  //eventually used for different  AliQAChecker::Instance()->Run
}

//____________________________________________________________________________ 
void AliITSQADataMakerRec::InitRaws()
{  
  // Initialization for RAW data 
	if(fSubDetector == 0 || fSubDetector == 1) {
	  AliDebug(1,"AliITSQADM:: SPD InitRaws\n");
	  fSPDDataMaker->InitRaws();
	}
	if(fSubDetector == 0 || fSubDetector == 2) {
 	  AliDebug(1,"AliITSQADM:: SDD InitRaws\n");
	  fSDDDataMaker->InitRaws();
	}
	if(fSubDetector == 0 || fSubDetector == 3) {
	  AliDebug(1,"AliITSQADM:: SSD InitRaws\n");
	  fSSDDataMaker->InitRaws();
	}
}

//____________________________________________________________________________
void AliITSQADataMakerRec::MakeRaws(AliRawReader* rawReader)
{ 
  // Fill QA for RAW   
  if(fSubDetector == 0 || fSubDetector == 1) fSPDDataMaker->MakeRaws(rawReader);
  if(fSubDetector == 0 || fSubDetector == 2) fSDDDataMaker->MakeRaws(rawReader);
  if(fSubDetector == 0 || fSubDetector == 3) fSSDDataMaker->MakeRaws(rawReader);
}

//____________________________________________________________________________ 
void AliITSQADataMakerRec::InitRecPoints()
{
  // Initialization for RECPOINTS
  if(fSubDetector == 0 || fSubDetector == 1) {
	AliDebug(1,"AliITSQADM:: SPD InitRecPoints\n");
    fSPDDataMaker->InitRecPoints();
  }
  if(fSubDetector == 0 || fSubDetector == 2) {
	AliDebug(1,"AliITSQADM:: SDD InitRecPoints\n");
	fSDDDataMaker->InitRecPoints();
  }
  if(fSubDetector == 0 || fSubDetector == 3) {
	AliDebug(1,"AliITSQADM:: SSD InitRecPoints\n");
	fSSDDataMaker->InitRecPoints();
  }
}

//____________________________________________________________________________ 
void AliITSQADataMakerRec::MakeRecPoints(TTree * clustersTree)
{
  // Fill QA for recpoints
  if(fSubDetector == 0 || fSubDetector == 1) fSPDDataMaker->MakeRecPoints(clustersTree);
  if(fSubDetector == 0 || fSubDetector == 2) fSDDDataMaker->MakeRecPoints(clustersTree);
  if(fSubDetector == 0 || fSubDetector == 3) fSSDDataMaker->MakeRecPoints(clustersTree);
}


