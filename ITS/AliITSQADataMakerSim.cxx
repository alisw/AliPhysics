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
#include "AliITSQADataMakerSim.h"
#include "AliITSQASPDDataMakerSim.h"
#include "AliITSQASDDDataMakerSim.h"
#include "AliITSQASSDDataMakerSim.h"
#include "AliLog.h"
#include "AliQA.h"
#include "AliQAChecker.h"
#include "AliITSQAChecker.h"
#include "AliRawReader.h"

ClassImp(AliITSQADataMakerSim)

//____________________________________________________________________________ 
AliITSQADataMakerSim::AliITSQADataMakerSim(Short_t subDet) :
AliQADataMakerSim(AliQA::GetDetName(AliQA::kITS), "ITS Quality Assurance Data Maker"),
fSubDetector(subDet),
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
    AliDebug(1,"AliITSQADM::Create SPD DataMakerSim\n");
	fSPDDataMaker = new AliITSQASPDDataMakerSim(this);
  }
  if(fSubDetector == 0 || fSubDetector == 2) {
	AliDebug(1,"AliITSQADM::Create SDD DataMakerSim\n");
	//printf("AliITSQADM::Create SDD DataMakerSim\n");		    
	fSDDDataMaker = new AliITSQASDDDataMakerSim(this);
  }
  if(fSubDetector == 0 || fSubDetector == 3) {
	AliDebug(1,"AliITSQADM::Create SSD DataMakerSim\n");
	fSSDDataMaker = new AliITSQASSDDataMakerSim(this);
  }
}

//____________________________________________________________________________ 
AliITSQADataMakerSim::~AliITSQADataMakerSim(){
  // destructor
  if(fSPDDataMaker)delete fSPDDataMaker;
  if(fSDDDataMaker)delete fSDDDataMaker;
  if(fSSDDataMaker)delete fSSDDataMaker;
}

//____________________________________________________________________________ 
AliITSQADataMakerSim::AliITSQADataMakerSim(const AliITSQADataMakerSim& qadm) :
AliQADataMakerSim(),
fSubDetector(qadm.fSubDetector),
fSPDDataMaker(NULL),
fSDDDataMaker(NULL),
fSSDDataMaker(NULL)
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
void AliITSQADataMakerSim::StartOfDetectorCycle()
{
  //Detector specific actions at start of cycle
  AliDebug(1,"AliITSQADM::Start of ITS Cycle\n");

  if(fSubDetector == 0 || fSubDetector == 1) fSPDDataMaker->StartOfDetectorCycle();
  if(fSubDetector == 0 || fSubDetector == 2) fSDDDataMaker->StartOfDetectorCycle();
  if(fSubDetector == 0 || fSubDetector == 3) fSSDDataMaker->StartOfDetectorCycle();
}

//____________________________________________________________________________ 
void AliITSQADataMakerSim::EndOfDetectorCycle(AliQA::TASKINDEX_t task, TObjArray* list)
{
  // launch the QA checking
  AliDebug(1,"AliITSDM instantiates checker with Run(AliQA::kITS, task, list)\n"); 
  if(fSubDetector == 0 || fSubDetector == 1) fSPDDataMaker->EndOfDetectorCycle(task, list);
  if(fSubDetector == 0 || fSubDetector == 2) fSDDDataMaker->EndOfDetectorCycle(task, list);
  if(fSubDetector == 0 || fSubDetector == 3) fSSDDataMaker->EndOfDetectorCycle(task, list);
  
  AliQAChecker *qac = AliQAChecker::Instance();
  AliITSQAChecker *qacb = (AliITSQAChecker *) qac->GetDetQAChecker(0);
  if(fSubDetector == 0 ) {
                Int_t offsetSPD = fSPDDataMaker->GetOffset(task); //+ fSPDDataMaker->GetOffsetS() + fSPDDataMaker->GetOffsetD() ;  
		Int_t offsetSDD = fSDDDataMaker->GetOffset(task); 
		Int_t offsetSSD = fSSDDataMaker->GetOffset(task);// + fSSDDataMaker->GetOffsetS() + fSSDDataMaker->GetOffsetD() ; 
    qacb->SetTaskOffset(offsetSPD, offsetSDD, offsetSSD); //Setting the offset for the QAChecker list		
	}
	qac->Run( AliQA::kITS , task, list);  //temporary skipping the checking
}

//____________________________________________________________________________ 
void AliITSQADataMakerSim::InitDigits()
{  
  // Initialization for RAW data 
	if(fSubDetector == 0 || fSubDetector == 1) {
	  AliDebug(1,"AliITSQADM:: SPD InitDigits\n");
	  fSPDDataMaker->InitDigits();
	}
	if(fSubDetector == 0 || fSubDetector == 2) {
 	  AliDebug(1,"AliITSQADM:: SDD InitDigits\n");
	  fSDDDataMaker->InitDigits();
	}
	if(fSubDetector == 0 || fSubDetector == 3) {
	  AliDebug(1,"AliITSQADM:: SSD InitDigits\n");
	  fSSDDataMaker->InitDigits();
	}
}

//____________________________________________________________________________
void AliITSQADataMakerSim::MakeDigits(TClonesArray * digits)
{ 
  // Fill QA for RAW   
  if(fSubDetector == 0 || fSubDetector == 1) fSPDDataMaker->MakeDigits(digits);
  if(fSubDetector == 0 || fSubDetector == 2) fSDDDataMaker->MakeDigits(digits);
  if(fSubDetector == 0 || fSubDetector == 3) fSSDDataMaker->MakeDigits(digits);
}

//____________________________________________________________________________
void AliITSQADataMakerSim::MakeDigits(TTree * digits)
{ 
  // Fill QA for RAW   
  if(fSubDetector == 0 || fSubDetector == 1) fSPDDataMaker->MakeDigits(digits);
  if(fSubDetector == 0 || fSubDetector == 2) fSDDDataMaker->MakeDigits(digits);
  if(fSubDetector == 0 || fSubDetector == 3) fSSDDataMaker->MakeDigits(digits);
}

//____________________________________________________________________________ 
void AliITSQADataMakerSim::InitSDigits()
{
  // Initialization for RECPOINTS
  if(fSubDetector == 0 || fSubDetector == 1) {
	AliDebug(1,"AliITSQADM:: SPD InitSDigits\n");
    fSPDDataMaker->InitSDigits();
  }
  if(fSubDetector == 0 || fSubDetector == 2) {
	AliDebug(1,"AliITSQADM:: SDD InitSDigits\n");
	fSDDDataMaker->InitSDigits();
  }
  if(fSubDetector == 0 || fSubDetector == 3) {
	AliDebug(1,"AliITSQADM:: SSD InitSDigits\n");
	fSSDDataMaker->InitSDigits();
  }
}

//____________________________________________________________________________ 
void AliITSQADataMakerSim::MakeSDigits(TClonesArray * sdigits)
{
  // Fill QA for recpoints
  if(fSubDetector == 0 || fSubDetector == 1) fSPDDataMaker->MakeSDigits(sdigits);
  if(fSubDetector == 0 || fSubDetector == 2) fSDDDataMaker->MakeSDigits(sdigits);
  if(fSubDetector == 0 || fSubDetector == 3) fSSDDataMaker->MakeSDigits(sdigits);
}

//____________________________________________________________________________ 
void AliITSQADataMakerSim::MakeSDigits(TTree * sdigits)
{
  // Fill QA for recpoints
  if(fSubDetector == 0 || fSubDetector == 1) fSPDDataMaker->MakeSDigits(sdigits);
  if(fSubDetector == 0 || fSubDetector == 2) fSDDDataMaker->MakeSDigits(sdigits);
  if(fSubDetector == 0 || fSubDetector == 3) fSSDDataMaker->MakeSDigits(sdigits);
}

//____________________________________________________________________________ 
void AliITSQADataMakerSim::InitHits()
{
  // Initialization for RECPOINTS
  if(fSubDetector == 0 || fSubDetector == 1) {
	AliDebug(1,"AliITSQADM:: SPD InitHits\n");
    fSPDDataMaker->InitHits();
  }
  if(fSubDetector == 0 || fSubDetector == 2) {
	AliDebug(1,"AliITSQADM:: SDD InitHits\n");
	fSDDDataMaker->InitHits();
  }
  if(fSubDetector == 0 || fSubDetector == 3) {
	AliDebug(1,"AliITSQADM:: SSD InitHits\n");
	fSSDDataMaker->InitHits();
  }
}

//____________________________________________________________________________ 
void AliITSQADataMakerSim::MakeHits(TClonesArray * hits)
{
  // Fill QA for recpoints
  if(fSubDetector == 0 || fSubDetector == 1) fSPDDataMaker->MakeHits(hits);
  if(fSubDetector == 0 || fSubDetector == 2) fSDDDataMaker->MakeHits(hits);
  if(fSubDetector == 0 || fSubDetector == 3) fSSDDataMaker->MakeHits(hits);
}

//____________________________________________________________________________ 
void AliITSQADataMakerSim::MakeHits(TTree * hits)
{
  // Fill QA for recpoints
  if(fSubDetector == 0 || fSubDetector == 1) fSPDDataMaker->MakeHits(hits);
  if(fSubDetector == 0 || fSubDetector == 2) fSDDDataMaker->MakeHits(hits);
  if(fSubDetector == 0 || fSubDetector == 3) fSSDDataMaker->MakeHits(hits);
}
