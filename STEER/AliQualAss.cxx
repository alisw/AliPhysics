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

//////////////////////////////////////////////////////////////////////////////
//
// Quality Assurance Object//_________________________________________________________________________
// Quality Assurance object. The QA status is held in one word per detector,
// each bit corresponds to a different status.
// bit 0-3  : QA raised during simulation      (SIM)
// bit 4-7  : QA raised during reconstruction  (REC)
// bit 8-11 : QA raised during ESD checking    (ESD)
// bit 12-15: QA raised during analysis        (ANA)
// Each of the 4 bits corresponds to a severity level of increasing importance
// from lower to higher bit (INFO, WARNING, ERROR, FATAL)
//
//*-- Yves Schutz CERN, July 2007 
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system ---
#include <TFile.h>
#include <TSystem.h>

// --- Standard library ---

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliQualAss.h"


ClassImp(AliQualAss)

  AliQualAss * AliQualAss::fgQA          = 0x0 ;
  TFile      * AliQualAss::fgOutput      = 0x0 ;   
  TString      AliQualAss::fgOutputName  = "QA.root" ;   

//____________________________________________________________________________
AliQualAss::AliQualAss() : 
  TNamed("", ""), 
  fNdet(12), 
  fQA(0x0), 
  fDet(kNULLDET),
  fTask(kNULLTASK)
{
  // default constructor
  // beware singleton: not to be used
}

//____________________________________________________________________________
AliQualAss::AliQualAss(const AliQualAss& qa) :
  TNamed(qa),
  fNdet(qa.fNdet),
  fQA(qa.fQA), 
  fDet(qa.fDet),
  fTask(qa.fTask)
{ 
  // cpy ctor
}

//_____________________________________________________________________________
AliQualAss& AliQualAss::operator = (const AliQualAss& qa)
{
// assignment operator

  this->~AliQualAss();
  new(this) AliQualAss(qa);
  return *this;
}

//_______________________________________________________________
AliQualAss::AliQualAss(ALITASK tsk) :
  TNamed("QA", "Quality Assurance status"), 
  fNdet(12), 
  fQA(0x0), 
  fDet(kNULLDET),
  fTask(tsk)
{
  // constructor to be used in the AliRoot module (SIM, REC, ESD or ANA)
  if (! CheckRange(tsk) ) {
    fTask = kNULLTASK ; 
    return ;
  } else {
    fQA = new ULong_t[fNdet] ;
    Int_t index ;
    for ( index = 0 ; index <= fNdet ; index++)
      ResetStatus(DETECTORINDEX(index)) ;
  }
}

//____________________________________________________________________________
AliQualAss::~AliQualAss() 
{
  // dtor  
  delete[] fQA ;
}

//_______________________________________________________________
const Bool_t AliQualAss::CheckRange(DETECTORINDEX det) const
{ 
  // check if detector is in given detector range: 0-fNdet

  Bool_t rv = ( det < 0 || det > fNdet )  ? kFALSE : kTRUE ;
  if (!rv)
    AliFatal(Form("Detector index %d is out of range: 0 <= index <= %d", det, kNDET)) ;
  return rv ;
}

//_______________________________________________________________
const Bool_t AliQualAss::CheckRange(ALITASK task) const
{ 
  // check if task is given taskk range: 0:kNTASK
  Bool_t rv = ( task < kSIM || task > kNTASK )  ? kFALSE : kTRUE ;
  if (!rv)
    AliFatal(Form("Module index %d is out of range: 0 <= index <= %d", task, kNTASK)) ;
  return rv ;
}

//_______________________________________________________________
const Bool_t AliQualAss::CheckRange(QABIT bit) const
{ 
  // check if bit is in given bit range: 0-kNBit

  Bool_t rv = ( bit < 0 || bit > kNBIT )  ? kFALSE : kTRUE ;
  if (!rv)
    AliFatal(Form("Status bit %d is out of range: 0 <= bit <= %d", bit, kNBIT)) ;
  return rv ;
}

//_______________________________________________________________
const char * AliQualAss::GetDetectorName(DETECTORINDEX det) const
{
  // returns the char name corresponding to detector index

  char * detName = "";
  switch (det) {
  case kNULLDET:
    break ; 
  case kITS:
    detName = "ITS" ;
    break ;
  case kTPC:
    detName = "TPC" ;
    break ;
  case kTRD:
    detName = "TRD" ;
    break ;
  case kTOF:
    detName = "TOF" ;
    break ;
  case kPHOS:
    detName = "PHOS" ;
    break ;
  case kHMPID:
    detName = "HMPID" ;
    break ;
  case kEMCAL:
    detName = "EMCAL" ;
    break ;
  case kMUON:
    detName = "MUON" ;
    break ;
  case kFMD:
    detName = "FMD" ;
    break ;
  case kZDC:
    detName = "ZDC" ;
    break ;
  case kPMD:
    detName = "PMD" ;
    break ;
  case kT0:
    detName = "TO" ;
    break ;
  case kVZERO:
    detName = "VZERO" ;
    break ;
  case kACORDE:
    detName = "ACORDE" ;
    break ;
  case kHLT:
    detName = "HLT" ;
    break ;
  default:
    AliError(Form("%d is not a valid detector index %d <= index <= %d\n", det, 0, kNDET-1)) ;
    break ;
  }
  return detName ;
}

//_______________________________________________________________
TFile * AliQualAss::GetQADMOutFile() 
{
  // opens the file to store the detectors Quality Assurance Data Maker results

  if (! fgOutput ) {     
	 char opt[6] ; 
     if  (gSystem->AccessPathName(fgOutputName.Data()))
	  sprintf(opt, "%s", "NEW") ;
     else 
      sprintf(opt, "%s", "UPDATE") ; 
    
     fgOutput = TFile::Open(fgOutputName.Data(), opt) ;
  }
  return fgOutput ; 
} 

//_______________________________________________________________
const char * AliQualAss::GetTaskName(ALITASK tsk) const
{
  // returns the char name corresponding to module index
  char * tskName = "" ;
  switch (tsk) {
  case kNULLTASK:
    break ; 
  case kSIM:
    tskName = "SIM" ;
    break ;
  case kREC:
    tskName = "REC" ;
    break ;
  case kESD:
    tskName = "ESD" ;
    break ;
  case kANA:
    tskName = "ANA" ;
    break ;
  default:
    AliError(Form("%d is not a valid module index %d <= index <= %d\n", tsk, 0, kNTASK-1)) ;
    break ;
  }
  return tskName ;
}

//_______________________________________________________________
const Bool_t AliQualAss::CheckFatal() const
{
  // check if any FATAL status is set
  Bool_t rv = kFALSE ;
  Int_t index ;
  for (index = 0; index < kNDET ; index++)
    rv = rv || IsSet(DETECTORINDEX(index), fTask, kFATAL) ;
  return rv ;
}

//_______________________________________________________________
const Bool_t AliQualAss::IsSet(DETECTORINDEX det, ALITASK tsk, QABIT bit) const
{
  // Checks is the requested bit is set

  CheckRange(det) ; 
  CheckRange(tsk) ;
  CheckRange(bit) ;

  ULong_t offset = Offset(tsk) ;
  ULong_t status = GetStatus(det) ;
  offset+= bit ;
  status = (status & 1 << offset) != 0 ;
  return status ;
}

//_______________________________________________________________
AliQualAss * AliQualAss::Instance()
{
  // Get an instance of the singleton.
  // Object must have been instantiated with Instance(ALITASK) first

  return fgQA ;
}

//_______________________________________________________________
AliQualAss * AliQualAss::Instance(DETECTORINDEX det)
{
  // Get an instance of the singleton. The only authorized way to call the ctor
  
  fgQA->Set(det) ;
  return fgQA ;
}

//_______________________________________________________________
AliQualAss * AliQualAss::Instance(ALITASK tsk)
{
  // get an instance of the singleton.

  if ( ! fgQA)
    switch (tsk) {
    case kNULLTASK:
      break ;
    case kSIM:
      fgQA = new AliQualAss(tsk) ;
      break ;
    case kREC:
      printf("fgQA = gAlice->GetQA()") ;
      break ;
    case kESD:
      printf("fgQA = dynamic_cast<AliQA *> (esdFile->Get(\"QA\")") ;
      break ;
    case kANA:
      printf("fgQA = dynamic_cast<AliQA *> (esdFile->Get(\"QA\")") ;
      break ;
    case kNTASK:
      break ;
    }
  if (fgQA) 
    fgQA->Set(tsk) ;
  return fgQA ;
}

//_______________________________________________________________
const ULong_t AliQualAss::Offset(ALITASK tsk) const
{
  // Calculates the bit offset for a given module (SIM, REC, ESD, ANA)

  CheckRange(tsk) ; 

  ULong_t offset = 0 ;
  switch (tsk) {
  case kNULLTASK:
    break ;
  case kSIM:
    offset+= 0 ;
    break ;
  case kREC:
    offset+= 4 ;
    break ;
  case kESD:
    offset+= 8 ;
    break ;
  case kANA:
    offset+= 12 ;
    break ;
  case kNTASK:
    break ;
  }

  return offset ;
}

//_______________________________________________________________
void AliQualAss::Set(QABIT bit)
{
  // Set the status bit of the current detector in the current module
  
  SetStatusBit(fDet, fTask, bit) ;
}

//_______________________________________________________________
void AliQualAss::SetStatusBit(DETECTORINDEX det, ALITASK tsk, QABIT bit)
{
 // Set the status bit for a given detector and a given task

  CheckRange(det) ;
  CheckRange(tsk) ;
  CheckRange(bit) ;

  ULong_t offset = Offset(tsk) ;
  ULong_t status = GetStatus(det) ;
  offset+= bit ;
  status = status | 1 << offset ;
  SetStatus(det, status) ;
}

//_______________________________________________________________
void AliQualAss::ShowAll() const
{
  // dispplay the QA status word
  Int_t index ;
  for (index = 0 ; index < kNDET ; index++)
    ShowStatus(DETECTORINDEX(index)) ;
}

//_______________________________________________________________
void AliQualAss::ShowStatus(DETECTORINDEX det) const
{
  // Prints the full QA status of a given detector
  CheckRange(det) ;
  ULong_t status = GetStatus(det) ;
  ULong_t simStatus = status & 0x000f ;
  ULong_t recStatus = status & 0x00f0 ;
  ULong_t esdStatus = status & 0x0f00 ;
  ULong_t anaStatus = status & 0xf000 ;

  AliInfo(Form("QA Status for %s sim=0x%x, rec=0x%x, esd=0x%x, ana=0x%x\n", GetDetectorName(det), simStatus, recStatus, esdStatus, anaStatus )) ;
}

