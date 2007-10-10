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
  AliQualAss * AliQualAss::fgQA        = 0x0 ;
  TFile      * AliQualAss::fgDataFile  = 0x0 ;   
  TString      AliQualAss::fgDataName  = "QA" ;   
  TString      AliQualAss::fgDetNames[]  = {"ITS", "TPC", "TRD", "TOF", "PHOS", "HMPID", "EMCAL", "MUON", "FMD",
					"ZDC", "PMD", "T0", "VZERO", "ACORDE", "HLT"} ;   
  TString      AliQualAss::fgTaskNames[]  = {"Raws", "Hits", "SDigits", "Digits", "RecPoints", "TrackSegments", "RecParticles", "ESDs"} ;   

//____________________________________________________________________________
AliQualAss::AliQualAss() : 
  TNamed("", ""), 
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
AliQualAss::AliQualAss(const DETECTORINDEX det) :
  TNamed("QA", "Quality Assurance status"), 
  fQA(new ULong_t[kNDET]), 
  fDet(det),
  fTask(kNULLTASK)
{
  // constructor to be used
  if (! CheckRange(det) ) {
    fDet = kNULLDET ; 
    return ;
  } 
  Int_t index ; 
  for (index = 0; index < kNDET; index++) 
    fQA[index] = 0 ; 
}
  
//_______________________________________________________________
AliQualAss::AliQualAss(const ALITASK tsk) :
  TNamed("QA", "Quality Assurance status"), 
  fQA(new ULong_t[kNDET]), 
  fDet(kNULLDET),
  fTask(tsk)
{
  // constructor to be used in the AliRoot module (SIM, REC, ESD or ANA)
  if (! CheckRange(tsk) ) {
    fTask = kNULLTASK ; 
    return ;
  } 
  Int_t index ; 
  for (index = 0; index < kNDET; index++) 
    fQA[index] = 0 ; 
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
  // check if detector is in given detector range: 0-kNDET

  Bool_t rv = ( det < 0 || det > kNDET )  ? kFALSE : kTRUE ;
  if (!rv)
    AliFatal(Form("Detector index %d is out of range: 0 <= index <= %d", det, kNDET)) ;
  return rv ;
}

//_______________________________________________________________
const Bool_t AliQualAss::CheckRange(ALITASK task) const
{ 
  // check if task is given taskk range: 0:kNTASK
  Bool_t rv = ( task < kRAW || task > kNTASK )  ? kFALSE : kTRUE ;
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
TFile * AliQualAss::GetQADMOutFile(const char * name, const Int_t run, const Int_t cycle) 
{
  // opens the file to store the detectors Quality Assurance Data Maker results
  char temp[100] ; 
  sprintf(temp, "%s.%s.%d.%d.root", name, fgDataName.Data(), run, cycle) ; 
  TString opt ; 
  if (! fgDataFile ) {     
    if  (gSystem->AccessPathName(temp))
      opt = "NEW" ;
    else 
      opt = "UPDATE" ; 
    fgDataFile = TFile::Open(temp, opt.Data()) ;
  } else {
   if ( (strcmp(temp, fgDataFile->GetName()) != 0) ) {
     if  (gSystem->AccessPathName(temp))
      opt = "NEW" ;
    else 
      opt = "UPDATE" ; 
    fgDataFile = TFile::Open(temp, opt.Data()) ;
   }
  }
  return fgDataFile ; 
} 

//_______________________________________________________________
const char * AliQualAss::GetDetName(Int_t det) 
{
	// returns the detector name corresponding to a given index (needed in a loop)

	if ( det >= 0 &&  det < kNDET) 
		return (fgDetNames[det]).Data() ; 
	else 
		return NULL ; 
}

//_______________________________________________________________
const char * AliQualAss::GetAliTaskName(ALITASK tsk)
{
  // returns the char name corresponding to module index
  TString tskName ;
  switch (tsk) {
  case kNULLTASK:
    break ; 
  case kRAW:
    tskName = "RAW" ;
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
    tsk = kNULLTASK ; 
    break ;
  }
  return tskName.Data() ;
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
AliQualAss * AliQualAss::Instance(const DETECTORINDEX det)
{
  // Get an instance of the singleton. The only authorized way to call the ctor
  
  if ( ! fgQA)
    fgQA = new AliQualAss(det) ;
  fgQA->Set(det) ;
  return fgQA ;
}

//_______________________________________________________________
AliQualAss * AliQualAss::Instance(const ALITASK tsk)
{
  // get an instance of the singleton.

  if ( ! fgQA)
    switch (tsk) {
    case kNULLTASK:
      break ;
	case kRAW:
      fgQA = new AliQualAss(tsk) ;
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
  case kRAW:
    offset+= 0 ;
    break ;
  case kSIM:
    offset+= 4 ;
    break ;
  case kREC:
    offset+= 8 ;
    break ;
  case kESD:
    offset+= 12 ;
    break ;
  case kANA:
    offset+= 16 ;
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
  ULong_t rawStatus = status & 0x0000f ;
  ULong_t simStatus = status & 0x000f0 ;
  ULong_t recStatus = status & 0x00f00 ;
  ULong_t esdStatus = status & 0x0f000 ;
  ULong_t anaStatus = status & 0xf0000 ;

  AliInfo(Form("QA Status for %s raw =0x%x, sim=0x%x, rec=0x%x, esd=0x%x, ana=0x%x\n", GetDetName(det).Data(), rawStatus, simStatus, recStatus, esdStatus, anaStatus )) ;
}

