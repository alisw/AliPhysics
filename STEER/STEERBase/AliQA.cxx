
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
// bit 0-3   : QA raised during simulation      (RAW)
// bit 4-7   : QA raised during simulation      (SIM)
// bit 8-11  : QA raised during reconstruction  (REC)
// bit 12-15 : QA raised during ESD checking    (ESD)
// bit 16-19 : QA raised during analysis        (ANA)
// Each of the 4 bits corresponds to a severity level of increasing importance
// from lower to higher bit (INFO, WARNING, ERROR, FATAL)
//
//*-- Yves Schutz CERN, July 2007 
//////////////////////////////////////////////////////////////////////////////


#include <cstdlib>
// --- ROOT system ---
#include <TClass.h>
#include <TFile.h>
#include <TSystem.h>
#include <TROOT.h>

// --- Standard library ---

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliQA.h"


ClassImp(AliQA)
AliQA    * AliQA::fgQA                   = 0x0 ;
TFile    * AliQA::fgQADataFile           = 0x0 ;   
TString    AliQA::fgQADataFileName       = "QA" ;  // will transform into Det.QA.run.cycle.root  
TFile    * AliQA::fgQARefFile            = 0x0 ;   
TString    AliQA::fgQARefDirName	       = "" ; 
TString    AliQA::fgQARefFileName        = "QA.root" ;
TFile    * AliQA::fgQAResultFile         = 0x0 ;  
TString    AliQA::fgQAResultDirName      = "" ;  
TString    AliQA::fgQAResultFileName     = "QA.root" ; 
TString    AliQA::fgDetNames[]           = {"ITS", "TPC", "TRD", "TOF", "PHOS", "HMPID", "EMCAL", "MUON", "FMD",
                                            "ZDC", "PMD", "T0", "VZERO", "ACORDE", "HLT", "Global", "CORR"
											#ifdef MFT_UPGRADE
	                                        , "MFT"
											#endif 
											};   
TString    AliQA::fgGRPPath              = "GRP/GRP/Data" ; 
TString    AliQA::fgRTNames[]            = {"UNKNOWN", "AUTO_TEST", "CALIBRATION", "CALIBRATION_PULSER", "CHANNEL_DELAY_TUNING", "COSMIC", 
																						"COSMICS", "DAQ_FO_UNIF_SCAN", "DAQ_GEN_DAC_SCAN", "DAQ_MEAN_TH_SCAN", "DAQ_MIN_TH_SCAN", 
																						"DAQ_NOISY_PIX_SCAN", "DAQ_PIX_DELAY_SCAN", "DAQ_UNIFORMITY_SCAN", "DCS_FO_UNIF_SCAN", 
																						"DCS_MEAN_TH_SCAN", "DCS_MIN_TH_SCAN", "DCS_PIX_DELAY_SCAN", "DCS_UNIFORMITY_SCAN", 
																						"DDL_TEST", "GAIN", "PEDESTAL", "INJECTOR",  "LASER", "MONTECARLO", "NOISE", "NOISY_PIX_SCAN", 
																						"PHYSICS", "PULSER", "STANDALONE", "STANDALONE_BC", "STANDALONE_CENTRAL", "STANDALONE_COSMIC", 
																						"STANDALONE_EMD", "STANDALONE_LASER", "STANDALONE_MB", "STANDALONE_PEDESTAL", 
																						"STANDALONE_SEMICENTRAL", "STANDALONE_PULSER" } ;   
TString       AliQA::fgTaskNames[]       = {"Raws", "Hits", "SDigits", "Digits", "RecPoints", "TrackSegments", "RecParticles", "ESDs"} ;   
const TString AliQA::fgkLabLocalFile     = "file://"  ; 
const TString AliQA::fgkLabLocalOCDB     = "local://" ;  
const TString AliQA::fgkLabAliEnOCDB     = "alien://" ;  
const TString AliQA::fgkRefFileName      = "QA.root" ; 
const TString AliQA::fgkQAName           = "QA"  ; 
const TString AliQA::fgkQACorrNtName     = "CorrQA" ;  
const TString AliQA::fgkRefOCDBDirName   = "Ref"  ; 
TString AliQA::fgkRefDataDirName	       = ""  ; 
const TString AliQA::fgkQARefOCDBDefault = "alien://folder=/alice/QA/20"  ; 
const TString AliQA::fgkExpert           = "Expert" ; 
const UInt_t  AliQA::fgkExpertBit        = 16 ; 
const UInt_t  AliQA::fgkQABit            = 17 ; 

//____________________________________________________________________________
AliQA::AliQA() : 
  TNamed("", ""), 
  fNdet(kNDET), 
  fQA(new ULong_t[fNdet]), 
  fDet(kNULLDET),
  fTask(kNULLTASK)
  	
{
  // default constructor
  // beware singleton: not to be used
  
  for (Int_t index = 0 ; index < fNdet ; index++) 
	fQA[index] = 0 ; 
}

//____________________________________________________________________________
AliQA::AliQA(const AliQA& qa) :
  TNamed(qa),
  fNdet(qa.fNdet), 
  fQA(new ULong_t[qa.fNdet]), 
  fDet(qa.fDet),
  fTask(qa.fTask)
{ 
  // cpy ctor
  for (Int_t index = 0 ; index < fNdet ; index++) 
	fQA[index] = qa.fQA[index] ; 
}

//_____________________________________________________________________________
AliQA& AliQA::operator = (const AliQA& qa)
{
// assignment operator

  this->~AliQA();
  new(this) AliQA(qa);
  return *this;
}

//_______________________________________________________________
AliQA::AliQA(const DETECTORINDEX_t det) :
  TNamed("QA", "Quality Assurance status"),
  fNdet(kNDET),  
  fQA(new ULong_t[fNdet]), 
  fDet(det),
  fTask(kNULLTASK) 
{
  // constructor to be used
  if (! CheckRange(det) ) {
    fDet = kNULLDET ; 
    return ;
  } 
  Int_t index ; 
  for (index = 0; index < fNdet; index++) 
    fQA[index] = 0 ; 
}
  
//_______________________________________________________________
AliQA::AliQA(const ALITASK_t tsk) :
  TNamed("QA", "Quality Assurance status"), 
  fNdet(kNDET),
  fQA(new ULong_t[fNdet]), 
  fDet(kNULLDET),
  fTask(tsk)
{
  // constructor to be used in the AliRoot module (SIM, REC, ESD or ANA)
  if (! CheckRange(tsk) ) {
    fTask = kNULLTASK ; 
    return ;
  } 
  Int_t index ; 
  for (index = 0; index < fNdet; index++) 
    fQA[index] = 0 ; 
}

//____________________________________________________________________________
AliQA::~AliQA() 
{
  // dtor  
  delete[] fQA ;
}

//_______________________________________________________________
void AliQA::Close() 
{
	// close the open files
	if (fgQADataFile) 
		if (fgQADataFile->IsOpen())
			fgQADataFile->Close() ; 
	if (fgQAResultFile) 
		if (fgQAResultFile->IsOpen()) 
			fgQAResultFile->Close() ;
	if (fgQARefFile)
		if (fgQARefFile->IsOpen())
			fgQARefFile->Close() ; 
} 

//_______________________________________________________________
Bool_t AliQA::CheckFatal() const
{
  // check if any FATAL status is set
  Bool_t rv = kFALSE ;
  Int_t index ;
  for (index = 0; index < kNDET ; index++)
    rv = rv || IsSet(DETECTORINDEX_t(index), fTask, kFATAL) ;
  return rv ;
}

//_______________________________________________________________
Bool_t AliQA::CheckRange(DETECTORINDEX_t det) const
{ 
  // check if detector is in given detector range: 0-kNDET

  Bool_t rv = ( det < 0 || det > kNDET )  ? kFALSE : kTRUE ;
  if (!rv)
    AliFatal(Form("Detector index %d is out of range: 0 <= index <= %d", det, kNDET)) ;
  return rv ;
}

//_______________________________________________________________
Bool_t AliQA::CheckRange(ALITASK_t task) const
{ 
  // check if task is given taskk range: 0:kNTASK
  Bool_t rv = ( task < kRAW || task > kNTASK )  ? kFALSE : kTRUE ;
  if (!rv)
    AliFatal(Form("Module index %d is out of range: 0 <= index <= %d", task, kNTASK)) ;
  return rv ;
}

//_______________________________________________________________
Bool_t AliQA::CheckRange(QABIT_t bit) const
{ 
  // check if bit is in given bit range: 0-kNBit

  Bool_t rv = ( bit < 0 || bit > kNBIT )  ? kFALSE : kTRUE ;
  if (!rv)
    AliFatal(Form("Status bit %d is out of range: 0 <= bit <= %d", bit, kNBIT)) ;
  return rv ;
}



//_______________________________________________________________
char * AliQA::GetAliTaskName(ALITASK_t tsk)
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
	return Form("%s", tskName.Data()) ;
}

//_______________________________________________________________
char * AliQA::GetBitName(QABIT_t bit) const
{
	// returns the char name corresponding to bit 
	TString bitName ;
	switch (bit) {
		case kNULLBit:
			break ; 
		case kINFO:
			bitName = "INFO" ;
			break ;  
		case kWARNING:
			bitName = "WARNING" ;
			break ;
		case kERROR:
			bitName = "ERROR" ;
			break ;
		case kFATAL:
			bitName = "FATAL" ;
			break ;
		default:
			bit = kNULLBit ; 
			break ;
	}
	return Form("%s", bitName.Data()) ; 
}

//_______________________________________________________________
AliQA::DETECTORINDEX_t AliQA::GetDetIndex(const char * name) 
{
	// returns the detector index corresponding to a given name
	TString sname(name) ; 
	DETECTORINDEX_t rv = kNULLDET ; 
	for (Int_t det = 0; det < kNDET ; det++) {
		if ( GetDetName(det) == sname ) {
			rv = DETECTORINDEX_t(det) ; 
			break ; 
		}
	}
	return rv ; 		
}

//_______________________________________________________________
const char * AliQA::GetDetName(Int_t det) 
{
	// returns the detector name corresponding to a given index (needed in a loop)
	
	if ( det >= 0 &&  det < kNDET) 
		return (fgDetNames[det]).Data() ; 
	else 
		return NULL ; 
}

//_______________________________________________________________
TFile * AliQA::GetQADataFile(const char * name, const Int_t run) 
{
  // opens the file to store the detectors Quality Assurance Data Maker results
	const char * temp = Form("%s.%s.%d.root", name, fgQADataFileName.Data(), run) ; 
	TString opt ; 
	if (! fgQADataFile ) {     
		if  (gSystem->AccessPathName(temp))
			opt = "NEW" ;
		else 
			opt = "UPDATE" ; 
		fgQADataFile = TFile::Open(temp, opt.Data()) ;
	} else {
		if ( strcmp(temp, fgQADataFile->GetName()) != 0 ) {
			fgQADataFile = dynamic_cast<TFile *>(gROOT->FindObject(temp)) ; 
			if ( !fgQADataFile ) {
				if  (gSystem->AccessPathName(temp))
					opt = "NEW" ;
				else 
					opt = "UPDATE" ; 
				fgQADataFile = TFile::Open(temp, opt.Data()) ;
			}
		}
  }
	return fgQADataFile ;
} 

//_____________________________________________________________________________
TFile * AliQA::GetQADataFile(const char * fileName)
{
  // Open if necessary the Data file and return its pointer

  if (!fgQADataFile) 
	if (!fileName) 
		fileName = AliQA::GetQADataFileName() ; 
	if  (!gSystem->AccessPathName(fileName)) {
		fgQADataFile =  TFile::Open(fileName) ;
	} else {
		AliFatalClass(Form("File %s not found", fileName)) ;
	}
  return fgQADataFile ; 
}

//_______________________________________________________________
TFile * AliQA::GetQAResultFile() 
{
  // opens the file to store the  Quality Assurance Data Checker results
	if (fgQAResultFile) 
		fgQAResultFile->Close() ; 
	fgQAResultFile = 0x0 ; 
//	if (!fgQAResultFile) { 
		TString dirName(fgQAResultDirName) ; 
		if ( dirName.Contains(fgkLabLocalFile)) 
			dirName.ReplaceAll(fgkLabLocalFile, "") ;
		TString fileName(dirName + fgQAResultFileName) ; 
		TString opt("") ; 
		if ( !gSystem->AccessPathName(fileName) )
			opt = "UPDATE" ; 
		else { 
			if ( gSystem->AccessPathName(dirName) )
				gSystem->mkdir(dirName) ; 
			opt = "NEW" ; 
		}
		fgQAResultFile = TFile::Open(fileName, opt) ;   
//	}
	
	return fgQAResultFile ;
}

//_______________________________________________________________
const TString AliQA::GetRunTypeName(RUNTYPE_t rt) 
{ 
  // returns the name of a run type
	TString rv("Invalid Run Type") ; 	
	if ( rt == kNULLTYPE ) {
		rv = "Known RUN_TYPE are: \n" ; 
		for (Int_t index = 0 ; index < kNTYPE; index++) {
			rv += Form("%2d -- %s\n", index, fgRTNames[index].Data()) ;  
		}
		AliErrorClass(Form("%s", rv.Data())) ; 
		return "" ; 
	}
	else {
		if ( rt > kNULLTYPE && rt < kNTYPE )
			rv = fgRTNames[rt] ; 
	}
	return rv ; 
}

//_______________________________________________________________
AliQA::TASKINDEX_t AliQA::GetTaskIndex(const char * name) 
{
	// returns the detector index corresponding to a given name
	TString sname(name) ; 
	TASKINDEX_t rv = kNULLTASKINDEX ; 
	for (Int_t tsk = 0; tsk < kNTASKINDEX ; tsk++) {
		if ( GetTaskName(tsk) == sname ) {
			rv = TASKINDEX_t(tsk) ; 
			break ; 
		}
	}
	return rv ; 		
}

//_______________________________________________________________
Bool_t AliQA::IsSet(DETECTORINDEX_t det, ALITASK_t tsk, QABIT_t bit) const
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
Bool_t AliQA::IsSetAny(DETECTORINDEX_t det, ALITASK_t tsk) const
{
  // Checks is the requested bit is set
	
  CheckRange(det) ; 
  CheckRange(tsk) ;
	
  ULong_t offset = Offset(tsk) ;
  ULong_t status = GetStatus(det) ;
	UShort_t st = 0 ; 
	for ( Int_t bit = 0 ; bit < kNBIT ; bit++) {
		offset+= bit ;
		st += (status & 1 << offset) != 0 ;		
	}
	if ( st == 0 ) 
		return kFALSE ; 
	else 
		return kTRUE ;
}
//_______________________________________________________________
Bool_t AliQA::IsSetAny(DETECTORINDEX_t det) const
{
  // Checks is the requested bit is set
	
  CheckRange(det) ; 
	
	ULong_t status = GetStatus(det) ;
	UShort_t st = 0 ; 
	for ( Int_t tsk = 0 ; tsk < kNTASK ; tsk++) {
		ULong_t offset = Offset(ALITASK_t(tsk)) ;
		for ( Int_t bit = 0 ; bit < kNBIT ; bit++) {
			offset+= bit ;
			st += (status & 1 << offset) != 0 ;		
		}
	}
	if ( st == 0 ) 
		return kFALSE ; 
	else 
		return kTRUE ;
}

//_______________________________________________________________
AliQA * AliQA::Instance()
{
  // Get an instance of the singleton.
  // Object must have been instantiated with Instance(ALITASK) first

  return fgQA ;
}

//_______________________________________________________________
AliQA * AliQA::Instance(const DETECTORINDEX_t det)
{
  // Get an instance of the singleton. The only authorized way to call the ctor
  
  if ( ! fgQA) {
    TFile * f = GetQAResultFile() ; 
	fgQA = dynamic_cast<AliQA *>(f->Get("QA")) ; 
    if ( ! fgQA ) 
		fgQA = new AliQA(det) ;
  }		
  fgQA->Set(det) ;
  return fgQA ;
}

//_______________________________________________________________
AliQA * AliQA::Instance(const ALITASK_t tsk)
{
  // get an instance of the singleton.

  if ( ! fgQA)
    switch (tsk) {
    case kNULLTASK:
      break ;
	case kRAW:
      fgQA = new AliQA(tsk) ;
      break ;
	case kSIM:
      fgQA = new AliQA(tsk) ;
      break ;
    case kREC:
      AliInfoClass("fgQA = gAlice->GetQA()") ;
      break ;
    case kESD:
      AliInfoClass("fgQA = dynamic_cast<AliQA *> (esdFile->Get(\"QA\")") ;
      break ;
    case kANA:
      AliInfoClass("fgQA = dynamic_cast<AliQA *> (esdFile->Get(\"QA\")") ;
      break ;
    case kNTASK:
      break ;
    }
  if (fgQA) 
    fgQA->Set(tsk) ;
  return fgQA ;
}

//_______________________________________________________________
AliQA *  AliQA::Instance(const TASKINDEX_t tsk) 
{
	// get an instance of the singleton.
	
	ALITASK_t index = kNULLTASK ; 

	if ( tsk == kRAWS )
		index = kRAW ;
	else if (tsk < kDIGITS)
		index = kSIM ;
	else if (tsk < kRECPARTICLES)
		index = kREC ; 
	else if (tsk == kESDS) 
		index = kESD ; 

	return Instance(index) ; 
}

//_______________________________________________________________
Long64_t AliQA::Merge(const TCollection * list) {
  // Merge the QA resuls in the list into this single AliQA object
  
  Long64_t nmerge=0;
  for (Int_t det = 0 ; det < kNDET ; det++) {
    Set(DETECTORINDEX_t(det)) ; 
    for (Int_t task = 0 ; task < kNTASK ; task++) {
      Set(ALITASK_t(task)) ; 
      for (Int_t bit = 0 ; bit < kNBIT ; bit++) {
	TIter next(list) ;
	AliQA * qa ; 
	while ( (qa = (AliQA*)next() ) ) {
	  ++nmerge;
	  qa->IsSet(DETECTORINDEX_t(det), ALITASK_t(task), QABIT_t(bit)) ;
	  Set(QABIT_t(bit)) ; 
	} // qa list
      } // bit
    } // task
  } // detector
  return nmerge;
}

//_______________________________________________________________
ULong_t AliQA::Offset(ALITASK_t tsk) const
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
void AliQA::Set(QABIT_t bit)
{
  // Set the status bit of the current detector in the current module
  
  SetStatusBit(fDet, fTask, bit) ;
}

//_____________________________________________________________________________
void AliQA::SetQARefStorage(const char * name)
{
	// Set the root directory where the QA reference data are stored

	fgQARefDirName = name ; 
	if ( fgQARefDirName.Contains(fgkLabLocalFile) )
		fgQARefFileName =  fgkRefFileName ; 
	else if ( fgQARefDirName.Contains(fgkLabLocalOCDB) )
		fgQARefFileName =  fgkQAName ; 
	else if ( fgQARefDirName.Contains(fgkLabAliEnOCDB) )
		fgQARefFileName =  fgkQAName ; 

  else {
	  AliErrorClass(Form("%s is an invalid storage definition\n", name)) ; 
	  fgQARefDirName  = "" ; 
	  fgQARefFileName = "" ; 
  }	
	TString tmp(fgQARefDirName) ; // + fgQARefFileName) ;
	AliInfoClass(Form("AliQA::SetQARefDir: QA references are in  %s\n", tmp.Data()) ) ;
}

//_____________________________________________________________________________
void AliQA::SetQARefDataDirName(const char * name) 
{
  // Set the lower level directory name where reference data are found
  TString test(name) ; 
  RUNTYPE_t rt = kNULLTYPE ; 
  for (Int_t index = 0; index < kNTYPE; index++) {
    if (test.CompareTo(fgRTNames[index]) == 0) {
      rt = (RUNTYPE_t) index ; 
      break ; 
    } 
	}
	
	if (rt == kNULLTYPE) {
      AliWarningClass(Form("AliQA::SetQARefDataDirName: %s is an unknown RUN TYPE name\n", name)) ; 
      return ; 
	}
 
	SetQARefDataDirName(rt) ; 
}

//_____________________________________________________________________________
void AliQA::SetQAResultDirName(const char * name)
{
  // Set the root directory where to store the QA status object

  fgQAResultDirName.Prepend(name) ; 
  AliInfoClass(Form("AliQA::SetQAResultDirName: QA results are in  %s\n", fgQAResultDirName.Data())) ;
  if ( fgQAResultDirName.Contains(fgkLabLocalFile)) 
    fgQAResultDirName.ReplaceAll(fgkLabLocalFile, "") ;
  fgQAResultFileName.Prepend(fgQAResultDirName) ;
}

//_______________________________________________________________
void AliQA::SetStatusBit(DETECTORINDEX_t det, ALITASK_t tsk, QABIT_t bit)
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
void AliQA::ShowAll() const
{
  // dispplay the QA status word
  Int_t index ;
  for (index = 0 ; index < kNDET ; index++) {
		for (Int_t tsk = kRAW ; tsk < kNTASK ; tsk++) {
			ShowStatus(DETECTORINDEX_t(index), ALITASK_t(tsk)) ;
		}
	}
}

//_______________________________________________________________
void AliQA::ShowStatus(DETECTORINDEX_t det, ALITASK_t tsk) const
{
	// Prints the full QA status of a given detector
	CheckRange(det) ;
	ULong_t status = GetStatus(det) ;
	ULong_t tskStatus[kNTASK] ; 
	tskStatus[kRAW] = status & 0x0000f ;
	tskStatus[kSIM] = status & 0x000f0 ;
	tskStatus[kREC] = status & 0x00f00 ;
	tskStatus[kESD] = status & 0x0f000 ;
	tskStatus[kANA] = status & 0xf0000 ;

	AliInfo(Form("====> QA Status for %8s raw =0x%lx, sim=0x%lx, rec=0x%lx, esd=0x%lx, ana=0x%lx", GetDetName(det).Data(), 
				 tskStatus[kRAW], tskStatus[kSIM], tskStatus[kREC], tskStatus[kESD], tskStatus[kANA] )) ;
	if (tsk == kNULLTASK) {
		for (Int_t itsk = kRAW ; itsk < kNTASK ; itsk++) {
			ShowASCIIStatus(det, ALITASK_t(itsk), tskStatus[itsk]) ; 
		} 
	} else {
			ShowASCIIStatus(det, tsk, tskStatus[tsk]) ; 
	}
}

//_______________________________________________________________
void AliQA::ShowASCIIStatus(DETECTORINDEX_t det, ALITASK_t tsk, const ULong_t status) const 
{
	// print the QA status in human readable format
	TString text; 
	for (Int_t bit = kINFO ; bit < kNBIT ; bit++) {
		if (IsSet(det, tsk, QABIT_t(bit))) {
			text = GetBitName(QABIT_t(bit)) ; 
			text += " " ; 
		}
	}
	if (! text.IsNull())
		AliWarningClass(Form("           %8s %4s 0x%4lx, Problem signalled: %8s \n", GetDetName(det).Data(), GetAliTaskName(tsk), status, text.Data())) ; 
}

//_______________________________________________________________
void AliQA::UnSet(QABIT_t bit)
{
	// UnSet the status bit of the current detector in the current module
	
	UnSetStatusBit(fDet, fTask, bit) ;
}

//_______________________________________________________________
void AliQA::UnSetStatusBit(DETECTORINDEX_t det, ALITASK_t tsk, QABIT_t bit)
{
	// UnSet the status bit for a given detector and a given task
	
	CheckRange(det) ;
	CheckRange(tsk) ;
	CheckRange(bit) ;
	
	ULong_t offset = Offset(tsk) ;
	ULong_t status = GetStatus(det) ;
	offset+= bit ;
	status = status & 0 << offset ;
	SetStatus(det, status) ;
}
