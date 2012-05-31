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
/* $Id: AliQAv1.cxx 31503 2009-03-16 11:01:16Z schutz $ */

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
#include <TH1.h>
#include <TSystem.h>
#include <TROOT.h>

// --- Standard library ---

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliQAv1.h"


ClassImp(AliQAv1)
AliQAv1  *     AliQAv1::fgQA                   = 0x0 ;
TFile    *     AliQAv1::fgQADataFile           = 0x0 ;   
TString        AliQAv1::fgQADataFileName       = "QA" ;  // will transform into Det.QA.run.root  
TFile    *     AliQAv1::fgQARefFile            = 0x0 ;   
TString        AliQAv1::fgQARefDirName	       = "" ; 
TString        AliQAv1::fgQARefFileName        = "QA.root" ;
TFile    *     AliQAv1::fgQAResultFile         = 0x0 ;  
TString        AliQAv1::fgQAResultDirName      = "" ;  
TString        AliQAv1::fgQAResultFileName     = "QA.root" ; 
TString        AliQAv1::fgDetNames[]           = {"ITS", "TPC", "TRD", "TOF", "PHOS", "HMPID", "EMCAL", "MUON", "FMD",
                                                  "ZDC", "PMD", "T0", "VZERO", "ACORDE", "HLT", "Global", "CORR"
												  #ifdef MFT_UPGRADE
												  , "MFT" 
												  #endif 
												  } ;   
TString        AliQAv1::fgGRPPath              = "GRP/GRP/Data" ; 
TString        AliQAv1::fgTaskNames[]          = {"Raws", "Hits", "SDigits", "Digits", "DigitsR", "RecPoints", "TrackSegments", "RecParticles", "ESDs"} ;   
TString        AliQAv1::fgModeNames[]          = {"", "Sim", "Rec", "QA"} ;   
const TString  AliQAv1::fgkLabLocalFile        = "file://"  ; 
const TString  AliQAv1::fgkLabLocalOCDB        = "local://" ;  
const TString  AliQAv1::fgkLabAliEnOCDB        = "alien://" ;  
const TString  AliQAv1::fgkRefFileName         = "QA.root" ; 
const TString  AliQAv1::fgkQAName              = "QA"  ; 
const TString  AliQAv1::fgkQACorrNtName        = "CorrQA" ;  
const TString  AliQAv1::fgkRefOCDBDirName      = "QA"  ; 
TString AliQAv1::fgRefDataDirName	             = ""  ; 
const TString  AliQAv1::fgkQARefOCDBDefault    = "alien://folder=/alice/QA/20"  ; 
const TString  AliQAv1::fgkExpert              = "Expert" ; 
//
const UInt_t   AliQAv1::fgkQAClonedBit         = BIT(14); // RS: bits set to the histos or arrays of cloned histos
const UInt_t   AliQAv1::fgkForbidCloningBit    = BIT(15); // to flag the state (cloned or not) of the 
const UInt_t   AliQAv1::fgkOrigHistoKeptBit    = BIT(16); // object at the given slot of fXXXQAList[]
//
const UInt_t   AliQAv1::fgkExpertBit           = BIT(19);  
const UInt_t   AliQAv1::fgkQABit               = BIT(20) ; 
const UInt_t   AliQAv1::fgkImageBit            = BIT(21) ; 
const Int_t    AliQAv1::fgkQADebugLevel        = 99 ; 
const TString  AliQAv1::fgkImageFileName       = "QAImage" ; 
const TString  AliQAv1::fImageFileFormat       = "ps" ; 
const UShort_t AliQAv1::fgkMaxQAObjects        = 10000 ; 

//____________________________________________________________________________
AliQAv1::AliQAv1() : 
  TNamed("", ""), 
  fNdet(kNDET), 
  fNEventSpecies(AliRecoParam::kNSpecies), 
  fLengthQA(fNdet*fNEventSpecies),
  fQA(new ULong_t[fLengthQA]), 
  fDet(kNULLDET),
  fTask(kNULLTASK), 
  fEventSpecie(AliRecoParam::kDefault), 
  fEventSpecies(new Bool_t[fNEventSpecies])
{
  // default constructor
  memset(fQA,0,fLengthQA*sizeof(ULong_t));
  memset(fEventSpecies,kFALSE,fNEventSpecies*sizeof(Bool_t));
}

//____________________________________________________________________________
AliQAv1::AliQAv1(const AliQAv1& qa) :
  TNamed(qa),
  fNdet(qa.fNdet), 
  fNEventSpecies(qa.fNEventSpecies), 
  fLengthQA(qa.fLengthQA),
  fQA(new ULong_t[fLengthQA]), 
  fDet(qa.fDet),
  fTask(qa.fTask), 
  fEventSpecie(qa.fEventSpecie), 
  fEventSpecies(new Bool_t[fNEventSpecies])
{ 
  // cpy ctor
  memcpy(fQA,qa.fQA,fLengthQA*sizeof(ULong_t));
  memcpy(fEventSpecies,qa.fEventSpecies,fNEventSpecies*sizeof(Bool_t));
}

//_____________________________________________________________________________
AliQAv1& AliQAv1::operator = (const AliQAv1& qa)
{
  // assignment operator
  if(&qa != this) {
    TNamed::operator=(qa);
    fNdet          = qa.fNdet;
    fNEventSpecies = qa.fNEventSpecies; 
    fLengthQA      = qa.fLengthQA;

    if(fQA) delete [] fQA;
    fQA = new ULong_t[fLengthQA];
    memcpy(fQA,qa.fQA,fLengthQA*sizeof(ULong_t));

    fDet = qa.fDet;
    fTask = qa.fTask;
    fEventSpecie = qa.fEventSpecie; 
    if(fEventSpecies) delete [] fEventSpecies;
    fEventSpecies = new Bool_t[fNEventSpecies];
    memcpy(fEventSpecies,qa.fEventSpecies,fNEventSpecies*sizeof(Bool_t));
  }  
  return *this;
}

//_______________________________________________________________
AliQAv1::AliQAv1(const Int_t qalength, ULong_t * qa, const Int_t eslength, Bool_t * es) :
TNamed("QA", "Quality Assurance status"),
fNdet(kNDET), 
fNEventSpecies(eslength), 
fLengthQA(qalength),
fQA(new ULong_t[fLengthQA]), 
fDet(kNULLDET),
fTask(kNULLTASK), 
fEventSpecie(AliRecoParam::kDefault), 
fEventSpecies(new Bool_t[fNEventSpecies])
{
  // constructor to be used
  memcpy(fQA, qa, fLengthQA*sizeof(ULong_t));
  memcpy(fEventSpecies, es, fNEventSpecies*sizeof(Bool_t));
}

//_______________________________________________________________
AliQAv1::AliQAv1(const DETECTORINDEX_t det) :
  TNamed("QA", "Quality Assurance status"),
  fNdet(kNDET), 
  fNEventSpecies(AliRecoParam::kNSpecies), 
  fLengthQA(fNdet*fNEventSpecies),
  fQA(new ULong_t[fLengthQA]), 
  fDet(det),
  fTask(kNULLTASK), 
  fEventSpecie(AliRecoParam::kDefault), 
  fEventSpecies(new Bool_t[fNEventSpecies])
{
  // constructor to be used
  if (! CheckRange(det) ) fDet = kNULLDET ; 
  memset(fQA,0,fLengthQA*sizeof(ULong_t));
  memset(fEventSpecies,kFALSE,fNEventSpecies*sizeof(Bool_t));
}
  
//_______________________________________________________________
AliQAv1::AliQAv1(const ALITASK_t tsk) :
  TNamed("QA", "Quality Assurance status"),
  fNdet(kNDET), 
  fNEventSpecies(AliRecoParam::kNSpecies), 
  fLengthQA(fNdet*fNEventSpecies),
  fQA(new ULong_t[fLengthQA]), 
  fDet(kNULLDET),
  fTask(tsk), 
  fEventSpecie(AliRecoParam::kDefault), 
  fEventSpecies(new Bool_t[fNEventSpecies])
{
  // constructor to be used in the AliRoot module (SIM, REC, ESD or ANA)
  if (! CheckRange(tsk) ) fTask = kNULLTASK ; 
  memset(fQA,0,fLengthQA*sizeof(ULong_t));
  memset(fEventSpecies,kFALSE,fNEventSpecies*sizeof(Bool_t));
}

//____________________________________________________________________________
AliQAv1::~AliQAv1() 
{
  // dtor  
  delete [] fQA;
  delete [] fEventSpecies;
}

//_______________________________________________________________
void AliQAv1::Close() 
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
Bool_t AliQAv1::CheckFatal() const
{
  // check if any FATAL status is set
  Bool_t rv = kFALSE ;
  Int_t index ;
  for (index = 0; index < kNDET ; index++)
    rv = rv || IsSet(DETECTORINDEX_t(index), fTask, fEventSpecie, kFATAL) ;
  return rv ;
}

//_______________________________________________________________
Bool_t AliQAv1::CheckRange(DETECTORINDEX_t det) const
{ 
  // check if detector is in given detector range: 0-kNDET

  Bool_t rv = ( det < 0 || det > kNDET )  ? kFALSE : kTRUE ;
  if (!rv)
    AliFatal(Form("Detector index %d is out of range: 0 <= index <= %d", det, kNDET)) ;
  return rv ;
}

//_______________________________________________________________
Bool_t AliQAv1::CheckRange(ALITASK_t task) const
{ 
  // check if task is given taskk range: 0:kNTASK
  Bool_t rv = ( task < kRAW || task > kNTASK )  ? kFALSE : kTRUE ;
  if (!rv)
    AliFatal(Form("Module index %d is out of range: 0 <= index <= %d", task, kNTASK)) ;
  return rv ;
}

//_______________________________________________________________
Bool_t AliQAv1::CheckRange(QABIT_t bit) const
{ 
  // check if bit is in given bit range: 0-kNBit

  Bool_t rv = ( bit < 0 || bit > kNBIT )  ? kFALSE : kTRUE ;
  if (!rv)
    AliFatal(Form("Status bit %d is out of range: 0 <= bit <= %d", bit, kNBIT)) ;
  return rv ;
}

//_______________________________________________________________
Bool_t AliQAv1::CheckRange(AliRecoParam::EventSpecie_t es) const
{ 
  // check if bit is in given bit range: 0-kNBit
  Bool_t rv = kFALSE ; 
  switch (es) {
    case AliRecoParam::kDefault: 
      rv = kTRUE ; 
      break ; 
    case AliRecoParam::kLowMult: 
      rv = kTRUE ; 
      break ; 
    case AliRecoParam::kHighMult: 
      rv = kTRUE ; 
      break ; 
    case AliRecoParam::kCosmic: 
      rv = kTRUE ; 
      break ; 
    case AliRecoParam::kCalib: 
      rv = kTRUE ; 
      break ; 
  }
  if (!rv)
    AliFatal(Form("Event Specie %d is not valid", es)) ;
  return rv ;
}

//_______________________________________________________________
const char * AliQAv1::GetAliTaskName(ALITASK_t tsk)
{
	// returns the char name corresponding to module index
	switch (tsk) {
		case kNULLTASK:
			break ; 
		case kRAW:
			return "RAW" ;
			break ;  
		case kSIM:
			return "SIM" ;
			break ;
		case kREC:
			return "REC" ;
			break ;
		case kESD:
			return "ESD" ;
			break ;
		case kANA:
			return "ANA" ;
			break ;
		default:
      return "" ; 
			break ;
	}
  return "" ;
}

//_______________________________________________________________
const char * AliQAv1::GetBitName(QABIT_t bit) const
{
	// returns the char name corresponding to bit 
	TString bitName ;
	switch (bit) {
		case kNULLBit:
			bitName = "NONE" ;
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
	return bitName.Data() ;
}

//_______________________________________________________________
TH1 * AliQAv1::GetData(TObjArray** list, Int_t index, AliRecoParam::EventSpecie_t eventSpecie)
{
    // retrieve QA data from the list at a given index and for a given event specie 
  TH1 * rv = NULL ; 
  Int_t esindex = AliRecoParam::AConvert(eventSpecie) ; 
  TObjArray * arr = list[esindex] ;
  if (arr) {
    if ( index > AliQAv1::GetMaxQAObj() ) {
			AliErrorClass(Form("Max number of authorized QA objects is %d", AliQAv1::GetMaxQAObj())) ; 
		} else {
      if ( arr->At(index) )  {
        rv = static_cast<TH1*>(arr->At(index)) ; 
      } 	
    }
  }  
  return rv ; 
}

//_______________________________________________________________
AliQAv1::DETECTORINDEX_t AliQAv1::GetDetIndex(const char * name) 
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
const char * AliQAv1::GetDetName(Int_t det) 
{
	// returns the detector name corresponding to a given index (needed in a loop)
	
	if ( det >= 0 &&  det < kNDET) 
		return (fgDetNames[det]).Data() ; 
	else 
		return NULL ; 
}

//_______________________________________________________________
TFile * AliQAv1::GetQADataFile(const char * name, Int_t run) 
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
			fgQADataFile = static_cast<TFile *>(gROOT->FindObject(temp)) ;
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
AliQAv1::MODE_t AliQAv1::Mode(TASKINDEX_t task) {
  // return "rec" or "sim" depending on the task
  
  switch (task) {
    case AliQAv1::kRAWS:
      return kRECMODE ; 
      break;
    case AliQAv1::kHITS:
      return kSIMMODE ; 
      break;
    case AliQAv1::kSDIGITS:
      return kSIMMODE ; 
      break;
    case AliQAv1::kDIGITS:
      return kSIMMODE ; 
      break;
    case AliQAv1::kDIGITSR:
      return kRECMODE ; 
      break;
    case AliQAv1::kRECPOINTS:
      return kRECMODE ; 
      break ; 
    case AliQAv1::kTRACKSEGMENTS:
      return kRECMODE ; 
      break;
    case AliQAv1::kRECPARTICLES:
      return kRECMODE ; 
      break;
    case AliQAv1::kESDS:
      return kRECMODE ; 
      break;
    default:
      break;
  }
  return AliQAv1::kNULLMODE;
}

//_____________________________________________________________________________
TFile * AliQAv1::GetQADataFile(const char * fileName)
{
  // Open if necessary the Data file and return its pointer

  if (!fgQADataFile) {
    if (!fileName) 
      fileName = AliQAv1::GetQADataFileName() ; 
    if  (!gSystem->AccessPathName(fileName)) {
      fgQADataFile =  TFile::Open(fileName) ;
    } else {
      AliFatalClass(Form("File %s not found", fileName)) ;
    }
  }
  return fgQADataFile ; 
}

//_______________________________________________________________
TFile * AliQAv1::GetQAResultFile() 
{
  // opens the file to store the  Quality Assurance Data Checker results	
  if (fgQAResultFile) {
   if (fgQAResultFile->IsOpen()) 
    fgQAResultFile->Close();
   delete fgQAResultFile;
  }
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
  return fgQAResultFile ; 
}

//_______________________________________________________________
AliQAv1::QABIT_t AliQAv1::GetQAStatusBit(AliRecoParam::EventSpecie_t es, DETECTORINDEX_t det, ALITASK_t tsk) const
{
    // returns the QA bit set
  QABIT_t rv = kNULLBit ; 
  if ( es == AliRecoParam::kDefault) 
    es = fEventSpecie ; 
  if ( det == kNULLDET ) 
    det = fDet ; 
  if ( tsk == kNULLTASK ) 
    tsk = fTask ; 
  for (Int_t bit = kINFO ; bit < kNBIT ; bit++) {
		if (IsSet(det, tsk, es, QABIT_t(bit))) 
      rv = QABIT_t(bit) ;
	}
  return rv ; 
}

//_______________________________________________________________
AliQAv1::TASKINDEX_t AliQAv1::GetTaskIndex(const char * name) 
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
Bool_t AliQAv1::IsSet(DETECTORINDEX_t det, ALITASK_t tsk, Int_t ies, QABIT_t bit) const 
{
  // Checks is the requested bit is set
   
  const AliRecoParam::EventSpecie_t es = AliRecoParam::Convert(ies) ; 
  return IsSet(det, tsk, es, bit) ; 
  
}  

//_______________________________________________________________
Bool_t AliQAv1::IsSet(DETECTORINDEX_t det, ALITASK_t tsk, AliRecoParam::EventSpecie_t es, QABIT_t bit) const
{
  // Checks is the requested bit is set
	
  CheckRange(det) ; 
  CheckRange(tsk) ;
  CheckRange(bit) ;
  CheckRange(es) ;
	
  ULong_t offset = Offset(tsk) ;
  ULong_t status = GetStatus(det, es) ;
  offset+= bit ;
  status = (status & 1 << offset) != 0 ;
  return status ;
}

//_______________________________________________________________
Bool_t AliQAv1::IsSetAny(DETECTORINDEX_t det, ALITASK_t tsk, AliRecoParam::EventSpecie_t es) const
{
  // Checks is the requested bit is set
	
  CheckRange(det) ; 
  CheckRange(tsk) ;
  CheckRange(es) ;
	
  ULong_t offset = Offset(tsk) ;
  ULong_t status = GetStatus(det, es) ;
	ULong_t st = 0 ; 
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
Bool_t AliQAv1::IsSetAny(DETECTORINDEX_t det, AliRecoParam::EventSpecie_t es) const
{
  // Checks is the requested bit is set
	
  CheckRange(det) ; 
  CheckRange(es) ; 
	
	ULong_t status = GetStatus(det, es) ;
	ULong_t st = 0 ; 
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
AliQAv1 * AliQAv1::Instance()
{
  // Get an instance of the singleton. The only authorized way to call the ctor

  if ( ! fgQA) {
    GetQAResultFile() ; 
    fgQA = static_cast<AliQAv1 *>(fgQAResultFile->Get("QA")) ; 
    if ( ! fgQA ) 
      fgQA = new AliQAv1() ;
  }
  return fgQA ;
}

//_______________________________________________________________
AliQAv1 * AliQAv1::Instance(const Int_t qalength, ULong_t * qa, const Int_t eslength, Bool_t * es)
{
  // Get an instance of the singleton. The only authorized way to call the ctor
  
  if ( ! fgQA) 
    fgQA = new AliQAv1(qalength, qa, eslength, es) ;
  return fgQA ;
}

//_______________________________________________________________
AliQAv1 * AliQAv1::Instance(const DETECTORINDEX_t det)
{
  // Get an instance of the singleton. The only authorized way to call the ctor
  
  if ( ! fgQA) {
    GetQAResultFile() ; 
    fgQA = static_cast<AliQAv1 *>(fgQAResultFile->Get(GetQAName())) ; 
    if ( ! fgQA ) 
      fgQA = new AliQAv1(det) ;
  }		
  fgQA->Set(det) ;
  return fgQA ;
}

//_______________________________________________________________
AliQAv1 * AliQAv1::Instance(const ALITASK_t tsk)
{
  // Get an instance of the singleton. The only authorized way to call the ctor

  if ( ! fgQA)
    switch (tsk) {
    case kNULLTASK:
      break ;
	case kRAW:
      fgQA = new AliQAv1(tsk) ;
      break ;
	case kSIM:
      fgQA = new AliQAv1(tsk) ;
      break ;
    case kREC:
      AliInfoClass("fgQA = gAlice->GetQA()") ;
      break ;
    case kESD:
      AliInfoClass("fgQA = static_cast<AliQAv1 *> (esdFile->Get(\"QA\")") ;
      break ;
    case kANA:
      AliInfoClass("fgQA = static_cast<AliQAv1 *> (esdFile->Get(\"QA\")") ;
      break ;
    case kNTASK:
      break ;
    }
  if (fgQA) 
    fgQA->Set(tsk) ;
  return fgQA ;
}

//_______________________________________________________________
AliQAv1 *  AliQAv1::Instance(const TASKINDEX_t tsk) 
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
Long64_t AliQAv1::Merge(const TCollection * list) {
	// Merge the QA resuls in the list into this single AliQAv1 object
	
  Long64_t nmerge = 0;
  for (Int_t det = 0 ; det < kNDET ; det++) {
    Set(DETECTORINDEX_t(det)) ; 
    for (Int_t task = 0 ; task < kNTASK ; task++) {
      Set(ALITASK_t(task)) ; 
      for (Int_t bit = 0 ; bit < kNBIT ; bit++) {
	TIter next(list) ;
	AliQAv1 * qa ; 
	while ( (qa = (AliQAv1*)next() ) ) {
          for (Int_t es = 0 ; es < fNEventSpecies ; es++) {
	    ++nmerge;
            if (qa->IsSet(DETECTORINDEX_t(det), ALITASK_t(task), es, QABIT_t(bit)))
              Set(QABIT_t(bit), es) ; 
          }
	} // qa list
      } // bit
    } // task
  } // detector
  return nmerge;
}

//_______________________________________________________________
ULong_t AliQAv1::Offset(ALITASK_t tsk) const
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
void AliQAv1::Reset(DETECTORINDEX_t det, ALITASK_t tsk, AliRecoParam::EventSpecie_t es)
{
    // resets all the bits for a given tsk and event specie
  CheckRange(det) ; 
  CheckRange(tsk) ;
  CheckRange(es) ;
	
  ULong_t offset = Offset(tsk) ;
  ULong_t status = GetStatus(det, es) ;
	for ( Int_t bit = 0 ; bit < kNBIT ; bit++) {
		offset+= bit ;
		status = status & ~1 << offset ;		
	}  
  SetStatus(det, es, status) ;
}

//_______________________________________________________________
void AliQAv1::ResetStatus(DETECTORINDEX_t det) 
{ 
  // reset the status of det for all event specie
  for (Int_t es = 0 ; es < fNEventSpecies ; es++)
    fQA[det*fNdet+es] = 0 ; 
}

//_______________________________________________________________
void AliQAv1::Set(QABIT_t bit, Int_t ies)
{
  // Set the status bit of the current detector in the current module and for the current event specie 
   Set(bit, AliRecoParam::Convert(ies)) ;
}

//_______________________________________________________________
void AliQAv1::Set(QABIT_t bit, AliRecoParam::EventSpecie_t es)
{
  // Set the status bit of the current detector in the current module and for the current event specie 
  
  SetStatusBit(fDet, fTask, es, bit) ;
}

//_____________________________________________________________________________
void AliQAv1::SetQARefStorage(const char * name)
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
	  AliErrorClass(Form("ERROR: %s is an invalid storage definition\n", name)) ; 
	  fgQARefDirName  = "" ; 
	  fgQARefFileName = "" ; 
  }	
	TString tmp(fgQARefDirName) ; // + fgQARefFileName) ;
	AliInfoClass(Form("AliQAv1::SetQARefDir: QA references are in  %s\n", tmp.Data() )) ;
}

//_____________________________________________________________________________
void AliQAv1::SetQAResultDirName(const char * name)
{
  // Set the root directory where to store the QA status object

  fgQAResultDirName.Prepend(name) ; 
  AliInfoClass(Form("AliQAv1::SetQAResultDirName: QA results are in  %s\n", fgQAResultDirName.Data())) ;
  if ( fgQAResultDirName.Contains(fgkLabLocalFile)) 
    fgQAResultDirName.ReplaceAll(fgkLabLocalFile, "") ;
  fgQAResultFileName.Prepend(fgQAResultDirName) ;
}

//_______________________________________________________________
void AliQAv1::SetStatusBit(DETECTORINDEX_t det, ALITASK_t tsk, AliRecoParam::EventSpecie_t es, QABIT_t bit)
{
 // Set the status bit for a given detector and a given task

  CheckRange(det) ;
  CheckRange(tsk) ;
  CheckRange(bit) ;
  CheckRange(es) ;
  
  Reset(det, tsk, es) ; 
  
  ULong_t offset = Offset(tsk) ;
  ULong_t status = GetStatus(det, es) ;
  offset+= bit ;
  status = status | 1 << offset ;
  SetStatus(det, es, status) ;
}

//_______________________________________________________________
void AliQAv1::Show(DETECTORINDEX_t det) const 
{ 
  // dispplay the QA status word
  if ( det == kNULLDET) 
    det = fDet ;  
  for (Int_t ies = 0 ; ies < fNEventSpecies ; ies++) {
    if ( IsEventSpecieSet(ies) )
      ShowStatus(det, kNULLTASK, AliRecoParam::ConvertIndex(ies)) ; 
  }
}

//_______________________________________________________________
void AliQAv1::ShowAll() const 
{
  // dispplay the QA status word
  Int_t index ;
  for (index = 0 ; index < kNDET ; index++) {
		for (Int_t tsk = kRAW ; tsk < kNTASK ; tsk++) {
      for (Int_t ies = 0 ; ies < fNEventSpecies ; ies++) {
        if ( IsEventSpecieSet(ies) )
          ShowStatus(DETECTORINDEX_t(index), ALITASK_t(tsk), AliRecoParam::ConvertIndex(ies)) ;
      }
    }
	}
}

//_______________________________________________________________
void AliQAv1::ShowStatus(DETECTORINDEX_t det, ALITASK_t tsk, AliRecoParam::EventSpecie_t es) const
{
	// Prints the full QA status of a given detector
	CheckRange(det) ;
	CheckRange(es) ;
	ULong_t status = GetStatus(det, es) ;
	ULong_t tskStatus[kNTASK] ; 
	tskStatus[kRAW] = status & 0x0000f ;
	tskStatus[kSIM] = status & 0x000f0 ;
	tskStatus[kREC] = status & 0x00f00 ;
	tskStatus[kESD] = status & 0x0f000 ;
	tskStatus[kANA] = status & 0xf0000 ;

	AliInfo(Form("====> QA Status for %8s %8s raw =0x%lx, sim=0x%lx, rec=0x%lx, esd=0x%lx, ana=0x%lx", GetDetName(det).Data(), AliRecoParam::GetEventSpecieName(es), 
				 tskStatus[kRAW], tskStatus[kSIM], tskStatus[kREC], tskStatus[kESD], tskStatus[kANA] )) ;
	if (tsk == kNULLTASK) {
		for (Int_t itsk = kRAW ; itsk < kNTASK ; itsk++) {
			ShowASCIIStatus(es, det, ALITASK_t(itsk), tskStatus[itsk]) ; 
		} 
	} else {
			ShowASCIIStatus(es, det, tsk, tskStatus[tsk]) ; 
	}
}

//_______________________________________________________________
void AliQAv1::ShowASCIIStatus(AliRecoParam::EventSpecie_t es, DETECTORINDEX_t det, ALITASK_t tsk, const ULong_t status) const 
{
	// print the QA status in human readable format
  const QABIT_t bit = GetQAStatusBit(es, det, tsk) ; 
  if ( bit != kNULLBit ) {
          AliInfoClass(Form("           %8s %8s %4s  \n", AliRecoParam::GetEventSpecieName(es), GetDetName(det).Data(), GetAliTaskName(tsk))) ; 
          AliInfoClass(Form("           %8s %8s %4s 0x%4lx, Problem signalled: %8s \n", AliRecoParam::GetEventSpecieName(es), GetDetName(det).Data(), GetAliTaskName(tsk), status, GetBitName(bit))) ; 
  }
}

//_______________________________________________________________
void AliQAv1::UnSet(QABIT_t bit, Int_t ies)
{
	// UnSet the status bit of the current detector in the current module
		UnSet(bit, AliRecoParam::Convert(ies)) ;
}

//_______________________________________________________________
void AliQAv1::UnSet(QABIT_t bit, AliRecoParam::EventSpecie_t es)
{
	// UnSet the status bit of the current detector in the current module
	
	UnSetStatusBit(fDet, fTask, es, bit) ;
}

//_______________________________________________________________
void AliQAv1::UnSetStatusBit(DETECTORINDEX_t det, ALITASK_t tsk, AliRecoParam::EventSpecie_t es, QABIT_t bit)
{
	// UnSet the status bit for a given detector and a given task
	
	CheckRange(det) ;
	CheckRange(tsk) ;
	CheckRange(bit) ;
	CheckRange(es) ;
	
	ULong_t offset = Offset(tsk) ;
	ULong_t status = GetStatus(det, es) ;
	offset+= bit ;
	status = status & 0 << offset ;
	SetStatus(det, es, status) ;
}
