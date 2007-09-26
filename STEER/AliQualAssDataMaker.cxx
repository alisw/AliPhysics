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

/*
  Base Class
  Produces the data needed to calculate the quality assurance. 
  All data must be mergeable objects.
  Y. Schutz CERN July 2007
*/

// --- ROOT system ---
#include <TSystem.h> 
#include <TFile.h>
#include <TClonesArray.h> 
#include <TTree.h>

// --- Standard library ---

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliQualAssDataMaker.h"
#include "AliESDEvent.h"

ClassImp(AliQualAssDataMaker)
  
TString AliQualAssDataMaker::fDetectorDirName("") ;

           
//____________________________________________________________________________ 
AliQualAssDataMaker::AliQualAssDataMaker(const char * name, const char * title) : 
  TNamed(name, title), 
  fOutput(0x0),
  fDetectorDir(0x0)
{
  // ctor
  TString tmp(GetName()) ; 
  tmp.Append("QA") ; 
  SetName(tmp.Data()) ; 
  fOutput = AliQualAss::GetQADMOutFile() ; 	
  fDetectorDirName = GetName() ; 
}

//____________________________________________________________________________ 
AliQualAssDataMaker::AliQualAssDataMaker(const AliQualAssDataMaker& qadm) :
  TNamed(qadm.GetName(), qadm.GetTitle()),
  fOutput(qadm.fOutput),
  fDetectorDir(qadm.fDetectorDir)
{
  //copy ctor
  fDetectorDirName = GetName() ; 
}

//__________________________________________________________________
AliQualAssDataMaker& AliQualAssDataMaker::operator = (const AliQualAssDataMaker& qadm )
{
  // Equal operator.
  this->~AliQualAssDataMaker();
  new(this) AliQualAssDataMaker(qadm);
  return *this;
}

//____________________________________________________________________________
void AliQualAssDataMaker::Exec(AliQualAss::TASKINDEX task, TObject * data) 
{ 
  // creates the quality assurance data for the various tasks (Hits, SDigits, Digits, ESDs)
 
  fDetectorDir = fOutput->GetDirectory(GetDetectorDirName()) ; 
  if (!fDetectorDir)
   fDetectorDir = fOutput->mkdir(GetDetectorDirName()) ; 
   
  switch (task) { 
  case AliQualAss::kHITS:
    AliInfo("Processing Hits QA") ; 
    MakeHits(data) ;
    break ; 

  case AliQualAss::kSDIGITS:
    AliInfo("Processing SDigits QA") ; 
    MakeSDigits(data) ;
    break ; 
    
  case AliQualAss::kDIGITS:
    MakeDigits(data) ;
    break ;  
 
   case AliQualAss::kRECPOINTS:
     AliInfo("Processing RecPoints QA") ; 
     {
       TTree * recpoints = dynamic_cast<TTree *>(data) ; 
       if (recpoints) 
	 MakeRecPoints(recpoints) ;
       else 
	 AliError("Wrong type of recpoints container") ; 
     }
    break ;  

   case AliQualAss::kTRACKSEGMENTS:
    AliInfo("Processing Track Segments QA: not existing anymore") ; 
//     TTree * ts = dynamic_cast<TTree *>(data) ; 
//     if (ts) 
//       MakeTrackSegments(ts) ;
//     else 
//       AliError("Wrong type of track segments container") ; 
    break ;  
  
  case AliQualAss::kRECPARTICLES:
    AliInfo("Processing RecParticles QA: not existing anymore") ; 
//     TTree * recpar = dynamic_cast<TTree *>(data) ; 
//     if (recpar) 
//       MakeRecParticles(recpar) ;
//     else 
//       AliError("Wrong type of recparticles container") ; 
    break ;  
    
  case AliQualAss::kESDS:
    AliInfo("Processing ESDs QA") ; 
    AliESDEvent * esd = dynamic_cast<AliESDEvent *>(data) ; 
    if (esd) 
      MakeESDs(esd) ;
    else 
      AliError("Wrong type of esd container") ; 
    break ;  
  }	
}

//____________________________________________________________________________ 
void AliQualAssDataMaker::Finish(AliQualAss::TASKINDEX task) const 
{ 
  // write to the output File 

  fDetectorDir->cd() ; 
  TDirectory * subDir = 0x0 ; 
  
//   switch (task) { 
//   case AliQualAss::kHITS:
  subDir = fDetectorDir->GetDirectory(AliQualAss::GetTaskName(task)) ; 
//     break ; 

//    case AliQualAss::kSDIGITS:
// 	subDir = fDetectorDir->GetDirectory("SDigits") ; 
//     break ; 

//    case AliQualAss::kDIGITS:
//     subDir = fDetectorDir->GetDirectory("Digits") ; 
//     break ;  
 
//    case AliQualAss::kRECPOINTS:
//     subDir = fDetectorDir->GetDirectory("RecPoints") ; 
//     break ;  

//    case AliQualAss::kTRACKSEGMENTS:
//     subDir = fDetectorDir->GetDirectory("TrackSegments") ; 
//     break ;  
  
//    case AliQualAss::kRECPARTICLES:
//     subDir = fDetectorDir->GetDirectory("RecParticles") ; 
// 	break ;  
     
//   case AliQualAss::kESDS:
//     subDir = fDetectorDir->GetDirectory("ESDs") ; 
//     break ;  
//   }	
  subDir->Write() ; 
} 

//____________________________________________________________________________ 
void AliQualAssDataMaker::Init(AliQualAss::TASKINDEX task)
{
  // general intialisation
  
  TDirectory * subDir = 0x0 ; 
  
  switch (task) { 
  case AliQualAss::kHITS: 
    subDir = fDetectorDir->GetDirectory("Hits") ; 
	if (!subDir)
      subDir = fDetectorDir->mkdir("Hits") ; 
	subDir->cd() ; 
    InitHits() ;
    break ; 
  
  case AliQualAss::kSDIGITS: 
	subDir = fDetectorDir->GetDirectory("SDigits") ; 
	if (!subDir)
		subDir = fDetectorDir->mkdir("SDigits") ; 
	subDir->cd() ; 
    InitSDigits() ;
    break ; 

  case AliQualAss::kDIGITS: 
	subDir = fDetectorDir->GetDirectory("Digits") ; 
	if (!subDir)
		subDir = fDetectorDir->mkdir("Digits") ; 
	subDir->cd() ; 
	InitDigits() ;
	break ; 
	  
  case AliQualAss::kRECPOINTS: 
	subDir = fDetectorDir->GetDirectory("RecPoints") ; 
	if(!subDir)
		subDir = fDetectorDir->mkdir("RecPoints") ; 
	subDir->cd() ; 
    InitRecPoints() ;
    break ; 

  case AliQualAss::kTRACKSEGMENTS: 
	subDir = fDetectorDir->GetDirectory("TrackSegments") ; 
	if (!subDir)
		subDir = fDetectorDir->mkdir("TrackSegments") ; 
	subDir->cd() ; 
    InitTrackSegments() ;
    break ; 
    
  case AliQualAss::kRECPARTICLES: 
	subDir = fDetectorDir->GetDirectory("RecParticles") ; 
	if (!subDir)
		subDir = fDetectorDir->mkdir("RecParticles") ; 
	subDir->cd() ; 
    InitRecParticles() ;
    break ; 
    
  case AliQualAss::kESDS: 
	subDir = fDetectorDir->GetDirectory("ESDs") ; 
	if (!subDir) 
		subDir = fDetectorDir->mkdir("ESDs") ;
	subDir->cd() ; 
    InitESDs() ;
    break ; 
  }  
}
