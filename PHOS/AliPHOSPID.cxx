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

//_________________________________________________________________________
//  Algorithm class for the identification of particles detected in PHOS        
//  base  class  of identified particle  
//  Why should I put meaningless comments
//  just to satisfy
//  the code checker                
                         
//                  
//*-- Author: Yves Schutz (SUBATECH) & Dmitri Peressounko


// --- ROOT system ---

// --- Standard library ---

// --- AliRoot header files ---
#include "AliConfig.h"
#include "AliPHOSPID.h"
#include "AliPHOSGetter.h"
#include "AliPHOSQualAssDataMaker.h" 

ClassImp(AliPHOSPID)

//____________________________________________________________________________
AliPHOSPID::AliPHOSPID():
  TTask("",""),
  fEventFolderName(""),
  fFirstEvent(0),
  fLastEvent(-1),
  fESD(0x0), 
  fQADM(0x0)
{
  // ctor
}


//____________________________________________________________________________
AliPHOSPID::AliPHOSPID(const TString alirunFileName, const TString eventFolderName):
  TTask("PHOS"+AliConfig::Instance()->GetPIDTaskName(), alirunFileName), 
  fEventFolderName(eventFolderName),
  fFirstEvent(0),
  fLastEvent(-1), 
  fESD(0x0), 
  fQADM(0x0)
{
  // ctor
  fQADM = new  AliPHOSQualAssDataMaker() ; //!Quality Assurance Data Maker
  GetQualAssDataMaker()->Init(AliQualAss::kRECPARTICLES) ;    
}

//____________________________________________________________________________
AliPHOSPID::AliPHOSPID(const AliPHOSPID & pid) :
  TTask(pid),fEventFolderName(pid.GetEventFolderName()),
  fFirstEvent(pid.GetFirstEvent()),fLastEvent(pid.GetLastEvent()), 
  fESD(pid.fESD), 
  fQADM(pid.fQADM)
{
  // Copy constructor
}
//____________________________________________________________________________
AliPHOSPID::~AliPHOSPID()
{
  // dtor
 //Remove this from the parental task before destroying
  if(AliPHOSGetter::Instance()->PhosLoader())
    AliPHOSGetter::Instance()->PhosLoader()->CleanPIDTask();
  delete fQADM ; 
}

