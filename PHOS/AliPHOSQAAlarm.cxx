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
// An alarm object that is instanciated by a AliPHOSQACheckable in response to
// a AliPHOSQAChecker
//*-- Author :  Yves Schutz (SUBATECH) 
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---

// --- Standard library ---

// --- AliRoot header files ---

#include "AliPHOSQAAlarm.h"
#include "AliRun.h"

ClassImp(AliPHOSQAAlarm)


//____________________________________________________________________________ 
  AliPHOSQAAlarm::AliPHOSQAAlarm(TString time, TString checked, TString checker, TString  message)  
{
  // ctor
 fTime    = time ; 
 fCable   = checked ; 
 fCer     = checker ; 
 fMessage = message ; 
 fEvent   = gAlice->GetEvNumber() ;
}

//____________________________________________________________________________ 
  AliPHOSQAAlarm::~AliPHOSQAAlarm()
{
  // dtor
}

//____________________________________________________________________________ 
  void AliPHOSQAAlarm::Print()
{
  // print the message 

  Info("Print", "Event# %d %s", fEvent, fMessage.Data()) ;  
}
