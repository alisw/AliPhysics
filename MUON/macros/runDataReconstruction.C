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

/* $Id: runReconstruction.C 23207 2007-12-20 09:59:20Z ivana $ */

/// \ingroup macros
/// \file runDataReconstruction.C
/// \brief Macro for running reconstruction
///
/// Macro for running reconstruction on the cosmics run data.
///
/// \author Laurent Aphecetche, Nicole Bastid, Bogdan Vulpescu, ...

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliCDBManager.h"
#include "AliReconstruction.h"
#include <TGrid.h>
#include <TSystem.h>
#endif

void runDataReconstruction(const char* input = "raw://run124360",
                           const char* ocdbPath = "raw://",
                           const char* recoptions="SAVEDIGITS",
                           Int_t numberOfEvents=1000)
{ 
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage(ocdbPath);

  AliReconstruction rec;
  
  rec.SetRunReconstruction("MUON");

  rec.SetRunQA("MUON:ALL");

  rec.SetQARefDefaultStorage("local://$ALICE_ROOT/QAref") ;

  rec.SetWriteESDfriend(kTRUE);
  rec.SetWriteAlignmentData();

  rec.SetInput(gSystem->ExpandPathName(input));

  rec.SetUseTrackingErrorsForAlignment("ITS");

  rec.SetCleanESD(kFALSE);

  rec.SetStopOnError(kFALSE);

  rec.SetOption("MUON",recoptions);  

  rec.SetQAWriteExpert(AliQAv1::kMUON);

  if ( numberOfEvents > 0 )
  {
    rec.SetEventRange(0,numberOfEvents);
  }

  AliLog::Flush();
  rec.Run();
  
}

