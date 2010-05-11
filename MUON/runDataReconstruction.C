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

void runDataReconstruction(const char* input = "raw://run117099",
                           const char* ocdbPath = "raw://",
                           const char* recoptions="SAVEDIGITS",
                           Int_t numberOfEvents=10000)
{ 
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage(ocdbPath);

  AliReconstruction MuonRec;
  
  MuonRec.SetInput(gSystem->ExpandPathName(input));
  MuonRec.SetRunVertexFinder(kFALSE);
  MuonRec.SetRunLocalReconstruction("MUON");
  MuonRec.SetRunTracking("MUON");
  MuonRec.SetFillESD(" ");
  MuonRec.SetLoadAlignData("MUON");
  MuonRec.SetNumberOfEventsPerFile(500); // must set a limit otherwise time per event increases with event number (this is a "bug" of the loaders)
  MuonRec.SetOption("MUON",recoptions);  
  MuonRec.SetRunQA("MUON:ALL");
  MuonRec.SetQAWriteExpert(AliQAv1::kMUON);

  if ( numberOfEvents > 0 )
  {
    MuonRec.SetEventRange(0,numberOfEvents);
  }
  MuonRec.Run();
  
}

