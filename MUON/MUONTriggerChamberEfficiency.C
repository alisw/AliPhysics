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

#if !defined(__CINT__) || defined(__MAKECINT__)
// ROOT includes
#include "TGrid.h"
#include "TString.h"

// MUON includes
#include "AliMUONCDB.h"
#include "AliMUONCalibrationData.h"
#include "AliMUONTriggerEfficiencyCells.h"
#include "AliMUONTriggerChamberEfficiency.h"
#include "AliCDBManager.h"
#include "AliCDBRunRange.h"
#include "Riostream.h"

#endif

/// \ingroup macros
/// \file MUONTriggerChamberEfficiency.C
/// \brief Macro to view and save the trigger chamber efficiency map 
/// calculated during reconstruction.
///
/// Efficiency map can be made available for next simulation.
///
/// \author Diego Stocco, Subatech, Nantes

void MUONTriggerChamberEfficiency(const Char_t* inputFile="./MUON.TriggerEfficiencyMap.root",
				  Bool_t addMapInLocalOCDB = kFALSE,
				  Int_t firstRun=0, Int_t lastRun = AliCDBRunRange::Infinity())
{
/// \param inputFile (default "./MUON.TriggerEfficiencyMaps.root")
///     File with the numerator and denominator histos for efficiency calculation
///     (It is the output of the PWG3/muon/AliAnalysisTaskTrigChEff analysis
/// \param addMapInLocalOCDB (default kFALSE)
///     if the map is OK, add it in local OCDB
/// \param firstRun (default 0)
///     first run of validity for CDB object
/// \param lastRun (default AliCDBRunRange::Infinity())
///     last run of validity for CDB Object

  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");

  AliMUONTriggerEfficiencyCells* effMap = new AliMUONTriggerEfficiencyCells(inputFile);
  if ( ! addMapInLocalOCDB ){
    AliCDBManager::Instance()->SetRun(firstRun);
    AliMUONTriggerChamberEfficiency trigChEff(effMap);
 
    trigChEff.DisplayEfficiency();
  }
  else
    AliMUONCDB::WriteToCDB("MUON/Calib/TriggerEfficiency", effMap, firstRun, lastRun,true);
}

//____________________________________________________________
void ShowOCDBmap(Int_t runNumber = 0, TString ocdbPath="local://$ALICE_ROOT/OCDB")
{
/// \param runNumber (default 0)
///     run number
/// \param ocdbPath(default "local://$ALICE_ROOT/OCDB")
///     path to OCDB
  if ( ocdbPath.BeginsWith("alien://") || ocdbPath.BeginsWith("raw://"))
    TGrid::Connect("alien://");

  AliCDBManager::Instance()->SetDefaultStorage(ocdbPath.Data());
  AliCDBManager::Instance()->SetRun(runNumber);
  AliMUONCalibrationData calib(runNumber);

  AliMUONTriggerChamberEfficiency trigChEff(calib.TriggerEfficiency());
  trigChEff.DisplayEfficiency();
}

