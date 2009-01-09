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
#include "TFile.h"
#include "TTree.h"

// MUON includes
#include "AliMUONCDB.h"
#include "AliMUONTriggerEfficiencyCells.h"
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
/// \author Diego Stocco, INFN Torino

void MUONTriggerChamberEfficiency(Bool_t addMapInSimulation=kFALSE,
				  const char *inputDir=".")
{
/// \param addMapInSimulation (default kFALSE);
///     kTRUE: creates file MUON/Calib/TriggerEfficiency/Run0_99999999_v0_s?.root
///            with calculated chamber efficiency which can be used in the next simulation
/// \param inputDir 
///     path to AliESDs.root (default ".")

    Char_t filename[150], *className = "AliMUONTriggerEfficiencyCells";
    sprintf(filename,"%s/AliESDs.root",inputDir);

    TFile *file = new TFile(filename,"read");
    if(!file){
	cerr << "Cannot find " << filename << "\nExit!" << endl;
	return;
    }
    
    TTree *esdTree = (TTree*)file->Get("esdTree");
    if(!esdTree){
	cerr << "Cannot find esdTree in " << filename << "\nExit!" << endl;
	return;
    }
    
    AliMUONTriggerEfficiencyCells *effMap = 
      (AliMUONTriggerEfficiencyCells*)esdTree->GetUserInfo()->FindObject(className);
    if(!effMap){
	cerr << "Cannot find " << className << " in esdTree.\nExit!" << endl;
	return;
    }

    effMap->DisplayEfficiency();

    if(!addMapInSimulation) return;

    AliMUONCDB muonCDB;
  muonCDB.WriteToCDB("MUON/Calib/TriggerEfficiency",effMap,0,AliCDBRunRange::Infinity(),true);
}
