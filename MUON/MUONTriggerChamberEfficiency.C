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
#include "AliMUONTriggerChamberEff.h"
#include "AliMUONTriggerEfficiencyCells.h"

#include "Riostream.h"

#endif

// Macro to view and save the trigger chamber efficiency map 
// calculated during reconstruction.
// Efficiency map can be made available for next simulation.

// Arguments:
//
// addMapInSimulation (default kFALSE):
//    kTRUE: creates file MUON/Calib/TriggerEfficiency/Run0_99999999_v0_s?.root
//           with calculated chamber efficiency which can be used in the next simulation
//
// inputDir (default "."):
//    path to AliESDs.root
// outDir (default 0x0)
//    directory where to store output file "MUON.TriggerEfficiencyMap.root"
//    (if not set outDir=inputDir is assumed)

void MUONTriggerChamberEfficiency(Bool_t addMapInSimulation=kFALSE,
				  const char *inputDir=".",
				  const char *outDir=0x0)
{
    char filename[150], outfileDir[150];
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
    
    AliMUONTriggerChamberEff *chEff = (AliMUONTriggerChamberEff*)esdTree->GetUserInfo()->FindObject("AliMUONTriggerChamberEff");
    if(!chEff){
	cerr << "Cannot find AliMUONTriggerChamberEff in esdTree.\nExit!" << endl;
	return;
    }

    if(!outDir) sprintf(outfileDir,"%s",inputDir);
    else sprintf(outfileDir,"%s",outDir);

    chEff->WriteEfficiencyMap(outfileDir);
    chEff->DisplayEfficiency();

    if(!addMapInSimulation) return;

    char outFileName[150];
    sprintf(outFileName,"%s/MUON.TriggerEfficiencyMap.root",outfileDir);
    
    AliMUONTriggerEfficiencyCells *effCells = new AliMUONTriggerEfficiencyCells(outFileName);
    AliMUONCDB muonCDB;
    muonCDB.WriteToCDB("MUON/Calib/TriggerEfficiency",effCells,0,99999999,true);
}
