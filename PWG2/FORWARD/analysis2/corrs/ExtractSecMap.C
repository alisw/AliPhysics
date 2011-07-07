/**
 * Script to draw the energy loss fits 
 * 
 * @ingroup pwg2_forward_analysis_scripts
 */
#ifndef __CINT__
#include <TFile.h>
#include <TList.h>
#include <TError.h>
#include "AliFMDCorrSecondaryMap.h"
#include "AliCentralCorrSecondaryMap.h"
#include "AliForwardCorrectionManager.h"
#endif 

//____________________________________________________________________
/** 
 * Extract the energy loss correction object from file and rename it 
 * according to the settings 
 * 
 * @param fname  File to extract from 
 * @param sys    Collision system (pp, PbPb)
 * @param sNN    Center of mass energy (in GeV) per nucleon
 * @param field  L3 magnetic field (-5,0,5) in kGaus
 * 
 * @ingroup pwg2_forward_analysis_scripts
 */
void
ExtractSecMap(const char* fname,
	     UShort_t sys=1, UShort_t sNN=900, Short_t field=5)
{
#ifdef __CINT__
  gROOT->Macro("$ALICE_ROOT/PWG2/FORWARD/analysis2/scripts/LoadLibs.C");
#endif

  TFile* file = TFile::Open(fname, "READ");
  if (!file) {
    Error("ExtractSecMap", "Couldn't open %s", fname);
    return;
  }

  ExtractFMDSecMap(file, sys, sNN, field);
  ExtractSPDSecMap(file, sys, sNN, field);
}

//____________________________________________________________________
/** 
 * Extract and copy FMD secondary map to file 
 * 
 * @param file  Input file
 * @param sys   Collision system (1:pp, 2:PbPb)
 * @param sNN   Center of mass energy (GeV) per nucleon
 * @param field L3 magnetic field
 *
 * @ingroup pwg2_forward_analysis_scripts
 */
void     
ExtractFMDSecMap(TFile* file, UShort_t sys, UShort_t sNN, Short_t field)
{
  TList* forward = static_cast<TList*>(file->Get("ForwardResults"));
  // static_cast<TList*>(file->Get("PWG2forwardDnDeta/Forward"));
  if (!forward) { 
    Error("ExtractSecMap", "Couldn't get forward list from %s", fname);
    return;
  }
  
  TString n(AliFMDCorrSecondaryMap::Class()->GetName());
  TObject* fmdCorr = forward->FindObject(n);
  if (!fmdCorr) { 
    Error("ExtractSecMap", "Couldn't get forward correction object %s", 
	  n.Data());
    return;
  }

  AliForwardCorrectionManager& mgr = AliForwardCorrectionManager::Instance();
  mgr.WriteFile(AliForwardCorrectionManager::kSecondaryMap, 
		sys, sNN, field, false, fmdCorr, false);
}

//____________________________________________________________________
/** 
 * Extract and copy SPD secondary map to file 
 * 
 * @param file  Input file
 * @param sys   Collision system (1:pp, 2:PbPb)
 * @param sNN   Center of mass energy (GeV) per nucleon
 * @param field L3 magnetic field
 *
 * @ingroup pwg2_forward_analysis_scripts
 */
void     
ExtractSPDSecMap(TFile* file, UShort_t sys, UShort_t sNN, Short_t field)
{
  TList* central = static_cast<TList*>(file->Get("CentralResults"));
  // static_cast<TList*>(file->Get("PWG2centralDnDeta/Central"));
  if (!central) { 
    Error("ExtractSecMap", "Couldn't get central list from %s", fname);
    return;
  }
  
  TString n(AliCentralCorrSecondaryMap::Class()->GetName());
  TObject* spdCorr  = central->FindObject(n);
  if (!spdCorr) { 
    Error("ExtractSecMap", "Couldn't get central correction object %s", 
	  n.Data());
    return;
  }


  AliCentralMultiplicityTask::Manager* mgr = new 
    AliCentralMultiplicityTask::Manager;
  // mgr->Dump();
  mgr->Print();
  mgr->WriteFile(0, sys, sNN, field, spdCorr, false);
}
    
  
//____________________________________________________________________
/** 
 * Extract the energy loss correction object from file and rename it 
 * according to the settings 
 * 
 * @param fname  File to extract from 
 * @param sys    Collision system (pp, PbPb)
 * @param sNN    Center of mass energy (in GeV) per nucleon
 * @param field  L3 magnetic field (-5,0,5) in kGaus
 * @param mc     Whether this is from MC data or not 
 * 
 * @ingroup pwg2_forward_analysis_scripts
 */
void
ExtractSecMap(const char* fname="forward_mccorr.root", 
	      const char* sys="p-p", 
	      Float_t     sNN=900, 
	      Float_t     field=5)
{
  gROOT->Macro("$ALICE_ROOT/PWG2/FORWARD/analysis2/scripts/LoadLibs.C");

  UShort_t uSys   = AliForwardUtil::ParseCollisionSystem(sys);
  UShort_t usNN   = AliForwardUtil::ParseCenterOfMassEnergy(uSys,sNN);
  Short_t  sField = AliForwardUtil::ParseMagneticField(field);

  ExtractSecMap(fname, uSys, usNN, sField);
}

//____________________________________________________________________
//
// EOF
//
