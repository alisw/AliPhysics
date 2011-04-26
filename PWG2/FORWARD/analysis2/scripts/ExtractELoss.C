/**
 * Script to draw the energy loss fits 
 * 
 * @ingroup pwg2_forward_analysis_scripts
 */
#ifndef __CINT__
#include <TFile.h>
#include <TList.h>
#include <TError.h>
#include "AliFMDCorrELossFit.h"
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
 * @param mc     Whether this is from MC data or not 
 * 
 * @ingroup pwg2_forward_analysis_scripts
 */
void
ExtractELoss(const char* fname="forward_eloss.root", 
	     UShort_t sys=1, UShort_t sNN=900, Short_t field=5, Bool_t mc=false)
{
#ifdef __CINT__
  gROOT->Macro("$ALICE_ROOT/PWG2/FORWARD/analysis2/scripts/LoadLibs.C");
#endif

  TFile* file = TFile::Open(fname, "READ");
  if (!file) {
    Error("ExtractELoss", "Couldn't open %s", fname);
    return;
  }
    
  TList* forward = static_cast<TList*>(file->Get("Forward"));
  // static_cast<TList*>(file->Get("PWG2forwardDnDeta/Forward"));
  if (!forward) { 
    Error("ExtractELoss", "Couldn't get forward list from %s", fname);
    return;
  }
  
  TList* fitter = static_cast<TList*>(forward->FindObject("fmdEnergyFitter"));
  if (!fitter) { 
    Error("ExtractELoss", "Couldn't get fitter folder");
    return;
  }

  TString cName(AliFMDCorrELossFit::Class()->GetName());

  AliFMDCorrELossFit* obj = 
    static_cast<AliFMDCorrELossFit*>(fitter->FindObject(cName));
  if (!obj) {
    Error("ExtractELoss", "Couldn't get %s correction object", cName.Data());
    return;
  }

  AliForwardCorrectionManager& mgr = AliForwardCorrectionManager::Instance();
  mgr.WriteFile(AliForwardCorrectionManager::kELossFits, 
		sys, sNN, field, mc, obj, false);
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
ExtractELoss(const char* fname="energyFits.root", 
	     const char* sys="p-p", 
	     Float_t     sNN=900, 
	     Float_t     field=5,
	     Bool_t      mc=false)
{
  UShort_t uSys   = AliForwardUtil::ParseCollisionSystem(sys);
  UShort_t usNN   = AliForwardUtil::ParseCenterOfMassEnergy(uSys,sNN);
  Short_t  sField = AliForwardUtil::ParseMagneticField(field);

  ExtractELoss(fname, uSys, usNN, sField, mc);
}

//____________________________________________________________________
//
// EOF
//
