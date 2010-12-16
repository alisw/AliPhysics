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
void
ExtractELoss(const char* fname="energyFits.root", 
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

  TString ofName(mgr.GetFileName(AliForwardCorrectionManager::kELossFits,
				 sys, sNN, field, mc));
  TFile* output = TFile::Open(ofName.Data(), "RECREATE");
  if (!output) { 
    Error("ExtractELoss", "Failed to open file %s", ofName.Data());
    return;
  }

  TString oName(mgr.GetObjectName(AliForwardCorrectionManager::kELossFits));
  obj->Write(oName);

  output->Write();
  output->ls();
  output->Close();
  
  TString dName(mgr.GetFileDir(AliForwardCorrectionManager::kELossFits));
  Info("ExtractELoss", "Wrote %s object %s to %s",cName.Data(),oName.Data(), 
       ofName.Data());
  Info("ExtractELoss", "%s should be copied to %s",ofName.Data(),dName.Data());
  Info("ExtractELoss", "Do for example\n\t"
       "aliroot $ALICE_ROOT/PWG2/FORWARD/analysis2/scripts/MoveCorrections.C\(0,0,0,0,1\)\n\t"
       "cp %s %s/", ofName.Data(), gSystem->ExpandPathName(dName.Data()));
}

    
  
//____________________________________________________________________
void
ExtractELoss(const char* fname="energyFits.root", 
	     const char* sys="p-p", 
	     Float_t     sNN=900, 
	     Float_t     field=5)
{
  UShort_t uSys   = AliForwardUtil::ParseCollisionSystem(sys);
  UShort_t usNN   = AliForwardUtil::ParseCenterOfMassEnergy(uSys,sNN);
  Short_t  sField = AliForwardUtil::ParseMagneticField(field);

  ExtractELoss(fname, uSys, usNN, sField);
}

//____________________________________________________________________
//
// EOF
//
