// $Id$
//
// Task to setup emcal related global objects.
//
// Author: C.Loizides

#include "AliEmcalSetupTask.h"
#include <TClonesArray.h>
#include <TGeoGlobalMagField.h>
#include <TGeoManager.h>
#include "AliAODEvent.h"
#include "AliAnalysisManager.h"
#include "AliCDBManager.h"
#include "AliEMCALGeometry.h"
#include "AliESDEvent.h"
#include "AliGRPManager.h"
#include "AliGeomManager.h"
#include "AliMagF.h"
#include "AliOADBContainer.h"

ClassImp(AliEmcalSetupTask)

//________________________________________________________________________
AliEmcalSetupTask::AliEmcalSetupTask() : 
  AliAnalysisTaskSE(),
  fOcdbPath(),
  fOadbPath("$ALICE_ROOT/OADB/EMCAL"),
  fGeoPath("."),
  fIsInit(kFALSE)
{
  // Constructor.
}

//________________________________________________________________________
AliEmcalSetupTask::AliEmcalSetupTask(const char *name) : 
  AliAnalysisTaskSE(name),
  fOcdbPath(),
  fOadbPath("$ALICE_ROOT/OADB/EMCAL"),
  fGeoPath("$ALICE_ROOT/OADB/EMCAL"),
  fIsInit(kFALSE)
{
  // Constructor.
  fBranchNames = "ESD:AliESDHeader.,AliESDRun.";
}

//________________________________________________________________________
AliEmcalSetupTask::~AliEmcalSetupTask()
{
  // Destructor.
}

//________________________________________________________________________
void AliEmcalSetupTask::UserExec(Option_t *) 
{
  // Main loop, called for each event.

  if (fIsInit)
    return;

  AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
  if (!am) {
    AliError("Manager zero, returning");
    return;
  }
  am->LoadBranch("AliESDRun.");
  am->LoadBranch("AliESDHeader.");

  Int_t runno = InputEvent()->GetRunNumber();
  TString geoname("EMCAL_COMPLETE12SMV1");
  Int_t year = 2013;
  if (runno<=139517) {
    year = 2010;
    geoname = "EMCAL_FIRSTYEARV1";
  } else if ((runno>139517) && (runno<=170593)) {
    year = 2011;
    geoname = "EMCAL_COMPLETEV1";
  } else if ((runno>170593) && (runno<=193766)) {
    year = 2012;
  }

  AliEMCALGeometry *geom = AliEMCALGeometry::GetInstance(geoname);
  if (!geom) {
    AliFatal(Form("Can not create geometry: %s",geoname.Data()));
    return;
  }

  AliCDBManager *man = 0;
  man = AliCDBManager::Instance();
  if (!man->IsDefaultStorageSet()) {
    if (fOcdbPath.Length()>0) {
      AliInfo(Form("Setting up OCDB"));
      man->SetDefaultStorage(fOcdbPath);
      man->SetRun(runno);
    } else {
      man = 0;
    }
  } else {
    if (man->GetRun()!=runno)
      man->SetRun(runno);
  }
  
  if (man) {
    AliInfo(Form("Loading grp data from OCDB"));
    AliGRPManager GRPManager;
    GRPManager.ReadGRPEntry();
    GRPManager.SetMagField();
    AliInfo(Form("Loading geometry from OCDB"));
    AliGeomManager::LoadGeometry();
    AliGeomManager::ApplyAlignObjsFromCDB("GRP ITS TPC TRD EMCAL");
  }

  TGeoManager *geo = AliGeomManager::GetGeometry();
  if (!geo) {
    TString fname(gSystem->ExpandPathName(Form("%s/geometry_%d.root", fGeoPath.Data(), year)));
    if (gSystem->AccessPathName(fname)==0) {
      AliInfo(Form("Loading geometry from %s", fname.Data()));
      AliGeomManager::LoadGeometry(fname);
      geo = AliGeomManager::GetGeometry();
    }
  }
  if (geo) {
    AliInfo(Form("Locking geometry"));
    geo->LockGeometry();
  }

  if (!TGeoGlobalMagField::Instance()->GetField()) { // construct field map
    InputEvent()->InitMagneticField();
  }

  if (fOadbPath.Length()>0) {
    AliOADBContainer emcalgeoCont(Form("emcal"));
    emcalgeoCont.InitFromFile(Form("%s/EMCALlocal2master.root",fOadbPath.Data()),
                              Form("AliEMCALgeo"));
    TObjArray *mobj=dynamic_cast<TObjArray*>(emcalgeoCont.GetObject(runno,"EmcalMatrices"));
    if (mobj) {
      for(Int_t mod=0; mod < (geom->GetEMCGeometry())->GetNumberOfSuperModules(); mod++){
        //AliInfo(Form("Misalignment matrix %d", mod));
        geom->SetMisalMatrix((TGeoHMatrix*) mobj->At(mod),mod);
      } 
    }
  }

  fIsInit = kTRUE;
}
