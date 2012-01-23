// $Id$
//
// Task to setup emcal related global objects
//
//

#include <TClonesArray.h>
#include <TGeoManager.h>
#include <TGeoGlobalMagField.h>
#include "AliAnalysisManager.h"
#include "AliCDBManager.h"
#include "AliEMCALGeometry.h"
#include "AliESDEvent.h"
#include "AliEmcalSetupTask.h"
#include "AliGeomManager.h"
#include "AliOADBContainer.h"

ClassImp(AliEmcalSetupTask)

//________________________________________________________________________
AliEmcalSetupTask::AliEmcalSetupTask() : 
  AliAnalysisTaskSE(),
  fOcdbPath(),
  fOadbPath("$ALICE_ROOT/OADB/EMCAL"),
  fGeoPath("."),
  fEsdEv(0),
  fIsInit(kFALSE)
{
  // Constructor.
}

//________________________________________________________________________
AliEmcalSetupTask::AliEmcalSetupTask(const char *name) : 
  AliAnalysisTaskSE(name),
  fOcdbPath(),
  fOadbPath("$ALICE_ROOT/OADB/EMCAL"),
  fGeoPath("."),
  fEsdEv(0),
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

  fEsdEv = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!fEsdEv) {
    AliError("Task works only on ESD events, returning");
    return;
  }

  AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
  if (!am) {
    AliError("Manager zero, returning");
    return;
  }
  am->LoadBranch("AliESDRun.");
  am->LoadBranch("AliESDHeader.");

  Int_t runno = fEsdEv->GetRunNumber();
  Bool_t is2010 = kTRUE;
  if (runno>139517) {
    is2010 = kFALSE;
  }

  TString geoname("EMCAL_FIRSTYEARV1");
  if (!is2010)
    geoname = "EMCAL_COMPLETEV1";
  AliEMCALGeometry *geom = AliEMCALGeometry::GetInstance(geoname);
  if (!geom) {
    AliFatal(Form("Can not create geometry: %s",geoname.Data()));
    return;
  }

  AliCDBManager *man = 0;
  if (fOcdbPath.Length()>0) {
    AliInfo(Form("Setting up OCDB"));
    man = AliCDBManager::Instance();
    man->SetDefaultStorage(fOcdbPath);
    man->SetRun(runno);
  }

  TGeoManager *geo = AliGeomManager::GetGeometry();
  if (!geo) {
    TString fname(Form("%s/geometry_2010.root", fGeoPath.Data()));
    if (!is2010)
      fname = Form("%s/geometry_2011.root", fGeoPath.Data());
    if (gSystem->AccessPathName(fname)!=0)
      fname = Form("%s/geometry.root", fGeoPath.Data());
    if (gSystem->AccessPathName(fname)==0) {
      AliInfo(Form("Loading geometry from %s", fname.Data()));
      AliGeomManager::LoadGeometry(fname);
    } else if (man) {
      AliInfo(Form("Loading geometry from OCDB"));
      AliGeomManager::LoadGeometry();
    }
  }
  if (geo) {
    AliGeomManager::ApplyAlignObjsFromCDB("EMCAL");
    AliInfo(Form("Locking geometry"));
    geo->LockGeometry();
  }

  if (!TGeoGlobalMagField::Instance()->GetField()) { // construct field map
    AliInfo("Constructing field map from ESD run info");
    fEsdEv->InitMagneticField();
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
