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
  fGeoPath("."),
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
  TString geoname("EMCAL_FIRSTYEARV1");
  Int_t year = 2010;
  if (runno>139517) {
    year = 2011;
    geoname = "EMCAL_COMPLETEV1";
  } 
  if (runno>170593) {
    year = 2012;
    geoname = "EMCAL_COMPLETE12SMV1";
  }

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
    TString fname(gSystem->ExpandPathName(Form("%s/geometry_%d.root", fGeoPath.Data(), year)));
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
    AliESDEvent *esdEv = dynamic_cast<AliESDEvent*>(InputEvent());
    if (esdEv) {
      AliInfo("Constructing field map from ESD run info");
      esdEv->InitMagneticField();
    } else {
      AliAODEvent *aodEv = dynamic_cast<AliAODEvent*>(InputEvent());
      if (aodEv) {
        Double_t curSol = 30000*aodEv->GetMagneticField()/5.00668;
        Double_t curDip = 6000 *aodEv->GetMuonMagFieldScale();
        AliMagF *field  = AliMagF::CreateFieldMap(curSol,curDip);
        TGeoGlobalMagField::Instance()->SetField(field);
      }
    }
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
