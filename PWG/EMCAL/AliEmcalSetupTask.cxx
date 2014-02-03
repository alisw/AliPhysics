// $Id$
//
// Task to setup emcal related global objects.
//
// Author: C.Loizides

#include "AliEmcalSetupTask.h"
#include <TClonesArray.h>
#include <TGeoGlobalMagField.h>
#include <TGeoManager.h>
#include <TRandom.h>
#include "AliAODEvent.h"
#include "AliAnalysisManager.h"
#include "AliCDBManager.h"
#include "AliEMCALGeometry.h"
#include "AliESDEvent.h"
#include "AliGRPManager.h"
#include "AliGeomManager.h"
#include "AliMagF.h"
#include "AliOADBContainer.h"
#include "AliTender.h"

ClassImp(AliEmcalSetupTask)

//________________________________________________________________________
AliEmcalSetupTask::AliEmcalSetupTask() : 
  AliAnalysisTaskSE(),
  fOcdbPath("uselocal"),
  fOadbPath("$ALICE_ROOT/OADB/EMCAL"),
  fGeoPath("$ALICE_ROOT/OADB/EMCAL"),
  fObjs("GRP ITS TPC TRD EMCAL"),
  fIsInit(kFALSE),
  fLocalOcdb(),
  fLocalOcdbStor()
{
  // Constructor.
}

//________________________________________________________________________
AliEmcalSetupTask::AliEmcalSetupTask(const char *name) : 
  AliAnalysisTaskSE(name),
  fOcdbPath("uselocal"),
  fOadbPath("$ALICE_ROOT/OADB/EMCAL"),
  fGeoPath("$ALICE_ROOT/OADB/EMCAL"),
  fObjs("GRP ITS TPC TRD EMCAL"),
  fIsInit(kFALSE),
  fLocalOcdb(),
  fLocalOcdbStor()
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
void AliEmcalSetupTask::ConnectInputData(Option_t *option)
{
  // Connect input data

  AliAnalysisTaskSE::ConnectInputData(option);

  if (fOcdbPath.Length()==0)
    return;

  AliCDBManager *man = AliCDBManager::Instance();
  if (man->IsDefaultStorageSet()) 
    return;

  if (fIsInit)
    return;

  AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
  if (!am)
    return;

  TObjArray *tasks = am->GetTasks();
  if (!tasks)
    return;

  AliTender *tender = 0;
  for (Int_t i=0; i<tasks->GetEntries(); ++i) {
    tender = dynamic_cast<AliTender*>(tasks->At(i));
    if (tender)
      break;
  }

  if (!tender)
    return;

  if (fOcdbPath != "uselocal") {
    tender->SetDefaultCDBStorage(fOcdbPath);
    return;
  }

  Int_t runno = AliAnalysisManager::GetAnalysisManager()->GetRunFromPath();
  if (runno<=0) {
    AliWarning(Form("Disabling tender, ignore possible message from tender below"));
    tender->SetDefaultCDBStorage("donotuse");
    return;
  }

  AliWarning(Form("Intercepting tender for run %d, ignore possible message from tender below", runno));
  Setup(runno);
  tender->SetDefaultCDBStorage(fLocalOcdbStor);
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
  Setup(runno);
}

//________________________________________________________________________
void AliEmcalSetupTask::Setup(Int_t runno) 
{
  // Setup everything

  // Setup AliEMCALGeometry corresponding to year
  TString geoname("EMCAL_COMPLETE12SMV1");
  Int_t year = 2013;
  if (runno>0 && runno<=139517) {
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

  if (runno<=0)
    return;

  // Setup CDB manager
  AliCDBManager *man = 0;
  man = AliCDBManager::Instance();
  if (man->IsDefaultStorageSet()) {
    AliInfo(Form("Default OCDB storage already set"));
  } else {
    if (fOcdbPath.Length()==0) {
      man = 0; // do not use OCDB
    } else if (fOcdbPath != "uselocal") {
      AliInfo(Form("Setting up OCDB to point to %s",fOcdbPath.Data()));
      man->SetDefaultStorage(fOcdbPath);
    } else { // use local copy of OCDB
      TString tmpdir=gSystem->WorkingDirectory();
      if (gSystem->AccessPathName(tmpdir))
	tmpdir = "/tmp";
      tmpdir+="/";
      tmpdir+=gSystem->GetUid();
      tmpdir+="-";
      TDatime t;
      tmpdir+=t.Get();
      tmpdir+="-";
      Int_t counter = 0;
      fLocalOcdb = tmpdir;
      fLocalOcdb += Form("%d%d%d",gRandom->Integer(999999999),gRandom->Integer(999999999),gRandom->Integer(999999999));
      while (!gSystem->AccessPathName(fLocalOcdb)) {
	fLocalOcdb = tmpdir;
	fLocalOcdb += Form("%d%d%d",gRandom->Integer(999999999),gRandom->Integer(999999999),gRandom->Integer(999999999));
	counter++;
	if (counter>100) {
	  AliFatal(Form("Could not create local directory for OCDB at %s",tmpdir.Data()));
	}
      }
      gSystem->MakeDirectory(fLocalOcdb);
      TString filename(Form("$ALICE_ROOT/PWG/EMCAL/data/%d.dat",year));
      TString cmd(Form("cd %s && tar -xf %s",fLocalOcdb.Data(),filename.Data()));
      Int_t ret = gSystem->Exec(cmd);
      if (ret==0) {
	TString locdb("local://");
	locdb+=fLocalOcdb;
	locdb+="/";
	locdb+=year;
	AliInfo(Form("Setting up local OCDB at %s",locdb.Data()));
	man->SetDefaultStorage(locdb);
	fLocalOcdbStor = locdb;
      } else {
	AliFatal(Form("Could not set up local OCDB at %s",fLocalOcdb.Data()));
      }
    }
  }
  
  // Load geometry from OCDB 
  if (man) {
    if (man->GetRun()!=runno)
      man->SetRun(runno);
    AliInfo(Form("Loading grp data from OCDB for run %d", runno));
    AliGRPManager GRPManager;
    GRPManager.ReadGRPEntry();
    GRPManager.SetMagField();
    AliInfo(Form("Loading geometry from OCDB"));
    AliGeomManager::LoadGeometry();
    if (!fObjs.IsNull())
      AliGeomManager::ApplyAlignObjsFromCDB(fObjs);
  }

  // Load geometry from file (does not use misalignment of ITS/TPC!)
  TGeoManager *geo = AliGeomManager::GetGeometry();
  if (!geo && fGeoPath.Length()>0) {
    TString fname(gSystem->ExpandPathName(Form("%s/geometry_%d.root", fGeoPath.Data(), year)));
    if (gSystem->AccessPathName(fname)==0) {
      AliInfo(Form("Loading geometry from file %s (should be avoided!)", fname.Data()));
      AliGeomManager::LoadGeometry(fname);
      geo = AliGeomManager::GetGeometry();
    }
  }

  // Lock geometry
  if (geo) {
    AliInfo(Form("Locking geometry"));
    geo->LockGeometry();
  }

  // Construct field map
  if (!TGeoGlobalMagField::Instance()->GetField()) { 
    InputEvent()->InitMagneticField();
  }

  // Apply mis-alignment matrices from OADB
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

//________________________________________________________________________
void AliEmcalSetupTask::Terminate(Option_t *) 
{
  // Called at the end.

  if (fLocalOcdb.Length()>0) {
    TString cmd(Form("rm -rf %s", fLocalOcdb.Data()));
    gSystem->Exec(cmd);
  }
}
