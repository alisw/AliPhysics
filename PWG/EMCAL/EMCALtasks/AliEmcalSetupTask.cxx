/**************************************************************************
 * Copyright(c) 1998-2013, ALICE Experiment at CERN, All rights reserved. *
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

/// \cond CLASSIMP
ClassImp(AliEmcalSetupTask)
/// \endcond

/**
 * Constructor. Setting default values:
 * -# OCDB:           Local snapshot
 * -# OCDB objects:   GRP ITS TPC TRD EMCAL
 * -# OADB:           $ALICE_PHYSICS/OADB/EMCAL
 * -# Geometry:       $ALICE_PHYSICS/OADB/EMCAL
 */
AliEmcalSetupTask::AliEmcalSetupTask() : 
       AliAnalysisTaskSE(),
       fOcdbPath("uselocal"),
       fOadbPath("$ALICE_PHYSICS/OADB/EMCAL"),
       fGeoPath("$ALICE_PHYSICS/OADB/EMCAL"),
       fObjs("GRP ITS TPC TRD EMCAL"),
       fNoOCDB(kFALSE),
       fIsInit(kFALSE),
       fLocalOcdb(),
       fLocalOcdbStor()
{
}

/**
 * Named constructor. Setting default values:
 * -# OCDB:           Local snapshot
 * -# OCDB objects:   GRP ITS TPC TRD EMCAL
 * -# OADB:           $ALICE_PHYSICS/OADB/EMCAL
 * -# Geometry:       $ALICE_PHYSICS/OADB/EMCAL
 * @param name Name of the setup task
 */
AliEmcalSetupTask::AliEmcalSetupTask(const char *name) : 
      AliAnalysisTaskSE(name),
      fOcdbPath("uselocal"),
      fOadbPath("$ALICE_PHYSICS/OADB/EMCAL"),
      fGeoPath("$ALICE_PHYSICS/OADB/EMCAL"),
      fObjs("GRP ITS TPC TRD EMCAL"),
      fNoOCDB(kFALSE),
      fIsInit(kFALSE),
      fLocalOcdb(),
      fLocalOcdbStor()
{
  fBranchNames = "ESD:AliESDHeader.,AliESDRun.";
}

/**
 * Destructor
 */
AliEmcalSetupTask::~AliEmcalSetupTask()
{
}

/**
 *
 * @param option
 */
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

/**
 * Main loop, called for each event. Executed only for the first event.
 * In case databases are not initialized, run Setup.
 * Attention: The task relies cannot handle run changes.
 * @param
 */
void AliEmcalSetupTask::UserExec(Option_t *) 
{
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
/**
 * Setup databases:
 * -# Setting up geometry by run number
 * -# Setting up OCDB on request: If OCDB is local, extract the snapshot file and set the OCDB path to this, otherwise
 *    initialize the OCDB Manager with the OCDB path provided
 * -# Setting up EMCAL OABD containers from the path specified.
 * @param runno Run number obtained from the input event
 */
void AliEmcalSetupTask::Setup(Int_t runno) 
{
  // Setup everything

  if ( runno < 0 )  // Run number 0 can occur for MC
    return;
  
  AliEMCALGeometry *geom = AliEMCALGeometry::GetInstanceFromRunNumber(runno);
  
  if (!geom) {
    AliFatal("Can not create geometry!");
    return;
  }
  
  // Setup AliEMCALGeometry corresponding to year
  Int_t year = 2013;
  if (runno>0 && runno<=139517) {
    year = 2010;
  } else if ((runno > 139517) && (runno <= 170593)) {
    year = 2011;
  } else if ((runno > 170593) && (runno <= 193766)) {
    year = 2012;
  } else if ((runno > 193766) && (runno <= 199161)) {
    year = 2013;
  } else if ( runno > 199161) { //MV: is this the last run of run 1?
    year = 2015;
  }
  
  // Setup CDB manager
  AliCDBManager *man = 0;
  if (!fNoOCDB) {
    man = AliCDBManager::Instance();
    if (!man)
      AliFatal(Form("Did not get pointer to CDB manager"));

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
        TString filename(Form("$ALICE_PHYSICS/PWG/EMCAL/data/%d.dat",year));
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
  }

  // Load geometry from OCDB 
  TGeoManager *geo(nullptr);
  if (man) {
    if (man->GetRun()!=runno)
      man->SetRun(runno);
    AliInfo(Form("Loading grp data from OCDB for run %d", runno));
    AliGRPManager GRPManager;
    GRPManager.ReadGRPEntry();
    GRPManager.SetMagField();
    if(!AliGeomManager::GetGeometry() && man){
      AliInfo(Form("Loading geometry from OCDB"));
      AliGeomManager::LoadGeometry();
      geo = AliGeomManager::GetGeometry();
      if (!fObjs.IsNull())
        AliGeomManager::ApplyAlignObjsFromCDB(fObjs);
      // Lock geometry
      if (geo) {
        AliInfo(Form("Locking geometry"));
        geo->LockGeometry();
      }
    }
  }

  // Load geometry from file (does not use misalignment of ITS/TPC!)
  geo = AliGeomManager::GetGeometry();
  if (!geo && fGeoPath.Length()>0) {
    TString fname(gSystem->ExpandPathName(Form("%s/geometry_%d.root", fGeoPath.Data(), year)));
    if (gSystem->AccessPathName(fname)==0) {
      AliInfo(Form("Loading geometry from file %s (should be avoided!)", fname.Data()));
      AliGeomManager::LoadGeometry(fname);
      geo = AliGeomManager::GetGeometry();
      // Lock geometry
      if (geo) {
        AliInfo(Form("Locking geometry"));
        geo->LockGeometry();
      }
    }
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
        if(mobj->At(mod))
        {
          geom->SetMisalMatrix((TGeoHMatrix*) mobj->At(mod),mod);
        }
        else if(gGeoManager)
        {
          AliWarning(Form("Set matrix for SM %d from gGeoManager\n",mod));
          geom->SetMisalMatrix(geom->GetMatrixForSuperModuleFromGeoManager(mod),mod);
        }
        else
        {
          AliError(Form("Matrix for SM %d is not available\n",mod));
        }
      } 
    }
  }

  fIsInit = kTRUE;
}

/**
 * Terminate function, called at the end of the analysis. Cleaning up the local OCDB snapshot (if created).
 * @param
 */
void AliEmcalSetupTask::Terminate(Option_t *) 
{
  if (fLocalOcdb.Length()>0) {
    TString cmd(Form("rm -rf %s", fLocalOcdb.Data()));
    gSystem->Exec(cmd);
  }
}
