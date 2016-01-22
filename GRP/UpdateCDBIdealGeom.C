#if !defined(__CINT__) || defined(__MAKECINT__)
#include "ARVersion.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBId.h"
#include "AliCDBMetaData.h"
#include "AliGeomManager.h"
#include "AliMC.h"
#include <TROOT.h>
#include "AliRun.h"
#include <TGeoManager.h>
#include <TString.h>
#include <TInterpreter.h>
#endif

void UpdateCDBIdealGeom(const char* cdbUri, const char* cfgFile){
	// Produce the ideal geometry and store it in the specified CDB
	// The second argument allows to specify the config file to be used
	// in particular for giving the choice to generate either a full or
	// a partial geometry.
	//

	AliCDBManager* cdb = AliCDBManager::Instance();
	// we set the default storage to the repository because some dets require
	// already at the time of geometry creation to find calibration objects in the cdb
	AliCDBStorage* storage = 0;
	if(!cdb->IsDefaultStorageSet()) cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
	storage = cdb->GetStorage(cdbUri);
	if(!storage) 
	{
		Printf("unable to create valid storage from: %s", cdbUri);
		return;
	}
	cdb->SetRun(0);
	AliCDBId id("GRP/Geometry/Data",0,AliCDBRunRange::Infinity());
	AliCDBMetaData *md= new AliCDBMetaData();

	// Get root and AliRoot versions
	const char* rootv = gROOT->GetVersion();
	TString av(ALIROOT_VERSION);
	TString revnum(ALIROOT_REVISION);

	Printf("root version: %s.  AliRoot %s, revision number %s",rootv,av.Data(),revnum.Data());

	md->SetAliRootVersion(av.Data());
	md->SetComment(Form("Geometry produced with root version %s and AliRoot %s, revision number %s",rootv,av.Data(),revnum.Data()));
	md->AddDateToComment();

	//gSystem->Exec("if [ -e geometry.root ]; then \necho deleting existing geometry.root \nrm -rf geometry.root \nfi");
	if(!gSystem->AccessPathName("geometry.root")){
		Printf("Deleting existing \"geometry.root\"");
		gSystem->Exec("rm -rf geometry.root");
	}
	gSystem->Load("libgeant321");
	gSystem->Load("libqpythia");
	gSystem->Load("libAliPythia6");
	gROOT->LoadMacro(cfgFile);
	gInterpreter->ProcessLine(gAlice->GetConfigFunction());

	gAlice->GetMCApp()->Init();

	if(!gGeoManager){
		Printf("Unable to produce a valid geometry to be put in the CDB!");
		return;
	}

	/*
	if (gSystem->Exec("if [ ! -e geometry.root ]; then \n return 1  \nfi")) {
		Printf("Did not find freshly written geometry.root file");
		return;
	} */
	if(gSystem->AccessPathName("geometry.root")){
		Printf("Did not find freshly written \"geometry.root\" file. Exiting ...");
		return;
	}

	Printf("Reloading freshly written geometry.root file");
	if (TGeoManager::IsLocked()) TGeoManager::UnlockGeometry();
	AliGeomManager::LoadGeometry("geometry.root");

	gGeoManager->DefaultColors(); // assign default colors according to Z of material
	// (many colors turned into dark gray nuances some time ago, when the root palette was changed)

	Printf("Storing in CDB geometry produced with root version %s and AliRoot version %s",rootv,av.Data());
	storage->Put(gGeoManager,id,md);
	// This is to allow macros lauched after this one in the same session to find the
	// newly produced geometry.
	storage->QueryCDB(cdb->GetRun());

}


