#include "ARVersion.h"
#if !defined(__CINT__) || defined(__MAKECINT__)
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
#include "iostream"
#include "fstream"
using namespace std;
#endif

Bool_t DiffGeomBeforeTagging(const char* recipient, const char* cdbUri="local://$ALICE_ROOT/OCDB", const char* cfgFile="$ALICE_ROOT/macros/Config.C"){
	// Compare the geometry created in the current AliRoot with the one
	// in the OCDB directory of the current aliroot. If they differ,
        // warn the recipients to update the geometry in $ALICE_ROOT/OCDB
        // before tagging and to consider updating it in the raw OCDB.
        //
        // Running this macro always before creating a release tag will
        // allow to have in release tags OCDBs the geometry corresponding
        // to their code.
        //
        // The recipient string is supposed to contain comma-separated 
        // valid email addresses, at least the one of the librarian and
        // probably also the one of the guy in charge of the geometry object.
        //

        // Uplaod the geometry from $ALICE_ROOT/OCDB and get the number of nodes
	AliCDBManager* cdb = AliCDBManager::Instance();
	cdb->SetDefaultStorage(cdbUri);
        cdb->Get("GRP/Geometry/Data",0);
        Int_t svnocdb_nnodes = gGeoManager->GetNNodes();
        Printf("The geometry in OCDB has %d nodes",svnocdb_nnodes);
	cdb->SetRun(0);

	gGeoManager = 0;

        // Create the geometry and get the number of nodes
	if(!gSystem->AccessPathName("geometry.root")){
		Printf("Deleting existing \"geometry.root\"");
		gSystem->Exec("rm -rf geometry.root");
	}

	gROOT->LoadMacro(cfgFile);
	gInterpreter->ProcessLine(gAlice->GetConfigFunction());
	gAlice->GetMCApp()->Init();

	if(!gGeoManager){
		Printf("Unable to produce a valid geometry to be put in the CDB!");
		return kFALSE;
	}

	if(gSystem->AccessPathName("geometry.root")){
		Printf("Did not find freshly written \"geometry.root\" file. Exiting ...");
		return kFALSE;
	}

	Printf("Reloading freshly written geometry.root file");
        if (TGeoManager::IsLocked()) TGeoManager::UnlockGeometry();
        AliGeomManager::LoadGeometry("geometry.root");

        Int_t svncode_nnodes = gGeoManager->GetNNodes();
        Printf("The generated geometry has %d nodes",svncode_nnodes);

        // If the number of nodes differs, remove the geometry file from the local OCDB
        // and replace it by calling the UpdateCDBIdealGeom.C macro
        if (svncode_nnodes != svnocdb_nnodes){
                Printf("The geometry generated and the one in the current OCDB differ.");
                Printf("Remove \"$ALICE_ROOT/OCDB/GRP/Geometry/Data/Run0_999999999_v0_s0.root\".");
                Printf("and run $ALICE_ROOT/GRP/UpdateCDBIdealGeom.C(\"local://$ALICE_ROOT/OCDB\",%s)",cfgFile);

		// Get root and AliRoot versions
		const char* rootv = gROOT->GetVersion();
		TString av(ALIROOT_SVN_BRANCH);
		Int_t revnum = ALIROOT_SVN_REVISION;
		Printf("root version: %s.  AliRoot %s, revision number %d",rootv,av.Data(),revnum);

		// create and send mail
		TString subject(Form("AliRoot revision %d is about to be tagged and the geometry changed!",revnum));
                // mail body 
                TString bodyFileName("mailbody.txt");
               // bodyFileName.Form("%s/mail.body", GetShuttleLogDir());
                gSystem->ExpandPathName(bodyFileName);

                ofstream mailBody;
                mailBody.open(bodyFileName, ofstream::out);

                if (!mailBody.is_open())
                {
                        Printf("Could not open mail body file %s", bodyFileName.Data());
                        return kFALSE;
                }

                TString recipients(recipient);

                TString body;
                body = Form("Dear AliRoot librarian or geometry expert,\n\n");
                body += Form("while tagging the current revision - r%d - the geometry produced by the code\n"
                                "appears to have changed w.r.t. to the one saved at the previous tag.\n\n"
				"You should remove \"$ALICE_ROOT/OCDB/GRP/Geometry/Data/Run0_999999999_v0_s0.root\"\n"
				"and run $ALICE_ROOT/GRP/UpdateCDBIdealGeom.C(\"local://$ALICE_ROOT/OCDB\",%s) "
				"before tagging.\n\n"
                                "Also consider updating the raw geometry on the raw OCDB.\n\n"
				"Je vous prie de bien vouloir croire en l'assurance de mes respectueuses et honorables salutations.",revnum,cfgFile);
		mailBody << body.Data();
		mailBody.close();

		Printf("Sending mail to \"%s\"",recipients.Data());
		TString mailCommand("");
                if(recipients.CountChar(',')==0){
                        mailCommand = Form("mail -s \"%s\" %s < %s",
                                        subject.Data(),
                                        recipients.Data(),
                                        bodyFileName.Data());
                }else{
                        TString cc(recipients);
                        recipients.Remove(recipients.First(','));
                        cc.Replace(0,cc.First(',')+1,"");
                        mailCommand = Form("mail -s \"%s\" -c %s %s < %s",
                                        subject.Data(),
                                        cc.Data(),
                                        recipients.Data(),
                                        bodyFileName.Data());
                }

                Bool_t result = gSystem->Exec(mailCommand.Data());
		return result;
        }else{
                Printf("There are no changes between the geometry generated with current code and the one in the current OCDB.");
        }
	return kTRUE;
}

