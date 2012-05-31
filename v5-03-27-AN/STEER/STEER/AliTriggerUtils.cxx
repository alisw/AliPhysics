/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

#include "AliTriggerUtils.h"
#include "AliTriggerConfiguration.h"
#include "AliTriggerClass.h"
#include "AliTriggerInput.h"
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliCDBManager.h"
#include "AliPDG.h"
#include "AliMC.h"
#include "AliModule.h"
#include "AliLog.h"
#include <TROOT.h>
#include <TInterpreter.h>
#include <TString.h>

ClassImp(AliTriggerUtils)

//_____________________________________________________________________________
Bool_t AliTriggerUtils::CheckConfiguration( TString& configfile, AliTriggerConfiguration * cfg )
{
   // To be used on the pre-creation of Configurations to check if the
   // conditions have valid inputs names.
   //
   // Initiate detectors modules from a Config file
   // Ask to each active module present in the fDetectorCluster
   // to create a Trigger detector and retrive the inputs from it
   // to create a list of inputs.
   // Each condition in the configuration is then checked agains 
   // the list of inputs


   if (!gAlice) {
      AliError( "no gAlice object. Restart aliroot and try again." );
      return kFALSE;
   }
   if (gAlice->Modules()->GetEntries() > 0) {
      AliError( "gAlice was already run. Restart aliroot and try again." );
      return kFALSE;
   }

   AliInfo( Form( "initializing gAlice with config file %s",
            configfile.Data() ) );
//_______________________________________________________________________
   gAlice->Announce();
   
   gROOT->LoadMacro(configfile.Data());
   gInterpreter->ProcessLine(gAlice->GetConfigFunction());
   
   if(AliCDBManager::Instance()->GetRun() >= 0) { 
     AliRunLoader::Instance()->SetRunNumber(AliCDBManager::Instance()->GetRun());
   } else {
     AliWarning("Run number not initialized!!");
   }
  
   AliRunLoader::Instance()->CdGAFile();
    
   AliPDG::AddParticlesToPdgDataBase();  

   gAlice->GetMCApp()->Init();
   
   //Must be here because some MCs (G4) adds detectors here and not in Config.C
   gAlice->InitLoaders();
   AliRunLoader::Instance()->MakeTree("E");
   AliRunLoader::Instance()->LoadKinematics("RECREATE");
   AliRunLoader::Instance()->LoadTrackRefs("RECREATE");
   AliRunLoader::Instance()->LoadHits("all","RECREATE");
   //
   // Save stuff at the beginning of the file to avoid file corruption
   AliRunLoader::Instance()->CdGAFile();
   gAlice->Write();

   AliRunLoader* runLoader = AliRunLoader::Instance();
   if( !runLoader ) {
      AliError( Form( "gAlice has no run loader object. "
                      "Check your config file: %s", configfile.Data() ) );
      return kFALSE;
   }

   // get the possible inputs to check the condition
   TObjArray inputs;
   TObjArray* detArray = runLoader->GetAliRun()->Detectors();

   TString detStr = cfg->GetTriggeringModules();

   for( Int_t iDet = 0; iDet < detArray->GetEntriesFast(); iDet++ ) {
      AliModule* det = (AliModule*) detArray->At(iDet);
      if( !det || !det->IsActive() ) continue;
      if( cfg->IsSelected( det->GetName(), detStr ) ) {
         AliInfo( Form( "Creating inputs for %s", det->GetName() ) );
         AliTriggerDetector* dtrg = det->CreateTriggerDetector();
         dtrg->AssignInputs(cfg->GetInputs());
         TObjArray* detInp = dtrg->GetInputs();
         for( Int_t i=0; i<detInp->GetEntriesFast(); i++ ) {
            AliInfo( Form( "Adding input %s", ((AliTriggerInput*)detInp->At(i))->GetName() ) );
            inputs.AddLast( detInp->At(i) );
         }
      }
   }

   // check if the condition is compatible with the triggers inputs
   Int_t ndesc = cfg->GetClasses().GetEntriesFast();
   Bool_t check = kTRUE;
   ULong64_t mask = 0L;
   for( Int_t j=0; j<ndesc; j++ ) {
     AliTriggerClass *trclass = (AliTriggerClass*)cfg->GetClasses().At( j );
     if( !(trclass->CheckClass( cfg )) ) check = kFALSE;
     else {
       if (trclass->IsActive(cfg->GetInputs(),cfg->GetFunctions())) {
	 AliInfo( Form( "Trigger Class (%s) OK, class mask (0x%llx)",
			trclass->GetName(), trclass->GetMask( ) ) );
       }
       else {
	 AliWarning( Form( "Trigger Class (%s) is NOT active, class mask (0x%llx)",
			   trclass->GetName(), trclass->GetMask( ) ) );
       }
     }
     // check if condition mask is duplicated
     if( mask & trclass->GetMask() ) {
       AliError( Form("Class (%s). The class mask (0x%llx) is ambiguous. It was already defined",
		      trclass->GetName(), trclass->GetMask()  ) );
       check = kFALSE;
     }
     mask |= trclass->GetMask();
   }

   return check;
}
